
args<-commandArgs(trailingOnly = TRUE)

amp_file <-  args[1] ; stopifnot(!is.na(amp_file))
bam_path <- args[2] ; stopifnot(!is.na(bam_path))
output_dir <- args[3] ; stopifnot(!is.na(output_dir))
bedfile <- args[4] ; stopifnot(!is.na(bedfile))
bin <- as.numeric(args[5]) ; stopifnot(!is.na(bin))
ref_genome <- args[6] ; stopifnot(!is.na(ref_genome))
assembly <- args[7] ; stopifnot(!is.na(assembly))
target_file <- args[8]; stopifnot(!is.na(target_file))
name <- ifelse(is.na(args[9]),"sample",args[9])
ploidy <- ifelse(is.na(args[10]),2,as.numeric(args[10])) 
chr_list <- args[11] #17 or 17,18 by example

library("BiocManager")
library("ChIPpeakAnno")
library("ACE")
library("GenomicRanges")
library("QDNAseq")
library("QDNAseq.hg19")
library("Rsamtools")
library("ggplot2")
library("plyr")
library("ggrepel")
library("ggpubr")

print(sessionInfo())
print(R.Version()[c("major","minor")])

set.seed(1)
options(scipen=999)

nano_cnvplot <- function(template, bedGR, cnv_conf,cellularity,ploidy,title,type=NULL)
{
    gene_alter_color <- c("grey60","darkorange","seagreen3","tomato3")
    names(gene_alter_color) <- c("19","2","3","segment")
    max_y <- 20
    max_y_epic <- 6

    if(!is.null(type) & type=="EPIC"){max_y <- max_y_epic}

    template[which(template$adjustedcopynumbers<0),"adjustedcopynumbers"] <- 0
    template[which(template$adjustedsegments<0),"adjustedsegments"] <- 0

    template$chr <- gsub("chr", "", template$chr, ignore.case = TRUE)
    template$chr <- gsub("X","23",template$chr)
    template <- template[which(!template$chr%in%c("Y")),]
    template <- template[order(template$chr,template$start),]

    #Overlap segment with genes

    tabGR <- GRanges(seqnames=template[,"chr"],IRanges(start=template[,"start"],end=template[,"end"]))
    ov <- suppressMessages(findOverlapsOfPeaks(tabGR, bedGR))

    tab2 <- data.frame(ov$overlappingPeaks)
    tab2$merge <- paste(tab2[,2],tab2[,3],tab2[,4],sep="_")
    tabg <- tapply(tab2[,13],tab2$merge,function(x){paste(unique(x),collapse="/")})

    template$merge <- paste(template[,"chr"],template[,"start"],template[,"end"],sep="_")
    template$Genes <- mapvalues(template$merge,names(tabg),tabg)
    template[grep("_",template$Genes),"Genes"] <- NA

    template$chr <- factor(template$chr,levels=sort(as.numeric(unique(template$chr))))
    template <- template[order(template$chr),]
    template$bin <- 1:nrow(template)

    chrcol <- 1 + as.numeric(template$chr)- 2 * floor(as.numeric(template$chr)/2) #num.mark = number of snps in a segment
    nn <- cumsum(table(template$chr))
    segbdry <- cumsum(c(0, table(template$chr)))
    segstart <- segbdry[-length(segbdry)]
    segend <- segbdry[-1]

    chromlevels_labels <- gsub(23,"X",unique(template$chr))

    # Segments

    segment <- Reduce(function(x,y) rbind(x,y),lapply(as.character(unique(template$chr)),function(chrom){
        tab <- template[which(as.character(template$chr)==chrom),]
        tabna <- tab[which(!is.na(tab$adjustedsegments)),]
        # tab[which(is.na(tab$adjustedsegments_max)),"adjustedsegments_max"] <- -1
        rl <- rle(as.vector(na.omit(tab$adjustedsegments)))
        u <- rl$lengths
        dd <- data.frame(chr=chrom,s=c(1,cumsum(u)[-1]-u[-1]+1),e=cumsum(u),seg=rl$values)
        dd <- as.data.frame(cbind(dd,t(apply(dd,1,function(x){
        start <- tabna[as.numeric(x["s"]),"bin"]
        end <- tabna[as.numeric(x["e"]),"bin"]
        return(c(start,end))}))))
        return(dd)}),accumulate=FALSE)
    colnames(segment) <- c("chr","startx","endx","adjustedsegments","start","end")

    # Cap the y axis

    template$adjustedcopynumbers_max <- template$adjustedcopynumbers
    template[which(template$adjustedcopynumbers>max_y),"adjustedcopynumbers_max"] <- max_y

    template$adjustedsegments_max <- template$adjustedsegments
    template[which(template$adjustedsegments>max_y),"adjustedsegments_max"] <- max_y

    segment$adjustedsegments_max <- segment$adjustedsegments
    segment[which(segment$adjustedsegments>max_y),"adjustedsegments_max"] <- max_y

    # Alternatively color genes

    if(length(table(template$Genes))!=0){
        template$chrcol <- mapvalues(template$Genes, unique(template[which(!is.na(template$Genes)),"Genes"]) , rep(c(2,3),ceiling(length(unique(template[which(!is.na(template$Genes)),"Genes"]))/2))[1:length(unique(template[which(!is.na(template$Genes)),"Genes"]))])
    }

    template[which(!template$chrcol%in%c(2,3)),"chrcol"] <- 19
    template$chrcol <- as.character(template$chrcol)

    dk <- data.frame(s=c(0,cumsum(table(template$chr))[-length(cumsum(table(template$chr)))]),e=cumsum(table(template$chr))-1)
    dk$color <- "white"
    dk[seq(1,nrow(dk),by=2),"color"] <- "grey85"

    g <- ggplot(template)+
            geom_rect(data=dk,aes(xmin=s,xmax=e,ymin=-Inf,ymax=+Inf,alpha=0.315,fill=color))+
            scale_fill_manual(values=unique(dk$color)) +
            geom_point(size=0.75,aes(x=bin,y=adjustedcopynumbers_max,color=chrcol))+ theme_minimal() +
            # geom_point(size=point_size,aes(x=bin,y=adjustedcopynumbers_max,color=chrcol))+ theme_minimal() +
            theme(legend.position="none",panel.grid=element_blank(),panel.border = element_rect(colour = "black", fill=NA, size=0.25), plot.margin=unit(c(1,1,1,1), "cm"),axis.text.y=element_text(size=6)) +
            labs(x="",y="Copy Number") +
            scale_x_continuous(breaks=(nn + c(0, nn[-length(nn)]))/2, labels=chromlevels_labels,expand=c(0.0015,0.0015))    +
            geom_hline(yintercept=0,col="grey85",lty=1) +
            geom_vline(xintercept=c(0,cumsum(table(template$chr))),col="grey80",lty=2,size=0.5) +
            coord_cartesian(ylim=c(0,round(1.05*max_y))) +
            geom_hline(yintercept=ploidy,lty=6,color="grey60")+
            geom_segment(data=segment,aes(x=start,xend=end,y=adjustedsegments_max,yend=adjustedsegments_max,colour="segment"),alpha=0.8,lwd=1.5) +
            scale_y_continuous(breaks=seq(0,round(1.05*max_y,1)))+
            scale_color_manual(values=c(gene_alter_color)) +
            geom_vline(xintercept=template[which(template$chrcol==2),"bin"],col=gene_alter_color["2"],alpha=0.4) +
            geom_vline(xintercept=template[which(template$chrcol==3),"bin"],col=gene_alter_color["3"],alpha=0.4)+
            ggtitle(paste(title," (cellularity=",cellularity,", ploidy=",ploidy,")",sep=""))

    # genes, plot at the average of genes points the label
    template$index_lab <- as.numeric(as.character(mapvalues(template$Genes,names(round(tapply(template$bin,template$Genes,median))),round(tapply(template$bin,template$Genes,median)))))

    template$size <- 0.01
    template[which(template$chrcol!="19"),"size"] <- 0.02

    jseg_lab1 <- template[which(template$chrcol=="2"),] ; jseg_lab1 <- jseg_lab1[which(duplicated(jseg_lab1$Genes)==FALSE),]

    jseg_lab2 <- template[which(template$chrcol=="3"),] ; jseg_lab2 <- jseg_lab2[which(duplicated(jseg_lab2$Genes)==FALSE),]

    glab <- ggplot(template,aes(x=bin,color=chrcol))+ geom_hline(yintercept=0,lty=2,col="grey80")  +
            geom_vline(xintercept=c(0,cumsum(table(template$chr))),col="grey80",lty=2,size=0.5)

        jtempN <- template[which(template$chrcol!="19"),]
        if(nrow(jtempN)!=0){glab <- glab + geom_point(data=jtempN,aes(x=bin,y=0,color=chrcol,size=size))}
        if(nrow(jseg_lab1)!=0){ glab <- glab + geom_label_repel(data=jseg_lab1,aes(label=Genes,y=1,x=index_lab),nudge_y=3,direction="both",size=2.75,segment.alpha=0.5,segment.size = 0.3,fontface = 'bold')}
        if(nrow(jseg_lab2)!=0){glab <- glab + geom_label_repel(data=jseg_lab2,aes(label=Genes,y=-1,x=index_lab),nudge_y=-3,direction="both",size=2.75,segment.alpha=0.5,segment.size = 0.3,fontface = 'bold') }
        #if(nrow(jseg_lab1)!=0){ glab <- glab + geom_label_repel(data=jseg_lab1,aes(label=Genes,y=1,x=index_lab),nudge_y=nudge,direction=segdir,size=size_label,segment.alpha=alphaval,segment.size = segmentval,fontface = 'bold')}
        #if(nrow(jseg_lab2)!=0){glab <- glab + geom_label_repel(data=jseg_lab2,aes(label=Genes,y=-1,x=index_lab),nudge_y=-nudge,direction=segdir,size=size_label,segment.alpha=alphaval,segment.size = segmentval,fontface = 'bold') }

        glab <- glab +theme_minimal() + labs(x="",y="") +
            coord_cartesian(xlim=c(0,max(template$bin)),ylim=c(-10,10)) +
            scale_x_continuous(expand=c(0.0015,0.0015)) +
            scale_color_manual(values=gene_alter_color) +
            theme(legend.position="none",axis.text.y=element_blank(),axis.text.x=element_blank(),panel.grid=element_blank(),  plot.margin=unit(c(0,1,-0.5,1), "cm")) +
            scale_size_continuous(range=c(0,2))

    pComb <- ggarrange(g,glab,ncol=1,align="v",nrow=2,heights=c(0.8,0.2))
    return(pComb)
}

## read length & gc content

ref_fasta <- scanFa(ref_genome)
gc <- letterFrequency(ref_fasta, letters="CG") / width(ref_fasta)
ref <- data.frame(quality=assembly, GC_content = round(as.numeric(gc)*100,2), readlength = unique(width(ref_fasta)))
ref <- ref[which(ref$GC_content>0.05),]

bam <- scanBam(bam_path, what=bamWhat("mapq","seq"))
nb_reads <- length(unique(bam[[1]]$qname))

bam_gc <- letterFrequency(bam[[1]]$seq, letters="CG") / width(bam[[1]]$seq)
sample <- data.frame(quality=name,GC_content = round(as.numeric(bam_gc)*100,2), readlength = width(bam[[1]]$seq))

all <- rbind(sample[!is.nan(sample$GC_content),],ref)

sum_length <- data.frame(name=name,nb_reads=nb_reads,min=min(sample$readlength),mean=round(mean(sample$readlength)),max=max(sample$readlength))

write.table(sum_length,file.path(output_dir, paste0(name, "_readlength_summary.txt")),quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)

### plot GC content histogram
GC_plot <- ggplot(all, aes(x=GC_content,group=quality,colour=quality)) +
  geom_freqpoly(aes(y=..ncount..)) +
  theme_classic() +
  scale_colour_brewer(type = "qual", palette = "Set1") +
  theme(aspect.ratio=1) +
  xlab("GC content") +
  ylab("no. of mapped reads")

ggsave(plot=GC_plot,filename=file.path(output_dir, paste0(name, "_gc.pdf")), width=4, height=3)

### plot read length distribution

readlength_plot <- ggplot(sample, aes(x=readlength)) +
  geom_freqpoly(aes(y=..ncount..)) +
  theme_classic() +
  scale_x_log10() +
  theme(aspect.ratio=1) +
  xlab("read length (bp)") +
  ylab("no. of mapped reads")

ggsave(plot=readlength_plot,filename=file.path(output_dir, paste0(name, "_readlength.pdf")), width=3, height=3)

## genomic profile
bins <- getBinAnnotations(binSize = bin, genome = assembly)

readCounts <- binReadCounts(bins, bamfiles = bam_path)
readCountsFiltered <- applyFilters(readCounts,residual = TRUE, blacklist = TRUE,chromosomes="Y")
readCountsFiltered <- estimateCorrection(readCountsFiltered)
copyNumbers <- correctBins(readCountsFiltered)
copyNumbers <- normalizeBins(copyNumbers)
copyNumbers <- smoothOutlierBins(copyNumbers)
copyNumbersSegmented <- segmentBins(copyNumbers, transformFun = "sqrt")

copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
sample_nb <- 1

# extract model (cellularity + error)
model1 <- singlemodel(copyNumbersSegmented, QDNAseqobjectsample = sample_nb, exclude=c())
cellularity <- model1$minima[tail(which(model1$rerror==min(model1$rerror)), 1)] #lowest error, not necessarely the most likely fit
error <- min(model1$rerror)
standard <- model1$standard

# import genes bed

bed <- read.table(bedfile,sep="\t",stringsAsFactors=FALSE)
colnames(bed) <- c("chr","start","end","gene")
bed[,1] <- gsub("chr","",bed[,1])
bed[,1] <- gsub("X","23",bed[,1])


## plot + tab

title <-  copyNumbersSegmented@phenoData@data[sample_nb,"name"]

template <- objectsampletotemplate(copyNumbersSegmented, sample_nb)

if (all(is.na(template$adjustedcopynumbers)) != TRUE){
    # adjusted copy numbers by standard which is the median bin segment value of a sample and times the ploidy
    template$adjustedcopynumbers <- ploidy + ((template$copynumbers - standard) * (cellularity * (ploidy - 2) + 2))/(cellularity * standard)
    template$adjustedsegments <- ploidy + ((template$segments - standard) * (cellularity * (ploidy - 2) + 2))/(cellularity * standard)
    template <- template[which(!template$chr%in%c("Y")),]

    segment <- Reduce(function(x,y) rbind(x,y),lapply(as.character(unique(template$chr)),function(chrom){
        tab <- template[which(as.character(template$chr)==chrom),]
        tabna <- tab[which(!is.na(tab$adjustedsegments)),]
        # tab[which(is.na(tab$adjustedsegments_max)),"adjustedsegments_max"] <- -1
        rl <- rle(as.vector(na.omit(tab$adjustedsegments)))
        u <- rl$lengths
        dd <- data.frame(chr=chrom,s=c(1,cumsum(u)[-1]-u[-1]+1),e=cumsum(u),seg=rl$values)
        dd <- as.data.frame(cbind(dd,t(apply(dd,1,function(x){
        start <- tabna[as.numeric(x["s"]),"start"]
        end <- tabna[as.numeric(x["e"]),"end"]
        return(c(start,end))}))))
        return(dd)}),accumulate=FALSE)
    colnames(segment) <- c("chr","startx","endx","adjustedsegments","start","end")

    write.table(segment[,c("chr","start","end","adjustedsegments")],paste(output_dir,"/",title,"_cnv_plot_segments.txt",sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

    bedGR <- GRanges(seqnames=bed[,1],IRanges(start=bed[,2],end=bed[,3]),name=bed[,"gene"])

    cnv_nano <- nano_cnvplot(template,bedGR,cnv_conf,cellularity,ploidy,title,"Nano")

    pdf(paste(output_dir,"/",title,"_cnv_plot.pdf",sep=""),width=20,height=6)
    print(cnv_nano)
    dev.off()

    ## amplification plot

    amp <- read.table(amp_file,stringsAsFactors=FALSE)

    if(!target_file == ""){
        targeted_genes <- read.table(target_file,stringsAsFactors=FALSE)
        amp[which(amp$V1%in%targeted_genes[,1]),"adaptive"] <- "Adaptive"
        amp[which(!(amp$V1%in%targeted_genes[,1])),"adaptive"] <- "Not Targeted"
    } else {
        amp[,"adaptive"] <- "Not Targeted"
    }

    p <- ggplot(amp,aes(x=V1,y=V5,fill=adaptive)) +
            geom_bar(stat="identity",position="dodge") +
            theme_minimal()+theme(legend.position="none",axis.title.y=element_text(size=6),axis.title.x=element_text(size=6),axis.text.x=element_text(angle=45,size=6,hjust=1),panel.grid.major=element_blank(),panel.background=element_rect(fill="grey98")) +
            # scale_x_discrete(expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0)) +
            labs(x="Genes", y="Average signal around gene") +
            ggtitle(title) + facet_wrap(~adaptive, scales="free_x") +
            scale_fill_manual(values=c("Adaptive"="tomato3","Not Targeted"="steelblue"))

    pdf(paste(output_dir,"/",title,"_amp_plot.pdf",sep=""),width=10,height=4)
    p
    dev.off()
}else{
    file.create(paste(output_dir,"/",title,"_cnv_plot_segments.txt",sep=""))
    file.create(paste(output_dir,"/",title,"_cnv_plot.pdf", sep=""))
    file.create(paste(output_dir,"/",title,"_amp_plot.pdf",sep=""))
}
