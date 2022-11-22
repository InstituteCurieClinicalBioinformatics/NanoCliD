# Create a cleaned up version of the annotSV table, keeping only columns deemed relevant for biological interpretation of somatic structural variants
library(readr)

arg <-
	commandArgs(trailingOnly = TRUE)

input.path <-
	arg[1]

output.path <-
	arg[2]

cols_to_keep <-
	c("AnnotSV_ID",
		"SV_chrom",
		"SV_start",
		"SV_end",
		"SV_length",
		"SV_type",
		"Nb_callers",
		"Sniffles",
		"cuteSV",
		"SVIM",
		"NanoVar",
		"Annotation_mode",
		"Gene_name",
		"Gene_count",
		"Tx",
		"Overlapped_tx_length",
		"Overlapped_CDS_length",
		"Overlapped_CDS_percent",
		"Frameshift",
		"Exon_count",
		"Location",
		"Location2",
		"Dist_nearest_SS",
		"Nearest_SS_type",
		"RE_gene",
#		"COSMIC_ID",
#		"COSMIC_MUT_TYP",
		"B_gain_source",
		"B_gain_coord",
		"B_loss_source",
		"B_loss_coord",
		"B_ins_source",
		"B_ins_coord",
		"B_inv_source",
		"B_inv_coord",
		"ENCODE_experiment",
		"Repeat_coord_left",
		"Repeat_type_left",
		"Repeat_coord_right",
		"Repeat_type_right",
		"SegDup_left",
		"SegDup_right",
		"ENCODE_blacklist_left",
		"ENCODE_blacklist_characteristics_left",
		"ENCODE_blacklist_right",
		"ENCODE_blacklist_characteristics_right",
		"OMIM_ID",
		"OMIM_phenotype",
		"OMIM_inheritance",
		"OMIM_morbid",
		"OMIM_morbid_candidate")

annotSV <-read_tsv(input.path)

annotSV$Nb_callers <-  substring(annotSV$INFO,6,6)
annotSV$Sniffles <-  substring(annotSV$INFO,17,17)
annotSV$cuteSV <- substring(annotSV$INFO,18,18)
annotSV$SVIM <-  substring(annotSV$INFO,19,19)
annotSV$NanoVar <-  substring(annotSV$INFO,20,20)

annotSV.cleaned <- annotSV[,cols_to_keep]

write.table(annotSV.cleaned,
	    file = output.path,
	    row.names = FALSE,
            sep = "\t",
            quote = FALSE)
