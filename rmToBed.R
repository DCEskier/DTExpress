args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
	stop("At least one argument must be supplied (input file).n", call. = FALSE)
} else if (length(args) == 1) {
	# default output file
	args[2] = "output.bed"
}

library(biomartr)
library(dplyr)
library(stringr)

hg<-read_rm(args[1])

#hg<-read_rm("~/Documents/MBG/Arabidopsis/GCF_000001735.4_TAIR10.1_rm.out")
last<-as.data.frame(str_split_fixed(hg$matching_class, "/", 2))
new_hg<-data.frame("chr"=hg$qry_id, "start"=hg$qry_start, "end"=hg$qry_end, "strand"=hg$matching_repeat, "repeat"=hg$repeat_id, "repeat_class"=last$V1, "repeat_family"=last$V2)
new_hg$strand<-as.character(new_hg$strand)
new_hg$strand <- replace(new_hg$strand, new_hg$strand=="C", "-")
new_hg$repeat_family[which(new_hg$repeat_family == "")] <- NA
new_hg$score <- rep(1, nrow(new_hg))
new_hg <- new_hg[, c("chr", "start", "end", "repeat.", "score", "strand", "repeat_class", "repeat_family")]

write.table(new_hg, args[2], quote = F, row.names = F, col.names = F, sep = "\t")