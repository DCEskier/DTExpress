args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
	stop("At least one argument must be supplied (input file).n", call. = FALSE)
} else if (length(args) == 1) {
	# default output file
	args[2] = "output.bed"
}

options(stringsAsFactors = F)

#read.csv("GCF_000001735.4_TAIR10.1_feature_table.txt", sep = "\t") -> features

read.csv(args[1], sep = "\t") -> features

features.subset <- subset(features, X..feature == "gene")

which(features.subset$symbol == "") -> symbolMissing

features.subset$symbol[symbolMissing] <-
	features.subset$locus_tag[symbolMissing]

features.bed <- features.subset[, c(7, 8, 9, 15, 16, 10)]

write.table(
	features.bed,
	args[2],
	quote = F,
	col.names = F,
	row.names = F,
	sep = "\t"
)

