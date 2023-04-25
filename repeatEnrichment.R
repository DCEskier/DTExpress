args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
	stop("At least one argument must be supplied (input file).n", call. = FALSE)
} else if (length(args) == 1) {
	# default output file
	args[5] = "output.txt"
}

library(GenomicRanges)
library(dplyr)

options(stringsAsFactors = F)

makeGRangeObj <- function(df) {
	gr <- with(df, GRanges(chr, IRanges(start, end), strand = strand))
	values(gr) <- df[, 4:5]
	return(gr)
}


toFindOverlaps <- function(gr_repeats, gr_genome) {
	# Find overlaps
	names(values(gr_genome)) <- c("gName", "gScore")
	names(values(gr_repeats)) <- c("rName", "rScore")
	m <- findOverlaps(gr_genome, gr_repeats, ignore.strand = T)
	gr_genome.matched <- gr_genome[queryHits(m)]
	
	# Add the metadata from gr2
	mcols(gr_genome.matched) <-
		cbind.data.frame(mcols(gr_genome.matched),
										 mcols(gr_repeats[subjectHits(m)]))
	
	return(gr_genome.matched)
}

getUpstream <- function(df, length) {
	dftemp <- df
	dfnew <- df
	
	indexP <- dftemp$strand == "+"
	indexN <- dftemp$strand == "-"
	
	dfnew$start[indexP] <- c(dftemp$start[indexP] - length)
	dfnew$end[indexN] <- c(dftemp$end[indexN] + length)
	
	
	return(dfnew)
}


allGenes <- read.delim(args[1], header = F)[, 1:6]
nrow(allGenes) -> sumG
names(allGenes) <-
	c("chr", "start", "end", "name", "score", "strand")
allGenes <- getUpstream(allGenes, as.numeric(args[4]))
allGenesGenomicRanges <- makeGRangeObj(allGenes)


allRepeats <-
	read.delim(args[2], header = F)
names(allRepeats) <-
	c("chr", "start", "end", "name", "score", "strand", "class", "family")
allRepeatsFiltered <- subset(allRepeats, class != "Low_complexity" & class != "Simple_repeat")[, 1:6]
names(allRepeatsFiltered) <-
	c("chr", "start", "end", "name", "score", "strand")
nrow(allRepeatsFiltered) -> sumR
allRepeatsGenomicRanges <- makeGRangeObj(allRepeatsFiltered)

toFindOverlaps(allRepeatsGenomicRanges, allGenesGenomicRanges) -> overlaps

genesDE <- read.delim(args[3], header = F)[, 1:6]
names(genesDE) <-
	c("chr", "start", "end", "name", "score", "strand")
genesDE <- getUpstream(genesDE, as.numeric(args[4]))
genesDEGenomicRanges <- makeGRangeObj(genesDE)

toFindOverlaps(allRepeatsGenomicRanges, genesDEGenomicRanges) -> overlapsDE

metadataOverlapsDE <- elementMetadata(overlapsDE)

list <-
	data.frame(table(metadataOverlapsDE$rName))  #to find count of the frame fields
header <- c('Repeats', 'Frequency')
names(list) <- header
orderList <- list %>% arrange(list$Frequency) #to order list asc
repeatFrequencyDE <- orderList[orderList$Frequency != 0,]

metadataOverlapsAll <- elementMetadata(overlaps)

list <-
	data.frame(table(metadataOverlapsAll$rName))  #to find count of the frame fields
header <- c('Repeats', 'Frequency')
names(list) <- header
orderList <- list %>% arrange(list$Frequency) #to order list asc
repeatFrequencyAll <- orderList[orderList$Frequency != 0,]

calculateESx_upGene <- function(calc_r, calc_g) {
	total <- right_join(calc_r, calc_g, by = "Repeats")
	names(total) <- c("Repeats", "r", "g")
	total$r[is.na(total$r)] <- 0
	sumR <- sum(total$r)
	sumG <- sum(total$g)
	total$R <- c(sumR)
	total$G <- c(sumG)
	total$ESx <- c((total$r / total$R) / (total$g / total$G))
	p <-
		apply(as.matrix(total[, 2:5]), 1, function(x)
			# to calculate pvalue with using fisher exact test
			fisher.test(matrix(round(x), ncol = 2), workspace = 1e9)$p.value)
	total$pValue <- signif(p, digits = 6)
	pAdjust <-
		p.adjust(p, method = "BH", n = length(p))     # to calculate p adjust value with using binom test
	total$pAdjust <- signif(pAdjust, digits = 6)
	return(total)
}

calculateESx_upGene(repeatFrequencyDE, repeatFrequencyAll) -> calculatedESx

calculatedESx <- subset(calculatedESx, r != 0)

write.table(
	calculatedESx,
	args[5],
	sep = "\t",
	quote = F,
	row.names = F
)
