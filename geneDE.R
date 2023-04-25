args = commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
	stop("At least one argument must be supplied (input file).n", call. = FALSE)
} else if (length(args) == 1) {
	# default output file
	args[4] = "output.txt"
}

options(stringsAsFactors = F)
library(edgeR)

read.csv(args[1], row.names = 1) -> counts

apply(counts, 1, sum) -> sums

counts[which(sums != 0), ] -> counts

pData <- read.delim(args[2], header = F, row.names = NULL)

meta <- data.frame(
	row.names = pData[, 1],
	condition = pData[, 2],
	libsize = pData[, 3]
)

libsize <- meta$libsize
condition <- factor(meta$condition)
design <- model.matrix( ~ 0 + condition)
colnames(design) <- levels(as.factor(meta$condition))

y <- DGEList(counts = counts, lib.size = libsize)

# Normalize the data
y <- calcNormFactors(y)
y$samples
plotMDS(y)

# Estimate the variance
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
plotBCV(y)

# Build an object to contain the normalized read abundance
logcpm <- cpm(y, log = TRUE, lib.size = libsize)
logcpm <- as.data.frame(logcpm)
colnames(logcpm) <- factor(meta$condition)

# Conduct fitting of the GLM
yfit <- glmFit(y, design)

# Initialize result matrices to contain the results of the GLM
results <- matrix(nrow = dim(counts)[1], ncol = 0)
logfc <- matrix(nrow = dim(counts)[1], ncol = 0)

read.delim(
	args[3],
	header = F,
	row.names = NULL
) -> contrasts.df

contrasts.vector <- as.vector(contrasts.df[, 2])

contrasts.names <- as.vector(contrasts.df[, 1])

my.contrasts <-
	makeContrasts(contrasts = contrasts.vector, levels = design)

attr(my.contrasts, "dimnames")$Contrasts <- contrasts.names

# Define the contrasts used in the comparisons
allcontrasts = contrasts.names

# Conduct a for loop that will do the fitting of the GLM for each comparison
# Put the results into the results objects
for (current_contrast in allcontrasts) {
	lrt <- glmLRT(yfit, contrast = my.contrasts[, current_contrast])
	plotSmear(lrt, de.tags = rownames(y))
	title(current_contrast)
	res <- topTags(lrt, n = dim(c)[1], sort.by = "none")$table
	colnames(res) <- paste(colnames(res), current_contrast, sep = ".")
	results <- cbind(results, res[, c(1, 5)])
	logfc <- cbind(logfc, res[c(1)])
}

# Sort the results table by the logFC
#results <- results[with(results, order(-abs(as.name(paste("logFC.", allcontrasts[1], sep = ""))))),]

# Save the results
write.table(results, args[4], quote = FALSE, sep = "\t")