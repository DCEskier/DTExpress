args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	stop(
		"At least one repeat count matrix, one phenotype data table, and one contrast table must be supplied (input file).n",
		call. = FALSE
	)
} else if (length(args) < 4) {
	# default output file
	args[4] = "results.csv"
}

options(stringsAsFactors = F)
library(edgeR)

read.delim(args[1]) -> counts

read.delim(args[2], row.names = NULL, header = F) -> pData

nrow(pData) -> sampleNum

names(counts)[1:8] <- c("Chromosome",
												"Start",
												"End",
												"Name",
												"Score",
												"Strand",
												"Type",
												"Family")

names(counts)[9:(8 + sampleNum)] <- pData[, 1]

b <-
	aggregate(
		list(counts[9:(8 + sampleNum)]),
		by = list(
			geneName = counts$Name,
			repeatClass = counts$Type ,
			repeatFamily = counts$Family
		),
		FUN = sum
	)

counts <- b

counts$Sum <- apply(counts, 1, function(x)
	sum(as.numeric(x[4:(3 + sampleNum)])))

counts <-
	subset(counts,
				 repeatClass != "Simple_repeat" &
				 	repeatFamily != "Low_complexity" & Sum != 0)

row.names(counts) <- make.names(counts$geneName, unique = T)

counts.annotation <- counts[, 1:3]

counts <- counts[, 4:(3 + sampleNum)]

libsize <- pData[, 3]
condition <- factor(pData[, 2])
design <- model.matrix( ~ 0 + condition)
colnames(design) <- levels(as.factor(pData[, 2]))

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
colnames(logcpm) <- factor(pData[, 2])

# Conduct fitting of the GLM
yfit <- glmFit(y, design)

# Initialize result matrices to contain the results of the GLM
results <- matrix(nrow = dim(counts)[1], ncol = 0)
logfc <- matrix(nrow = dim(counts)[1], ncol = 0)

read.delim(
	args[3],
	header = F,
	sep = "\t",
	row.names = NULL
) -> contrasts.df

contrasts.vector <- as.vector(contrasts.df[, 2])

contrasts.names <- as.vector(contrasts.df[, 1])

my.contrasts <-
	makeContrasts(contrasts = contrasts.vector, levels = design)

attr(my.contrasts, "dimnames")$Contrasts <- contrasts.names

# Define the contrasts used in the comparisons
allcontrasts = c(contrasts.names)

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

results$class <- counts.annotation$repeatClass
results$type <- counts.annotation$repeatFamily

# Sort the results table by the logFC
results <- results[with(results, order(-abs(logFC.salt))),]

# Save the results
write.table(results, 'results.txt', quote = FALSE, sep = "\t")

# Plot Fold Changes for repeat classes and types
for (current_contrast in allcontrasts) {
	logFC <- results[, paste0("logFC.", current_contrast)]
	# Plot the repeat classes
	classes <- with(results, reorder(class, -logFC, median))
	par(mar = c(6, 10, 4, 1))
	boxplot(
		logFC ~ classes,
		data = results,
		outline = FALSE,
		horizontal = TRUE,
		las = 2,
		xlab = "log(Fold Change)",
		main = current_contrast
	)
	abline(v = 0)
	# Plot the repeat types
	types <- with(results, reorder(type, -logFC, median))
	boxplot(
		logFC ~ types,
		data = results,
		outline = FALSE,
		horizontal = TRUE,
		las = 2,
		xlab = "log(Fold Change)",
		main = current_contrast
	)
	abline(v = 0)
}

write.csv(results, args[4])