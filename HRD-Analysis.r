#devtools::install_github("buschlab/sequenza", build_vignettes = FALSE) sequenza installation!
#install_github('sztup/scarHRD',build_vignettes = TRUE) scar installation
#BiocManager::install("igordot/copynumber") copynumber installation

# Load required libraries
library(sequenza)
library(scarHRD)
library(copynumber)
library(knitr)
library(devtools)

# Specify the path to the .seqz file generated from sequenza-utils
data.file <- "~/Desktop/output-sequenza_files/9HRD/small.out.seqz"


#~/Downloads/scarHRD Analysis/out.seqz.txt"

# Extract Sequenza data from the .seqz file
test <- sequenza.extract(data.file, verbose = FALSE)


# Use test as input for scar_score
#scar_score(test, reference = "grch37", seqz = TRUE, chr.in.names = FALSE)

# Fit the Sequenza model
CP <- sequenza.fit(test)

# Generate results and save them to the specified output directory
sequenza.results(sequenza.extract = test,
                 cp.table = CP,
                 sample.id = "Test",
                 out.dir = "TEST")

# Results and description list
res_list <- c("alternative_fit.pdf", "alternative_solutions.txt",
              "chromosome_depths.pdf", "chromosome_view.pdf",
              "CN_bars.pdf", "confints_CP.txt",
              "CP_contours.pdf", "gc_plots.pdf",
              "genome_view.pdf", "model_fit.pdf",
              "mutations.txt", "segments.txt",
              "sequenza_cp_table.RData", "sequenza_extract.RData",
              "sequenza_log.txt")
res_list <- paste("Test", res_list, sep = "_")

description_list <- c(
  "Alternative solution fit to the segments. One solution per slide.",
  "List of all ploidy/cellularity alternative solutions.",
  "Coverage visualization in normal and tumor samples, before and after normalization.",
  "Chromosome view of depth ratio, B-allele frequency, and mutations. One chromosome per slide.",
  "Bar plot showing the percentage of genome in detected copy number states.",
  "Confidence interval table for the best solution from the model.",
  "Likelihood density for cellularity/ploidy solution with local maxima.",
  "GC correction visualization in normal and tumor samples.",
  "Genome-wide visualization of allele-specific and absolute copy number results.",
  "Model fit diagnostic plot.",
  "Table with mutation data and estimated number of mutated alleles (Mt).",
  "Table listing detected segments with copy number state estimates.",
  "RData file with maxima a posteriori computation.",
  "RData file of all sample information.",
  "Log with version and time information."
)
knitr::kable(data.frame(Files = res_list, Description = description_list))

# Load and display segmentation data, showing only segments with copy number <= 4
seg.tab <- read.table("TEST/Test_segments.txt", header = TRUE, sep = "\t")
alt_res <- read.table("TEST/Test_alternative_solutions.txt", header = TRUE, sep = "\t", fileEncoding = "UTF-8")
seg.tab <- seg.tab[seg.tab$CNt <= 4, ]
is.num <- sapply(seg.tab, is.numeric)
seg.tab[is.num] <- lapply(seg.tab[is.num], round, 3)
knitr::kable(head(seg.tab))

# Plot genome-wide views
sequenza:::genome.view(seg.tab)
sequenza:::genome.view(seg.tab, info.type = "CNt")

# Raw genome view
sequenza:::plotRawGenome(test, cellularity = alt_res$cellularity[1], ploidy = alt_res$ploidy[1])

# Cellularity and ploidy likelihood contours
cp.plot(CP)
cp.plot.contours(CP, add = TRUE, likThresh = c(0.999, 0.95), col = c("lightsalmon", "red"), pch = 20)

# Chromosome view with depth ratio, BAF, and segments
chromosome.view(mut.tab = test$mutations[[1]], baf.windows = test$BAF[[1]],
                ratio.windows = test$ratio[[1]], min.N.ratio = 1,
                segments = test$segments[[1]], main = test$chromosomes[1],
                cellularity = 0.89, ploidy = 1.9, avg.depth.ratio = 1)

# Calculate the HRD score using scarHRD with grch38 or grch37 as the reference genome
scar_score(data.file, reference = "grch37", seqz = TRUE, chr.in.names = TRUE) #chr.in.names = FALSE when chr are numerical(1234...22 )
