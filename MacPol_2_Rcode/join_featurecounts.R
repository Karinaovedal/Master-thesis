setwd("/Volumes/home/Karina/masterthesis/RNAseq_PGE2/featureCounts")

file_list <- list.files(pattern = "*.counts.txt$")

### Read in file names
datalist <- lapply(file_list, function(x) {
    dat <- read.table(file = x, header = TRUE, sep = "\t")
    names(dat)[2] <- x
    return(dat)
})

# Brew version of r studio
library(plyr, lib.loc = "/opt/homebrew/lib/R/4.3/site-library")
library(dplyr, lib.loc = "/opt/homebrew/lib/R/4.3/site-library")
# Web downloaded version of r studio
#library(plyr, lib.loc = "/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library")
#library(dplyr, lib.loc = "/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library")
joined <- join_all(dfs = datalist, by = "Geneid", type = "full")

tail(joined, 5)
head(joined, 5)
dim(joined)
colnames(joined)
joined2 <- dplyr::select(joined, contains("bam"))
head(joined2)
dim(joined2)
rownames(joined2) <- joined$Geneid
rownames(joined2) <- gsub(" ", "_", rownames(joined2))
for (col in 1:ncol(joined2)) {
    colnames(joined2)[col] <- sub("..", "", colnames(joined2)[col])
} # remove X from column names
for (col in 1:ncol(joined2)) {
    colnames(joined2)[col] <- sub(".bam", "", colnames(joined2)[col])
}

write.table(joined2, file = "Count_DataFrame.csv", quote = F, sep = ",", row.names = TRUE, col.names = NA)
