#load mafs using filepaths
myMaf <- read.csv("/Users/amarsandhu/Documents/R_Documents/kirc_maf_merged.txt", sep = "\t", stringsAsFactors = F)
library("maftools")
kirc_maf = read.maf(myMaf)

myMaf1 <- read.csv("/Users/amarsandhu/Documents/R_Documents/tRCC_TCGA.maf", sep = "\t", stringsAsFactors = F)
trcc_maf = read.maf(myMaf1)

#build trinucleotidematrix
kirc_tnm = trinucleotideMatrix(kirc_maf, ref_genome = "/Users/amarsandhu/Documents/R_Documents/b37.fasta")
trcc_tnm = trinucleotideMatrix(trcc_maf, ref_genome = "/Users/amarsandhu/Documents/R_Documents/b37.fasta")

#remove rows of all zeros
kirc_fm = kirc_tnm$nmf_matrix[-which(rowSums(kirc_tnm$nmf_matrix) == 0), ]

#build exposures matrices using quadratic programming method
library(SignatureEstimation)
library(deconstructSigs)
trcc_exposures = findSigExposures(t(trcc_tnm$nmf_matrix), t(signatures.cosmic), decomposition.method = decomposeQP)
kirc_exposures = findSigExposures(t(kirc_fm), t(signatures.cosmic), decomposition.method = decomposeQP)

#construct matrices with at least greatest 95% of signature contributions
temp = trcc_exposures
for(i in 1:12){
  for(j in 1:30){
    colmin = apply(temp$exposures, 2, which.min)
    if (colSums(temp$exposures, na.rm = T)[[i]] > 0.95 & colSums(temp$exposures, na.rm = T)[[i]] - min(temp$exposures[, i], na.rm = T) >= 0.95) {
      temp$exposures[colmin[[i]], i] = NA
    } else {
      next()
    }
  }
}
temp2 = kirc_exposures
for(i in 1:416){
  for(j in 1:30){
    colmin2 = apply(temp2$exposures, 2, which.min)
    if (colSums(temp2$exposures, na.rm = T)[[i]] > 0.95 & colSums(temp2$exposures, na.rm = T)[[i]] - min(temp2$exposures[, i], na.rm = T) >= 0.95) {
      temp2$exposures[colmin2[[i]], i] = NA
    } else {
      next()
    }
  }
}

#test statistical significance
cosineSim <- function(x, y){
  cos_mat <- (t(x)%*%y)/(sqrt((rowSums(x^2) %>% as.numeric() %>% t()) %*% (rowSums(y^2) %>% as.numeric())))[1,1]
  return(rowMeans(cos_mat))
}
mean_cos = cosineSim(trcc_exposures$exposures, kirc_exposures$exposures)

tmat = numeric()
for (i in 1:30){
  tmat = c(tmat, t.test(trcc_exposures$exposures[i, ], kirc_exposures$exposures[i, ])$p.value)
}

umat = numeric()
for (i in 1:30){
  umat = c(umat, wilcox.test(trcc_exposures$exposures[i, ], kirc_exposures$exposures[i, ])$p.value)
}
# trcc3 = melt(trcc_exposures$exposures, id.vars = 30)
# kirc3 = melt(kirc_exposures$exposures, id.vars = 30)
# utest = wilcox.test(x = trcc3$value, y = kirc3$value)

#plot data - do this in new window or click zoom - the default r window might remain blank
library(reshape2)
trcc2 = melt(temp, id.vars = 30)
kirc2 = melt(temp2, id.vars = 30)

library(ggplot2)
g1 = ggplot(trcc2[which(!is.na(trcc2$value) & trcc2$L1 == "exposures"), ], aes(x = Var2, y = value, fill = Var1)) +
      geom_bar(stat = "identity") +
      xlab("Mean Sample Cosine Similarity to KIRC") +
      ylab("Contributions") +
      ggtitle(paste("tRCC Signature Contributions(≥95%) per Sample, Mean Cosine Similarity for All Samples = ", round(mean(mean_cos), 3))) +
      scale_x_discrete(labels = as.character(round(mean_cos, 3))) +
      theme_bw()
g2 = ggplot(kirc2[which(!is.na(kirc2$value) & kirc2$L1 == "exposures"), ], aes(x = Var2, y = value, fill = Var1)) +
      geom_bar(stat = "identity") +
      xlab("Sample") +
      ylab("Contributions") +
      ggtitle("KIRC Signature Contributions(≥95%) per Sample") +
      scale_x_discrete(labels = NULL) +
      theme_bw()

spltrcc = as.data.frame(as.character(trcc2[which(!is.na(trcc2$value) & trcc2$L1 == "exposures"), ]$Var1) %>% strsplit(split = '.', fixed = T))
spltrcc = as.character.numeric_version(spltrcc[2, ])
uqspltrcc = as.numeric(unique(spltrcc))

g3 = ggplot(trcc2[which(!is.na(trcc2$value) & trcc2$L1 == "exposures"), ], aes(x = Var1, y = value)) +
      geom_boxplot() +
      xlab("Signature") +
      ylab("Contribution") +
      scale_x_discrete(labels = sort(uqspltrcc)) +
      scale_y_continuous(limits = c(0.0, 0.7)) +
      ggtitle("tRCC Signature Contributions")
g4 = ggplot(kirc2[which(!is.na(kirc2$value) & kirc2$L1 == "exposures"), ], aes(x = Var1, y = value)) +
      geom_boxplot() +
      xlab("Signature") +
      ylab("Contribution") +
      scale_x_discrete(labels = as.character(1:30)) +
      scale_y_continuous(limits = c(0.0, 0.7))
      ggtitle("KIRC Signature Contributions")

library(gridExtra)
grid.arrange(g1, g3, g2, g4, nrow = 2)

#plot t- and u-tests
tmat2 = data.frame(tmat, signature = 1:30)
ggplot(tmat2, aes(x = signature, y = 1 - tmat, fill = "magenta")) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks = 1:30) +
  theme_bw() +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "blue") +
  guides(fill = F)
umat2 = data.frame(umat, signature = 1:30)
ggplot(umat2, aes(x = signature, y = 1 - umat, fill = "magenta")) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks = 1:30) +
  theme_bw() +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "blue") +
  guides(fill = F)
