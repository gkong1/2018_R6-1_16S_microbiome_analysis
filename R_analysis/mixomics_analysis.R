setwd("~/phd/Kong et al. 2018. Gut dysbiosis in transgenic mouse model of HD/")
source("~/phd/scripts/common_functions.R")
# Import data
data <- read.delim("analysis_files/count_table_raw.txt", header = T, row.names = 1)
data <- data[,order(colnames(data))]

# Import meta
meta <- read.delim("analysis_files/Sample_metadata.tsv", header = T, row.names = 1)
meta <- meta[order(rownames(meta)),]
meta$Sex = factor(meta$Sex, levels = c("M","F"))
meta$Genotype = factor(meta$Genotype, levels = c("WT","HD"))

# Check whether sample names match
colnames(data)
rownames(meta)
colnames(data) = rownames(meta)

data.TSS <- data_cleanup(data, samples_as_rows = FALSE, filter_percent = 0.01, 
                         offset = TRUE, offset_number = 0.01)
# Mixomics analysis
library(mixOmics)
pca.res <- pca(data.TSS, ncomp = 10)

# Plot PCA
plot(pca.res)
plotIndiv(pca.res, group = meta$Genotype, pch = meta$Sex, legend = T)

# SPLS-DA analysis
set.seed(33)  # for reproducible results for this code
data.tune.splsda = tune.splsda(data.TSS, 
                               Y = meta$Genotype, 
                               ncomp = 3, 
                               multilevel = NULL, 
                               logratio = 'none',
                               test.keepX = c(seq(5,150, 5)), 
                               validation = c('loo'), 
                               dist = 'max.dist', # prediction distance can be chosen according to tune.plsda results
                               progressBar = FALSE)
plot(data.tune.splsda)

#select.keepX = data.tune.splsda.1$choice.keepX[1:2]
select.keepX = c(20, 20) # originally c(50,35)

data.splsda = splsda(X = data.TSS,  Y = meta$Genotype, 
                     ncomp = 2, keepX = select.keepX, 
                     logratio= "CLR")

jpeg('Images/sPLSDA_tp12.jpg', units = "in", width = 6, height = 5,
     res = 300)
plotIndiv(data.splsda, 
          ind.names = F, 
          #col.per.group = color.mixo(1:2), 
          comp = c(1,2),
          pch = meta$Genotype, 
          ellipse = TRUE,
          legend = TRUE,
          title = 'sPLS-DA (tp12): Combined', 
          X.label = "Component 1", 
          Y.label = "Component 2")
dev.off()

set.seed(34)  # for reproducible results for this code
data.perf.splsda = perf(data.splsda, validation = 'loo', 
                        progressBar = FALSE, dist = 'max.dist')

data.perf.splsda$error.rate


plot(data.perf.splsda)

head(selectVar(data.splsda, comp = 1)$value) 
varcontrib = selectVar(data.splsda, comp = 1)$value
var.comp1 = selectVar(data.splsda, comp = 1)$name

# stability of OTUs selected on comp 1
varcontrib$stability = data.perf.splsda.1$features$stable[[1]][var.comp1.metab]

# Loading plot
jpeg('Images/Loading plot_12wks.jpg', units = "in", width = 6, height = 5,
     res = 300)
plotLoadings(data.splsda, comp = 1, method = 'mean', contrib = 'max',
             size.title = 1, ndisplay = 20, size.name = 0.65, size.legend = 0.6, title = "Contribution on comp 1: Combined M&F")
dev.off()

#-----------------------------------------------------------------
# Male sPLS-DA
index <- which(meta$Sex == "M")
meta.m <- meta[index,]
data.TSS.m <- data.TSS[index,]
set.seed(33)  # for reproducible results for this code
data.tune.splsda.m = tune.splsda(data.TSS.m, 
                               Y = meta.m$Genotype, 
                               ncomp = 3, 
                               multilevel = NULL, 
                               logratio = 'none',
                               test.keepX = c(seq(5,150, 5)), 
                               validation = c('loo'), 
                               dist = 'max.dist', # prediction distance can be chosen according to tune.plsda results
                               progressBar = FALSE)
plot(data.tune.splsda.m)

#select.keepX = data.tune.splsda.1$choice.keepX[1:2]
select.keepX = c(10, 20) # originally c(50,35)

data.splsda.m = splsda(X = data.TSS.m,  Y = meta.m$Genotype, 
                     ncomp = 2, keepX = select.keepX, 
                     logratio= "CLR")

jpeg('Images/sPLSDA_tp12_male.jpg', units = "in", width = 6, height = 5,
     res = 300)
plotIndiv(data.splsda.m, 
          ind.names = F, 
          #col.per.group = color.mixo(1:2), 
          comp = c(1,2),
          pch = meta.m$Genotype, 
          ellipse = TRUE,
          legend = TRUE,
          title = 'Genes (tp12) sPLS-DA', 
          X.label = "Component 1", 
          Y.label = "Component 2")
dev.off()

set.seed(34)  # for reproducible results for this code
data.perf.splsda.m = perf(data.splsda.m, validation = 'loo', 
                        progressBar = FALSE, dist = 'max.dist')

data.perf.splsda.m$error.rate


plot(data.perf.splsda.m)

head(selectVar(data.splsda.m, comp = 1)$value) 
varcontrib.m = selectVar(data.splsda.m, comp = 1)$value
var.comp1.m = selectVar(data.splsda.m, comp = 1)$name

# stability of OTUs selected on comp 1
varcontrib.m$stability = data.perf.splsda.m$features$stable[[1]][var.comp1.metab]

# Loading plot
jpeg('Images/Loading plot_12wks_male.jpg', units = "in", width = 6, height = 5,
     res = 300)
plotLoadings(data.splsda.m, comp = 1, method = 'mean', contrib = 'max',
             size.title = 1, ndisplay = 5, size.name = 0.65, size.legend = 0.6, 
             title = "Contribution comp 1 - Male")
dev.off()

taxa <- read.delim(file = "analysis_files/tax_id.tsv", header = T, row.names = 1)
taxa.m.sig <- taxa[which(rownames(taxa) %in% var.comp1.m[1:5]),]
#-----------------------------------------------------
# Female sPLS-DA
index <- which(meta$Sex == "F")
meta.f <- meta[index,]
data.TSS.f <- data.TSS[index,]

set.seed(33)  # for reproducible results for this code
data.tune.splsda.f = tune.splsda(data.TSS.f, 
                               Y = meta.f$Genotype, 
                               ncomp = 3, 
                               multilevel = NULL, 
                               logratio = 'none',
                               test.keepX = c(seq(5,150, 5)), 
                               validation = c('loo'), 
                               dist = 'max.dist', # prediction distance can be chosen according to tune.plsda results
                               progressBar = FALSE)
plot(data.tune.splsda.f)

#select.keepX = data.tune.splsda.1$choice.keepX[1:2]
select.keepX = c(10, 20) # originally c(50,35)

data.splsda.f = splsda(X = data.TSS.f,  Y = meta.f$Genotype, 
                     ncomp = 2, keepX = select.keepX, 
                     logratio= "CLR")

jpeg('Images/sPLSDA_tp12_female.jpg', units = "in", width = 6, height = 5,
     res = 300)
plotIndiv(data.splsda.f, 
          group = meta.f$Genotype,
          ind.names = F, 
          col.per.group = color.mixo(1:2), 
          comp = c(1,2),
          pch = meta.f$Genotype, 
          ellipse = TRUE,
          legend = TRUE,
          title = 'Genes (tp12) sPLS-DA: Female data', 
          X.label = "Component 1", 
          Y.label = "Component 2")
dev.off()

set.seed(34)  # for reproducible results for this code
data.perf.splsda.f = perf(data.splsda.f, validation = 'loo', 
                        progressBar = FALSE, dist = 'max.dist')

data.perf.splsda.f$error.rate


plot(data.perf.splsda.f)

head(selectVar(data.splsda.f, comp = 1)$value) 
varcontrib.f = selectVar(data.splsda.f, comp = 1)$value
var.comp1.f = selectVar(data.splsda.f, comp = 1)$name

# stability of OTUs selected on comp 1
varcontrib.f$stability = data.perf.splsda.f$features$stable[[1]][var.comp1.f]

jpeg('Images/Loading plot_12wks_female.jpg', units = "in", width = 6, height = 5,
     res = 300)
plotLoadings(data.splsda.f, comp = 1, method = 'mean', contrib = 'max',
             size.title = 1, ndisplay = 10, size.name = 0.65, size.legend = 0.6, title = "Contribution comp 1 - Female")
dev.off()

#-----------------------------------------------------
# Loading plot






