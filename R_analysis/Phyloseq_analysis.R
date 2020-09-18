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

# Import phylogenetic tree
tree <- ape::read.tree("analysis_files/")

# Alpha-div analysis
library(phyloseq)
phy.obj <- phyloseq(otu_table(data, taxa_are_rows = TRUE), sample_data(meta))

## Rarefy
set.seed(123)
phy.obj.rare = rarefy_even_depth(phy.obj)

### Alpha div analysis
alpha_div = plot_richness(phy.obj.rare, color = "Genotype", x = "Genotype",
                          measures = c("Observed", "Shannon", "InvSimpson","Fisher"), 
                          title = "Alpha Diversity of ASVs")
jpeg('Images/Alpha_div_12wks.jpg', units = "in", width = 5, height = 5,
     res = 300)
alpha_div + geom_point(size = 2, alpha = 0.8) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()

### Data frame for kruskal wallis test
alphadt <- alpha_div$data
obs = alphadt[which(alphadt$variable == "Observed"),]
invS = alphadt[which(alphadt$variable == "InvSimpson"),]
shan = alphadt[which(alphadt$variable == "Shannon"),]
fisher = alphadt[which(alphadt$variable == "Fisher"),]

kruskal.test(value ~ Genotype, data = obs) #p = 0.351
kruskal.test(value ~ Genotype, data = invS) #p = 0.008
kruskal.test(value ~ Genotype, data = shan) #p = 0.05
kruskal.test(value ~ Genotype, data = fisher) #p = 0.351

#### Testing for male diff in Shannon index
shan.m = shan[which(shan$Sex == "M"),]
kruskal.test(value ~ Genotype, data = shan.m) #p = 0.021

### Barplots 
obs_plot <- plotIndiv_bar(obs, x = Genotype, y = value, title = "Number of species observed") + 
  facet_grid(~Sex) +
  scale_fill_manual(values = c('#388ECC','#F68B33'))
shan_plot <- plotIndiv_bar(shan, x = Genotype, y = value, title = "Species evenness") + 
  facet_grid(~Sex) +
  scale_fill_manual(values = c('#388ECC','#F68B33'))
plotIndiv_bar(fisher, x = Genotype, y = value, title = "Species fisher") + 
  facet_grid(~Sex) +
  scale_fill_manual(values = c('#388ECC','#F68B33'))
plotIndiv_bar(invS, x = Genotype, y = value, title = "Species evenness") + 
  facet_grid(~Sex) +
  scale_fill_manual(values = c('#388ECC','#F68B33'))

jpeg(filename = "Images/Observed_Shannon_barplot.jpg", width = 15, height = 10, units = "cm", res = 300)
gridExtra::grid.arrange(obs_plot, shan_plot, nrow = 1)
dev.off()

# Beta div analysis
## Filtering low counts
data.filter <- low.count.removal(data, percent = 0.01)
data.filter <- data.filter$data.filter
dim(data.filter)

## Transform to relative abundance
data.TSS <- apply(data.filter, 2, TSS.divide)

## To phyloseq
phy.obj.ra <- phyloseq(otu_table(data.TSS, taxa_are_rows = TRUE), sample_data(meta))

## Jaccard distance
distjac = distance(phy.obj.ra, method = "jaccard")
vegan::adonis(distjac ~ Genotype, data.frame(sample_data(phy.obj.ra)))
vegan::adonis(distjac ~ Genotype, meta)

## Bray Curtis dissimilarity
distBC = distance(phy.obj.ra, method = "bray")
vegan::adonis(distBC ~ Genotype, meta)
ordBC = ordinate(phy.obj.ra, method = "PCoA", distance = distBC)
pcoa_bc = plot_ordination(phy.obj.ra, ordBC, color = "Sex", shape = "Genotype") + 
  geom_point(size = 2, alpha = 0.8) + ggtitle("Principal Components Analysis plot") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = c('#388ECC','#F68B33')) + 
  labs(y="PCoA 2 [11.3%]", x = "PCoA 1 [50.1%]")
jpeg(filename = "Images/PCoA_Bray_curtis.jpg", width = 12, height = 11, units = "cm", res = 300)
pcoa_bc
dev.off()


