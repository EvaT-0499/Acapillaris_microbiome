# Script for analysis after phyloseq object
# Eva Tanneau

# load packages
library("DECIPHER") 
library("ape")
library("DESeq2")
library("ggplot2")
library("phyloseq")
library("plotly")
library("vegan")
library("philr")
library("tidyverse")
library("adespatial")
library("devtools")
library("qiime2R")
library("MicrobeR")
library("microbiome")
library("microbiomeSeq")
library("pander")
library("ranacapa")
library("grid")
library("gridExtra")
library("knitr")
library("png")
library("ggpubr")
library("RColorBrewer")
library("remotes")
library("microbiomeutilities")
library("DAtest")
library("MicEco")
library("ggrepel")
library("SpiecEasi")
library("NetCoMi")
library("dplyr")
library("gridExtra")
library("usethis")
library("devtools")
library("cowplot")
library("wesanderson")


# load the pruned ps object - according the right year
ps.ITS <- readRDS("/ps_ITS_2022.rds")
ps.16S <- readRDS("/ps_16S_2022_.rds")


# load rarefied ps object - according the right year
ps.rar.ITS <- readRDS("/ps_ITS_2022_rarefied.rds")
ps.rar.16S <- readRDS("/ps_16S_2022_rarefied.rds")

# Set text size for figures
My_Theme = theme(axis.title.x = element_text(size = 20),
                 axis.text.x = element_text(angle = 0, size = 20, hjust = 0.5),
                 axis.title.y = element_text(size = 20),
                 text = element_text(size = 20))


# Subset leaves or roots
ps.16S.leaf <- subset_samples(ps.16S, Material == "Leaf")
ps.ITS.leaf <- subset_samples(ps.ITS, Material == "Leaf")

ps.rar.16S.leaf <- subset_samples(ps.rar.16S, Material == "Leaf")
ps.rar.ITS.leaf <- subset_samples(ps.rar.ITS, Material == "Leaf")

ps.rar.16S.root <- subset_samples(ps.rar.16S, Material == "Root")
ps.rar.ITS.root <- subset_samples(ps.rar.ITS, Material == "Root")

# Remove soil samples
ps.rar.16S <- subset_samples(ps.rar.16S, Material=="Leaf"|Material=="Root")
ps.rar.ITS <- subset_samples(ps.rar.ITS, Material=="Leaf"|Material=="Root")


################################################################

# Blocks (verify that the blocks are comparable)
# Supp Figure 2

# Bacteria - Full dataset
plot_diversity_stats(ps.16S, group = "Block", 
                     index = "diversity_shannon", 
                     group.order = c("61","62"),                      
                     #group.colors = mycols,
                     label.format="p.format",
                     stats = TRUE,
) + ylab("Shannon Diversity") + xlab("")+ scale_fill_brewer( palette = "Dark2") + scale_color_brewer(palette = "Dark2") + My_Theme

# Fungi - Full dataset
plot_diversity_stats(ps.ITS, group = "Block", 
                     index = "diversity_shannon", 
                     group.order = c("61","62"),                      
                     #group.colors = mycols,
                     label.format="p.format",
                     stats = TRUE,
) + ylab("Shannon Diversity") + xlab("")+ scale_fill_brewer( palette = "Dark2") + scale_color_brewer(palette = "Dark2") + My_Theme


################################################################

# Plant compartment comparison
# Figure 2

# Bacteria - Alpha diversity
plot_diversity_stats(ps.rar.16S, group = "Material", 
                     index = "diversity_shannon", 
                     group.order = c("Leaf","Root", "Soil"),                      
                     #group.colors = mycols,
                     label.format="p.format",
                     stats = TRUE,
) + ylab("Shannon Diversity") + xlab("")+ scale_fill_brewer( palette = "Accent") + scale_color_brewer(palette = "Accent") + My_Theme


# Fungi - Alpha diversity
plot_diversity_stats(ps.rar.ITS, group = "Material", 
                     index = "diversity_shannon", 
                     group.order = c("Leaf","Root", "Soil"),                      
                     #group.colors = mycols,
                     label.format="p.format",
                     stats = TRUE,
) + ylab("Shannon Diversity") + xlab("")+ scale_fill_brewer( palette = "Accent") + scale_color_brewer(palette = "Accent") + My_Theme


# Bacteria - Beta diversity
physeq.ord <- ordinate(ps.rar.16S, "PCoA", "bray")
b.div.bray <- plot_ordination(ps.rar.16S, physeq.ord, type= "samples", color= "Material") + geom_point(size=4)
b.div.bray <- b.div.bray + theme_classic() + scale_color_brewer("Material", palette = "Accent", name="Plant compartment" ) + My_Theme
print(b.div.bray)

# Fungi - Beta diversity
physeq.ord <- ordinate(ps.rar.ITS, "PCoA", "bray")
b.div.bray <- plot_ordination(ps.rar.ITS, physeq.ord, type= "samples", color= "Material") + geom_point(size=4)
b.div.bray <- b.div.bray + theme_classic() + scale_color_brewer("Material", palette = "Accent", name="Plant compartment") + My_Theme
print(b.div.bray)

#Permanova Bacteria
ps.prop.16S <- transform_sample_counts(ps.rar.16S, function(otu) otu/sum(otu))
bray.dist <- phyloseq::distance(ps.prop.16S, method = "bray")
metadata <- as(sample_data(ps.prop.16S), "data.frame")
adonis2(bray.dist ~ Material*Treatment,
        data = metadata)

#Permanova Fungi
ps.prop.ITS <- transform_sample_counts(ps.rar.ITS, function(otu) otu/sum(otu))
bray.dist <- phyloseq::distance(ps.prop.ITS, method = "bray")
metadata <- as(sample_data(ps.prop.ITS), "data.frame")
adonis2(bray.dist ~ Material*Treatment,
        data = metadata)


# Bacteria - Relative abundance at Family level
merged_ps = merge_samples(ps.rar.16S, "Material")
normalized_ps_phylum = tax_glom(merged_ps, "Family", NArm = TRUE)

normalized_ps_phylum_relabun = transform_sample_counts(normalized_ps_phylum, function(OTU) OTU/sum(OTU))
taxaSums = data.frame(tax_table(normalized_ps_phylum_relabun)[,"Family"],
                      taxa_sums = taxa_sums(normalized_ps_phylum_relabun)) %>%
  arrange(desc(taxa_sums)) # reverse sort

N <- 30
top20 <- names(sort(taxa_sums(normalized_ps_phylum), decreasing = TRUE))[1:N]
GP.genus.prop <- transform_sample_counts(normalized_ps_phylum, function(x) x / sum(x) )
GP.genus.prop.top <- prune_taxa(top20, GP.genus.prop)
y = GP.genus.prop.top


df1 = data.frame(ID = c(taxa_names(y), "Other"), Family = c(tax_table(y)[,"Family"], "Other"))
df2 = t(cbind(otu_table(y), data.frame(Other = 1 - sample_sums(y))))

df = cbind(df1, df2) %>%
  pivot_longer(-c(ID, Family), names_to = "Material", values_to = "Abundance") %>%
  as.data.frame
df$Family = as.factor(df$Family)

# Save the data
write_csv(df,"/rel_ab_16S_2022_top30.csv")

colourCount = 31
getPalette = colorRampPalette(brewer.pal(12,"Set3"))

ggplot(df, aes(Material, Abundance, fill = Family)) +
  geom_col(color = "black") + scale_fill_manual(values = getPalette(colourCount)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Material", y = "Average relative abundance") + My_Theme + guides(fill = guide_legend(ncol = 1))



# Fungi - Relative abundance at Family level
merged_ps = merge_samples(ps.rar.ITS, "Material")
normalized_ps_phylum = tax_glom(merged_ps, "Family", NArm = TRUE)

normalized_ps_phylum_relabun = transform_sample_counts(normalized_ps_phylum, function(OTU) OTU/sum(OTU))
taxaSums = data.frame(tax_table(normalized_ps_phylum_relabun)[,"Family"],
                      taxa_sums = taxa_sums(normalized_ps_phylum_relabun)) %>%
  arrange(desc(taxa_sums)) # reverse sort

N <- 30
top20 <- names(sort(taxa_sums(normalized_ps_phylum), decreasing = TRUE))[1:N]
GP.genus.prop <- transform_sample_counts(normalized_ps_phylum, function(x) x / sum(x) )
GP.genus.prop.top <- prune_taxa(top20, GP.genus.prop)
y = GP.genus.prop.top


df1 = data.frame(ID = c(taxa_names(y), "Other"), Family = c(tax_table(y)[,"Family"], "Other"))
df2 = t(cbind(otu_table(y), data.frame(Other = 1 - sample_sums(y))))

df = cbind(df1, df2) %>%
  pivot_longer(-c(ID, Family), names_to = "Material", values_to = "Abundance") %>%
  as.data.frame
df$Family = as.factor(df$Family)

# Save the data
write_csv(df,"/rel_ab_ITS_2022_top30.csv")

colourCount = 31
getPalette = colorRampPalette(brewer.pal(12,"Set3"))

ggplot(df, aes(Material, Abundance, fill = Family)) +
  geom_col(color = "black") + scale_fill_manual(values = getPalette(colourCount)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Material", y = "Average relative abundance") + My_Theme + guides(fill = guide_legend(ncol = 1))

########################################################################


# Year comparison
# Figure 3

# Load data from each years
ps.16S.2020 <- readRDS("/ps_16S_2020_rarefied.rds")
ps.16S.2021 <- readRDS("/ps_16S_2021_rarefied.rds")
ps.16S.2022 <- readRDS("/ps_16S_2022_rarefied.rds")

ps.ITS.2020 <- readRDS("/ps_ITS_2020_rarefied.rds")
ps.ITS.2021 <- readRDS("/ps_ITS_2021_rarefied.rds")
ps.ITS.2022 <- readRDS("/ps_ITS_2022_rarefied.rds")
# Add the year in the metadata
# bacteria
sample_data(ps.16S.2020)$Year <- "2020"
sample_data(ps.16S.2021)$Year <- "2021"
sample_data(ps.16S.2022)$Year <- "2022"

# fungi
sample_data(ps.ITS.2020)$Year <- "2020"
sample_data(ps.ITS.2021)$Year <- "2021"
sample_data(ps.ITS.2022)$Year <- "2022"

# Merge data from each year
ps.16S.all <- merge_phyloseq(ps.16S.2021, ps.16S.2020, ps.16S.2022)
ps.ITS.all <- merge_phyloseq(ps.ITS.2021, ps.ITS.2020, ps.ITS.2022)

# Keep only the leaf samples
ps.16S.leaf <- subset_samples(ps.16S.all, Material == "Leaf")
ps.ITS.leaf <- subset_samples(ps.ITS.all, Material == "Leaf")

# Subset by management practices
ps.16S.mow <- subset_samples(ps.16S.leaf, Treatment == "Mowing")
ps.16S.graz <- subset_samples(ps.16S.leaf, Treatment == "Perm_access")
ps.16S.excl <- subset_samples(ps.16S.leaf, Treatment == "Perm_exclosure")

ps.ITS.mow <- subset_samples(ps.ITS.leaf, Treatment == "Mowing")
ps.ITS.graz <- subset_samples(ps.ITS.leaf, Treatment == "Perm_access")
ps.ITS.excl <- subset_samples(ps.ITS.leaf, Treatment == "Perm_exclosure")


#Alpha diversity
# Bacteria
plot_diversity_stats(ps.16S.root, group = "Year", 
                     index = "diversity_shannon", 
                     group.order = c("2020","2021","2022"),
                     label.format="p.format",
                     stats = TRUE) + ylab("Shannon Diversity") + xlab("")+ scale_fill_brewer( palette = "Pastel1")+ scale_color_brewer( palette = "Pastel1")+ My_Theme


# Fungi
plot_diversity_stats(ps.ITS.root, group = "Year", 
                     index = "diversity_shannon", 
                     group.order = c("2020","2021","2022"),
                     label.format="p.format",
                     stats = TRUE) + ylab("Shannon Diversity") + xlab("")+ scale_fill_brewer( palette = "Pastel1")+ scale_color_brewer( palette = "Pastel1")+ My_Theme



#Beta diversity
ps.prop.16S <- transform_sample_counts(ps.16S.root, function(otu) otu/sum(otu))
ps.prop.ITS <- transform_sample_counts(ps.ITS.root, function(otu) otu/sum(otu))

# Run PermANOVA
# bacteria
bray.dist <- phyloseq::distance(ps.prop.16S, method = "bray")
metadata <- as(sample_data(ps.prop.16S), "data.frame")
adonis2(bray.dist ~ Year,
        data = metadata)

# fungi
bray.dist <- phyloseq::distance(ps.prop.ITS, method = "bray")
metadata <- as(sample_data(ps.prop.ITS), "data.frame")
adonis2(bray.dist ~ Year,
        data = metadata)


# Beta diversity graphs
# bacteria
physeq.ord <- ordinate(ps.16S.root, "PCoA", "bray")
b.div.bray <- plot_ordination(ps.16S.root, physeq.ord, type= "samples", color= "Year") + geom_point(size=4)
b.div.bray <- b.div.bray + stat_ellipse()  + scale_color_brewer("Year", palette = "Pastel1") + My_Theme + theme_classic(base_size = 20)
print(b.div.bray)

#fungi
physeq.ord <- ordinate(ps.ITS.root, "PCoA", "bray")
b.div.bray <- plot_ordination(ps.ITS.root, physeq.ord, type= "samples", color= "Year") + geom_point(size=4)
b.div.bray <- b.div.bray + stat_ellipse()  + scale_color_brewer("Year", palette = "Pastel1") + My_Theme + theme_classic(base_size = 20)
print(b.div.bray)




# Relative abundance
# Fungi
merged_ps = merge_samples(ps.ITS.root, "Year")
normalized_ps_phylum = tax_glom(merged_ps, "Family", NArm = TRUE)

normalized_ps_phylum_relabun = transform_sample_counts(normalized_ps_phylum, function(OTU) OTU/sum(OTU))
taxaSums = data.frame(tax_table(normalized_ps_phylum_relabun)[,"Family"],
                      taxa_sums = taxa_sums(normalized_ps_phylum_relabun)) %>%
  arrange(desc(taxa_sums)) # reverse sort

N <- 20
top20 <- names(sort(taxa_sums(normalized_ps_phylum), decreasing = TRUE))[1:N]
GP.genus.prop <- transform_sample_counts(normalized_ps_phylum, function(x) x / sum(x) )
GP.genus.prop.top <- prune_taxa(top20, GP.genus.prop)
y = GP.genus.prop.top


df1 = data.frame(ID = c(taxa_names(y), "Other"), Family = c(tax_table(y)[,"Family"], "Other"))
df2 = t(cbind(otu_table(y), data.frame(Other = 1 - sample_sums(y))))

df = cbind(df1, df2) %>%
  pivot_longer(-c(ID, Family), names_to = "Year", values_to = "Abundance") %>%
  as.data.frame
df$Family = as.factor(df$Family)

colourCount = 21
getPalette = colorRampPalette(brewer.pal(8,"Pastel1"))

ggplot(df, aes(Year, Abundance, fill = Family)) +
  geom_col(color = "black") + scale_fill_manual(values = getPalette(colourCount)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Year", y = "Average relative abundance") + My_Theme + guides(fill = guide_legend(ncol = 1))



# Bacteria
merged_ps = merge_samples(ps.16S.root, "Year")
normalized_ps_phylum = tax_glom(merged_ps, "Family", NArm = TRUE)

normalized_ps_phylum_relabun = transform_sample_counts(normalized_ps_phylum, function(OTU) OTU/sum(OTU))
taxaSums = data.frame(tax_table(normalized_ps_phylum_relabun)[,"Family"],
                      taxa_sums = taxa_sums(normalized_ps_phylum_relabun)) %>%
  arrange(desc(taxa_sums)) # reverse sort

N <- 20
top20 <- names(sort(taxa_sums(normalized_ps_phylum), decreasing = TRUE))[1:N]
GP.genus.prop <- transform_sample_counts(normalized_ps_phylum, function(x) x / sum(x) )
GP.genus.prop.top <- prune_taxa(top20, GP.genus.prop)
y = GP.genus.prop.top


df1 = data.frame(ID = c(taxa_names(y), "Other"), Family = c(tax_table(y)[,"Family"], "Other"))
df2 = t(cbind(otu_table(y), data.frame(Other = 1 - sample_sums(y))))

df = cbind(df1, df2) %>%
  pivot_longer(-c(ID, Family), names_to = "Year", values_to = "Abundance") %>%
  as.data.frame
df$Family = as.factor(df$Family)

colourCount = 21
getPalette = colorRampPalette(brewer.pal(8,"Pastel1"))

ggplot(df, aes(Year, Abundance, fill = Family)) +
  geom_col(color = "black") + scale_fill_manual(values = getPalette(colourCount)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Year", y = "Average relative abundance") + My_Theme + guides(fill = guide_legend(ncol = 1))


################################################################


# Focus on leaves
# Managment practices differences
# Figure 4

# Bacteria - Alpha diversity
plot_diversity_stats(ps.rar.16S.root, group = "Treatment", 
                     index = "diversity_shannon", 
                     group.order = c("Mowing","Perm_access", "Perm_exclosure"),                      
                     #group.colors = mycols,
                     label.format="p.format",
                     stats = TRUE,
) + ylab("Shannon Diversity") + xlab("")+ scale_fill_brewer( palette = "Dark2") + scale_color_brewer(palette = "Dark2") + My_Theme + 
  scale_x_discrete(labels = c('Mowing','Grazing','Undisturbed vegetation'))


# Fungi - Alpha diversity
plot_diversity_stats(ps.rar.ITS.root, group = "Treatment", 
                     index = "diversity_shannon", 
                     group.order = c("Mowing","Perm_access", "Perm_exclosure"),                      
                     #group.colors = mycols,
                     label.format="p.format",
                     stats = TRUE,
) + ylab("Shannon Diversity") + xlab("")+ scale_fill_brewer( palette = "Dark2") + scale_color_brewer(palette = "Dark2") + My_Theme + 
  scale_x_discrete(labels = c('Mowing','Grazing','Undisturbed vegetation'))


# Bacteria - Beta diversity
physeq.ord <- ordinate(ps.rar.16S.root, "PCoA", "bray")
b.div.bray <- plot_ordination(ps.rar.16S.root, physeq.ord, type= "samples", color= "Treatment") + geom_point(size=4)
b.div.bray <- b.div.bray + theme_classic() + scale_color_brewer("Treatment", palette = "Dark2", name="Management practices", labels=c("Mowing","Grazing","Undisturbed vegetation")) + My_Theme
print(b.div.bray)

# Fungi - Beta diversity
physeq.ord <- ordinate(ps.rar.ITS.root, "PCoA", "bray")
b.div.bray <- plot_ordination(ps.rar.ITS.root, physeq.ord, type= "samples", color= "Treatment") + geom_point(size=4)
b.div.bray <- b.div.bray + theme_classic() + scale_color_brewer("Treatment", palette = "Dark2", name="Management practices", labels=c("Mowing","Grazing","Undisturbed vegetation")) + My_Theme
print(b.div.bray)

#Permanova Bacteria
ps.prop.16S <- transform_sample_counts(ps.rar.16S.leaf, function(otu) otu/sum(otu))
bray.dist <- phyloseq::distance(ps.prop.16S, method = "bray")
metadata <- as(sample_data(ps.prop.16S), "data.frame")
adonis2(bray.dist ~ Treatment,
        data = metadata)

#Permanova Fungi
ps.prop.ITS <- transform_sample_counts(ps.rar.ITS.leaf, function(otu) otu/sum(otu))
bray.dist <- phyloseq::distance(ps.prop.ITS, method = "bray")
metadata <- as(sample_data(ps.prop.ITS), "data.frame")
adonis2(bray.dist ~ Treatment,
        data = metadata)


# Differential abundance analysis - Leaves - Family level

# Group comparison
# bacteria
ps.leafMP.16S <- subset_samples(ps.rar.16S.leaf, Treatment=="Mowing" | Treatment=="Perm_access")
ps.leafEP.16S <- subset_samples(ps.rar.16S.leaf, Treatment=="Perm_exclosure" | Treatment=="Perm_access")
ps.leafEM.16S <- subset_samples(ps.rar.16S.leaf, Treatment=="Perm_exclosure" | Treatment=="Mowing")

# fungi
ps.leafMP.ITS <- subset_samples(ps.rar.ITS.leaf, Treatment=="Mowing" | Treatment=="Perm_access")
ps.leafEP.ITS <- subset_samples(ps.rar.ITS.leaf, Treatment=="Perm_exclosure" | Treatment=="Perm_access")
ps.leafEM.ITS <- subset_samples(ps.rar.ITS.leaf, Treatment=="Perm_exclosure" | Treatment=="Mowing")



# Comparison Grazing - Undisturbed vegetation
# Agglomerate to Genus/Family level
phy_genus <- tax_glom(ps.leafEP.ITS, "Family")

# Only keep genera present in at least 10 samples
phy_genus_pre <- ps_prune(phy_genus, min.samples = 10)

# Include covariate
res_covar <- DA.ds2(phy_genus_pre, predictor = "Treatment", covars = "Block")
res_covar[res_covar$pval.adj <= 0.05, ]

# Subset only highly significant
res_covar_sig <- res_covar[res_covar$pval.adj <= 0.01, ]

tab1 <- data.frame(res_covar_sig$Family, res_covar_sig$log2FoldChange)
names <- c("Family","log2FoldChange")
colnames(tab1)<- names

# change the direction of the sign to have Undisturbed vegetation as the base
tab1$log2FoldChange<--1*(tab1$log2FoldChange)

# add the kingdom
tab1$Kingdom <- "Fungi"


# Agglomerate to Genus/Family level
phy_genus <- tax_glom(ps.leafEP.16S, "Family")

# Only keep genera present in at least 10 samples
phy_genus_pre <- ps_prune(phy_genus, min.samples = 10)

# Include covariate
res_covar <- DA.ds2(phy_genus_pre, predictor = "Treatment", covars = "Block")
res_covar[res_covar$pval.adj <= 0.05, ]

# Subset only highly significant
res_covar_sig <- res_covar[res_covar$pval.adj <= 0.01, ]

tab2 <- data.frame(res_covar_sig$Family, res_covar_sig$log2FoldChange)
names <- c("Family","log2FoldChange")
colnames(tab2)<- names

# change the direction of the sign to have Undisturbed vegetation as the base
tab2$log2FoldChange<--1*(tab2$log2FoldChange)

# add the kingdom
tab2$Kingdom <- "Bacteria"

# Merge bacteria and fungi table
tab <- rbind(tab1,tab2)
tab$Kingdom <- factor(tab$Kingdom, levels = c("Fungi", "Bacteria"))

# Save the table
write_csv(tab, "/DSeq_grazing_leaf_2022_family.csv")


# Create a new .csv document after checking manually ANCOM BC results (for both mowing and grazing):
# selecting only logFC <-1 and >1 + significant adjusted p value
#Open this new document:
tab <- read_csv("/mowing_dseq_genus_pruned_2022.csv")
my_colors <- c("#FF6666","#0066CC")


My_Theme = theme(axis.title.x = element_text(size = 20),
                 axis.text.x = element_text(angle = 0, size = 20, hjust = 0.5),
                 axis.title.y = element_text(size = 20),
                 text = element_text(size = 15))

plt <- ggplot(graz) +
  geom_col(aes(logFC, Genus, fill = Kingdom), width = 0.6) + 
  theme_minimal() + My_Theme + scale_fill_manual(values = my_colors)

plt





# Comparison Mowing - Undisturbed vegetation
# Agglomerate to Genus/Family level
phy_genus <- tax_glom(ps.leafEM.ITS, "Family")

# Only keep genera present in at least 10 samples
phy_genus_pre <- ps_prune(phy_genus, min.samples = 10)

# Include covariate
res_covar <- DA.ds2(phy_genus_pre, predictor = "Treatment", covars = "Block")
res_covar[res_covar$pval.adj <= 0.05, ]

# Subset only highly significant
res_covar_sig <- res_covar[res_covar$pval.adj <= 0.01, ]

tab1 <- data.frame(res_covar_sig$Family, res_covar_sig$log2FoldChange)
names <- c("Family","log2FoldChange")
colnames(tab1)<- names

# change the direction of the sign to have Undisturbed vegetation as the base
tab1$log2FoldChange<--1*(tab1$log2FoldChange)

# add the kingdom
tab1$Kingdom <- "Fungi"


# Agglomerate to Genus/Family level
phy_genus <- tax_glom(ps.leafEM.16S, "Family")

# Only keep genera present in at least 10 samples
phy_genus_pre <- ps_prune(phy_genus, min.samples = 10)

# Include covariate
res_covar <- DA.ds2(phy_genus_pre, predictor = "Treatment", covars = "Block")
res_covar[res_covar$pval.adj <= 0.05, ]

# Subset only highly significant
res_covar_sig <- res_covar[res_covar$pval.adj <= 0.01, ]

tab2 <- data.frame(res_covar_sig$Family, res_covar_sig$log2FoldChange)
names <- c("Family","log2FoldChange")
colnames(tab2)<- names

# change the direction of the sign to have Undisturbed vegetation as the base
tab2$log2FoldChange<--1*(tab2$log2FoldChange)

# add the kingdom
tab2$Kingdom <- "Bacteria"

# Merge bacteria and fungi table
tab <- rbind(tab1,tab2)
tab$Kingdom <- factor(tab$Kingdom, levels = c("Fungi", "Bacteria"))

write_csv(tab, "")

tab <- read.csv("/media/tanneau/601C-60FB/Doctorat/Publications/1st publication/Supp. data/ANCOM-BC/Final_sheets/2020_mowing_ANCOMBC_logFC_padj.csv")

my_colors <- c("#FF6666","#0066CC")
my_colors <- c("#0066CC")

plt <- ggplot(tab) +
  geom_col(aes(logFC, forcats::fct_reorder(Genus, logFC), fill = Kingdom), width = 0.6) + 
  theme_minimal() + My_Theme + scale_fill_manual(values = my_colors)+labs(y= "Genus", x = "log fold-change")

plt



###########################################################


# Network building
# Figure 5

# Merge ITS and 16S data (or open the ps. object with both fungi and bacteria data)
ps <- readRDS("/ps_2022_16S_ITS.rds")

# Subset to keep only the data from the leaves
ps.leaf <- subset_samples(ps, Material=="Leaf")

# Subset by management practices
ps.mow.leaf <- subset_samples(ps.leaf, Treatment=="Mowing")
ps.graz.leaf <- subset_samples(ps.leaf, Treatment=="Perm_access")
ps.con.leaf <- subset_samples(ps.leaf, Treatment=="Perm_exclosure")

# Keeping taxa with taxa sum >5
ps.mow.leaf <- prune_taxa(taxa_sums(ps.mow.leaf) > 5, ps.mow.leaf)
ps.graz.leaf <- prune_taxa(taxa_sums(ps.graz.leaf) > 5, ps.graz.leaf)
ps.con.leaf <- prune_taxa(taxa_sums(ps.con.leaf) > 5, ps.con.leaf)





# Construction network for Mowing treatment
amgut_genus <- tax_glom(ps.con.leaf, taxrank = "Genus")

#Taxonomic table
taxtab <- as(tax_table(amgut_genus), "matrix")

amgut_genus_renamed <- renameTaxa(amgut_genus, 
                                  pat = "<name>", 
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus")

net_genus <- netConstruct(amgut_genus_renamed,
                          taxRank = "Genus",
                          measure = "spieceasi",
                          zeroMethod = "multRepl",
                          normMethod = "clr",
                          sparsMethod = "threshold",
                          thresh = 0.3,
                          verbose = 3)


props_genus <- netAnalyze(net_genus, clustMethod = "cluster_fast_greedy")
# Metrics information
summary(props_genus, numbNodes = 50L)

#phylcol <- c("#984EA3","#FF7F00")
phylcol <- c("#FF6666","#0066CC")



#library(DAAG)
# Get phyla names
taxtab <- as(tax_table(amgut_genus_renamed), "matrix")
phyla <- as.factor(gsub("k__", "", taxtab[, "Kingdom"]))
names(phyla) <- taxtab[, "Genus"]

# Generate the plot
plot(props_genus,
     edgeInvisFilter = "threshold",
     edgeInvisPar = 0.07, #edge filtering
     layout = "spring",
     repulsion = 1, # shape with closest node
     shortenLabels = "none",
     #charToRm = "g__",
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr", #size of nodes based on number of counts
     nodeSizeSpread = 4, #size of big nodes
     nodeColor = "feature", 
     featVecCol = phyla, 
     colorVec =  phylcol,
     posCol = "#66CC99", 
     negCol = "#CC3366",
     edgeTranspLow = 1,
     edgeTranspHigh = 40, #transparency
     cexNodes = 2, # size nodes
     cexLabels = 0, # 0 = no ; 1=name of every node
     # cexLabels as 1 = all node names
     cexHubLabels = NULL, # name of the hubs - NULL or size 0.9
     title1 = "Network on genus level with SpiecEasi - Undisturbed vegetation - Leaf", 
     showTitle = TRUE,
     cexTitle = 1.2,
     highlightHubs = TRUE,
     hubLabelFont = NULL,
     hubBorderCol = "black")

# Add legend
phylcol_transp <- colToTransp(phylcol, 40)

legend(-1, 1, cex = 0.85, pt.cex = 3, title = "Kingdom", 
       legend=levels(phyla), col = phylcol_transp, bty = "n", pch = 16) 

legend(0.60, 1, cex = 0.9, title = "estimated correlation:",
       legend = c("+","-"), lty = 1.9, lwd = 2, col = c("#66CC99","#CC3366"), 
       bty = "n", horiz = TRUE)







# Construction network for Undisturbed vegetation treatment
amgut_genus <- tax_glom(ps.con.leaf, taxrank = "Genus")

#Taxonomic table
taxtab <- as(tax_table(amgut_genus), "matrix")

amgut_genus_renamed <- renameTaxa(amgut_genus, 
                                  pat = "<name>", 
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus")

net_genus <- netConstruct(amgut_genus_renamed,
                          taxRank = "Genus",
                          measure = "spieceasi",
                          zeroMethod = "multRepl",
                          normMethod = "clr",
                          sparsMethod = "threshold",
                          thresh = 0.3,
                          verbose = 3)

props_genus <- netAnalyze(net_genus, clustMethod = "cluster_fast_greedy")
# Metrics information
summary(props_genus, numbNodes = 50L)

# Get phyla names
taxtab <- as(tax_table(amgut_genus_renamed), "matrix")
phyla <- as.factor(gsub("k__", "", taxtab[, "Kingdom"]))
names(phyla) <- taxtab[, "Genus"]

# Generate the plot
plot(props_genus,
     edgeInvisFilter = "threshold",
     edgeInvisPar = 0.07, #edge filtering
     layout = "spring",
     repulsion = 1, # shape with closest node
     shortenLabels = "none",
     #charToRm = "g__",
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr", #size of nodes based on number of counts
     nodeSizeSpread = 4, #size of big nodes
     nodeColor = "feature", 
     featVecCol = phyla, 
     colorVec =  phylcol,
     posCol = "#66CC99", 
     negCol = "#CC3366",
     edgeTranspLow = 1,
     edgeTranspHigh = 40, #transparency
     cexNodes = 2, # size nodes
     cexLabels = 0, # 0 = no ; 1=name of every node
     # cexLabels as 1 = all node names
     cexHubLabels = NULL, # name of the hubs - NULL or size 0.9
     title1 = "Network on genus level with SpiecEasi - Undisturbed vegetation - Leaf", 
     showTitle = TRUE,
     cexTitle = 1.2,
     highlightHubs = TRUE,
     hubLabelFont = NULL,
     hubBorderCol = "black")


# Add legend
phylcol_transp <- colToTransp(phylcol, 40)

legend(-1, 1, cex = 0.85, pt.cex = 3, title = "Kingdom", 
       legend=levels(phyla), col = phylcol_transp, bty = "n", pch = 16) 

legend(0.60, 1.05, cex = 0.9, title = "estimated correlation:",
       legend = c("+","-"), lty = 1.9, lwd = 2, col = c("#66CC99","#CC3366"), 
       bty = "n", horiz = TRUE)





# Construction network for Grazing treatment
amgut_genus <- tax_glom(ps.graz.leaf, taxrank = "Genus")

#Taxonomic table
taxtab <- as(tax_table(amgut_genus), "matrix")

amgut_genus_renamed <- renameTaxa(amgut_genus, 
                                  pat = "<name>", 
                                  substPat = "<name>_<subst_name>(<subst_R>)",
                                  numDupli = "Genus")

net_genus <- netConstruct(amgut_genus_renamed,
                          taxRank = "Genus",
                          measure = "spieceasi",
                          zeroMethod = "multRepl",
                          normMethod = "clr",
                          sparsMethod = "threshold",
                          thresh = 0.3,
                          verbose = 3)

props_genus <- netAnalyze(net_genus, clustMethod = "cluster_fast_greedy")
# Metrics information
summary(props_genus, numbNodes = 50L)

# Get phyla names
taxtab <- as(tax_table(amgut_genus_renamed), "matrix")
phyla <- as.factor(gsub("k__", "", taxtab[, "Kingdom"]))
names(phyla) <- taxtab[, "Genus"]

# Generate the plot
plot(props_genus,
     edgeInvisFilter = "threshold",
     edgeInvisPar = 0.07, #edge filtering
     layout = "spring",
     repulsion = 1, # shape with closest node
     shortenLabels = "none",
     #charToRm = "g__",
     labelScale = FALSE,
     rmSingles = TRUE,
     nodeSize = "clr", #size of nodes based on number of counts
     nodeSizeSpread = 4, #size of big nodes
     nodeColor = "feature", 
     featVecCol = phyla, 
     colorVec =  phylcol,
     posCol = "#66CC99", 
     negCol = "#CC3366",
     edgeTranspLow = 1,
     edgeTranspHigh = 40, #transparency
     cexNodes = 2, # size nodes
     cexLabels = 0, # 0 = no ; 1=name of every node
     # cexLabels as 1 = all node names
     cexHubLabels = NULL, # name of the hubs - NULL or size 0.9
     title1 = "Network on genus level with SpiecEasi - Grazing - Leaf", 
     showTitle = TRUE,
     cexTitle = 1.2,
     highlightHubs = TRUE,
     hubLabelFont = NULL,
     hubBorderCol = "black")


# Add legend
phylcol_transp <- colToTransp(phylcol, 40)

legend(-1, 1, cex = 0.85, pt.cex = 3, title = "Kingdom", 
       legend=levels(phyla), col = phylcol_transp, bty = "n", pch = 16) 

legend(0.60, 1, cex = 0.9, title = "estimated correlation:",
       legend = c("+","-"), lty = 1.9, lwd = 2, col = c("#66CC99","#CC3366"), 
       bty = "n", horiz = TRUE)

