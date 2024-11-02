
setwd("C:/Users/szd_z/OneDrive/Desktop/WGCNA/webtutorial")

#Setting string not as factor
options(stringsAsFactors = FALSE)

#Loading WGCNA
library(WGCNA)
options(stringsAsFactors = FALSE)
#Enable multithread
enableWGCNAThreads()
#Reading the raw data (rows are the sample and columns the genes)
expressiondata = read.csv("counts.csv", header = TRUE, row.names = 1)
expression <- t(expressiondata)

#Group data in a dendogram to check outliers
sampleTree = hclust(dist(expression), method = "average")
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

#Plot a line showing the cut-off
abline(h = 2e+06, col = "red") #This value of 31000 was chosen based on my data, you need to check the best value to your data
######################################

traitData = read.csv("metadata.csv", header = TRUE, row.names = 1)
#########################################
gsg <- goodSamplesGenes(expression)
summary(gsg)
gsg$allOK
#Removing outliers
table(gsg$goodGenes)
table(gsg$goodSamples)
data <- expressiondata[gsg$goodGenes == TRUE,]
#detecting outliers Samples hierachical clusstring method1

htree <- hclust(dist(t(data)), method = 'average')
plot(htree)

#####

############## part 2 ##################################################
#revoming outliers
head(data)
samples.to.be.exculed <- c('danaus_leg6')
data.subst <- data[,!(colnames(data)%in% samples.to.be.exculed)]
#=====================================================
library(tidyverse)
library(DESeq2)
library(WGCNA)
metadata <- traitData %>%
  filter(!row.names(.) %in% samples.to.be.exculed)
################# metadata and subset #### done

#making the rownames and colnames identical
all(rownames(metadata) %in% colnames(data.subst))
dim(data.subst)
dim(metadata)
all(rownames(metadata) == colnames(data.subst))

htree <- hclust(dist(t(data.subst)), method = 'average')
plot(htree)
#Plot a line showing the cut-off
abline(h = 2e+06, col = "red") #This value of 31000 was chosen based on my data, you need to check the best value to your data
########################################################
###########Removing NA values
data.subst <- data.subst %>%
  mutate_all(~ ifelse(is.na(.), 1, .))


########################################
#normalization with deseq2
dds <- DESeqDataSetFromMatrix(countData = data.subst,
                              colData = metadata,
                              design = ~ 1)

dds75 <- dds[rowSums(counts(dds) >= 10) >= 06,]
nrow(dds75)
dds_norm <- vst(dds75)
norm_counts <- assay(dds_norm) %>%
  t()
################################
#Network Construction.
Power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

sft <- pickSoftThreshold(norm_counts,
                         powerVector = Power,
                         networkType = "signed",
                         verbose = 5)
#####################
#######################
sft.data <- sft$fitIndices
library(ggplot2)
#visualize to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red')+
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()
#installing grid
#install.packages("gridExtra")
library(gridExtra)

grid.arrange(a1, a2, nrow = 2)
#grid.arrange(a1, a2, ncol = 2)
##########################################value is 18

#convert matrix to numeric
norm_counts[] <- sapply(norm_counts, as.numeric)
#norm_counts <-  t(norm_counts)
soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor


#memory estimates  block size

bwnet <- blockwiseModules(norm_counts,
                          maxBlockSize = 6000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

#module eigen genes pressent in 
module_eigengenes <- bwnet$MEs
head(module_eigengenes)
#get numbers of genes in each module
color_table <- table(bwnet$colors)
color_table <- data.frame(color_table)
write.csv(color_table, "color_table.csv")

#plot the dendogram and module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)

##########################################
#################################
#####################

#creating trait file
head(metadata)

meta_bin <- metadata %>%
  mutate(Leg = ifelse(Tissues == "leg", 1, 0),
         Antenae = ifelse(Tissues == "Antenae", 1, 0),
         Tribolium_castaneum = ifelse(Species == "Tribolium castaneum", 1, 0),
         Danaus_plexipus = ifelse(Species == "Danaus plexipus", 1, 0),
         Heliconius_melpomene = ifelse(Species == "heliconius melpomene", 1, 0))

head(meta_bin)
test <- meta_bin[-1]
test1 <- test[-1]
test2 <- test1[-1]



nSamples <- nrow(norm_counts)
nGenes <- ncol(norm_counts)

nrow(module_eigengenes)
nrow(test2)

module.trait.cor <- cor(module_eigengenes, test2, use = 'p')
module.trait.cor.pvals <- corPvalueStudent(module.trait.cor, nSamples)
# visualize module-trait association as a heat map

# Install the package if you haven't already
#install.packages("CorLevelPlot")

# Load the package
library(CorLevelPlot)


heatmap.data <- merge(module_eigengenes, test2, by = 'row.names')
head(heatmap.data)

heatmap.data <- heatmap.data %>%
  column_to_rownames(var = 'Row.names')

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[09:10],
             y = names(heatmap.data)[01:08],
             col = c("blue1", "skyblue", "white", "pink", "red"))

# extracting associated modules
module.gene.mapping <- as.data.frame(bwnet$colors)

#extracting red module
red.module <- module.gene.mapping %>%
  filter(bwnet$colors == 'red') %>%
  rownames()

red.module
red.module <- data.frame(red.module)
write.csv(red.module, "redgene.csv")
##############done#################Intra-modular analysis ########################



# geting hub genes
# getting association between eigengenes and trait
#this will quanitfy the similarity of all genes in array to each module

module.membership.measure <- cor(module_eigengenes, norm_counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

#we can look at which genes have high pvalue  for each module 
#this is needed to know hub gene in interesting module show this command will show you 1to 7 genes
module.membership.measure.pvals[1:07,1:07]

#calculate the gene significance and associated p-values gene significance
gene.signf.corr <- cor(norm_counts, test2$Leg, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

#you are now looking at top 30 genes here associated with leg
gene.signf.corr.pvals %>%
  as.data.frame() %>%
  arrange(V1) %>%
  head(30)
#top 30 genes significantly acoociated with leg
# now can convert them into gene symbols and analyze downstream 
##### making red module plot for both tissues
awon <- merge(module_eigengenes, metadata, by = 0)

ggplot(
  awon,
  aes(
    x = Tissues,
    y = MEred,
    color = Tissues
  )
) +
  # a boxplot with outlier points hidden (they will be in the sina plot)
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  # A sina plot to show all of the individual data points
  ggforce::geom_sina(maxwidth = 0.3) +
  theme_classic()
##########################

modNames = substring(names(MEs), 3)

###### #### next training


module_df <- data.frame(
  gene_id = names(bwnet$colors),
  colors = bwnet$colors
)





# pick out a few modules of interest here
modules_of_interest = c("red")


# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
norm_counts[1:5,1:10]
norm_counts_T <- t(norm_counts)

subexpr = norm_counts_T[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )
###

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")

################
##visualizaiton of network
geneTree = bwnet$dendrograms[[1]]
moduleColors = bwnet$colors

## Calculate topological overlap anew: this could be done more efficiently by saving the TOM
## calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(norm_counts, power = 18);
## Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
## Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
## Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

nSelect = 400
## For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
## There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
## Open a graphical window
sizeGrWindow(9,9)
## Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
## the color palette; setting the diagonal to NA also improves the clarity of the plot
##png('123.png')
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
### done #### network plotting

####### exproting files for red module 
# pick out a few modules of interest here only red
modules_of_interest = c("red")
genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = norm_counts_T[genes_of_interest$gene_id,]
expr_of_interest[1:5,1:5]
#>                       B-3      B-4      B-5      L-3      L-4
#> AC186512.3_FG001 6.901539 7.389644 6.975945 6.859593 7.370816
#> AC186512.3_FG007 7.919688 7.754506 7.670946 7.417760 7.988427
#> AC190623.3_FG001 6.575155 7.170788 7.438024 8.223261 8.008850
#> AC196475.3_FG004 6.054319 6.439899 6.424540 5.815344 6.565299
#> AC196475.3_FG005 6.194406 5.872273 6.207174 6.499828 6.314952

# Only recalculate TOM for modules of interest (faster, altho there's some online discussion if this will be slightly off)
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = soft_power)

# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

head(edge_list)
#> # A tibble: 6 x 5
#>   gene1            gene2            correlation module1   module2  
#>   <chr>            <chr>                  <dbl> <chr>     <chr>    
#> 1 AC186512.3_FG001 AC186512.3_FG007      0.0238 turquoise turquoise
#> 2 AC186512.3_FG001 AC190623.3_FG001      0.0719 turquoise turquoise
#> 3 AC186512.3_FG001 AC196475.3_FG004      0.143  turquoise turquoise
#> 4 AC186512.3_FG001 AC196475.3_FG005      0.0117 turquoise turquoise
#> 5 AC186512.3_FG001 AC196489.3_FG002      0.0181 turquoise turquoise
#> 6 AC186512.3_FG001 AC198481.3_FG004      0.0240 turquoise turquoise

# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list,
            file = "edgelist.tsv",
            delim = "\t")








