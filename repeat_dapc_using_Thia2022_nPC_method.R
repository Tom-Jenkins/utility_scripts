# -------------------------- #
#
# Repeat DAPC using Thia 2022 Recommendations
#
# Author: Tom Jenkins
# R Version: 4.0.2
# Last updated: November 2022
#
# -------------------------- #

# Thia 2022 paper
# https://doi.org/10.1111/1755-0998.13706

# Import packages
library(tidyverse) # v1.3.1
library(adegenet) # v2.1.5
library(poppr) # v2.9.3

# Download lobster SNP genotypes from Jenkins et al. 2019 to script directory
# https://doi.org/10.1111/eva.12849
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.2v1kr38
url = "https://tomjenkins.netlify.app/files/Lobster_SNP_Genotypes.csv"
download.file(url, "./Lobster_SNP_Genotypes.csv")


#--------------#
#
# Import Data & Filtering
#
#--------------#

# Apply the same code and filters as in the paper
# https://tomjenkins.netlify.app/tutorials/r-popgen-getting-started/#3

lobster = read.csv("Lobster_SNP_Genotypes.csv")
lobster_wide = reshape(lobster, idvar = c("ID","Site"), timevar = "Locus", direction = "wide", sep = "")
colnames(lobster_wide) = gsub("Genotype", "", colnames(lobster_wide))
snpgeno = lobster_wide[ , 3:ncol(lobster_wide)]
snps_to_remove = c("25580","32362","41521","53889","65376","8953","21197","15531","22740","28357","33066","51507","53052","53263","21880","22323","22365")
snpgeno = snpgeno[ , !colnames(snpgeno) %in% snps_to_remove]
ind = as.character(lobster_wide$ID) # individual ID
site = as.character(lobster_wide$Site) # site ID
lobster_gen = df2genind(snpgeno, ploidy = 2, ind.names = ind, pop = site, sep = "")
lobster_gen = missingno(lobster_gen, type = "geno", cutoff = 0.20)
lob_dups = c("Laz4","Eye15","Eye16","Eye01","Laz2","Eye08","Gul101","Eye25","Iom02","Hel07","Eye27","Eye05","Eye06","Eye23","Eye22","Eye11","Cro08","Tar1","Eye14","Tar3","Lyn04","Lyn15","Eye07","Eye02","Eye20")
lob_Nodups = indNames(lobster_gen)[! indNames(lobster_gen) %in% lob_dups]
lobster_gen = lobster_gen[lob_Nodups, ]
isPoly(lobster_gen) %>% summary
lobster_gen


#--------------#
#
# Plotting Parameters
#
#--------------#

# Colours
cols = c("#377EB8","#FDB462","#E31A1C","#FC8D59")
scales::show_col(cols)

# Custom theme for ggplot2
ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                plot.title = element_text(size=15),
                plot.subtitle = element_text(size=10)
)

# Function to add region labels to data frame
add_region_labs = function(df){
  newDF = df %>% 
    mutate(
      Region = case_when(
        Site == "Ale"|Site == "The"|Site == "Tor"|Site == "Sky" ~ " Aegean Sea ",
        Site == "Sar13"|Site == "Sar17"|Site == "Tar"|Site == "Laz" ~ " Central Mediterranean ",
        Site == "Jer"|Site == "Idr16"|Site == "Idr17"|Site == "Hel" ~ " Atlantic ",
        Site == "Cor"|Site == "Hoo"|Site == "Kil"|Site == "Mul" ~ " Atlantic ",
        Site == "Ven"|Site == "Tro"|Site == "Ber"|Site == "Flo" ~ " Atlantic ",
        Site == "Sin"|Site == "Gul"|Site == "Kav"|Site == "Lys" ~ " Atlantic ",
        Site == "Vig"|Site == "Brd"|Site == "Cro"|Site == "Eye" ~ " Atlantic ",
        Site == "Heb"|Site == "Iom"|Site == "Ios"|Site == "Loo" ~ " Atlantic ",
        Site == "Lyn"|Site == "Ork"|Site == "Pad"|Site == "Pem" ~ " Atlantic ",
        Site == "She"|Site == "Sbs"|Site == "Sul" ~ " Atlantic ",
        Site == "Oos" ~ " Oosterschelde "
      )
    )
  return(newDF)
}

# Order of region labels
region_order = c(" Atlantic ", " Central Mediterranean ", " Aegean Sea ", " Oosterschelde ")


#--------------#
#
# Standard PCA
#
#--------------#

# Replace missing data with the mean allele frequencies
x = tab(lobster_gen, NA.method = "mean")

# Perform PCA
pca1 = dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 3)

# Analyse how much percent of genetic variance is explained by each axis
percent = pca1$eig/sum(pca1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,12),
        names.arg = round(percent, 1))

# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(pca1$li)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")

# Add a column containing individuals
ind_coords$Ind = indNames(lobster_gen)

# Add a column with the site IDs
ind_coords$Site = lobster_gen$pop

# Add regional labels to data frame
ind_coords = add_region_labs(ind_coords)

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)
centroid = add_region_labs(centroid)
centroid$Region = factor(centroid$Region, levels = region_order)

# Add centroid coordinates to ind_coords data frame
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))
ind_coords$Region = factor(ind_coords$Region, levels = region_order)

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Scatter plot axis 1 vs. 2
pca_plt = ggplot(data = ind_coords, aes(x = Axis1*-1, y = Axis2*-1))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen*-1, yend = Axis2.cen*-1, colour = Region), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Region), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Region), size = 4, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("European lobster PCA", subtitle = "1,271 individuals genotyped at 79 biallelic SNPs")+
  # custom theme
  ggtheme
pca_plt

# Export plot
ggsave(plot = pca_plt, filename = "Lobster_PCA.png", width = 12, height = 8, dpi = 600)


#--------------#
#
# DAPC: Original Method in Jenkins et al. 2019
#
#--------------#

# xvalDAPC() was used to assess the number of principal components to use in DA
# 70 PCs

# Cross validation to find the optimal number of PCs to retain in dapc
# x = tab(lobster_gen, NA.method = "mean")
# set.seed(123)
# xval = xvalDapc(x, lobster_gen$pop, n.pca.max=300, training.set=0.9,
#                 result="groupMean", center=TRUE, scale=FALSE,
#                 n.rep=30, n.pca=NULL)

# Number of PCs with best stats
# xval$`Number of PCs Achieving Highest Mean Success`
# xval$`Number of PCs Achieving Lowest MSE`
# xval$`Root Mean Squared Error by Number of PCs of PCA` # lower score = better
# 70 PCs to retain

# Run a DAPC using site IDs as priors
dapc1 = dapc(lobster_gen, lobster_gen$pop, n.pca = 70, n.da = 3)

# Analyse how much percent of genetic variance is explained by each axis
percent = dapc1$eig/sum(dapc1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,60),
        names.arg = round(percent, 1))

# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(dapc1$ind.coord)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")

# Add a column containing individuals
ind_coords$Ind = indNames(lobster_gen)

# Add a column with the site IDs
ind_coords$Site = lobster_gen$pop

# Add regional labels to data frame
ind_coords = add_region_labs(ind_coords)

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)
centroid = add_region_labs(centroid)
centroid$Region = factor(centroid$Region, levels = region_order)

# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))
ind_coords$Region = factor(ind_coords$Region, levels = region_order)

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Scatter plot axis 1 vs. 2
dapc_plt1 = ggplot(data = ind_coords, aes(x = Axis1*-1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen*-1, yend = Axis2.cen, colour = Region), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Region), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Region), size = 4, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("European lobster DAPC: Retaining 70 PCs as in Jenkins et al. 2019 (xvalDapc)", subtitle = "1,271 individuals genotyped at 79 biallelic SNPs")+
  # custom theme
  ggtheme
dapc_plt1

# Export plot
ggsave(plot = dapc_plt1, filename = "Lobster_DAPC_70_PCs.png", width = 12, height = 8, dpi = 600)



#--------------#
#
# DAPC: K-1 Method in Thia 2022
#
#--------------#

# Because of the genetic cline in the Atlantic, setting K was difficult
# However, based on the PCA and on previous studies and known information about lobster,
# K was set to 5. Therefore, 4 principal components were retained in the DAPC.

# Run a DAPC using site IDs as priors
dapc2 = dapc(lobster_gen, lobster_gen$pop, n.pca = 4, n.da = 3)

# Analyse how much percent of genetic variance is explained by each axis
percent = dapc2$eig/sum(dapc2$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,60),
        names.arg = round(percent, 1))

# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(dapc2$ind.coord)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")

# Add a column containing individuals
ind_coords$Ind = indNames(lobster_gen)

# Add a column with the site IDs
ind_coords$Site = lobster_gen$pop

# Add regional labels to data frame
ind_coords = add_region_labs(ind_coords)

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Site, data = ind_coords, FUN = mean)
centroid = add_region_labs(centroid)
centroid$Region = factor(centroid$Region, levels = region_order)

# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Site", suffix = c("",".cen"))
ind_coords$Region = factor(ind_coords$Region, levels = region_order)

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Scatter plot axis 1 vs. 2
dapc_plt2 = ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Region), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Region), shape = 21, size = 3, show.legend = FALSE)+
  # centroids
  geom_label(data = centroid, aes(label = Site, fill = Region), size = 4, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  ggtitle("European lobster DAPC: Retaining 4 PCs using Thia 2022 K-1 Method", subtitle = "1,271 individuals genotyped at 79 biallelic SNPs")+
  # custom theme
  ggtheme
dapc_plt2

# Export plot
ggsave(plot = dapc_plt2, filename = "Lobster_DAPC_4_PCs.png", width = 12, height = 8, dpi = 600)


# Combine plots
library(patchwork) # v1.1.1
(dapc_fig = dapc_plt1 + dapc_plt2)
ggsave(plot = dapc_fig, filename = "Lobster_DAPC.png", width = 18, height = 8, dpi = 600)
