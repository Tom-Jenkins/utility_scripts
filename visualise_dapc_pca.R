# --------------------------- #
#
# Description:
# Visualise DAPC / PCA results
#
# DAPC: discriminant analysis of principal components
# PCA: principal components analysis
#
# --------------------------- #

# Load libraries
library(adegenet)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(scales)

# Import example dataset
data("dapcIllus")
geno = dapcIllus$a

# Data summary
geno
summary(geno$pop)

# Change population labels
popNames(geno) = c("Pop1","Pop2","Pop3","Pop4","Pop5","Pop6")
summary(geno$pop)


#--------------#
#
# Perform DAPC
#
#--------------#

# Perform cross validation to find the optimal number of PCs to retain in DAPC
set.seed(123)
x = tab(geno, NA.method = "mean")
crossval = xvalDapc(x, geno$pop, result = "groupMean", xval.plot = TRUE)

# Number of PCs with best stats (lower score = better)
crossval$`Root Mean Squared Error by Number of PCs of PCA`
crossval$`Number of PCs Achieving Highest Mean Success`
crossval$`Number of PCs Achieving Lowest MSE`
numPCs = as.numeric(crossval$`Number of PCs Achieving Lowest MSE`)

# Run a DAPC using population IDs as priors
dapc1 = dapc(geno, geno$pop, n.pca = numPCs, n.da = 3)
dapc1

# Analyse how much percent of genetic variance is explained by each axis
percent = dapc1$eig/sum(dapc1$eig)*100
percent
barplot(percent, ylab = "Percent of genetic variance explained by eigenvectors", 
        names.arg = round(percent, 2))

# Create a dataframe containing individual coordinates
ind_coords = as.data.frame(dapc1$ind.coord)


#--------------#
#
# Perform PCA
#
#--------------#

# Replace missing data with the mean allele frequencies
x = tab(geno, NA.method = "mean")
x

# Perform PCA and plot results
pca1 = dudi.pca(x, scannf = FALSE, scale = FALSE, nf = 3)
pca1

# Analyse how much percent of genetic variance is explained by each axis
percent = pca1$eig/sum(pca1$eig)*100
percent
barplot(percent, ylab = "Percent of genetic variance explained by eigenvectors",
        names.arg = round(percent, 2))

# Create a dataframe containing individual coordinates
ind_coords = as.data.frame(pca1$li)


#--------------#
#
# Visualise DAPC results
#
#--------------#

# Rename column of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")

# Add a column containing individuals
ind_coords$Ind = indNames(geno)

# Add a column with the population IDs
ind_coords$Pop = geno$pop

# Conditional function that adds other labels to dataframe (optional)
addlabel = function(x){
  # If population label x is present function will output y
  if(x=="Pop1"|x=="Pop2"|x=="Pop3") y = "Region1"
  if(x=="Pop4"|x=="Pop5"|x=="Pop6") y = "Region2"
return(y)
}
ind_coords$Region = sapply(ind_coords$Pop, addlabel)  
head(ind_coords)

# Reorder labels for plotting (optional)
# ind_coords$Region = factor(ind_coords$Region, levels = c("Region2","Region1"))

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Pop,
                     data = ind_coords,
                     FUN = mean)

# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Pop", suffix = c("",".cen"))
head(ind_coords)

# Add region labels to centroid dataframe
centroid$Region = sapply(centroid$Pop, addlabel)  
centroid

# Define colour palette
show_col(brewer.pal(9, "Set1"))
cols = brewer.pal(nPop(geno), "Set1") ; cols

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Custom ggplot2 theme
ggtheme = theme(legend.title = element_blank(),
                axis.text.y = element_text(colour="black", size=14),
                axis.text.x = element_text(colour="black", size=14),
                axis.title = element_text(colour="black", size=14),
                legend.position = "right",
                legend.text = element_text(size=15),
                legend.key = element_rect(fill = NA),
                legend.key.size = unit(0.7, "cm"),
                legend.box.spacing = unit(0, "cm"),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                # title centered
                plot.title = element_text(hjust=0.5, size=25) 
)

# Scatter plot axis 1 vs. 2 [colour by population]
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Pop), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Pop), shape = 21, size = 4, show.legend = TRUE)+
  # centroids
  geom_label(data = centroid, aes(label = Pop, fill = Pop), size = 5, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  # custom theme
  ggtheme
  
# Export plot
ggsave("Figure1.png", width = 12, height = 8, dpi = 600)

# Scatter plot axis 1 vs. 2 [colour by region]
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # spider segments
  geom_segment(aes(xend = Axis1.cen, yend = Axis2.cen, colour = Region), show.legend = FALSE)+
  # points
  geom_point(aes(fill = Region), shape = 21, size = 4, show.legend = TRUE)+
  # centroids
  geom_label(data = centroid, aes(label = Pop, fill = Region), size = 5, show.legend = FALSE)+
  # colouring
  scale_fill_manual(values = cols)+
  scale_colour_manual(values = cols)+
  # custom labels
  labs(x = xlab, y = ylab)+
  # custom theme
  ggtheme

# Export plot
ggsave("Figure2.png", width = 12, height = 8, dpi = 600)