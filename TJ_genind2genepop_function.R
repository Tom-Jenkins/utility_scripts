#======================================
# R function to export a genind object in genepop format
# 
# Tom Jenkins t.l.jenkins@exeter.ac.uk
#
# July 2018
#
#======================================

## data: genind object
## file: file name to write
## Function only works on biallelic SNP datasets.
## Missing data must be recorded as NA.
## Must have adegenet and miscTools installed

## Example use: 
# library(adegenet)
# ind = as.character(paste("ind_", seq(1:100), sep=""))
# pop = as.character(c(rep("Pop1",25), rep("pop2",25), rep("pop3",25), rep("pop4",25)))
# loci = list(c("AA","AC","CC"), c("GG","GC","CC"), c("TT","TA","AA"), c("CC","CT","TT"))
# loci = sample(loci, 100, replace=T)
# loci = lapply(loci, sample, size=100, replace=TRUE)
# geno = as.data.frame(loci, col.names= .genlab("loc",100))
# gen = df2genind(geno, ploidy=2, ind.names=ind, pop=pop, sep="")
# genind2genepop(gen, file="example_genepop.gen")


genind2genepop = function(data, file=""){
  
  ## Check input file a genind object
  if(!"genind" %in% class(data)){
    warning("Function was designed for genind objects.")
  }
  
  ## Check adegenet and miscTools are installed
  if(!require(adegenet)){install.packages("adegenet")}
  if(!require(miscTools)){install.packages("miscTools")}
  
  
  # ---------------- #
  #
  # Preamble
  #
  # ---------------- #
  
  ## Convert genind to dataframe object
  df = genind2df(data, usepop=FALSE)
  
  ## Convert A-01, C-02, G-03, T-04
  mat = as.matrix(df)
  mat = apply(mat, FUN=gsub, MARGIN=2, pattern="A", replacement="01") 
  mat = apply(mat, FUN=gsub, MARGIN=2, pattern="C", replacement="02")
  mat = apply(mat, FUN=gsub, MARGIN=2, pattern="G", replacement="03")
  mat = apply(mat, FUN=gsub, MARGIN=2, pattern="T", replacement="04")
  
  ## Convert NAs to 0000
  mat[is.na(mat)] = "0000"
  
  ## Add a column containing individual names with a comma afterwards
  ind = paste(indNames(data),",", sep="")
  mat = cbind(ind, mat)
  
  
  # ---------------- #
  #
  # Insert a Pop row
  #
  # ---------------- #
  
  ## Pop label that will separate each population
  popline = c("Pop", rep("", ncol(mat)-1 ))
  
  ## Count the number of individuals in each population
  pop_counts = data.frame(Counts = summary(data$pop))
  
  ## Add a column totalling the cumulative sum 
  pop_counts$Sum = cumsum(pop_counts$Counts)
  
  ## Insert a Pop row between each population
  for (i in 1:nrow(pop_counts)){
    
    # i is the row number and increases by 1 after each interation to compensate 
    # for the extra row being inserted each run through the loop
    pop.row = rep(NA, nrow(pop_counts))
    pop.row[i] = pop_counts$Sum[i] + i
    mat = insertRow(mat, pop.row[i], popline) 
  }
  
  # Remove the last Pop row
  mat = mat[-nrow(mat), ] 
  
  
  # ---------------- #
  #
  # Construct the Genepop file
  #
  # ---------------- #
  
  ## Genepop header
  file_date = format(Sys.time(), "%Y%m%d@%H%M") # date and time
  header = paste("Genepop file format", file_date)
  
  ## List of loci names separated by commas
  loc_names = paste(locNames(data), collapse=",")
  
  # Insert title, locus and pop rows at the beginning
  mat = insertRow(mat, 1, c(header, rep("", ncol(mat)-1 )))
  mat = insertRow(mat, 2, c(loc_names, rep("", ncol(mat)-1 )))
  mat = insertRow(mat, 3, popline)
  
  # Export file
  write.table(mat, file=file, quote=FALSE, col.names=F, row.names=F)

}
