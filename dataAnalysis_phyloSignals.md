Temporal dynamics of alien flora accumulation
================
Banerjee AK, Bhowmick AR
July 27, 2023

- [Preparation](#preparation)
- [Making the phylogenetic tree](#making-the-phylogenetic-tree)
- [Phylogenetic signal](#phylogenetic-signal)
  - [Biotic variables](#biotic-variables)
  - [A biotic variables](#a-biotic-variables)

**Objective: To ascertain presence of phylogenetic signals in the
variables**

# Preparation

``` r
#Load the required libraries
devtools::install_github ("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)
library(phytools)
library(picante)
library(ade4)
library(ape)
library(nlme)
library(geiger)
library(phangorn)
library(phylobase)

#Set the working directory
setwd("~/regression_analysis")
```

# Making the phylogenetic tree

``` r
data<-read.csv("species_list_phylo.csv")
tree.a<-phylo.maker(sp.list = data,tree = GBOTB.extended,nodes = nodes.info.1,
                    scenarios = "S3")
write.tree(tree.a$scenario.3,"tree1.tre")
```

# Phylogenetic signal

## Biotic variables

``` r
#GROWTH FORM#
tree <-read.tree("tree2.tre") #species names without "_"
trait<-read.csv("species_data_gf.csv")
trait<-as.data.frame(trait)
tmpTr <- drop.tip(tree, tree$tip.label[! tree$tip.label %in% trait[, 1]])
tr4 <- phylo4d(tmpTr, trait, label.type="column")
pruned <- as(extractTree(tr4), "phylo")
ldata <- tdata(tr4, "tip")
ldata <- ldata[pruned$tip.label, 1]
names(ldata) <- pruned$tip.label
#run the function
phylo.signal.disc(ldata, pruned) ->result
write(result$.Randomization.Results, file="gf.out")

#HEIGHT#
tree<-read.tree("tree2.tre") #species names without "_"
comm<-read.table("comm_height.txt") #species names arranged in columns
tree1<-prune.sample(comm, tree) #to match tree tips with species names in comm
comm<-comm[,tree1$tip.label] # to order species in tree and comm
comm #visualize
write.csv(comm,"comm_2.csv") #write the ordered species
trait<-read.csv("species_data_height.csv") #read trait data (ordered in excel)
trait<-as.data.frame(trait) #making the data frame
k.trait<-phylosig(tree1,trait$height_scaled,test = TRUE)
print(k.trait)
plot(k.trait)
lambda.trait<-phylosig(tree1,trait$height_scaled,method = "lambda",test = TRUE)
print(lambda.trait)
plot(lambda.trait)

#SEED#
tree<-read.tree("tree2.tre")
comm<-read.table("comm_seed.txt")
tree1<-prune.sample(comm, tree)
comm<-comm[,tree1$tip.label]
comm
write.csv(comm,"comm_2.csv")
trait<-read.csv("species_data_seed.csv")
trait<-as.data.frame(trait)
k.trait<-phylosig(tree1,trait$seed_scaled,test = TRUE)
print(k.trait)
plot(k.trait)
lambda.trait<-phylosig(tree1,trait$seed_scaled,method = "lambda",test = TRUE)
print(lambda.trait)
plot(lambda.trait)

#SLA#
tree<-read.tree("tree2.tre")
comm<-read.table("comm_sla.txt")
tree1<-prune.sample(comm, tree)
comm<-comm[,tree1$tip.label]
comm
write.csv(comm,"comm_2.csv")
trait<-read.csv("species_data_sla.csv")
trait<-as.data.frame(trait)
k.trait<-phylosig(tree1,trait$sla,test = TRUE)
print(k.trait)
plot(k.trait)
lambda.trait<-phylosig(tree1,trait$sla,method = "lambda",test = TRUE)
print(lambda.trait)
plot(lambda.trait)

#NATIVE CONGENERS#
tree<-read.tree("tree2.tre") 
comm<-read.table("comm_natcon.txt") 
tree1<-prune.sample(comm, tree)
comm<-comm[,tree1$tip.label]
comm #visualize
write.csv(comm,"comm_2.csv")
trait<-read.csv("species_data_natcon.csv")
trait<-as.data.frame(trait)
k.trait<-phylosig(tree1,trait$Nat.con,test = TRUE)
print(k.trait)
plot(k.trait)
lambda.trait<-phylosig(tree1,trait$Nat.con,method = "lambda",test = TRUE) 
print(lambda.trait)
plot(lambda.trait)
```

## A biotic variables

``` r
#USE#
tree<-read.tree("tree2.tre") #species names without "_"
comm<-read.table("comm_use.txt") #species names arranged in columns
tree1<-prune.sample(comm, tree) #to match tree tips with species names in comm
comm<-comm[,tree1$tip.label] # to order species in tree and comm
comm #visualize
write.csv(comm,"comm_2.csv") #write the ordered species
trait<-read.csv("species_data_use.csv") #read trait data (ordered in excel)
trait<-as.data.frame(trait) #making the data frame
k.trait<-phylosig(tree1,trait$uses_nu,test = TRUE)
print(k.trait)
plot(k.trait)
lambda.trait<-phylosig(tree1,trait$uses_nu,method = "lambda",test = TRUE)
print(lambda.trait)
plot(lambda.trait)

#INTRODUCTION PATHWAY#
tree<-read.tree("tree2.tre") 
comm<-read.table("comm_int.txt") 
tree1<-prune.sample(comm, tree) 
comm<-comm[,tree1$tip.label] 
comm #visualize
write.csv(comm,"comm_2.csv") 
trait<-read.csv("species_data_int.csv") 
trait<-as.data.frame(trait) 
k.trait<-phylosig(tree1,trait$int_nu,test = TRUE)
print(k.trait)
plot(k.trait)
lambda.trait<-phylosig(tree1,trait$int_nu,method = "lambda",test = TRUE) 
print(lambda.trait)
plot(lambda.trait)

#HABITAT#
tree<-read.tree("tree2.tre") 
comm<-read.table("comm_habitat.txt") 
tree1<-prune.sample(comm, tree) 
comm<-comm[,tree1$tip.label] 
comm #visualize
write.csv(comm,"comm_2.csv") 
trait<-read.csv("species_data_habitat.csv") 
trait<-as.data.frame(trait)
k.trait<-phylosig(tree1,trait$habitat,test = TRUE) 
print(k.trait)
plot(k.trait)
lambda.trait<-phylosig(tree1,trait$habitat,method = "lambda",test = TRUE) 
print(lambda.trait)
plot(lambda.trait)

#NATIVE RANGE#
tree<-read.tree("tree2.tre") 
comm<-read.table("comm_native.txt") 
tree1<-prune.sample(comm, tree)
comm<-comm[,tree1$tip.label]
comm #visualize
write.csv(comm,"comm_2.csv")
trait<-read.csv("species_data_native.csv") 
trait<-as.data.frame(trait) 
k.trait<-phylosig(tree1,trait$natRng,test = TRUE) 
print(k.trait)
plot(k.trait)
lambda.trait<-phylosig(tree1,trait$natRng,method = "lambda",test = TRUE) 
print(lambda.trait)
plot(lambda.trait)

#NATURALIZED RANGE#
tree<-read.tree("tree2.tre")
comm<-read.table("comm_naturalized.txt")
tree1<-prune.sample(comm, tree) 
comm<-comm[,tree1$tip.label] 
comm #visualize
write.csv(comm,"comm_2.csv") 
trait<-read.csv("species_data_naturalized.csv")
trait<-as.data.frame(trait) 
k.trait<-phylosig(tree1,trait$natuRng,test = TRUE) 
print(k.trait)
plot(k.trait)
lambda.trait<-phylosig(tree1,trait$natuRng,method = "lambda",test = TRUE) 
print(lambda.trait)
plot(lambda.trait)

########################################################
####### The above code chunk is fitted for the three invasion categories.
####### The necessary data files are in the same folder.
####### The text files are named as: 
## "comm_<variable name>_in (for invasive)
## "comm_<variable name>_nt (for naturalized)
## "comm_<variable name>_int (for introduced)
####### The CSV files are named as:
## "species_data_<variable name>_in (for invasive)
## "species_data_<variable name>_nt (for naturalized)
## "species_data_<variable name>_int (for introduced)
########################################################
########################################################
```

**END**
