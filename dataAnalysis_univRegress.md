Temporal dynamics of alien flora accumulation
================
Banerjee AK, Bhowmick AR
July 27, 2023

- [Preparation](#preparation)
- [Analysis for variables with phylogenetic signals: Phylogenetic
  Generalized Least Square
  (PGLS)](#analysis-for-variables-with-phylogenetic-signals-phylogenetic-generalized-least-square-pgls)
  - [Biotic variables](#biotic-variables)
  - [Abiotic variables](#abiotic-variables)
  - [Plot](#plot)
- [Analysis for variables without phylogenetic
  signals](#analysis-for-variables-without-phylogenetic-signals)
  - [Plot](#plot-1)

**Objective: To identify direct influence of the biotic and abiotic
variables on minimum residence time (MRT) after considering the
phylogenetic relationship between the species**

# Preparation

``` r
#Load the required libraries
devtools::install_github ("jinyizju/V.PhyloMaker")
library(V.PhyloMaker)
library(phytools)
library(picante)
library(ade4) # source of example data and tree
library(ape) # tree handling
library(nlme) # regression modelling
library(geiger)
library(phangorn)
library(phylobase)
library(tidyverse)

#Set the working directory
setwd("~/regression_analysis")
```

# Analysis for variables with phylogenetic signals: Phylogenetic Generalized Least Square (PGLS)

## Biotic variables

``` r
#GROWTH FORM#
tree<-read.tree("tree2.tre") #species names without "-"
comm<-read.table("comm_gf.txt") #species names arranged in columns
tree1<-prune.sample(comm, tree) #to match tree tips with species names in comm
comm<-comm[,tree1$tip.label] #to order species in tree and comm
data<-read.table("species_data_gf.csv",header=TRUE, sep=",") #species names arranged in columns
data<-as.data.frame(data)
head(data)
data_frame <- data %>% column_to_rownames(var="species")
#Using corSymm correlation structure
mat <- vcv(tree1, corr=TRUE) #construct matrix
fit <- gls(mrt ~ gf, correlation=corSymm(mat[lower.tri(mat)],fixed=TRUE), data=data_frame)
summary(fit)
#Using corBrownian correlation structure
fit2 <- gls(mrt ~ gf, correlation=corBrownian(phy=tree1), data=data_frame)
summary(fit2)
anova(fit, fit2) #models are the same

#HEIGHT#
tree<-read.tree("tree2.tre") 
comm<-read.table("comm_height.txt") 
tree1<-prune.sample(comm, tree) 
comm<-comm[,tree1$tip.label] 
data<-read.table("species_data_height.csv",header=TRUE, sep=",") 
data<-as.data.frame(data)
head(data)
data_frame <- data %>% column_to_rownames(var="species")
#Using corSymm correlation structure
mat <- vcv(tree1, corr=TRUE) 
fit <- gls(mrt ~ height_scaled, correlation=corSymm(mat[lower.tri(mat)],fixed=TRUE), data=data_frame)
summary(fit)
#Using corBrownian correlation structure
fit2 <- gls(mrt ~ height_scaled, correlation=corBrownian(phy=tree1), data=data_frame)
summary(fit2)
anova(fit, fit2) 

#SEED#
tree<-read.tree("tree2.tre") 
comm<-read.table("comm_seed.txt") 
tree1<-prune.sample(comm, tree) 
comm<-comm[,tree1$tip.label] 
data<-read.table("species_data_seed.csv",header=TRUE, sep=",") 
data<-as.data.frame(data)
head(data)
data_frame <- data %>% column_to_rownames(var="species")
#Using corSymm correlation structure
mat <- vcv(tree1, corr=TRUE) # construct matrix
fit <- gls(mrt ~ seed_scaled, correlation=corSymm(mat[lower.tri(mat)],fixed=TRUE), data=data_frame)
summary(fit)
#Using corBrownian correlation structure
fit2 <- gls(mrt ~ seed_scaled, correlation=corBrownian(phy=tree1), data=data_frame)
summary(fit2)
anova(fit, fit2) 

#SLA#
tree<-read.tree("tree2.tre") 
comm<-read.table("comm_sla.txt") 
tree1<-prune.sample(comm, tree) 
comm<-comm[,tree1$tip.label] 
data<-read.table("species_data_sla.csv",header=TRUE, sep=",") 
data<-as.data.frame(data)
head(data)
data_frame <- data %>% column_to_rownames(var="species")
#Using corSymm correlation structure
mat <- vcv(tree1, corr=TRUE) 
fit <- gls(mrt ~ sla, correlation=corSymm(mat[lower.tri(mat)],fixed=TRUE), data=data_frame)
summary(fit)
#Using corBrownian correlation structure
fit2 <- gls(mrt ~ sla, correlation=corBrownian(phy=tree1), data=data_frame)
summary(fit2)
anova(fit, fit2) 

#NATIVE CONGENER#
tree<-read.tree("tree2.tre") 
comm<-read.table("comm_natcon.txt") 
tree1<-prune.sample(comm, tree) 
comm<-comm[,tree1$tip.label] 
data<-read.table("species_data_natcon.csv",header=TRUE, sep=",") 
data<-as.data.frame(data)
head(data)
data_frame <- data %>% column_to_rownames(var="species")
#Using corSymm correlation structure
mat <- vcv(tree1, corr=TRUE) 
fit <- gls(mrt ~ Nat.con, correlation=corSymm(mat[lower.tri(mat)],fixed=TRUE), data=data_frame)
summary(fit)
#Using corBrownian correlation structure
fit2 <- gls(mrt ~ Nat.con, correlation=corBrownian(phy=tree1), data=data_frame)
summary(fit2)
anova(fit, fit2)
```

## Abiotic variables

``` r
#USE#
tree<-read.tree("tree2.tre") #species names without "-"
comm<-read.table("comm_use.txt") #species names arranged in columns
tree1<-prune.sample(comm, tree) #to match tree tips with species names in comm
comm<-comm[,tree1$tip.label] # to order species in tree and comm
data<-read.table("species_data_use.csv",header=TRUE, sep=",") #species names arranged in columns
data<-as.data.frame(data)
head(data)
data_frame <- data %>% column_to_rownames(var="species")
#Using corSymm correlation structure
mat <- vcv(tree1, corr=TRUE) # construct matrix
fit <- gls(mrt ~ uses_nu, correlation=corSymm(mat[lower.tri(mat)],fixed=TRUE), data=data_frame)
summary(fit)
#Using corBrownian correlation structure
fit2 <- gls(mrt ~ uses_nu, correlation=corBrownian(phy=tree1), data=data_frame)
summary(fit2)
anova(fit, fit2) # models are the same

#INTRODUCTION PATHWAY#
tree<-read.tree("tree2.tre") #species names without "-"
comm<-read.table("comm_int.txt") #species names arranged in columns
tree1<-prune.sample(comm, tree) #to match tree tips with species names in comm
comm<-comm[,tree1$tip.label] # to order species in tree and comm
data<-read.table("species_data_int.csv",header=TRUE, sep=",") #species names arranged in columns
data<-as.data.frame(data)
head(data)
data_frame <- data %>% column_to_rownames(var="species")
logistic_model <- glm(mrt~int_nu, family = gaussian, data, maxit = 100)
summary(logistic_model)

#HABITAT#
tree<-read.tree("tree2.tre") 
comm<-read.table("comm_habitat.txt") 
tree1<-prune.sample(comm, tree) 
comm<-comm[,tree1$tip.label] 
data<-read.table("species_data_habitat.csv",header=TRUE, sep=",") 
data<-as.data.frame(data)
head(data)
data_frame <- data %>% column_to_rownames(var="species")
#Using corSymm correlation structure
mat <- vcv(tree1, corr=TRUE) 
fit <- gls(mrt ~ habitat, correlation=corSymm(mat[lower.tri(mat)],fixed=TRUE), data=data_frame)
summary(fit)
#Using corBrownian correlation structure
fit2 <- gls(mrt ~ habitat, correlation=corBrownian(phy=tree1), data=data_frame)
summary(fit2)
anova(fit, fit2) 

#NATIVE RANGE#
tree<-read.tree("tree2.tre") 
comm<-read.table("comm_native.txt") 
tree1<-prune.sample(comm, tree) 
comm<-comm[,tree1$tip.label] 
data<-read.table("species_data_native.csv",header=TRUE, sep=",") 
data<-as.data.frame(data)
head(data)
data_frame <- data %>% column_to_rownames(var="species")
#Using corSymm correlation structure
mat <- vcv(tree1, corr=TRUE) 
fit <- gls(mrt ~ natRng, correlation=corSymm(mat[lower.tri(mat)],fixed=TRUE), data=data_frame)
summary(fit)
#Using corBrownian correlation structure
fit2 <- gls(mrt ~ natRng, correlation=corBrownian(phy=tree1), data=data_frame)
summary(fit2)
anova(fit, fit2) 

#NATURALIZED RANGE#
tree<-read.tree("tree2.tre")
comm<-read.table("comm_naturalized.txt") 
tree1<-prune.sample(comm, tree) 
comm<-comm[,tree1$tip.label]
data<-read.table("species_data_naturalized.csv",header=TRUE, sep=",")
data<-as.data.frame(data)
head(data)
data_frame <- data %>% column_to_rownames(var="species")
#Using corSymm correlation structure
mat <- vcv(tree1, corr=TRUE) 
fit <- gls(mrt ~ natuRng, correlation=corSymm(mat[lower.tri(mat)],fixed=TRUE), data=data_frame)
summary(fit)
#Using corBrownian correlation structure
fit2 <- gls(mrt ~ natuRng, correlation=corBrownian(phy=tree1), data=data_frame)
summary(fit2)
anova(fit, fit2) 
```

## Plot

``` r
#For one variable: NATURALIZED RANGE

library(ggplot2)
library(ggsci)
library(ggeffects)
mydf <- ggpredict(fit2, terms = "naturalized")

p<-ggplot(mydf, aes(x, predicted)) +
  geom_line(color="salmon",linewidth=1.0) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1,fill="red")+
  geom_point()+
  scale_color_jco()+
  scale_fill_jco()
p
p + labs(x = "SLA", y = "MRT") + 
  theme(axis.title = element_text(size = 16,color = "black",face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 14,color = "black"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
  )
```

# Analysis for variables without phylogenetic signals

``` r
########################################################
####### All variables have phylogenetic signals when three invasion categories
####### were considered together (see Table S7).
####### For individual invasion categories, some variables do not have
####### phylogenetic signals. Logistic regression models are fitted for
####### these variables. One example (SLA for invasive alien category) 
####### is given below.
####### The above code chunk is fitted for the three invasion categories
####### separately for those variables having phylogenetic signals.
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

#SLA of invasive aliens#
tree<-read.tree("tree2.tre") #species names without "-"
comm<-read.table("comm_sla_in.txt") #species names arranged in columns
tree1<-prune.sample(comm, tree) #to match tree tips with species names in comm
comm<-comm[,tree1$tip.label] # to order species in tree and comm
data<-read.table("species_data_sla_in.csv",header=TRUE, sep=",") #species names arranged in columns
data<-as.data.frame(data)
head(data)
data_frame <- data %>% column_to_rownames(var="species")
logistic_model <- glm(mrt~sla, family = gaussian, data, maxit = 100)
summary(logistic_model)
```

## Plot

``` r
library(ggplot2)
library(ggsci)
library(ggeffects)
mydf <- ggpredict(logistic_model, terms = "sla")

p<-ggplot(mydf, aes(x, predicted)) +
  geom_line(color="salmon",linewidth=1.0) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .1,fill="red")+
  geom_point()+
  scale_color_jco()+
  scale_fill_jco()
p
p + labs(x = "SLA", y = "MRT") + 
  theme(axis.title = element_text(size = 16,color = "black",face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 14,color = "black"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
  )
```

**END**
