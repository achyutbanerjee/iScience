Temporal dynamics of alien flora accumulation
================
Banerjee AK, Bhowmick AR
July 27, 2023

- [Preparation](#preparation)
- [Phylogenetic path analysis](#phylogenetic-path-analysis)

**Objective: To decipher the direct and indirect relationships among
variables and their influences on minimum residence time (MRT) after
considering the phylogenetic relationship between the species**

# Preparation

``` r
#Load the required libraries
library(phylopath)
library(picante)
library(ggraph)

#Set the working directory
setwd("~/regression_analysis")
```

# Phylogenetic path analysis

``` r
tree<-read.tree("tree2.tre") #species names without "-"(e.g.Opuntia ficus-indica)
comm<-read.table("comm_species.txt") #species names arranged in columns
tree1<-prune.sample(comm, tree) #to match tree tips with species names in comm

#Read data
data<-read.csv("data1.csv") #imputed data
rownames(data) <- data$Species
head(data)

#Define models
models <- define_model_set(
  one   = c(use ~ nat+natu_l+nc+gf+hab),
  two = c(int ~ use+nat+natu_l+nc+gf+hab),
  three  = c(nat ~ gf+hab+use+seed_s+height_s),
  four  = c(natu_l ~ gf+hab+nat+use+seed_s+height_s),
  five   = c(hab ~ gf+seed_s+height_s),
  .common = c(mrt ~ gf+nc+use+int+nat+natu_l+hab+seed_s+height_s)
)
result <- phylo_path(models, data = data, tree = tree1, model = 'lambda',order = NULL)
result
(s <- summary(result))

#To view the best ranked model
#This returns a DAG with standardized regression coefficients, as well as a matrix of standard errors.
(best_model <- best(result))
plot(best_model)

#confidence in the regression coefficients
#coef_plot can visualize the estimates and their approximate confidence intervals
#need to use the boot parameter when fitting the DAG if you are not model averaging
best_model1 <- best(result,boot=500)
coef_plot(best_model1)

#Path averaging is only done for the models that actually contain that path. 
#This facilitates the detection of weak effects, but also biases coefficients away from zero.
average_model <- average(result)
average_model
plot(average_model, algorithm = 'mds', curvature = 0.5)
coef_plot(average_model)

#Alternatively, we can assume the coefficients (and their variance) for absent paths to be zero by setting
#avg_method=full
average_model_full <- average(result, avg_method = "full")
average_model_full
plot(average_model_full, algorithm = 'mds', curvature = 0.5)
coef_plot(average_model_full)

########################################################
####### The above code chunk is fitted for the three invasion categories.
####### The necessary data files are in the same folder.
####### The imputed data files (.csv) are named as: 
## "data1_in (for invasive)
## "data1_nt (for naturalized)
## "data1_int (for introduced)
####### The species files (.txt) are named as:
## "comm_species_in (for invasive)
## "comm_species_nt (for naturalized)
## "comm_species_int (for introduced)
########################################################
########################################################
```

**END**
