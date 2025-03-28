Temporal dynamics of alien flora accumulation
================
Banerjee AK, Bhowmick AR
July 27, 2023

- [Preparation](#preparation)
- [Variable 1: Growth form](#variable-1-growth-form)
- [Variable 2: Use](#variable-2-use)
- [Variable 3: Introduction pathway](#variable-3-introduction-pathway)
- [Variable 4: Habitat](#variable-4-habitat)
- [Variable 5: Native range](#variable-5-native-range)
- [Variable 6: Naturalized range](#variable-6-naturalized-range)

**Objective: To compare the minimum residence time (MRT) values between
the three invasion categories of alien species (invasive, naturalized
and introduced), with and without considering the categorical
variables**

# Preparation

``` r
#Load the required libraries
library(dplyr)
library(ggpubr)
library(car)
library(tidyverse)
library(rstatix)
library(datarium)
library(FSA)
library(Rmisc)
library(ARTool)

#Set the working directory
setwd("~/comparative_analysis")
```

# Variable 1: Growth form

``` r
data<-read.csv("data_gf.csv")
head(data)
data$group<-as.factor(data$Category)
data$group1<-as.factor(data$invasion_status)
tgc <- summarySE(data, measurevar="mrt", groupvars=c("group","group1"))
tgc
#Shapiro wilk test (normality check)
data %>%
  group_by(group) %>%
  shapiro_test(mrt)
data %>%
  group_by(group1) %>%
  shapiro_test(mrt)
#Levene test (homogeneity of variance check)
leveneTest(mrt ~ group, data = data)
leveneTest(mrt ~ group1, data = data)
#comparison
m = art(mrt ~ group * group1, data=data)
anova(m)
art.con(m, "group", adjust="holm") %>% 
  summary() %>%  # add significance stars to the output
  mutate(sig. = symnum(p.value, corr=FALSE, na=FALSE,
                       cutpoints = c(0, .001, .01, .05, .10, 1),
                       symbols = c("***", "**", "*", ".", " ")))
```

# Variable 2: Use

``` r
data<-read.csv("data_use.csv")
head(data)
data$group<-as.factor(data$Category)
data$group1<-as.factor(data$invasion_status)
tgc <- summarySE(data, measurevar="mrt", groupvars=c("group","group1"))
tgc
#Shapiro wilk test (normality check)
data %>%
  group_by(group) %>%
  shapiro_test(mrt)
data %>%
  group_by(group1) %>%
  shapiro_test(mrt)
#Levene test (homogeneity of variance check)
leveneTest(mrt ~ group, data = data)
leveneTest(mrt ~ group1, data = data)
#comparison
m = art(mrt ~ group * group1, data=data)
anova(m)
art.con(m, "group", adjust="holm") %>% 
  summary() %>%  # add significance stars to the output
  mutate(sig. = symnum(p.value, corr=FALSE, na=FALSE,
                       cutpoints = c(0, .001, .01, .05, .10, 1),
                       symbols = c("***", "**", "*", ".", " ")))
art.con(m, "group1", adjust="holm") %>% 
  summary() %>%  # add significance stars to the output
  mutate(sig. = symnum(p.value, corr=FALSE, na=FALSE,
                       cutpoints = c(0, .001, .01, .05, .10, 1),
                       symbols = c("***", "**", "*", ".", " ")))
```

# Variable 3: Introduction pathway

``` r
data<-read.csv("data_int.csv")
head(data)
data$group<-as.factor(data$Category)
data$group1<-as.factor(data$invasion_status)
tgc <- summarySE(data, measurevar="mrt", groupvars=c("group","group1"))
tgc
#Shapiro wilk test (normality check)
data %>%
  group_by(group) %>%
  shapiro_test(mrt)
data %>%
  group_by(group1) %>%
  shapiro_test(mrt)
#Levene test (homogeneity of variance check)
leveneTest(mrt ~ group, data = data)
leveneTest(mrt ~ group1, data = data)
#comparison
m <- aov(mrt ~ group * group1, data=data)
summary(m)
kruskal.test(mrt ~ group1, data = data)
dunnTest(mrt~group1,data=data)
```

# Variable 4: Habitat

``` r
data<-read.csv("data_hab.csv")
head(data)
data$group<-as.factor(data$Category)
data$group1<-as.factor(data$invasion_status)
tgc <- summarySE(data, measurevar="mrt", groupvars=c("group","group1"))
tgc
#Shapiro wilk test (normality check)
data %>%
  group_by(group) %>%
  shapiro_test(mrt)
data %>%
  group_by(group1) %>%
  shapiro_test(mrt)
#Levene test (homogeneity of variance check)
leveneTest(mrt ~ group, data = data)
leveneTest(mrt ~ group1, data = data)
#comparison
m = art(mrt ~ group * group1, data=data)
anova(m)
art.con(m, "group", adjust="holm") %>% 
  summary() %>%  # add significance stars to the output
  mutate(sig. = symnum(p.value, corr=FALSE, na=FALSE,
                       cutpoints = c(0, .001, .01, .05, .10, 1),
                       symbols = c("***", "**", "*", ".", " ")))
```

# Variable 5: Native range

``` r
data<-read.csv("data_native.csv")
head(data)
data$group<-as.factor(data$Category)
data$group1<-as.factor(data$invasion_status)
tgc <- summarySE(data, measurevar="mrt", groupvars=c("group","group1"))
tgc
#Shapiro wilk test (normality check)
data %>%
  group_by(group) %>%
  shapiro_test(mrt)
data %>%
  group_by(group1) %>%
  shapiro_test(mrt)
#Levene test (homogeneity of variance check)
leveneTest(mrt ~ group, data = data)
leveneTest(mrt ~ group1, data = data)
#comparison
m = art(mrt ~ group * group1, data=data)
anova(m)
art.con(m, "group", adjust="holm") %>% 
  summary() %>%  # add significance stars to the output
  mutate(sig. = symnum(p.value, corr=FALSE, na=FALSE,
                       cutpoints = c(0, .001, .01, .05, .10, 1),
                       symbols = c("***", "**", "*", ".", " ")))
```

# Variable 6: Naturalized range

``` r
data<-read.csv("data_naturalized.csv")
head(data)
data$group<-as.factor(data$Category)
data$group1<-as.factor(data$invasion_status)
tgc <- summarySE(data, measurevar="mrt", groupvars=c("group","group1"))
tgc
#Shapiro wilk test (normality check)
data %>%
  group_by(group) %>%
  shapiro_test(mrt)
data %>%
  group_by(group1) %>%
  shapiro_test(mrt)
#Levene test (homogeneity of variance check)
leveneTest(mrt ~ group, data = data)
leveneTest(mrt ~ group1, data = data)
#comparison
m = art(mrt ~ group * group1, data=data)
anova(m)
art.con(m, "group", adjust="holm") %>% 
  summary() %>%  # add significance stars to the output
  mutate(sig. = symnum(p.value, corr=FALSE, na=FALSE,
                       cutpoints = c(0, .001, .01, .05, .10, 1),
                       symbols = c("***", "**", "*", ".", " ")))
art.con(m, "group1", adjust="holm") %>% 
  summary() %>%  # add significance stars to the output
  mutate(sig. = symnum(p.value, corr=FALSE, na=FALSE,
                       cutpoints = c(0, .001, .01, .05, .10, 1),
                       symbols = c("***", "**", "*", ".", " ")))
```

**END**
