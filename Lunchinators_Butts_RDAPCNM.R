rm(list=ls())


# RDA-PCNM Lunchinators

# Example of how PCNM works - From EEB698-Species Composition Analysis ####

# If not installed uncomment this line - install.packages('ade4')
library(ade4)
library(vegan)

# Thanks Dr. Dixon for this example code! 

# generate a 10 x 6 regular grid and compute pcnm scores

xy <- expand.grid(x=1:10, y = 1:6)
head(xy)
tail(xy)
xy

xy.pcnm <- pcnm(dist(xy)) #Default distance matrix = Eucl.
dim(xy.pcnm$vectors)

# visualize various PCNM vectors using image plots
par(mfrow=c(3,2), mar=c(3,3,0,0)+0.2, mgp=c(2,0.8,0))

axis <- 1;
image(1:10, 1:6, matrix(xy.pcnm$vectors[,axis], nrow=10))
axis <- 2;
image(1:10, 1:6, matrix(xy.pcnm$vectors[,axis], nrow=10))
axis <- 3;
image(1:10, 1:6, matrix(xy.pcnm$vectors[,axis], nrow=10))
axis <- 4;
image(1:10, 1:6, matrix(xy.pcnm$vectors[,axis], nrow=10))
axis <- 5;
image(1:10, 1:6, matrix(xy.pcnm$vectors[,axis], nrow=10))
axis <- 6;
image(1:10, 1:6, matrix(xy.pcnm$vectors[,axis], nrow=10))

axis <- 10;
image(1:10, 1:6, matrix(xy.pcnm$vectors[,axis], nrow=10))
axis <- 20;
image(1:10, 1:6, matrix(xy.pcnm$vectors[,axis], nrow=10))
axis <- 30;
image(1:10, 1:6, matrix(xy.pcnm$vectors[,axis], nrow=10))
axis <- 35;
image(1:10, 1:6, matrix(xy.pcnm$vectors[,axis], nrow=10))
axis <- 40;
image(1:10, 1:6, matrix(xy.pcnm$vectors[,axis], nrow=10))
axis <- 42;
image(1:10, 1:6, matrix(xy.pcnm$vectors[,axis], nrow=10))


# Further extension - Variation partitioning ####
  ## Paritition variation in fungal communities associated with tree 
  ## leaves to three categories: leaf type, leaf chemistry, 
  ## and spatial relationships on the landscape 

# read in tables containing species, environmental variables, 
  # and geographic distances

# Set working directory to wherever you downloaded these files! 

endo.spp <- read.csv('endophytes.csv')  # column names represent OTUs
endo.env <- read.csv('endophytes_env.csv') # Environmental data 
endo.dist <- read.csv('endophytes_dist.csv') # Spatial Data 
par(mfrow=c(1,1))

dim(endo.spp)
str(endo.env)
str(endo.dist)

# look at the spatial arrangment of samples in two dimentions
plot(endo.dist)

# represent spatial patterns through PCNMs
endo.pcnm <- pcnm(dist(endo.dist))

# loadings for each PCNM axis can be extracted using scores()
str(scores(endo.pcnm))

# Two-dimensional sampling, especially when samples are not evenly spaced,
  ## results in variables with patterns that are not as easy to interpret,
  ## but the first axes tend to represent large scale patterns 
  ## while the later axes tend to represent smaller scale patterns.
  ## Here we plot the first six PCNM axes to visualise the patterns 
  ## that they represent.
# ordisurf() fits a smooth surface model to the data, 
  ## but the output is suppressed here
# as a biproduct, ordisurf returns the plot that we want
par(mfrow=c(3, 2))
ordisurf(endo.dist,scores(endo.pcnm,choi=1),bubble=4,
         main='PCNM 1')
ordisurf(endo.dist,scores(endo.pcnm,choi=2),bubble=4,
         main='PCNM 2')
ordisurf(endo.dist,scores(endo.pcnm,choi=3),bubble=4,
         main='PCNM 3')
ordisurf(endo.dist,scores(endo.pcnm,choi=4),bubble=4,
         main='PCNM 4')
ordisurf(endo.dist,scores(endo.pcnm,choi=5),bubble=4,
         main='PCNM 5')
ordisurf(endo.dist,scores(endo.pcnm,choi=6),bubble=4,
         main='PCNM 6')

# Now there are two matrices 
  ## one representing variation in measured environmental variables 
  ## and other representing spatial distributions of samples 
# We can now estimate the amount of variation in community comp. that each 
  ## explains 


# do predictor matrices affect community composition?

# environmental matrix as predictor
# Capsscale function is a distance-based redundancy analysis 
  ## Same as a normal RDA, but allows different dissimilarity matrices
cap.env <- capscale(endo.spp ~ ., data=endo.env, 
                    distance='bray')
cap.env

# spatial matrix as predictor
cap.pcnm <- capscale(endo.spp ~ ., 
                     data=as.data.frame(scores(endo.pcnm)), 
                     distance='bray')
cap.pcnm

# select particular variables to proceed with 
  ## (here we use both forward and backward selection 
  ## but could use either one separately)

# set up the null cases with no predictors
mod0.env <- capscale(endo.spp ~ 1, data=endo.env, distance='bray')
mod0.pcnm <- capscale(endo.spp ~ 1, data=as.data.frame(scores(endo.pcnm)), distance='bray')

# select variables in each predictor table
step.env <- ordistep(mod0.env, scope=formula(cap.env))

step.pcnm <- ordistep(mod0.pcnm, scope=formula(cap.pcnm))

# Now can we can look at the result for each explanatory matrix including
  ## only those variables that best explain variation between communities 
  ## can also assess significance level 

# species, tissue type, tissue CN ratio and N concentration 
  ## predict variation in community composition
step.env
step.env$anova  # presents results in an ANOVA-like table

# only six/seven of the PCNM axes appear to predict variation 
  ## in community composition
  ## significance of PCNM11 varies each time because it is based on permutations
step.pcnm
step.pcnm$anova  # presents results in an ANOVA-like table

# Now we can perform variation partitioning by including both types of 
  ## predictor matrices as separate arguments 
# This is done using a species-sample matrix as the response matrix 
  ## Basically a redundancy analysis 

# create pcnm table with only significant axes
endo.pcnm.sub <- scores(endo.pcnm, 
                        choices=c(1:4, 6, 11, 14))

# partition variation among four predictor tables:
#   1) leaf species 
#   2) leaf type (canopy/litter)
#   3) leaf chemistry
#   4) spatial gradients

## varpart function partitions variation in community data or community 
  ## dissimilarity with respect to 2, 3, or 4 explanatory tables using 
  ## adjusted R-squared in RDA or db-RDA 
endo.var <- varpart(endo.spp, 
                    ~ species, 
                    ~ type, 
                    ~ CNratio + percentN, 
                    endo.pcnm.sub, data=endo.env)
#endo.var
par(mfrow=c(1, 2))
showvarparts(4)
plot(endo.var)

#The plot shows the adjusted R2 values associated with each partition 
  ## or for overlapping partitions. Although the plot does not compare labels, 
  ## we can print the object containing the result of the call to 
  ## and get the relevant information from the 'Explanatory tables' and 
  ## 'Partition table' (under 'Individual fractions') blocks.

# Then we can test the significance of each individual component. 
  ## Here we use 'rda' function since that is how we performed variation partitioning.
  ## We test for the effect of each individual component after already accounting 
  ## for the effect of the other components using the function 
  ## (note that we need to specify the name of the 
  ## dataframe for variables named within ).  

# significance of partition X1
anova(rda(endo.spp  ~ species + Condition(endo.env$type) +  
            + Condition(endo.env$CNratio) + Condition(endo.env$percentN) 
          + Condition(endo.pcnm.sub), 
          data=endo.env))

# significance of partition X2
anova(rda(endo.spp  ~ type + Condition(endo.env$species) +  
            + Condition(endo.env$CNratio) + Condition(endo.env$percentN) 
          + Condition(endo.pcnm.sub), 
          data=endo.env))

# significance of partition X3
anova(rda(endo.spp  ~ CNratio + percentN 
          + Condition(endo.env$species) + Condition(endo.env$type) 
          + Condition(endo.pcnm.sub), 
          data=endo.env))

# significance of partition X4
anova(rda(endo.spp  ~ endo.pcnm.sub 
          + Condition(endo.env$species) + Condition(endo.env$type) 
          + Condition(endo.env$CNratio) + Condition(endo.env$percentN)))

# Leaf species, type (litter/fresh), and chemistry are all statistically significant 
  ## in their contributions to fungal community composition,
  ## even though the amount of variation explained by some of these variables is small.
  ## There was no statistically significant spatial pattern in fungal community composition 
  ## independent of these other measured drivers.

