#######################################################################
#
# GAM_Intro.R
#
# Me following along with Wieling (2018) on using GAMs for linguistics
#
# UC Santa Cruz * M. Brinkerhoff * 2022-11-28 (M)
#
#######################################################################

# install packages if not yet installed
packages <- c("mgcv","itsadug","lme4","sp")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))  
}

# load required packages
library(mgcv)
library(itsadug)
library(sp) # for colors which also print well in grayscale
library(lme4)

# version information
R.version.string

# loading dataset
if (!file.exists("full.rda")) { 
  download.file("http://www.let.rug.nl/wieling/Tutorial/full.rda",
                "full.rda") # 2 MB
}

load("full.rda")

# investigating the dataset
str(full)
head(full)

# data wrangling 
words = c("tent","tenth")
dat = droplevels(full[full$Word %in% words,]) # removes unwanted lines
dat$Word = relevel(dat$Word,"tent") # tent is set as the reference level

# track first observation per trial (to correct for autocorrelation later on)
# and sort per trajectory
dat$start.event = (dat$Time == 0) 
dat = dat[order(dat$Speaker,dat$Trial,dat$Time),]
# dat <- start_event(dat,event=c("Speaker","Trial")) # bug in itsadug 2.4

# fitting linear models
# bam and gam do essentially the same thing. 
# bam is used for very large datasets
m1 <- bam(Pos ~ Word, # formula
          data=dat, # dataset
          method="fREML" #method used for GAM
          )
(smry1 <- summary(m1)) # print the variable you define

# fitting non-linear models
m2 <- bam(Pos ~ Word + # formula looking at position of T1 sensor by words
            s(Time, # what you are smoothing over
              by=Word, # for which levels 
              bs="tp", # type of smooth (tp = thin plate regression spline)
              k=10 # size of basis function
              ), # smooth over time 
          data=dat # where the data comes from
          )
(smry2 <- summary(m2)) 

# checking the fitting procedure and results
gam.check(m2)

# The first lines show that the model converged on a solution. 
# The bottom lines are associated with the smooths. It shows 
# the edf values together with k' (i.e. k - 1). If the value of 
# k-index is lower than 1 and the associated p-value is low, this 
# suggests that the basis dimension has been restricted too much. 
# In that case, it is good practice to refit the model with the 
# value of k doubled. In this case, there is no reason to do so, 
# as the value of k-index is not smaller than 1 and the p-value is 
# relatively high. 
#
# In principle, the k-parameter can be set as high as the number of 
# unique values in the data, as penalization will result in the 
# appropriate shape. However, allowing for more complexity negatively 
# impacts computation time.

# visualizing the results
par(mfrow=c(1,2)) # setting the parameters for the graphic
plot(m2, # mgcv function for ploting GAMs. 
     ylim=c(-1,2), # the y limits for each plot
     select=1 # which part of the parameter it gets printed in
     )
abline(h=0) # adds a straight line to the plot at zero
plot(m2, ylim=c(-1,2), select=2)
abline(h=0)

# This is ok but only visualizes the two non-linear patterns without 
# taking into account anything else in the model. This means that 
# only the partial effects are visualized. 

# get around this limitation by using plot_smooth
two.colors = bpy.colors(n = 2, 
                        cutoff.tails = 0.7, 
                        alpha = 1.0
                        ) # spec colors that look great in greyscale 
plot_smooth(m2, 
            view="Time", 
            plot_all="Word",
            main="",
            rug=FALSE, 
            ylim=c(-1,2), 
            col=two.colors
            ) 

# plotting the difference between the 
plot_diff(m2, 
          view="Time", 
          comp=list(Word=c("tenth","tent")),
          ylim=c(-1,2)
          )

# do we need more complexity in our model
# Even tho we see that there is a clear difference it is 
# important to see if this is true statistically. 
# There are three ways that we can investigate this:
# 1. model comparison
# 2. Refitting the model with a binary difference smooth
# 3. Refitting the model with an ordered factor difference smooth

# model comparison
# a model comparison is only possible when using 
m2a.ml <- bam(Pos ~ Word + # formula
              s(Time), # smoothing over time
              data=dat, # dataset
              method="ML" # maximum likelihood (ML) estimation method
              )
m2b.ml <- bam(Pos ~ Word + s(Time, by=Word), data=dat, method="ML")
compareML(m2a.ml, m2b.ml)