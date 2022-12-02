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
            main="m2",
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

# binary difference smooth
dat$IsTenth <- (dat$Word == "tenth")*1

m2.bin <- bam(Pos ~ s(Time) + s(Time,by=IsTenth), data=dat)
(smry2bin <- summary(m2.bin))

# ordered factor difference smooth
dat$WordO <- as.ordered(dat$Word)
contrasts(dat$WordO) <- "contr.treatment"

m2.ord <- bam(Pos ~ s(Time) + s(Time,by=WordO) + WordO, data=dat)
(smrym2.ord <- summary(m2.ord))

plot(m2.ord, select=2, shade=TRUE, ylim=c(-1,2))
abline(h=0)

# model criticism
gam.check(m2)

# mixed-effects
# random intercepts 
m3 <- bam(Pos ~ Word +
            s(Time, by=Word) + 
            s(Speaker,bs="re"),
          data=dat)
(smrym3 <- summary(m3))

par(mfrow=c(1,2))
plot_smooth(m3, view="Time", 
            plot_all="Word", 
            main="m3", 
            rug=FALSE, 
            rm.ranef=T, 
            ylim=c(-1,2), 
            col=two.colors
            )

plot_diff(m3, 
          view="Time",
          comp=list(Word=c("tenth","tent")), 
          rm.ranef=T, 
          ylim=c(-1,2)
          )

# Random slopes
m4 <- bam(Pos ~ Word + 
            s(Time, by=Word) + 
            s(Speaker,bs="re") + 
            s(Speaker,Word,bs="re"), 
          data=dat
          )
(smrym4 <- summary(m4))

# comparing m3 to m4
compareML(m3,m4)

# plotting the smooth and difference for m4
par(mfrow=c(1,2))
plot_smooth(m4, view="Time", 
            plot_all="Word", 
            main="m4", 
            rug=FALSE, 
            rm.ranef=T, 
            ylim=c(-1,2), 
            col=two.colors
            ) 

plot_diff(m4, 
          view="Time", 
          comp=list(Word=c("tenth",
                           "tent")), 
          rm.ranef=T, 
          ylim=c(-1,2)
          )

# Factor smooths
m5 <- bam(Pos ~ Word + 
            s(Time, by=Word) + 
            s(Speaker,Word,bs="re") + 
            s(Time,Speaker,bs="fs",m=1), 
          data=dat
          )
(smrym5 <- summary(m5))

# plotting m5 simple
plot(m5, select=4)

# comparing m4 and m5
compareML(m4,m5)

# plotting the smooth and difference for m5
par(mfrow=c(1,2))
plot_smooth(m5, 
            view="Time", 
            plot_all="Word",
            main="m5",
            rug=FALSE, 
            rm.ranef=T, 
            ylim=c(-1,2), 
            col=two.colors
            ) 

plot_diff(m5, 
          view="Time",
          comp=list(Word=c("tenth","tent")),
          rm.ranef=T, 
          ylim=c(-1,2)
          )

# generating another model for comparing
m6 <- bam(Pos ~ Word + 
            s(Time, by=Word) + 
            s(Time,
              Speaker,
              by=Word,
              bs="fs",
              m=1),
          data=dat
          )
(smrym6 <- summary(m6))

compareML(m5,m6)

par(mfrow=c(1,2))
plot_smooth(m6, 
            view="Time", 
            plot_all="Word", 
            main="m6", 
            rug=FALSE, 
            rm.ranef=T,
            ylim=c(-1,2), 
            col=two.colors
            ) 

plot_diff(m6, 
          view="Time", 
          comp=list(Word=c("tenth","tent")), 
          rm.ranef=T, 
          ylim=c(-1,2)
          )

# autocorrection
m6acf <- acf_resid(m6)

(rhoval <- m6acf[2]) # autocorrelation at lag 1

m7 <- bam(Pos ~ Word +
            s(Time, by=Word) +
            s(Time,Speaker,by=Word,bs="fs",m=1), 
          data=dat, 
          rho=rhoval, 
          AR.start=dat$start.event) 

acf_resid(m7)

summary(m7)

par(mfrow=c(1,2))
plot_smooth(m7, 
            view="Time", 
            plot_all="Word", 
            main="m7", 
            rug=FALSE, 
            rm.ranef=T, 
            ylim=c(-1,2), 
            col=two.colors
            ) 

plot_diff(m7,
          view="Time",
          comp=list(Word=c("tenth","tent")), 
          rm.ranef=T, 
          ylim=c(-1,2)
          )

m7.alt <- bam(Pos ~ Word + s(Time, by=Word) + s(Time,Speaker,by=Word,bs="fs",m=1), data=dat, 
              rho=rhoval - 0.1, AR.start=dat$start.event) 

acf_resid(m7.alt)

m7.alt2 <- bam(Pos ~ Word + s(Time, by=Word) + s(Time,Speaker,by=Word,bs="fs",m=1), data=dat, 
               rho=0.99, AR.start=dat$start.event) 

acf_resid(m7.alt2)

par(mfrow=c(1,2))
plot_smooth(m7.alt2, view="Time", plot_all="Word", main="m7.alt2", rug=FALSE, rm.ranef=T, ylim=c(-1,2), 
            col=two.colors) 

plot_diff(m7.alt2, view="Time", comp=list(Word=c("tenth","tent")), rm.ranef=T, ylim=c(-1,2))

# Tensor product interaction
m8 <- bam(Pos ~ Word + te(Time, Trial, by=Word) + s(Time,Speaker,by=Word,bs="fs",m=1), data=dat, 
          rho=rhoval, AR.start=dat$start.event)
(smrym8 <- summary(m8))

par(mfrow=c(2,2))
fvisgam(m8, view=c("Time","Trial"), cond=list(Word=c("tent")), main="m8: tent", rm.ranef=T,
        zlim=c(-0.9,1.6), ylim=c(0,600), color="gray")

fvisgam(m8, view=c("Time","Trial"), cond=list(Word=c("tenth")), main="m8: tenth", rm.ranef=T,
        zlim=c(-0.9,1.6), ylim=c(0,600), color="gray") 

plot_diff2(m8, view=c("Time","Trial"), comp=list(Word=c("tenth","tent")), rm.ranef=T,  se=0,
           main="Difference tenth - tent", zlim=c(-0.1,1.2), ylim=c(0,600), color="gray")

par(mfcol=c(3,2))
plot_diff2(m8, view=c("Time","Trial"), comp=list(Word=c("tenth","tent")), rm.ranef=T,  se=0,
           main="Difference tenth - tent", zlim=c(-0.1,1.2), ylim=c(0,600), color="gray",
           cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5) 

abline(h=500,lty=2,lwd=2,col="white")

plot_diff2(m8, view=c("Time","Trial"), comp=list(Word=c("tenth","tent")), rm.ranef=T,  se=0,
           main="Difference tenth - tent", zlim=c(-0.1,1.2), ylim=c(0,600), color="gray",
           cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5) 

abline(h=300,lty=2,lwd=2,col="white")

plot_diff2(m8, view=c("Time","Trial"), comp=list(Word=c("tenth","tent")), rm.ranef=T,  se=0,
           main="Difference tenth - tent", zlim=c(-0.1,1.2), ylim=c(0,600), color="gray",
           cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5) 

abline(h=100,lty=2,lwd=2,col="white")

plot_diff(m8, view="Time", comp=list(Word=c("tenth","tent")), cond=list(Trial=500), rm.ranef=T, 
          main="Difference for trial 500", ylim=c(-1,2), cex.lab=1.5, cex.axis=1.5, 
          cex.main=1.5, cex.sub=1.5) 

abline(h=1.035,lty=3)
abline(v=0.745,lty=3)

plot_diff(m8, view="Time", comp=list(Word=c("tenth","tent")), cond=list(Trial=300), rm.ranef=T, 
          main="Difference for trial 300", ylim=c(-1,2), cex.lab=1.5, cex.axis=1.5, 
          cex.main=1.5, cex.sub=1.5)

abline(h=1.035,lty=3)
abline(v=0.745,lty=3)

plot_diff(m8, view="Time", comp=list(Word=c("tenth","tent")), cond=list(Trial=100), rm.ranef=T, 
          main="Difference for trial 100", ylim=c(-1,2), cex.lab=1.5, cex.axis=1.5, 
          cex.main=1.5, cex.sub=1.5)

abline(h=1.035,lty=3)
abline(v=0.745,lty=3)

# Decomposition of tensor
m8.dc <- bam(Pos ~ Word + s(Time, by=Word) + s(Trial, by=Word) + ti(Time, Trial, by=Word) + 
               s(Time,Speaker,by=Word,bs="fs",m=1), data=dat, rho=rhoval, 
             AR.start=dat$start.event)
(smrym8.dc <- summary(m8.dc))

par(mfrow=c(2,2))
plot(m8.dc, select=1, shade=T, rug=F, cex.lab=1.35, cex.axis=1.35, cex.main=1.35, cex.sub=1.35,
     ylim=c(-1.5,1.5))
abline(h=0)
plot(m8.dc, select=2, shade=T, rug=F, cex.lab=1.35, cex.axis=1.35, cex.main=1.35, cex.sub=1.35,
     ylim=c(-1.5,1.5))
abline(h=0)
plot(m8.dc, select=3, shade=T, rug=F, cex.lab=1.35, cex.axis=1.35, cex.main=1.35, cex.sub=1.35,
     ylim=c(-1.5,1.5))
abline(h=0)
plot(m8.dc, select=4, shade=T, rug=F, cex.lab=1.35, cex.axis=1.35, cex.main=1.35, cex.sub=1.35,
     ylim=c(-1.5,1.5))
abline(h=0)