##################################################################
# Statistical Analysis of Gowning/Gloving Behavior Change Models #
##################################################################

## Import Packages and Import Data
library(vioplot)
library(ggplot2)

mrsa <- as.data.frame(read.csv("ContactChangeOutcomes_mrsa.csv"))

# Plot baseline vs. targeted distribution
# Convert to patient days: 12 patients * 365 = 4380 pd
par(mfrow=c(1,1))
par(mar=c(5.1,5.1,4.1,2.1))
pd_den <- density((mrsa$Baseline/(18*365))*1000, adjust=5,kernel="gaussian",from=0)
plot(pd_den,main="",ylab="Density",xlab="Cases per 1000 Patient-Days",cex.lab=1.5,lwd=5,lty=1)
abline(v=mean((mrsa$Baseline/(18*365))*1000),lwd=3,col="red")
abline(v=5.94,lwd=3,lty=2)
legend("topright", c("Smoothed Density","Median Incidence Rate","Target Incidence Rate"), lwd=3, col=c("black","red","black"),
       lty=c(1,1,2), pch=c(NA,NA), bty='n',inset=0.02, cex=1.3)

# Reorder data as a single column with a scenario variable
bl <- mrsa[1]
colnames(bl) <- "Incident"
bl$Scenario <- factor("Baseline")
bl$rate <- bl$Incident/(18*365)*1000

lowcontact <- mrsa[2]
colnames(lowcontact) <- "Incident"
lowcontact$Scenario <- factor("Reduced Contact")
lowcontact$rate <- lowcontact$Incident/(18*365)*1000

efficient <- mrsa[3]
colnames(efficient) <- "Incident"
efficient$Scenario <- factor("Efficient Contact")
efficient$rate <- efficient$Incident/(18*365)*1000

redhh <- mrsa[4]
colnames(redhh) <- "Incident"
redhh$Scenario <- factor("Reduced Contact + HH")
redhh$rate <- redhh$Incident/(18*365)*1000

effhh <- mrsa[5]
colnames(effhh) <- "Incident"
effhh$Scenario <- factor("Efficient Contact + HH")
effhh$rate <- effhh$Incident/(18*365)*1000

collected <- rbind(bl,lowcontact,efficient)

# ANOVA and Tukey HSD for the differences in the main scenarios
IncidentTest <- aov(collected$rate ~ collected$Scenario)

print(summary(IncidentTest))
print(TukeyHSD(IncidentTest))

print(mean(bl$rate))
print(sd(bl$rate))
print(mean(lowcontact$rate))
print(sd(lowcontact$rate))
print(mean(efficient$rate))
print(sd(efficient$rate))
print(mean(redhh$rate))
print(sd(redhh$rate))
print(mean(effhh$rate))
print(sd(effhh$rate))

## Figure 4 - Cumultive incident and recurrent cases
par(mfrow=c(1,1))
vioplot(bl$rate,lowcontact$rate,efficient$rate,names=c("Baseline","Reduced Contact","Efficient Contact"),col="Grey90",drawRect=FALSE)
title(ylab="MRSA Acquisitions (per 1,000 patient-days)",cex.lab=1.5)
title(xlab="Scenario",cex.lab=1.5)
segments(0.75,mean(bl$rate),1.25,mean(bl$rate),col="black",lwd=3)
segments(1.75,mean(lowcontact$rate),2.25,mean(lowcontact$rate),col="black",lwd=3)
segments(2.75,mean(efficient$rate),3.25,mean(efficient$rate),col="black",lwd=3)
legend("topleft", c("Scenario Mean"), lwd=3,col=c("black"),lty=1, bty='n', cex=1.25)

## Analysis of the sweep data
sweep <- read.csv("ContactChangeOutcomes_MRSA_Sweep.csv")
reducedsweep <- glm(ReducedCases ~ ReducedPer,family="poisson", data=sweep)
summary(reducedsweep)

efficientsweep <- glm(EfficientCases ~ EfficientPer,family="poisson", data=sweep)
summary(efficientsweep)

# Plot
par(mfrow=c(2,1))
predict_percent <- as.data.frame(seq(from=0,to=50,by=0.5))
colnames(predict_percent) <- "ReducedPer"
predict_percent$EfficientPer <- seq(from=0, to=50, by=0.50)

reducedpredict <- predict(reducedsweep,newdata=predict_percent,type="response")
plot(sweep$ReducedPer,sweep$ReducedCases/(18*365)*1000,pch=21,col="grey90",bg="grey90",
     xlab="Percent Reduction in Contact Rate",ylab="MRSA Acquisitions (per 1,000 patient-days)",cex.lab=1.25)
lines(predict_percent$ReducedPer,reducedpredict/(18*365)*1000,col="black",lwd=5)
legend("topleft", c("Poisson Fit"), lwd=3,col=c("black"),lty=1, bty='n', cex=1.1)

effpredict <- predict(efficientsweep,newdata=predict_percent,type="response")
plot(sweep$EfficientPer,sweep$EfficientCases/(18*365)*1000,pch=21,col="grey90",bg="grey90",
     xlab="Percent Reduction in Contact Rate",ylab="MRSA Acquisitions (per 1,000 patient-days)",cex.lab=1.25)
lines(predict_percent$EfficientPer,effpredict/(18*365)*1000,col="black",lwd=5)
legend("topleft", c("Poisson Fit"), lwd=3,col=c("black"),lty=1, bty='n', cex=1.1)

