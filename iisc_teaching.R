

library(ggplot2)
library(dplyr)

out<-read.csv(file = "Extracted_Data_from_Graph.csv")

plot(out$Rotifers..mL..1.)
plot(out$Rotifers..mL..1..1)


Algal_d<-c(14,12, 9, 8, 8.5, 4, 3.8, 3.2,4,4.2,3.8,3.1, 3, 2.9,2.8,2.6, 1.4, 1.1, 1, 0.9, 0.5, 0.1,0.0,0)*10^6
rotifer_d<-c(0.1,1,3,5,13.5,10,6,7,12,8,13,15.6,14.6,12,7.2,3,2,2.2,1,0.2,0.1,0.0,0,0)

time<-1:24
data.frame(Algae= Algal_d,time,rotifer_d)

plot(time, Algal_d, cex=4, pch=20, col="black", ylim=c(0,17),xlab = "Time")
points( rotifer_d, cex=4,pch=21,col="firebrick",xlab="Time")

plot(time, Algal_d, typ="l", lwd=3, col="black", ylim=c(0,17),xlab = "Time", ylab="Species Density")
lines(rotifer_d,col="firebrick",lwd=4)



plot(time, Algal_d, type="l", col="blue", lwd=2, 
     ylab="Algal Density (cells/mL)", xlab="Time (days)",
     ylim=c(0, max(Algal_d))) # Adjust y-axis limit accordingly

# Add the second dataset (Rotifer_d) with a secondary y-axis
par(new=TRUE)  # Allow for a new plot on top
plot(time, rotifer_d, type="l", col="red", lwd=2, 
     axes=FALSE, xlab="", ylab="", ylim=c(0, max(rotifer_d))) 

# Add the right y-axis
axis(4)
mtext("Rotifer Density (individuals/mL)", side=4, line=3)

# Add a legend
legend("topright", legend=c("Algal Density", "Rotifer Density"), 
       col=c("blue", "red"), lty=1, lwd=2)

