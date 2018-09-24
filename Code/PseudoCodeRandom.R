# Takes approx 1min 14s to produce 1 simulation
st <- proc.time()
PseudoSim( DirName = "Test2", Pars= c(0.0015, 0.07, 0.0028, 0.013, 0.045, 0.16, 7e-5, 0.1, 0.025, 0.02))
proc.time() - st

# Read in parameter configurations for LHD
# Each parameter is one column in $theta
# order of parameters: ACTIN FK GI_SPEED bB dB diffB KM sA BASAL_A diffA dA 
LHDtheta <- read.mat("LHC/MS_Thetas.mat")
ParNames <- c("ACTIN_FK", "GI_SPEED", "bB", "dB", "diffB", "KM", "sA", "BASAL_A", "diffA", "dA")



Par######## SUMMARY STATISTICS
#### GI activator


## GIs matrix containing GI data for each LHD parameter
GIs <- matrix( nrow=1000, ncol=100)

for ( i in 1:100)
{
    GI <- PseudoRead( Path = paste("MS", i, sep=""), Output = "GI" )
    GIs[,i] <- GI$V1
}

# Plotting GIs 


# time index - what times to plot at

# flat vs 
tind <- c(3:25)

# beginning - no overlap in the first time interval (clear pattern)
tind <- c(1:3)

# 



plot( GIs[tind, 1], x = tind, type='l', ylab="Global inhibitor", xlab ="GI_speed", main="Time 0-1000",  ylim = c(0, max(GIs)), col=alpha("black", 0.3))

for (i in 2:100)
{
    lines( GIs[tind, i], x = tind, col = alpha(i, 0.3))
}

# boxplots of GIs sorted by gi_speed
#gi_speed <- LHDtheta$theta[,3]
#ind <- rank(gi_speed)
for(k in 1:10) {
    ptest <- LHDtheta$theta[,k]
    ind <- rank(ptest)
    
    GIs2 <- matrix( nrow=1000, ncol=100)
    for (i in 1:100)
        GIs2[, i] <-  GIs [, which(ind==i)]
    
    png(paste("GI", ParNames[k], ".png", sep=""), width = 1200, height = 600)
   
    par(mfrow=c(2,1), mar=c(2,3,2,1))
    boxplot(GIs2, cex=0.2, main= "Global inhibitor")
    plot(sort(LHDtheta$theta[,k]), x=c(1:100), col="red", type="l", main = ParNames[k])
    
    dev.off()
    
}





    
    