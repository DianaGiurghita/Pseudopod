# Takes approx 1min 14s to produce 1 simulation
st <- proc.time()
PseudoSim( DirName = "Test2", Pars= c(0.0015, 0.07, 0.0028, 0.013, 0.045, 0.16, 7e-5, 0.1, 0.025, 0.02))
proc.time() - st

# Read in parameter configurations for LHD
# Each parameter is one column in $theta
# order of parameters: ACTIN FK GI_SPEED bB dB diffB KM sA BASAL_A diffA dA 
LHDtheta <- read.mat("LHC/MS_Thetas.mat")
ParNames <- c("ACTIN_FK", "GI_SPEED", "bB", "dB", "diffB", "KM", "sA", "BASAL_A", "diffA", "dA")



######## SUMMARY STATISTICS
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
   
    par( mfrow=c(2,1), mar=c(2,3,2,1))
    boxplot( GIs2, cex=0.2, main= "Global inhibitor")
    plot( sort( LHDtheta$theta[,k]), x=c(1:100), col="red", type="l", main = ParNames[k])
    
    dev.off()
    
}





##### PEAK ANALYSIS
##### Detecting number of maximum of LI and LA

## Plotting settings for restoring later
op <- par(no.readonly = TRUE)

# Time index 
tind <- c( 100, 300, 500, 1000)

for ( j in 1:100) {
        
    path <- paste("MS", j, sep="")  
    
    # Read in data
    theta1_LI <- PseudoRead( Path = path, Output = "LI" )
    theta1_LA <- PseudoRead( Path = path, Output = "LA" )
    theta1_X <- PseudoRead( Path = path, Output = "X" )
    theta1_Y <- PseudoRead( Path = path, Output = "Y" )


   # pdf( paste("Theta", j, ".pdf", sep=""))
    
    # Plotting settings for LI/LA peaks plots
    par( mfrow = c(4,3), mar = c(4,4,2,0.5), oma = c(0, 0, 3, 0))
    val <- 0.15
    
    # Plot peak locations of LI and LA on cell membrane
    for (i in 1:4){
        
        ## get index of positions without NAs
        theta1_indLI <- ! is.na( theta1_LI [tind[i], ])
        theta1_indLA <- ! is.na( theta1_LA [tind[i], ])
        
        vecLI <- as.numeric (matrix ( theta1_LI[tind[i], theta1_indLI ]))
        vecLA <- as.numeric (matrix ( theta1_LA[tind[i], theta1_indLA ]))
        
        # where the peaks of LI are, col 1 = LI value, col2 = node location
        # peaksLI <- findpeaks( vecLI, threshold = 1.1 * min(vecLI), ndowns=5, nups=5 )
        # peaksLA <- findpeaks( vecLA, threshold = 1.1 * min(vecLA), ndowns=5, nups=5 )
        
        peaksLI <- findpeaks( vecLI )
        peaksLA <- findpeaks( vecLA )

        # peaksLI <- findpeaks( vecLI, minpeakheight = val* (max(vecLI) - min(vecLI)) )
        # peaksLA <- findpeaks( vecLA, minpeakheight = val* (max(vecLA) - min(vecLA)) )
        # 
        
        adjLI <- 0.1 * max(vecLI) 
        adjLA <- 0.1 * max(vecLA)  
        
        # Plot Cell
        xc <- matrix( theta1_X[tind[i], ])
        yc <- matrix( theta1_Y[tind[i], ])
        
        CellTitle = paste("Cell contour Time", c(tind[i]))
        plot(x=xc, y=yc, type="l", xlab = "x coordinate", ylab="y coordinate", main = CellTitle )
            points ( xc[peaksLA[,2]], yc[peaksLA[,2]], pch=1, lwd=2, col="#41b6c4")
            
            points ( xc[peaksLI[,2]], yc[peaksLI[,2]], pch=20, col="#c51b8a")
            text (peaksLI[,2]+adjLI, peaksLI[,1]+adjLA, pch= paste( c(1:nrow(peaksLI))),  col="#c51b8a")
            
            grid()
        
        # Plot LI
        plot( vecLI , type = "l", ylab = "LI", main= "LI", col="navy", xlab="Node number")
            grid()
            points (peaksLI[,2], peaksLI[,1], pch = 20, col="#c51b8a")
            text (peaksLI[,2], peaksLI[,1]+adjLI, pch= paste( c(1:nrow(peaksLI))),  col="#c51b8a")
        #abline( h= val* (max(vecLI) - min(vecLI)), col = "gray" )
        
        # Plot LA
        plot( vecLA , type = "l", ylab = "LA", main= "LA", col="navy", xlab="Node number")
            grid()
            points (peaksLA[,2], peaksLA[,1], pch = 1, col="#41b6c4")
            text (peaksLA[,2], peaksLA[,1]+adjLA, pch= paste( c(1:nrow(peaksLA))),  col="#41b6c4")
        #abline( h= val* (max(vecLA) - min(vecLA)), col = "gray" )
        
        mtext( paste( "Theta", j , ": Peak location of LI and LA on cell membrane", sep=""), outer = TRUE, side = 3, cex = 1.2, line = 1, adj = 0 )
    
        npeakLA[i] <- nrow( peaksLA > val* (max(vecLA) - min(vecLA)))
        npeakLI[i] <- nrow( peaksLI > val* (max(vecLI) - min(vecLI))) 
        }
    
   # dev.off()
    
    pdf( paste("Theta", j, ".pdf", sep=""))
    par(mfrow = c(1,2))
    plot(y=npeakLA, x=tind,  type='b', main =  "LA", xlab = "Time", ylab = "Number of peaks")
    plot(y=npeakLI, x=tind,  type='b', main =  "LI", xlab = "Time", ylab = "Number of peaks")
    mtext( paste( "Theta", j , " number of peaks over time", sep=""), outer = TRUE, side = 3, cex = 1.2, line = 1, adj = 0 )
    dev.off()
    

}



## resotre plot settings
par( op)
