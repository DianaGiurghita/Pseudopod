###### Procrustes analysis


### Step 1: building matrices with the same dimensions in time (x,y - coord). 
### First check whether there is substantial loss of information from scaling down - the obvious method.

## Plots for every theta that compare the scaled down version of x,y time serie with original one
for (j in 1:100) {

    # Define paths where cell theta j is located
    path <- paste("MS", j, sep="")  
    fpath <- paste("~/ownCloud/Phd Diana/Cside 2018/Simulator/Init_LHS/", path, "/Cell1/cellX.csv", sep="")
    
    
    # Get minimum number of nodes for cell theta j from the csv file 
    minlen <- min( count.fields( fpath, sep = ","))
    
    
    # Read in x and y coords for theta js
    theta1_X <- PseudoRead( Path = path, Output = "X" )
    theta1_Y <- PseudoRead( Path = path, Output = "Y" )
    
    # Define times for plotting cell membrane at
    Ktimes <- c(1, 200, 400, 600, 800, 1000)
    
    
    # Open pdf file to write plot in 
    pdf( paste("Theta", j, "_thin.pdf", sep=""), width = 7, height = 4)
    
    # Set plotting area
    par( mfrow = c(1,3), mar = c(4,4,2,0.5), oma = c(0, 0, 3, 0))
    
    
    ### Plot j cell over time 
    {
      
        K <- length( Ktimes)    
        
            
        # plot x,y
        plot ( x = c( min(theta1_X, na.rm = TRUE), max( theta1_X, na.rm = TRUE)), 
               y = c( min(theta1_Y, na.rm = TRUE), max( theta1_Y, na.rm = TRUE)), 
               type="n", ylab = "Y coordinate", xlab = "X coordinate")
        grid()
        
        for ( k in 1:K ) {
            
            x <- theta1_X[ Ktimes[k], ]
            y <- theta1_Y[ Ktimes[k], ]
            
            
            # Get number of nodes for each timepoint Ktimes[k]
            timelen <- length( na.omit( theta1_X [Ktimes[k] ,] ))
            
            # Index of equally spaced positions of length = minlen
            ind <- round( seq( 1, timelen, len = minlen))
            
            # Thinned series as lines
            lines ( x[ind], y[ind], col = alpha(brewer.pal(n=9, name="Oranges")[k+1], 0.7), lwd=2  )
            
            # Original series as points
            points( x, y, col = brewer.pal(n=9, name="Blues")[k+1], lwd=1, pch=20, cex=0.3 ) }
        
        
        # plot x
        plot ( x = c( 1, dim(theta1_X)[2] ), 
               y = c( min(theta1_X, na.rm = TRUE), max( theta1_X, na.rm = TRUE)), 
               type="n", ylab = "X coordinate", xlab = "Node index")
        grid()
        
        for ( k in 1:K ) {
            
            x <- theta1_X[ Ktimes[k], ]
            
            # Get number of nodes for each timepoint Ktimes[k]
            timelen <- length( na.omit( theta1_X [Ktimes[k] ,] ))
            
            # Index of equally spaced positions of length = minlen
            ind <- round( seq( 1, timelen, len = minlen))
            
            # Thinned series as lines
            lines ( x[ind], x=ind, col = alpha( brewer.pal(n=9, name="Oranges")[k+1], 0.7), lwd=2  )
            
            # Original series as points
            points( x~ c(1:length(x)), col = brewer.pal(n=9, name="Blues")[k+1], lwd=1, pch=20, cex=0.3 ) }
        
        
        # plot y
        plot ( x = c( 1, dim(theta1_Y)[2] ), 
               y = c( min(theta1_Y, na.rm = TRUE), max( theta1_Y, na.rm = TRUE)), 
               type="n", ylab = "Y coordinate", xlab = "Node index")
        grid()
        
        for ( k in 1:K ) {
            
            x <- theta1_Y[ Ktimes[k], ]
            
            # Get number of nodes for each timepoint Ktimes[k]
            timelen <- length( na.omit( theta1_X [Ktimes[k] ,] ))
            
            # Index of equally spaced positions of length = minlen
            ind <- round( seq( 1, timelen, len = minlen))
            
            # Thinned series as lines
            lines ( x[ind], x=ind, col = alpha( brewer.pal(n=9, name="Oranges")[k+1], 0.7), lwd=2  )
            
            # Original series as points
            points( x~ c(1:length(x)), col = brewer.pal(n=9, name="Blues")[k+1], lwd=1, pch=20, cex=0.3 ) }
        
        mtext( paste( "Theta", j , ": Thinning of x, y coordinate series", sep=""), outer = TRUE, side = 3, cex = 1.2, line = 1, adj =  ) }
        
    
    # Add legend
    legend("bottomright", c("Original", "Thinned"), pch = c( 20, NA), lwd=c(NA,2), col = c(brewer.pal(n=9, name="Blues")[4], brewer.pal(n=9, name="Oranges")[4]), 
           bty ='n', inset=c(0,1), xpd=TRUE, horiz = T)
        
    dev.off()
}

### Step 2: apply Procrustes analysis 

library( shapes)
# demo(shapes)


# Cell 1

j <- 1
res <- thinXY(j=j)

X <- res$X
Y <- res$Y


# Rearrange data in an array of size: 79 x 2 x 1000 (i.e.: nodes x dimensions x cell_at_timepoint)

cell1 <- array (data = NA, dim = c( ncol(X), 2, 1000  ))

for (i in 1:1000){
    cell1[ , 1, i] <- X[i, ]
    cell1[ , 2, i] <- Y[i, ] }


# Procrustes analysis WITH scaling
proc_cell1 <- procGPA(cell1)


# Procrustes analysis WITHOUT scaling
proc_cell1_ns <- procGPA(cell1, scale = FALSE)


# Plotting
# Open pdf file to write plot in 
pdf( paste("Theta", j, "_thin.pdf", sep=""), width = 7, height = 7)

par(mfrow = c(2,1))

plotshapes( proc_cell1$rotated[,, Ktimes], joinline=c(1:79, 1), symbol = 20, color = alpha( "blue", 0.7))
            
mtext( paste("Theta" , j, ": Procrustes analysis with scaling", sep="") )


plotshapes( proc_cell1_ns$rotated[,, Ktimes], joinline=c(1:79, 1), symbol = 20, color = alpha( "red", 0.7))

mtext( paste("Theta" , j, ": Procrustes analysis without scaling", sep="") )
dev.off()


