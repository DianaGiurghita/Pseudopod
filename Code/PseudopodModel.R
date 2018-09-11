### Columns names needed 

# detects maximum number of elements in the files (in each row) and takes the max of that
Rdim <- max(count.fields("cellX.csv", sep = ","))

# forms column names N1, ... N Rnum - this is needed to force import missing values at the end of 
# the beginning rows.
CN <- paste0("N", seq_len(Rdim))




### Global inhibitor value
### cellGI.csv 
# 1 value at each timepoint 

GI <- read.csv("cellGI.csv", header = F)

    par( mfrow = c(1,2))
    plot(GI$V1, type = 'l', ylab="Global inhibitor", xlab ="Time", main="Time 0-1000")
    plot(GI$V1[500:999], x=c(500: 999), type = 'l', ylab="Global inhibitor", xlab ="Time", main="Time 500-1000")
    

    
### Local activator concentration at each node around the boundary 
### cellLA.csv

LA <- read.csv("cellLA.csv", header = F, col.names = CN, fill = T )
  
    par( mfrow = c(2,2))
    plot(matrix(LA[1,]), type = 'l', ylab="Local activator", main="t=1")
    plot(matrix(LA[300,]), type = 'l', ylab="Local activator", main="t=300")
    plot(matrix(LA[600,]), type = 'l', ylab="Local activator", main="t=600")
    plot(matrix(LA[1000,]), type = 'l', ylab="Local activator", main="t=1000")
    
    

### Local inhibitor concentration at each node around the boundary 
### cellLI.csv

LI <- read.csv("cellLI.csv", header = F, col.names = CN, fill = T )

    par( mfrow = c(2,2))
    plot(matrix(LI[1,]), type = 'l', ylab="Local inhibitor", main="t=1")
    plot(matrix(LI[300,]), type = 'l', ylab="Local inhibitor", main="t=300")
    plot(matrix(LI[600,]), type = 'l', ylab="Local inhibitor", main="t=600")
    plot(matrix(LI[1000,]), type = 'l', ylab="Local inhibitor", main="t=1000")


    
### Local values of s(x,t) at each node on the boundary
### cellS.csv
    
S <- read.csv("cellS.csv", header = F, col.names = CN, fill = T )
 
    par( mfrow = c(2,2))
    plot(matrix(S[1,]), type = 'l', ylab="Stimulus strength", main="t=1")
    plot(matrix(S[300,]), type = 'l', ylab="Stimulus strength", main="t=300")
    plot(matrix(S[600,]), type = 'l', ylab="Stimulus strength", main="t=600")
    plot(matrix(S[1000,]), type = 'l', ylab="Stimulus strength", main="t=1000")



   
### Cell x coordinate and y coordinates 
### cellX.csv 
### cellY.csv
    
cellX <- read.csv("cellX.csv", header = F, col.names = CN, fill = T )
cellY <- read.csv("cellY.csv", header = F, col.names = CN, fill = T )


    par( mfrow = c(2,2))
    xc <- matrix( cellX[1,])
    yc <- matrix( cellY[1,])
    plot(x=xc, y=yc, type="l", xlab = "x coordinate", ylab="y coordinate", main = "Timepoint 1")

    xc <- matrix( cellX[300,1:Rdim])
    yc <- matrix( cellY[300,1:Rdim])
    plot(x=xc, y=yc, type="l", xlab = "x coordinate", ylab="y coordinate", main = "Timepoint 300")

    xc <- matrix( cellX[600,1:Rdim])
    yc <- matrix( cellY[600,1:Rdim])
    plot(x=xc, y=yc, type="l", xlab = "x coordinate", ylab="y coordinate", main = "Timepoint 600")

    xc <- matrix( cellX[1000,1:Rdim])
    yc <- matrix( cellY[1000,1:Rdim])
    plot(x=xc, y=yc, type="l", xlab = "x coordinate", ylab="y coordinate", main = "Timepoint 1000")

## x vs t
    
    par( mfrow = c(2,4))
    
    plot(cellX$N1, type='l', xlab = "Time", ylab = "x", main = "Node 1")
    plot(cellX$N30, type='l', xlab = "Time", ylab = "x", main = "Node 30")
    plot(cellX$N60, type='l', xlab = "Time", ylab = "x", main = "Node 60")
    plot(cellX$N80, type='l', xlab = "Time", ylab = "x", main = "Node 80")
    
    plot(y=cellX[1,], x=c(1:86), type='l', xlab = "Node", ylab = "x", main = "Time=1")
    plot(y=cellX[300,], x=c(1:86), type='l', xlab = "Node", ylab = "x", main = "Time=300")
    plot(y=cellX[600,], x=c(1:86), type='l', xlab = "Node", ylab = "x", main = "Time=600")
    plot(y=cellX[1000,], x=c(1:86), type='l', xlab = "Node", ylab = "x", main = "Time=1000")
    

## y vs t    
    plot(cellY$N1, type='l', xlab = "Time", ylab = "y", main = "Node 1")
    plot(cellY$N30, type='l', xlab = "Time", ylab = "y", main = "Node 30")
    plot(cellY$N60, type='l', xlab = "Time", ylab = "y", main = "Node 60")
    plot(cellY$N80, type='l', xlab = "Time", ylab = "y", main = "Node 80")
    
    plot(y=cellY[1,], x=c(1:86), type='l', xlab = "Node", ylab = "y", main = "Time=1")
    plot(y=cellY[300,], x=c(1:86), type='l', xlab = "Node", ylab = "y", main = "Time=300")
    plot(y=cellY[600,], x=c(1:86), type='l', xlab = "Node", ylab = "y", main = "Time=600")
    plot(y=cellY[1000,], x=c(1:86), type='l', xlab = "Node", ylab = "y", main = "Time=1000")
    