### Fitting a cubic splines model on x 


par( mfrow = c(2,4))

plot(cellX$N1, type='l', xlab = "Time", ylab = "x", main = "Node 1")
plot(cellX$N30, type='l', xlab = "Time", ylab = "x", main = "Node 30")
plot(cellX$N60, type='l', xlab = "Time", ylab = "x", main = "Node 60")
plot(cellX$N80, type='l', xlab = "Time", ylab = "x", main = "Node 80")

plot(y=cellX[1,], x=c(1:86), type='l', xlab = "Node", ylab = "x", main = "Time=1")
plot(y=cellX[300,], x=c(1:86), type='l', xlab = "Node", ylab = "x", main = "Time=300")
plot(y=cellX[600,], x=c(1:86), type='l', xlab = "Node", ylab = "x", main = "Time=600")
plot(y=cellX[1000,], x=c(1:86), type='l', xlab = "Node", ylab = "x", main = "Time=1000")

x <- as.numeric(c(1:86))
y1 <- as.numeric(cellX[1, 1:86]) 
y2 <- as.numeric(cellY[1, 1:86])
my1 <- lm( y1 ~ bs(x, df = 5 )) 
fy1 <- fitted(my1)

my2 <- lm( y2 ~ bs(x, df = 5 )) 
fy2 <- fitted(my2)

par (mfrow = c(3,1))
plot(y=cellX[1,], x=c(1:86), type='l', xlab = "Node", ylab = "x", main = "Time=1", lwd =2)
    lines( fy1, col="red", lty =2, lwd =2)
   
plot(y=cellY[1,], x=c(1:86), type='l', xlab = "Node", ylab = "y", main = "Time=1")
    lines( fy2, col="red", lty =2, lwd =2)

    xc <- matrix( cellX[1,])
    yc <- matrix( cellY[1,])
    plot(x=xc, y=yc, type="l", xlab = "x coordinate", ylab="y coordinate", main = "Timepoint 1", lwd=2)
    lines( x=fy1, y=fy2, col="red", lty =2, lwd =2)

mt1 <- lm()
          