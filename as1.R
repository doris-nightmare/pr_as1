library(MASS)

mu1<-matrix(2,2)
mu2 <-matrix(2,-2)
mu3 <-matrix(c(-2,0))

sigma <-matrix(c(0.25,0,0,0.25),nrow=2)

omega1 <- solve(sigma) %*% mu1
omega2 <- solve(sigma) %*% mu2
omega3 <- solve(sigma) %*% mu3

ln13 <- log(0.333333, base = exp(1))
omega10 <- ln13 - 0.5 * t(mu1) %*% solve(sigma) %*% mu1
omega20 <- ln13 - 0.5 * t(mu2) %*% solve(sigma) %*% mu2
omega30 <- ln13 - 0.5 * t(mu3) %*% solve(sigma) %*% mu3

#generate random points
p1 <- mvrnorm(100, mu1, sigma)
p2 <- mvrnorm(100, mu2, sigma)
p3 <- mvrnorm(100, mu3, sigma)

#plot the points
plot.new()
plot.window(xlim=c(-4,3), ylim=c(-4,3))
axis(1)
axis(2)
title(main="Three Classes")
title(xlab="x1")
title(ylab="x2")
points(p1, pch = 1, col = 1)
points(p2, pch = 2, col = 2)
points(p3, pch = 3, col = 3)

#decision hyperplance
abline(h = 0, col = 1) # 12
text(0,0,"l1",col=1)
abline(a=-1, b= 2, col = 2) # 23
text(-1,-3,"l2",col=2)
abline(a=1, b= -2, col = 3) # 13
text(-1,3,"l3",col=3)
legend("bottomleft",c("1",'2','3'),pch=c(1,2,3),col=c(1,2,3))

#calculate the estimates
mu1.hat = as.matrix(apply(p1, 2, mean))
mu2.hat = as.matrix(apply(p2, 2, mean))
mu3.hat = as.matrix(apply(p3, 2, mean))

calc.sigma <- function(points, mu){
    n <- nrow(points)
    col <- apply(t(points),2,function(x) {x-mu})
    sigma <- col %*% t(col) / n
    return(sigma)
}

sigma1.hat = calc.sigma(p1, mu1.hat)
sigma2.hat = calc.sigma(p2, mu2.hat)
sigma3.hat = calc.sigma(p3, mu3.hat)

W1 <- 0.5 * solve(sigma1.hat)
W2 <- 0.5 * solve(sigma2.hat)
W3 <- 0.5 * solve(sigma3.hat)

w1 <- solve(sigma1.hat)  %*% mu1.hat
w2 <- solve(sigma2.hat)  %*% mu2.hat
w3 <- solve(sigma3.hat)  %*% mu3.hat

calc.omega0 <- function(mu, sigma){
    omega0 <- -0.5 * t(mu) %*% solve(sigma) %*% mu 
        - 0.5 * log(det(sigma),base = exp(1)) +log(0.3333,base=exp(1))
    return(omega0)
}

omega10 <- calc.omega0(mu1.hat, sigma1.hat)
omega20 <- calc.omega0(mu2.hat, sigma2.hat)
omega30 <- calc.omega0(mu3.hat, sigma3.hat)

#plot the points
plot.new()
plot.window(xlim=c(-4,3), ylim=c(-4,3))
axis(1)
axis(2)
title(main="Three Classes")
title(xlab="x1")
title(ylab="x2")
points(p1, pch = 1, col = 1)
points(p2, pch = 2, col = 2)
points(p3, pch = 3, col = 3)
text(0,0,"l1'",col=1)
text(-1,-3,"l2'",col=2)
text(-1,3,"l3'",col=3)
legend("bottomleft",c("1",'2','3'),pch=c(1,2,3),col=c(1,2,3))

#plot implict function
bound1 <- function(x,y){-3*x^2+0.58*x*y+0.182*y^2-3.203*x+18.062*y+4.929}
bound2 <- function(x,y){-0.013*x^2-1.12*x*y+0.298*y^2+19.016*x-9.603*y-11.070}
bound3 <- function(x,y){-0.312*x^2-0.536*x*y+0.48*y^2+15.813*x+8.459*y-6.142}
x<-seq(-4,3,length=1000)
y<-seq(-4,3,length=1000)
contour(x,y,outer(x,y,bound1),level=0, add=T)
contour(x,y,outer(x,y,bound2), col = 2, level=0, add=T)
contour(x,y,outer(x,y,bound3), col = 3, level=0, add=T)


#generate new data points
p1.new <- mvrnorm(100, mu1, sigma)
p2.new <- mvrnorm(100, mu2, sigma)
p3.new <- mvrnorm(100, mu3, sigma)

# new estimates
#calculate the estimates
mu1.new.hat = as.matrix(apply(p1.new, 2, mean))
mu2.new.hat = as.matrix(apply(p2.new, 2, mean))
mu3.new.hat = as.matrix(apply(p3.new, 2, mean))

sigma1.new.hat = calc.sigma(p1.new, mu1.new.hat)
sigma2.new.hat = calc.sigma(p2.new, mu2.new.hat)
sigma3.new.hat = calc.sigma(p3.new, mu3.new.hat)

W1.new <- 0.5 * solve(sigma1.new.hat)
W2.new <- 0.5 * solve(sigma2.new.hat)
W3.new <- 0.5 * solve(sigma3.new.hat)

w1.new <- solve(sigma1.new.hat)  %*% mu1.new.hat
w2.new <- solve(sigma2.new.hat)  %*% mu2.new.hat
w3.new <- solve(sigma3.new.hat)  %*% mu3.new.hat

omega10.new <- calc.omega0(mu1.new.hat, sigma1.new.hat)
omega20.new <- calc.omega0(mu2.new.hat, sigma2.new.hat)
omega30.new <- calc.omega0(mu3.new.hat, sigma3.new.hat)

#plot new points
plot.new()
plot.window(xlim=c(-4,3), ylim=c(-4,3))
axis(1)
axis(2)
title(main="Another Plot for Three Classes")
title(xlab="x1")
title(ylab="x2")
points(p1.new, pch = 1, col = 1)
points(p2.new, pch = 2, col = 2)
points(p3.new, pch = 3, col = 3)
text(0,0,"l1''",col=1)
text(-1,-3,"l2''",col=2)
text(-1,3,"l3''",col=3)
legend("bottomleft",c("1",'2','3'),pch=c(1,2,3),col=c(1,2,3))

#plot implict function
bound1.new <- function(x,y){0.557*x^2-0.648*x*y+0.174*y^2+1.036*x+16.155*y-0.294}
bound2.new <- function(x,y){0.128*x^2+0.46*x*y+0.076*y^2+15.522*x-9.542*y-9.387}
bound3.new <- function(x,y){0.686*x^2-1.864*x*y+0.249*y^2+16.558*x+6.614*y-9.68}
x<-seq(-4,3,length=1000)
y<-seq(-4,3,length=1000)
contour(x,y,outer(x,y,bound1.new),level=0, add=T)
contour(x,y,outer(x,y,bound2.new), col = 2, level=0, add=T)
contour(x,y,outer(x,y,bound3.new), col = 3, level=0, add=T)

save(list=ls(),file='as1.data')
