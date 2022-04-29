rm(list=ls())
library(rgl)
require("parallel")
require("Rcpp")

# function
e = function(x, y, k = 120) {
  ifelse (test=((abs(x) > k) | (abs(y) > k)),
          yes=0,
          no=(1 - x/k)*(1 - y/k)*(1 + x/k)*(1 + y/k)*
            (-(y + 47) * sin(sqrt(abs(y + x/2 + 47))) - x*sin(sqrt(abs(x - (y + 47)))))
  )
}

#Complicated derivative is found using:
#https://www.derivative-calculator.net/
e.prime = function(par,k=120){
  x=par[1]
  y=par[2]
  dy = -((1-x/k)*(x/k+1)*(y/k+1)*((-y-47)*sin(sqrt(abs(y+x/2+47)))-x*sin(sqrt(abs(y-x+47)))))/k+((1-x/k)*(x/k+1)*(1-y/k)*((-y-47)*sin(sqrt(abs(y+x/2+47)))-x*sin(sqrt(abs(y-x+47)))))/k+(1-x/k)*(x/k+1)*(1-y/k)*(y/k+1)*(-sin(sqrt(abs(y+x/2+47)))+((-y-47)*(y+x/2+47)*cos(sqrt(abs(y+x/2+47))))/(2*abs(y+x/2+47)^(3/2))-(x*(y-x+47)*cos(sqrt(abs(y-x+47))))/(2*abs(y-x+47)^(3/2)))
  dx = -((1-y/k)*(y/k+1)*(x/k+1)*((-y-47)*sin(sqrt(abs(x/2+y+47)))-x*sin(sqrt(abs(x-y-47)))))/k+((1-y/k)*(y/k+1)*(1-x/k)*((-y-47)*sin(sqrt(abs(x/2+y+47)))-x*sin(sqrt(abs(x-y-47)))))/k+(1-y/k)*(y/k+1)*(1-x/k)*(x/k+1)*(-sin(sqrt(abs(x-y-47)))-(x*(x-y-47)*cos(sqrt(abs(x-y-47))))/(2*abs(x-y-47)^(3/2))+((-y-47)*cos(sqrt(abs(x/2+y+47)))*(x/2+y+47))/(4*abs(x/2+y+47)^(3/2)))
  return(c(dx,dy))
}

gradient.descent <- function(starting.x, starting.y, gr=e.prime,  gamma = .225, epsilon = .001, n = 5000, verbose = FALSE, graph=FALSE, k = 120, ...) {
  e = function(x, y, k = 120) {
    ifelse (test=((abs(x) > k) | (abs(y) > k)),
            yes=0,
            no=(1 - x/k)*(1 - y/k)*(1 + x/k)*(1 + y/k)*
              (-(y + 47) * sin(sqrt(abs(y + x/2 + 47))) - x*sin(sqrt(abs(x - (y + 47)))))
    )
  }
  e.prime = function(par,k=120){
    x=par[1]
    y=par[2]
    dy = -((1-x/k)*(x/k+1)*(y/k+1)*((-y-47)*sin(sqrt(abs(y+x/2+47)))-x*sin(sqrt(abs(y-x+47)))))/k+((1-x/k)*(x/k+1)*(1-y/k)*((-y-47)*sin(sqrt(abs(y+x/2+47)))-x*sin(sqrt(abs(y-x+47)))))/k+(1-x/k)*(x/k+1)*(1-y/k)*(y/k+1)*(-sin(sqrt(abs(y+x/2+47)))+((-y-47)*(y+x/2+47)*cos(sqrt(abs(y+x/2+47))))/(2*abs(y+x/2+47)^(3/2))-(x*(y-x+47)*cos(sqrt(abs(y-x+47))))/(2*abs(y-x+47)^(3/2)))
    dx = -((1-y/k)*(y/k+1)*(x/k+1)*((-y-47)*sin(sqrt(abs(x/2+y+47)))-x*sin(sqrt(abs(x-y-47)))))/k+((1-y/k)*(y/k+1)*(1-x/k)*((-y-47)*sin(sqrt(abs(x/2+y+47)))-x*sin(sqrt(abs(x-y-47)))))/k+(1-y/k)*(y/k+1)*(1-x/k)*(x/k+1)*(-sin(sqrt(abs(x-y-47)))-(x*(x-y-47)*cos(sqrt(abs(x-y-47))))/(2*abs(x-y-47)^(3/2))+((-y-47)*cos(sqrt(abs(x/2+y+47)))*(x/2+y+47))/(4*abs(x/2+y+47)^(3/2)))
    return(c(dx,dy))
  }
  if(graph){
    x=seq(-k, k)
    y=seq(-k, k)
    z=outer(x,y,e)
    persp.out=persp3d(x, y, z, aspect = c(1, 1, 1), col = "grey", xlab = "X", ylab = "Y", zlab = "Z", polygon_offset = 1)
  }
  par = c(starting.x,starting.y)
  z.old = e(par[1],par[2])
  for (i in seq_len(n)) {
    gradient = gr(par, ...)
    par = par + gamma * gradient
    z = e(par[1], par[2])
    gradient.size = sum(abs(gradient))
    if (verbose) {
      cat(sep = "", "i = ", i, ", par = c(", paste(signif(par, 4), collapse = ","),
          "), gradient = c(", paste(signif(gradient, 4), collapse = ","),
          "), size = ", signif(gradient.size, 4), "\n")
    }
    if (graph) { # This block is not part of the algorithm. It just draws.
      points3d(x=par[1],y=par[2],z=0,col="blue",pch=16,add=T)
      points3d(x=par[1],y=par[2],z=z,col="red",pch=16,add=T)
      lines3d(x=c(par[1],par[1]),y=c(par[2],par[2]),z=c(0,z),col="orange",add=T)
    }
    if (gradient.size < epsilon || z.old > z) {
      break
    }
    z.old=z
  }
  return(par)
}

#A local Maximum
peak = gradient.descent(starting.x=0, starting.y = 0, gr = e.prime, gamma = 0.2, n = 1000, verbose = FALSE, graph=TRUE)

#All local Maximums, global maximum, and number of local maximums
peak.finder = function(step.size=10, gradient = e.prime, k = 120, unique = TRUE, display = FALSE){
  all.peaks = matrix(nrow=0,ncol=3)
  #Finding all peaks starts here (slowest part of the function)
  for(x in seq(-k+5,k-5,step.size)){
    for(y in seq(-k+5,k-5,step.size)){
      p = gradient.descent(starting.x = x, starting.y = y, gr=gradient)
      all.peaks = rbind(all.peaks, c(p,e(p[1],p[2])))
    }
  }
  #Finding all peaks ends here
  if(unique){
    all.peaks = round(all.peaks,digit=0)
    all.peaks = unique(all.peaks)
    all.peaks = all.peaks[which(abs(all.peaks[,1])<120 & abs(all.peaks[,2])<120),]
    peaks = matrix(nrow=1,ncol=3)
    peaks[1,]=all.peaks[1,]
    for(i in 2:nrow(all.peaks)){
      close.to.peaks = FALSE
      for(j in 1:nrow(peaks)){
        if(sqrt((all.peaks[i,1]-peaks[j,1])**2+(all.peaks[i,2]-peaks[j,2])**2) < 3){
          close.to.peaks = TRUE
        }
      }
      if(!close.to.peaks){
        peaks = rbind(peaks,all.peaks[i,])
      }
    }
    peaks = peaks[which(peaks[,3]>0),]
  }
  if(display){
    k = 120
    x = seq(-k, k)
    y = seq(-k, k)
    z = outer(x, y, e)
    global.maximum = peaks[which.max(peaks[,3]),]
    persp3d(x, y, z, aspect = c(1, 1, 1), col = "grey", xlab = "X", ylab = "Y", zlab = "Z", polygon_offset = 1)
    points3d(x = global.maximum[1], y = global.maximum[2], z = global.maximum[3], col = "green",size=10,add=T)
    points3d(x = peaks[,1], y = peaks[,2], z = peaks[,3], col = "red",size=10,add=T)
    n = nrow(peaks)
    print(paste0("Number of distinct local maximums is: ",n))
  }
}
peak.finder(display=TRUE)

print("Time elapsed for displaying, removing duplicates peaks, and finding all peaks (total work)")
print(system.time(peak.finder(display=TRUE)))

print("Time elapsed for removing duplicates peaks, and finding all peaks")
print(system.time(peak.finder(display=FALSE)))

print("Time elapsed for finding all peaks is:")
print(system.time(peak.finder(unique = FALSE, display=FALSE)))
#It is apparent that time elapsed finding all the peaks is the bulk of the total time taken

#Method 2 (using mapply at bottleneck)
peak.finder.M2 = function(step.size=10, gradient = e.prime, k = 120, unique = TRUE, display = FALSE){
  starting = matrix(nrow=0,ncol=2)
  for(i in seq(-k+5,k-5,step.size)){
    for(j in seq(-k+5,k-5,step.size)){
      starting = rbind(starting,c(i, j))
    }
  }
  all.peaks = mapply(FUN = gradient.descent, starting[,1], starting[,2])
  all.peaks = t(all.peaks)
  all.peaks = cbind(all.peaks,vector(length = nrow(all.peaks)))
  for(i in 1:nrow(all.peaks)){
    all.peaks[i,3]=e(all.peaks[i,1],all.peaks[i,2])
  }
  if(unique){
    all.peaks = round(all.peaks,digit=0)
    all.peaks = unique(all.peaks)
    all.peaks = all.peaks[which(abs(all.peaks[,1])<120 & abs(all.peaks[,2])<120),]
    peaks = matrix(nrow=1,ncol=3)
    peaks[1,]=all.peaks[1,]
    for(i in 2:nrow(all.peaks)){
      close.to.peaks = FALSE
      for(j in 1:nrow(peaks)){
        if(sqrt((all.peaks[i,1]-peaks[j,1])**2+(all.peaks[i,2]-peaks[j,2])**2) < 3){
          close.to.peaks = TRUE
        }
      }
      if(!close.to.peaks){
        peaks = rbind(peaks,all.peaks[i,])
      }
    }
    peaks = peaks[which(peaks[,3]>0),]
  }
  if(display){
    k = 120
    x = seq(-k, k)
    y = seq(-k, k)
    z = outer(x, y, e)
    global.maximum = peaks[which.max(peaks[,3]),]
    persp3d(x, y, z, aspect = c(1, 1, 1), col = "grey", xlab = "X", ylab = "Y", zlab = "Z", polygon_offset = 1)
    points3d(x = global.maximum[1], y = global.maximum[2], z = global.maximum[3], col = "green",size=10,add=T)
    points3d(x = peaks[,1], y = peaks[,2], z = peaks[,3], col = "red",size=10,add=T)
    n = nrow(peaks)
    print(paste0("Number of distinct local maximums is: ",n))
  }
}
#to verify method 2 works use:
peak.finder.M2(display=TRUE)

print("Time elapsed for finding all peaks for method 2 is:")
print(system.time(peak.finder.M2(unique = FALSE, display=FALSE)))

#Method 3 (using parallel computing at bottleneck)
peak.finder.M3 = function(step.size=10, gradient = e.prime, k = 120, unique = TRUE, display = FALSE){
  starting = matrix(nrow=0,ncol=2)
  for(i in seq(-k+5,k-5,step.size)){
    for(j in seq(-k+5,k-5,step.size)){
      starting = rbind(starting,c(i, j))
    }
  }
  
  n.cores = detectCores()
  if (.Platform$OS.type == "windows") {
    cluster = makePSOCKcluster(names = n.cores)
    all.peaks = clusterMap(cl = cluster, fun = gradient.descent, starting[,1], starting[,2])
    stopCluster(cluster)
  } else {
    all.peaks =  mcmapply(FUN = gradient.descent, starting[,1], starting[,2], mc.cores = n.cores)
  }
  temp = matrix(nrow = length(all.peaks),ncol=3)
  for(i in 1:length(all.peaks)){
    temp[i,]=c(all.peaks[[i]][1],all.peaks[[i]][2],e(all.peaks[[i]][1],all.peaks[[i]][2]))
  }
  all.peaks = temp
  if(unique){
    all.peaks = round(all.peaks,digit=0)
    all.peaks = unique(all.peaks)
    all.peaks = all.peaks[which(abs(all.peaks[,1])<120 & abs(all.peaks[,2])<120),]
    peaks = matrix(nrow=1,ncol=3)
    peaks[1,]=all.peaks[1,]
    for(i in 2:nrow(all.peaks)){
      close.to.peaks = FALSE
      for(j in 1:nrow(peaks)){
        if(sqrt((all.peaks[i,1]-peaks[j,1])**2+(all.peaks[i,2]-peaks[j,2])**2) < 3){
          close.to.peaks = TRUE
        }
      }
      if(!close.to.peaks){
        peaks = rbind(peaks,all.peaks[i,])
      }
    }
    peaks = peaks[which(peaks[,3]>0),]
  }
  if(display){
    k = 120
    x = seq(-k, k)
    y = seq(-k, k)
    z = outer(x, y, e)
    global.maximum = peaks[which.max(peaks[,3]),]
    persp3d(x, y, z, aspect = c(1, 1, 1), col = "grey", xlab = "X", ylab = "Y", zlab = "Z", polygon_offset = 1)
    points3d(x = global.maximum[1], y = global.maximum[2], z = global.maximum[3], col = "green",size=10,add=T)
    points3d(x = peaks[,1], y = peaks[,2], z = peaks[,3], col = "red",size=10,add=T)
    n = nrow(peaks)
    print(paste0("Number of distinct local maximums is: ",n))
  }
}
#to verify method 3 works use:
peak.finder.M3(display=TRUE)

print("Time elapsed for finding all peaks for method 3 is:")
print(system.time(peak.finder.M3(unique = FALSE, display=FALSE)))

cppFunction("
#include <Rcpp.h>
#include <cstdlib>
#include <cmath>
#include <list>
#include<math.h> 
using namespace Rcpp;

double e(double x, double y, int k = 120){
  return (1 - x/k)*(1 - y/k)*(1 + x/k)*(1 + y/k)*(-(y + 47) * sin(sqrt(abs(y + x/2 + 47))) - x*sin(sqrt(abs(x - (y + 47)))));
}

std::vector<double> e_prime(std::vector<double> par, int k = 120){
  double x = par[0];
  double y = par[1];
  double dy = -((1-x/k)*(x/k+1)*(y/k+1)*((-y-47)*sin(sqrt(abs(y+x/2+47)))-x*sin(sqrt(abs(y-x+47)))))/k+((1-x/k)*(x/k+1)*(1-y/k)*((-y-47)*sin(sqrt(abs(y+x/2+47)))-x*sin(sqrt(abs(y-x+47)))))/k+(1-x/k)*(x/k+1)*(1-y/k)*(y/k+1)*(-sin(sqrt(abs(y+x/2+47)))+((-y-47)*(y+x/2+47)*cos(sqrt(abs(y+x/2+47))))/(2*pow(abs(y+x/2+47),1.5))-(x*(y-x+47)*cos(sqrt(abs(y-x+47))))/(2*pow(abs(y-x+47),1.5)));
  double dx = -((1-y/k)*(y/k+1)*(x/k+1)*((-y-47)*sin(sqrt(abs(x/2+y+47)))-x*sin(sqrt(abs(x-y-47)))))/k+((1-y/k)*(y/k+1)*(1-x/k)*((-y-47)*sin(sqrt(abs(x/2+y+47)))-x*sin(sqrt(abs(x-y-47)))))/k+(1-y/k)*(y/k+1)*(1-x/k)*(x/k+1)*(-sin(sqrt(abs(x-y-47)))-(x*(x-y-47)*cos(sqrt(abs(x-y-47))))/(2*pow(abs(x-y-47),1.5))+((-y-47)*cos(sqrt(abs(x/2+y+47)))*(x/2+y+47))/(4*pow(abs(x/2+y+47),1.5)));
  std::vector<double> grad;
  grad.push_back(dx);
  grad.push_back(dy);
  return grad;
}

double sum_abs(std::vector<double> gradient){
  return sqrt(pow(gradient[0],2)) + sqrt(pow(gradient[1],2));
}

std::vector<double> gradientDescent(double starting_x, double starting_y, double gamma = 0.225, double epsilon = 0.001, int n = 5000) {
  double z_old = e(starting_x,starting_y);
  std::vector<double> par;
  par.push_back(starting_x);
  par.push_back(starting_y);
  for(int i =0; i <5000; i++){
    std::vector<double> gradient = e_prime(par);
    par[0] = par[0] + gamma * gradient[0];
    par[1] = par[1] + gamma * gradient[1];
    double z = e(par[0],par[1]);
    if(sum_abs(gradient) < epsilon || z_old > z){
      break;
    }
    z_old = z;
  }
  return par;
}
")

#Method 4 (using c++ at bottleneck)
peak.finder.M4 = function(step.size=10, gradient = e.prime, k = 120, unique = TRUE, display = FALSE){
  starting = matrix(nrow=0,ncol=2)
  for(i in seq(-k+5,k-5,step.size)){
    for(j in seq(-k+5,k-5,step.size)){
      starting = rbind(starting,c(i, j))
    }
  }
  #the c++ gradient function does not have graph or verbose ability
  #Note for myself: I can also call the external cpp function using sourceCpp("GradientDescentOptimizationHelper.cpp")
  all.peaks = mapply(FUN = gradientDescent, starting[,1], starting[,2])
  all.peaks = t(all.peaks)
  all.peaks = cbind(all.peaks,vector(length = nrow(all.peaks)))
  for(i in 1:nrow(all.peaks)){
    all.peaks[i,3]=e(all.peaks[i,1],all.peaks[i,2])
  }
  if(unique){
    all.peaks = round(all.peaks,digit=0)
    all.peaks = unique(all.peaks)
    all.peaks = all.peaks[which(abs(all.peaks[,1])<120 & abs(all.peaks[,2])<120),]
    peaks = matrix(nrow=1,ncol=3)
    peaks[1,]=all.peaks[1,]
    for(i in 2:nrow(all.peaks)){
      close.to.peaks = FALSE
      for(j in 1:nrow(peaks)){
        if(sqrt((all.peaks[i,1]-peaks[j,1])**2+(all.peaks[i,2]-peaks[j,2])**2) < 3){
          close.to.peaks = TRUE
        }
      }
      if(!close.to.peaks){
        peaks = rbind(peaks,all.peaks[i,])
      }
    }
    peaks = peaks[which(peaks[,3]>0),]
  }
  if(display){
    k = 120
    x = seq(-k, k)
    y = seq(-k, k)
    z = outer(x, y, e)
    global.maximum = peaks[which.max(peaks[,3]),]
    persp3d(x, y, z, aspect = c(1, 1, 1), col = "grey", xlab = "X", ylab = "Y", zlab = "Z", polygon_offset = 1)
    points3d(x = global.maximum[1], y = global.maximum[2], z = global.maximum[3], col = "green",size=10,add=T)
    points3d(x = peaks[,1], y = peaks[,2], z = peaks[,3], col = "red",size=10,add=T)
    n = nrow(peaks)
    print(paste0("Number of distinct local maximums is: ",n))
  }
}
#to verify method 4 works use:
peak.finder.M4(display=TRUE)

print("Time elapsed for finding all peaks for method 4 is:")
print(system.time(peak.finder.M4(unique = FALSE, display=FALSE)))
#Method 4 has best improvement