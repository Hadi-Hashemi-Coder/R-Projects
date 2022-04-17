#initializes
rm(list=ls())
load("cleaned.sales.Rda")
load("exchangeRate.Rda")

#Creates vector of sales for a given brand at a branch
#For missing data points, an average of sales before and after months are taken
prep = function(brand,branch,data=cleaned.sales){
  sales = c()
  for(month in unique(data$Date)){
    for(row in 1:nrow(data)){
      if(data[row,"Date"] == month & data[row,"Brand"] == brand){
        sales = c(sales, data[row,branch])
      }
    }
  }
  sales=as.numeric(sales)
  
  for(m in 1:(length(sales)-1)){
    if(is.na(sales[m])){
      if(is.na(sales[m+1])){
        sales[m]=mean(c(sales[m-1],sales[m+2]))
      } else {
        sales[m]=mean(c(sales[m-1],sales[m+1]))
      }
    }
  }
  return(sales)
}

#Exponential smoothing (of past month's sales)on Fusen gum sales in Dubai

exp.smoothing = function(brand,branch,data=cleaned.sales,seasons = FALSE){
  sales = prep(brand=brand,branch=branch,data=data)
  ts.sales = ts(data=sales,frequency=12,start=c(2012,7))
  exp.s=HoltWinters(ts.sales, beta = FALSE, gamma = seasons)
  plot(exp.s,main=paste0("Exponential Smoothing for ",brand," in",branch),ylab=paste0("Sales for ",brand))
  return(exp.s)
}

exp.s.s = exp.smoothing(brand="Fusen Gum",branch="Dubai",seasons=TRUE)
exp.s.ns = exp.smoothing(brand="Fusen Gum",branch="Dubai",seasons=FALSE)

exp.s.s$SSE
exp.s.ns$SSE


#Exponential smoothing (of past month's sales) and weighted addition of current month's exchange rate on Fusen gum sales in Genaveh
#You need to have a variable brand and branch before calling this

mean.forecast.error = function(par){
  B0=par[1]
  B1=par[2]
  sum.forecast.error = 0
  forecasts = c(data[1])
  for(month in 2:length(data)){
    Y = B0*data[month-1]+(1-B0)*forecasts[month-1]+B1*exchange.rate[month,3]
    forecasts = c(forecasts, Y)
    sum.forecast.error = sum.forecast.error + (data[month]-forecasts[month])**2
  }
  return((1/length(data))*sum.forecast.error)
}

optimizer = function(brand, branch, par){
  data=prep(brand=brand,branch=branch,data=cleaned.sales)
  B.optimal=optim(par=c(par),mean.forecast.error)
  return(B.optimal)
}

#At the time of making this I didn't have Jan, Feb, March, and April of 2022 for cleaned.sales
exchange.rate = exchange.rate[-c(115,116,117),]
B.optimal = optimizer(brand = "Fusen Gum",branch = "Genaveh", par=c(0.5,0.5))

plotter = function(B,brand,branch){
  B0=B[1]
  B1=B[2]
  forecasts = c(data[1])
  for(month in 2:length(data)){
    Y = B0*data[month-1]+(1-B0)*forecasts[month-1]+B1*exchange.rate[month,3]
    forecasts = c(forecasts, Y)
  }
  plot(ts(data,frequency=12,start=c(2012,7)),main = paste0("Exponential Smoothing with weighted exchange rate addition for ",brand," at  ",branch),ylab="sales")
  lines(ts(forecasts,frequency=12,start=c(2012,7)),col="red")
}

plotter(B.optimal$par,"Fusen Gum","Genaveh")
