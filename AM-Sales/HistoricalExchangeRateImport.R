#Note: Other rates can be found on this website too: https://www.exchangerates.org.uk/AED-IRR-spot-exchange-rates-history-2012.html

#initializes
rm(list=ls())
load("cleaned.sales.Rda")

exchange.rate = data.frame(matrix(ncol = 3, nrow = 0))
colnames(exchange.rate)=c("Year","Month","Rate")

years = seq(2012,2022)
years = as.character(years)
for(year in years){
  lines = readLines(paste0("https://www.exchangerates.org.uk/AED-IRR-spot-exchange-rates-history-",year,".html"))
  months.table.index = grep(pattern="\\s\\s\\d{1,2}\\s\\w*\\s\\d\\d\\d\\d",x=lines)
  for(month.table.index in months.table.index){
    month.table = readHTMLTable(lines[(month.table.index-6):(month.table.index+2)])
    last.row = length(month.table$hist[,1])
    row.to.keep = month.table$hist[last.row,1]
    curr.row = nrow(exchange.rate)+1
    exchange.rate[curr.row,"Rate"] = as.numeric(sub(pattern=".*(\\s\\d+.\\d*$)",x=row.to.keep,replacement="\\1"))
    exchange.rate[curr.row,"Month"] = sub(pattern=".*\\s(\\w*)(\\s\\d*).(\\s\\d*\\.\\d*)$",x=row.to.keep,replacement="\\1")
    exchange.rate[curr.row,"Year"] = as.numeric(year)
  }
}

#Exclude this current month due to it html not being like the others
exchange.rate = exchange.rate[-nrow(exchange.rate),]
#Exclude January to June of 2012
exchange.rate = exchange.rate[-seq(1,6),]
rownames(exchange.rate)=NULL

save(exchange.rate,file="exchangeRate.Rda")

#plot(ts(exchange.rate[,3],frequency=12,start=c(1,2012)),main="AED to Iranian Rials exchange rate",ylab="rate")
