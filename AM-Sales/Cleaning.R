#Notes: if columns are added then this has to be reviewed

#initializes
rm(list=ls())
setwd("/cloud/project/Sales")

Sales=read_excel("AlMahasenSales.xlsx")
#ignoring unwanted columns
Sales = Sales[,-c(11)]

#Collecting Headers
headers = c()
for(i in 1:nrow(Sales)){
  if(!is.na(Sales[i,1])){
    for(j in 2:length(Sales[i,])){
      if(!(Sales[i,j] %in% headers) & !is.na(Sales[i,j])){
        headers=c(headers, Sales[i,j])
      }
    }
  }
}
headers=c(headers)

#Creating new empty data frame
names = c("Date","Brand",headers[-c(4,5,8)],headers[c(5,8)])
cleaned.sales = data.frame(matrix(ncol = 11, nrow = 0))
colnames(cleaned.sales)=names

#Entering Values
original.headers.rows = grep(pattern="*\\w*,\\s\\d*",x=unlist(Sales[,1]),ignore.case = TRUE)
original.headers.rows = c(original.headers.rows, nrow(Sales)+1)
for(header.row.index in 1:(length(original.headers.rows)-1)){
  curr.H = original.headers.rows[header.row.index]
  next.H = original.headers.rows[header.row.index+1]
  curr.row = nrow(cleaned.sales)
  cleaned.sales[(curr.row+1):(curr.row+next.H-curr.H-1),"Date"] = Sales[curr.H,1]
  cleaned.sales[(curr.row+1):(curr.row+next.H-curr.H-1),"Brand"] = Sales[(curr.H+1):(next.H-1),9]
  cleaned.sales[(curr.row+1):(curr.row+next.H-curr.H-1),"Approximate profit"] = Sales[(curr.H+1):(next.H-1),10]
  for(colnum in 1:length(Sales[curr.H,])){
    if(Sales[curr.H,colnum] %in% headers[-c(4,8)]){
      cleaned.sales[(curr.row+1):(curr.row+next.H-curr.H-1),unlist(Sales[curr.H,colnum])] = Sales[(curr.H+1):(next.H-1),colnum]
    }
  }
}

#Remove dashes
for(Col in colnames(cleaned.sales)){
  cleaned.sales[,Col]=gsub(pattern="-{3,}",x=cleaned.sales[,Col],replacement="")
}

#July, 2017 & August, 2017 is missing so replacing it with empty values for now
July = data.frame(matrix(ncol=11,nrow=16))
July[,2]=cleaned.sales[cleaned.sales[,1]=="June, 2017",2]
July[,1]="July, 2017"
colnames(July)=names
August = data.frame(matrix(ncol=11,nrow=16))
August[,2]=cleaned.sales[cleaned.sales[,1]=="June, 2017",2]
August[,1]="August, 2017"
colnames(August)=names
cleaned.sales = rbind(cleaned.sales[1:890,],July,August,cleaned.sales[891:nrow(cleaned.sales),])

#February, 2016 is misspelled
cleaned.sales[cleaned.sales[,1]=="Febuary, 2016",1] = "February, 2016"

View(cleaned.sales)
save(cleaned.sales,file="cleaned.sales.Rda")
#can load with load("cleaned.sales.Rda")