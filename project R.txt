
#install packages
#install.packages("BiocManager")
#BiocManager::install("GEOquery")

#call packages
library(GEOquery)


#get data
Gset<- getGEO('GSE180594', destdir="." , GSEMatrix = TRUE ,AnnotGPL = T)
Gset <- Gset[[1]]
Data<-exprs(Gset)
myfinaldata<-head(Data[ ,1:20])
print(myfinaldata)
View((myfinaldata))
print(Gset)


#classification of data to alive and dead
colnames1<-c("G1","G2","G3","G4","G5","G6","G7","G9","G10","G11","G12","G13","G14","G15","G16","G18","G19")
colnames2<-c("G8","G17","G20")
DataColnames<-c("G1","G2","G3","G4","G5","G6","G7","G8","G9","G10","G11","G12","G13","G14","G15","G16","G17","G18","G19","G20")
rownames<-c(1:249)
data<-matrix(c(Data[,1:20]),nrow = 249,ncol = 20,byrow = F,list(rownames,DataColnames))
alive<-matrix(c(data[,c(-8,-17,-20)]),nrow = 249,ncol = 17,byrow = F,list(rownames,colnames1))
dead<-matrix(c(data[,c(8,17,20)]),nrow = 249,ncol = 3,byrow = F,list(rownames,colnames2))
print(alive)
print(dead)



#t test
Pvalue<-c()
Tvalue<-c()
GenIndex<-c()
counter=1
i=1
while(counter<=20){
  test<-t.test(alive[i,],dead[i,],var.equal = TRUE)
  if((!is.nan(test$p.value)) & (test$p.value<=0.05)){
    Pvalue<-c(Pvalue,test$p.value)
    Tvalue<-c(Tvalue,test$statistic)
    GenIndex<-c(GenIndex,i)
    counter=counter+1
  }
  i=i+1
}
DataFram<-data.frame(GenIndex=GenIndex,T.value=Tvalue,P.value=Pvalue)
print(DataFram)


#visualization
setsname<-c("dead","alive")
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))

color=c("green","green","green","green","green","green","green","yellow","green","green","green","green","green","green","green","green","yellow","green","green","yellow")
status=c("alive","alive","alive","alive","alive","alive","alive","dead","alive","alive","alive","alive","alive","alive","alive","alive","dead","alive","alive","dead")


#boxplot of all data before filtering
boxplot(Data,boxwex=0.7,xlab = "Name of Gene", ylab = "Mileage Data",main="GSE180594 before filtering",notch=F,outline=FALSE,las=1,col=palette())

# boxplot of 20 gene after filtering
datdata <- data.frame(dead,alive)
boxplot( datdata,xlab = "Name of Gene", ylab = "Mileage Data",main="alive and dead boxplot",col = color)

# boxplot of alive status
boxplot(alive, notch=F,varwidth=T,xlab = "Name of Gene", ylab = "Mileage Data",main="alive boxplot", outline=TRUE,las=2,col="#8EEC9C")

# boxplot of dead status
boxplot(dead, notch=F,varwidth=T,xlab = "Name of Gene", ylab = "Mileage Data",main="dead boxplot", outline=TRUE,las=2,col="#F08080")

#scatter of 20 gene data after filtering
plot(x=data[GenIndex,],xlab="samples",main = "GSE180594")

#scatter of two variable (GenIndex , Pvalue)
attach(DataFram)
plot(GenIndex,Pvalue , main="Scatterplot Example1",
     xlab="GenIndex ", ylab="Pvalue ", pch=19)

DataFrama<-data.frame(alive)
print(DataFrama)
#scatter of two variable (G1 , G3)
attach(DataFrama)
plot(G1,G3 , main="Scatterplot Example2",
     xlab="G1 ", ylab="G3 ", pch=19)

#heatmap of 20 gene data after filteringS
heatmap(cbind(dead, alive),col=color, main="heatmap of 20 gene")

#heatmap of alive status
heatmap(alive, main="heatmap of alive status")

#heatmap of head status
heatmap(dead, main="heatmap of dead status")

# pie chart of 20 gene data
pie(GenIndex, labels = status, main = "GSE180594 pie chart",col = color)

# Bar Charts of 20 gene data 
M <- c("G1","G2","G3","G4","G5","G6","G7","G8","G9","G10","G11","G12","G13","G14","G15","G16","G17","G18","G19","G20")
barplot(GenIndex,names.arg=M,xlab="name of Gene",ylab = "Mileage Data",col=color, main="GSE180594 chart", border="red")

# histogram of 20 gene data 
hist(GenIndex,xlab = "name of Gene",col = "green",border = "red", xlim = c(0,20), ylim = c(0,3),  breaks = 10)

#line graph of GSE180594 
plot(GenIndex,type = "o", col = color, xlab = "name of Gene", ylab = "classification", main = "line graph of GSE180594 ")



