library(data.table)
setwd("C:/Users/asus/Desktop/碩論/data")

data0519_mes=fread("DIS_2017-05-19_34200000_57600000_message_30.csv",header=F)
data0519_lob=fread("DIS_2017-05-19_34200000_57600000_orderbook_30.csv",header=F)
colnames(data0519_mes)=c("Time(sec)","Event Type","Order ID","Size","Price","Direction")
colnames(data0519_lob)=rep(c("Ask price","Volume","Bid price","Volume"),30)
data0519=cbind(data0519_mes,data0519_lob)
ask_0519=data0519[data0519$Direction==-1,]
ask_0519=as.data.frame(ask_0519)
ask_0519$spread=ask_0519[,7]-ask_0519[,9]

#取出十點到十二點的資料
da=ask_0519[ask_0519$`Time(sec)`>36000 & ask_0519$`Time(sec)`<43200,]
da$index=NA
index=seq(36000,43200,30)
for(i in c(1:(length(index)-1))){
    da[da$`Time(sec)`>=index[i] & da$`Time(sec)`<index[i+1],]$index=i=i
}

#算ask price 前五的量平均
ask_30=aggregate(da[,c(8,12,16,20,24,127)],by = list(index=da$index), mean)
ask_30$Volume_five=apply(ask_30[,c(2:6)],MARGIN=1,sum)
hist(ask_30$Volume_five)
quantile(ask_30$Volume_five,prob=c(0,1/3,2/3,1))
quantile(ask_30$Volume_five,prob=c(0,0.25,0.75,1))

#找出有執行的量平均
da_execution=da[da$`Event Type`==4,]
execution_30=aggregate(da_execution$Size,by=list(index=da_execution$index),mean)
da_index=data.frame(index=seq(1,length(index)))
execution_30=merge(da_index,execution_30,all.x = T)
execution_30[is.na(execution_30$x),2]=0
execution_30=execution_30[-241,]
names(execution_30)[2]="average_execution"

#把兩筆資料合併
da_30=cbind(ask_30,execution_30)
da_low=da_30[da_30$Volume_five<quantile(ask_30$Volume_five,prob=c(0,1/3,2/3,1))[2],]
da_mid=da_30[da_30$Volume_five>=quantile(ask_30$Volume_five,prob=c(0,1/3,2/3,1))[2] &
             da_30$Volume_five<quantile(ask_30$Volume_five,prob=c(0,1/3,2/3,1))[3],]
da_high=da_30[da_30$Volume_five>=quantile(ask_30$Volume_five,prob=c(0,1/3,2/3,1))[3],]

da_low$state=round(mean(da_low$Volume_five),digits=-2)
da_mid$state=round(mean(da_mid$Volume_five),digits=-2)
da_high$state=round(mean(da_high$Volume_five),digits=-2)
da_30=rbind(da_low,da_mid,da_high)
da_30=da_30[order(da_30$index),]

da_30$spread=da_30$spread/10000
#估計rho

da_30$spread_next=c(da_30$spread[-1],NA)
da_30=da_30[-240,]
da_30$x=da_30$spread+(da_30$average_execution/da_30$state)
model=lm(log(spread_next)~log(x),data=da_30)
summary(model)

table(da_30$state)
par(mfrow=c(1,1))
plot(log(da_30$spread_next)~log(da_30$x),pch=16,xlab="log(x)",ylab="log(spread_next)")
abline(coefficients(model)[1],coefficients(model)[2],col=2)
# plot(model)
# return()

par(mfrow=c(1,1))
hist(da_30$spread,breaks = 8, main="The frequency of spread",xlab="spread")
summary(da_30$spread)

