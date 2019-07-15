library(ggplot2)
cbPalette <-c("slateblue","#009E73","#D55E00")
#two state
data=read.csv("spread01_100shares_twostate.csv")
#data=data[data$method!="Naive",]

state=paste(data$state1,data$state2,data$state3,sep=",")
state=paste("(",state,")")
data$state=state

ggplot(data,aes(x=state,y=gap_rate,group=method,color=method))+geom_line()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_point(mapping = aes(x = state, y = gap_rate))+
  ylim(-0.3,0.3)+scale_colour_manual(values=cbPalette)

#three state
data1=read.csv("spread01_100shares_threestate.csv")
#data1=data1[data1$method!="Naive",]

state=paste(data1$state1,data1$state2,data1$state3,sep=",")
state=paste("(",state,")")
data1$state=state

ggplot(data1,aes(x=state,y=gap_rate,group=method,color=method))+geom_line()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1))+
  geom_point(mapping = aes(x = state, y = gap_rate))+
  ylim(-0.3,0.3)+scale_colour_manual(values=cbPalette)

