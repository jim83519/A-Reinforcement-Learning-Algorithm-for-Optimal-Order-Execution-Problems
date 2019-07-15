#Given index
n=2
delta_t=1/n
rho=6
index=exp(-rho*delta_t)


da=read.csv("spread01_100shares_threestate.csv")
da$spread2=index*da$spread1+index*(1/da$state1)*da$first
da$spread3=index*da$spread2+index*(1/da$state2)*da$second
da$cost=da$first*(da$spread1+da$first/(2*da$state1))+da$second*(da$spread2+da$second/(2*da$state2))+
        da$third*(da$spread3+da$third/(2*da$state3))

da$method=c(rep("Naive",27),rep("Partition",27),rep("Q-learning",27))
da$gap_rate=(da$cost-da$cost[da$method=="Naive"])/da$cost[da$method=="Naive"]

#write.csv(da,"spread01_100shares_threestate.csv",row.names = F)

naive=da$cost[da$method=="Naive"]
partition=da$cost[da$method=="Partition"]
q_learning=da$cost[da$method=="Q-learning"]

which(naive<partition)
which(naive<q_learning)

test=cbind(naive,partition)
test[which(naive<partition),]
