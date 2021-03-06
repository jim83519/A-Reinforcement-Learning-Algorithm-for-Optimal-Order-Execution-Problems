coef=function(coef1,coef2){
  ma=matrix(rep(NA,9),ncol=3)
  ma[1,1]=coef1[1]*coef2[1]
  ma[1,2]=coef1[1]*coef2[2]
  ma[1,3]=coef1[1]*coef2[3]
  ma[2,1]=coef1[2]*coef2[1]
  ma[2,2]=coef1[2]*coef2[2]
  ma[2,3]=coef1[2]*coef2[3]
  ma[3,1]=coef1[3]*coef2[1]
  ma[3,2]=coef1[3]*coef2[2]
  ma[3,3]=coef1[3]*coef2[3]
  return(ma)
}
coef_sq=function(coef1,coef2){
  ma=matrix(rep(NA,4),ncol=2)
  ma[1,1]=coef1[1]*coef2[1]
  ma[1,2]=coef1[1]*coef2[2]
  ma[2,1]=coef1[2]*coef2[1]
  ma[2,2]=coef1[2]*coef2[2]
  return(ma)
}

################
s=c(100,200) #status
pr=matrix(rep(0.5,4),nrow=2,byrow=T) #transition probabilities (p11,p12,p21,p22)

# theta=5
# 
# for (i in 1:2){
#   for(j in 1:2){
#     pr[i,j]=exp(-(s[j]-s[i]-theta*(u[]-s[i])*delta_t)^2/(2*sigma^2*s[i]^2*delta_t)) #u is affected by which time step
#   }
# }
# pr[1,1]
# 
# 
n=2
delta_t=1/n
rho=6
index=exp(-rho*delta_t)
#index=1/2
################

matA=matrix(c(1,0,0,index),nrow = 2,ncol = 2,byrow =TRUE)
matF1=t(matrix(c(0,0),nrow=1,ncol=2,byrow=TRUE))
matF2=t(matrix(c(1,0),nrow=1,ncol=2,byrow=TRUE))
matG=t(matrix(c(1,0),nrow=1))

j=1  
matB_s1=t(matrix(c(-1,index*(1/s[j])),nrow=1,ncol=2,byrow=TRUE))###
matC_s1=matrix(c(1/(2*s[j]),0,1/2,0,0,0,1/2,0,0),nrow=3,ncol=3,byrow=TRUE)

temp1_s1=c(matB_s1[1,1],matA[1,1],matA[1,2]) #u,x,d coef (x_(i+1)=xi-ui)
temp2_s1=c(matB_s1[2,1],matA[2,1],matA[2,2]) #u,x,d coef (d_(i+1)=0.5*di+0.5*ui)

j=2
matB_s2=t(matrix(c(-1,index*(1/s[j])),nrow=1,ncol=2,byrow=TRUE))###
matC_s2=matrix(c(1/(2*s[j]),0,1/2,0,0,0,1/2,0,0),nrow=3,ncol=3,byrow=TRUE)

temp1_s2=c(matB_s2[1,1],matA[1,1],matA[1,2])  #u,x,d coef (x_(i+1)=xi-ui)
temp2_s2=c(matB_s2[2,1],matA[2,1],matA[2,2])  #u,x,d coef (d_(i+1)=0.5*di+0.25*ui)

#i=2
#####
#i=2,s1
cost_s1_2=matC_s1
cost_s1_2
#i=2,s2
cost_s2_2=matC_s2
cost_s2_2

#s1個別的cost function(不同的u,在限制下)(i=2)
  u_s1_2=matrix(rep(NA,6),ncol=2)
  a_s1_2=cost_s1_2[1,1]
  b_s1_2=c(cost_s1_2[1,2]+cost_s1_2[2,1],cost_s1_2[1,3]+cost_s1_2[3,1])
  #u<=Fx
  u_s1_2[1,]=t(matF2)
  temp1=cost_s1_2[1,1]*coef_sq(u_s1_2[1,],u_s1_2[1,])+(cost_s1_2[1,3]+cost_s1_2[3,1])*coef_sq(u_s1_2[1,],c(0,1))
  #Fx<u<Gx
  u_s1_2[2,]=-b_s1_2/(2*a_s1_2)
  temp2=cost_s1_2[1,1]*coef_sq(u_s1_2[1,],u_s1_2[1,])+(cost_s1_2[1,3]+cost_s1_2[3,1])*coef_sq(u_s1_2[1,],c(0,1))
  #u>=Gx
  u_s1_2[3,]=t(matG)
  temp3=cost_s1_2[1,1]*coef_sq(u_s1_2[1,],u_s1_2[1,])+(cost_s1_2[1,3]+cost_s1_2[3,1])*coef_sq(u_s1_2[1,],c(0,1))
  
  value_s1_2=cbind(temp1,temp2,temp3) #value function (x,d)

#s2個別的cost function(不同的u,在限制下)(i=2)
  u_s2_2=matrix(rep(NA,6),ncol=2)
  a_s2_2=cost_s2_2[1,1]
  b_s2_2=c(cost_s2_2[1,2]+cost_s2_2[2,1],cost_s2_2[1,3]+cost_s2_2[3,1])
  #u<=Fx
  u_s2_2[1,]=t(matF2)
  temp1=cost_s2_2[1,1]*coef_sq(u_s2_2[1,],u_s2_2[1,])+(cost_s2_2[1,3]+cost_s2_2[3,1])*coef_sq(u_s2_2[1,],c(0,1))
  #Fx<u<Gx
  u_s2_2[2,]=-b_s2_2/(2*a_s2_2)
  temp2=cost_s2_2[1,1]*coef_sq(u_s2_2[1,],u_s2_2[1,])+(cost_s2_2[1,3]+cost_s2_2[3,1])*coef_sq(u_s2_2[1,],c(0,1))
  #u>=Gx
  u_s2_2[3,]=t(matG)
  temp3=cost_s2_2[1,1]*coef_sq(u_s2_2[1,],u_s2_2[1,])+(cost_s2_2[1,3]+cost_s2_2[3,1])*coef_sq(u_s2_2[1,],c(0,1))
  
  value_s2_2=cbind(temp1,temp2,temp3)

  
#here since the cost function for region is the smae 
#i=1,s1
cost_s1_1=matC_s1+
          pr[1,1]*(value_s1_2[1,1]*coef(temp1_s1,temp1_s1)+(value_s1_2[1,2]+value_s1_2[2,1])*coef(temp1_s1,temp2_s1))+
          pr[1,2]*(value_s2_2[1,1]*coef(temp1_s1,temp1_s1)+(value_s2_2[1,2]+value_s2_2[2,1])*coef(temp1_s1,temp2_s1))
cost_s1_1
#i=1,s2
cost_s2_1=matC_s2+
          pr[2,1]*(value_s1_2[1,1]*coef(temp1_s2,temp1_s2)+(value_s1_2[1,2]+value_s1_2[2,1])*coef(temp1_s2,temp2_s2))+
          pr[2,2]*(value_s2_2[1,1]*coef(temp1_s2,temp1_s2)+(value_s2_2[1,2]+value_s2_2[2,1])*coef(temp1_s2,temp2_s2))
cost_s2_1

#s1個別的cost function(不同的u,在限制下)(i=1)
  u_s1_1=matrix(rep(NA,6),ncol=2)
  a_s1_1=cost_s1_1[1,1]
  b_s1_1=c(cost_s1_1[1,2]+cost_s1_1[2,1],cost_s1_1[1,3]+cost_s1_1[3,1])
  #u<Fx
  u_s1_1[1,]=t(matF1)
  temp1=cost_s1_1[1,1]*coef_sq(u_s1_1[1,],u_s1_1[1,])+(cost_s1_1[1,2]+cost_s1_1[2,1])*coef_sq(u_s1_1[1,],c(1,0))+
        (cost_s1_1[1,3]+cost_s1_1[3,1])*coef_sq(u_s1_1[1,],c(0,1))+cost_s1_1[2:3,2:3]
  #Fx<u<Gx
  u_s1_1[2,]=-b_s1_1/(2*a_s1_1)
  temp2=cost_s1_1[1,1]*coef_sq(u_s1_1[2,],u_s1_1[2,])+(cost_s1_1[1,2]+cost_s1_1[2,1])*coef_sq(u_s1_1[2,],c(1,0))+
        (cost_s1_1[1,3]+cost_s1_1[3,1])*coef_sq(u_s1_1[2,],c(0,1))+cost_s1_1[2:3,2:3]
  #u>Gx
  u_s1_1[3,]=t(matG)
  temp3=cost_s1_1[1,1]*coef_sq(u_s1_1[3,],u_s1_1[3,])+(cost_s1_1[1,2]+cost_s1_1[2,1])*coef_sq(u_s1_1[3,],c(1,0))+
        (cost_s1_1[1,3]+cost_s1_1[3,1])*coef_sq(u_s1_1[3,],c(0,1))+cost_s1_1[2:3,2:3]
  
  value_s1_1=cbind(temp1,temp2,temp3)

#s2個別的cost function(不同的u,在限制下)(i=1)
  u_s2_1=matrix(rep(NA,6),ncol=2)
  a_s2_1=cost_s2_1[1,1]
  b_s2_1=c(cost_s2_1[1,2]+cost_s2_1[2,1],cost_s2_1[1,3]+cost_s2_1[3,1])
  #u<Fx
  u_s2_1[1,]=t(matF1)
  temp1=cost_s2_1[1,1]*coef_sq(u_s2_1[1,],u_s2_1[1,])+(cost_s2_1[1,2]+cost_s2_1[2,1])*coef_sq(u_s2_1[1,],c(1,0))+
        (cost_s2_1[1,3]+cost_s2_1[3,1])*coef_sq(u_s2_1[1,],c(0,1))+cost_s2_1[2:3,2:3]
  #Fx<u<Gx
  u_s2_1[2,]=-b_s2_1/(2*a_s2_1)
  temp2=cost_s2_1[1,1]*coef_sq(u_s2_1[2,],u_s2_1[2,])+(cost_s2_1[1,2]+cost_s2_1[2,1])*coef_sq(u_s2_1[2,],c(1,0))+
        (cost_s2_1[1,3]+cost_s2_1[3,1])*coef_sq(u_s2_1[2,],c(0,1))+cost_s2_1[2:3,2:3]
  #u>Gx
  u_s2_1[3,]=t(matG)
  temp3=cost_s2_1[1,1]*coef_sq(u_s2_1[3,],u_s2_1[3,])+(cost_s2_1[1,2]+cost_s2_1[2,1])*coef_sq(u_s2_1[3,],c(1,0))+
        (cost_s2_1[1,3]+cost_s2_1[3,1])*coef_sq(u_s2_1[3,],c(0,1))+cost_s2_1[2:3,2:3]
  
  value_s2_1=cbind(temp1,temp2,temp3)

#####

#i=0,s1
#####
region_A_s1=pr[1,1]*(value_s1_1[1:2,5:6])+pr[1,2]*(value_s2_1[1:2,5:6])
region_B_s1=pr[1,1]*(value_s1_1[1:2,3:4])+pr[1,2]*(value_s2_1[1:2,5:6])
region_C_s1=pr[1,1]*(value_s1_1[1:2,3:4])+pr[1,2]*(value_s2_1[1:2,3:4])
region_D_s1=pr[1,1]*(value_s1_1[1:2,1:2])+pr[1,2]*(value_s2_1[1:2,3:4])
region_E_s1=pr[1,1]*(value_s1_1[1:2,1:2])+pr[1,2]*(value_s2_1[1:2,1:2])

#s1個別的cost function(不同的u,在限制下)(i=0)
costA_s1=matC_s1+region_A_s1[1,1]*coef(temp1_s1,temp1_s1)+(region_A_s1[1,2]+region_A_s1[2,1])*coef(temp1_s1,temp2_s1)+
         region_A_s1[2,2]*coef(temp2_s1,temp2_s1)

  u_s1_0_A=matrix(rep(NA,6),ncol=2)
  a_s1_0=costA_s1[1,1]
  b_s1_0=c(costA_s1[1,2]+costA_s1[2,1],costA_s1[1,3]+costA_s1[3,1])
  #u<Fx
  u_s1_0_A[1,]=t(matF1)
  temp1=costA_s1[1,1]*coef_sq(u_s1_0_A[1,],u_s1_0_A[1,])+(costA_s1[1,2]+costA_s1[2,1])*coef_sq(u_s1_0_A[1,],c(1,0))+
        (costA_s1[1,3]+costA_s1[3,1])*coef_sq(u_s1_0_A[1,],c(0,1))+costA_s1[2:3,2:3]
  #Fx<u<Gx
  u_s1_0_A[2,]=-b_s1_0/(2*a_s1_0)
  temp2=costA_s1[1,1]*coef_sq(u_s1_0_A[2,],u_s1_0_A[2,])+(costA_s1[1,2]+costA_s1[2,1])*coef_sq(u_s1_0_A[2,],c(1,0))+
        (costA_s1[1,3]+costA_s1[3,1])*coef_sq(u_s1_0_A[2,],c(0,1))+costA_s1[2:3,2:3]
  #u>Gx
  u_s1_0_A[3,]=t(matG)
  temp3=costA_s1[1,1]*coef_sq(u_s1_0_A[3,],u_s1_0_A[3,])+(costA_s1[1,2]+costA_s1[2,1])*coef_sq(u_s1_0_A[3,],c(1,0))+
        (costA_s1[1,3]+costA_s1[3,1])*coef_sq(u_s1_0_A[3,],c(0,1))+costA_s1[2:3,2:3]

  value_s1_0_A=cbind(temp1,temp2,temp3)


costB_s1=matC_s1+region_B_s1[1,1]*coef(temp1_s1,temp1_s1)+(region_B_s1[1,2]+region_B_s1[2,1])*coef(temp1_s1,temp2_s1)+
         region_B_s1[2,2]*coef(temp2_s1,temp2_s1)

  u_s1_0_B=matrix(rep(NA,6),ncol=2)
  a_s1_0=costB_s1[1,1]
  b_s1_0=c(costB_s1[1,2]+costB_s1[2,1],costB_s1[1,3]+costB_s1[3,1])
  #u<Fx
  u_s1_0_B[1,]=t(matF1)
  temp1=costB_s1[1,1]*coef_sq(u_s1_0_B[1,],u_s1_0_B[1,])+(costB_s1[1,2]+costB_s1[2,1])*coef_sq(u_s1_0_B[1,],c(1,0))+
        (costB_s1[1,3]+costB_s1[3,1])*coef_sq(u_s1_0_B[1,],c(0,1))+costB_s1[2:3,2:3]
  #Fx<u<Gx
  u_s1_0_B[2,]=-b_s1_0/(2*a_s1_0)
  temp2=costB_s1[1,1]*coef_sq(u_s1_0_B[2,],u_s1_0_B[2,])+(costB_s1[1,2]+costB_s1[2,1])*coef_sq(u_s1_0_B[2,],c(1,0))+
        (costB_s1[1,3]+costB_s1[3,1])*coef_sq(u_s1_0_B[2,],c(0,1))+costB_s1[2:3,2:3]
  #u>Gx
  u_s1_0_B[3,]=t(matG)
  temp3=costB_s1[1,1]*coef_sq(u_s1_0_B[3,],u_s1_0_B[3,])+(costB_s1[1,2]+costB_s1[2,1])*coef_sq(u_s1_0_B[3,],c(1,0))+
        (costB_s1[1,3]+costB_s1[3,1])*coef_sq(u_s1_0_B[3,],c(0,1))+costB_s1[2:3,2:3]

  value_s1_0_B=cbind(temp1,temp2,temp3)

costC_s1=matC_s1+region_C_s1[1,1]*coef(temp1_s1,temp1_s1)+(region_C_s1[1,2]+region_C_s1[2,1])*coef(temp1_s1,temp2_s1)+
         region_C_s1[2,2]*coef(temp2_s1,temp2_s1)

  u_s1_0_C=matrix(rep(NA,6),ncol=2)
  a_s1_0=costC_s1[1,1]
  b_s1_0=c(costC_s1[1,2]+costC_s1[2,1],costC_s1[1,3]+costC_s1[3,1])
  #u<Fx
  u_s1_0_C[1,]=t(matF1)
  temp1=costC_s1[1,1]*coef_sq(u_s1_0_C[1,],u_s1_0_C[1,])+(costC_s1[1,2]+costC_s1[2,1])*coef_sq(u_s1_0_C[1,],c(1,0))+
        (costC_s1[1,3]+costC_s1[3,1])*coef_sq(u_s1_0_C[1,],c(0,1))+costC_s1[2:3,2:3]
  #Fx<u<Gx
  u_s1_0_C[2,]=-b_s1_0/(2*a_s1_0)
  temp2=costC_s1[1,1]*coef_sq(u_s1_0_C[2,],u_s1_0_C[2,])+(costC_s1[1,2]+costC_s1[2,1])*coef_sq(u_s1_0_C[2,],c(1,0))+
        (costC_s1[1,3]+costC_s1[3,1])*coef_sq(u_s1_0_C[2,],c(0,1))+costC_s1[2:3,2:3]
  #u>Gx
  u_s1_0_C[3,]=t(matG)
  temp3=costC_s1[1,1]*coef_sq(u_s1_0_C[3,],u_s1_0_C[3,])+(costC_s1[1,2]+costC_s1[2,1])*coef_sq(u_s1_0_C[3,],c(1,0))+
        (costC_s1[1,3]+costC_s1[3,1])*coef_sq(u_s1_0_C[3,],c(0,1))+costC_s1[2:3,2:3]
  
  value_s1_0_C=cbind(temp1,temp2,temp3)

costD_s1=matC_s1+region_D_s1[1,1]*coef(temp1_s1,temp1_s1)+(region_D_s1[1,2]+region_D_s1[2,1])*coef(temp1_s1,temp2_s1)+
         region_D_s1[2,2]*coef(temp2_s1,temp2_s1)

  u_s1_0_D=matrix(rep(NA,6),ncol=2)
  a_s1_0=costD_s1[1,1]
  b_s1_0=c(costD_s1[1,2]+costD_s1[2,1],costD_s1[1,3]+costD_s1[3,1])
  #u<Fx
  u_s1_0_D[1,]=t(matF1)
  temp1=costD_s1[1,1]*coef_sq(u_s1_0_D[1,],u_s1_0_D[1,])+(costD_s1[1,2]+costD_s1[2,1])*coef_sq(u_s1_0_D[1,],c(1,0))+
        (costD_s1[1,3]+costD_s1[3,1])*coef_sq(u_s1_0_D[1,],c(0,1))+costD_s1[2:3,2:3]
  #Fx<u<Gx
  u_s1_0_D[2,]=-b_s1_0/(2*a_s1_0)
  temp2=costD_s1[1,1]*coef_sq(u_s1_0_D[2,],u_s1_0_D[2,])+(costD_s1[1,2]+costD_s1[2,1])*coef_sq(u_s1_0_D[2,],c(1,0))+
        (costD_s1[1,3]+costD_s1[3,1])*coef_sq(u_s1_0_D[2,],c(0,1))+costD_s1[2:3,2:3]
  #u>Gx
  u_s1_0_D[3,]=t(matG)
  temp3=costD_s1[1,1]*coef_sq(u_s1_0_D[3,],u_s1_0_D[3,])+(costD_s1[1,2]+costD_s1[2,1])*coef_sq(u_s1_0_D[3,],c(1,0))+
        (costD_s1[1,3]+costD_s1[3,1])*coef_sq(u_s1_0_D[3,],c(0,1))+costD_s1[2:3,2:3]

  value_s1_0_D=cbind(temp1,temp2,temp3)

costE_s1=matC_s1+region_E_s1[1,1]*coef(temp1_s1,temp1_s1)+(region_E_s1[1,2]+region_E_s1[2,1])*coef(temp1_s1,temp2_s1)+
         region_E_s1[2,2]*coef(temp2_s1,temp2_s1)

  u_s1_0_E=matrix(rep(NA,6),ncol=2)
  a_s1_0=costE_s1[1,1]
  b_s1_0=c(costE_s1[1,2]+costE_s1[2,1],costE_s1[1,3]+costE_s1[3,1])
  #u<Fx
  u_s1_0_E[1,]=t(matF1)
  temp1=costE_s1[1,1]*coef_sq(u_s1_0_E[1,],u_s1_0_E[1,])+(costE_s1[1,2]+costE_s1[2,1])*coef_sq(u_s1_0_E[1,],c(1,0))+
        (costE_s1[1,3]+costE_s1[3,1])*coef_sq(u_s1_0_E[1,],c(0,1))+costE_s1[2:3,2:3]
  #Fx<u<Gx
  u_s1_0_E[2,]=-b_s1_0/(2*a_s1_0)
  temp2=costE_s1[1,1]*coef_sq(u_s1_0_E[2,],u_s1_0_E[2,])+(costE_s1[1,2]+costE_s1[2,1])*coef_sq(u_s1_0_E[2,],c(1,0))+
        (costE_s1[1,3]+costE_s1[3,1])*coef_sq(u_s1_0_E[2,],c(0,1))+costE_s1[2:3,2:3]
  #u>Gx
  u_s1_0_E[3,]=t(matG)
  temp3=costE_s1[1,1]*coef_sq(u_s1_0_E[3,],u_s1_0_E[3,])+(costE_s1[1,2]+costE_s1[2,1])*coef_sq(u_s1_0_E[3,],c(1,0))+
        (costE_s1[1,3]+costE_s1[3,1])*coef_sq(u_s1_0_E[3,],c(0,1))+costE_s1[2:3,2:3]

  value_s1_0_E=cbind(temp1,temp2,temp3)



#####

#i=0,s2
#####
region_A_s2=pr[2,1]*(value_s1_1[1:2,5:6])+pr[2,2]*(value_s2_1[1:2,5:6])
region_B_s2=pr[2,1]*(value_s1_1[1:2,3:4])+pr[2,2]*(value_s2_1[1:2,5:6])
region_C_s2=pr[2,1]*(value_s1_1[1:2,3:4])+pr[2,2]*(value_s2_1[1:2,3:4])
region_D_s2=pr[2,1]*(value_s1_1[1:2,1:2])+pr[2,2]*(value_s2_1[1:2,3:4])
region_E_s2=pr[2,1]*(value_s1_1[1:2,1:2])+pr[2,2]*(value_s2_1[1:2,1:2])

#s2個別的cost function(不同的u,在限制下)(i=0)
costA_s2=matC_s2+region_A_s2[1,1]*coef(temp1_s2,temp1_s2)+(region_A_s2[1,2]+region_A_s2[2,1])*coef(temp1_s2,temp2_s2)+
         region_A_s2[2,2]*coef(temp2_s2,temp2_s2)

  u_s2_0_A=matrix(rep(NA,6),ncol=2)
  a_s2_0=costA_s2[1,1]
  b_s2_0=c(costA_s2[1,2]+costA_s2[2,1],costA_s2[1,3]+costA_s2[3,1])
  #u<Fx
  u_s2_0_A[1,]=t(matF1)
  temp1=costA_s2[1,1]*coef_sq(u_s2_0_A[1,],u_s2_0_A[1,])+(costA_s2[1,2]+costA_s2[2,1])*coef_sq(u_s2_0_A[1,],c(1,0))+
        (costA_s2[1,3]+costA_s2[3,1])*coef_sq(u_s2_0_A[1,],c(0,1))+costA_s2[2:3,2:3]
  #Fx<u<Gx
  u_s2_0_A[2,]=-b_s2_0/(2*a_s2_0)
  temp2=costA_s2[1,1]*coef_sq(u_s2_0_A[2,],u_s2_0_A[2,])+(costA_s2[1,2]+costA_s2[2,1])*coef_sq(u_s2_0_A[2,],c(1,0))+
        (costA_s2[1,3]+costA_s2[3,1])*coef_sq(u_s2_0_A[2,],c(0,1))+costA_s2[2:3,2:3]
  #u>Gx
  u_s2_0_A[3,]=t(matG)
  temp3=costA_s2[1,1]*coef_sq(u_s2_0_A[3,],u_s2_0_A[3,])+(costA_s2[1,2]+costA_s2[2,1])*coef_sq(u_s2_0_A[3,],c(1,0))+
        (costA_s2[1,3]+costA_s2[3,1])*coef_sq(u_s2_0_A[3,],c(0,1))+costA_s2[2:3,2:3]

  value_s2_0_A=cbind(temp1,temp2,temp3)


costB_s2=matC_s2+region_B_s2[1,1]*coef(temp1_s2,temp1_s2)+(region_B_s2[1,2]+region_B_s2[2,1])*coef(temp1_s2,temp2_s2)+
         region_B_s2[2,2]*coef(temp2_s2,temp2_s2)

  u_s2_0_B=matrix(rep(NA,6),ncol=2)
  a_s2_0=costB_s2[1,1]
  b_s2_0=c(costB_s2[1,2]+costB_s2[2,1],costB_s2[1,3]+costB_s2[3,1])
  #u<Fx
  u_s2_0_B[1,]=t(matF1)
  temp1=costB_s2[1,1]*coef_sq(u_s2_0_B[1,],u_s2_0_B[1,])+(costB_s2[1,2]+costB_s2[2,1])*coef_sq(u_s2_0_B[1,],c(1,0))+
        (costB_s2[1,3]+costB_s2[3,1])*coef_sq(u_s2_0_B[1,],c(0,1))+costB_s2[2:3,2:3]
  #Fx<u<Gx
  u_s2_0_B[2,]=-b_s2_0/(2*a_s2_0)
  temp2=costB_s2[1,1]*coef_sq(u_s2_0_B[2,],u_s2_0_B[2,])+(costB_s2[1,2]+costB_s2[2,1])*coef_sq(u_s2_0_B[2,],c(1,0))+
        (costB_s2[1,3]+costB_s2[3,1])*coef_sq(u_s2_0_B[2,],c(0,1))+costB_s2[2:3,2:3]
  #u>Gx
  u_s2_0_B[3,]=t(matG)
  temp3=costB_s2[1,1]*coef_sq(u_s2_0_B[3,],u_s2_0_B[3,])+(costB_s2[1,2]+costB_s2[2,1])*coef_sq(u_s2_0_B[3,],c(1,0))+
        (costB_s2[1,3]+costB_s2[3,1])*coef_sq(u_s2_0_B[3,],c(0,1))+costB_s2[2:3,2:3]
  
  value_s2_0_B=cbind(temp1,temp2,temp3)
  
costC_s2=matC_s2+region_C_s2[1,1]*coef(temp1_s2,temp1_s2)+(region_C_s2[1,2]+region_C_s2[2,1])*coef(temp1_s2,temp2_s2)+
         region_C_s2[2,2]*coef(temp2_s2,temp2_s2)
  
  u_s2_0_C=matrix(rep(NA,6),ncol=2)
  a_s2_0=costC_s2[1,1]
  b_s2_0=c(costC_s2[1,2]+costC_s2[2,1],costC_s2[1,3]+costC_s2[3,1])
  #u<Fx
  u_s2_0_C[1,]=t(matF1)
  temp1=costC_s2[1,1]*coef_sq(u_s2_0_C[1,],u_s2_0_C[1,])+(costC_s2[1,2]+costC_s2[2,1])*coef_sq(u_s2_0_C[1,],c(1,0))+
        (costC_s2[1,3]+costC_s2[3,1])*coef_sq(u_s2_0_C[1,],c(0,1))+costC_s2[2:3,2:3]
  #Fx<u<Gx
  u_s2_0_C[2,]=-b_s2_0/(2*a_s2_0)
  temp2=costC_s2[1,1]*coef_sq(u_s2_0_C[2,],u_s2_0_C[2,])+(costC_s2[1,2]+costC_s2[2,1])*coef_sq(u_s2_0_C[2,],c(1,0))+
        (costC_s2[1,3]+costC_s2[3,1])*coef_sq(u_s2_0_C[2,],c(0,1))+costC_s2[2:3,2:3]
  #u>Gx
  u_s2_0_C[3,]=t(matG)
  temp3=costC_s2[1,1]*coef_sq(u_s2_0_C[3,],u_s2_0_C[3,])+(costC_s2[1,2]+costC_s2[2,1])*coef_sq(u_s2_0_C[3,],c(1,0))+
        (costC_s2[1,3]+costC_s2[3,1])*coef_sq(u_s2_0_C[3,],c(0,1))+costC_s2[2:3,2:3]
  
  value_s2_0_C=cbind(temp1,temp2,temp3)
  
costD_s2=matC_s2+region_D_s2[1,1]*coef(temp1_s2,temp1_s2)+(region_D_s2[1,2]+region_D_s2[2,1])*coef(temp1_s2,temp2_s2)+
         region_D_s2[2,2]*coef(temp2_s2,temp2_s2)
  
  u_s2_0_D=matrix(rep(NA,6),ncol=2)
  a_s2_0=costD_s2[1,1]
  b_s2_0=c(costD_s2[1,2]+costD_s2[2,1],costD_s2[1,3]+costD_s2[3,1])
  #u<Fx
  u_s2_0_D[1,]=t(matF1)
  temp1=costD_s2[1,1]*coef_sq(u_s2_0_D[1,],u_s2_0_D[1,])+(costD_s2[1,2]+costD_s2[2,1])*coef_sq(u_s2_0_D[1,],c(1,0))+
        (costD_s2[1,3]+costD_s2[3,1])*coef_sq(u_s2_0_D[1,],c(0,1))+costD_s2[2:3,2:3]
  #Fx<u<Gx
  u_s2_0_D[2,]=-b_s2_0/(2*a_s2_0)
  temp2=costD_s2[1,1]*coef_sq(u_s2_0_D[2,],u_s2_0_D[2,])+(costD_s2[1,2]+costD_s2[2,1])*coef_sq(u_s2_0_D[2,],c(1,0))+
        (costD_s2[1,3]+costD_s2[3,1])*coef_sq(u_s2_0_D[2,],c(0,1))+costD_s2[2:3,2:3]
  #u>Gx
  u_s2_0_D[3,]=t(matG)
  temp3=costD_s2[1,1]*coef_sq(u_s2_0_D[3,],u_s2_0_D[3,])+(costD_s2[1,2]+costD_s2[2,1])*coef_sq(u_s2_0_D[3,],c(1,0))+
        (costD_s2[1,3]+costD_s2[3,1])*coef_sq(u_s2_0_D[3,],c(0,1))+costD_s2[2:3,2:3]
  
  value_s2_0_D=cbind(temp1,temp2,temp3)
  
costE_s2=matC_s2+region_E_s2[1,1]*coef(temp1_s2,temp1_s2)+(region_E_s2[1,2]+region_E_s2[2,1])*coef(temp1_s2,temp2_s2)+
         region_E_s2[2,2]*coef(temp2_s2,temp2_s2)
  
  u_s2_0_E=matrix(rep(NA,6),ncol=2)
  a_s2_0=costE_s2[1,1]
  b_s2_0=c(costE_s2[1,2]+costE_s2[2,1],costE_s2[1,3]+costE_s2[3,1])
  #u<Fx
  u_s2_0_E[1,]=t(matF1)
  temp1=costE_s2[1,1]*coef_sq(u_s2_0_E[1,],u_s2_0_E[1,])+(costE_s2[1,2]+costE_s2[2,1])*coef_sq(u_s2_0_E[1,],c(1,0))+
        (costE_s2[1,3]+costE_s2[3,1])*coef_sq(u_s2_0_E[1,],c(0,1))+costE_s2[2:3,2:3]
  #Fx<u<Gx
  u_s2_0_E[2,]=-b_s2_0/(2*a_s2_0)
  temp2=costE_s2[1,1]*coef_sq(u_s2_0_E[2,],u_s2_0_E[2,])+(costE_s2[1,2]+costE_s2[2,1])*coef_sq(u_s2_0_E[2,],c(1,0))+
        (costE_s2[1,3]+costE_s2[3,1])*coef_sq(u_s2_0_E[2,],c(0,1))+costE_s2[2:3,2:3]
  #u>Gx
  u_s2_0_E[3,]=t(matG)
  temp3=costE_s2[1,1]*coef_sq(u_s2_0_E[3,],u_s2_0_E[3,])+(costE_s2[1,2]+costE_s2[2,1])*coef_sq(u_s2_0_E[3,],c(1,0))+
        (costE_s2[1,3]+costE_s2[3,1])*coef_sq(u_s2_0_E[3,],c(0,1))+costE_s2[2:3,2:3]
  
  value_s2_0_E=cbind(temp1,temp2,temp3)

#####  

#各階段斜率係數
#####
u_s1_1
m_s1_2=matrix(rep(NA,4),ncol=2)
m_s1_2[1,]=c(u_s1_2[1,]-u_s1_2[2,])
m_s1_2[2,]=c(u_s1_2[3,]-u_s1_2[2,])
m_s1_2

m_s2_2=matrix(rep(NA,4),ncol=2)
m_s2_2[1,]=c(u_s2_2[1,]-u_s2_2[2,])
m_s2_2[2,]=c(u_s2_2[3,]-u_s2_2[2,])
m_s2_2

m_s1_1=matrix(rep(NA,4),ncol=2)
m_s1_1[1,]=-(u_s1_1[1,]-u_s1_1[2,])
m_s1_1[2,]=c(u_s1_1[3,]-u_s1_1[2,])
m_s1_1

m_s2_1=matrix(rep(NA,4),ncol=2)
m_s2_1[1,]=-c(u_s2_1[1,]-u_s2_1[2,])
m_s2_1[2,]=c(u_s2_1[3,]-u_s2_1[2,])
m_s2_1

m_s1_A_0=matrix(c(u_s1_0_A[1,]-u_s1_0_A[2,],u_s1_0_A[3,]-u_s1_0_A[2,]),ncol=2,byrow =TRUE)
m_s1_B_0=matrix(c(u_s1_0_B[1,]-u_s1_0_B[2,],u_s1_0_B[3,]-u_s1_0_B[2,]),ncol=2,byrow =TRUE)
m_s1_C_0=matrix(c(u_s1_0_C[1,]-u_s1_0_C[2,],u_s1_0_C[3,]-u_s1_0_C[2,]),ncol=2,byrow =TRUE)
m_s1_D_0=matrix(c(u_s1_0_D[1,]-u_s1_0_D[2,],u_s1_0_D[3,]-u_s1_0_D[2,]),ncol=2,byrow =TRUE)
m_s1_E_0=matrix(c(u_s1_0_E[1,]-u_s1_0_E[2,],u_s1_0_E[3,]-u_s1_0_E[2,]),ncol=2,byrow =TRUE)

m_s2_A_0=matrix(c(u_s2_0_A[1,]-u_s2_0_A[2,],u_s2_0_A[3,]-u_s2_0_A[2,]),ncol=2,byrow =TRUE)
m_s2_B_0=matrix(c(u_s2_0_B[1,]-u_s2_0_B[2,],u_s2_0_B[3,]-u_s2_0_B[2,]),ncol=2,byrow =TRUE)
m_s2_C_0=matrix(c(u_s2_0_C[1,]-u_s2_0_C[2,],u_s2_0_C[3,]-u_s2_0_C[2,]),ncol=2,byrow =TRUE)
m_s2_D_0=matrix(c(u_s2_0_D[1,]-u_s2_0_D[2,],u_s2_0_D[3,]-u_s2_0_D[2,]),ncol=2,byrow =TRUE)
m_s2_E_0=matrix(c(u_s2_0_E[1,]-u_s2_0_E[2,],u_s2_0_E[3,]-u_s2_0_E[2,]),ncol=2,byrow =TRUE)
#####

#aggregrate coefficients & best ans for each region
#####
m_aggregrate=matrix(c(m_s1_1[2,],m_s2_1[2,],-m_s1_1[1,],-m_s2_1[1,]),ncol=2,byrow=TRUE)
m_aggregrate

temp1=(m_aggregrate[1,1]*temp1_s1+m_aggregrate[1,2]*temp2_s1)[2:3]+(m_aggregrate[1,1]*temp1_s1+m_aggregrate[1,2]*temp2_s1)[1]*u_s1_0_A[3,]
temp2=-m_s1_C_0[1,] #因為交集限制
temp3=(m_aggregrate[3,1]*temp1_s1+m_aggregrate[3,2]*temp2_s1)[2:3]+(m_aggregrate[3,1]*temp1_s1+m_aggregrate[3,2]*temp2_s1)[1]*u_s1_0_D[1,]
temp4=(m_aggregrate[4,1]*temp1_s1+m_aggregrate[4,2]*temp2_s1)[2:3]+(m_aggregrate[4,1]*temp1_s1+m_aggregrate[4,2]*temp2_s1)[1]*u_s1_0_E[1,]

m_s1_0=rbind(temp1,temp2,temp3,temp4)
u_s1_0=rbind(u_s1_0_A[3,],u_s1_0_C[2,],u_s1_0_C[1,],u_s1_0_D[1,],u_s1_0_E[1,]) #since region B 與其linear boundary 的無交集


temp1=(m_aggregrate[1,1]*temp1_s2+m_aggregrate[1,2]*temp2_s2)[2:3]+(m_aggregrate[1,1]*temp1_s2+m_aggregrate[1,2]*temp2_s2)[1]*u_s2_0_A[3,]
temp2=-m_s2_C_0[1,]#因為交集限制
temp3=(m_aggregrate[3,1]*temp1_s2+m_aggregrate[3,2]*temp2_s2)[2:3]+(m_aggregrate[3,1]*temp1_s2+m_aggregrate[3,2]*temp2_s2)[1]*u_s2_0_D[1,]
temp4=(m_aggregrate[4,1]*temp1_s2+m_aggregrate[4,2]*temp2_s2)[2:3]+(m_aggregrate[4,1]*temp1_s2+m_aggregrate[4,2]*temp2_s2)[1]*u_s2_0_E[1,]

m_s2_0=rbind(temp1,temp2,temp3,temp4)
u_s2_0=rbind(u_s2_0_A[3,],u_s2_0_C[2,],u_s2_0_C[1,],u_s2_0_D[1,],u_s2_0_E[1,])
#####
rm(temp1,temp2,temp3,temp4)


partALG=function(x,d,pr,init){ #張數、價差、轉移機率(p11,p12,p21,p22)、初始狀態
  best=NULL
  state=init
  when=NULL #紀錄狀態
  spread=d #紀錄價差(for cost)
  
  if(state==s[1]){
    if(x>0 & as.numeric(m_s1_0[1,1])*x+as.numeric(m_s1_0[1,2])*d<0){
      u=as.numeric(u_s1_0[1,1]*x+u_s1_0[1,2]*d)
    }else if(as.numeric(m_s1_0[1,1])*x+as.numeric(m_s1_0[1,2])*d>0 & as.numeric(m_s1_0[2,1])*x+as.numeric(m_s1_0[2,2])*d>0){
      u=as.numeric(u_s1_0[2,1]*x+u_s1_0[2,2]*d)
    }else if(as.numeric(m_s1_0[2,1])*x+as.numeric(m_s1_0[2,2])*d<0 & as.numeric(m_s1_0[3,1])*x+as.numeric(m_s1_0[3,2])*d>0){
      u=as.numeric(u_s1_0[3,1]*x+u_s1_0[3,2]*d)
    }else if(as.numeric(m_s1_0[3,1])*x+as.numeric(m_s1_0[3,2])*d<0 & as.numeric(m_s1_0[4,1])*x+as.numeric(m_s1_0[4,2])*d>0){
      u=as.numeric(u_s1_0[4,1]*x+u_s1_0[4,2]*d)
    }else{u=as.numeric(u_s1_0[5,1]*x-u_s1_0[5,2]*d)}
    
    u=round(u)
    u=max(u,0)
    x=temp1_s1[2]*x+temp1_s1[1]*u
    d=temp2_s1[3]*d+temp2_s1[1]*u
    
  }else{  #s=s[2]
    if(x>0 & as.numeric(m_s2_0[1,1])*x+as.numeric(m_s2_0[1,2])*d<0){
      u=as.numeric(u_s2_0[1,1]*x+u_s2_0[1,2]*d)
    }else if(as.numeric(m_s2_0[1,1])*x+as.numeric(m_s2_0[1,2])*d>0 & as.numeric(m_s2_0[2,1])*x+as.numeric(m_s2_0[2,2])*d>0){
      u=as.numeric(u_s2_0[2,1]*x+u_s2_0[2,2]*d)
    }else if(as.numeric(m_s2_0[2,1])*x+as.numeric(m_s2_0[2,2])*d<0 & as.numeric(m_s2_0[3,1])*x+as.numeric(m_s2_0[3,2])*d>0){
      u=as.numeric(u_s2_0[3,1]*x+u_s2_0[3,2]*d)
    }else if(as.numeric(m_s2_0[3,1])*x+as.numeric(m_s2_0[3,2])*d<0 & as.numeric(m_s2_0[4,1])*x+as.numeric(m_s2_0[4,2])*d>0){
      u=as.numeric(u_s2_0[4,1]*x+u_s2_0[4,2]*d)
    }else{u=as.numeric(u_s2_0[5,1]*x+u_s2_0[5,2]*d)}
    
    u=round(u)
    u=max(u,0)
    x=temp1_s2[2]*x+temp1_s2[1]*u
    d=temp2_s2[3]*d+temp2_s2[1]*u
  }
  
  best=append(best,u)
  when=append(when,state)
  spread=append(spread,d)
  if(state==s[1]){prob=pr[1,]
  }else{prob=pr[2,]}
  state=sample(s,1,replace=TRUE,prob=prob) #狀態隨機
  
  if(state==s[1]){
    if(x>0 & m_s1_1[2,1]*x+m_s1_1[2,2]*d<0){
      u=as.numeric(u_s1_1[3,1]*x+u_s1_1[3,2]*d)
    }else if(m_s1_1[2,1]*x+m_s1_1[2,2]*d>0 & m_s1_1[1,1]*x+m_s1_1[1,2]*d>0){
      u=as.numeric(u_s1_1[2,1]*x+u_s1_1[2,2]*d)
    }else{u=as.numeric(u_s1_1[1,1]*x+u_s1_1[1,2]*d)}
    
    u=round(u)
    u=max(u,0)
    x=temp1_s1[2]*x+temp1_s1[1]*u
    d=temp2_s1[3]*d+temp2_s1[1]*u
    
  }else{ #s=s[2]
    if(x>0 & m_s2_1[2,1]*x+m_s2_1[2,2]*d<0){
      u=as.numeric(u_s2_1[3,1]*x+u_s2_1[3,2]*d)
    }else if(m_s2_1[2,1]*x+m_s2_1[2,2]*d>0 & m_s2_1[1,1]*x+m_s2_1[1,2]*d>0){
      u=as.numeric(u_s2_1[2,1]*x+u_s2_1[2,2]*d)
    }else{u=as.numeric(u_s2_1[1,1]*x+u_s2_1[1,2]*d)}
    
    u=round(u)
    u=max(u,0)
    x=temp1_s2[2]*x+temp1_s2[1]*u
    d=temp2_s2[3]*d+temp2_s2[1]*u
  }
  
  when=append(when,state)
  best=append(best,u)
  spread=append(spread,d)
  if(state==s[1]){prob=pr[1,]
  }else{prob=pr[2,]}
  state=sample(s,1,replace=TRUE,prob=prob) #狀態隨機
  when=append(when,state)
  best=append(best,x) #final 全買
  ls=list(best_order=best,state=when,spread=spread)
  return(ls)
}

init=100
x=100;d=0.1

set.seed(4)
output1=partALG(x,d,pr,init) #(100,100,300)
set.seed(2)
output2=partALG(x,d,pr,init) #(100,300,100)

init=200
x=100;d=0.1

set.seed(4)
output3=partALG(x,d,pr,init) #(300,100,300)
set.seed(2)
output4=partALG(x,d,pr,init) #(300,300,100)

output1
