m<-1.74545454545455
sd<-0.596284793999944
N<-207

y3<-m
y2<-m-sd/2
y4<-m+sd/2
y1<-m-1.5*sd
y5<-m+1.5*sd

y<-c(0.845885315,1.477492265,1.749343516,2.049552035,2.696683833)
#y<-c(y1,y2,y3,y4,y5)

A<-matrix(0,5,5)

A[1,]<-y
A[2,]<-N/(N-1)*(y-m)^2
A[3,]<-1
A[4,]<-c(y2-y1, y1-y3, y3-y2,0,0)
A[5,]<-c(0,0,y4-y3,y3-y5,y5-y4)

b<-c(m,sd^2,1,0,0)

x<-solve(A,b)

x
y

results<-read.csv2("hist_results.csv")
colnames(results)<-c(Par_names,"y")

y_traj<-matrix(0,ncol = length(f_model_traj(mu)),nrow = nrow(results))

for(i in 1:nrow(results)){
  print(i)
  pars_tmp<-results[i,1:length(Par_names)]
  names(pars_tmp)<-Par_names
  y_traj[i,]<-f_model_traj(pars_tmp)
}
pars_tmp

write.csv2(y_traj,"hist_result_traj.csv",row.names = F)

