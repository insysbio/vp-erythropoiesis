setwd("")
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library('nleqslv')
library('VGAM')
library('EnvStats')


calc_delta<-function(x,a,b,m,sd){
  
  mu = x["mu"]
  sigma = x["sigma"]
  
  a1 = (1/2)*(sqrt(2)*(-sigma^2+log(a)-mu)/sigma)
  a2 = (1/2)*(sqrt(2)*(sigma^2+mu-log(b))/sigma)
  a3 = (1/2)*(sqrt(2)*(log(a)-mu)/sigma) 
  a4 = (1/2)*(sqrt(2)*(-log(b)+mu)/sigma)
  a5 = (1/2)*(sqrt(2)*(-2*sigma^2+log(a)-mu)/sigma)
  a6 = (1/2)*(sqrt(2)*(2*sigma^2+mu-log(b))/sigma)  
  a7 = (1/2)*(sqrt(2)*(-2*sigma^2+log(b)-mu)/sigma)
  a8 = (1/2)*(sqrt(2)*(log(b)-mu)/sigma)  
  
  N = exp(sigma^2)*(erf(a3)*erf(a5)*exp(sigma^2)+
                      erf(a3)*erf(a6)*exp(sigma^2)+
                      erf(a4)*erf(a5)*exp(sigma^2)+
                      erf(a4)*erf(a6)*exp(sigma^2)-
                      (erf(a1))^2-
                      2*erf(a1)*erf(a2)-
                      (erf(a2))^2)*exp(2*mu)   
  
  Et = exp(1/2*(sigma^2))*exp(mu)*(erf(a1)+erf(a2))/(erf(a3)+erf(a4))
  
  Vt = N/(erf(a3)+erf(a4))^2  
  
  eq1=(Et - m)
  eq2=Vt-sd^2
  
  eqs=c(eq1,eq2)
  names(eqs)=c("delta_m","delta_sd")
  
  return(eqs)
}


#Par_names<<-c("syn_factor_th1","syn_factor_th2","syn_factor_th17","syn_factor_th22","syn_factor_treg","syn_factor_pbm","syn_factor_idec","syn_factor_m1","syn_factor_mc","syn_factor_kc","kmax_0_il4_pro_kcss_le","Imax_0_il4_syn_pfss_le","Imax_0_il4_syn_pfsg_le","Imax_0_il4_syn_plbss_le","Imax_0_il4_tr_plbsg_le","kmax_tnfa_il22_pro_kcss_le","Imax_0_il22_syn_pfss_le","Imax_0_il22_syn_pfsg_le","Imax_0_il22_tr_pfsg_le","K_portion_nls","k0base_influx_agene","kbase_dif_pbm_lc_pl_ld","Imax_0_il4_dif_pbm_lc_pl", "Placebo_max_Pruritus",	"Placebo_max_SCORAD",	"Placebo_max_EASI",	"Placebo_max_Sleep_loss")
#Par_list<-param(model)

#mu<<-as.numeric(Par_list[Par_names])
#names(mu)<-Par_names

Npars<<-39

m_list=c(1.05818181818182,	1.12363636363636,	0.92969696969697,	1.55151515151515,	1.74545454545455,	1.74545454545455,	2.3030303030303,	1.55151515151515,	3.05454545454545,	1.80606060606061,	1.49090909090909,	1.11272727272727,	0.995151515151515,	1.43030303030303,	0.561212121212121,	0.678787878787879)
sd_list=c(0.696569054809025,	0.441792824645413,	0.393005886954508,	0.514973231181769,	0.704700211090843,	0.596284793999944,	0.786011773909017,	0.92153104527264,	1.92437365336346,	1.84306209054528,	1.46360813072713,	0.872744107581736,	0.593574408572671,	0.650492502545393,	0.417399355799961,	0.317115094990879)


Times<<-c(30024,	30048,	30072,	30096,	30120,	30144,	30168,	30192,	30216,	30240,	30264,	30288,	30312,	30336,	30504,	30672)
ColNames<-c("X30024", "X30048", "X30072", "X30096", "X30120", "X30144", "X30168", "X30192",	"X30216",	"X30240",	"X30264",	"X30288",	"X30312",	"X30336",	"X30504",	"X30672")

t_fit = 6



rand_sample=read.table("RET_pl_norm_sample_5_wt_nan.csv",header = T,sep=";",dec = ",")

Npatients<-207

y_dist_matrix<-matrix(0,ncol = length(Times), nrow = Npatients)

mu<-1
sigma<-1


for(i in 1:length(Times)){
  t<-Times[i]
  print(t)
  a=min(rand_sample[,ColNames[i]])
  b=max(rand_sample[,ColNames[i]])
  
  m=m_list[i]
  sd=sd_list[i]
  
  m1=m-a
  
  
  x0=c(mu,sigma)
  names(x0)=c("mu","sigma")
  
  result=nleqslv(x0,calc_delta,a=0,b=b-a,m=m1,sd=sd)
  
  #print(result)
  mu<-result$x["mu"]
  sigma<-result$x["sigma"]
  
  y_dist<-rlnormTrunc(Npatients,
                      meanlog = result$x["mu"], 
                      sdlog = result$x["sigma"], 
                      min = 0, max = b-a) + a
  
  print(paste("mean of y_dist = ", mean(y_dist),", m =", m, sep = ""))
  print(paste("sd of y_dist = ",sd(y_dist),", sd =", sd, sep = ""))
  
  
  y_dist<-sort(y_dist)
  
  y_dist_matrix[,i]<-y_dist
  
}

rs<-rand_sample[,paste("X",Times[t_fit],sep="")]
names(rs)<-seq(length(rs))
rs_sorted<-sort(rs)



y_sample<-c()
y_sample_numbers<-c()

for(y in y_dist_matrix[,t_fit]){
  for(i in 1:length(rs_sorted)){
    y_point<-rs_sorted[i]
    if(y_point<y){
      less_point<-y_point
    }else{
      more_point<-y_point
      if(abs(less_point-y)<abs(more_point-y)){
        y_sample<-c(y_sample,less_point)
        y_sample_numbers<-c(y_sample_numbers,names(less_point))
      }else{
        y_sample<-c(y_sample,more_point)
        y_sample_numbers<-c(y_sample_numbers,names(more_point))
      }
      break
      
    }
  }
}

y_sample_matrix=matrix(0, nrow = length(y_sample_numbers), ncol = length(Times))

for(i in 1:length(y_sample_numbers)){
  y_sample_matrix[i,]<-as.numeric(rand_sample[y_sample_numbers[i],])
}

m_sample<-apply(y_sample_matrix,2,mean)
m_matrix<-cbind(m_sample,m_list)

print(m_matrix)
plot(Times,m_matrix[,"m_list"],main="N_patients = 5",ylab="RET_pl_norm",xlab="Time, hours")
lines(Times,m_matrix[,"m_sample"],col="red")


sd_sample<-apply(y_sample_matrix,2,sd)
sd_matrix<-cbind(sd_sample,sd_list)
print(sd_matrix)

print(paste("mean of y_sample = ", mean(y_sample)))
print(paste("sd of y_sample = ", sd(y_sample)))

print("selected rows:")
print(as.numeric(y_sample_numbers))


######t-test#####
Data_sample_selected<-data.frame(y_sample_matrix)
write.csv2(Data_sample_selected, file = "ERY_outputs_four_method_207_sample.csv")

t.test(Data_sample_selected$X1, mu = 1.05818181818182)

t.test(Data_sample_selected$X2, mu = 1.12363636363636)

t.test(Data_sample_selected$X3, mu = 0.92969696969697)

t.test(Data_sample_selected$X4, mu = 1.55151515151515)

t.test(Data_sample_selected$X5, mu = 1.74545454545455)

t.test(Data_sample_selected$X6, mu = 1.74545454545455)

t.test(Data_sample_selected$X7, mu = 2.3030303030303)

t.test(Data_sample_selected$X8, mu = 1.55151515151515)

t.test(Data_sample_selected$X9, mu = 3.05454545454545)

t.test(Data_sample_selected$X10, mu = 1.80606060606061)

t.test(Data_sample_selected$X11, mu = 1.49090909090909)

t.test(Data_sample_selected$X12, mu = 1.11272727272727)

t.test(Data_sample_selected$X13, mu = 0.995151515151515)

t.test(Data_sample_selected$X14, mu = 1.43030303030303)

t.test(Data_sample_selected$X15, mu = 0.561212121212121)

t.test(Data_sample_selected$X16, mu = 0.678787878787879)

