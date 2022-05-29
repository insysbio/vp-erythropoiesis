### set working directory:
setwd("")
###############################
#### libraries ################
###############################

library("mrgsolve")
library("TruncatedNormal")

input_filename<-'Ery_parameters_ranges.csv'

###############################
###### Model simulation #######
###############################
f_model<-function(pars){
  model_sim <- tryCatch((model %>% 
    data_set(data_switch)%>%
    param(pars) %>%
    mrgsim(
      #delta=1, #10,
      tgrid = tg,
      maxsteps=1e9, 
      hmax = 0,# hmax = 1, 
      # hmin = 0,
      atol = 1e-2, #1e-7
      rtol = 1e-3, #1e-4
      end = 3e4+672
      #end = round(30*24*50.4609)
      #end=1e7+37000
    )),error=function(e) NA)
    
  func<-model_sim$RET_pl_norm[which(model_sim$time %in% Times)]
  return(func)
}


###############################
###### Likelihood #############
###############################
f_likelihood<-function(output,m,sd){
  for(i in 1:length(Times)){
    f_list<-dnorm((output[i]-m[i])/sd[i])
    F_a<-pnorm((0-m[i])/sd[i])
    F_b<-1
    f_list<-f_list/(sd[i]*(F_b-F_a))
    L_list<-log(f_list)
    L=sum(L_list)
  }
  return(L)
}


###############################
###### Poster #################
###############################
P_poster<-function(pars){
  f_list=dnorm((pars-mu_prior)/sigma_prior)
  F_a<-pnorm((0-mu_prior)/sigma_prior)
  F_b<-1
  f_list<-f_list/(sigma_prior*(F_b-F_a))
  P_list<-log(f_list)
  P=sum(P_list)
  return(P)
}


###############################
###### MCMC ###################
###############################
Step_MCMC<-function(pars_cur,m,sd,Likelihood,beta=1){
  Step_made<-FALSE
  p_cur<-Likelihood
  iter<-0
  while(!Step_made){
    pars_prop<-pars_cur+t(mvrandn(l=-pars_cur,u=rep(Inf,Npars),Sig = diag(width^2),n=1))
    colnames(pars_prop)<-Par_names
    if(prod(ifelse(pars_prop>0,1,0))){
      iter<-iter+1
      print(iter)
      output<-f_model(pars_prop[1,])
      L<-sum(f_likelihood(output,m,sd))
      P<-P_poster(pars_prop)
      p_prop<-L+P
      p_accept<-exp((p_prop-p_cur)/beta)
      rnd<-runif(1,0,1)
      if(! is.nan(p_accept)){
        if(rnd<p_accept){
          Step_made<-TRUE
        }
      }
      if((iter>10000)&(!Step_made)){
        break
      }
    }
  }
  result<-list(pars_prop,p_prop,output,Step_made)
  names(result)<-c("pars","L","output","Step_made")
  return(result)
}


################################
##### model compilation ########
################################
tg0 = tgrid(0, 3e4, 3e3)
tg1 = tgrid(3e4, 3e4+672, 1)

tg <- c(tg0, tg1)

time_ev_1 = 30000
time_ev_2 = 30000.01

data_switch<-as_data_set(ev(ID = 1, cmt = 'switch_inj', time = time_ev_1, amt = 1, ii = 0, evid = 1, addl = 0),
                         ev(ID = 1, cmt = 'switch_inj', time = time_ev_2, amt = -1, ii = 0, evid = 1, addl = 0))
data_switch$ID[2]<-1

model_full<-mread(model="200513_mrgsolve_output_dose_ev",file = "200513_mrgsolve_output_dose_ev.cpp")

model<<-model_full

parameters_table<-read.table(input_filename, sep=';', header=T, stringsAsFactors = F)
Par_names<<-parameters_table[,'parameter_name']
Par_list<-param(model)

mu<<-as.numeric(Par_list[Par_names])
names(mu)<-Par_names

Npars<<-length(Par_names)

Times<<-c(3e4+c(24,	48,	72,	96,	120,	144,	168,	192,	216,	240,	264,	288,	312,	336,	504,	672))


###############################
##### set prior ###############
###############################
mu_prior<<-mu
sigma_prior<<-mu/2
width<<-mu/2

################################
##### experimental data ########
################################
m <- c(1.05818181818182,	1.12363636363636,	0.92969696969697,	1.55151515151515,	1.74545454545455,	1.74545454545455,	2.3030303030303,	1.55151515151515,	3.05454545454545,	1.80606060606061,	1.49090909090909,	1.11272727272727,	0.995151515151515,	1.43030303030303,	0.561212121212121,	0.678787878787879)
sd <- c(0.696569054809025,	0.441792824645413,	0.393005886954508,	0.514973231181769,	0.704700211090843,	0.596284793999944,	0.786011773909017,	0.92153104527264,	1.92437365336346,	1.84306209054528,	1.46360813072713,	0.872744107581736,	0.593574408572671,	0.650492502545393,	0.417399355799961,	0.317115094990879)


###############################
####### Main ##################
###############################
N_pars_sample<-207

Likelihood<-rep(0,N_pars_sample)

Matr_pars<-matrix(0,N_pars_sample,Npars)

colnames(Matr_pars)<-Par_names

Matr_pars[1,]<-mu_prior+t(mvrandn(l=-mu_prior,u=rep(Inf,Npars),Sig=diag(sigma_prior^2),n=1))


output_0<-f_model(Matr_pars[1,])



L<-f_likelihood(output_0,m,sd)
P<-P_poster(Matr_pars[1,])
Likelihood[1]<-sum(L)+P
Likelihood[1]
Output<-matrix(0,nrow = N_pars_sample,ncol = length(Times))
Output[1,]<-c(output_0)

Output

for(i in 2:N_pars_sample){
  print(paste("iteration =",i,sep = " "))
  result_MCMC<-Step_MCMC(Matr_pars[i-1,],m,sd,Likelihood[i-1],beta = 1)
  if(!result_MCMC[["Step_made"]]){
    break
  }
  Matr_pars[i,]<-result_MCMC[["pars"]]
  Likelihood[i]<-result_MCMC[["L"]]
  Output[i,]<-result_MCMC[["output"]]
  print(Likelihood[i])
  write.csv2(Output,"MCMC_output_207.csv",row.names = F)
}

m
apply(Output,2,mean)
sd
apply(Output,2,sd)
    
write.csv2(Matr_pars,"MCMC_matr_pars_207.csv",row.names = F)
