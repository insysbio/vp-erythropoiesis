### set working directory:
setwd("")
###############################
#### libraries ################
###############################
library("mrgsolve")
library("randtoolbox")
library("rngWELL")
library("TruncatedNormal")
library("plotly")

input_filename<-'Ery_parameters_ranges.csv'
###############################
###### Model simulation #######
###############################
f_model<-function(pars){
  model_sim <- model %>% 
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
    )
  
  func<-model_sim$RET_pl_norm[which(model_sim$time %in% Times)]
  return(func)
}


################################
######## F_fit_params ##########
################################
f_fit_params_exact<-function(param,y_series){
  y_model<-f_model(param)
  return(sum((y_model-y_series)^2))
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

lower<<-parameters_table[,'min']
#names(l)<-Par_names
upper<<-parameters_table[,'max']
#names(u)<-Par_names



Npars<<-length(Par_names)

Times<<-c(3e4+c(24,	48,	72,	96,	120,	144,	168,	192,	216,	240,	264,	288,	312,	336,	504,	672))


################################
######## Fitting ###############
################################

m <- c(1.05818181818182,	1.12363636363636,	0.92969696969697,	1.55151515151515,	1.74545454545455,	1.74545454545455,	2.3030303030303,	1.55151515151515,	3.05454545454545,	1.80606060606061,	1.49090909090909,	1.11272727272727,	0.995151515151515,	1.43030303030303,	0.561212121212121,	0.678787878787879)
sd <- c(0.696569054809025,	0.441792824645413,	0.393005886954508,	0.514973231181769,	0.704700211090843,	0.596284793999944,	0.786011773909017,	0.92153104527264,	1.92437365336346,	1.84306209054528,	1.46360813072713,	0.872744107581736,	0.593574408572671,	0.650492502545393,	0.417399355799961,	0.317115094990879)


sigma_ln<-sqrt(log(sd^2/m^2+1))
mu_ln<-log(m)-sigma_ln^2/2


N_f_sample<-5 #207

f_matrix<-matrix(0,nrow = N_f_sample, ncol = length(Times))

for(i in 1:length(Times)){
  f_matrix[,i]<-sort(rlnorm(N_f_sample,mu_ln[i],sigma_ln[i]))
}

f_matrix

par_matrix<-matrix(0,ncol = Npars,nrow = N_f_sample)

for (i in 1:N_f_sample){
  print(paste("iteration =",i,sep = " "))
  
  fitting_result<-tryCatch(optim(par=mu,fn=f_fit_params_exact,
                        y_series=f_matrix[i,],
                        method = "L-BFGS-B",
                        lower = lower,
                        upper = upper), error = function(e) list(par=rep(0,Npars)))
  #fitting_result<-optim(par=mu,fn=f_fit_params_exact,
  #                     y_series=f_matrix[i,])
  #print("par:")
  #print(fitting_result$par)
  #print("conv:")
  #print(fitting_result$convergence)
  par_matrix[i,]<-fitting_result$par
  write.csv2(par_matrix,"Ery_par_fitting_6.csv",row.names = F)
}

#input_filename_par<-'Ery_par_fitting_point.csv'
#par_matrix<-read.table(input_filename_par, sep=';', header=T, stringsAsFactors = F)
#names(par_matrix)<-Par_names

Output_exact<-matrix(0,ncol = length(Times),nrow(par_matrix))

for(i in 1:nrow(par_matrix)){
  print(paste("Iteration",i))
  pars<-par_matrix[i,]
  names(pars)<-Par_names
  print(pars)
  if(sum(pars)>0){
    y_new<-f_model(pars)
    print(y_new)
    Output_exact[i,]<-y_new
  }
}
write.csv2(Output_exact,"Ery_output_fitting.csv",row.names = F)


m
apply(Output_exact[Output_exact[,1]!=0,],2,mean)
plot(m)
lines(apply(Output_exact[Output_exact[,1]!=0,],2,mean), col="red")

sd
apply(Output_exact[Output_exact[,1]!=0,],2,sd)
