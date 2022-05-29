setwd("")
library(mrgsolve)

IniDataMatrix<-read.table("ERY_param_sel_4000_sample.csv", header = TRUE, sep=';', dec = ",")

#pars<-IniDataMatrix[i,1:2]
#Npar<-24

model_full<-mread(model="200513_mrgsolve_output_dose_ev",file = "200513_mrgsolve_output_dose_ev.cpp")
model<<-model_full

tg0 = tgrid(0, 3e4, 3e3)
tg1 = tgrid(3e4, 3e4+672, 1)

tg <- c(tg0, tg1)

time_ev_1 = 30000
time_ev_2 = 30000.01

data_switch<-as_data_set(ev(ID = 1, cmt = 'switch_inj', time = time_ev_1, amt = 1, ii = 0, evid = 1, addl = 0),
                         ev(ID = 1, cmt = 'switch_inj', time = time_ev_2, amt = -1, ii = 0, evid = 1, addl = 0))
data_switch$ID[2]<-1

Times<-3e4+c(24,	48,	72,	96,	120,	144,	168,	192,	216,	240,	264,	288,	312,	336,	504,	672)
#Times <- 100*round(Times/100)
#Times <- c(10,100)

res_1<-matrix(NA,nrow = nrow(IniDataMatrix),ncol = length(Times))
save_numbers = c()
#res_2<-matrix(0,nrow = nrow(IniDataMatrix),ncol = length(Times)+1)

for(i in 1:nrow(IniDataMatrix)){
  #for(i in 1:1){  
  data<-IniDataMatrix[i,] 
  
  model_sim <- tryCatch((model %>% 
    data_set(data_switch)%>%
    param(data) %>%
    # idata_set(idata) %>%
    #param(switch_dupi = 1) %>%
    #param(switch_dupi = 1 , K_portion_nls_adl = 0.25) %>%
    #amt = 46666.667; 233333.33; 933333.33
    #ev(ID = 1, cmt = 'AAb_sc_amt_', time = 10000000, amt = 46666.667, ii = 672, evid = 1, addl = 20) %>%
    #ev(ID = 1, cmt = 'EPO_sc', time = 3e4, amt = 6084750, ii = 0, evid = 1, addl = 0, end = 3e4+0.01) %>%
    mrgsim(
      #delta=1, #10,
      tgrid = tg,
      maxsteps=1e9, 
      hmax = 0,# hmax = 1, 
      # hmin = 0,
      atol = 1e-9, #1e-7
      rtol = 1e-4, #1e-4
      end = 3e4+672
      #end = round(30*24*50.4609)
      #end=1e7+37000
    )),error=function(e) NA)
  
  if(! is.na(model_sim)){
    res_1[i,]<-model_sim$RET_pl_norm[which(model_sim$time %in% Times)]
    save_numbers = c(save_numbers,i)
  }
  
  #res_2[i,]<-model_sim$Pruritus_score[which(model_sim$time %in% Times)]
  
}

res_1_full=res_1


res_1 = res_1[save_numbers,]
res_1_nan <- na.omit(res_1)

#res_1_1<-res_1[,1:ncol(res_1)]
#res_1<-res_1[,2:ncol(res_1)]

write.csv2(res_1,"RET_pl_norm_sample_5.csv",row.names = F)
write.csv2(res_1_nan,"RET_pl_norm_sample_5_wt_nan.csv",row.names = F)
write.csv2(save_numbers,"Save_numbers_sample_5.csv",row.names = F)
write.csv2(res_1_full,"RET_pl_norm_sample_5_full.csv",row.names = F)

###Selection for parameters###
Full_frame_par_output<-cbind(IniDataMatrix,res_1_full)
Save_numbers_frame_par_output<-Full_frame_par_output[save_numbers,]
Save_numbers_frame_par_output_nan<-na.omit(Save_numbers_frame_par_output)
Save_numbers_frame_par_output_final<-Save_numbers_frame_par_output_nan[,2:56]
Save_numbers_frame_par_output_final<-Save_numbers_frame_par_output_final[,1:39]
write.csv2(Save_numbers_frame_par_output_final,"Parameters_sample_success_simulation_SBS.csv",row.names = F)

###Testing plot###
model_sim$ACC_nle
plot(model_sim$time[11:236],model_sim$EASI_rel[11:236])
plot(model_sim, ~EASI)
