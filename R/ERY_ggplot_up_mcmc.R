setwd("/Users/galinakolesova/Documents/Work_in_home/ERY/Work/R_methodics/Sample_207/MCMC")
library('ggplot2')
library("RColorBrewer")
library('tidyverse')
library('magrittr')
library('readr')
library('dplyr')
file1<-'Ery_ret_exp_day.txt'
file2<-'Ery_ret_pred_day_mcmc.txt'
ylabels_file<-'ylabels_Ery.csv'
target_functions<-c('RET_pl_norm')
argument_name<-'Time_day'
argument_max_value<-28
title_x_axis<-"Time, days" 
#plot_title<-'Reticulocyte numbers versus time profiles \n S_patients=5, Epoetin alpha dose = 300 IU/kg, single dose'
plot_title<-'Reticulocyte numbers versus time profiles'
plot_subtitle<-'S_patients=207, Epoetin alpha dose = 300 IU/kg, single dose' 
ylabels_type<-'from_table'
set_SD_model_title1<- 'Experimental data: SD
Reticulocytes, PMID:15317827'
set_SD_model_title2<- 'Model simulation: SD
Reticulocytes'

set_mean_model_title1<- 'Experimental data: mean                         
Reticulocytes, PMID:15317827'
set_mean_model_title2<- 'Model simulation: mean                         
Reticulocytes'

#caption_title<- 'References: \n 1.	Alexander Stepanov, Galina Lebedeva, ACoP 10 (2019) From stem cell to erythrocyte and platelet: QSP model of erythropoiesis and thrombopoiesis for assessing the impact of pharmacological interventions. Poster. \n 2.	Galina Kolesova, Oleg Demin, PAGE 28 (2019) Abstr 9008 [www.page-meeting.org/?abstract=9008]. \n 3. Experimental data sourse: PMID:15317827'#title describing conditions of simulations for file1 (drug doses,for example)
#dataset_label_exp<-'ibr420QD_3'
medicine_in_set1 <-'zan160BID' # 'ibr560QD', 'aca100BID', 'zan160BID', 'zan320QD' for target PK functions automatic choice, 'other1', 'other2'  - to plot other medicines with manually target setting
medicine_in_set2 <-'aca100BID'
medicine_in_set3 <- 'zan320QD'
width<-6 #width of the plot when saving pdf
height<-4 #height of the plot when saving pdf


Data_parset_1 <- read_tsv(file1)
dataset_label<-rep(medicine_in_set1, times=length(Data_parset_1[,1]))
Data_parset_1<-cbind(Data_parset_1, dataset_label)

Data_parset_2 <- read_tsv(file2)
dataset_label<-rep(medicine_in_set2, times=length(Data_parset_2[,1]))
Data_parset_2<-cbind(Data_parset_2, dataset_label)

Data_parset<-rbind(Data_parset_1[,c(target_functions, argument_name,'y_sd_min','y_sd_max','dataset_label')], Data_parset_2[,c(target_functions,argument_name,'y_sd_min','y_sd_max','dataset_label')])

if (ylabels_type=='from_table'){
  ylabels<-read.table(ylabels_file, sep=';', header = T, stringsAsFactors = F)
}

#plotting
fill_colours<-c("#E41A1C","#525252","#525252","#525252","#525252","#377EB8","#377EB8","#377EB8","#4DAF4A","#F781BF")
names(fill_colours)<-c('zan160BID','zan320QD','zan40QD','zan80QD','zan160QD','ibr560QD','ibr420QD','ibr350QD','aca100BID','aca200QD')

for (model_function in target_functions){
  if (ylabels_type=='from_model'){
    ylabel=model_function
  }
  if (ylabels_type=='default'){
    ylabel=ylabel_default
  }
  if (ylabels_type=='from_table'){
    ylabel=ylabels[which(ylabels$model_name==model_function),'yaxis_name']
    ymin=ylabels[which(ylabels$model_name==model_function),'ymin']
    ymax=ylabels[which(ylabels$model_name==model_function),'ymax']
  }
  pdf_name<-paste0(model_function,'.png')
  
  df <- data.frame(
    x = Data_parset[,"Time_day"],
    y = Data_parset[,model_function],
    curve = Data_parset$dataset_label,
    ymin = Data_parset$y_sd_min,
    ymax = Data_parset$y_sd_max
  )
  
  ggplot()+
    geom_ribbon(data=df, mapping = aes(x=x, ymin=ymin, ymax=ymax, fill=curve), alpha=0.4)+
    geom_line(data=df, mapping = aes(x=x,y=y, group = curve, color = curve))+
    #geom_errorbar(data = experimental_data_for_function, mapping = aes(x=time, ymin=lower, ymax=upper, width=0.3, colour=medicine), size=0.3)+
    #geom_point(data = experimental_data_for_function, mapping = aes(x=time, y=median, shape=medicine), size=2)+
    #ggtitle(plot_title)+
    labs(title = plot_title, subtitle = plot_subtitle)+
    scale_alpha_manual(guide=F)+
    scale_fill_manual(values=fill_colours,
                      breaks=c(medicine_in_set1,medicine_in_set2),
                      labels=c(set_SD_model_title1,set_SD_model_title2))+
    scale_colour_manual(values=fill_colours,
                        breaks=c(medicine_in_set1,medicine_in_set2),
                        labels=c(set_mean_model_title1,set_mean_model_title2))+
    #scale_shape_manual(values=c(16),
    #breaks=c(medicine_in_set),
    #labels=c(set_experimental_title))+
    scale_x_continuous(breaks=c(0, 1,	4,	7,	10,	13,	16,	19,	22,	25,	28))+
    scale_y_continuous(limits=c(ymin, ymax))+
    #scale_y_continuous(trans = 'log10')+
    #      ylim(ymin,max(experimental_data_for_function$upper))+
    theme(
      panel.grid.major = element_line(size = 0.1, linetype = 'solid',colour = "gray"),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill='white', colour='black'),
      legend.position="bottom",
      legend.title=element_blank(),
      legend.text = element_text(color = "black", size = 4.4),
      legend.box.margin = margin(0.2, 0.2, 0.2, 0.2, "cm"),
      legend.text.align = 0,
      legend.title.align = 0,
      legend.box.just = 'left'
    )+
    xlab(title_x_axis)+
    ylab(ylabel)
  ggsave(pdf_name, width = width, height = height)
  
}


