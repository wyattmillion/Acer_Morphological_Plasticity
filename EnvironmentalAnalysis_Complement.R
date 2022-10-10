#Separate script for running environmental analysis
#will require some things from MorphologicalPlasticityAnalysis.R

#This script will 
  #1) make maps for manuscript
  #2) produce thermal and water quality metrics for sites
  #3) run Bayesian analysis for effect of environment on trait values

##################
library(tidyverse)
library(lubridate)
library(rworldmap)
library(devtools)
require(gstudio) 
require(ggmap)
library(ggfortify)
library(naniar)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(maps)
library(ggrepel)
library(ggsn)
library(mapdata)
library(acs)
library(choroplethr)
library(choroplethrMaps)
library(RColorBrewer)
library(ggplot2)

setwd("/Volumes/GoogleDrive/My Drive/SCHOOL/LAB/CRRPTransplant") #change path (stuff within "") to the appropriate location
coords<-read.csv("outplant site coordinates.csv",header=T) #change path (stuff within "") to the appropriate location of datasheet if not in current wd
coords<- coords %>% arrange(outplant.site = factor(outplant.site, levels=c("E. Sambo","Marker 32","W. Sambo","Big Pine Shoals","Dave's Ledge","Looe Key","Eastern Dry Rocks","Maryland Shoals","Bahia Honda")))

######1. Making maps for manuscript
#########
#read in the outplant site coordinates to make map
sites<-coords[,1:3]
stations<-read.csv("serc.csv") #change path (stuff within "") to the appropriate location

register_google(key = "AIzaSyDaUEsZAbllpuQvFajCUpTeh1fBxuTbrMg")
api.key.install("208638c72972ed33a5808ab2e6dab52c5878ce37")

#qmplot(Lon, Lat, data=sites)
geocode("Saddlebunch Keys")
bw_map<-get_googlemap(center=c(-81.55,24.6),zoom=10, 
                      color = "bw",
                      style = "feature:road|visibility:off&style=feature:administrative|visibility:off")

#order based on risk scores from MorphologicalPlasticityAnalysis.R
#Genet survival rank order: "36","1","44","50","3","7","31","13","41","62"
colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(10)
colorG<-c("#246BAE","#E1EDF3","#FAE7DC","#DC6F58","#F7B799","#053061","#B6202E","#549EC9","#A6CFE3","#67001F")

#Site survival rank order: "E. Sambo","Dave's Ledge","Marker 32","Big Pine","Looe Key","W. Sambo","EDR","Maryland Shoals","Bahia Honda"
colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(9)
colorS <-c("#67001F","#C1DDEB","#F7F7F7","#2971B1","#BB2A33","#053061","#FACDB5","#6AABD0","#E58267","#666666")

#Site map
ggmap(bw_map,extent="device",legend="bottomleft") +
  geom_point(data=sites,aes(x=Lon,y=Lat),size=4,color=colorS,alpha=0.7)+
  geom_point(data=stations,aes(x=Lon,y=Lat),color="white",pch=18,size=2) +
  geom_text(data=sites,aes(x=Lon+.035, y=Lat-.015, label=outplant.site),size=2)+
  xlab("Longitude") + ylab("Latitude") + 
  scalebar(x.min=-81.9,x.max=-81.2,y.min=24.45,y.max=24.9,dist=8,dist_unit="km",transform=TRUE,model="WGS84",st.bottom=TRUE,st.dist=0.04,st.size=3)
#######

#####2a. Temperature data
#######
#read in the individual site data or just skip this to read in full file with all data together
BH<-read.csv("EnvironmentalData/BHTempTotal.csv", header = T)
BH$Date<-mdy_hm(BH$Date)
BH$Site<-"Bahia Honda"

BP<-read.csv("EnvironmentalData/BPTempTotal.csv", header = T)
BP$Date<-mdy_hm(BP$Date)
BP$Site<-"Big Pine"

LK<-read.csv("EnvironmentalData/LKTempTotal.csv", header = T)
LK$Date<-mdy_hm(LK$Date) #Looe Key missing data from April 2018-October 2018
LK$Site<-"Looe Key"
names(LK)[1]<-"number"

DL<-read.csv("EnvironmentalData/DLTempTotal.csv", header = T)
DL$Date<-mdy_hm(DL$Date)
DL$Site<-"Dave's Ledge"

MS<-read.csv("EnvironmentalData/MSTempTotal.csv", header = T)
MS$Date<-mdy_hm(MS$Date)
MS$Site<-"Maryland Shoals"

ES<-read.csv("EnvironmentalData/ESTempTotal.csv", header = T)
ES$Date<-mdy_hm(ES$Date)
ES$Site<-"E. Sambo"

WS<-read.csv("EnvironmentalData/WSTempTotal.csv", header = T)
WS$Date<-mdy_hm(WS$Date)
WS$Site<-"W. Sambo"

M32<-read.csv("EnvironmentalData/M32TempTotal.csv", header = T)
M32$Date<-mdy_hm(M32$Date)
M32$Site<-"Marker 32"

EDR<-read.csv("EnvironmentalData/EDRTempTotal.csv", header = T)
EDR$Date<-mdy_hm(EDR$Date)
EDR$Site<-"EDR"

temp<-rbind(BH,BP,LK,DL,MS,ES,WS,M32,EDR)

temp<-temp[complete.cases(temp[,3]),]

unique(temp$Site) #should be 9 sites
temp<-temp %>% mutate(Site = factor(Site, levels=c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","Looe Key","EDR","Maryland Shoals","Bahia Honda")))
temp <- temp %>% filter(Date>=as.Date("2018-05-1") & Date <= as.Date("2019-04-30")) #isolating data to only experimental series
#write.csv(temp,"EnvironmentalData/OutplantTempData.csv",row.names = FALSE)

####BEGIN HERE IF USING GITHUB DATA###
temp<-read.csv("EnvironmentalData/OutplantTempData.csv",header = T) #change path (stuff within "") to the appropriate location
temp<-temp %>% mutate(Site = factor(Site, levels=c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","Looe Key","EDR","Maryland Shoals","Bahia Honda")))

summary(temp$Site)
#**Important to note that Looe Key have much fewer observations. That is because the logger flooded from Apr to Oct 2018

#keeping if a temp is greater than certain threshold
temp$DH30.5<-ifelse(temp$Temp>=30.5,1,0) #noting the hours above 30.5
temp$DH32<-ifelse(temp$Temp>=32,1,0) #noting the hours above 32

####Visualizing the temperature day as a time series
quartz()
temp %>% #group_by(Site,day=floor_date(Date,"day")) %>% 
  #summarise(DailyMean=mean(Temp),
  #          SD=sd(Temp)) %>%
  ggplot(aes(x=Date,y=Temp,color=Site))+
  #geom_point(aes(group=Site),position= position_jitter(w=100))+
  geom_line(aes(group=Site,color=Site),position=position_jitter(w=.2))+
  #geom_errorbar(aes(ymin=MonthlyMean-SD,ymax=MonthlyMean+SD),width=.3)+
  scale_color_manual(values = colorS)+ylab("Hourly Temperature")+
  #geom_vline(xintercept = temp[1465,2])+ #to help visualize start of summer
  #geom_vline(xintercept = temp[3672,2])+ #to help visualize end of summer
  theme_classic()+facet_wrap(~Site)+theme(legend.position = 'none')

#Figure 1 will include the time series data with an insert showing hourly summer temps
p<-temp %>% group_by(Site,day=floor_date(Date,"day")) %>% 
  summarise(DailyAvg=mean(Temp),
            SD=sd(Temp)) %>%
  ggplot(aes(x=day,y=DailyAvg))+
  geom_line(aes(group=Site,color=Site),position=position_jitter(w=.2))+ #sites colored from best to worst survival
  #stat_summary(fun=mean,geom = 'line',color="black")+ #if you want to add a black line of average
  #stat_summary(fun.data = 'mean_sdl',fun.args = list(mult=1), geom = 'smooth',se=T)+
  theme_classic()+xlab(NULL)+ylab("Average Daily Temperature")+
  theme(legend.position="none")+
  scale_color_manual(values = colorS)+
  geom_hline(yintercept = 30.5,linetype="dashed") #showing the bleaching threshold of 30.5
p

inP<-temp %>% filter(Date>=as.Date("2018-07-01") & Date <= as.Date("2018-09-30")) %>% 
  group_by(Site,day=floor_date(Date,"day")) %>% 
  summarise(DailyAvg=mean(Temp),
            SD=sd(Temp)) %>%
  ggplot(aes(day,DailyAvg))+
  geom_line(aes(group=Site,color=Site))+
  #scale_color_manual(values = colorS)+
  #stat_summary(fun=mean,geom = 'line',color="black",size=1)+
  theme_classic()+
  theme(legend.position ="none" )+xlab(NULL)+ylab(NULL)+
  scale_color_manual(values = colorS[-6])+
  geom_hline(yintercept = 30.5,linetype="dashed")
inP

#plot the timeseries with the inset
library(cowplot)
library(grid)
quartz()
ggdraw()+
  draw_plot(p)+
  draw_plot(inP,.09, .05, .4, .3)
#draw_plot(rectGrob(x=.3,y=.86,width=.03,height=.18,gp = gpar(fill = "skyblue2", alpha = 0.3)))
#####

#####generating site specific temperature characteristics and saving the all in one dataframe
####################
#maximum monthly mean
coords$MaxMonthlyMean<-temp %>% group_by(Site,month=floor_date(Date,"month")) %>%
  summarise(MonthlyMean=mean(Temp)) %>%
  group_by(Site)%>%
  summarize(max(MonthlyMean)) %>% pull(`max(MonthlyMean)`)

#annual mean
coords$AnnualMean<- temp %>% group_by(Site) %>%
  summarise(Mean=mean(Temp)) %>% pull(`Mean`)

#annual range
coords$AnnualRange<-temp %>%group_by(Site) %>%
  summarize(Range=max(Temp)-min(Temp))%>% pull(Range)

#average daily range
coords$AvgDailyRange<-temp %>% group_by(Site,day=floor_date(Date,"day")) %>%
  summarize(dayRange=max(Temp)-min(Temp)) %>%
  group_by(Site) %>%
  summarize(mean(dayRange)) %>% pull(`mean(dayRange)`)

#site survival risk scores
m2$coefficients[1:8] #remember the best site is considered 1
coords$RiskScore<-c(1,1.070873,1.378437,1.308932,1.413481,2.124245,2.067483,2.152316,2.564903)

#summing the number of days above 30.5
#used as metric of cumuluative stress per year in Manzello et al 2007 https://doi.org/10.1016/j.marpolbul.2007.08.009
# They also found days above 30.5C were a useful indicator of bleaching among reefs in the FKRT
HotDays<-aggregate(temp$DH30.5~Site+floor_date(Date,'month')+floor_date(Date,"day"), data = temp, sum)
HotDays$DD30.5<-ifelse(HotDays$`temp$DH30.5`>0,1,0)

#Can do the same for 32C 
HotDays$DH32<-aggregate(temp$DH32~Site+floor_date(Date,'month')+floor_date(Date,"day"), data = temp, sum) %>% pull(`temp$DH32`)
HotDays$DD32<-ifelse(HotDays$DH32>0,1,0)

#Summing the hours and days above certain threshold, later on you see most are these are correlated with eachother
TotalHots<-HotDays %>% group_by(Site) %>%
  summarise(TDD30.5=sum(DD30.5),
            TDH30.5=sum(`temp$DH30.5`),
            TDD32=sum(DD32),
            TDH32=sum(DH32))
coords<-cbind(coords,TotalHots[,2:5]) 

####estimating the autocorrelation in temperature

#Using an "area under the curve" approach. It is hard to decide a relevant time scale to look at thermal predictability,
#So instead, plot autocorrelation across numerous scales (lags). 
#Autocor tends to declines as the lag gets bigger (harder to predict increasingly more future temps from current temps)
#Sites that are less predictable will have sharper declines in autocor as lag increases 
#Quantifying this predictablilty as the sum of autocor before autocor drops below 0 (i.e. no information about future temps from current temps)
#acf() function calculates autocor over any lag you want so use this to make a curve and find the area under it that is positive autocor
#focusing this on the time when temp is hottest and most variable (June through Sep- see time series fig)

#use this as a test to see what acf() produces. First filter by dates and site of interest
ts<-temp %>% filter(Date>=as.Date("2018-07-01") & Date <= as.Date("2018-09-30")) %>% filter(Site=="Bahia Honda")
acf(ts[,3],lag.max =2000,plot=T,main="Autocor over Lags in FL Summer")

sites<-c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","EDR","Maryland Shoals","Bahia Honda") #looe key is missing data for the summer
a<-matrix(ncol=2,c("Site","Predictability"))
for (i in sites) {
  ts<-
    temp %>% filter(Date>=as.Date("2018-07-01") & Date <= as.Date("2018-09-30")) %>% #isolates time series to only FL summer, on the time series, these dates are after the summmer rise and before fall fall
    #group_by(Site) %>% 
    #summarise(WeeklyAvg=mean(Temp)) %>%
    filter(Site==i)
  lag<-acf(ts[,3],lag.max =2000,plot=F,main=i) #2000 is roughly 2.5 months,ACF reaches 0 well before this lag.max  
  auc<-sum(lag$acf[1:((which(lag$acf<0)[1])-1)]) #to find area under curve to first cross of x axis
  a<-rbind(a,c(i,auc))
  #od<-lm(lag$acf~lag$lag)  #to find slope 
  #a<-rbind(a,c(i,mod$coefficients[2])) #to find slope
}

plot(a[,2]) #seems like predictability goes down with increasing risk
a<-as.data.frame(a[-1,])
a$V2
coords$Predictability<-c(209.141, 195.195, 214.138, 134.525, 138.193,NA, 186.321,182.886,105.005)
#####

#######2b. SERC Reef water quality data
################

serc<-read.csv("EnvironmentalData/SERC_WQ_AUG2020_uM.csv",header=T, stringsAsFactors = T)
summary(serc) #its a bit messy, especially because sites are not labeled the same way as the outplant experiment names
unique(serc$SITE)
unique(serc$STATION)
serc$STATION<-as.factor(serc$STATION)
serc$STATION=factor(serc$STATION,levels=c("256","259","263","267","270","273","403","276","280"))
#American Shoal=Dave's Ledge, Western Sambo #2= W. Sambo, Western Sambo = Marker 32
#make those changes and others to make downstream steps easier
serc$SITE<-ifelse(serc$SITE=="Bahia Honda Offshore","Bahia Honda",
                  ifelse(serc$SITE=="Big Pine Shoal","Big Pine",
                         ifelse(serc$SITE=="Looe Key","Looe Key",
                                ifelse(serc$SITE=="American Shoal","Dave's Ledge",
                                       ifelse(serc$SITE=="Maryla Shoal","Maryland Shoals",
                                              ifelse(serc$SITE=="Western Sambo #2","W. Sambo",
                                                     ifelse(serc$SITE=="Western Sambo","Marker 32",
                                                            ifelse(serc$SITE=="Eastern Dry Rocks","EDR",
                                                                   ifelse(serc$SITE=="Eastern Sambo Offshore","E. Sambo",NA)))))))))
#Reorder sites from best to worst survival
serc<-serc %>% mutate(SITE = factor(SITE, levels=c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","Looe Key","EDR","Maryland Shoals","Bahia Honda")))

serc$DATE<-lubridate::mdy(serc$DATE)
serc<-dplyr::arrange(serc,DATE)
sercH<-filter(serc,DATE<=as.Date("2019-6-1")) #historical serc data is everything from start of SERC to end of exp.
#pulling out only outplant time period
sercC<-filter(serc,DATE>=as.Date("2018-3-1")& DATE<=as.Date("2019-6-1")) #contemporary SERC data from only experimental period

#individual ANOVAs for each water quality  parameter in both historical and contemporary sets
#change dataset as needed           Cont  Hist
summary(aov(lm(NOX.B~SITE,sercH))) #.813 #0.0193 *
summary(aov(lm(NO2.B~SITE,sercH))) #.723 #0.147
summary(aov(lm(NO3_B~SITE,sercH))) #.741 #0.0145*
summary(aov(lm(NH4.B~SITE,sercH))) # .579 #0.414
summary(aov(lm(TN.B~SITE,sercH))) #0.37 #0.597
summary(aov(lm(DIN.B~SITE,sercH))) # 0.444 #0.346
summary(aov(lm(TON.B~SITE,sercH))) # 0.379 #0.63
summary(aov(lm(TP.B~SITE,sercH))) #0.991 #0.99
summary(aov(lm(SRP.B~SITE,sercH))) #.267 #0.798
summary(aov(lm(TURB.B~SITE,sercH))) # 0.0824 #0.697
summary(aov(lm(TOC.B~SITE,sercH))) # 0.54 #0.488
summary(aov(lm(SiO2.B~SITE,sercH))) #.856 #0.0405 *
summary(aov(lm(SAL.B~SITE,sercH))) #.993 #0.955
summary(aov(lm(DO.B~SITE,sercH))) #.77 #0.997
summary(aov(lm(Kd~SITE,sercH))) #.101 #0.661
summary(aov(lm(X.SAT_B~SITE,sercH))) #.768 #0.997
summary(aov(lm(Temp~Site,temp))) #2e-16

model<-aov(Temp~Site,temp) 
summary(model) 
TukeyHSD(x=model, 'SITE', conf.level=0.95)

#Summarizing historical and contemporary averages and ranges for all parameters at each site
sercH.sum <-
  sercH %>% #just change to sercC to sum contemporary values
  group_by(SITE) %>%
  summarize(mean_NO3.B=mean(NO3_B, na.rm=T),
            range_NO3.B=(max(NO3_B, na.rm = T)-(min(NO3_B,na.rm = T))),
            mean_NO2.B=mean(NO2.B, na.rm=T),
            range_NO2.B=(max(NO2.B, na.rm = T)-(min(NO2.B,na.rm = T))),
            mean_NH4.B=mean(NH4.B, na.rm=T),
            range_NH4.B=(max(NH4.B, na.rm = T)-(min(NH4.B,na.rm = T))),
            mean_TN.B=mean(TN.B, na.rm=T),
            range_TN.B=(max(TN.B, na.rm = T)-(min(TN.B,na.rm = T))),
            mean_TP.B=mean(TP.B, na.rm=T),
            range_TP.B=(max(TP.B, na.rm = T)-(min(TP.B,na.rm = T))),
            mean_TON.B=mean(TON.B, na.rm=T),
            range_TON.B=(max(TON.B, na.rm = T)-(min(TON.B,na.rm = T))),
            mean_TOC.B=mean(TOC.B, na.rm=T),
            range_TOC.B=(max(TOC.B, na.rm = T)-(min(TOC.B,na.rm = T))),
            mean_Kd=mean(Kd, na.rm=T), #light attenuation coefficient 
            range_Kd=(max(Kd, na.rm = T)-(min(Kd,na.rm = T))),
            mean_SiO2.B=mean(SiO2.B, na.rm=T),
            range_SiO2.B=(max(SiO2.B, na.rm = T)-(min(SiO2.B,na.rm = T))),
            mean_TURB.B=mean(TURB.B, na.rm=T),
            range_TURB.B=(max(TURB.B, na.rm = T)-(min(TURB.B,na.rm = T))),
            mean_DO.B=mean(DO.B, na.rm=T),
            range_DO.B=(max(DO.B, na.rm = T)-(min(DO.B,na.rm = T))),
            mean_X.SAT_B=mean(X.SAT_B, na.rm=T),
            range_X.SAT_B=(max(X.SAT_B, na.rm = T)-(min(X.SAT_B,na.rm = T))))

####SAVE THIS IF REVIEWERS ASK - otherwise continue on with analysis using overall range in historical data
#for historical data instead of overall range, calculating average annual ranges 
sercH.sum <-
  sercH %>%
  group_by(SITE) %>%
  summarize(mean_NO3.B=mean(NO3_B, na.rm=T),mean_NO2.B=mean(NO2.B, na.rm=T),
            mean_NH4.B=mean(NH4.B, na.rm=T),mean_TN.B=mean(TN.B, na.rm=T),
            mean_TP.B=mean(TP.B, na.rm=T), mean_TON.B=mean(TON.B, na.rm=T),
            mean_TOC.B=mean(TOC.B, na.rm=T),mean_Kd=mean(Kd, na.rm=T), #light attenuation coefficient 
            mean_SiO2.B=mean(SiO2.B, na.rm=T),mean_TURB.B=mean(TURB.B, na.rm=T),
            mean_DO.B=mean(DO.B, na.rm=T),mean_X.SAT_B=mean(X.SAT_B, na.rm=T))
#Finding average annual range, this is only needed for historical SERC data since contemporary only covers a singe year
sercH.sum$AvgRange_NO3_B<-aggregate(sercH$NO3_B~SITE+floor_date(DATE,'year'), data = sercH, FUN = function(x) max(x)-min(x)) %>% group_by(SITE) %>% summarize(AvgRange_NO3_B=mean(`sercH$NO3_B`,na.rm = T)) %>% pull(`AvgRange_NO3_B`)
sercH.sum$AvgRange_NO2.B<-aggregate(sercH$NO2.B~SITE+floor_date(DATE,'year'), data = sercH, FUN = function(x) max(x)-min(x)) %>% group_by(SITE) %>% summarize(AvgRange_NO2.B=mean(`sercH$NO2.B`,na.rm = T)) %>% pull(`AvgRange_NO2.B`)
sercH.sum$AvgRange_NH4.B<-aggregate(sercH$NH4.B~SITE+floor_date(DATE,'year'), data = sercH, FUN = function(x) max(x)-min(x)) %>% group_by(SITE) %>% summarize(AvgRange_NH4.B=mean(`sercH$NH4.B`,na.rm = T)) %>% pull(`AvgRange_NH4.B`)
sercH.sum$AvgRange_TN.B<-aggregate(sercH$TN.B~SITE+floor_date(DATE,'year'), data = sercH, FUN = function(x) max(x)-min(x)) %>% group_by(SITE) %>% summarize(AvgRange_TN.B=mean(`sercH$TN.B`,na.rm = T)) %>% pull(`AvgRange_TN.B`)
sercH.sum$AvgRange_TP.B<-aggregate(sercH$TP.B~SITE+floor_date(DATE,'year'), data = sercH, FUN = function(x) max(x)-min(x)) %>% group_by(SITE) %>% summarize(AvgRange_TP.B=mean(`sercH$TP.B`,na.rm = T)) %>% pull(`AvgRange_TP.B`)
sercH.sum$AvgRange_TON.B<-aggregate(sercH$TON.B~SITE+floor_date(DATE,'year'), data = sercH, FUN = function(x) max(x)-min(x)) %>% group_by(SITE) %>% summarize(AvgRange_TON.B=mean(`sercH$TON.B`,na.rm = T)) %>% pull(`AvgRange_TON.B`)
sercH.sum$AvgRange_TOC.B<-aggregate(sercH$TOC.B~SITE+floor_date(DATE,'year'), data = sercH, FUN = function(x) max(x)-min(x)) %>% group_by(SITE) %>% summarize(AvgRange_TOC.B=mean(`sercH$TOC.B`,na.rm = T)) %>% pull(`AvgRange_TOC.B`)
sercH.sum$AvgRange_Kd<-aggregate(sercH$Kd~SITE+floor_date(DATE,'year'), data = sercH, FUN = function(x) max(x)-min(x)) %>% group_by(SITE) %>% summarize(AvgRange_Kd=mean(`sercH$Kd`,na.rm = T)) %>% pull(`AvgRange_Kd`)
sercH.sum$AvgRange_SiO2.B<-aggregate(sercH$SiO2.B~SITE+floor_date(DATE,'year'), data = sercH, FUN = function(x) max(x)-min(x)) %>% group_by(SITE) %>% summarize(AvgRange_SiO2.B=mean(`sercH$SiO2.B`,na.rm = T)) %>% pull(`AvgRange_SiO2.B`)
sercH.sum$AvgRange_TURB.B<-aggregate(sercH$TURB.B~SITE+floor_date(DATE,'year'), data = sercH, FUN = function(x) max(x)-min(x)) %>% group_by(SITE) %>% summarize(AvgRange_TURB.B=mean(`sercH$TURB.B`,na.rm = T)) %>% pull(`AvgRange_TURB.B`)
sercH.sum$AvgRange_DO.B<-aggregate(sercH$DO.B~SITE+floor_date(DATE,'year'), data = sercH, FUN = function(x) max(x)-min(x)) %>% group_by(SITE) %>% summarize(AvgRange_DO.B=mean(`sercH$DO.B`,na.rm = T)) %>% pull(`AvgRange_DO.B`)
sercH.sum$AvgRange_X.SAT_B<-aggregate(sercH$X.SAT_B~SITE+floor_date(DATE,'year'), data = sercH, FUN = function(x) max(x)-min(x)) %>% group_by(SITE) %>% summarize(AvgRange_X.SAT_B=mean(`sercH$X.SAT_B`,na.rm = T)) %>% pull(`AvgRange_X.SAT_B`)

#####

#####2c. Principal components analysis that combines SERC and temp data
#########
#Pick back up here with which ever version of historical SERC summary
envH<-cbind(sercH.sum[,1],coords[,2:13],sercH.sum[,2:25])
envC<-cbind(sercC.sum[,1],coords[,2:13],sercC.sum[,2:25])

#PCA with either historical SERC+thermal & contemporary SERC+thermal
pcaH<-prcomp(envH[-6,c(2:7,9:37)],scale=T) #have to remove LK because thermal characteristics are off due to missing data
pcaC<-prcomp(envC[-6,c(2:7,9:37)],scale=T)

quartz()
autoplot(pcaH,data=envH[-6,],loadings=F,loadings.label=F,colour='SITE',size=4,
         loadings.label.size=3,loadings.label.repel=F, loading.label.color="darkgrey", loadings.colour='lightgrey',
         label.size=3)+
  #geom_text(vjust=1.5,label=envC[-6,1],size=3)+
  #geom_text_repel(label=envC[-6,1])+
  theme_classic()+
  #theme(legend.position = "none")+
  scale_color_manual(values = c("#053061", "#2971B1", "#6AABD0", "#C1DDEB", "#B5B5B5", "#E58267", "#BB2A33", "#67001F"))

########

##Pulling in trait data from the MorphologicalPlasticityAnalysis.R script
grow<-read.csv("/Volumes/GoogleDrive/My Drive/SCHOOL/LAB/CRRPTransplant/AcerMorphologyData_Imputed2.csv")
grow$Tag<-as.factor(grow$Tag)
grow$Array<-as.factor(grow$Array)
grow$Genotype<-as.factor(grow$Genotype)
grow$FragID<-as.factor(grow$FragID)
grow$NurseyGroup<-as.factor(grow$NurseyGroup)
#relevel sites from East to West
grow= grow %>% 
  mutate(Site = factor(Site, levels=c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","Looe Key","EDR","Maryland Shoals","Bahia Honda")))
grow= grow %>%
  mutate(Genotype=factor(Genotype,levels=c("36","1","50","3","44","7","31","13","62","41")))

#Adding up breakage but allowing some potential error which is calculated from the CV for the average size of corals at each time point.
grow$T3_Break<-ifelse(grow$T3_TLE-grow$T0_TLE<(-.2),1,0) 
grow$T6_Break<-ifelse(grow$T6_TLE-grow$T3_TLE<(-.25),1,0)
grow$T9_Break<-ifelse(grow$T9_TLE-grow$T6_TLE<(-.35),1,0)
grow$T12_Break<-ifelse(grow$T12_TLE-grow$T9_TLE<(-.52),1,0)
grow$CulumativeBreaks<-rowSums(grow[,c("T3_Break","T6_Break","T9_Break","T12_Break")],na.rm = T)
grow$CulumativeBreaks<-ifelse(grow$T3_Status=="M",NA,grow$CulumativeBreaks) #if a coral is missing from the start it should not be included


#3a. Start of preparing a dataframe for Bayesian model
#this is a process because instead of taking a single value of survival, growth, temp, water quality variable for every site, every site will have 4 values
#this is possible because growth was measure over 4 time periods, SERC has a value for each of those periods, and temp characteristics can be calculated within each period
# the follow is X steps to creating this dataframe for bayesian analysis
#############################
#Step 1: create site averages for change in size in each trait at each time point 
  sSize<-grow
  sSize<-
    sSize %>% 
    mutate(TLE_0.3=(T3_TLE-T0_TLE),
           TLE_3.6=(T6_TLE-T3_TLE),
           TLE_6.9=(T9_TLE-T6_TLE),
           TLE_9.12=(T12_TLE-T9_TLE),
           SA_0.3=(T3_SA-T0_SA),
           SA_3.6=(T6_SA-T3_SA),
           SA_6.9=(T9_SA-T6_SA),
           SA_9.12=(T12_SA-T9_SA),
           V_0.3=(T3_V-T0_V),
           V_3.6=(T6_V-T3_V),
           V_6.9=(T9_V-T6_V),
           V_9.12=(T12_V-T9_V),
           Vinter_0.3=(T3_Vinter-T0_Vinter),
           Vinter_3.6=(T6_Vinter-T3_Vinter),
           Vinter_6.9=(T9_Vinter-T6_Vinter),
           Vinter_9.12=(T12_Vinter-T9_Vinter))
  
  #make sure there are NAs in all following time points if dead (some instances where we measured a dead coral)
  sSize[which(sSize$T3_Status!="A"),c(58:73)]<-NA #if growth is negative or a coral is dead, replace all following GRs with NA
  sSize[which(sSize$T6_Status!="A"),c(59:61,63:65,67:69,71:73)]<-NA
  sSize[which(sSize$T9_Status!="A"),c(60,61,64,65,68,69,72,73)]<-NA
  sSize[which(sSize$T12_Status!="A"),c(61,65,69,73)]<-NA
  
  sSize<-sSize %>% group_by(Site) %>% summarize_at(c(names(.)[58:73]),mean, na.rm=T)
  #should have a dataframe that is 9 rows with 17 columns that contain the change in size for each trait/time
  sSize2<-sSize %>% pivot_longer(cols=2:17,names_to = c('trait','time'),names_sep='_',values_to="ChangeInSize")
  #reorder the two dfs to match
  sSize2<- sSize2 %>% pivot_wider(names_from = "trait",values_from="ChangeInSize")

#Step 2: adding time-point specific Risk Score values
  #Essentially we will preform the coxph model as if each time point is the last
  library(survival)
  library(coxme)
  library(survminer)
  surv<-grow
  str(surv)
  
  #This will require individual columns for each time describing status and time of death
  #in status 1= dead 0=alive NA=missing
  #in TOD, time point of death will be listed, if still alive, we will indicate that ramet survived to the next "month"
  surv<- surv %>%
    mutate(T3_Dead = case_when(
      T3_Status=="A" ~ "0",
      T3_Status=="D" ~ "1", 
      T3_Status=="M" ~ "NA" ),
      T3_TOD= case_when(
        T3_Status=="A" ~ "4",
        T3_Status=="D" ~ "3", 
        T3_Status=="M" ~ "NA"),
      T6_Dead = case_when(
        T6_Status=="A" ~ "0",
        T6_Status=="D" ~ "1", 
        T6_Status=="M" ~ "NA" ),
      T6_TOD= case_when(
        T6_Status=="A" ~ "7",
        T3_Status=="D" & T6_Status=="D" ~ "3",
        T3_Status=="A" & T6_Status=="D" ~ "6",
        T3_Status=="A" & T6_Status=="M" ~ "NA",
        T6_Status=="M" ~ "NA"),
      T9_Dead = case_when(
        T9_Status=="A" ~ "0",
        T9_Status=="D" ~ "1", 
        T9_Status=="M" ~ "NA" ),
      T9_TOD= case_when(
        T9_Status=="A" ~ "10",
        T3_Status=="D" & T6_Status=="D" & T9_Status=="D" ~ "3",
        T3_Status=="A" & T6_Status=="D" & T9_Status=="D" ~ "6",
        T3_Status=="A" & T6_Status=="A" & T9_Status=="D" ~ "9",
        T3_Status=="A" & T6_Status=="M" & T9_Status=="M" ~ "NA",
        T3_Status=="A" & T6_Status=="A" & T9_Status=="M" ~ "NA",
        T9_Status=="M" ~ "NA")
      )
  
  surv$T3_TOD<-as.numeric(surv$T3_TOD)
  surv$T6_TOD<-as.numeric(surv$T6_TOD)
  surv$T9_TOD<-as.numeric(surv$T9_TOD)
  surv$T12_TOD<-as.numeric(surv$T12_TOD)
  
  surv$T3_Dead<-as.numeric(surv$T3_Dead)
  surv$T6_Dead<-as.numeric(surv$T6_Dead)
  surv$T9_Dead<-as.numeric(surv$T9_Dead)
  surv$T12_Dead<-as.numeric(surv$T12_Dead)
  
  surv$Site<-relevel(surv$Site,ref="E. Sambo")## picked site with best survival as reference (want to know how much worse others are)
  surv$Genotype<-relevel(surv$Genotype,ref="36")## picked genet with best survival as reference
  
  scoreT3<-surv[complete.cases(surv$T3_TOD),]
  scoreT6<-surv[complete.cases(surv$T6_TOD),]
  scoreT9<-surv[complete.cases(surv$T9_TOD),]
  scoreT12<-surv[complete.cases(surv$T12_TOD),]
  
  #creates a object to be used in the model
  died3=Surv(scoreT3$T3_TOD, scoreT3$T3_Dead)
  died6=Surv(scoreT6$T6_TOD, scoreT6$T6_Dead)
  died9=Surv(scoreT9$T9_TOD, scoreT9$T9_Dead)
  died12=Surv(scoreT12$T12_TOD, scoreT12$T12_Dead)
  
  #running the models
  control<-coxme.control(iter.max=30)
  mT3<-coxme(died3~Site+(1|CulumativeBreaks), data=scoreT3) #model did not coverge so can't trust coefficients
  mT6<-coxme(died6~Site+Genotype+(1|CulumativeBreaks), data=scoreT6)
  mT9<-coxme(died9~Site+Genotype+(1|CulumativeBreaks), data=scoreT9)
  mT12<-coxme(died12~Site+Genotype+(1|CulumativeBreaks), data=scoreT12)
  
  #So now we need to pull the coefficients from each model to describe site risk score
  mT3$coefficients[1:8] #these coefficients are quite wild and is most likely because model didn't converge
  
  #For t3 only, I am going to take the proportion surviving at each site (0<x<1). 
  #then to set the E. Sambo (the reference at other time points) to 1 and to make worse surviving sites have a higher value, subtract proportion surviving from (1+E.Sambo value)
  #now reference is 1 and higher values = greater risk and values are on a similar scale with other time points
  S=survfit(died3~Site, data=scoreT3)
  S # ES .966; M32 .931; WS 1; BP 1, DL .931, LK .833, EDR .852; MS .967; BH .733  
  1.966+(c(.966,.931,1,1,.931,.833,.852,.967,.733)*-1) 
  #1.000 1.035 0.966 0.966 1.035 1.133 1.114 0.999 1.233
  
  #risk scores over T3 to T6
  mT6$coefficients[1:8]
  #SiteMarker 32        SiteW. Sambo        SiteBig Pine    SiteDave's Ledge        SiteLooe Key    SiteEDR SiteMaryland Shoals     SiteBahia Honda
  #   0.8764844           0.6484101           0.9694663           1.0815462           1.7422409     1.8215727           1.2932383           2.3766294
  
  #risk scores over T6 to T9
  mT9$coefficients[1:8]
  #SiteMarker 32        SiteW. Sambo        SiteBig Pine    SiteDave's Ledge        SiteLooe Key          SiteEDR   SiteMaryland Shoals SiteBahia Honda 
  #   1.005173            1.087909            1.289528            1.084357            2.136075            2.063799            2.028440      2.579155 
  
  #risk scores over T9 to T12
  mT12$coefficients[1:8]
  #SiteMarker 32        SiteW. Sambo        SiteBig Pine    SiteDave's Ledge        SiteLooe Key             SiteEDR SiteMaryland Shoals     SiteBahia Honda 
  #   1.070873            1.378437            1.308932            1.413481            2.124245            2.067483            2.152316            2.564903 
  
  #now put it all together. I am sure there is a better way to do it but I like the control when doing it manually
  rs<-data.frame(rep(c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","Looe Key","EDR","Maryland Shoals","Bahia Honda"),4),#making sure sites align
             c(rep("0.3",9),rep("3.6",9),rep("6.9",9),rep("9.12",9)),#adding timepoints to keep track
             c(c(1.000, 1.035, 0.966, 0.966, 1.035, 1.133, 1.114, 0.999, 1.233),#T3 scores
             c(1,0.8764844 ,0.6484101,0.9694663,1.0815462,1.7422409,1.8215727,1.2932383,2.3766294), #t6 scores
             c(1,1.005173,1.087909,1.289528,1.084357,2.136075, 2.063799,2.028440,2.579155), #T9 scores
             c(1,1.070873,1.378437,1.308932,1.413481,2.124245,2.067483,2.152316,2.564903))) 
  names(rs)<-c("SITE","TimePoint","RiskScore")

#Step 3: add time-point specific SERC values to the size DF - fairly easy because there is only one value measured by SERC over each time interval
  serc3<- sercC %>%
    mutate(TimePoint = case_when(
      DATE>=as.Date("2018-5-1") & DATE<=as.Date("2018-6-10") ~ "0.3",#Apr to Jul
      DATE>=as.Date("2018-9-1") & DATE<=as.Date("2018-9-30") ~ "3.6", #Jul to Oct
      DATE>=as.Date("2018-11-1") & DATE<=as.Date("2018-12-10") ~ "6.9", #Oct to Jan
      DATE>=as.Date("2019-5-1") & DATE<=as.Date("2019-5-30") ~ "9.12")) #Jan to Apr
  serc3<-serc3[complete.cases(serc3$TimePoint),]
  serc3<- serc3 %>% arrange(factor(SITE,levels=c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","Looe Key","EDR","Maryland Shoals","Bahia Honda")))

#Step 4: Calculate time interval specific thermal parameters (pain i the ass)
  #Average over time intervals
  spring<-temp %>% filter(Date>=as.Date("2018-05-01") & Date <= as.Date("2018-07-10"))%>%
    group_by(Site) %>% summarize(MeanTemp=mean(Temp))
  summer<-temp %>% filter(Date>=as.Date("2018-07-11") & Date <= as.Date("2018-10-20"))%>%
    group_by(Site) %>% summarize(MeanTemp=mean(Temp))
  fall<-temp %>% filter(Date>=as.Date("2018-10-24") & Date <= as.Date("2019-01-09"))%>%
    group_by(Site) %>% summarize(MeanTemp=mean(Temp))
  winter<-temp %>% filter(Date>=as.Date("2019-01-10") & Date <= as.Date("2019-04-30"))%>%
    group_by(Site) %>% summarize(MeanTemp=mean(Temp))
  
  #Range
  spring$TempRange<-temp %>% filter(Date>=as.Date("2018-05-01") & Date <= as.Date("2018-07-10"))%>%
    group_by(Site) %>% summarize(Range=max(Temp)-mean(Temp)) %>% pull(`Range`)
  summer$TempRange<-temp %>% filter(Date>=as.Date("2018-07-11") & Date <= as.Date("2018-10-20"))%>%
    group_by(Site) %>% summarize(Range=max(Temp)-mean(Temp)) %>% pull(`Range`)
  fall$TempRange<-temp %>% filter(Date>=as.Date("2018-10-24") & Date <= as.Date("2019-01-09"))%>%
    group_by(Site) %>% summarize(Range=max(Temp)-mean(Temp)) %>% pull(`Range`)
  winter$TempRange<-temp %>% filter(Date>=as.Date("2019-01-10") & Date <= as.Date("2019-04-30"))%>%
    group_by(Site) %>% summarize(Range=max(Temp)-mean(Temp)) %>% pull(`Range`)
  
  #Avg Daily range
  spring$AvgDailyRange<-temp %>% group_by(Site,day=floor_date(Date,"day")) %>% filter(Date>=as.Date("2018-05-01") & Date <= as.Date("2018-07-10")) %>%
    summarize(dayRange=max(Temp)-min(Temp)) %>% group_by(Site) %>%
    summarize(mean(dayRange)) %>% pull(`mean(dayRange)`)
  summer$AvgDailyRange<-temp %>% group_by(Site,day=floor_date(Date,"day")) %>% filter(Date>=as.Date("2018-07-11") & Date <= as.Date("2018-10-20")) %>%
    summarize(dayRange=max(Temp)-min(Temp)) %>% group_by(Site) %>%
    summarize(mean(dayRange)) %>% pull(`mean(dayRange)`)
  fall$AvgDailyRange<-temp %>% group_by(Site,day=floor_date(Date,"day")) %>% filter(Date>=as.Date("2018-10-24") & Date <= as.Date("2019-01-09"))%>%
    summarize(dayRange=max(Temp)-min(Temp)) %>% group_by(Site) %>%
    summarize(mean(dayRange)) %>% pull(`mean(dayRange)`)
  winter$AvgDailyRange<-temp %>% group_by(Site,day=floor_date(Date,"day")) %>% filter(Date>=as.Date("2019-01-10") & Date <= as.Date("2019-04-30")) %>%
    summarize(dayRange=max(Temp)-min(Temp)) %>% group_by(Site) %>%
    summarize(mean(dayRange)) %>% pull(`mean(dayRange)`)
  
  #Days above 30.5
  spring$DD30.5<-temp %>% filter(Date>=as.Date("2018-05-01") & Date <= as.Date("2018-07-10")) %>% group_by(Site,day=floor_date(Date,"day"))%>%
    summarize(DH=mean(DH30.5)) %>% count(DH>0) %>% filter(`DH > 0`==TRUE) %>% pull(n)
  summer$DD30.5<-temp %>% filter(Date>=as.Date("2018-07-11") & Date <= as.Date("2018-10-20")) %>% group_by(Site,day=floor_date(Date,"day"))%>%
    summarize(DH=mean(DH30.5)) %>% count(DH>0) %>% filter(`DH > 0`==TRUE) %>% pull(n)
  fall$DD30.5<-rep(0,9) #no days where the temperature reached above 30.5 in fall or winter, above code won't work so just manually add in zeros
  winter$DD30.5<-rep(0,9)
  
  #predictability ***Don't just ctrl +enter through, you need to adjust dates as necessary*****
  sites<-c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge",
           "Looe Key",
           "EDR","Maryland Shoals","Bahia Honda") #looe key is missing data for the spring/summer, but not fall/winter
  
  #change the date as needed to do spring, summer, fall, winter
  #Spring (no Looe Key): filter(Date>=as.Date("2018-05-01") & Date <= as.Date("2018-07-10")) 
  #Summer (no Looe Key): filter(Date>=as.Date("2018-07-11") & Date <= as.Date("2018-10-20"))
  #Fall: filter(Date>=as.Date("2018-10-24") & Date <= as.Date("2019-01-09"))
  #Winter: filter(Date>=as.Date("2019-01-10") & Date <= as.Date("2019-04-30"))
  a<-matrix(ncol=2,c("Site","Predictability"))
  for (i in sites) {
    ts<-
      temp %>% filter(Date>=as.Date("2019-01-10") & Date <= as.Date("2019-04-30")) %>% 
      filter(Site==i)
    lag<-acf(ts[,3],lag.max =2000,plot=F,main=i) #2000 is roughly 2.5 months,ACF reaches 0 well before this lag.max  
    auc<-sum(lag$acf[1:((which(lag$acf<0)[1])-1)]) #to find area under curve to first cross of x axis
    a<-rbind(a,c(i,auc))
  }
  spring$Predictability<-as.numeric(a[2:9,2])
  summer$Predictability<-as.numeric(a[2:9,2])
  fall$Predictability<-as.numeric(a[2:10,2])
  winter$Predictability<-as.numeric(a[2:10,2])
  #summer and spring both have missing Looe Key data so add in row of all NAs
  x<-c("Looe Key", rep(NA,5))
  spring<-rbind(spring[1:5,],x,spring[6:8,])
  summer<-rbind(summer[1:5,],x,summer[6:8,])
  
  spring$TimePoint<-rep("0.3",9)
  summer$TimePoint<-rep("3.6",9)
  fall$TimePoint<-rep("6.9",9)
  winter$TimePoint<-rep("9.12",9)
  
  tempBay<-rbind(spring,summer,fall,winter)

#Step 5: Make sure survival, size,temp and SERC data are in the same order before combining them
  sSize2<- sSize2 %>% arrange(factor(Site,levels=c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","Looe Key","EDR","Maryland Shoals","Bahia Honda")),
                              factor(time, levels = c("0.3", "3.6", "6.9","9.12")))
  serc3<- serc3 %>% arrange(factor(SITE,levels=c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","Looe Key","EDR","Maryland Shoals","Bahia Honda")),
                            factor(TimePoint, levels = c("0.3", "3.6", "6.9","9.12")))
  tempBay<-tempBay %>% arrange(factor(Site,levels=c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","Looe Key","EDR","Maryland Shoals","Bahia Honda")),
                               factor(TimePoint, levels = c("0.3", "3.6", "6.9","9.12")))
  rs<- rs %>% arrange(factor(SITE,levels=c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","Looe Key","EDR","Maryland Shoals","Bahia Honda")),
                            factor(TimePoint, levels = c("0.3", "3.6", "6.9","9.12")))
  envBay<-cbind(serc3,tempBay,sSize2,rs)
  envBay<-envBay[,-c(1:7,9:12,14,16,18,20,22,24,26,28,30:34,36,38,40:44,47,52,56,57,64,65,70,71)] #removing columns we aren't interested in
  envBay$MeanTemp<-as.numeric(envBay$MeanTemp)
  envBay$TempRange<-as.numeric(envBay$TempRange)
  envBay$AvgDailyRange<-as.numeric(envBay$AvgDailyRange)
  envBay$DD30.5<-as.numeric(envBay$DD30.5)
  envBay$Predictability<-as.numeric(envBay$Predictability)
  #write.csv(envBay,"Bayesian_Environmental_Dataframe.csv",row.names = F)
  
  envBay<-read.csv("/Volumes/GoogleDrive-109655319948019787656/My Drive/SCHOOL/LAB/CRRPTransplant/Bayesian_Environmental_Dataframe.csv",header = T)
  str(envBay)

#Step 6: Remove highly correlated variables
  #considering many variables are highly correlated when looking at the averages, would be best for bayesian model to remove these characters
  library(corrplot)
  res<-cor(envBay[-6,2:27], use = "complete.obs") #use the summarized data set to identify highly correlated variables
  round(res,2)  
  quartz()
  corrplot(res, type="upper",method='number',size=.5)
  
  #it'll be really big, messy plot so might just be better to look at coef in a list which you can make below 
  df<-as.data.frame(res)
  df$WQ1<-row.names(df)
  df<-df %>% pivot_longer(col=1:21,names_to = "WQ2",values_to="coef")
  df<- df %>% filter (coef>0.7 | coef <(-.7))
  #I took this DF into excel and manually went through and identified correlated variables to be removed, prioritized variables with fewer missing values
  #I tended to keep variables if the relationship wasn't obvious chemcially/physically 
  
  #remove: X.IO [21], TN.B [6], TON.B [8], TN.TP [16], NO2.B [4], SiO2.B [12]
  #DIN.B [7], NOX.B [2], X.SAT.B [20], DSIGT [22], NO3_B [3], Si.DIN [19],SRP [10]
  #N.P[17], DIN.TP[18]
  
  #however RJags wont work if there are NAs in the data set
  #if two variables are highly correlated pic the one that has the least missing values
  colSums(is.na(envBay)) 
  #NH4 has no NAs as is highly correlated with all other Nitrogen components
  #SiO2 will remove 9 values but is highly correlatated with T0C.B (no NAs) 
  #N.P will remove 7 values, not highly correlated with other WQ but is also not sig different between sites
  #Kd and DO have identical missing values so if you keep one you might as well keep the other
  #DIN.TP have 2 missing values, it is highly correlated with multiple thermal componets so remove
  
  envBay2<-cbind(envBay[,-c(21,6,8,16,3,4,7,2,20,22,17,18,19,10,12)])
  envBay2 <- subset(envBay2, !is.na(Predictability))
  envBay2 <- subset(envBay2, !is.na(Kd))
  colSums(is.na(envBay2)) #should have no NAs in the final df
  str(envBay2)
  envBay2<-envBay2 %>% relocate(TimePoint, .after=last_col())
#####
  
####3b. Now you can actually start running the Bayesian code

#Standardize covariates
MyStd <- function(x) { (x - mean(x,na.rm=T)) / sd(x,na.rm = T)}

for (i in 2:12) {
  envBay2[i]<-MyStd(envBay2[,i])
}

install.packages("AER")
library(R2jags)
require(lattice)
library(AER)
#JAGS model for Acer morpholog
# Zuur's code for graphs
source("/Volumes/GoogleDrive-109655319948019787656//My Drive/SCHOOL/LAB/RCNtraining/MCMCSupportHighstatV4.R")
# Prepare model data for JAGS
colnames(envBay2[2:12])

#This has to match exact to the variables you wish to include as predictors in the model
X <- model.matrix(~ NH4.B+TP.B+TOC.B+TURB.B+DO.B+Kd+MeanTemp+TempRange+AvgDailyRange+DD30.5+Predictability,data = envBay2)
K <- ncol(X)

### new part for GLMM
Site <- as.numeric(as.factor(envBay2$SITE))  #just give it this order because Looe Key was  removed
Nre<- length(unique(Site))  
### End of New part

envBay2$RiskScore<-as.factor(envBay2$RiskScore) #response variable (change as needed)
#JAGS.data <- list(Y = disease_data$WPL_P_A,
JAGS.data <- list(Y = envBay2$RiskScore, 
                  X = X,
                  K = K,    # New part
                  Site=Site, #Random effects identification 
                  Nre = Nre, # Number of random effects
                  #End of New part
                  N = nrow(envBay2))

# JAGS code for a Negative binomial Generalized linear model from Zuur et al 2016
#Zero-inflated models
sink("AcerMorphPlast_Bayes.txt")
cat("
    model{
    #1A. Priors beta 
    for (i in 1:K) { beta[i] ~ dnorm(0, 0.0001)}
    
    #1B Priors for k (size)
    size ~ dunif(0, 20)
    
    #New part for Priors for random intercepts
    for (i in 1:Nre) {a[i] ~ dnorm(0, tau_Site)}
    
    #Prior for tau_Site
    sigma_Site ~ dgamma(0.001, 0.001)
    tau_Site <- 1 / (sigma_Site * sigma_Site)
    #End of new part
    
    #2. Likelihood
    for (i in 1:N) {
    #This is just how JAGS implements a NB distribution
    #It is based on general mathematical rules for a NB 	
    Y[i] ~  dnegbin(p[i], size)
    p[i] <- size / (size + mu[i])  
    
    log(mu[i]) <- eta[i]
    eta[i]     <- inprod(beta[], X[i,] + a[Site[i]]) # New part 
    
    #3. Discrepancy measures 
    #Pearson residuals
    Exp[i] <- mu[i] 
    Var[i] <- mu[i] + mu[i] * mu[i] / size
    E[i]   <- (Y[i]  - Exp[i]) / sqrt(Var[i])    
    
    #Simulated data with mean/variance taken from the fitted model 
    #See text under block B, below.    
    YNew[i] ~  dnegbin(p[i], size)                      
    
    #Pearson residual for predicted data      
    ENew[i] <- (YNew[i] - Exp[i]) / sqrt(Var[i]) 
    
    #Squared residuals
    D[i]    <- pow(E[i], 2)                      
    DNew[i] <- pow(ENew[i], 2)   
    }          
    
    #Sum of squared Pearson residuals:
    Fit     <- sum(D[1:N])      
    #Sum of squared predicted Pearson residuals:     
    FitNew  <- sum(DNew[1:N])   
    }
    ",fill = TRUE)
sink()

# Initial values & parameters to save
inits  <- function () {
  list(
    beta = rnorm(K, 0, 0.1),
    size = runif(1, 0, 20 ) )  }

params <- c("beta")

# JAGS
G1 <- jags(data       = JAGS.data,
           inits      = inits,
           parameters = params,
           model      = "AcerMorphPlast_Bayes.txt",
           n.thin     = 10,
           n.chains   = 1,
           n.burnin   = 4000,
           n.iter     = 5000)
G2   <- update(G1, n.iter = 100000, n.thin = 10)
#G3   <- update(G2, n.iter = 100000, n.thin = 10)
outGrisk<-G2$BUGSoutput
#outGv<-G2$BUGSoutput
#outGsa<-G2$BUGSoutput
#outGtle
#outGvi 
save(outGtle,outGsa,outGv,outGvi,outGrisk, file = "BayesOutputs.RData")
#Save the model results as you make them so you don't have to keep rerunning the model
#outTLE<-outG
#outSA<-outG
#outV<-outG
#outVi<-outG
outRS <- outG

outG <- outRS

#identify the variables
vars<- c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "beta[6]","beta[7]","beta[8]","beta[9]","beta[10]",
         "beta[11]","beta[12]")

OUT1 <- MyBUGSOutput(outGrisk, vars) 
print(OUT1, digits =3)
summary(G2)

#Check for overdispersion
E    <- G2$BUGSoutput$mean$E
N    <- nrow(envBay2)
p    <- K
sum(E^2) / (N - p)
#Looks good

library(mcmcplots)
#colored caterpillar plot
MyNames <- colnames(X)

# Fancy plot
G2_df=data.frame(variableG2=MyNames,MeanG2=OUT1[,1],Down=OUT1[,3],Up=OUT1[,4])

G2_df$color <- ("#E1EDF3")
G2_df$color[(G2_df$MeanG2 > 0) & (G2_df$Down>0)] <- '#DC6F58'
G2_df$color[(G2_df$MeanG2 < 0) & (G2_df$Up<0)] <- '#549EC9'
G2_df$color[(G2_df$variableG2=="(Intercept)")]<-'#E1EDF3' #just uncoloring the intercept point since it does not make sense to color like other variables
#tiff(file=file.path(home,'output','Figure_2.tif'),height=2000,width=2700,res=300)

quartz()
ggplot(G2_df,aes(x=reorder(variableG2, MeanG2), MeanG2)) +
  geom_errorbar(aes(ymax=Up, ymin=Down), width=0) +
  geom_point(pch=23, size=3, fill=G2_df$color, color="black") +
  coord_flip() +
  theme_grey(base_size=15) +
  guides(colour="none")+
  geom_hline(yintercept=0, linetype="dashed", color="gray") +
  labs(y=expression(paste("Estimated ",gamma," coefficients")), x="")+theme_classic()

#simpler plot
ggplot(G2_df,aes(x=reorder(variableG2, MeanG2), MeanG2)) +
  geom_errorbar(aes(ymax=G2_df$Up, ymin=G2_df$Down), width=0) +
  geom_point(pch=21, size=3, fill=G2_df$color, color="black") +
  coord_flip() +
  theme_grey(base_size=15) +
  guides(colour=FALSE)+
  geom_hline(yintercept=0, linetype="dashed", color="gray") +
  labs(y=expression(paste("Estimated ",gamma," coefficients")), x="")

#Histograms
MyBUGSHist(outG, vars, PanelNames=MyNames)
MyBUGSACF(outG, vars)
MyBUGSHist(outG, vars)
MyBUGSChains(outG, vars) #always get an error

MyBUGSChains(outG, c(uNames("beta", K)), PanelNames = MyNames)
MyBUGSACF(outG, c(uNames("beta", K)), PanelNames = MyNames)
OUT1 <- MyBUGSOutput(outG, c(uNames("beta", K)),VarNames = MyNames)
print(OUT1, digits = 4)

library(mcmcplots)
denplot(outG, "beta", col="black")
MyNames <- colnames(X)

#caterplot(outG, "beta", collapse=TRUE, labels.loc='above', lwd= c(3,6), col="black", pch=16)
quartz()
caterplot(outG, "beta", labels = MyNames, labels.loc='axis', cex.axis=3,lwd= c(1,2), col="black", pch=5, cex.label=.55)+abline(v=0,lty=2)
#caterplot(outG, "beta", labels.loc="axis", labels = newnames, denstrip=TRUE, lwd= c(3,6), col="black", pch=16)


#dev.off()
plot(env2$AUC,env2$RiskScore)

####End of Environmental Analysis