#Morphological Plasticity Analysis

#install.packages("reshape2")
library(tidyr)
library(dplyr)
options(dplyr.summarise.inform=F)
library(ggplot2)
library(lme4)
library(lmerTest)
library(multcomp)
library(pbkrtest)
library(ggfortify)
library(gridExtra)
library(ggpubr)
library(lsmeans)
library(tidyverse)
library(reshape2)
library(naniar)
library(piecewiseSEM)
library(foreach)
library(car)
library(agricolae)

#Start here if you need to impute the data, go to line 224 if not
grow<-read.csv("/Volumes/Google Drive/My Drive/SCHOOL/LAB/CRRPTransplant/Acer_Growth_Data.csv", header = T)
str(grow)

grow$Tag<-as.factor(grow$Tag)
grow$Array<-as.factor(grow$Array)
grow$Genotype<-as.factor(grow$Genotype)
grow$FragID<-as.factor(grow$FragID)
grow$NurseyGroup<-as.factor(grow$NurseyGroup)

#relevel sites from East to West (can choose to order in anyway, eventually they will be ordered by survival)
grow= grow %>% 
  mutate(Site = factor(Site, levels=c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","Looe Key","EDR","Maryland Shoals","Bahia Honda")))
grow= grow %>%
  mutate(Genotype=factor(Genotype,levels=c("36","1","50","3","44","7","31","13","62","41")))

#generating volume of interstitial space for each time point
grow<-
  grow %>% 
  mutate(T0_Vinter=T0_Vcx-T0_V,
         T3_Vinter=T3_Vcx-T3_V,
         T6_Vinter=T6_Vcx-T6_V,
         T9_Vinter=T9_Vcx-T9_V,
         T12_Vinter=T12_Vcx-T12_V)

#how many ramets do we have for each time point
sum(!is.na(grow$T0_V)) #213 (should be 270 because no death has occurred, this is because some 3D models did not work)
sum(!is.na(grow$T3_V)) #212
sum(!is.na(grow$T6_V)) #231
sum(!is.na(grow$T9_V)) #200
sum(!is.na(grow$T12_V)) #204

#A common mistake when phenotyping is recording V as the SA of a mesh with closed holes. 
#If this is the case, V will be greater than SA but this shouldn't happen just based on these coral sizes.
#this should be 0, just as a check
sum(grow$T12_SA<grow$T12_V,na.rm = T) #looks good for all time points

######Imputing some of the missing data for the pre-outplant time point#######
#This is because those 3D models didn't turn out for T0, but by hand TLE can be used to impute data
#First, create a data frame of all time point's data concatenated, this will be used to determine if we use only T0 data or all data for imputing data
tregs<-cbind(data.frame(TLE=c(grow$T0_TLE,grow$T3_TLE,grow$T6_TLE,grow$T9_TLE,grow$T12_TLE)),
             data.frame(SA=c(grow$T0_SA,grow$T3_SA,grow$T6_SA,grow$T9_SA,grow$T12_SA)),
             data.frame(V=c(grow$T0_V,grow$T3_V,grow$T6_V,grow$T9_V,grow$T12_V)),
             data.frame(Vcx=c(grow$T0_Vcx,grow$T3_Vcx,grow$T6_Vcx,grow$T9_Vcx,grow$T12_Vcx)),
             data.frame(SAcx=c(grow$T0_Vcx,grow$T3_SAcx,grow$T6_SAcx,grow$T9_SAcx,grow$T12_SAcx)),
             data.frame(Vinter=c(grow$T0_Vinter,grow$T3_Vinter,grow$T6_Vinter,grow$T9_Vinter,grow$T12_Vinter)))

##Relationship between SA and TLE
  #equation from just T0 data
  summary(lm(grow$T0_SA~grow$T0_TLE)) # r2= .7208,inter=5.3295,slope=3.1943
  plot(grow$T0_TLE,grow$T0_SA)

  #Just manually testing this relationship to see if it produces something that makes sense - for my own sanity
  grow[2,c(1,2,8,10)] #actual values: T0_TLE=9.069	T0_SA=35.7
  3.1943*grow[2,8]+5.3295 #=34.29861 so the regressed value is pretty close to the actual value

  #Now the equation from all time points data
  summary(lm(tregs$SA~tregs$TLE)) #r=0.9387, inter = 4.65103, slope=3.45284
  
  #comparing the two equations to each other
  T0_SA.T0<-(3.1943*grow$T0_TLE+5.3295)
  T0_SA.All<-(3.45284*grow$T0_TLE+4.65103)
  plot(T0_SA.All~T0_SA.T0)+abline(0,1) #if the points align on a 1 to 1 line, then there isn't much difference between the two regressions
  
  #visualizing the difference in the imputation by each regression
  T0_SA.All<-ifelse(is.na(grow$T0_SA),(3.45284*grow$T0_TLE+4.65103),grow$T0_SA)
  plot(x=grow$T0_TLE,y=T0_SA.All, col = ifelse(is.na(grow$T0_SA),"red","black"))+
    points(y=T0_SA.T0,x=grow$T0_TLE,col="green") # black points are actual data, red are imputed using full data regression, green is only from T0 data
  
  #So for SA, use the full data set to create an equation to calculate missing SA values from TLE since data from the two types of regressions are very similar

##Relationship between V and TLE  
  #equation for V from only T0 data
  summary(lm(grow$T0_V~grow$T0_TLE)) #r=.5217, inter=.20034, slope = .92925
  plot(grow$T0_TLE,grow$T0_V)
  
  #how does a polynomial look? doesn't improve the fit any more
  summary(lm(T0_V ~ poly(T0_TLE,2,raw=T),data=grow)) #r=0.5239, x2=-.01688, x=1.25922, inter=-1.2163
  
  #equation for full data
  summary(lm(tregs$V~tregs$TLE)) # r= 0.8114, inter=0.57084, slope=1.05577
  plot(tregs$TLE,tregs$V)#looks like a fairly linear relationship and polynomial doesn't improve r-squared
  
  #comparing the two equations
  T0_V.T0<-ifelse(is.na(grow$T0_V),(.92925*grow$T0_TLE+.20034),NA) 
  T0_V.All<-ifelse(is.na(grow$T0_V),(1.05577*grow$T0_TLE+0.57084),NA)
  plot(T0_V.All,T0_V.T0)+abline(0,1) # pretty consistent offset between the two equations with the full data regression estimating larger V 
  
  T0_V.All<-ifelse(is.na(grow$T0_V),(1.05577*grow$T0_TLE+0.57084),grow$T0_V)
  plot(x=grow$T0_TLE,y=T0_V.All, col = ifelse(is.na(grow$T0_SA),"red","black"))+points(y=T0_V.T0,x=grow$T0_TLE,col="green") #green are from T0 regression
  T0_V.All-T0_V.T0
  
  #for reference in the entire data set
  plot(x=tregs$TLE,y=tregs$V)+
    points(y=T0_V.All,x=grow$T0_TLE, col="red")+
    points(y=T0_V.T0,x=grow$T0_TLE,col="green")
  
  
  #the full data regression overestimates, by ~1.5 cm^3 for the 57 colonies without 3D models
  #for reference, 1cm^3 is 1mL, and all colonies at T0 are >2cm^3 in volume and of the missing colonies, most are generally larger in TLE
  #Considering the disconnect between TLE and V that can occur as colonies get bigger, its best to use the T0 data for imputation with a linear model

##Relationship between V of convex hull and TLE  
  #equation for Volume of convex hull using T0 only
  plot(grow$T0_Vcx~grow$T0_TLE) 
  summary(lm(grow$T0_Vcx~grow$T0_TLE)) #r=0.7516, inter= -6.0113, slope= 2.8592
  #polynomial for T0 data doesn't improve
  summary(lm(T0_Vcx ~ poly(T0_TLE,2,raw=T),data=grow)) #r=0.7505 intercept = -6.30927, x2=-.003554, x=2.928696
  
  #equation for full data
  plot(tregs$TLE,tregs$Vcx) #nonlinear, will use polynomial model
  summary(lm(Vcx~poly(TLE,2,raw=T),data=tregs)) # r= .9229, inter = -10.277617, x=3.956962, x2=0.040357
  
  
  T0_Vc.T0<-ifelse(is.na(grow$T0_Vcx),(2.8592*grow$T0_TLE-6.0113),NA)
  T0_Vc.All<-ifelse(is.na(grow$T0_Vcx),(0.040357*(grow$T0_TLE)^2+3.956962*grow$T0_TLE-10.277617),NA)
  plot(T0_Vc.All,T0_Vc.T0)+abline(0,1) #clear difference due to polynomial model
  
  T0_Vc.All<-ifelse(is.na(grow$T0_Vcx),(0.040357*(grow$T0_TLE)^2+3.956962*grow$T0_TLE-10.277617),grow$T0_Vcx)
  plot(x=grow$T0_TLE,y=T0_Vc.All, col = ifelse(is.na(grow$T0_Vcx),"red","black"))+points(y=T0_Vc.T0,x=grow$T0_TLE,col="green") #green points are from T0 regression
  
  #Use the linear regression from T0 data, full data set is too influenced by larger data points

#Relationship between SA of convex hull and TLE
  #SAcx from T0 data only (SAcx will be used for invariant traits but only for T12)
  plot(grow$T0_SAcx~grow$T0_TLE)
  summary(lm(grow$T0_SAcx~grow$T0_TLE)) #r=0.7531, inter= 6.4704, slope= 4.0134
    
  #using full data set
  plot(tregs$TLE,tregs$SAcx)
  summary(lm(tregs$SAcx~tregs$TLE)) #r=0.9102 , inter=3.76496, slope= 4.66994
  
  #comparing the two sets of imputed data
  T0_SAc.T0<-ifelse(is.na(grow$T0_SAcx),(4.0134*grow$T0_TLE+6.4704),NA)
  T0_SAc.All<-ifelse(is.na(grow$T0_SAcx),(4.66994*grow$T0_TLE+3.76496),NA)
  plot(T0_SAc.All,T0_SAc.T0)+abline(0,1)#seems to overestimate values when using full data but is more aligned at using the polynomial at least when coral are small
  
  T0_SAc.All<-ifelse(is.na(grow$T0_SAcx),(4.66994*grow$T0_TLE+3.76496),grow$T0_SAcx)
  plot(x=grow$T0_TLE,y=T0_SAc.All, col = ifelse(is.na(grow$T0_SAcx),"red","black"))+points(y=T0_SAc.T0,x=grow$T0_TLE,col="green")
  
  #The two equations give fairly similar imputations when colonies are small but there is some differentiation as they get bigger
  #best to go with the T0 linear regression

####Applying those equations to fill in the data
####
grow$T0_SA<-ifelse(is.na(grow$T0_SA),(3.1943*grow$T0_TLE+5.3295),grow$T0_SA)
grow$T0_V<-ifelse(is.na(grow$T0_V),(.92925*grow$T0_TLE+.20034),grow$T0_V)
grow$T0_Vcx<-ifelse(is.na(grow$T0_Vcx),(2.8592*grow$T0_TLE-6.0113),grow$T0_Vcx)
grow$T0_SAcx<-ifelse(is.na(grow$T0_SAcx),(4.0134*grow$T0_TLE+6.4704),grow$T0_SAcx)
sum(is.na(grow$T0_SA)) #now there should be no missing values for that inital T0 time point!

#regenerating volume of interstitial space to update the T0 data
grow<-
  grow %>% 
  mutate(T0_Vinter=T0_Vcx-T0_V)

##Can we impute for some of the T3 corals since there is a large chunk of missing models for EDR 1-20
  summary(lm(grow$T3_SA~grow$T3_TLE)) # r2= .908,inter=9.30429,slope=2.90982
  plot(grow$T3_TLE,grow$T3_SA)
  
  summary(lm(filter(grow,Site=="EDR")$T3_SA~filter(grow,Site=="EDR")$T3_TLE)) # r2=.9729, inter=4.2166, slope=2.7014
  SA_edr<-ifelse(is.na(grow$T3_SA),(2.7014*grow$T3_TLE+4.2166),NA)
  SA.All<-ifelse(is.na(grow$T3_SA),(2.90982*grow$T3_TLE+9.30429),NA)
  plot(filter(grow,Site=="EDR")$T3_TLE,filter(grow,Site=="EDR")$T3_SA)+points(y=SA_edr,x=grow$T3_TLE,col="green")+points(y=SA.All,x=grow$T3_TLE,col="red")
  #using the EDR based values
  
  #Volume
  summary(lm(grow$T3_V~grow$T3_TLE)) # r2= .7855,slope=0.7094,inter=2.6366 
  plot(grow$T3_TLE,grow$T3_V)
  
  summary(lm(filter(grow,Site=="EDR")$T3_V~filter(grow,Site=="EDR")$T3_TLE)) # r2=.9729, slope=.79728, inter=2.212866
  V_edr<-ifelse(is.na(grow$T3_V),(.79728*grow$T3_TLE+2.212866),NA)
  V.All<-ifelse(is.na(grow$T3_V),(.7094*grow$T3_TLE+2.6366),NA)
  plot(filter(grow,Site=="EDR")$T3_TLE,filter(grow,Site=="EDR")$T3_V)+points(y=V_edr,x=grow$T3_TLE,col="green")+points(y=V.All,x=grow$T3_TLE,col="red")
  
  #Vcx
  summary(lm(grow$T3_Vcx~grow$T3_TLE)) # r2= .8333,inter=-5.4428,slope=3.7225
  plot(grow$T3_TLE,grow$T3_V)
  summary(lm(filter(grow,Site=="EDR")$T3_Vcx~filter(grow,Site=="EDR")$T3_TLE)) # r2=.949, inter=-19.6043, slope=4.1235
  summary(lm(T3_Vcx~poly(T3_TLE,2,raw=T),data=filter(grow,Site=="EDR"))) #r2=.9486, inter=-6.94964,x=2.57763, x2=.03640
  
  #Vcx_edr<-ifelse(is.na(grow$T3_Vcx),(4.1235*grow$T3_TLE-19.6043),NA)
  Vcx_edr<-ifelse(is.na(grow$T3_Vcx),(.03640*(grow$T3_TLE)^2+2.57763*grow$T3_TLE-6.94964),NA)
  Vcx.All<-ifelse(is.na(grow$T3_Vcx),(3.7225*grow$T3_TLE-5.4419),NA)
  plot(filter(grow,Site=="EDR")$T3_TLE,filter(grow,Site=="EDR")$T3_Vcx)+points(y=Vcx_edr,x=grow$T3_TLE,col="green")+points(y=Vcx.All,x=grow$T3_TLE,col="red")

  #SAcx
  plot(filter(grow,Site=="EDR")$T3_SAcx~filter(grow,Site=="EDR")$T3_TLE)
  summary(lm(filter(grow,Site=="EDR")$T3_SAcx~filter(grow,Site=="EDR")$T3_TLE)) #r= .9723, inter= 7.3565, slope=3.8744
  SAcx_edr<-ifelse(is.na(grow$T3_SAcx),(3.8744*grow$T3_TLE+7.3565),NA)
  plot(filter(grow,Site=="EDR")$T3_TLE,filter(grow,Site=="EDR")$T3_SAcx)+points(y=SAcx_edr,x=grow$T3_TLE,col="green")
  
####Applying those equations to fill in the data
####
grow$T3_SA<-ifelse(is.na(grow$T3_SA)&grow$Site=="EDR",(2.7014*grow$T3_TLE+4.2166),grow$T3_SA)
grow$T3_V<-ifelse(is.na(grow$T3_V)&grow$Site=="EDR",(.79728*grow$T3_TLE+2.212866),grow$T3_V)
grow$T3_Vcx<-ifelse(is.na(grow$T3_Vcx)&grow$Site=="EDR",(.03640*(grow$T3_TLE)^2+2.57763*grow$T3_TLE-6.94964),grow$T3_Vcx)
grow$T3_SAcx<-ifelse(is.na(grow$T3_SAcx)&grow$Site=="EDR",(3.8744*grow$T3_TLE+7.3565),grow$T3_SAcx)

#recalculating volume of interstitial space to update the T3 data
grow<-
  grow %>% 
  mutate(T0_Vinter=T0_Vcx-T0_V)
write.csv(grow,"AcerMorphologyData_Imputed2.csv",row.names = FALSE)
######End of imputation  step########

#Start here if you don't have to change the imputation step and are using the Imputed2 data set from github
grow<-read.csv("/Volumes/Google Drive/My Drive/SCHOOL/LAB/CRRPTransplant/AcerMorphologyData_Imputed2.csv")

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

#######Survival Analysis#####
#using coxPH ratios
#install.packages("nationalparkcolors")
library(survival) # for survival analysis 
library(vegan) # for plotting ellipses in principal coordinates analysis plots
library(corr) # for correlations, not available for R v 3.6.3
library(MCMCglmm) # for mcmcglmm stats
library(nlme) #for lme
library(MASS) # for stepAIC
library(coxme) # for mixed effects cox model
library(ggridges) # for ridge plots
library(corrplot) # for correlations
library(summarytools) # for dfSummary
library(RColorBrewer)
library(survminer)
surv<-grow
str(surv)

surv$T12_TOD<-as.numeric(surv$T12_TOD)

####Making a survival object where TOD is time of death event, and Dead is whether or not a fragment was dead
score<-surv[complete.cases(surv$T12_TOD),] #removing ramets Missing from outplant 

score$Site<-relevel(score$Site,ref="E. Sambo")## picked site with best survival as reference (want to know how much worse others are)
score$Genotype<-relevel(score$Genotype,ref="36")## picked genet with best survival as reference

died=Surv(score$T12_TOD, score$T12_Dead)

####stats
scurveS=survfit(died~Site, data=score)
scurveS
#n events median 0.95LCL 0.95UCL
#Site=E. Sambo        28      1     NA      NA      NA#96.4
#Site=Marker 32       27      2     NA      NA      NA #92.6
#Site=W. Sambo        28      4     NA      NA      NA #85.7
#Site=Big Pine        27      4     NA      NA      NA #85.2
#Site=Dave's Ledge    26      4     NA      NA      NA #84.6
#Site=Looe Key        25      6     NA      NA      NA #80
#Site=EDR             25      6     NA      NA      NA #80
#Site=Maryland Shoals 30      9     NA      NA      NA #70
#Site=Bahia Honda     30     12     NA       9      NA #60
score$Site<-factor(score$Site,levels=c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","Looe Key","EDR","Maryland Shoals","Bahia Honda"))## picked site with best survival as reference (want to know how much worse others are)

scurveG<- survfit(died~Genotype,data=score)
scurveG

#n events median 0.95LCL 0.95UCL
#Genotype=36 25      1     NA      NA      NA #96 
#Genotype=1  20      2     NA      NA      NA #90
#Genotype=50 26      3     NA      NA      NA #88.5
#Genotype=3  26      3     NA      NA      NA #88.5
#Genotype=44 25      3     NA      NA      NA #88
#Genotype=7  23      4     NA      NA      NA #82.6
#Genotype=31 26      5     NA      NA      NA #80.8 
#Genotype=13 25      8     NA      NA      NA #68
#Genotype=62 24      9     NA       9      NA #62.5
#Genotype=41 26     10     NA       9      NA #61.5
score$Genotype<-factor(score$Genotype,levels=c("36","1","50","3","44","7","31","13","62","41"))## picked genet with best survival as reference

#to see which genotypes survive best at each site
scurveGxE<- survfit(died~Site+Genotype,data=score)
scurveGxE
##- identify model that accounts for full effects and tests for potential interaction
control<-coxme.control(iter.max=30)
m1<-coxme(died~Site*Genotype+(1|CulumativeBreaks), data=score,control=control)
m1$loglik
-2*-204.8711  #409.7422
m2<-coxme(died~Site+Genotype+(1|CulumativeBreaks), data=score)
m2$loglik
-2*-234.114 #468.228

468.228-409.7422 # large improvement in log likelihood with interaction (want lower value), even after accounting for breakage

#but cannot trust coefficients bc of lack of convergence. Will need to use additive model
m2
anova(m1)
anova(m2)


m3<-coxme(died~Site+(1|CulumativeBreaks), data=score)
m3$loglik
-2*-245.8340 #491.668

491.668-468.228 #23.44 - including genotype in a model containing site improves model

m4<-coxme(died~Genotype+(1|CulumativeBreaks), data=score)
m4$loglik
-2*-241.8249   #483.65

483.65-468.228 #15.422 - including site in a model containing genotype improves the model


# Visualize
#++++++++++++++++++++++++++++++++++++

#brewer.pal(10,"Paired")
#geno<-c("36","50","44","7","3","31","1","13","41","62")
colorG <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(10)
colorS <- colorRampPalette(rev(brewer.pal(n = 11, name ="RdBu")))(9)
#The middle colors are hard to see so adjusting the 5th color ro a darker grey
colorS <- c("#053061", "#2971B1" ,"#6AABD0" ,"#C1DDEB", "#B5B5B5", "#FACDB5" ,"#E58267", "#BB2A33", "#67001F")

quartz()
ggsurvplot(scurveG, data = score,
           palette=colorG,fun = "pct", conf.int = FALSE,legend=c(.136,.361),legend.title="Genotype",
           legend.labs=c("36","1","50","3","44","7","31","13","62","41"),
           ggtheme = theme_bw(),ylim=c(30,100))

quartz()
ggsurvplot(scurveS, data = score,
           palette=colorS,fun = "pct", conf.int = FALSE,legend=c(.2,.361),legend.title="Site",
           legend.labs=c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","Looe Key","EDR","Maryland Shoals","Bahia Honda"),
           ggtheme = theme_bw(),ylim=c(30,100))

### Now, use correlogram to show correlations among genet performance across sites
scurve <- survfit(died~Site+Genotype, data=score)

res<-data.frame(summary(scurve)$table)
site<-str_split(rownames(res),",",simplify=TRUE)[,1]
site<-gsub('Site=','',site)
genet<-str_split(rownames(res),",",simplify=TRUE)[,2]
genet<-gsub('Genotype=','',genet)
df<-data.frame(res,site,genet,stringsAsFactors=TRUE)

df$n.surv<-df$n.start-df$events
df$PercSurv<-df$n.surv/df$n.start
rownames(df)<-seq_len(nrow(df))
df.sort<-df[order(df$genet),]
df.short<-df.sort[,c(10,11,13)]
df.short$site = factor(df.short$site, levels=c("E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","Looe Key","EDR","Maryland Shoals","Bahia Honda"))
# Rank order "E. Sambo","Marker 32","W. Sambo","Big Pine","Dave's Ledge","Looe Key","EDR","Maryland Shoals","Bahia Honda"
# Geographic order "Bahia Honda","Big Pine", "Looe Key","Dave's Ledge","Maryland Shoals","E. Sambo","W. Sambo","Marker 32", "EDR"
GxSMatrix<-spread(df.short,site,PercSurv)
GcorrBySite<-cor(GxSMatrix[,c(2:10)])
quartz()
corrplot(GcorrBySite,type="upper",tl.col="black",tl.srt=45,method='ellipse')
corrplot.mixed(GcorrBySite,lower = "number",upper="ellipse",tl.col="black")

plotlist<-c(p1,p2)
plot_grid(plotlist = plotlist,ncol=2,labels = c("A","B"))

####
#non mixed-effects version, same story, interaction is significant
mInt<-coxph(died~Site*Genotype, data=score)
#note, model does not converge, coefficients should be ignored but log likelihood score remains valid see ?coxph
mInt$loglik
#first value is for null model, second is for full model
-2*-205.4487 #calculate the -2LogL = 410.8974

mAdd<-coxph(died~Site+Genotype, data=score)
mAdd$loglik
-2*-237.4866 #calculate the -2LogL = 474.9732

#calculate log likelihood ratio statistic = (-2LogL mAdd) - (-2LogL mInt)
474.9732-410.8974 #64.0758-> large improvement in log likelihood, should retain interaction term in model but cannot trust coefficients...
###

########End of survival analysis#####

##########Breakage Analysis
#########
#1. Is fragmentation more likely at certain sites or in genotypes
install.packages("RVAideMemoire")
library(ordinal)
str(grow$CulumativeBreaks)

#Visualizing some trends
p4<-grow %>% group_by(Genotype) %>% 
  summarize(TotalFrag=sum(CulumativeBreaks,na.rm=T)) %>% 
  ggplot(aes(x=Genotype,y=TotalFrag))+geom_bar(stat = "identity")+
  ylab("# of Breaks")

p5<-grow %>% group_by(Site) %>% 
  summarize(TotalFrag=sum(CulumativeBreaks,na.rm=T)) %>% 
  ggplot(aes(x=Site,y=TotalFrag))+geom_bar(stat = "identity")+
  ylab("# of Breaks")

quartz()
ggarrange(p4,p5,nrow=2,labels = c("A","B"))

#In order to compare size and broken (0 vs 1):
test <-grow %>% pivot_longer(cols=c(8,15,22,29),names_to = c('time','trait'),names_sep='_',values_to="size",values_drop_na=F) #sizes by time point
test2<-grow %>% pivot_longer(cols=c(42:45),names_to="BreakTime",values_to="broken",values_drop_na=F) #getting breakage
test3<-grow %>% pivot_longer(cols=c(14,21,28,35),names_to="trait",values_to="Status",values_drop_na=F)
test<-cbind(test,test2[55],test3[55]) #may have to change this depending on the position of last columns 
test<-test[complete.cases(test$Status),]
test<-filter(test,Status!="M")

#GLM for effect of size on potential to break
output <- glm(broken~size, family="binomial", data=test)
summary(output) #size does not impact likelihood of breakage but this may not be the best way of looking at it
anova(output)

#Fisher's exact test for breakage, first set the size classes
test %>% filter(size<5 & broken==0) %>% summarize(Broken=n())  #173
test %>% filter(size<5 & broken==1) %>% summarize(Broken=n()) #15
test %>% filter(size>=5 & size<10 & broken==0) %>% summarize(Broken=n()) #197
test %>% filter(size>=5 & size<10 & broken==1) %>% summarize(Broken=n()) #77
test %>% filter(size>=10 & size<15 & broken==0) %>% summarize(Broken=n()) #131
test %>% filter(size>=10 & size<15 & broken==1) %>% summarize(Broken=n()) #32
test %>% filter(size>=15 & size<20 & broken==0) %>% summarize(Broken=n()) #69
test %>% filter(size>=15 & size<20 & broken==1) %>% summarize(Broken=n()) #20
test %>% filter(size>=20 & size<30 & broken==0) %>% summarize(Broken=n()) #89
test %>% filter(size>=20 & size<30 & broken==1) %>% summarize(Broken=n()) #22
test %>% filter(size>=30 & broken==0) %>% summarize(Broken=n()) #61
test %>% filter(size>=30 & broken==1) %>% summarize(Broken=n()) #11

#making a dataframe that can test the size classes (probably a cleaner way of doing this)
df<-data.frame("Broken"=c(15,77,32,20,22,11), "Unbroken"=c(173,197,131,69,89,61),row.names = c(">5","5-10","10-15","15-20","20-30","30+"))
mosaicplot(df, color = TRUE, main="Size Classes")
library(stats)

fisher.test(df,simulate.p.value = T, B=100)
library(rstatix)
row_wise_fisher_test(df, p.adjust.method = "holm",detailed = F)
# A tibble: 6 × 5
#group     n          p      p.adj p.adj.signif
#* <chr> <dbl>      <dbl>      <dbl> <chr>       
# 1 >5      897 0.00000156 0.00000936 ****        
# 2 5-10    897 0.0000552  0.000276   ***         
# 3 10-15   897 1          1          ns          
# 4 15-20   897 0.485      1          ns          
# 5 20-30   897 1          1          ns          
# 6 30+     897 0.358      1          ns       
#beeswarm plot to show distribution of size ~broken vs. not broken
library(ggbeeswarm)
ggplot(test)+geom_quasirandom(aes(x=as.factor(test$broken),y=size), size = 0.5,na.rm=T,col=as.numeric(test$Status))

#running the Cumulative Link mixed model
broke<-grow
broke$CulumativeBreaks<-factor(broke$CulumativeBreaks,ordered = T,levels = c("0","1","2"))
broke$Genotype<-relevel(broke$Genotype,ref="1")
broke$Site<-relevel(broke$Site,ref="Marker 32")
frag=clm(CulumativeBreaks~Site+Genotype,data=broke,Hess=TRUE,link="probit") 
#warning message about singularity and non convergence goes away when interaction is removed
frag
summary(frag)
logLik(frag)


#1. How often does fragmentation immediately precede death?
sum(grow$CulumativeBreaks,na.rm=T) #177 individual breakage events
grow %>% filter(T3_Break==1 & T3_Status=="D") %>% dim() #19 instances
grow %>% filter(T6_Break==1 & T6_Status=="D") %>% dim() #6
grow %>% filter(T9_Break==1 & T9_Status=="D") %>% dim() #1
grow %>% filter(T12_Break==1 & T12_Status=="D") %>% dim() #2

sum(grow$T12_Dead==1,na.rm=T) #48 dead at T12
sum(grow$T12_Dead==0,na.rm=T) #198 alive at T12
sum(is.na(grow$T12_Dead)) #24 missing

sum(grow$T12_TOD=="13",na.rm=T) #198
sum(grow$T12_TOD=="12",na.rm=T) #3
sum(grow$T12_TOD=="9",na.rm=T) #10
sum(grow$T12_TOD=="6",na.rm=T) #12
sum(grow$T12_TOD=="3",na.rm=T) #23

sum(grow$T12_Status=="M",na.rm=T) #24 missing colonies (puck and frag missing)

28/177 # 15.8% of breakage events were associated with the death of the fragment
#which means 84.2% of breaks were non-fatal

28/48 #58.3% of the mortality events had breakage occuring during the same time point
#42.7% of the mortality events didn't have any breakage associated with it

sum(grow$T3_Break,na.rm = T) #86
sum(grow$T6_Break,na.rm = T) #53
sum(grow$T9_Break,na.rm = T) #24
sum(grow$T12_Break,na.rm = T) #14

##############
#End of breakage analysis

#Just calculating the overall averages for sites and genotype
grow %>% group_by(Genotype) %>%
  summarise(mean(T12_TLE,na.rm=T),
            mean(T12_SA,na.rm=T),
            mean(T12_V,na.rm=T),
            mean(T12_Vinter,na.rm=T))

grow %>% group_by(Site) %>%
  summarise(mean(T12_TLE,na.rm=T),
            mean(T12_SA,na.rm=T),
            mean(T12_V,na.rm=T),
            mean(T12_Vinter,na.rm=T))

####### Mixed models for absolute size 
####################

#1. Effects on absolute size in order to have largest data set
#   - basically considering all of things that happen to a coral, what impacts final size
#   - going to include all breakage no matter how severe

#removing dead colonies (once you die, all future measures are not included)
grow2<-grow
grow2[which(grow2$T3_Status!="A"),c(15:20,22:27,29:34,36:41,54:57)]<-NA #if a coral is dead, replace all following measures with NA
grow2[which(grow2$T6_Status!="A"),c(22:27,29:34,36:41,55:57)]<-NA
grow2[which(grow2$T9_Status!="A"),c(29:34,36:41,56:57)]<-NA
grow2[which(grow2$T12_Status!="A"),c(36:41,57)]<-NA

#pivoting data set to have a column for trait value and time and trait
size<-grow2 %>% pivot_longer(cols=c(15,17,18,22,24,25,29,31,32,36,38,39,54:57),names_to = c('time','trait'),names_sep='_',values_to="size",values_drop_na=T)
size$time<-factor(size$time,ordered = T,levels=c("T3","T6","T9","T12")) #making time ordinal

#mixed model for TLE size
  TLEmod<-lmer(log(size+1)~Genotype*Site*time+T0_TLE+(1|Site:Array)+(1|CulumativeBreaks),data=size[which(size$trait=="TLE"),],REML = T)
  qqPlot(residuals(TLEmod),xlab="Theoretical Quantiles",ylab="Observed Quantiles",main = "TLEmod") #try a transform if data are non-normal, log(x+1) looks really good
  
  plot(TLEmod,which=1,main="TLEmod") #want random scatter; no apparent trend line, looks better after log(x+1)
  
  #Looking at results of MEM
  anova(TLEmod,ddf="Kenward-Roger") # genotype, time, initial size, and GxS are sig
  summary(TLEmod) #
  rand(TLEmod) #both random effects are very sig

#SA
  SAmod<-lmer(sqrt(size)~Genotype*Site*time+T0_SA+(1|Site:Array)+(1|CulumativeBreaks),data=size[which(size$trait=="SA"),],REML = T)
  qqPlot(residuals(SAmod),xlab="Theoretical Quantiles",ylab="Observed Quantiles",main = "SAmod") #try a transform if data are non-normal, sqrt() looks best, log(x+1) is good too
  
  plot(SAmod,which=1,main="SAmod") #want random scatter; no apparent trendline, looks better after sqrt()
  
  min(size$size[which(size$trait=="SA")]) #no real need to log (x+1) because there are no negative values
  
  anova(SAmod,ddf="Kenward-Roger") # genotype, time, initial Size, GxS, SxT are sig
  summary(SAmod) #doesnt look like any G or GxS combinations are sig from reference
  rand(SAmod) #both are significant

#V
  Vmod<-lmer(log(size+1)~Genotype*Site*time+T0_V+(1|Site:Array)+(1|CulumativeBreaks),data=size[which(size$trait=="V"),],REML = T)
  qqPlot(residuals(Vmod),xlab="Theoretical Quantiles",ylab="Observed Quantiles",main = "Vmod") #try a transform if data are non-normal, log(x+1) is good
  
  plot(Vmod,which=1,main="Vmod") #want random scatter; no apparent trendline, looks better after log(x+1)
  min(size$size[which(size$trait=="V")]) # log (x+1) because there are 0 values
  
  anova(Vmod,ddf="Kenward-Roger") # genotype, time,initial size, GxS, SxT are sig
  summary(Vmod) #doesnt look like any G or GxS combinations are sig from reference
  rand(Vmod) #both are significant

#Vinter
  Vimod<-lmer(log(size+1)~Genotype*Site*time+T0_Vinter+(1|Site:Array)+(1|CulumativeBreaks),data=size[which(size$trait=="Vinter"),],REML = T)
  qqPlot(residuals(Vimod),xlab="Theoretical Quantiles",ylab="Observed Quantiles",main = "Vimod") #try a transform if data are non-normal, log(x+1) is good
  
  plot(Vimod,which=1,main="Vimod") #want random scatter; no apparent trendline, looks better after log(x+1)
  min(size$size[which(size$trait=="Vinter")]) # log (x+1) because there are 0 values
  
  anova(Vimod,ddf="Kenward-Roger") # genotype, time,initial size, GxS, SxT are sig
  summary(Vimod) #doesnt look like any G or GxS combinations are sig from reference
  rand(Vimod) #both are significant

#figure to show how initial size predicts final size
ggplot(size[which(size$trait=="TLE"),],aes(x=T0_TLE,y=size))+
  geom_point()+
  labs(y='TLE  (cm)')+
  geom_smooth(method='lm')+
  stat_cor(method = "pearson", label.x = NULL, label.y = NULL)+
  facet_wrap(~time, scales = "free")

#####
#End of mixed modeling for absolute size

####Mixed models for growth rate
#####
#   - perhaps makes more sennse to be used for plasticity
#   - removing any coral with a break, only want positive growth
#   - using growth rate per 3 month period as a trait

##Calculating the growth metrics, this will be growth over each 3 month time point recorded as growth rate per month
#could consider using something like 'standardized productivity' where new growth is standardized by old growth 
#   final-initial/time/initial size  <- does this have issues with the repeat of initial size?
#   drury et al. used inital number of branches to standardize by inital size
#TLE
grow<-
  grow %>% 
  mutate(TLE_0.3=(T3_TLE-T0_TLE)/3,
         TLE_3.6=(T6_TLE-T3_TLE)/3,
         TLE_6.9=(T9_TLE-T6_TLE)/3,
         TLE_9.12=(T12_TLE-T9_TLE)/3)
#SA
grow<-
  grow %>% 
  mutate(SA_0.3=(T3_SA-T0_SA)/3,
         SA_3.6=(T6_SA-T3_SA)/3,
         SA_6.9=(T9_SA-T6_SA)/3,
         SA_9.12=(T12_SA-T9_SA)/3)
#V
grow<-
  grow %>% 
  mutate(V_0.3=(T3_V-T0_V)/3,
         V_3.6=(T6_V-T3_V)/3,
         V_6.9=(T9_V-T6_V)/3,
         V_9.12=(T12_V-T9_V)/3)
#Vinter
grow<-
  grow %>% 
  mutate(Vinter_0.3=(T3_Vinter-T0_Vinter)/3,
         Vinter_3.6=(T6_Vinter-T3_Vinter)/3,
         Vinter_6.9=(T9_Vinter-T6_Vinter)/3,
         Vinter_9.12=(T12_Vinter-T9_Vinter)/3)

#Now we want to try to deal with breakage. 
#We want to remove every colony that had breakage. 
#Removing colonies that lost any of its TLE over any 3 month period.

#what if I completely removed a colony with any severe breakage -> leaves us with 76 corals and extremely poor replication across genotypes and sites.
#test<-grow
#grow2<-grow%>%
#  filter((TLE_0.3)>0) %>% 
#  filter((TLE_3.6)>0) %>% 
#  filter((TLE_6.9)>0) %>% 
#  filter((TLE_9.12)>0)  

grow3<-grow
#want to put NA in all following time points if dead (some instances where we measured a dead coral)
grow3[which(grow3$T3_Status!="A"),c(58:73)]<-NA #if growth is negative or a coral is dead, replace all following GRs with NA
grow3[which(grow3$T6_Status!="A"),c(59:61,63:65,67:69,71:73)]<-NA
grow3[which(grow3$T9_Status!="A"),c(60,61,64,65,68,69,72,73)]<-NA
grow3[which(grow3$T12_Status!="A"),c(61,65,69,73)]<-NA

#Removing the instance of breakage but still allowing for tubthumping (i.e. I get knocked down, but I get up again)
#but what is considered breakage, any negative growth, 
#some negative to allow for errors in models (if so, use numbers in comments that are based on CV for each trait)
grow3[which(grow3$TLE_0.3<(0)),c("TLE_0.3","Vinter_0.3")]<-NA #-.19
grow3[which(grow3$SA_0.3<(0)),c("SA_0.3")]<-NA #-1.5
grow3[which(grow3$V_0.3<(0)),c("V_0.3")]<-NA #-.6

grow3[which(grow3$TLE_3.6<(0)),c("TLE_3.6","Vinter_3.6")]<-NA #-.23
grow3[which(grow3$SA_3.6<(0)),c("SA_3.6")]<-NA #-1.98
grow3[which(grow3$V_3.6<(0)),c("V_3.6")]<-NA #-.89

grow3[which(grow3$TLE_6.9<(0)),c("TLE_6.9","Vinter_6.9")]<-NA #-.33
grow3[which(grow3$SA_6.9<(0)),c("SA_6.9")]<-NA #-2.67
grow3[which(grow3$V_6.9<(0)),c("V_6.9")]<-NA #-1.27

grow3[which(grow3$TLE_9.12<(0)),c("TLE_9.12","Vinter_9.12")]<-NA #-.52
grow3[which(grow3$SA_9.12<(0)),c("SA_9.12")]<-NA #-4.13
grow3[which(grow3$V_9.12<(0)),c("V_9.12")]<-NA #-1.98

#Lirman et al 2014, and many others show impact of size on growth rate,
#therefore including size at begining of a time period when modeling growth rate is necessary 
#pivoting data frame to make it easier to work with
grow4<-grow3 %>% pivot_longer(cols=58:73,names_to = c('trait','time'),names_sep='_',values_to="growth_rate") #removing only negative GR events instead of entire colonies that broken give 4x more phenotype measures
grow4<- grow4 %>% separate(time,c("start","end"),"([.])",convert=T,remove=F) #separate time in to start and end point

#now get the size of each fragment at the start of each time interval with this
#there is probably a faster way to do it but this works too
test<-grow3 %>% pivot_longer(col=c(8,10,11,15,17,18,22,24,25,29,31,32,53:56), names_to=c("iTime","trait"),names_sep="_", values_to="size")
test$iTime<-gsub("T","\\1",test$iTime)
#need to arrange the rows of test so that it follows the order of grow4 
test<-test %>% arrange(factor(Site,levels = c("Bahia Honda","Big Pine","Dave's Ledge","E. Sambo","EDR","Looe Key","Marker 32","Maryland Shoals","W. Sambo")),
                       Tag,
                       factor(trait,levels=c("TLE","SA","V","Vinter")),
                       iTime) 

#now add in the initial size for each frag for each trait at the start of each time interval
grow4$size<-ifelse(paste(test$Site,test$Tag,test$trait,test$iTime)==paste(grow4$Site,grow4$Tag,grow4$trait,grow4$start),test$size,NA) #ensuring correct replacement by making sure rows match before pasting value 
grow4<-drop_na(grow4, growth_rate)#there will be missing GR values because of dead/missing coral, however if GR is present, there will always be an associated size
sum(is.na(grow4$size)) #should be 0,  results in 2685 values. 

str(grow4$time)
#making time an ordered factor because these intervals are subsequent
grow4$time <- factor(grow4$time, order = TRUE, 
                                   levels = c("0.3", "3.6", "6.9","9.12"))

######Running the mixed effects models and apply transformations if needed
###
min(grow4$growth_rate[which(grow4$trait=="TLE")]) #lowest value is 0 so log(x+1) is necessary
TLEmod<-lmer(log(growth_rate+1)~Genotype*Site*time+size+(1|CulumativeBreaks)+(1|Site:Array),data=grow4[which(grow4$trait=="TLE"),],REML = T)
qqPlot(residuals(TLEmod),xlab="Theoretical Quantiles",ylab="Observed Quantiles",main = "TLEmod") #looks good with log(x+1) transformation
plot(TLEmod,which=1,main="TLEmod") #want random scatter; no apparent trendline, a little bit of a horizontal volcano

#Pulling out results
anova(TLEmod,ddf="Kenward-Roger") 
#including some negative growth in TLE, changes anova results dramatically
summary(TLEmod) #summary of all permutations of fixed effects
rand(TLEmod) #array is significant, break is not, test to see if should included in model

#fitting model without breaks
TLEmod2<-lmer(log(growth_rate+1)~Genotype*Site*time+size+(1|Site:Array),data=grow4[which(grow4$trait=="TLE"),],REML = T)
AIC(TLEmod,k=2) #1398.163 
AIC(TLEmod2,k=2)#1396.163, slightly better but not huge improvement, but removing breaks makes GxS sig

#SA
min(grow4$growth_rate[which(grow4$trait=="SA")]) #lowest value is .01133
SAmod<-lmer(sqrt(growth_rate)~Genotype*Site*time+size+(1|CulumativeBreaks)+(1|Site:Array),data=grow4[which(grow4$trait=="SA"),],REML = T)
qqPlot(residuals(SAmod),xlab="Theoretical Quantiles",ylab="Observed Quantiles",main="SAmod") #try a transform if data are non-normal
plot(SAmod,which=1) #want random scatter; no apparent trendline

anova(SAmod,ddf="Kenward-Roger") 
summary(SAmod)
rand(SAmod) #array sig but not breakage, should it be removed?

SAmod2<-lmer(sqrt(growth_rate)~Genotype*Site*time+size+(1|Site:Array),data=grow4[which(grow4$trait=="SA"),],REML = T)
AIC(SAmod)  #1954.756
AIC(SAmod2) #1952.756 again, not a huge improvement but removing breakage does not change pattern of results

#V
min(grow4$growth_rate[which(grow4$trait=="V")]) #lowest  value is 0 so log(x+1) is needed
Vmod<-lmer(log(growth_rate+1)~Genotype*Site*time+size+(1|CulumativeBreaks)+(1|Site:Array),data=grow4[which(grow4$trait=="V"),],REML = T)
qqPlot(residuals(Vmod),xlab="Theoretical Quantiles",ylab="Observed Quantiles",main = "Vmod")+ #try a transform if data are non-normal
plot(Vmod,which=1) #want random scatter; no apparent trendline

anova(Vmod,ddf="Kenward-Roger") 
summary(Vmod)
rand(Vmod) #breaks is sig but not array

#Vinter
min(grow4$growth_rate[which(grow4$trait=="Vinter")]) #lowest value is -6.8 so will have to do log(x+7) to successful transform
Vimod<-lmer(log(growth_rate+7)~Genotype*Site*time+size+(1|CulumativeBreaks)+(1|Site:Array),data=grow4[which(grow4$trait=="Vinter"),],REML = T)
qqPlot(residuals(Vimod),xlab="Theoretical Quantiles",ylab="Observed Quantiles") #try a transform if data are non-normal
plot(Vimod,which=1) #want random scatter; no apparent trendline

anova(Vimod,ddf="Kenward-Roger")
summary(Vimod)
rand(Vimod) #array sig, breakage is not sig but also not 1.

Vimod2<-lmer(log(growth_rate+7)~Genotype*Site*time+size+(1|Site:Array),data=grow4[which(grow4$trait=="Vinter"),],REML = T)
AIC(Vimod)  #1435.014
AIC(Vimod2) #1433.018 again, not a huge improvement but removing breakage makes GxS almost significant
########
##End of mixed modeling for growth rate


###Some visualizations for morphology
#####
####For visualizations only
size2<-grow2 %>% pivot_longer(cols=c(8,10,11,15,17,18,22,24,25,29,31,32,36,38,39,53:57),names_to = c('time','trait'),names_sep='_',values_to="size",values_drop_na=T)
size2<-size2 %>% mutate(time=factor(time,levels = c("T0","T3","T6","T9","T12")))

#Wanted to do a PCA for morphology but using absolute size isn't meaningful because the dominant signal is size change over time so instead looking at shape in traits not dependent on size
#Evaluating invariant shape of colonies at end of the experiment at each site
#using only nonbroken colonies to look at shapes of coral able to grow without any physical damage
shape<-grow2 %>% filter(T12_Status=="A")%>% filter(CulumativeBreaks==0) 
shape<-shape %>% mutate(SAtoV=T12_SA/T12_V,
                        SAtoTLE=T12_SA/T12_TLE,
                        Shericity=(4*pi*(((3*T12_V)/(4*pi))^(1/3))^2)/T12_SA,
                        Convexity=T12_V/T12_Vcx,
                        Packing=T12_SA/T12_SAcx)

#Morphology principal components 
pcaMor<-prcomp(shape[,c(58:62)],scale=T)
#if you want to color points by categories of good vs bad sites
shape$Quality<-ifelse(shape$Site %in% c("E. Sambo","Marker 32","W. Sambo"),"Top 3",ifelse(shape$Site %in% c("EDR","Maryland Shoals","Bahia Honda"),"Bottom 3","Mid 3"))
shape$Quality<-as.factor(shape$Quality)

#change components of autoplot to color by site, quality, genet, etc
quartz()
autoplot(pcaMor,data=shape,loadings=T,loadings.label=T,colour='Site',size=2,frame=T,frame.type=,
         loadings.label.size=3,loadings.label.repel=T, loading.label.color="darkgrey", loadings.colour='lightgrey',
         label.size=3)+
  #geom_text(vjust=1.5,label=envC[-6,1],size=3)+
  #geom_text_repel(label=envC[-6,1])+
  theme_classic()+
  theme(legend.position = "right")+
  scale_color_manual(values = colorG)+scale_fill_manual(values=colorS)

#plot of averaged genotype growth at a site
quartz()
size2[which(size2$trait=="TLE"),] %>% #change trait as needed
  group_by(Genotype,Site,time) %>%
  summarise(size=mean(size)) %>%
  mutate(time=factor(time,levels = c("T0","T3","T6","T9","T12"))) %>%
  ggplot(aes(x=time,y=size,group=Genotype))+
  geom_point(aes(color=Genotype))+
  geom_line(aes(color=Genotype,group=Genotype))+
  scale_color_manual(values = colorG)+
  ylab("TLE (cm)")+ xlab(label = NULL)+theme_bw()+
  facet_wrap(~Site,scales = "free")

#use this to plot individual fragments at a site, colored by array
quartz()
size2[which(size2$trait=="TLE"& size2$Site=="Maryland Shoals"),] %>%  #change trait and site as needed
  #group_by(Genotype,Site,time) %>%
  #summarise(size=mean(size)) %>%
  mutate(time=factor(time,levels = c("T0","T3","T6","T9","T12"))) %>%
  ggplot(aes(x=time,y=size,color=Array))+
  geom_point(aes(color=Array))+
  geom_line(aes(color=Array,group=Array))+
  ylab("TLE (cm)")+ xlab(label=NULL)+
  ggtitle("Maryland Shoals")+
  facet_wrap(~Genotype,scales = "free")

#Relationship between size and growth rate, which could be a standardized GR by size
ggplot(grow4[which(grow4$trait=="TLE"),],aes(x=size,y=growth_rate/size))+
  geom_point()+
  geom_smooth()+
  stat_cor(method = "pearson", label.x = NULL, label.y = NULL)+
  labs(y='Standardized TLE growth rate (cm/month)')+facet_wrap(~time,scales = "free") #this looks similar to what Lirman et al. 2014 - standardized growth decreases with size

#visualizing GR over time 
quartz()
grow4[which(grow4$trait=="TLE"),] %>% #change trait as needed
  group_by(Genotype,Site,time) %>%
  summarise(gr=mean(growth_rate)) %>%
  mutate(time=factor(time,levels = c("0.3","3.6","6.9","9.12"))) %>%
  ggplot(aes(x=time,y=gr,color=Genotype))+
  geom_point(aes(color=Genotype))+
  geom_line(aes(color=Genotype,group=Genotype))+
  scale_color_manual(values = colorG)+
  ylab("TLE GR (cm/month)")+ #change title as needed
  facet_wrap(~Site,scales = "free")

#plotting each individual coral at a site by array
grow4[which(grow4$trait=="TLE"& grow4$Site=="Maryland Shoals"),] %>% #change trait and site as needed
  #group_by(Genotype,Site,time) %>%
  #summarise(size=mean(size)) %>%
  #mutate(time=factor(time,levels = c())) %>%
  ggplot(aes(x=time,y=growth_rate,color=Array))+
  geom_point(aes(color=Array))+
  geom_line(aes(color=Array,group=Array))+
  ylab("TLE GR (cm/month)")+ xlab(label=NULL)+  #change title as needed
  ggtitle("Maryland Shoals")+ #change title as needed
  facet_wrap(~Genotype,scales = NULL)

pd=position_dodge(.2)
grow4[which(grow4$trait=="TLE" & grow4$time=="9.12"),] %>% #change trait as needed
  group_by(Genotype,Site) %>%
  summarise(gr=mean(growth_rate)) %>%
  ggplot(aes(x=Site,y=gr, color=Genotype))+
  geom_point(aes(color=Genotype), size = 2,position = pd)+
  geom_line(aes(color=Genotype,linetype=Genotype, group=Genotype),size=1,position = pd)+
  #geom_errorbar(aes(ymin=gr-se, ymax=gr+se, color=Genotype), width=.2,position = pd)+
  scale_color_manual(values = colorG)+
  stat_summary(fun=mean,geom = 'crossbar',color="black", width=0.3)+
  theme(axis.title = element_blank())

######

##Joint Regression Analysis
###########

#function to preform the joint regression analysis and provide the coeffienent, r.squared, and pvalue
#it first calculates environmental averages
#then calculates the genotype averages at each site
#then fits linear regression for a given genotype, trait, time; data is the full dataset with all trait values
jra<-function(grouping1,grouping2,genotype,trait,time,data){
  S<-data %>% group_by_at(grouping1) %>% 
    summarize(Savg=mean(size)) %>% ungroup() %>%
    complete(Site,trait,time)
  G<-data %>% group_by_at(grouping2) %>% 
    summarize(Gavg=mean(size))%>% ungroup() %>%
    complete(Genotype,Site,trait,time)
  mod<-lm(filter(G,Genotype==UQ(genotype)& trait==UQ(trait) &time==UQ(time))$Gavg~
       filter(S,trait==UQ(trait) &time==UQ(time))$Savg)
  return(c(mod$coefficients[[2]],summary(mod)$adj.r.squared,summary(mod)$coefficients[,4][[2]]))
}

group1<-c("Site","trait","time")
group2<-c("Genotype","Site","trait","time")

b<-matrix(nrow=1,c("Genotype","Trait","Time","Plasticity","adj.r.sq","p.val"))
for (ti in unique(size$time)){
  for (tr in unique(size$trait)){
    for(g in unique(size$Genotype)) {
      eop<-jra(grouping1=group1,grouping2=group2,genotype = g,trait = tr,time=ti,data=size)
      b<-rbind(b,c(g,tr,ti,as.numeric(eop)))
    }
  }
}

jr.pls<-as.data.frame(b[-1,]) #safe matrix as dataframe
names(jr.pls)<-c("Genotype", "Trait", "Time", "Plasticity", "adj.r.sq", "p.val") #rename columns
jr.pls$Plasticity<- as.numeric(as.character(jr.pls$Plasticity))  
jr.pls$adj.r.sq<- as.numeric(as.character(jr.pls$adj.r.sq))
jr.pls$p.val<- as.numeric(as.character(jr.pls$p.val))

##If you want to spot check the loop
jra(grouping1=group1,grouping2=group2,genotype = "1",trait = "SA",time="T3",data=size)

#or to check the entire function
  S<-size %>% group_by_at(group1) %>% 
    summarize(Savg=mean(size)) %>% ungroup() %>%
    complete(Site,trait,time)
  G<-size %>% group_by(Genotype, Site, trait, time) %>% 
    summarize(Gavg=mean(size,na.rm=T)) %>% ungroup() %>%
    complete(Genotype,Site,trait,time)
  
  summary(lm(filter(G,Genotype=="1", trait=="SA" &time=="T3")$Gavg~
           filter(S,trait=="SA" &time=="T3")$Savg))

##IF THE REVIEWERS ASK
#####EDR seems to be an outlier for some genotypes, so doing JR with all sites except EDR may be useful
  group1<-c("Site","trait","time")
  group2<-c("Genotype","Site","trait","time")
  
  b2<-matrix(nrow=1,c("Genotype","Trait","Time","Plasticity","adj.r.sq","p.val"))
  for (ti in unique(grow4[which(grow4$Site!="EDR"),]$time)){
    for (tr in unique(grow4[which(grow4$Site!="EDR"),]$trait)){
      for(g in unique(grow4[which(grow4$Site!="EDR"),]$Genotype)) {
        eop<-jra(grouping1=group1,grouping2=group2,genotype = g,trait = tr,time=ti,data=grow4[which(grow4$Site!="EDR"),])
        b2<-rbind(b2,c(g,tr,ti,as.numeric(eop)))
      }
    }
  }
  
  jr.pls2<-as.data.frame(b2[-1,]) #safe matrix as dataframe
  names(jr.pls2)<-c("Genotype", "Trait", "Time", "Plasticity", "adj.r.sq", "p.val") #rename columns
  jr.pls2$Plasticity<- as.numeric(as.character(jr.pls2$Plasticity))  
  jr.pls2$adj.r.sq<- as.numeric(as.character(jr.pls2$adj.r.sq))
  jr.pls2$p.val<- as.numeric(as.character(jr.pls2$p.val))

  #####Removing EDR, a location of very high growth for 2 genotypes, reduced the strength of the 2-way relationships between plasticity and survival.
  #Conclusions may not be very different, plasticity still tends to increase survival but is not significant
  #However,the plasticity vs growth relationships strengthen.
####END OF JOINT REGRESSION
######

#######Relationships between Plasticity, growth, and survival
########

#adding on the survival data taken from Cox models above
jr.pls<-jr.pls %>% arrange(factor(Genotype,levels = c("36","1","50","3","44","7","31","13","62","41"))) #reorder dfs to more easily add new columns
m2$coefficients[9:17] 
#cox risk scores, higher scores means greater chance of death
#Genotype1 Genotype50  Genotype3 Genotype44  Genotype7 Genotype31 Genotype13  Genotype62 Genotype41
#0.9524501  1.2405797  1.1978105  1.0984207  1.5670542  1.6818513  2.3255781  2.5519734  2.6081417
 
#probably a cleaner way of doing this
jr.pls$RiskScore<-c(rep(1,16), rep(0.952,16), rep(1.241,16), rep(1.198,16), rep(1.0984207,16), rep(1.5670542,16), rep(1.6818513,16) ,rep(2.3255781,16) ,rep(2.5519734,16), rep(2.6081417,16))
#pulling global growth rate for each genotype
globeGR<-grow4 %>% group_by(Genotype,trait,time) %>% summarize(meanGR=mean(growth_rate),sd=sd(growth_rate))

#reorder the two dfs to match
globeGR<- globeGR %>% arrange(factor(Genotype,levels=c("36","1","50","3","44","7","31","13","62","41")),
                              factor(trait,levels=c("TLE","SA","V","Vinter")),
                              factor(time, levels = c("0.3", "3.6", "6.9","9.12")))
jr.pls<- jr.pls %>% arrange(factor(Genotype,levels=c("36","1","50","3","44","7","31","13","62","41")),
                            factor(Trait,levels=c("TLE","SA","V","Vinter")),
                            factor(Time, levels = c("T3", "T6", "T9","T12")))

#double check to make sure the two df have identical rows
paste(jr.pls$Genotype,jr.pls$Trait,jr.pls$Time)==paste(globeGR$Genotype,globeGR$trait,globeGR$time)

jr.pls<-cbind(jr.pls,globeGR[,4:5]) #add growth rate to JR dataframe


#Visualizing in 2D space the relationship between plasticity and risk score / p vs GR / GR vs RS

plots<-foreach (x = unique(jr.pls$Trait)) %do% {
  foreach (p=unique(jr.pls$Time)) %do% {
    ggscatter(data=filter(jr.pls,Trait==x,Time==p),x="Plasticity",y="meanGR",
              add="reg.line", conf.int=TRUE,
              cor.coef=TRUE, cor.method="pearson",
              xlab=paste("Plasticity in",x ,"at",p), ylab="Mean GR")
  }
}

#pulling out just a single time point (for figure 4 in manuscript)
plots<-foreach (p=unique(jr.pls$Trait)) %do% {
    ggscatter(data=filter(jr.pls,Trait==p,Time=="T12"),x="Plasticity",y="meanGR",
              add="reg.line", conf.int=TRUE,
              cor.coef=TRUE, cor.method="pearson", color = "RiskScore",
              xlab=paste("Plasticity in",p), ylab="Mean Growth Rate")+
    gradient_color(c( "#0000CC", "#E1EDF3","#FAE7DC","#67001F"))+theme(legend.position = "none")
}

#have to pull out plots from for loop plot list before putting in ggarrange
flist<-c(plots[[1]],plots[[2]],plots[[3]],plots[[4]])
quartz()
ggarrange(plotlist =plots,common.legend = T) 


#write.csv(jr.pls,"/Volumes/GoogleDrive/My Drive/SCHOOL/LAB/CRRPTransplant/JointRegressionResults.csv")

###Visualize the relationship between Plasticity, risk score, and growth rate in 3D space

#simple plots of points in 3D space
library(rgl)
plot3d(x=filter(jr.pls,Trait=="TLE",Time=="9.12")$Plasticity,
       y=filter(jr.pls,Trait=="TLE",Time=="9.12")$RiskScore,
       z=filter(jr.pls,Trait=="TLE",Time=="9.12")$meanGR,
       type="s",radius=.1,
       xlab="Plasticity",ylab="Risk Score",zlab="Growth Rate")

#If you want a surface to give a fitness landscape you need to fill in the space around the points'
#To do this you need to make a matrix where rows and columns are the x and y (in this case plasticity and growth rate)
#Then extrapolated to fill the matrix based on existing points 
library(plotly)

#adding on a surface or "fitness landscape"
  dat=filter(jr.pls,Trait=="SA",Time=="T12") #pick the trait and/or time point to focus on
  fit<-lm(RiskScore~meanGR+Plasticity,data=dat) #fit a lm to base surface on
  
  graph_reso <- 60 #set the number of data points to have, increasing this will increase resolution of surface
  x <- seq(min(dat$Plasticity), max(dat$Plasticity), length.out = graph_reso) #samples x values from range of points
  y <- seq(min(dat$meanGR), max(dat$meanGR), length.out = graph_reso) #samples y values from range of points
  risk_surface <- expand.grid(Plasticity = x,meanGR =  y,KEEP.OUT.ATTRS = F) #create open data frame waiting for new z values
  
  risk_surface$RiskScore <- predict.lm(fit, newdata = risk_surface) # predicts z values from x and y
  risk_surface <- acast(risk_surface, Plasticity~meanGR, value.var = "RiskScore") #converts x,y,z values into xy-plane (makes x and y the rows/columns of matrix and z values for each pair)

#finally making the plot itself
p <- plot_ly(data = dat) %>%
  #adds scatter plot in 3D space
  add_trace(x = ~Plasticity, y = ~meanGR, z = ~RiskScore,
            type = "scatter3d", mode = "markers",
            text=~Genotype,textfont=list(size=7), #changes text characteristics
            #colors="black",#colors=colorG,
            opacity = 1,size=5) %>% #changes point characteristics
  layout(title="Relationships for V T6 to T9",scene=list(zaxis=list(autorange="reversed",title="Risk Score"),yaxis=list(title="Average Growth Rate"))) %>%  #general plot layout
  #adds plane (fitness landscape) from the predicted z values generated above
  add_surface(z = risk_surface, x=x,y=y,
            type = "surface",opacity=.5, colorscale = list(c(0,1), c("blue", "red")))

quartz() #running this in base R (no Rstudio) will pop out a html link which may be useful to you
p 
##########
###End of multi-trait comparisons
