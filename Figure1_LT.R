load("/mpp/lipid/dat_analysis.rdata")
load("/mpp/XinyueLi2/lipid/druglipid_MPP.RData")
library(reshape2)
library(ggplot2)
library(data.table)
library(plyr)
library(gridExtra)
################################
### Figure 1a
################################
calc <- function(df) {
  re <- round(sum(df$IS_TC>=240)/nrow(df),3) #TC>=240
  re <- c(re,round(sum(df$IS_LDL>=160)/nrow(df),3)) #LDL>=160
  re <- c(re,round(sum(df$IS_HDL<40)/nrow(df),3)) #HDL<40
  re <- c(re,round(sum(df$IS_TG>=200)/nrow(df),3)) #TG>=200
  re <- c(re,round(sum((df$med_wm_forget+df$med_tcm+df$medLipid)>0)/nrow(df),3))
  #re <- c(re,round(sum(df$IS_TG>=200 | df$IS_TC>=240 | df$IS_HDL<40 | df$IS_LDL>=160 |(df$med_wm_forget+df$med_tcm+df$medLipid)>0)/nrow(df),3))
  re <- c(re,round(sum(df$dyslipidemia_measure>0 |(df$med_wm_forget+df$med_tcm+df$medLipid)>0)/nrow(df),3))
  return(re)
}

## detect outlier
#################################
## detect TC outlier
myfunc_tc <- function(x) {
  if(x==99.9 | x==99.99){x=-999}
  if(x>=25.9*38.67 & x<=103.6*38.67){x=x/10}
  else if(x>=259*38.67){x=x/100}
  else {x=x}
}

IS_TC_new_str <- lapply(dat$IS_TC,myfunc_tc)
IS_TC_new <- as.data.frame(do.call(rbind.data.frame, IS_TC_new_str))
names(IS_TC_new) <- c("IS_TC_new")
dat <- cbind(dat, IS_TC_new)
dat$IS_TC_new[dat$IS_TC_new==-999] <- NA

## detect LDL outlier
myfunc_ldl <- function(x) {
  if(x==99.9 | x==99.99){x=-999}
  if(x>=200*38.67){x=x/100}
  else {x=x}
}

IS_LDL_new_str <- lapply(dat$IS_LDL,myfunc_ldl)
IS_LDL_new <- as.data.frame(do.call(rbind.data.frame, IS_LDL_new_str))
names(IS_LDL_new) <- c("IS_LDL_new")
dat <- cbind(dat, IS_LDL_new)
dat$IS_LDL_new[dat$IS_LDL_new==-999] <- NA

## detect HDL outlier
myfunc_hdl <- function(x) {
  if(x==99.9 | x==99.99 | x==9.99){x=-999}
  if(x>=3.9*38.67 & x<=25.9*38.67){x=x/10}
  else if(x>=39*38.67){x=x/100}
  else {x=x}
}

IS_HDL_new_str <- lapply(dat$IS_HDL,myfunc_hdl)
IS_HDL_new <- as.data.frame(do.call(rbind.data.frame, IS_HDL_new_str))
names(IS_HDL_new) <- c("IS_HDL_new")
dat <- cbind(dat, IS_HDL_new)
dat$IS_HDL_new[dat$IS_HDL_new==-999] <- NA



# detect tg outlier
myfunc_tg <- function(x) {
  if(x==99.9 | x==99.99 | x==9.99){x=-999}
  if(x>=5.7*88.57 & x<=56.5*88.57){x=x/10}
  else if(x>=57*88.57){x=x/100}
  else {x=x}
}

IS_TG_new_str <- lapply(dat$IS_TG,myfunc_tg)
IS_TG_new <- as.data.frame(do.call(rbind.data.frame, IS_TG_new_str))
names(IS_TG_new) <- c("IS_TG_new")
dat <- cbind(dat, IS_TG_new)
dat$IS_TG_new[dat$IS_TG_new==-999] <- NA

# clean data 
dat$IS_TC <- dat$IS_TC_new
dat$IS_HDL <- dat$IS_HDL_new
dat$IS_LDL <- dat$IS_LDL_new
dat$IS_TG <- dat$IS_TG_new


#Distributions of lipid profile overall
calc(dat)#=0.066 0.034 0.150 0.160 0.024 0.319
sum(dat$IS_HDL)/dim(dat)[1] #Mean concentration of HDL-C 56.31041
sum(dat$IS_TG)/dim(dat)[1]#Mean concentration of TG 136.828
sum(dat$IS_LDL>=160|(dat$med_wm_forget+dat$med_tcm+dat$medLipid)>0)/nrow(dat) #0.04698857

# Among people with dyslipidemia, XX (XX%) were treated with lipid-lowering medications
sum((dat$med_wm_forget+dat$med_tcm+dat$medLipid)>0)/sum(dat$IS_LDL>=160|(dat$med_wm_forget+dat$med_tcm+dat$medLipid)>0)
# XX (XX%) achieved control goals for all four biomarkers 
sum(dat$IS_LDL<160&dat$IS_TG<200&dat$IS_TC<240&dat$IS_HDL>=40&(dat$med_wm_forget+dat$med_tcm+dat$medLipid)>0)/sum((dat$med_wm_forget+dat$med_tcm+dat$medLipid)>0)


##by age group
#g_age <- c(35,45,55,65,75)
g_age <- c(35,40,45,50,55,60,65,70,75)

re_age_m <- re_age_f <- matrix(0,nrow=(length(g_age)-1),ncol=6)
colnames(re_age_m) <- colnames(re_age_f) <- c("TC>=240","LDL>=160","HDL<40","TG>=200",
                                              "Lipid-Lowering Medication","Total Dyslipidemia")
temp1 <- dat[dat$Sex==1,]
temp2 <- dat[dat$Sex==2,]
for (i in 1:(length(g_age)-1)) {
  if (i < length(g_age)-1) {
    temp <- temp1[temp1$Age>=g_age[i] & temp1$Age<g_age[i+1],]
    re_age_m[i,] <- calc(temp)*100
    temp <- temp2[temp2$Age>=g_age[i] & temp2$Age<g_age[i+1],] 
    re_age_f[i,] <- calc(temp)*100 
  } else {
    temp <- temp1[temp1$Age>=g_age[i] & temp1$Age<=g_age[i+1],]
    re_age_m[i,] <- calc(temp)*100
    temp <- temp2[temp2$Age>=g_age[i] & temp2$Age<=g_age[i+1],] 
    re_age_f[i,] <- calc(temp)*100 
  }#Why seperate???
}
re_age_m <- cbind(g_age[1:(length(g_age)-1)],re_age_m)
re_age_f <- cbind(g_age[1:(length(g_age)-1)],re_age_f)
colnames(re_age_m)[1] <- colnames(re_age_f)[1] <- "age"
save(re_age_m,re_age_f,file="fig1.RData")
save(re_age_m,re_age_f,file="/mpp/Lingyi Tan/lipidfig1.RData")
rm(temp1,temp2,temp)

###################################
### Male and Female
###################################
load("fig1.RData")
re_age_m <- melt(as.data.frame(re_age_m),id.vars="age") #Wide to long
re_age_f <- melt(as.data.frame(re_age_f),id.vars="age") #Wide to long
re_age <- rbind(re_age_m,re_age_f)
re_age$gender <- c(rep("Male",nrow(re_age_m)),rep("Female",nrow(re_age_f)))
re_age$gender <- ifelse(re_age$gender=="Male","Men","Women")
re_age$gender <- factor(re_age$gender,levels=c("Men","Women"))
re_age$Criteria <- re_age$variable



p1 <- ggplot(re_age,aes(x=age,y=value,fill=Criteria))+
geom_bar(stat = "identity",position="dodge")+
scale_x_continuous(name = "Age",
breaks=c(35,40,45,50,55,60,65,70),
labels=c("35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-75")) +
#breaks=c(35,45,55,65),
#labels=c("35-44","45-54","55-64","65-75")) +
scale_y_continuous(name = "Percentage (%)",
breaks=seq(0,45,by=5),
labels=seq(0,45,by=5),
limits=c(0,80)) +
#scale_fill_brewer() +
scale_fill_manual(
#Change colors here
values=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c"))+
facet_wrap(~gender)+
coord_cartesian(xlim = c(30, 75)) +
theme_bw()+
theme(panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
strip.background = element_rect(colour = "grey", fill = "white"),
panel.border = element_rect(colour = "grey", fill=NA),
legend.key = element_rect(colour = "white"),
axis.text.x  = element_text(angle=90, vjust=0.5))









####################################
## Figure 1b
## treatment and control
#####################################
#####################################
#####################################old calc
#calc <- function(df) {
#  pat <- sum(df$IS_LDL>=160 |(df$med_wm_forget+df$med_tcm+df$medLipid)>0)
#  re <- c(round(pat/nrow(df),3))#prev
#  re <- c(re,round(sum((df$med_wm_forget+df$med_tcm+df$medLipid)>0)/pat,3))#trt
#  re <- c(re,round(sum(df$IS_LDL<=160 & (df$med_wm_forget+df$med_tcm+df$medLipid)>0 )/pat,3))#control
#  return(re)
#}

#####################################
#####################################
#####################################new calc for control
#extremely high risk: control<70mg/dL
#high risk:control<100mg/dL
#low moderate risk:control<130mg/dL

#save(lipid_cleaned,file='/mpp/Lingyi Tan/lipid/lipid_cleaned.rdata')
#load('/mpp/Lingyi Tan/lipid/lipid_cleaned.rdata')
load("/mpp/XinyueLi2/bp/dat2.RData")
lipid_cleaned<-dat
lipid_cleaned$Patient_Id<-as.numeric(lipid_cleaned$Patient_Id)
new_df<-merge(x = lipid_cleaned, y = dat[,c("Patient_Id","TCM_HBP","WM_HBP")], by= "Patient_Id", all.x = TRUE)
s <- new_df$Systolic_Average>=140
d <- new_df$Diastolic_Average>=90
M1 <- new_df$TCM_HBP #chinese med
M2 <- new_df$WM_HBP #Western med
#simulation: 
# lipid_cleaned$TCM_HBP<-rbinom(n=dim(lipid_cleaned)[1],size=1,prob=0.1)
# lipid_cleaned$WM_HBP<-rbinom(n=dim(lipid_cleaned)[1],size=1,prob=0.2)
# new_df<-lipid_cleaned
# s=new_df$Systolic_Average>=140
# d=new_df$Diastolic_Average>=90
# new_df$HBP <- ifelse(s+d+M1+M2>0,1,0)

###Three risk factors:
# new_df$HDL_low<-ifelse(new_df$IS_HDL<40,1,0)
# new_df$age_risk<-ifelse(new_df$Sex==1&new_df$Age>=45|new_df$Sex==2&new_df$Age>=55,1,0)
# new_df$smoke[df$Smoking_Currently==1]<-1
# new_df$smoke[df$Smoking_Currently==2]<-0
#simulation:
# new_df$HDL<-rnorm(dim(new_df)[1],mean=35,sd=8)
# new_df$HDL_low<-ifelse(new_df$HDL<40,1,0)
# new_df$age_risk<-ifelse(new_df$Sex=="Male"&new_df$Age>=45|new_df$Sex=="Female"&new_df$Age>=55,1,0)
# new_df$smoke<-as.numeric(new_df$Smoke)-1
# new_df$IS_LDL<-rnorm(dim(new_df)[1], mean=90, sd=19) 
# new_df$IS_TC<-rnorm(dim(new_df)[1], mean=120, sd=23) 
# new_df$Systolic_Average<-rnorm(dim(new_df)[1], mean=110, sd=15)
# new_df$Diastolic_Average<-rnorm(dim(new_df)[1], mean=80, sd=13)
# new_df$Hx_MI<-rbinom(n=dim(new_df)[1],size=1,prob=0.13)
# new_df$Hx_Stroke<-rbinom(n=dim(new_df)[1],size=1,prob=0.23)
# new_df$med_wm_forget<-rbinom(n=dim(new_df)[1],size=1,prob=0.03)
# new_df$med_tcm<-rbinom(n=dim(new_df)[1],size=1,prob=0.05)
# new_df$medLipid<-rbinom(n=dim(new_df)[1],size=1,prob=0.06)


#Calculate: prevalence of people with different degree of risk
calc <- function(df) {
  pat <- sum(df$IS_LDL>=160 |(df$med_wm_forget+df$med_tcm+df$medLipid)>0)
  statin<- sum(dat$medStatin>0)#9387 for dat
  #extremely high risk
  df$ASCVD <-ifelse(df$Hx_MI>0|df$Hx_Stroke>0,1,0) 
  #high risk in cell 1
  df$HR1<-ifelse(df$IS_TC>=3.1*38.67 & df$IS_TC<=7.2*38.67 & df$Age >40 | df$IS_LDL<=4.9*38.67 & df$IS_LDL>=1.8*38.67 & df$Age >40 |df$IS_TC>=7.2*38.67 |df$IS_LDL>=4.9*38.67 ,1,0) 
  df$TCLDLrisk[df$IS_TC<4.1*38.67&df$IS_TC>=3.1*38.67 |df$IS_LDL<2.6*38.67&df$IS_LDL>=1.8*38.67]<-1
  df$TCLDLrisk[df$IS_TC<5.2*38.67&df$IS_TC>=4.1*38.67 |df$IS_LDL<3.4*38.67&df$IS_LDL>=2.6*38.67]<-2
  df$TCLDLrisk[df$IS_TC<7.2*38.67&df$IS_TC>=5.2*38.67 |df$IS_LDL<4.9*38.67&df$IS_LDL>=3.4*38.67]<-3 
  df$TCLDLrisk[is.na(df$TCLDLrisk)]<-0  
  #risk number in cell 2
  df$num_risk2<-df$smoke+df$HDL_low+df$age_risk
  #risk in cell 3
  df$risk1<-ifelse(df$Systolic_Average>=160 | df$Diastolic_Average>=100,1,0)
  df$risk2<-ifelse(df$IS_TC-df$IS_LDL>=200,1,0)
  df$non_HDL<-df$IS_TC-df$IS_LDL
  df$risk3<-ifelse(df$IS_HDL<=40,1,0)
  df$risk4<-ifelse(df$BMI>=28,1,0)
  df$num_risk3<- df$risk1+  df$risk2+  df$risk3+  df$risk4+df$smoke
  #extremely high risk:1
  #high risk:2
  #moderate risk:3
  #low risk:4
  df$risk<-0
  df$risk[df$ASCVD==1]<-1
  df$risk[df$risk==0 & df$HR1==1| df$risk==0&df$num_risk2==2 & df$TCLDLrisk>1 & df$HBP==1 |df$risk==0 & df$num_risk2==3 & df$TCLDLrisk>0 & df$HBP==1  ]<-2
  df$risk[df$risk==0 &df$num_risk2==2 & df$TCLDLrisk==3 & df$HBP==0 |df$risk==0 &df$num_risk2==1 & df$TCLDLrisk>1 & df$HBP==1 |df$risk==0 &df$num_risk2==2 & df$TCLDLrisk==1 & df$HBP==1 ]<-3
  df$risk[df$risk==0 &df$num_risk2<2 & df$HBP==0 & df$TCLDLrisk>0 |df$risk==0 &df$num_risk2==2 & df$TCLDLrisk <3 & 0< df$TCLDLrisk & df$HBP==0 | df$risk==0 &df$num_risk2==3 & df$HBP==0 & df$TCLDLrisk==1| df$risk==0 &df$num_risk2==0 & df$HBP==1 & df$TCLDLrisk > 0|df$risk==0 & df$num_risk2==1 & df$HBP==1 & df$TCLDLrisk==1 ]<-4
  df$risk[df$risk==3 & df$Age < 55 &df$num_risk3>=2]<-2
  #re <- c(round(pat/nrow(df),3))#prev
  re1<-data.frame(round(table(df$risk)/nrow(df),3)) #prev for diff risk
  colnames(re1)<-c("risk_level","prev")
  #re <- c(re,round(sum((df$med_wm_forget+df$med_tcm+df$medLipid)>0)/pat,3))#trt
  df$med<-ifelse(df$med_wm_forget+df$med_tcm+df$medLipid>0,1,0)
  trt<-NULL
  for (i in levels(re1[,1])){
    temp<-df[df$risk==i,]
    s<-sum((temp$med_wm_forget+temp$med_tcm+temp$medLipid)>0)/dim(temp)[1] #prob of trt/risk
    trt<-c(trt,round(s,3))
  }
  re1<-cbind(re1,trt)
  cont<-NULL
  for (i in levels(re1[,1])){
    temp<-df[df$risk==i,]
    if (i =="0"){s<-sum(temp$med>0)}
    else if (i=="1"){s<-sum(temp$IS_LDL<=70 & temp$med>0)}
    else if (i=="2"){s<-sum(temp$IS_LDL<=100 & temp$med >0)}
    else if (i=="3"){s<-sum(temp$IS_LDL<=130 & temp$med >0)}
    else if (i=="4"){s<-sum(temp$IS_LDL<=130 & temp$med >0)}
    # med<-sum(temp$med)
    s<-s/dim(temp)[1] #prob of control/risk
    cont<-c(cont,round(s,3))
  }
  re1<-cbind(re1,cont)
  #re <- c(re,round(sum(df$risk==0 & df$med >0 | df$risk==1 & df$IS_LDL<=70 & df$med >0 | df$risk==2 & df$IS_LDL<=100 & df$med >0 | df$risk==3 & df$IS_LDL<=130 & df$med >0| df$risk==4 & df$IS_LDL<=130 & df$med >0)/pat,3))#control
  colnames(re1)<-c("risk","Prevalence", "Treatment","Control")
  return(re1)
}


##by age group
#g_age <- c(35,45,55,65,75)
g_age <- c(35,40,45,50,55,60,65,70,75)

dat<-new_df

re_age_m <- re_age_f <- matrix(NA,ncol=5)
colnames(re_age_m) <- colnames(re_age_f) <- c("risk","Prevalence", "Treatment","Control","g_age[i]")
temp_m <- dat[dat$Sex==1,]
temp_f <- dat[dat$Sex==2,]
# Simulation
# temp_m <- dat[dat$Sex=="Male",]
# temp_f <- dat[dat$Sex=="Female",]

for (i in 1:(length(g_age)-1)) {
  if (i < length(g_age)-1) {
    temp <- temp_m[temp_m$Age>=g_age[i] & temp_m$Age<g_age[i+1],]
    temp1<-calc(temp)#*100
    temp1<-cbind(temp1,g_age[i])
    re_age_m<-rbind(re_age_m,temp1)
    
    ###### risk5=extremely high + high
    # risk5<-temp1[temp1$risk==1,]+temp1[temp1$risk==2,]
    # risk5[1]<-5
    # risk5[5]<-g_age[i]
    # re_age_m<-rbind(re_age_m,risk5)

    temp <- temp_f[temp_f$Age>=g_age[i] & temp_f$Age<g_age[i+1],] 
    temp2<-calc(temp)#*100
    temp2<-cbind(temp2,g_age[i])
    re_age_f<-rbind(re_age_f,temp2)
    
    ###### risk5=extremely high + high
    # risk5<-temp2[temp2$risk==1,]+temp2[temp2$risk==2,]
    # risk5[1]<-5
    # risk5[5]<-g_age[i]
    # re_age_m<-rbind(re_age_m,risk5)

} else {
    temp <- temp_m[temp_m$Age>=g_age[i] & temp_m$Age<=g_age[i+1],]
    temp1<-calc(temp)#*100
    temp1<-cbind(temp1,g_age[i])
    re_age_m<- rbind(re_age_m,temp1)#*100
    temp <- temp_f[temp_f$Age>=g_age[i] & temp_f$Age<=g_age[i+1],] 
    temp2<-calc(temp)#*100
    temp2<-cbind(temp2,g_age[i])
    re_age_f<-rbind(re_age_f,temp2)
  }
}

# re_age_m[re_age_m$risk!=c(1,2),]
# re_age_f[re_age_f$risk!=c(1,2),]

colnames(re_age_m)[5] <- colnames(re_age_f)[5] <- "age"
re_age_m<-re_age_m[-1,]
re_age_f<-re_age_f[-1,]
save(re_age_m,re_age_f,file="fig1b.RData")
# rm(temp1,temp2,temp)

calc(dat) # Overall
###################################
### Male and Female
###################################
load("fig1b.RData")
library(reshape2)
re_age_m <- melt(as.data.frame(re_age_m),id.vars=c("age","risk")) #wide-long
re_age_f <- melt(as.data.frame(re_age_f),id.vars=c("age","risk"))
re_age <- rbind(re_age_m,re_age_f)
re_age$risk<-as.factor(re_age$risk)
levels(re_age$risk)<-c("No risk","Extremely high risk", "High risk", "Moderate risk","Low risk")
re_age$risk<-ordered(re_age$risk,levels=c("No risk","Low risk","Moderate risk","High risk","Extremely high risk"))
re_age$gender <- c(rep("Men",nrow(re_age_m)),rep("Female",nrow(re_age_f)))

re_age$Criteria <- re_age$variable
re_age$value<-as.numeric(re_age$value)*100


re_age_tc <- re_age[re_age$variable!="Prevalence",]
re_age_p <- re_age[re_age$variable=="Prevalence",]


##################################
# Updated Figure1b (1)
# Classified by risk level

ggplot(re_age_p,aes(x=age,y=value,fill=risk))+
  geom_bar( stat = "identity",position="fill")+
  scale_x_continuous(name = "Age",
                     breaks=c(35,40,45,50,55,60,65,70),
                     labels=c("35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-75")) +
  #breaks=c(35,45,55,65),
  #labels=c("35-44","45-54","55-64","65-75")) +
  scale_y_continuous(name = "Percentage",
                      breaks=seq(0,1,by=0.1),
                      labels=seq(0,1,by=0.1)
  #                   limits=c(0,50)
                                    ) +
  #scale_fill_brewer() +
  #scale_fill_manual(name = element_blank(),
  #                  labels=c("TC>=240","LDL>=160","HDL<40","TG>=200"),
  #                  values=c("#99CCFF","#3399FF","#0066CC","#0000CC"))+
  facet_wrap(~gender)+
  coord_cartesian(xlim = c(30, 75)) +
  theme_bw()+
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "grey", fill = "white"),
        panel.border = element_rect(colour = "grey", fill=NA),
        legend.key = element_rect(colour = "white"),
        axis.text.x  = element_text(angle=90, vjust=0.5))+
  guides(fill=guide_legend(title=NULL))

##################################



##################################
# Updated Figure1b (2)
# Classified by risk level
# Treatment and control are distinguished by transparency
load("fig1b.RData")
library(reshape2)
re_age <- rbind(re_age_m,re_age_f)
re_age$risk<-as.factor(re_age$risk)
levels(re_age$risk)<-c("No risk","Extremely high risk", "High risk", "Moderate risk","Low risk")
re_age$risk<-ordered(re_age$risk,levels=c("No risk","Low risk","Moderate risk","High risk","Extremely high risk"))
re_age$gender <- c(rep("Men",nrow(re_age_m)),rep("Female",nrow(re_age_f)))
re_age$Criteria <- re_age$variable
re_age[,c(2:4)]<-re_age[,c(2:4)]*100

ggplot(re_age,aes(x=age,y=value))+
  geom_bar( stat = "identity",position="dodge",aes(y=Treatment,fill=risk,alpha=0.7))+
  geom_bar(stat="identity", position="dodge",aes(y=Control,fill=risk,alpha=0.7))+
#  geom_bar( stat = "identity",position="dodge",aes(group=risk),color="#000000")+
#  scale_fill_manual(values=c("#F8766D","#0fC3C7"))+
  scale_x_continuous(name = "Age",
                     breaks=c(35,40,45,50,55,60,65,70),
                     labels=c("35-39","40-44","45-49","50-54","55-59","60-64","65-69","70-75")) +
  
  scale_y_continuous(name = "Percentage (%)",
                     breaks=seq(0,45,by=5),
                     labels=seq(0,45,by=5),
                     limits=c(0,20)) +
  facet_wrap(~gender)+
  coord_cartesian(xlim = c(30, 75)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour = "grey", fill = "white"),
        panel.border = element_rect(colour = "grey", fill=NA),
        legend.key = element_rect(colour = "white"),
        axis.text.x  = element_text(angle=90, vjust=0.5))+
  guides(fill=guide_legend(title=NULL))


########################################################

####
#### Figure 3 Histogram
####
########################################################
re <- round(sum(df$IS_TC>=240)/nrow(df),3)
re <- c(re,round(sum(df$IS_LDL>=160)/nrow(df),3))
re <- c(re,round(sum(df$IS_HDL<40)/nrow(df),3))
re <- c(re,round(sum(df$IS_TG>=200)/nrow(df),3))
dat$dy_TC <- dat$IS_TC>=240
dat$dy_LDL <- dat$IS_LDL>=160
dat$dy_HDL <- dat$IS_HDL<40
dat$dy_TG <- dat$IS_TG>=200
dat$Age_group <- cut(dat$Age,c(34,seq(40,75,by=5)),cinclude.lowest=T,right=T)
par(mfrow=c(2,1))
countTC1 <- unlist(lapply(split(dat[dat$Sex==1,],dat$Age_group[dat$Sex==1]),function(x) sum(x$dy_TC)/nrow(x)))
countTC2 <- unlist(lapply(split(dat[dat$Sex==2,],dat$Age_group[dat$Sex==2]),function(x) sum(x$dy_TC)/nrow(x)))
dfc1 <- data.frame(TC=c(countTC1,countTC2),
                   Age=c(names(countTC1),names(countTC2)),
                   Gender=factor(c(rep("Male",length(countTC1)),rep("Female",length(countTC2)))))
countLDL1 <- unlist(lapply(split(dat[dat$Sex==1,],dat$Age_group[dat$Sex==1]),function(x) sum(x$dy_LDL)/nrow(x)))
countLDL2 <- unlist(lapply(split(dat[dat$Sex==2,],dat$Age_group[dat$Sex==2]),function(x) sum(x$dy_LDL)/nrow(x)))
dfc2 <- data.frame(LDL=c(countLDL1,countLDL2),
                   Age=c(names(countLDL1),names(countLDL2)),
                   Gender=factor(c(rep("Male",length(countLDL1)),rep("Female",length(countLDL2)))))
countHDL1 <- unlist(lapply(split(dat[dat$Sex==1,],dat$Age_group[dat$Sex==1]),function(x) sum(x$dy_HDL)/nrow(x)))
countHDL2 <- unlist(lapply(split(dat[dat$Sex==2,],dat$Age_group[dat$Sex==2]),function(x) sum(x$dy_HDL)/nrow(x)))
dfc3 <- data.frame(HDL=c(countHDL1,countHDL2),
                   Age=c(names(countHDL1),names(countHDL2)),
                   Gender=factor(c(rep("Male",length(countHDL1)),rep("Female",length(countHDL2)))))
countTG1 <- unlist(lapply(split(dat[dat$Sex==1,],dat$Age_group[dat$Sex==1]),function(x) sum(x$dy_TG)/nrow(x)))
countTG2 <- unlist(lapply(split(dat[dat$Sex==2,],dat$Age_group[dat$Sex==2]),function(x) sum(x$dy_TG)/nrow(x)))
dfc4 <- data.frame(TG=c(countTG1,countTG2),
                   Age=c(names(countTG1),names(countTG2)),
                   Gender=factor(c(rep("Male",length(countTG1)),rep("Female",length(countTG2)))))

f1 <- ggplot(dfc1, aes(x=Age, y=TC, fill=Gender)) + 
  geom_bar(stat="identity", position=position_dodge())+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
f2 <- ggplot(dfc2, aes(x=Age, y=LDL, fill=Gender)) + 
  geom_bar(stat="identity", position=position_dodge())+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
f3 <- ggplot(dfc3, aes(x=Age, y=HDL, fill=Gender)) + 
  geom_bar(stat="identity", position=position_dodge())+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
f4 <- ggplot(dfc4, aes(x=Age, y=TG, fill=Gender)) + 
  geom_bar(stat="identity", position=position_dodge())+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
grid.arrange(f1,f2,f3,f4,ncol=2)

s############################
## Appendix Figure 2
## TC, LDL-C, HDL-C, and TG 
## TC ≥240 mg/dL, LDL-C ≥160 mg/dL, HDL-C <40 mg/dL, triglycerides ≥200 mg/dL
############################
theme <-  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

s1_1 <- ggplot(dat, aes(x = IS_TC)) + 
  geom_histogram(position="identity")  + scale_x_continuous(name="Total Cholesterol (mg/dL)")+ theme
s1_2 <- ggplot(dat, aes(x = IS_TG)) + 
  geom_histogram(position="identity")  + scale_x_continuous(name="Triglycerides (mg/dL)")+ theme
s1_3 <- ggplot(dat, aes(x = IS_LDL)) + 
  geom_histogram(position="identity")  + scale_x_continuous(name="LDL (mg/dL)")+ theme
s1_4 <- ggplot(dat, aes(x = IS_HDL)) + 
  geom_histogram(position="identity")  + scale_x_continuous(name="HDL (mg/dL)")+ theme

grid.arrange(s1_1,s1_2,s1_3,s1_4,ncol=2)


################################
###
### Table 2: based on medication
###
################################
#Use of lipid-lowering medications (statin vs. non-statin) among adults who have LDL-C ≥ 160 mg/dL and who have LDL-C ≥ 190 mg/dL, by age 
load("/mpp/XinyueLi2/lipid/druglipid_western.RData")
load("/mpp/XinyueLi2/lipid/druglipid_MPP.RData")
druglipid$is_statin <- as.numeric(druglipid$class_name=="HMG")
druglipid$is_niacin <- as.numeric(druglipid$class_name=="yansuan")
druglipid$is_fibrate <- as.numeric(druglipid$class_name=="benyangsuan")
druglipid$is_ezetimibe <- druglipid$chem_id %in% c(2860,3384)
medStatin <- medWM[medWM$chem %in% druglipid$chem_id[druglipid$is_statin==1],]
medNonstatin <- medWM[medWM$chem %in% druglipid$chem_id[druglipid$is_statin==0],]
dat$medStatin <- as.numeric(dat$Patient_Id %in% medStatin$PID)
dat$medNonstatin <- as.numeric(dat$Patient_Id %in% medNonstatin$PID)
dat$medLipid <- dat$medStatin+dat$medNonstatin
dat$medNiacin <- as.numeric(dat$Patient_Id %in%medWM$PID[medWM$chem %in% druglipid$chem_id[druglipid$is_niacin==1]])
dat$medFibrate <- as.numeric(dat$Patient_Id %in%medWM$PID[medWM$chem %in% druglipid$chem_id[druglipid$is_fibrate==1]])
dat$medEzetimibe <- as.numeric(dat$Patient_Id %in%medWM$PID[medWM$chem %in% druglipid$chem_id[druglipid$is_ezetimibe==1]])
dat$medLipidAll <- dat$medStatin+dat$medNonstatin+dat$med_tcm
dat$medLipidAll2 <- dat$medStatin+dat$medNonstatin+dat$med_tcm+(dat$med_wm_forget>0)

##LDL population
ldl160 <- dat[dat$IS_LDL>=160 | dat$medLipidAll2>0,]
ldl190 <- dat[dat$IS_LDL>=190 | dat$medLipidAll2>0,]
##function
calc <- function(df) {
  ci <- function(p,n) {
    sd <- sqrt(p*(1-p)/n)
    rrr <- round(c(p-qnorm(0.975)*sd,p+qnorm(0.975)*sd),3)
    rrr[1] <- max(rrr[1],0)
    return(rrr)
  }
  form <- function(r,n) {
    rtemp <- ci(r/n,n)
    paste(r,"(",round(r/n,3)*100,"%,",rtemp[1]*100,"%-",rtemp[2]*100,"%)",sep="")
  }
  
  #statin
  statin <- sum(df$medStatin==1)
  statin1 <- sum(df$medStatin==1 & df$medLipidAll==1)
  statin2 <- sum(df$medStatin==1 & df$medLipidAll>1)
  #non-statin
  nonstatin <- sum(df$medNonstatin==1)
  nonstatin1 <- sum(df$medNonstatin==1 & df$medLipidAll==1)
  nonstatin2 <- sum(df$medNonstatin==1 & df$medLipidAll>1)
  #categories
  Niacin <- sum(df$medNiacin>0)
  Fibrate <- sum(df$medFibrate>0)
  Ezetimibe <- sum(df$medEzetimibe>0)
  #TCM
  tcm <- sum(df$med_tcm>0)
  tcm1 <- sum(df$med_tcm>0 & df$medLipidAll==1)
  tcm2 <- sum(df$med_tcm>0 & df$medLipidAll>1)
  #take but not remember
  forget <- sum(df$med_wm_forget>0)
  all <- nrow(df)
  return(c(form(statin,all),
           form(statin1,all),
           form(statin2,all),
           form(nonstatin,all),
           form(nonstatin1,all),
           form(nonstatin2,all),
           form(Niacin,all),
           form(Fibrate,all),
           form(Ezetimibe,all),
           form(tcm,all),
           form(tcm1,all),
           form(tcm2,all),
           form(forget,all),
           form(all,all)))
}

g_age <- c(35,45,55,65,75)
N_cat <- 14

re_med_m <- re_med_f <- matrix(0,ncol=length(g_age),nrow=N_cat)
temp1 <- ldl160
temp2 <- ldl190
re_med_m[,1] <- calc(temp1)
re_med_f[,1] <- calc(temp2)
for (i in 1:(length(g_age)-1)) {
  if (i < length(g_age)-1) {
    temp <- temp1[temp1$Age>=g_age[i] & temp1$Age<g_age[i+1],]
    re_med_m[,i+1] <- calc(temp)
    temp <- temp2[temp2$Age>=g_age[i] & temp2$Age<g_age[i+1],] 
    re_med_f[,i+1] <- calc(temp)
  } else {
    temp <- temp1[temp1$Age>=g_age[i] & temp1$Age<=g_age[i+1],]
    re_med_m[,i+1] <- calc(temp)
    temp <- temp2[temp2$Age>=g_age[i] & temp2$Age<=g_age[i+1],] 
    re_med_f[,i+1] <- calc(temp)
  }
}
write.table(re_med_m,file="table2_160.txt",sep="\t",quote=F,col.names=F)
write.table(re_med_f,file="table2_190.txt",sep="\t",quote=F,col.names=F)
