library(ggplot2)
library(viridis)
library(ggExtra)

sens_HCW <- read.csv("~/SPH Dropbox/HK_COVID/programs/martina/pooled_samples/data/SensPool_HCWscreen.csv")
sens_HCW <- sens_HCW[sens_HCW$Prevalence<5,]

se = sens_HCW$Sensitivity <- sens_HCW$Sensitivity/100
sens_low = sens_HCW$SensCIlow <- sens_HCW$SensCIlow/100
sens_upp = sens_HCW$SensCIupp <- sens_HCW$SensCIupp/100
npool <- sens_HCW$npool
prev = sens_HCW$Prevalence <- sens_HCW$Prevalence/100 # Set prevalence

datapool<-as.data.frame(sens_HCW)

## Plot pooled sens vs prevalence by pool size

datapoolsubn <- datapool[datapool$npool<=5,]
ggplot(datapoolsubn,aes(x=Prevalence,y=Sensitivity,group=as.factor(npool),color=as.factor(npool)))+
  geom_line()+geom_ribbon(aes(ymin=SensCIlow, ymax=SensCIupp,color=as.factor(npool)),color=NA,alpha=0.1)+
  ylab("Pooled sensitivity")+xlab("Prevalence")+ylim(c(0.7,1))+
    theme_bw()+scale_color_viridis(discrete=TRUE)+
  theme(legend.title = element_text(size = 11),legend.text = element_text(size = 10),
        axis.text=element_text(size=10),
        axis.title=element_text(size=11))+labs(color=expression(n[pool]))


## Plot pooled sens vs pool size by prevalence

datapoolsubp <- datapool[which(datapool$Prevalence<=0.01),]
label.sens <- c("0.01"="1","0.005"="0.5","2e-04"="0.02")

ggsens <- ggplot(datapoolsubp,aes(x=npool,y=Sensitivity,group=as.factor(Prevalence),color=as.factor(Prevalence)))+
  geom_line()+geom_ribbon(aes(ymin=SensCIlow, ymax=SensCIupp,color=as.factor(Prevalence)),color=NA,alpha=0.1)+
  ylab("Pooled sensitivity")+#geom_point(size=1)+
  xlab(expression(n[pool]))+ylim(c(0.60,1))+xlim(c(1,20))+
   theme_gray()+removeGridX()+#scale_color_brewer(palette = "Dark2")+
  theme(legend.title = element_text(size = 10),legend.text = element_text(size = 10),
        axis.text=element_text(size=9),
        axis.title=element_text(size=10))+labs(color="Prevalence (%)") + scale_color_viridis(discrete=TRUE,labels=label.sens)



#### Plot pooled sensitivity for 70%,80%,90%

data_plow <- datapoolsubp[which(datapoolsubp$Prevalence==2e-04),]
data_plow <- data_plow[which(data_plow$npool<=10),]

sensvec <- c(0.7,0.8,0.9)
sensM <- sapply(sensvec, function(x){data_plow$Sensitivity*x})
sensL <- sapply(sensvec, function(x){data_plow$SensCIlow*x})
sensH <- sapply(sensvec, function(x){data_plow$SensCIupp*x})

poolsens <- c(sensM[,1],sensM[,2],sensM[,3])
poolsensL <- c(sensL[,1],sensL[,2],sensL[,3])
poolsensH <- c(sensH[,1],sensH[,2],sensH[,3])

sensnopool <- c(rep("70",10),rep("80",10),rep("90",10))
npool <- c(rep(seq(1,10,1),3))
data_sens <- data.frame(poolsens,poolsensL,poolsensH,sensnopool,npool)


ggpool <- ggplot(data_sens,aes(x=npool,y=poolsens,group=as.factor(sensnopool),color=as.factor(sensnopool)))+
  geom_line()+geom_ribbon(aes(ymin=poolsensL, ymax=poolsensH,color=as.factor(sensnopool)),color=NA,alpha=0.1)+
  ylab("Pooled sensitivity")+scale_x_continuous(breaks=c(0,2,4,6,8,10))+#geom_point(size=1)+
  xlab(expression(n[pool]))+ylim(c(0.4,1))+
  theme_light()+removeGridX()+scale_color_brewer(palette = "Dark2")+
  theme(legend.title = element_text(size = 10),legend.position="bottom",legend.text = element_text(size = 10),
        axis.text=element_text(size=9),
        axis.title=element_text(size=10))+labs(color="Unpooled \nsensitivity (%)") #+ scale_color_viridis(discrete=TRUE,labels=label.sens)


