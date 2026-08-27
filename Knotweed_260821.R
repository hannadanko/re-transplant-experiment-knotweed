#Knotweed re-transplant
#2026-08-21
#Sophie Karrenberg

# 1.working directory####
setwd("~/Documents/PROJECTS/Knotweed//Knotweed re-transplant/Analysis retransplant 2024")

# 2.libraries####
library(lme4)
library(ggplot2)
library(lmerTest)
library(dplyr)
library(effects)
library(MASS)
library(lmerPerm)
library(ggpubr)
library(emmeans)
library(multcomp)
library(multcompView)

# 3.data prep####

BM.all <- read.table("CombinedBiomass_250717_SK.csv", sep=";", dec=",", header=T)

str(BM.all)

BM <- data.frame(
  popID=BM.all$X.Population, 
  Ind=BM.all$Individual, 
  pop.ind= paste(BM.all$X.Population, BM.all$Individual, sep="-"),
  Site_2022 = BM.all$GrownInGarden2022, 
  Site_2023 = BM.all$TransplantedInGarden2023,
  Block_2023 = BM.all$Block, 
  remove =rep(0,360),
  Rhiz_wgt_g = BM.all$InitialRhizomeWeight_g,
  Lat_origin = BM.all$Lat_Origin, 
  DW_total_g = BM.all$TotalBiomass_g,
  DW_below_g = BM.all$BelowDW_g,
  DW_above_g = BM.all$AboveDw_g,
  prop.below = BM.all$BelowDW_g/BM.all$TotalBiomass_g)

  str(BM)
  

BM$site.block <- paste(BM$Site_2023, BM$Block_2023, sep=".")
BM$site.block <- as.factor(BM$site.block)
  
BM$popID <- as.factor(BM$popID)
BM$Ind <- as.factor(BM$Ind)
BM$Site_2022 <- as.factor(BM$Site_2022)
BM$Site_2023 <- as.factor(BM$Site_2023)

BM$label <- paste(BM$popID, BM$Ind, sep="-")
BM$label <- as.factor(BM$label)
str(BM)

## 3.1 data checks####
### 3.1.1 data check:Rhizome weight, per harvest and transplant site####
par(mfrow=c(1,3))
interaction.plot(response=BM$Rhiz_wgt_g[BM$Site_2022=="Torino"],x.factor = BM$Site_2023[BM$Site_2022=="Torino"], 
                 trace.factor = BM$popID[BM$Site_2022=="Torino"], legend=F , main="Torino harvest 2022", ylim=c(0,12), xlab="Rhizome given to", las=1, ylab="Rhizome weight (g)")
interaction.plot(response=BM$Rhiz_wgt_g[BM$Site_2022=="Tübingen"],x.factor = BM$Site_2023[BM$Site_2022=="Tübingen"], 
                 trace.factor = BM$popID[BM$Site_2022=="Tübingen"], legend=F, main="Tübingen harvest 2022" , ylim=c(0,12),xlab="Rhizome given to", las=1, ylab="Rhizome weight (g)")
interaction.plot(response=BM$Rhiz_wgt_g[BM$Site_2022=="Uppsala"],x.factor = BM$Site_2023[BM$Site_2022=="Uppsala"], 
                 trace.factor = BM$popID[BM$Site_2022=="Uppsala"], legend=F, main="Uppsala harvest 2022" , ylim=c(0,12), xlab="Rhizome given to", las=1, ylab="Rhizome weight (g)")
# looks like random variation (well done Tübingen!!!), some very heavy rhizomes in all harvests
# more small values + fewer large values -> log-transform

par(mfrow=c(1,3))
interaction.plot(response=log(BM$Rhiz_wgt_g[BM$Site_2022=="Torino"]),x.factor = BM$Site_2023[BM$Site_2022=="Torino"], 
                 trace.factor = BM$popID[BM$Site_2022=="Torino"], legend=F , main="Torino harvest 2022", ylim=c(0,2.5), 
                 xlab="log(Rhizome given to", las=1, ylab="Rhizome weight (g))")
interaction.plot(response=log(BM$Rhiz_wgt_g[BM$Site_2022=="Tübingen"]),x.factor = BM$Site_2023[BM$Site_2022=="Tübingen"], 
                 trace.factor = BM$popID[BM$Site_2022=="Tübingen"], legend=F, main="Tübingen harvest 2022" , ylim=c(0,2.5),
                 xlab="Rhizome given to", las=1, ylab="log(Rhizome weight (g))")
interaction.plot(response=log(BM$Rhiz_wgt_g[BM$Site_2022=="Uppsala"]),x.factor = BM$Site_2023[BM$Site_2022=="Uppsala"], 
                 trace.factor = BM$popID[BM$Site_2022=="Uppsala"], legend=F, main="Uppsala harvest 2022" , ylim=c(0,2.5), 
                 xlab="Rhizome given to", las=1, ylab="log(Rhizome weight (g))")

# looks better, use log(rhizome weight)

### 3.1.2 data check: check rhizome weight and latitude of origin####
par(mfrow=c(1,1))

ggplot(data = BM, aes(Lat_origin, log(Rhiz_wgt_g))) +
  geom_point(show.legend=FALSE) + 
  facet_wrap(~ Site_2022*Site_2023, labeller = "label_both")

### 3.1.3 data check: check rhizome weight and T1 ####
par(mfrow=c(1,1))
boxplot(log(Rhiz_wgt_g)  ~ Site_2022*Site_2023, data=BM, col=c("red", "purple", "lightblue"), las=1)

### 3.1.4 data check: correlation below and above ground weights####
ggplot(data = BM, aes(DW_below_g, DW_above_g)) +
  geom_point(show.legend=FALSE) + 
  facet_wrap(~ Site_2022*Site_2023 ,labeller = "label_both")



### 3.1.5 data check:block effects####
ggplot(data = BM[BM$remove==0,], aes(DW_below_g, DW_above_g, color=site.block )) +
  geom_point(show.legend=FALSE) + 
  geom_label(data = . %>% group_by(Site_2022*Site_2023) %>% filter(remove == 1), aes(label = label, hjust = -0.3, vjust=0.3),size=3,show.legend=FALSE)+
  facet_wrap(~ Site_2022*Site_2023, labeller = "label_both" )

#some block effects, but not extreme


BM[BM$DW_total_g<20,]


# 4.lmer models ####--------------------------------------------
## 4.1 lmer model prep: standardization of rhizome weight and latitude####
BM$popID <- as.factor(BM$popID)
median(BM$DW_above_g, na.rm=T) #16
median(BM$DW_below_g, na.rm=T) #40
median(BM$DW_total_g, na.rm=T) #56
BM$RW.std <- scale(log(BM$Rhiz_wgt_g))

hist(BM$Lat_origin)
BM$Lat_origin.std <- BM$Lat_origin -mean(BM$Lat_origin)

hist(BM$Lat_origin.std)


## 4.2.1 Above ground DW ####
mod.above <- lmer(log(DW_above_g+20 )~  
                    RW.std + Lat_origin.std + Site_2022 + Site_2023+  RW.std:Site_2023 +Lat_origin.std:Site_2023 +Site_2022:Site_2023 + (1|site.block), data=BM ) 

(R.above <- anova(mod.above)[,1:6])

par(mfrow=c(1,2))
plot(mod.above)
qqnorm(resid(mod.above))
qqline(resid(mod.above))

pairs(emmeans(mod.above , "Site_2023"))
pairs(emmeans(mod.above , "Site_2022"))

cld(emmeans(mod.above , "Site_2023"))
cld(emmeans(mod.above , "Site_2022"))
exp(3.62)-20
exp(3.75)-20
(22.52108-12.78595)/12.78595
(22.52108-17.33757)/17.33757
write.table(R.above, "Results_260116.txt", dec=".", sep=",")


## 4.2.2 submodels above-ground BW ####
#Torino
mod.above.To <- lmer(log(DW_above_g+20)~  
                       RW.std + Lat_origin.std + Site_2022 + (1|site.block), data=BM[BM$Site_2023=="Torino",] ) 
anova.above.To <- anova(mod.above.To)[,1:6]
name <- "above.BM.To" 
write.table(name, "S_table.csv")
write.table(anova.above.To, "S_table.csv", append=T)

#Tü
mod.above.Tu <- lmer(log(DW_above_g+20)~  
                       RW.std + Lat_origin.std + Site_2022 + (1|site.block), data=BM[BM$Site_2023=="Tübingen",] ) 
anova.above.Tu <- anova(mod.above.Tu)[,1:6]

name <- "above.BM.Tu" 
write.table(name, "S_table.csv", append=T)
write.table(anova.above.Tu, "S_table.csv", append=T)

#Up
mod.above.Up <- lmer(log(DW_above_g+20)~  
                       RW.std + Lat_origin.std + Site_2022 + (1|site.block), data=BM[BM$Site_2023=="Uppsala",] ) 
anova.above.Up <- anova(mod.above.Up)[,1:6]

name <- "above.BM.Up" 
write.table(name, "S_table.csv", append=T)
write.table(anova.above.Up, "S_table.csv", append=T)




## 4.3.1 Below ground DW ####
mod.below <- lmer(log(DW_below_g+10)~  
                    RW.std + Lat_origin.std + Site_2022 + Site_2023+  RW.std:Site_2023 + Lat_origin.std:Site_2023 +Site_2022:Site_2023 + (1|site.block), data=BM ) 

(R.below <- anova(mod.below)[,1:6])
cld(emmeans(mod.below , "Site_2023"))

exp(3.7)-10
exp(3.85)-10
exp(4.3)-10
(63.69979-30.4473)/63.69979
(36.99306-30.4473)/36.99306
plot(mod.below)
qqnorm(resid(mod.below))
qqline(resid(mod.below))
write.table(R.below, "Results_260116.txt", dec=".", sep=",", append=T)

#Site_2022 : Site_2023 significant#Site_2022 : Site_2023 significanttotal
## 4.3.2 submodels below-ground BM####
#Torino
mod.below.To <- lmer(log(DW_below_g+10)~  
                    RW.std + Lat_origin.std + Site_2022 + (1|site.block), data=BM[BM$Site_2023=="Torino",] ) 

anova.below.To <- anova(mod.below.To)[,1:6]

name <- "below.BM.To" 
write.table(name, "S_table.csv", append=T)
write.table(anova.below.To, "S_table.csv", append=T)

plot(mod.below.To)
qqnorm(resid(mod.below.To))
qqline(resid(mod.below.To))

cld(emmeans(mod.below.To, "Site_2022"), Letters=LETTERS)
pairs(emmeans(mod.below.To, "Site_2022"), Letters=LETTERS)



#Tübingen
mod.below.Tu <- lmer(log(DW_below_g+10)~  
                       RW.std + Lat_origin.std + Site_2022 + (1|site.block), data=BM[BM$Site_2023=="Tübingen",] ) 

anova.below.Tu <- anova(mod.below.Tu)[,1:6]

name <- "below.BM.Tu" 
write.table(name, "S_table.csv", append=T)
write.table(anova.below.Tu, "S_table.csv", append=T)

plot(mod.below.Tu)
qqnorm(resid(mod.below.Tu))
qqline(resid(mod.below.Tu))

pairs(emmeans(mod.below.To, "Site_2022"))
# Tü: site_2022: ns! 

#Uppsala
mod.below.Up <- lmer(log(DW_below_g+10)~  
                       RW.std + Lat_origin.std + Site_2022 + (1|site.block), data=BM[BM$Site_2023=="Uppsala",] ) 

BM[BM$Site_2023=="Uppsala","site.block"]
table(BM[BM$Site_2023=="Uppsala","site.block"])
table(is.na(BM[BM$Site_2023=="Uppsala","DW_below_g"]), BM[BM$Site_2023=="Uppsala","site.block"])

anova.below.Up <- anova(mod.below.Up)[,1:6]

name <- "below.BM.Up" 
write.table(name, "S_table.csv", append=T)
write.table(anova.below.Up, "S_table.csv", append=T)

plot(mod.below.Up)
qqnorm(resid(mod.below.Up))
qqline(resid(mod.below.Up))

cld(emmeans(mod.below.Up, "Site_2022"), Letters=LETTERS)
pairs(emmeans(mod.below.Up, "Site_2022"), Letters=LETTERS)

## 4.4.1 Total BM ####
median(BM$DW_total_g, na.rm=T) #57

mod.total <- lmer(log(DW_total_g+50)~  
                    RW.std + Lat_origin.std + Site_2022 + Site_2023+  RW.std:Site_2023 + Lat_origin.std:Site_2023 +Site_2022:Site_2023 + (1|site.block), data=BM ) 


anova(mod.total)
(R.total <- anova(mod.total)[,1:6])

plot(mod.total)
qqnorm(resid(mod.total))
qqline(resid(mod.total))

cld(emmeans(mod.total, "Site_2022"), Letters=LETTERS)
pairs(emmeans(mod.total, "Site_2022"), Letters=LETTERS)

cld(emmeans(mod.total, "Site_2023"), Letters=LETTERS)
pairs(emmeans(mod.total, "Site_2023"), Letters=LETTERS)
write.table(R.total, "Results_260116.txt", dec=".", sep=",", append=T)

## 4.4.2 submodels total BM####
#Torino
mod.total.To <- lmer(log(DW_total_g+50)~  
                       RW.std + Lat_origin.std + Site_2022 + (1|site.block), data=BM[BM$Site_2023=="Torino",] ) 


anova.total.To <- anova(mod.total.To)[,1:6]

name <- "total.BM.To" 
write.table(name, "S_table.csv", append=T)
write.table(anova.total.To, "S_table.csv", append=T)

#Tübingen
mod.total.Tu <- lmer(log(DW_total_g+50)~  
                       RW.std + Lat_origin.std + Site_2022 + (1|site.block), data=BM[BM$Site_2023=="Tübingen",] ) 


anova.total.Tu <- anova(mod.total.Tu)[,1:6]

name <- "total.BM.Tu" 
write.table(name, "S_table.csv", append=T)
write.table(anova.total.Tu, "S_table.csv", append=T)

#Uppsala
mod.total.Up <- lmer(log(DW_total_g+50)~  
                       RW.std + Lat_origin.std + Site_2022 + (1|site.block), data=BM[BM$Site_2023=="Uppsala",] ) 

anova.total.Up <- anova(mod.total.Up)[,1:6]

name <- "total.BM.Up" 
write.table(name, "S_table.csv", append=T)
write.table(anova.total.Up, "S_table.csv", append=T)


## 4.5.1 proportion below####

mod.prop <- lmer(asin(sqrt(prop.below))~  
                    RW.std + Lat_origin.std + Site_2022 + Site_2023+  RW.std:Site_2023 + Lat_origin.std:Site_2023 +Site_2022:Site_2023 + (1|site.block), data=BM ) 

(R.prop <- anova(mod.prop)[,1:6])
cld(emmeans(mod.prop, "Site_2023"))

sin(0.875)^2
sin(1.042)^2
(0.7861073-0.589123)/0.7861073
(0.7454854-0.589123)/0.7454854
plot(mod.prop)
qqnorm(resid(mod.prop))
qqline(resid(mod.prop))

# 2022:2023 interaction significant

## 4.5.2 submodels proportion below ####
#Torino
prop.below.To <- lmer(asin(sqrt(prop.below))~  
                       RW.std + Lat_origin.std + Site_2022 + (1|site.block), data=BM[BM$Site_2023=="Torino",] ) 


anova.prop.To <- anova(prop.below.To)[,1:6]


name <- "prop.To" 
write.table(name, "S_table.csv", append=T)
write.table(anova.prop.To, "S_table.csv", append=T)
plot(prop.below.To)
qqnorm(resid(prop.below.To))
qqline(resid(prop.below.To))

# ns! 



#Tübingen
prop.below.Tu <- lmer(asin(sqrt(prop.below))~  
                       RW.std + Lat_origin.std + Site_2022 + (1|site.block), data=BM[BM$Site_2023=="Tübingen",] ) 


anova.prop.Tu <- anova(prop.below.Tu)[,1:6]

name <- "prop.Tu" 
write.table(name, "S_table.csv", append=T)
write.table(anova.prop.Tu, "S_table.csv", append=T)
plot(prop.below.Tu)
qqnorm(resid(prop.below.Tu))
qqline(resid(prop.below.Tu))

# Tü: site_2022: ns! 

#Uppsala
prop.below.Up <- lmer(asin(sqrt(prop.below))~  
                       RW.std + Lat_origin.std + Site_2022 + (1|site.block), data=BM[BM$Site_2023=="Uppsala",] ) 



anova.prop.Up <- anova(prop.below.Up)[,1:6]

name <- "prop.Up" 
write.table(name, "S_table.csv", append=T)
write.table(anova.prop.Up, "S_table.csv", append=T)

plot(prop.below.Up)
qqnorm(resid(prop.below.Up))
qqline(resid(prop.below.Up))

cld(emmeans(prop.below.Up, "Site_2022"), Letters=LETTERS)
pairs(emmeans(prop.below.Up, "Site_2022"), Letters=LETTERS)
cld(emmeans(mod.prop, "Site_2023"), Letters=LETTERS)
pairs(emmeans(prop.below.Up, "Site_2022"), Letters=LETTERS)




#
# 5. Plots for effects from lmers ----------------------------------------- 
## 5.1 Site effects plot (initial plot for Fig. 3) ####

p1 <- ggplot(x1, aes(x=Site_2023, y= exp(fit)-20 , fill=Site_2022)) + 
  geom_bar(position="dodge", width=0.7, alpha=1, stat="identity") + 
  scale_fill_manual(values = c("red", "purple", "blue"))+
  geom_errorbar(aes(ymin=exp(fit-se)-20, ymax=exp(fit+se)-20), 
                position =  position_dodge(width=0.7), width=0.4) + 
  theme_bw(base_size=12)+ labs(y="Dry biomass above ground (g)")


p2 <- ggplot(x2, aes(x=Site_2023, y= exp(fit)-10 , fill=Site_2022)) + 
  geom_bar(position="dodge", width=0.7, alpha=1, stat="identity") + 
  scale_fill_manual(values = c("red", "purple", "blue"))+
  geom_errorbar(aes(ymin=exp(fit-se)-10, ymax=exp(fit+se)-10), 
                position =  position_dodge(width=0.7), width=0.4) + 
  theme_bw(base_size=12)+ labs(y="Dry biomass below ground (g)")



p3 <- ggplot(x3, aes(x=Site_2023, y= exp(fit)-50 , fill=Site_2022)) + 
  geom_bar(position="dodge", width=0.7, alpha=1, stat="identity") + 
  scale_fill_manual(values = c("red", "purple", "blue"))+
  geom_errorbar(aes(ymin=exp(fit-se)-50, ymax=exp(fit+se)-50), 
                position =  position_dodge(width=0.7), width=0.4) + 
  theme_bw(base_size=12)+ labs(y="Total dry biomass ground (g)")


p3 <- ggplot(x3, aes(x=Site_2023, y= exp(fit)-50 , fill=Site_2022)) + 
  geom_bar(position="dodge", width=0.7, alpha=1, stat="identity") + 
  scale_fill_manual(values = c("red", "purple", "blue"))+
  geom_errorbar(aes(ymin=exp(fit-se)-50, ymax=exp(fit+se)-50), 
                position =  position_dodge(width=0.7), width=0.4) + 
  theme_bw(base_size=12)+ labs(y="Total dry biomass ground (g)")

p4 <- ggplot(x4, aes(x=Site_2023, y= sin(fit)^2 , fill=Site_2022)) + 
  geom_bar(position="dodge", width=0.7, alpha=1, stat="identity") + 
  scale_fill_manual(values = c("red", "purple", "blue"))+
  geom_errorbar(aes(ymin=sin(fit-se)^2, ymax=sin(fit+se)^2), 
                position =  position_dodge(width=0.7), width=0.4) + 
  theme_bw(base_size=12)+ labs(y="Prop. below-ground biomass")


ggarrange(p1, p2, p3, p4)


## 5.2 Plot rhizome weight effects (initial plot for Fig. 4) ####

p1 <- plot(Effect(c("RW.std", "Site_2023"), mod.above), 
           multiline = T, ci.style="bands", colors = c("red", "purple", "blue"), main="")
p2 <- plot(Effect(c("RW.std", "Site_2023"), mod.below), multiline = T, ci.style="bands", colors = c("red", "purple", "blue"), main="")
p3 <- plot(Effect(c("RW.std", "Site_2023"), mod.total), multiline = T, ci.style="bands", colors = c("red", "purple", "blue"), main="")
p4 <- plot(Effect(c("RW.std", "Site_2023"), mod.prop), multiline = T, ci.style="bands", colors = c("red", "purple", "blue"), main="")

ggarrange(p1, p2, p3, p4)

