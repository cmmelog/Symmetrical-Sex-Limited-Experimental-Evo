library(ggplot2)
library(plotrix) #for like confidence intervals
library(boot)#for bootstrapping
library(lme4)
library(rstatix) #Anova
library(ggpubr)
library(tidyverse)
#library(readxl)#read excel
library(plyr)#


# #columns #Line: LineID number #Poulation: Population ID (H for High
# temperature Treatment, B for Benign Temperature) #Treatment: The # temperature the selection was performed at, for Ancestors, the
# temperature of testing #Chromsome: whether the Female-Limited (FL) or # Male-Limited (ML) treatments #Sex: Male or Female meaning Male or Female
# fitness assay #Lay.Focal.Num/Lay.Competitor.Num: the number of flies of # the testing sex (See previous column) alive when flies are cleared. So
# the number of females laying eggs for female fitness, or males who made # it alive during the lay period. Also includes whether flies escaped. So
# 6.0 (or 6) indicates 6 alive females in laying vial (until the end) and # 0 escapees. 6.1 means 6 females alive when flipping females into laying
# vial, but one female escaped during flipping process (so did not make it #                                                       into vial due to human process error). Escappees were corrected for
# female fitness assay in csv file. All male fitness corrections made # here. No correction was made dead flies. There was no distinction of
# when flies died. #Lay.Ancestor.Num: the number of T4/Ancestor flies of # the non-testing sex present in laying vial. For male fitness assay this
# would be females. This was only counted for male fitness assay. Female # Fitness assays started with 12 males and only deaths that occured during
# interaction phase were counted. #Block: Block number. Female and Male # assays were each done in 2 blocks each.
# #Offspring.Focal/Offsring.Competitor: the total offspring count of both
# sexes. For female assay the entire vial with 12 laying females was # counted. For male fitness assay
# #Female.Offspring.Focal/Female.Offspring.Focal: The number of offspring # that were female from focal or competitor flies

####------ Functions I wrote section-----####
#for substracting mean Realtive Fitness of a group (chr/sex/treatment combo) from Relative.Ft2
groupNorm_fit <-function(data, chr, sex, means){
  Group.Norm.Fitness=subset(data,Focal.Sex==sex&Chromosome==chr)
  #the normalized fitness column
  Group.Norm.Fitness$meanGroupNorm.Fit=NA
  for (pop in means$Population){
    Group.Norm.Fitness[Group.Norm.Fitness$Population==pop,"meanGroupNorm.Fit"]=Group.Norm.Fitness[Group.Norm.Fitness$Population==pop,"Relative.Fit2"]-
      means[means$Focal.Sex==sex&means$Population==pop&means$Chromosome==chr,"Relative.Fit2"]
  }
  return(Group.Norm.Fitness)  
}

#I should rename since it's specific
boot_corr<-function(data,chr,pop,col1,col2){
  BFL.Adboot=replicate(1000,{
    #pull out replicates per population. The default size is n. So it might be uncessary to calculate length?
    samp.B1=sample(data[data$Chromosome==chr&data$Population==pop[1],"Line.ID"],size=length(data[data$Chromosome==chr&data$Population==pop[1],"Line.ID"]),replace=TRUE)#32
    samp.B2=sample(data[data$Chromosome==chr&data$Population==pop[2],"Line.ID"],size=length(data[data$Chromosome==chr&data$Population==pop[2],"Line.ID"]),replace=TRUE)#30
    samp.B3=sample(data[data$Chromosome==chr&data$Population==pop[3],"Line.ID"],size=length(data[data$Chromosome==chr&data$Population==pop[3],"Line.ID"]),replace=TRUE)#29
    
    samp.BFL=append(samp.B1,samp.B2)
    samp.BFL=append(samp.BFL,samp.B3)
    
    #rho is value four
    BFL.rmf=cor.test(data[data$Line.ID%in%samp.BFL,col1],data[data$Line.ID%in%samp.BFL,col2])
    #rho is value four
    #BFL.rmf[4]#this will keep that it's rho
    BFL.rmf[[4]]#this is just value
  })
  
  
  return(BFL.Adboot)
}

##this will groupnormalize any column by population-chr-treatment-Focal.Sex
GroupNorm_func<-function(data,new_col,old_col,sex){
  #the data.frame should contain the new column already
  counter=0
  #this gets the mean per each group
  if(length(sex)>1){
  new.GroupNorm=data_summary(data,old_col,c("Population","Chromosome","Treatment","Focal.Sex"))
  #order: B1-FL-F, B1-FL-M, B1-ML-F, B1-ML-F
  populations=unique(new.GroupNorm$Population)
  for(pop in populations){
    #if sex includes Females
    if("F"%in%sex){
    data[data$Population==pop&data$Chromosome=="FL"&data$Focal.Sex=="F",new_col]=
      data[data$Population==pop&data$Chromosome=="FL"&data$Focal.Sex=="F",old_col]-new.GroupNorm[new.GroupNorm$Focal.Sex=="F",old_col][counter+1]
    data[data$Population==pop&data$Chromosome=="ML"&data$Focal.Sex=="F",new_col]=
      data[data$Population==pop&data$Chromosome=="ML"&data$Focal.Sex=="F",old_col]-new.GroupNorm[new.GroupNorm$Focal.Sex=="F",old_col][counter+2]
    
    #if Males
    if("M"%in%sex){
    data[data$Population==pop&data$Chromosome=="FL"&data$Focal.Sex=="M",new_col]=
      data[data$Population==pop&data$Chromosome=="FL"&data$Focal.Sex=="M",old_col]-new.GroupNorm[new.GroupNorm$Focal.Sex=="M",old_col][counter+1]
    data[data$Population==pop&data$Chromosome=="ML"&data$Focal.Sex=="M",new_col]=
      data[data$Population==pop&data$Chromosome=="ML"&data$Focal.Sex=="M",old_col]-new.GroupNorm[new.GroupNorm$Focal.Sex=="M",old_col][counter+2]
    }
    counter=counter+2
    #print(pop)
  }
    
    
  
  }}
  
  if(length(sex)==1){
    #if there's only one sex won't matter what it is
    new.GroupNorm=data_summary(data,old_col,c("Population","Chromosome","Treatment"))
    #order: B1-FL-F, B1-FL-M, B1-ML-F, B1-ML-F
    populations=unique(new.GroupNorm$Population)
    for(pop in populations){
      
        data[data$Population==pop&data$Chromosome=="FL",new_col]=
          data[data$Population==pop&data$Chromosome=="FL",old_col]-new.GroupNorm[,old_col][counter+1]
        data[data$Population==pop&data$Chromosome=="ML",new_col]=
          data[data$Population==pop&data$Chromosome=="ML",old_col]-new.GroupNorm[,old_col][counter+2]
       
        counter=counter+2
        #print(pop)
      }}
      
  #maybe this is unneccesary
  return(data)
}
#The column name of the variable, the column name(s) of the grouping
summary_func <- function(x, col){
  c(mean = mean(x[[col]], na.rm=TRUE),
    sd = sd(x[[col]], na.rm=TRUE),
    err = sd(x[[col]], na.rm=TRUE)/length(x[[col]]))
}
bootmean.func <- function(x,i){mean(x[i])}
data_summary <- function(data, varname, groupnames_col){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      err = sd(x[[col]], na.rm=TRUE)/(length(x[[col]])-1))
  }
  data_sum<-ddply(data, groupnames_col, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
} 
graph.diff.density<-function(data,IDX,col,legends,leg_pos){
  x <-data[IDX,] #SBGE.chr2$totalExonLength
  x<-x[!is.na(x[[col]]),col]
  y <-data[!IDX,]
  y<-y[!is.na(y[[col]]),col]
  #generate main LINE density
  u<-range(c(x, y))
  dSig.chr2<-density(x,from=u[1], to=u[2]) #works
  dNSig.chr2<-density(y,from=u[1], to=u[2]) #works
  ##dd_xy<-dx$y -dy$y
  diff.SigNon <-dSig.chr2$y - dNSig.chr2$y
  
  #Generate Simple Kernel Density Estimate uing the default R function
  dSig.fit2 <- replicate(1000,{
    #Sample with replacement (for bootstrap from original dataset). Save the resample to x
    samp <- sample(x, replace=TRUE)                    
    #Generate the density from the resampled dataset, and extract y coordinates to generate variablity bands
    density(samp, from=min(dSig.chr2$x), to=max(dSig.chr2$x))$y}) 
  #note!!!.................................mac(X value data))
  #.........................................................^ end DENSITY
  #.........................................................fit stores y ONLY
  #density(x, from=u[1], to=u[2])}) 
  #Apply the quantile function to the y coordinates to get the bounds of the polygon to be drawn on the y axis?
  dSig.fit3 <- apply(dSig.fit2, 1, quantile, c(0.025,0.975) )
  
  dNSig.fit2 <- replicate(1000,{
    samp <- sample(y, replace=TRUE)  
    #NOTE, this is NOT the harmonized range!
    density(samp, from=min(dNSig.chr2$x), to=max(dNSig.chr2$x))$y}) 
  dNSig.fit3 <- apply(dNSig.fit2, 1, quantile, c(0.025,0.975) )
  
  #generate CI for differences
  test=dSig.fit2 - dNSig.fit2
  diff.CI=apply(dSig.fit2-dNSig.fit2, 1, quantile, c(0.025,0.975) )
  diff.CI.1=apply(test, 1, quantile, c(0.025,0.975) )
  lower_graph=min(diff.SigNon)*1.1
  max=range(c(dSig.fit2,dNSig.fit2))
  upper_graph=max[2]
  
  plot(dSig.chr2, col=2,ylab="Density",xlim=c(u[1],u[2]),xlab=col,
       main=paste('Difference in',col, "Distribution" ),ylim=c(lower_graph,upper_graph))#,xlim=c(-10,15))
  #plot(dNSig.chr2 , col=1, ylim=c(0, .20), main='Difference in SBGE Distribution', xlab='',ylab="Density",xlim=c(-10,15))
  polygon( c(dNSig.chr2$x, rev(dNSig.chr2$x)), c(dNSig.fit3[1,], rev(dNSig.fit3[2,])),
           col='lightgrey', density = -0.5, border=F)
  polygon( c(dSig.chr2$x, rev(dSig.chr2$x)), c(dSig.fit3[1,], rev(dSig.fit3[2,])),
           col='lightgreen', density = -0.5, border=F)
  #same x as normal plots but different y
  polygon( c(dSig.chr2$x, rev(dSig.chr2$x)), c(diff.CI[1,], rev(diff.CI[2,])),
           col=2, density = -0.5, border=F)
  abline(h=0, lty=3, col=8) # lty=2, lwd=2,
  lines(dNSig.chr2, col=1,lwd=2)
  lines(dSig.chr2, col=3,lwd=2)
  lines(dSig.chr2$x, diff.SigNon , col="darkred", lty=1, lwd=2)  ## <---------------- difference#
  #mtext(sprintf('N(x) = %s  Bandwidth(x) = %s', dSig.chr2$n, signif(dSig.chr2$bw, 3)), 1, 2)
  #mtext(sprintf('N(y) = %s  Bandwidth(y) = %s', dNSig.chr2$n, signif(dNSig.chr2$bw, 3)), 1, 3)
  legend(leg_pos, legend=legends, col=c(2,1,3), 
         lty=c(1, 1, 2), lwd=c(1, 1, 2))
}



#####-------------------------------------------IMPORT DATA-----
#   [1] "Line"                        "Population"                  "Treatment"                  
# [4] "Chromosome"                  "Focal.Sex"                   "Focal.Female.Death"         
# [7] "Competitor.Female.Death"     "Lay.Ancestor.Death"          "Block"                      
# [10] "Offspring.Focal"             "Offspring.Competitor"        "Offspring.Female.Focal"     
# [13] "Offspring.Female.Competitor" "Not.Total.Fecundity"         "Male.Death"                 
# [16] "Fit1"                        "Sex.Ratio.Focal"             "Sex.Ratio.Comp"             
# [19] "Fit2"                        "Relative.Fit2"               "Offspring.Total"            
# [22] "Line.ID"                     "Rep.ID"                     

SSL.fit=read.table("C:\\Users\\user\\Documents\\Documents\\SSL_FitnessAssay_Temp.csv",header=TRUE,sep=",")
#change types of some columns to factors
SSL.fit$Population=as.factor(SSL.fit$Population)#2
SSL.fit$Treatment=as.factor(SSL.fit$Treatment)#3
SSL.fit$Chromosome=as.factor(SSL.fit$Chromosome)#4
SSL.fit$Focal.Sex=as.factor(SSL.fit$Focal.Sex)#5
SSL.fit$Block=as.factor(SSL.fit$Block)#9



####---MAIN SORTING OF DATA----####
#SSL.fit=subset(SSL.fit.RAW)
#Calculating competitive fitness as ratio of Focal:Competitor Offspring
#SSL.fit$Fit1=SSL.fit$Offspring.Focal/SSL.fit$Offspring.Competitor
#Calculating Sex Ratio
SSL.fit$Sex.Ratio.Focal=SSL.fit$Offspring.Female.Focal/SSL.fit$Offspring.Focal
SSL.fit$Sex.Ratio.Comp=SSL.fit$Offspring.Female.Competitor/SSL.fit$Offspring.Competitor
#Calculating competitive fitness as percent of total flies that were focal
SSL.fit$Fit2=SSL.fit$Offspring.Focal/(SSL.fit$Offspring.Competitor+SSL.fit$Offspring.Focal)

#old way of doing Relative.Fit2 that didn't account for block
# #Derriving Ancestoral average fitness (I should probably do this by block as there are two blocks in each group)
# Ancestor.H.M=mean(SSL.fit[SSL.fit$Chromosome=="T4"&SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="M","Fit2"],na.rm=TRUE)
# Ancestor.H.F=mean(SSL.fit[SSL.fit$Chromosome=="T4"&SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="F","Fit2"],na.rm=TRUE)
# Ancestor.B.M=mean(SSL.fit[SSL.fit$Chromosome=="T4"&SSL.fit$Treatment=="BENIGN"&SSL.fit$Focal.Sex=="M","Fit2"],na.rm=TRUE)
# Ancestor.B.F=mean(SSL.fit[SSL.fit$Chromosome=="T4"&SSL.fit$Treatment=="BENIGN"&SSL.fit$Focal.Sex=="F","Fit2"],na.rm=TRUE)
# #Ancestor.means=c(Ancestor.H.M,Ancestor.B.M,Ancestor.H.F,Ancestor.B.F)
# 
# #Fitness Relative to Average Ancestor
# SSL.fit[SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="M","Relative.Fit2"]=SSL.fit[SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="M","Fit2"]/Ancestor.H.M
# SSL.fit[SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="F","Relative.Fit2"]=SSL.fit[SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="F","Fit2"]/Ancestor.H.F
# SSL.fit[SSL.fit$Treatment=="BENIGN"&SSL.fit$Focal.Sex=="M","Relative.Fit2"]=SSL.fit[SSL.fit$Treatment=="BENIGN"&SSL.fit$Focal.Sex=="M","Fit2"]/Ancestor.B.M
# SSL.fit[SSL.fit$Treatment=="BENIGN"&SSL.fit$Focal.Sex=="F","Relative.Fit2"]=SSL.fit[SSL.fit$Treatment=="BENIGN"&SSL.fit$Focal.Sex=="F","Fit2"]/Ancestor.B.F

#New way of doing Relative.Fit2 accounting for Block

SSL.T4.mean=data_summary(SSL.fit[SSL.fit$Chromosome=="T4",],"Fit2",c("Treatment","Focal.Sex","Block"))$Fit2
# ORDER OF T4 means: [1] "BENIGN T4 F 1" "BENIGN T4 F 2" "BENIGN T4 M 3" [4] "BENIGN T4 M 4" "HIGH T4 F 1"   "HIGH T4 F 2"  
# [7] "HIGH T4 M 3"   "HIGH T4 M 4"   "HIGH T4 M 5"
SSL.fit$Relative.Fit2=0
#B, F1. This will include all CHR, FL, ML and T4
SSL.fit[SSL.fit$Treatment=="BENIGN"&SSL.fit$Focal.Sex=="F"&SSL.fit$Block==1,"Relative.Fit2"]=
  SSL.fit[SSL.fit$Treatment=="BENIGN"&SSL.fit$Focal.Sex=="F"&SSL.fit$Block==1,"Fit2"]/SSL.T4.mean[1]
#B, F2
SSL.fit[SSL.fit$Treatment=="BENIGN"&SSL.fit$Focal.Sex=="F"&SSL.fit$Block==2,"Relative.Fit2"]=  
  SSL.fit[SSL.fit$Treatment=="BENIGN"&SSL.fit$Focal.Sex=="F"&SSL.fit$Block==2,"Fit2"]/SSL.T4.mean[2]
#B, M3
SSL.fit[SSL.fit$Treatment=="BENIGN"&SSL.fit$Focal.Sex=="M"&SSL.fit$Block==3,"Relative.Fit2"]=
  SSL.fit[SSL.fit$Treatment=="BENIGN"&SSL.fit$Focal.Sex=="M"&SSL.fit$Block==3,"Fit2"]/SSL.T4.mean[3]
#B, M4
SSL.fit[SSL.fit$Treatment=="BENIGN"&SSL.fit$Focal.Sex=="M"&SSL.fit$Block==4,"Relative.Fit2"]=   
  SSL.fit[SSL.fit$Treatment=="BENIGN"&SSL.fit$Focal.Sex=="M"&SSL.fit$Block==4,"Fit2"]/SSL.T4.mean[4]
#H, F1
SSL.fit[SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="F"&SSL.fit$Block==1,"Relative.Fit2"]=
  SSL.fit[SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="F"&SSL.fit$Block==1,"Fit2"]/SSL.T4.mean[5]
#H, F2
SSL.fit[SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="F"&SSL.fit$Block==2,"Relative.Fit2"]=
  SSL.fit[SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="F"&SSL.fit$Block==2,"Fit2"]/SSL.T4.mean[6]
#H, M3
SSL.fit[SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="M"&SSL.fit$Block==3,"Relative.Fit2"]=
  SSL.fit[SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="M"&SSL.fit$Block==3,"Fit2"]/SSL.T4.mean[7]
#H, M4
SSL.fit[SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="M"&SSL.fit$Block==4,"Relative.Fit2"]=
  SSL.fit[SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="M"&SSL.fit$Block==4,"Fit2"]/SSL.T4.mean[8]
#H, M5
SSL.fit[SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="M"&SSL.fit$Block==5,"Relative.Fit2"]=
  SSL.fit[SSL.fit$Treatment=="HIGH"&SSL.fit$Focal.Sex=="M"&SSL.fit$Block==5,"Fit2"]/SSL.T4.mean[9]





#Creating a Total Fecundity measure
SSL.fit$Offspring.Total=SSL.fit$Offspring.Focal+SSL.fit$Offspring.Competitor
#creating the various line and rep IDs
SSL.fit$Line.ID=paste0(SSL.fit$Population,SSL.fit$Chromosome,SSL.fit$Line)
SSL.fit[SSL.fit$Chromosome=="T4","Line.ID"]=paste(SSL.fit[SSL.fit$Chromosome=="T4","Line.ID"],SSL.fit[SSL.fit$Chromosome=="T4","Block"])
SSL.fit$Rep.ID=as.factor(paste(SSL.fit$Line.ID, SSL.fit$Focal.Sex))

SSL.fit$Surviving.Females=12
SSL.fit[SSL.fit$Focal.Sex=="M","Surviving.Females"]=12-SSL.fit[SSL.fit$Focal.Sex=="M","Lay.Ancestor.Death"]
SSL.fit[SSL.fit$Focal.Sex=="F","Surviving.Females"]=SSL.fit[SSL.fit$Focal.Sex=="F","Focal.Female.Death"]+ SSL.fit[SSL.fit$Focal.Sex=="F","Competitor.Female.Death"]
SSL.fit$Avg.Offspring.per.F=SSL.fit$Offspring.Total/SSL.fit$Surviving.Females

######creating subgroups as necessary
#No Ancestor
NoA=subset(SSL.fit, Chromosome!="T4"&Rep.ID!="H1FL12.5 F")
#Treatment Temperature
# SSL.B.fit=subset(SSL.fit,Treatment=="BENIGN")
# SSL.H.fit=subset(SSL.fit,Treatment=="HIGH")
# #Separating by sex
# SSL.female=subset(SSL.fit,Focal.Sex=="F"&Rep.ID!="H1FL12.5 F")
# NoA.fem=subset(SSL.female, Chromosome!="T4")
# SSL.male=subset(SSL.fit,Focal.Sex=="M")
# #These are outliers-->exclude both ML and FL outliers
# outliers=SSL.male$Fit2<0.5&SSL.male$Chromosome!="T4"
# SSL.male.clean=SSL.male[!outliers,]
# SSL.male.clean=SSL.male.clean[!is.na(SSL.male.clean$Block),]#only if you have Fit2 data
# NoA.m=subset(SSL.male.clean, Population!="T4"&Chromosome!="T4")

####========================================------->MAIN RELATIVE FITNESS Analysis (with subsample Ancestor)--(use)----####

##the old way that I performed fitness bootstrapping, using fecundity and mate harm as example
####---ALL BAR COMPARISON, SEX BY CHROMOSOME by TREATMENT (bc need for ancestor)----Fit2
###SUBSAMPLE and get the MEANS!!
Sex.Chr.Treat.boot=replicate(1000,{
  #pull out replicates (all)
  samp.ID<-sample(SSL.fit[,"Rep.ID"],replace=TRUE)
  #get the specific control for each group (there are 9)
  samp.T4.mean=data_summary(SSL.fit[SSL.fit$Rep.ID%in%samp.ID&SSL.fit$Chromosome=="T4",],"Fit2",c("Treatment","Focal.Sex","Block"))$Fit2
  # ORDER OF T4 means: [1] "BENIGN T4 F 1" "BENIGN T4 F 2" "BENIGN T4 M 3" [4] "BENIGN T4 M 4" "HIGH T4 F 1"   "HIGH T4 F 2"  
  # [7] "HIGH T4 M 3"   "HIGH T4 M 4"   "HIGH T4 M 5"
  samp.fit=SSL.fit[SSL.fit$Rep.ID%in%samp.ID,]
  samp.fit$Relative.Fitness=0
  #B, F1
  samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="F"&samp.fit$Block==1,"Relative.Fitness"]=
    samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="F"&samp.fit$Block==1,"Fit2"]/samp.T4.mean[1]
  #B, F2
  samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="F"&samp.fit$Block==2,"Relative.Fitness"]=  
    samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="F"&samp.fit$Block==2,"Fit2"]/samp.T4.mean[2]
  #B, M3
  samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="M"&samp.fit$Block==3,"Relative.Fitness"]=
    samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="M"&samp.fit$Block==3,"Fit2"]/samp.T4.mean[3]
  #B, M4
  samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="M"&samp.fit$Block==4,"Relative.Fitness"]=   
    samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="M"&samp.fit$Block==4,"Fit2"]/samp.T4.mean[4]
  #H, F1
  samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="F"&samp.fit$Block==1,"Relative.Fitness"]=
    samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="F"&samp.fit$Block==1,"Fit2"]/samp.T4.mean[5]
  #H, F2
  samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="F"&samp.fit$Block==2,"Relative.Fitness"]=
    samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="F"&samp.fit$Block==2,"Fit2"]/samp.T4.mean[6]
  #H, M3
  samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="M"&samp.fit$Block==3,"Relative.Fitness"]=
    samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="M"&samp.fit$Block==3,"Fit2"]/samp.T4.mean[7]
  #H, M4
  samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="M"&samp.fit$Block==4,"Relative.Fitness"]=
    samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="M"&samp.fit$Block==4,"Fit2"]/samp.T4.mean[8]
  #H, M5
  samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="M"&samp.fit$Block==5,"Relative.Fitness"]=
    samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="M"&samp.fit$Block==5,"Fit2"]/samp.T4.mean[9]
  
  #this will output results? The order doesn't matter here because necessary changes have already been made
  #Use samp.fit since subsample
  #get values for each population and then average the populations
  pops=data_summary(samp.fit,"Relative.Fitness",c("Chromosome","Treatment","Focal.Sex","Population"))
  data_summary(samp.fit,"Relative.Fitness",c("Chromosome","Treatment","Focal.Sex"))$Relative.Fitness
})
#GET ANCESTOR MEAN
#-------there are 27 combinations of CHR, Treatment, Sex and Block, corresponding to each row in the following order
#  [1] "BENIGN FL F 1" "BENIGN FL F 2" "BENIGN FL M 3" [4] "BENIGN FL M 4" "HIGH FL F 1"   "HIGH FL F 2"  
#  [7] "HIGH FL M 3"   "HIGH FL M 4"   "HIGH FL M 5"   [10] "BENIGN ML F 1" "BENIGN ML F 2" "BENIGN ML M 3"
# [13] "BENIGN ML M 4" "HIGH ML F 1"   "HIGH ML F 2"   [16] "HIGH ML M 3"   "HIGH ML M 4"   "HIGH ML M 5"  
# [19] "BENIGN T4 F 1" "BENIGN T4 F 2" "BENIGN T4 M 3" [22] "BENIGN T4 M 4" "HIGH T4 F 1"   "HIGH T4 F 2"  
# [25] "HIGH T4 M 3"   "HIGH T4 M 4"   "HIGH T4 M 5"  
#-----there are 9 Ancestors, and therefore 18 adjustments

#do this to get the mean Relative.Fit2 for each group which the above should be equivalent to. Make sure ORDER is same as above boot
Sex.Chr.Treat.plotting=data_summary(SSL.fit,"Relative.Fit2",c("Chromosome","Treatment","Focal.Sex"))
#calculate all CI including median, just in case
Sex.Chr.Treat.CI=apply(Sex.Chr.Treat.boot, 1, quantile, c(0.025,0.975,0.5))
Sex.Chr.Treat.plotting$lower_CI=Sex.Chr.Treat.CI[1,]
Sex.Chr.Treat.plotting$upper_CI=Sex.Chr.Treat.CI[2,]
Sex.Chr.Treat.plotting$median=Sex.Chr.Treat.CI[3,]
#change the name of Focal.Sex variables for graphing
#library(plyr)
Sex.Chr.Treat.plotting$Focal.Sex<-revalue(Sex.Chr.Treat.plotting$Focal.Sex,c("F"="Females","M"="Males"))
Sex.Chr.Treat.plotting$Treatment<-revalue(Sex.Chr.Treat.plotting$Treatment,c("HIGH"="NOVEL"))

##


#####-----GRAPH MAIN Relative FITNESS FIGURES: 2 x 2 x 2 comparison!---(use)------####

#ggplot(data=Sex.Chr.Treat.plotting[1:8,],aes(fill=Chromosome,y=median,x=Focal.Sex,colour=Chromosome))+
Fig1B=ggplot(data=Sex.Chr.Treat.plotting[1:8,],aes(x=Focal.Sex,y=Relative.Fit2,fill=Chromosome))+
  scale_fill_brewer(palette="Set1")+ylab("Mean Relative Fitness")+ggtitle("Relative Fit2, \nbootstrap CI subsample ancestor")+
  geom_bar(position=position_dodge(),stat = "summary", fun = "mean")+theme_classic() + xlab(" ")+
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.4,size=0.5,position=position_dodge(.9), colour="black")+
    facet_wrap(~Treatment)+
  geom_hline(yintercept=1.0,colour="black")+
  theme(text = element_text(size = 14,color="black"),legend.position="bottom")+
  theme(plot.title = element_text(size = 17,color="black")) + theme(axis.text.x = element_text(size = 12,color="black"))+
  theme(axis.text.y = element_text(size = 12,color="black"))


####----Combine treatment Sex by CHR---(use)----
Sex.Chr.boot=replicate(1000,{
  #pull out replicates (all)
  samp.ID<-sample(SSL.fit[,"Rep.ID"],replace=TRUE)
  #get the specific control for each group (there are 9)
  samp.T4.mean=data_summary(SSL.fit[SSL.fit$Rep.ID%in%samp.ID&SSL.fit$Chromosome=="T4",],"Fit2",c("Treatment","Focal.Sex","Block"))$Fit2
  # ORDER OF T4 means: [1] "BENIGN T4 F 1" "BENIGN T4 F 2" "BENIGN T4 M 3" [4] "BENIGN T4 M 4" "HIGH T4 F 1"   "HIGH T4 F 2"  
  # [7] "HIGH T4 M 3"   "HIGH T4 M 4"   "HIGH T4 M 5"
  samp.fit=SSL.fit[SSL.fit$Rep.ID%in%samp.ID,]
  samp.fit$Relative.Fitness=0
  #B, F1
  samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="F"&samp.fit$Block==1,"Relative.Fitness"]=
    samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="F"&samp.fit$Block==1,"Fit2"]/samp.T4.mean[1]
  #B, F2
  samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="F"&samp.fit$Block==2,"Relative.Fitness"]=  
    samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="F"&samp.fit$Block==2,"Fit2"]/samp.T4.mean[2]
  #B, M3
  samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="M"&samp.fit$Block==3,"Relative.Fitness"]=
    samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="M"&samp.fit$Block==3,"Fit2"]/samp.T4.mean[3]
  #B, M4
  samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="M"&samp.fit$Block==4,"Relative.Fitness"]=   
    samp.fit[samp.fit$Treatment=="BENIGN"&samp.fit$Focal.Sex=="M"&samp.fit$Block==4,"Fit2"]/samp.T4.mean[4]
  #H, F1
  samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="F"&samp.fit$Block==1,"Relative.Fitness"]=
    samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="F"&samp.fit$Block==1,"Fit2"]/samp.T4.mean[5]
  #H, F2
  samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="F"&samp.fit$Block==2,"Relative.Fitness"]=
    samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="F"&samp.fit$Block==2,"Fit2"]/samp.T4.mean[6]
  #H, M3
  samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="M"&samp.fit$Block==3,"Relative.Fitness"]=
    samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="M"&samp.fit$Block==3,"Fit2"]/samp.T4.mean[7]
  #H, M4
  samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="M"&samp.fit$Block==4,"Relative.Fitness"]=
    samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="M"&samp.fit$Block==4,"Fit2"]/samp.T4.mean[8]
  #H, M5
  samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="M"&samp.fit$Block==5,"Relative.Fitness"]=
    samp.fit[samp.fit$Treatment=="HIGH"&samp.fit$Focal.Sex=="M"&samp.fit$Block==5,"Fit2"]/samp.T4.mean[9]
  
  #this will output results? The order doesn't matter here because necessary changes have already been made
  #Use samp.fit since subsample
  #get values for each population and then average the populations
  pops=data_summary(samp.fit,"Relative.Fitness",c("Chromosome","Focal.Sex","Population"))
  data_summary(pops,"Relative.Fitness",c("Chromosome","Focal.Sex"))$Relative.Fitness
})
#GET ANCESTOR MEAN
#-------there are 27 combinations of CHR, Treatment, Sex and Block, corresponding to each row in the following order
#  [1] "BENIGN FL F 1" "BENIGN FL F 2" "BENIGN FL M 3" [4] "BENIGN FL M 4" "HIGH FL F 1"   "HIGH FL F 2"  
#  [7] "HIGH FL M 3"   "HIGH FL M 4"   "HIGH FL M 5"   [10] "BENIGN ML F 1" "BENIGN ML F 2" "BENIGN ML M 3"
# [13] "BENIGN ML M 4" "HIGH ML F 1"   "HIGH ML F 2"   [16] "HIGH ML M 3"   "HIGH ML M 4"   "HIGH ML M 5"  
# [19] "BENIGN T4 F 1" "BENIGN T4 F 2" "BENIGN T4 M 3" [22] "BENIGN T4 M 4" "HIGH T4 F 1"   "HIGH T4 F 2"  
# [25] "HIGH T4 M 3"   "HIGH T4 M 4"   "HIGH T4 M 5"  
#-----there are 9 Ancestors, and therefore 18 adjustments

#do this to get the mean Relative.Fit2 for each group which the above should be equivalent to. Make sure ORDER is same as above boot
Sex.Chr.plotting=data_summary(SSL.fit,"Relative.Fit2",c("Chromosome","Focal.Sex"))
#calculate all CI including median, just in case
Sex.Chr.CI=apply(Sex.Chr.boot, 1, quantile, c(0.025,0.975,0.5))
Sex.Chr.plotting$lower_CI=Sex.Chr.CI[1,]
Sex.Chr.plotting$upper_CI=Sex.Chr.CI[2,]
Sex.Chr.plotting$median=Sex.Chr.CI[3,]
#change the name of Focal.Sex variables for graphing
#library(plyr)
Sex.Chr.plotting$Focal.Sex<-revalue(Sex.Chr.plotting$Focal.Sex,c("F"="Females","M"="Males"))

##


#####-----GRAPH MAIN Relative FITNESS FIGURES: 2 x 2 comparison!----(use)------####
Fig1A=ggplot(data=Sex.Chr.plotting[1:4,],aes(fill=Chromosome,y=Relative.Fit2,x=Focal.Sex))+
  scale_fill_brewer(palette="Set1")+ylab("Mean Relative Fitness")+ggtitle("Relative Fit2, \nbootstrap CI subsample ancestor")+
  geom_bar(position=position_dodge(),stat = "summary", fun = "mean")+theme_classic() +
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.4,position=position_dodge(.9), colour="black",size=0.5)+
  geom_hline(yintercept=1,colour="black")+xlab(" ")+
  theme(text = element_text(size = 14,color="black"),legend.position="none")+
  theme(plot.title = element_text(size = 17,color="black")) + theme(axis.text.x = element_text(size = 12,color="black"))+
  theme(axis.text.y = element_text(size = 12,color="black"))

#####=====================------Relative FITNESS analysis with structure 2 x 2 x 2 x 3 (WIP)-----

Sex.Chr.Treat.Pop.boot=replicate(1000,{
   #Noting the the default sample is N, this should proportionally sample each subgroup
  #Benign Female Block 1
  B1.FFL1=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  B2.FFL1=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  B3.FFL1=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  B1.FML1=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  B2.FML1=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  B3.FML1=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  
  B.F1=append(append(append(B1.FFL1,B2.FFL1),append(B3.FFL1,B1.FML1)),append(B2.FML1,B3.FML1))

    #Beningn Female Block 2
  B1.FFL2=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  B2.FFL2=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  B3.FFL2=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  B1.FML2=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  B2.FML2=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  B3.FML2=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  
  B.F2=append(append(append(B1.FFL2,B2.FFL2),append(B3.FFL2,B1.FML2)),append(B2.FML2,B3.FML2))
  
  #BENIGN MALE BLOCK 3
  B1.MFL3=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  B2.MFL3=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  B3.MFL3=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  B1.MML3=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  B2.MML3=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  B3.MML3=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  
  B.M3=append(append(append(B1.MFL3,B2.MFL3),append(B3.MFL3,B1.MML3)),append(B2.MML3,B3.MML3))
  
  #BENIGN MALE BLOCK 3
  B1.MFL4=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  B2.MFL4=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  B3.MFL4=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  B1.MML4=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  B2.MML4=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  B3.MML4=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  
  B.M4=append(append(append(B1.MFL4,B2.MFL4),append(B3.MFL4,B1.MML4)),append(B2.MML4,B3.MML4))
  
  #HIGH Female Block 1
  H1.FFL1=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  H2.FFL1=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  H3.FFL1=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  H1.FML1=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  H2.FML1=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  H3.FML1=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  
  H.F1=append(append(append(H1.FFL1,H2.FFL1),append(H3.FFL1,H1.FML1)),append(H2.FML1,H3.FML1))
  
  #HIGH Female block 2
  H1.FFL2=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  H2.FFL2=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  H3.FFL2=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  H1.FML2=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  H2.FML2=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  H3.FML2=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  
  H.F2=append(append(append(H1.FFL2,H2.FFL2),append(H3.FFL2,H1.FML2)),append(H2.FML2,H3.FML2))
  
  #HIGH MALE BLOCK 3
  H1.MFL3=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  H2.MFL3=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  H3.MFL3=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  H1.MML3=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  H2.MML3=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  H3.MML3=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  
  H.M3=append(append(append(H1.MFL3,H2.MFL3),append(H3.MFL3,H1.MML3)),append(H2.MML3,H3.MML3))
  
  #HIGH MALE BLOCK 4
  H1.MFL4=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  H2.MFL4=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  H3.MFL4=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  H1.MML4=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  H2.MML4=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  H3.MML4=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  
  H.M4=append(append(append(H1.MFL4,H2.MFL4),append(H3.MFL4,H1.MML4)),append(H2.MML4,H3.MML4))
  
  #HIGH MALE BLOCK 5
  H1.MFL5=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="5","Rep.ID"],replace=TRUE)
  H2.MFL5=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="5","Rep.ID"],replace=TRUE)
  H3.MFL5=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="5","Rep.ID"],replace=TRUE)
  H1.MML5=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="5","Rep.ID"],replace=TRUE)
  H2.MML5=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="5","Rep.ID"],replace=TRUE)
  H3.MML5=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="5","Rep.ID"],replace=TRUE)
  
  H.M5=append(append(append(H1.MFL5,H2.MFL5),append(H3.MFL5,H1.MML5)),append(H2.MML5,H3.MML5))
  
  #saving these separate so easier to check if correct
  samp.B=append(append(B.F1,B.F2),append(B.M3,B.M4))
  samp.H=append(append(append(H.F1,H.F2),append(H.M3,H.M4)),H.M5)
  samp.BH=append(samp.B,samp.H)
  #put together the full sample in it's correct proportions
  samp=SSL.fit[SSL.fit$Rep.ID%in%samp.BH,]
  
  #calculate all the means for each block for T4. save them as variables to check
  #Benign Female Block 1
  T4.BF1=sample(SSL.fit[SSL.fit$Focal.Sex=="F"&SSL.fit$Treatment=="BENIGN"&SSL.fit$Chromosome=="T4"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  T4.BF1.mean=mean(SSL.fit[SSL.fit$Rep.ID%in%T4.BF1,"Fit2"])
  #Benign Female Block 2
  T4.BF2=sample(SSL.fit[SSL.fit$Focal.Sex=="F"&SSL.fit$Treatment=="BENIGN"&SSL.fit$Chromosome=="T4"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  T4.BF2.mean=mean(SSL.fit[SSL.fit$Rep.ID%in%T4.BF2,"Fit2"])
  #Benign Male Block 3
  T4.BM3=sample(SSL.fit[SSL.fit$Focal.Sex=="M"&SSL.fit$Treatment=="BENIGN"&SSL.fit$Chromosome=="T4"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  T4.BM3.mean=mean(SSL.fit[SSL.fit$Rep.ID%in%T4.BM3,"Fit2"])
  #Benign Male Block 4
  T4.BM4=sample(SSL.fit[SSL.fit$Focal.Sex=="M"&SSL.fit$Treatment=="BENIGN"&SSL.fit$Chromosome=="T4"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  T4.BM4.mean=mean(SSL.fit[SSL.fit$Rep.ID%in%T4.BM4,"Fit2"])
  
  #HIGH Blocks
  T4.HF1=sample(SSL.fit[SSL.fit$Focal.Sex=="F"&SSL.fit$Treatment=="HIGH"&SSL.fit$Chromosome=="T4"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  T4.HF1.mean=mean(SSL.fit[SSL.fit$Rep.ID%in%T4.HF1,"Fit2"])
  T4.HF2=sample(SSL.fit[SSL.fit$Focal.Sex=="F"&SSL.fit$Treatment=="HIGH"&SSL.fit$Chromosome=="T4"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  T4.HF2.mean=mean(SSL.fit[SSL.fit$Rep.ID%in%T4.HF2,"Fit2"])
  T4.HM3=sample(SSL.fit[SSL.fit$Focal.Sex=="M"&SSL.fit$Treatment=="HIGH"&SSL.fit$Chromosome=="T4"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  T4.HM3.mean=mean(SSL.fit[SSL.fit$Rep.ID%in%T4.HM3,"Fit2"])
  T4.HM4=sample(SSL.fit[SSL.fit$Focal.Sex=="M"&SSL.fit$Treatment=="HIGH"&SSL.fit$Chromosome=="T4"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  T4.HM4.mean=mean(SSL.fit[SSL.fit$Rep.ID%in%T4.HM4,"Fit2"])
  T4.HM5=sample(SSL.fit[SSL.fit$Focal.Sex=="M"&SSL.fit$Treatment=="HIGH"&SSL.fit$Chromosome=="T4"&SSL.fit$Block=="5","Rep.ID"],replace=TRUE)
  T4.HM5.mean=mean(SSL.fit[SSL.fit$Rep.ID%in%T4.HM5,"Fit2"])
  
  
  ##For each sample, make sure the appropriate Fit2 is turned into Relative.Fitness
  samp$Relative.Fitness=0
  samp[samp$Rep.ID%in%B.F1,"Relative.Fitness"]= samp[samp$Rep.ID%in%B.F1,"Fit2"]/T4.BF1.mean
  samp[samp$Rep.ID%in%B.F2,"Relative.Fitness"]= samp[samp$Rep.ID%in%B.F2,"Fit2"]/T4.BF2.mean
  samp[samp$Rep.ID%in%B.M3,"Relative.Fitness"]= samp[samp$Rep.ID%in%B.M3,"Fit2"]/T4.BM3.mean
  samp[samp$Rep.ID%in%B.M4,"Relative.Fitness"]= samp[samp$Rep.ID%in%B.M4,"Fit2"]/T4.BM4.mean
  samp[samp$Rep.ID%in%H.F1,"Relative.Fitness"]= samp[samp$Rep.ID%in%H.F1,"Fit2"]/T4.HF1.mean
  samp[samp$Rep.ID%in%H.F2,"Relative.Fitness"]= samp[samp$Rep.ID%in%H.F2,"Fit2"]/T4.HF2.mean
  samp[samp$Rep.ID%in%H.M3,"Relative.Fitness"]= samp[samp$Rep.ID%in%H.M3,"Fit2"]/T4.HM3.mean
  samp[samp$Rep.ID%in%H.M4,"Relative.Fitness"]= samp[samp$Rep.ID%in%H.M4,"Fit2"]/T4.HM4.mean
  samp[samp$Rep.ID%in%H.M5,"Relative.Fitness"]= samp[samp$Rep.ID%in%H.M5,"Fit2"]/T4.HM5.mean  
  #this will output results? The order doesn't matter here because necessary changes have already been made
  #All of the raw Fit2 has now been strandardized by appropriate control, so can summarize the MEAN
  by.pops.first=data_summary(samp,"Relative.Fitness",c("Chromosome","Treatment","Focal.Sex","Population"))
  data_summary(by.pops.first,"Relative.Fitness",c("Chromosome","Treatment","Focal.Sex"))$Relative.Fitness
})

#note that even though I standardize by pop I don't take the mean by pop

#do this to get the mean Relative.Fit2 for each group which the above should be equivalent to. Make sure ORDER is same as above boot
Sex.Chr.Treat.Pop.plotting=data_summary(SSL.fit[SSL.fit$Chromosome!="T4",],"Relative.Fit2",c("Chromosome","Treatment","Focal.Sex"))
#calculate all CI including median, just in case
Sex.Chr.Treat.Pop.CI=apply(Sex.Chr.Treat.Pop.boot, 1, quantile, c(0.025,0.975,0.5))
Sex.Chr.Treat.Pop.plotting$mean=apply(Sex.Chr.Treat.Pop.boot, 1, mean)
Sex.Chr.Treat.Pop.plotting$lower_CI=Sex.Chr.Treat.Pop.CI[1,]
Sex.Chr.Treat.Pop.plotting$upper_CI=Sex.Chr.Treat.Pop.CI[2,]
Sex.Chr.Treat.Pop.plotting$median=Sex.Chr.Treat.Pop.CI[3,]
#change the name of Focal.Sex variables for graphing
#library(plyr)
Sex.Chr.Treat.Pop.plotting$Focal.Sex<-revalue(Sex.Chr.Treat.Pop.plotting$Focal.Sex,c("F"="Females","M"="Males"))


#####-----GRAPH RelativeFITNESS FIGURES: 2 x 2 x 2 x 3 comparison (WIP)------####

#there is disagreement between the mean and CI here for MANY of the ones. it makes HMales more sim and Bmales lower
Fig1B=ggplot(data=Sex.Chr.Treat.Pop.plotting[1:8,],aes(fill=Chromosome,y=median,x=Focal.Sex,colour=Chromosome))+
  scale_fill_brewer(palette="Set1")+ylab("Relative Mean Fitness")+ggtitle("Relative Fit2, \nbootstrap CI subsample ancestor")+
  geom_bar(position=position_dodge(),stat = "summary", fun = "mean")+theme_classic() +
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,position=position_dodge(.9), colour="black")+
  facet_wrap(~Treatment)+
  geom_hline(yintercept=1.0,colour="black")+
  theme(text = element_text(size = 14,color="black"),legend.position="right")+
  theme(plot.title = element_text(size = 17,color="black")) + theme(axis.text.x = element_text(size = 12,color="black"))+
  theme(axis.text.y = element_text(size = 12,color="black"))



#####===========================================---->-Mortality and HARM (Ancestor Female Fecundity)--####
#-------ALTERNATIVE ANCESTRAL FECUNDITY MATE HARM (bootstrap together)----->(use)----####

fecundity.boot=replicate(1000,{
  #pull out replicates
  samp<-sample(SSL.fit[SSL.fit$Focal.Sex=="M","Rep.ID"],replace=TRUE)
  data_summary(SSL.fit[SSL.fit$Rep.ID%in%samp,],"Avg.Offspring.per.F",c("Chromosome","Treatment"))$Avg.Offspring.per.F
  #data_summary(SSL.fit[SSL.fit$Rep.ID%in%samp,],"Offspring.Total",c("Chromosome","Treatment"))$Offspring.Total
})
fecundity.CI=apply(fecundity.boot, 1, quantile, c(0.025,0.975) )
fecundity.error=data_summary(subset(SSL.fit,Focal.Sex=="M"),"Avg.Offspring.per.F",c("Chromosome","Treatment"))
#fecundity.error=data_summary(subset(SSL.fit,Focal.Sex=="M"),"Offspring.Total",c("Chromosome","Treatment"))
fecundity.error$lower_CI=fecundity.CI[1,]
fecundity.error$upper_CI=fecundity.CI[2,]

mate.harm.fec.plot=ggplot(data=fecundity.error[1:4,],aes(fill=Chromosome,y=Avg.Offspring.per.F,x=c("BENIGN","NOVEL","BENIGN","NOVEL"),colour=Chromosome))+
  scale_fill_brewer(palette="Set1")+ylab("Grand Ancestral Female Offspring per Female")+
  geom_bar(position=position_dodge(),stat = "summary", fun = "mean")+theme_classic() +xlab("Thermal Regime")+
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,position=position_dodge(.9), colour="black") 


#-----Alternative way Moraltiy (Harm) (bootstrap together)---------------->(use)----####
mortality.boot=replicate(1000,{
  #pull out replicates
  samp<-sample(NoA[NoA$Focal.Sex=="M","Rep.ID"],replace=TRUE)
  data_summary(NoA[NoA$Rep.ID%in%samp,],"Lay.Ancestor.Death",c("Chromosome","Treatment"))$Lay.Ancestor.Death
})
Mortality.error=data_summary(subset(NoA,Focal.Sex=="M"),"Lay.Ancestor.Death",c("Chromosome","Treatment"))
Mortality.error$Ancestor.Survival=12-Mortality.error$Lay.Ancestor.Death

mortality.CI=apply(mortality.boot, 1, quantile, c(0.025,0.975) )

Mortality.error$lower_CI=12-mortality.CI[1,]
Mortality.error$upper_CI=12-mortality.CI[2,]

#mate.harm.fec.plot
mate.harm.surv.plot=ggplot(data=Mortality.error,aes(fill=Chromosome,y=(Ancestor.Survival),x=c("BENIGN","NOVEL","BENIGN","NOVEL"),colour=Chromosome))+
  scale_fill_brewer(palette="Set1")+ylab("Female Grand Ancestor Survival")+xlab("Thermal Regime")+
  geom_bar(position=position_dodge(),stat = "summary", fun = "mean")+
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,position=position_dodge(.9), colour="black")+
  theme_classic()+theme(plot.title = element_text(size = 17,color="black"), axis.text.x = element_text(size = 12,color="black"))+
  theme(axis.text.y = element_text(size = 12,color="black"),legend.position="none")
#make the aesthetic match
mate.harm.fec.plot=mate.harm.fec.plot+theme(axis.text.x = element_text(size = 12,color="black"), 
axis.text.y = element_text(size = 12,color="black"),legend.margin=margin(-10,0,0,-15))#legend.key.size = unit(1, "lines"))

ggarrange(mate.harm.surv.plot,mate.harm.fec.plot,labels=c("A","B"),ncol=2,widths=c(0.8,1.1))

######------MATE HARM CORR WITH MALE SUCCESS CORR GRAPH---(second)--(use)----####

Per.Line.Fit$Offspring.Total.GroupNorm=0
#Per.Line.Fit=GroupNorm_func(Per.Line.Fit,"Offspring.Total.GroupNorm","Offspring.Total","M")
#doing this only for focal offspring but being too lazy to change the code in later sections
Per.Line.Fit=GroupNorm_func(Per.Line.Fit,"Offspring.Total.GroupNorm","Offspring.Total","M")

BFL.Adboot=replicate(1000,{
  #pull out replicates per population. The default size is n. So it might be uncessary to calculate length?
  samp.B1=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="B1","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="B1","Line.ID"]),replace=TRUE)#32
  samp.B2=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="B2","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="B2","Line.ID"]),replace=TRUE)#30
  samp.B3=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="B3","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="B3","Line.ID"]),replace=TRUE)#29
  
  samp.BFL=append(samp.B1,samp.B2)
  samp.BFL=append(samp.BFL,samp.B3)
  
  #rho is value four
  BFL.rmf=cor.test(Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.BFL,"GroupNorm.Fit.M"],Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.BFL,"Offspring.Total.GroupNorm"])
  #rho is value four
  #BFL.rmf[4]#this will keep that it's rho
  BFL.rmf[[4]]#this is just value
})
BML.Adboot=replicate(1000,{
  #pull out replicates (all)
  
  samp.B1=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="B1","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="B1","Line.ID"]),replace=TRUE)#32
  samp.B2=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="B2","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="B2","Line.ID"]),replace=TRUE)#30
  samp.B3=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="B3","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="B3","Line.ID"]),replace=TRUE)#29
  
  samp.BML=append(samp.B1,samp.B2)
  samp.BML=append(samp.BML,samp.B3)
  
  #rho is value four
  BML.rmf=cor.test(Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.BML,"GroupNorm.Fit.M"],Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.BML,"Offspring.Total.GroupNorm"])
  #rho is value four
  #BFL.rmf[4]#this will keep that it's rho
  BML.rmf[[4]]#this is just value
})

HFL.Adboot=replicate(1000,{
  #pull out replicates (all)
  
  samp.H1=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="H1","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="H1","Line.ID"]),replace=TRUE)#32
  samp.H2=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="H2","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="H2","Line.ID"]),replace=TRUE)#30
  samp.H3=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="H3","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="H3","Line.ID"]),replace=TRUE)#29
  
  samp.HFL=append(samp.H1,samp.H2)
  samp.HFL=append(samp.HFL,samp.H3)
  
  #rho is value four
  HFL.rmf=cor.test(Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.HFL,"GroupNorm.Fit.M"],Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.HFL,"Offspring.Total.GroupNorm"])
  #rho is value four
  #BFL.rmf[4]#this will keep that it's rho
  HFL.rmf[[4]]#this is just value
})
HML.Adboot=replicate(1000,{
  #pull out replicates (all)
  
  samp.H1=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="H1","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="H1","Line.ID"]),replace=TRUE)#32
  samp.H2=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="H2","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="H2","Line.ID"]),replace=TRUE)#30
  samp.H3=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="H3","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="H3","Line.ID"]),replace=TRUE)#29
  
  samp.HML=append(samp.H1,samp.H2)
  samp.HML=append(samp.HML,samp.H3)
  
  #rho is value four
  HML.rmf=cor.test(Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.HML,"GroupNorm.Fit.M"],Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.HML,"Offspring.Total.GroupNorm"])
  #rho is value four
  #BFL.rmf[4]#this will keep that it's rho
  HML.rmf[[4]]#this is just value
})

#compare the distributions
#can do a corr of corrs? or cor1-cor2?
##start with ML. ##does the pairing matter?????
##if the corrs are the same the difference will be 0. If there is no 0 diff, then SIG DIFF
diff.ML=HML.Adboot-BML.Adboot
diff.ML.CI=quantile(diff.ML,c(0.975,0.025))# interval(0.5441168 <-> 0.1124657) SIG
diff.FL=HFL.Adboot-BFL.Adboot
diff.FL.CI=quantile(diff.FL,c(0.975,0.025))#interval (0.414297915 <-> -0.003945277) NOT SIG

#now get the average per treatment
Ball.Adboot=(BFL.Adboot+BML.Adboot)/2
Hall.Adboot=(HFL.Adboot+HML.Adboot)/2

diff.all=Hall.Adboot-Ball.Adboot
diff.all.CI=quantile(diff.all,c(0.975,0.025))

plotting_values=data.frame(mean=c(mean(BFL.Adboot),mean(BML.Adboot,na.rm=TRUE),mean(HFL.Adboot),mean(HML.Adboot),mean(Hall.Adboot),mean(Ball.Adboot)),
                           Chromosome=c("BFL","BML","HFL","HML","Havg","Bavg"),
                           upper_CI=c(quantile(BFL.Adboot,0.975),quantile(BML.Adboot,0.975,na.rm=TRUE),quantile(HFL.Adboot,0.975),quantile(HML.Adboot,0.975),quantile(Hall.Adboot,0.975),quantile(Ball.Adboot,0.975)),
                           lower_CI=c(quantile(BFL.Adboot,0.025),quantile(BML.Adboot,0.025,na.rm=TRUE),quantile(HFL.Adboot,0.025),quantile(HML.Adboot,0.025),quantile(Hall.Adboot,0.025),quantile(Ball.Adboot,0.025)))
#overall=c(BFL.corr[[4]],BML.corr[[4]],HML.corr[[4]],HFL.corr[[4]]))

#examine if increase mortality confers success
mate.harm.success.plot=ggplot(plotting_values,aes(x=Chromosome,y=as.numeric(mean),fill=Chromosome,colour=Chromosome))+
  geom_point(size=3,alpha=.9,shape=21)+ 
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI,colour=Chromosome), width=.4,position=position_dodge(.9))+
  scale_fill_manual(values=c("black","red","#3e86ff","black","red","#3e86ff"))+
  scale_colour_manual(values=c("black","red","#3e86ff","black","red","#3e86ff"))+#"#f6564b
  geom_vline(xintercept=3.5,colour="grey",linewidth=0.5)+ # ggtitle("groupNorm Ancestor Fitness w/\n Male Competitive fit")
  ylab("correlation of Grand Ancestor # Focal individuals and \nFocal Male Fitness ")+xlab("Chromosome Pool")+
  theme_classic() +theme(text = element_text(size = 12,color="black"),legend.position="none",panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7))+
  theme(axis.text.x = element_text(size = 12,color="black"), axis.text.y = element_text(size = 12,color="black"))

#mate.harm.fec.plot=mate.harm.fec.plot+ylab("Grand Ancestor Offspring \n(per Fem)")

#ggarrange(ggarrange(mate.harm.surv.plot,mate.harm.fec.plot,labels=c("A","B"),ncol=2,widths=c(1,1),common.legend=TRUE),
#      mate.harm.success.plot,nrow=2)



########-----(use)randomize Corr Grand Ancestor Fitness and Focal Male Fitness within pop blocks (use)####
##so within a block scramble focal and calculate new total
SSL.fit$pop.block=paste(SSL.fit$Block, SSL.fit$Population)
#have to create the new column before calling the function
SSL.fit$meanGroupNorm.Fit=0
SSL.fit$Fake.Total.Offspring=0
#the function won't modify the original data.frame so save new->#this is the new normalized
#normalizing removes T4
SSL.fit$Real.meanGroupNorm.Fit=0
SSL.fit=GroupNorm_func(SSL.fit,"Real.meanGroupNorm.Fit","Relative.Fit2",c("F","M"))
#only care about Male fitness correlation
M.fit=subset(SSL.fit,Focal.Sex=="M")
#order is NML, NFL, BML, BFL #replicate 1000 times
perF=FALSE
Fake.cor=replicate(1000,{
  #sample per block and generate one randomized go through
  for(pop.block in unique(M.fit$pop.block)){
    #print(pop.block)
    #pull out ML and FL separately for each pop and block (pops are specific to treatment)
    temp.ML.ID=M.fit$pop.block==pop.block&M.fit$Chromosome=="ML"
    temp.FL.ID=M.fit$pop.block==pop.block&M.fit$Chromosome=="FL"
    temp.T4.ID=M.fit$Block==str_split(pop.block," ")[[1]][1]&M.fit$Chromosome=="T4"
    #I sample them separately to make sure I am swaping with ML and FL (per pop)
    temp.ML.focal=sample(M.fit[temp.ML.ID,"Offspring.Focal"])
    temp.FL.focal=sample(M.fit[temp.FL.ID,"Offspring.Focal"])
    #temp.T4.focal=sample(M.fit[temp.T4.ID,"Offspring.Focal"])
    temp.ML.comp=sample(M.fit[temp.ML.ID,"Offspring.Competitor"])
    temp.FL.comp=sample(M.fit[temp.FL.ID,"Offspring.Competitor"])
    #temp.T4.comp=sample(M.fit[temp.T4.ID,"Offspring.Competitor"])
    M.fit[temp.ML.ID,"Fake.Total.Offspring"]=temp.ML.comp+temp.ML.focal
    M.fit[temp.FL.ID,"Fake.Total.Offspring"]=temp.FL.comp+temp.FL.focal
    M.fit[temp.ML.ID,"Fake.Offspring.perF"]= M.fit[temp.ML.ID,"Fake.Total.Offspring"]/(12-M.fit[temp.ML.ID,"Lay.Ancestor.Death"])
    M.fit[temp.FL.ID,"Fake.Offspring.perF"]= M.fit[temp.FL.ID,"Fake.Total.Offspring"]/(12-M.fit[temp.FL.ID,"Lay.Ancestor.Death"])
    
    #M.fit[temp.T4.ID,"Fake.Total.Offspring"]=temp.T4.comp+temp.T4.focal
    M.fit[temp.ML.ID,"Fake.Fit2"]=temp.ML.focal/M.fit[temp.ML.ID,"Fake.Total.Offspring"]
    M.fit[temp.FL.ID,"Fake.Fit2"]=temp.FL.focal/M.fit[temp.FL.ID,"Fake.Total.Offspring"]
    #M.fit[temp.T4.ID,"Fake.Fit2"]=temp.T4.focal/M.fit[temp.T4.ID,"Fake.Total.Offspring"]
    #I should use Fake.Fit2 here as well
    Ancestor.mean.fit=mean(M.fit[temp.T4.ID,"Fit2"])
    #Ancestor.mean.fit=mean(M.fit[temp.T4.ID,"Fake.Fit2"])
  }
  #normalize the singular randomized go through
  M.fit=GroupNorm_func(M.fit,"normFake.nT","Fake.Total.Offspring",c("M"))
  if(perF){
    M.fit=GroupNorm_func(M.fit,"normFake.nT","Fake.Offspring.perF",c("M"))}
  #this will replace meanGroupNorm.Fit with the Randomized version for males (not changing column name for ease)
  M.fit=GroupNorm_func(M.fit,"meanGroupNorm.Fit","Fake.Fit2",c("M"))
  #M.fit=GroupNorm_func(M.fit,"meanGroupNorm.Fit","Relative.Fake.Fit2",c("M"))
  
  #find the correlation and save that as one value
  #note that the meanGroupNorm.Fit will always stay the same and it's the others that change
  
  c(as.numeric(cor.test(M.fit[M.fit$Chromosome=="ML"&M.fit$Treatment=="HIGH","meanGroupNorm.Fit"],
                        M.fit[M.fit$Chromosome=="ML"&M.fit$Treatment=="HIGH","normFake.nT"])[[4]]),
    as.numeric(cor.test(M.fit[M.fit$Chromosome=="FL"&M.fit$Treatment=="HIGH","meanGroupNorm.Fit"],
                        M.fit[M.fit$Chromosome=="FL"&M.fit$Treatment=="HIGH","normFake.nT"])[[4]]),
    as.numeric(cor.test(M.fit[M.fit$Chromosome=="ML"&M.fit$Treatment=="BENIGN","meanGroupNorm.Fit"],
                        M.fit[M.fit$Chromosome=="ML"&M.fit$Treatment=="BENIGN","normFake.nT"])[[4]]),
    as.numeric(cor.test(M.fit[M.fit$Chromosome=="FL"&M.fit$Treatment=="BENIGN","meanGroupNorm.Fit"],
                        M.fit[M.fit$Chromosome=="FL"&M.fit$Treatment=="BENIGN","normFake.nT"])[[4]]))
})

Bavg=(Fake.cor[3,]+Fake.cor[4,])/2
Navg=(Fake.cor[1,]+Fake.cor[2,])/2
Fake.cor=rbind(Fake.cor,Bavg)
Fake.cor=rbind(Fake.cor,Navg)
lower_CI=apply(Fake.cor,1,quantile,0.025)
upper_CI=apply(Fake.cor,1,quantile,0.975)
median=apply(Fake.cor,1,quantile,0.5)

###now normalize actual total offspring
M.fit=GroupNorm_func(M.fit,"norm.Offspring.Total","Offspring.Total",c("M"))
if(perF){
  M.fit=GroupNorm_func(M.fit,"norm.Offspring.Total","Avg.Offspring.per.F",c("M"))}
#M.fit=GroupNorm_func(SSL.fit,"Real.meanGroupNorm.Fit","Relative.Fit2",("M"))

real.means=c(as.numeric(cor.test(M.fit[M.fit$Chromosome=="ML"&M.fit$Treatment=="HIGH","Real.meanGroupNorm.Fit"],
                                 M.fit[M.fit$Chromosome=="ML"&M.fit$Treatment=="HIGH","norm.Offspring.Total"])[[4]]),
             as.numeric(cor.test(M.fit[M.fit$Chromosome=="FL"&M.fit$Treatment=="HIGH","Real.meanGroupNorm.Fit"],
                                 M.fit[M.fit$Chromosome=="FL"&M.fit$Treatment=="HIGH","norm.Offspring.Total"])[[4]]),
             as.numeric(cor.test(M.fit[M.fit$Chromosome=="ML"&M.fit$Treatment=="BENIGN","Real.meanGroupNorm.Fit"],
                                 M.fit[M.fit$Chromosome=="ML"&M.fit$Treatment=="BENIGN","norm.Offspring.Total"])[[4]]),
             as.numeric(cor.test(M.fit[M.fit$Chromosome=="FL"&M.fit$Treatment=="BENIGN","Real.meanGroupNorm.Fit"],
                                 M.fit[M.fit$Chromosome=="FL"&M.fit$Treatment=="BENIGN","norm.Offspring.Total"])[[4]]),
             #these are the Bavg and Navg which is just compiling them
             as.numeric(cor.test(M.fit[M.fit$Treatment=="BENIGN","Real.meanGroupNorm.Fit"],
                                 M.fit[M.fit$Treatment=="BENIGN","norm.Offspring.Total"])[[4]]),
             as.numeric(cor.test(M.fit[M.fit$Treatment=="HIGH","Real.meanGroupNorm.Fit"],
                                 M.fit[M.fit$Treatment=="HIGH","norm.Offspring.Total"])[[4]]))

plotting_values=data.frame(Chromosome=c("ML","FL","ML","FL","Z","Z"),median=median,lower_CI=lower_CI,upper_CI=upper_CI,real.means=real.means)
plotting_values$Treatment=c("NOVEL","NOVEL","BENIGN","BENIGN","BENIGN","NOVEL")

##way to graph as points
#examine if increase mortality confers success
Fig3C=ggplot(plotting_values[1:4,],aes(x=Chromosome,y=as.numeric(median)))+
  geom_point(size=3,alpha=.7,shape=21)+ 
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI,colour="darkgrey"), width=.4,position=position_dodge(.9))+
  facet_wrap(~Treatment,scales="free_x")+
  ylab("correlation of Female Mating Partner Fitness and \nFocal Male Competitive Fitness")+xlab("Chromosome Pool")+
  #ggtitle("random nT and wM w/in pop/block \n+ Fake.Fit2 + Fake.Offspring.perF")+
  geom_point(inherit.aes=FALSE,aes(y=real.means,x=Chromosome,colour=Chromosome,fill=Chromosome),size=2.5,shape=16)+
  scale_fill_manual(values=c("black","red","#3e86ff","darkgrey","red","#3e86ff"))+
  #scale_colour_manual(values=c("black","purple","red","darkgrey","#3e86ff","orange","green"))+
  scale_colour_manual(values=c("darkgrey","red","#3e86ff","darkgrey","red","#3e86ff"))+ #"#f6564b
  theme_classic() +theme(text = element_text(size = 12,color="black"),legend.position="none",panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7))+
  theme(axis.text.x = element_text(size = 12,color="black"), axis.text.y = element_text(size = 12,color="black"))

#Fig3C
##Now doing it with a box plot to look visually somewhere inbetween

plotting_values3=data.frame(Treatment=c(rep("NOVEL",2000),rep("BENIGN",2000)),
                            Chromosome=c(rep("ML",1000),rep("FL",1000),rep("ML",1000),rep("FL",1000)),
                            corrs=append(append(Fake.cor[1,],Fake.cor[2,]),append(Fake.cor[3,],Fake.cor[4,])))

AltFig3C=  ggplot(plotting_values3,aes(x=Chromosome,y=corrs,colour="grey",fill="lightgrey"))+
  geom_boxplot(notch=TRUE,outliers=FALSE,position=position_dodge2(.5))+ylim(-0.35,0.15)+
  ylab("Pearson correlation of \nFemale Mating Partner Fitness and \nFocal Male Competitive Fitness")+xlab("Chromosome Pool")+
  facet_wrap(~Treatment,scale="free_x") +
  geom_point(inherit.aes=FALSE,data=plotting_values[1:4,],aes(y=real.means,x=Chromosome,colour=Chromosome,fill=Chromosome),size=2.5,shape=16)+
  geom_point(inherit.aes=FALSE,data=plotting_values[1:4,],aes(y=real.means,x=Chromosome,colour=Chromosome,fill=Chromosome),size=5,shape=4)+
  scale_fill_manual(values=c(alpha("red",1),alpha("lightgrey",0.7),alpha("#3e86ff",1)))+
  scale_colour_manual(values=c("red","black","#3e86ff","darkgrey","orange","green"))+
  #scale_colour_manual(values=c("darkgrey","red","#3e86ff","darkgrey","red","#3e86ff"))+ #"#f6564b
  
  theme_classic() +theme(text = element_text(size = 12,color="black"),legend.position="none",panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7))+
  theme(axis.text.x = element_text(size = 12,color="black"), axis.text.y = element_text(size = 12,color="black"))



# ###alternative way to graph for mate harm
# plotting_values3=data.frame(Treatment=c(rep("NOVEL",2000),rep("BENIGN",2000)), 
#                 Chromosome=c(rep("ML",1000),rep("FL",1000),rep("ML",1000),rep("FL",1000)), 
#                 corrs=append(append(Fake.cor[1,],Fake.cor[2,]),append(Fake.cor[3,],Fake.cor[4,])))
# 
# ggplot(plotting_values3,aes(x=Chromosome,y=corrs,fill="darkgrey"))+
#   geom_violin(draw_quantiles=TRUE)+
#   facet_wrap(~Treatment,scale="free_x") +
#   geom_point(inherit.aes=FALSE,data=plotting_values[1:4,],aes(y=real.means,x=Chromosome,colour=Chromosome,fill=Chromosome),size=2.5,shape=16)+
#   geom_point(inherit.aes=FALSE,data=plotting_values[1:4,],aes(y=real.means,x=Chromosome,colour=Chromosome,fill=Chromosome),size=5,shape=4)+
#   scale_fill_manual(values=c("grey","red","#3e86ff","grey","red","#3e86ff"))+ #"#f6564b
#   scale_colour_manual(values=c("red","#3e86ff","red","#3e86ff"))+ #"#f6564b
#   
#   theme_classic() +theme(text = element_text(size = 12,color="black"),legend.position="none",panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7))+
#   theme(axis.text.x = element_text(size = 12,color="black"), axis.text.y = element_text(size = 12,color="black"))



####================---------->FOCAL FEMALE FECUNDITY----####
#bootstrap separately or together
#-----Alternative way for Focal Female Fecundity (bootstrap together)----->(use)----####
fecundity.boot=replicate(1000,{
  #pull out replicates
  samp<-sample(SSL.fit[SSL.fit$Focal.Sex=="F","Rep.ID"],replace=TRUE)
  data_summary(SSL.fit[SSL.fit$Rep.ID%in%samp,],"Offspring.Competitor",c("Chromosome","Treatment"))$Offspring.Competitor
})
fecundity.CI=apply(fecundity.boot, 1, quantile, c(0.025,0.975) )
fecundity.error=data_summary(subset(SSL.fit,Focal.Sex=="F"),"Offspring.Competitor",c("Chromosome","Treatment"))
fecundity.error$lower_CI=fecundity.CI[1,]
fecundity.error$upper_CI=fecundity.CI[2,]

ggplot(data=fecundity.error,aes(fill=Chromosome,y=Offspring.Competitor,x=Treatment,colour=Chromosome))+
  scale_fill_brewer(palette="Set1")+ylab("Competitor Female Offspring (F)")+
  geom_bar(position=position_dodge(),stat = "summary", fun = "mean")+theme_classic() +
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,position=position_dodge(.9), colour="black") 

mortality.boot=replicate(1000,{
  #pull out replicates
  samp<-sample(NoA[NoA$Focal.Sex=="F","Rep.ID"],replace=TRUE)
  data_summary(NoA[NoA$Rep.ID%in%samp,],"Focal.Female.Death",c("Chromosome","Treatment"))$Focal.Female.Death
})
Mortality.error=data_summary(subset(NoA,Focal.Sex=="F"),"Focal.Female.Death",c("Chromosome","Treatment"))
Mortality.error$Focal.Female.Survival=Mortality.error$Focal.Female.Death

mortality.CI=apply(mortality.boot, 1, quantile, c(0.025,0.975) )
#Note Focal.Female Death is NOT the number of dead but the number of living
Mortality.error$lower_CI=mortality.CI[1,]
Mortality.error$upper_CI=mortality.CI[2,]

ggplot(data=Mortality.error,aes(fill=Chromosome,y=(Focal.Female.Survival),x=Treatment,colour=Chromosome))+
  scale_fill_brewer(palette="Set1")+ylab("Focal Female Survival (F)")+
  geom_bar(position=position_dodge(),stat = "summary", fun = "mean")+
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,position=position_dodge(.9), colour="black")+
  theme_classic()+theme(plot.title = element_text(size = 17,color="black")) + theme(axis.text.x = element_text(size = 12,color="black"))+
  theme(axis.text.y = element_text(size = 12,color="black"))

comp.mortality.boot=replicate(1000,{
  #pull out replicates
  samp<-sample(NoA[NoA$Focal.Sex=="F","Rep.ID"],replace=TRUE)
  data_summary(NoA[NoA$Rep.ID%in%samp,],"Competitor.Female.Death",c("Chromosome","Treatment"))$Competitor.Female.Death
})
Mortality.error=data_summary(subset(NoA,Focal.Sex=="F"),"Competitor.Female.Death",c("Chromosome","Treatment"))
Mortality.error$Competitor.Female.Survival=Mortality.error$Competitor.Female.Death

mortality.CI=apply(comp.mortality.boot, 1, quantile, c(0.025,0.975) )

Mortality.error$lower_CI=mortality.CI[1,]
Mortality.error$upper_CI=mortality.CI[2,]

ggplot(data=Mortality.error,aes(fill=Chromosome,y=(Competitor.Female.Survival),x=Treatment,colour=Chromosome))+
  scale_fill_brewer(palette="Set1")+ylab("Competitor Female Survival (F)")+
  geom_bar(position=position_dodge(),stat = "summary", fun = "mean")+
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,position=position_dodge(.9), colour="black")+
  theme_classic()+theme(plot.title = element_text(size = 17,color="black")) + theme(axis.text.x = element_text(size = 12,color="black"))+
  theme(axis.text.y = element_text(size = 12,color="black"))

####===================================================------->MODEL---GENERAL COUNT BASED QUASI-BINOMIAL FITNESS MODEL----####
#Model to use
NoA.test=glmer(cbind(Offspring.Focal,Offspring.Competitor)~Chromosome*Focal.Sex*Treatment+(1|Block)+(1|Population)+(1|Rep.ID),data=subset(NoA[!is.na(NoA$Offspring.Focal),],Treatment!="T4"),family=binomial)#
Anova(NoA.test,type="III")

summary(NoA.test)
# NoA.test1=glmer(cbind(Offspring.Focal,Offspring.Competitor)~Chromosome*Focal.Sex*Treatment+Block+Population+(1|Rep.ID),data=subset(NoA[!is.na(NoA$Offspring.Focal),],Treatment!="T4"),family=binomial)#
# Anova(NoA.test1,type="III")# fixed effect model matrix is rank deficient so dropping 2 columns/co-eff. 
#this is because Populations and Block are combinations of other fixed effects which would require nesting

# NoA.test1=glmer(cbind(Offspring.Focal,Offspring.Competitor)~Chromosome*Focal.Sex*Treatment+Population+(1|Rep.ID),data=subset(NoA[!is.na(NoA$Offspring.Focal),],Treatment!="T4"),family=binomial)#
# Anova(NoA.test1,type="III")# fixed effect model matrix is rank deficient so dropping 1 columns/co-eff. 

#~~~~~~~Luke Holman Request Model~~~~~#
#Model to run as requested by Luke Holman
#According to him Model <- glmer(cbind(n_focal, n_non_focal) ~ sex * treatment +(1 | RANDOM), family = “quasibinomial”, data = my_data)
#I think treatment in above means Chromosome. He also mentions Block and replicate ID so I'm sure these are all fine
Luke.Holman=glmer(cbind(Offspring.Focal,Offspring.Competitor)~Chromosome*Focal.Sex+(1|Block)+(1|Population)+(1|Rep.ID),data=subset(NoA[!is.na(NoA$Offspring.Focal),],Treatment!="T4"),family=binomial)#
Anova(Luke.Holman,type="III")

library(emmeans)
Simple.LH=emmeans(Luke.Holman, ~ Chromosome * Focal.Sex)
Regular=emmeans(NoA.test, ~ Chromosome * Focal.Sex * Treatment)

contrast(Regular, method = "pairwise")
contrast(Simple.LH, method = "pairwise")


##----Submodels for subsampling

##Focal Female Fecundity
focal.fem.fecund=lmer(Offspring.Focal~Chromosome*Treatment+Focal.Female.Death+(1|Block)+(1|Population)+(1|Rep.ID),data=subset(NoA[!is.na(NoA$Offspring.Focal),],Focal.Sex=="F"&Population!="T4"))#
Anova(focal.fem.fecund,type="III")
#Compare change in competitor offspring but only between T4 and FL
comp.fem.fecund=lmer(Offspring.Competitor~Chromosome*Treatment+Competitor.Female.Death+(1|Block)+(1|Rep.ID),data=subset(SSL.fit[!is.na(SSL.fit$Offspring.Competitor),],Focal.Sex=="F"&Chromosome!="ML"))##
Anova(comp.fem.fecund,type="III")

#Benign only comparison FL to Ancestor
B.focal.fem.fecund=lmer(Offspring.Focal~Chromosome+Focal.Female.Death+(1|Block), data=subset(SSL.fit[!is.na(SSL.fit$Offspring.Focal),],Focal.Sex=="F"&Treatment=="BENIGN"&Chromosome!="ML"))#
Anova(B.focal.fem.fecund,type="III")

##FEMALE MORTALITY IN FEMALE FITNESS
#comparing dead vs alive and what effects that number (not percent to total)
focalF.mortality=glmer(cbind(Focal.Female.Death,6-Focal.Female.Death)~Chromosome*Treatment+ Block+(1|Population)+(1|Rep.ID), data=subset(SSL.fit[!is.na(SSL.fit$Offspring.Focal),],Focal.Sex=="F"&Chromosome!="T4"&Population!="T4"),family=binomial)#
Anova(focalF.mortality,type="III")

compF.mortality=glmer(cbind(Competitor.Female.Death,6-Competitor.Female.Death)~Chromosome*Treatment+ Block+(1|Population)+(1|Rep.ID), data=subset(SSL.fit[!is.na(SSL.fit$Offspring.Focal),],Focal.Sex=="F"&Chromosome!="T4"&Population!="T4"),family=binomial)#
Anova(compF.mortality,type="III")
#non-binomial trial
compF.mortality=lmer(Competitor.Female.Death ~Chromosome*Treatment+ Block+(1|Population)+(1|Rep.ID), data=subset(SSL.fit[!is.na(SSL.fit$Offspring.Focal),],Focal.Sex=="F"&Chromosome!="T4"&Population!="T4"))#
Anova(compF.mortality,type="III")

#the balance of who dies, similar to fecudunity
netF.mortality=glmer(cbind(Competitor.Female.Death,Focal.Female.Death)~Chromosome*Treatment+ Block+(1|Population)+(1|Rep.ID), data=subset(SSL.fit[!is.na(SSL.fit$Offspring.Focal),],Focal.Sex=="F"&Chromosome!="T4"&Population!="T4"),family=binomial)#
Anova(netF.mortality,type="III")
netF.mortality=glmer(cbind((Competitor.Female.Death+Focal.Female.Death),12-(Competitor.Female.Death+Focal.Female.Death))~Chromosome*Treatment+ Block+(1|Population)+(1|Rep.ID), data=subset(SSL.fit[!is.na(SSL.fit$Offspring.Focal),],Focal.Sex=="F"&Chromosome!="T4"&Population!="T4"),family=binomial)#
Anova(netF.mortality,type="III")

#
#mate harm female mortality as quasi-binomial
#not binomial. Chr p=0.01369 and Treatment=0.063
MH.mortality=lmer(Lay.Ancestor.Death~Chromosome*Treatment+(1|Block)+(1|Population), data=subset(SSL.fit[!is.na(SSL.fit$Offspring.Focal),],Focal.Sex=="M"&Chromosome!="T4"))#
# 
#to avoid "must be between 0 and 1 error do COUNTS. Compared to 12 as total or remainder 12-deaths (dead, survived)
#comparing counts of dead vs alive
MH.mortality=glmer(cbind(Lay.Ancestor.Death,12-Lay.Ancestor.Death)~Chromosome*Treatment+ Block+(1|Population)+(1|Rep.ID), data=subset(SSL.fit[!is.na(SSL.fit$Offspring.Focal),],Focal.Sex=="M"&Chromosome!="T4"&Population!="T4"),family=binomial)#
Anova(MH.mortality,type="III")
summary(MH.mortality)

##Mate harm for fecundity model
#MH.fecundity=lmer(Offspring.Focal~Chromosome*Treatment+Lay.Ancestor.Death+(1|Block)+(1|Population), data=subset(SSL.fit[!is.na(SSL.fit$Offspring.Focal),],Focal.Sex=="M"&Chromosome!="ML"))#
MH.fecundity=lmer(Offspring.Total~Chromosome*Treatment+Lay.Ancestor.Death+(1|Block)+(1|Population), data=subset(SSL.fit[!is.na(SSL.fit$Offspring.Focal),],Focal.Sex=="M"&Chromosome!="T4"))#
#MH.fecundity=lm(Offspring.Total~Chromosome*Treatment+Lay.Ancestor.Death+Block+Treatment:Population, data=subset(SSL.fit[!is.na(SSL.fit$Offspring.Focal),],Focal.Sex=="M"&Chromosome!="T4"))#

MH.fecundity=lmer(Offspring.Total~Chromosome*Treatment+Lay.Ancestor.Death+Block+(1|Population:Treatment), data=subset(SSL.fit[!is.na(SSL.fit$Offspring.Focal),],Focal.Sex=="M"&Chromosome!="T4"))#
#Offspring.Total
MH.fecundity=glmer(Offspring.Total~Chromosome*Treatment+Lay.Ancestor.Death+Block+(1|Population), data=subset(SSL.fit[!is.na(SSL.fit$Offspring.Focal),],Focal.Sex=="M"&Chromosome!="T4"))#
Anova(MH.fecundity,type="III")
summary(MH.fecundity)

#get fecundity per individual female


#see if fecundity per individual female is still lower
MH.fecundity=lmer(Avg.Offspring.per.F~Chromosome*Treatment*Lay.Ancestor.Death+Block+(1|Population), data=subset(SSL.fit[!is.na(SSL.fit$Offspring.Focal),],Focal.Sex=="M"&Chromosome!="T4"))#

Anova(MH.fecundity,type="III")

Per.Line.Death=merge(SSL.fit[SSL.fit$Focal.Sex=="M",c("Line.ID","Population","Treatment","Chromosome","Lay.Ancestor.Death")],
               SSL.fit[SSL.fit$Focal.Sex=="F",c("Line.ID","Competitor.Female.Death")],by="Line.ID")

#Remove Block from the model as these were not done the same per line
MH.mortality=glmer(cbind(Lay.Ancestor.Death,12-Lay.Ancestor.Death)~Chromosome*Treatment+Competitor.Female.Death+(1|Population)+(1|Line.ID), data=subset(Per.Line.Death,Chromosome!="T4"&Population!="T4"),family=binomial)#
Anova(MH.mortality,type="III")

#examination of correlation of surviv

####==============================================------------------->RMF-----####

###rank by TreatmentSex (2x2)###
#pulling these out so that I can merge and compare their fitness later
F.fit=subset(SSL.fit,Focal.Sex=="F"&Chromosome!="T4")
F.fit$rankedFitness=NA
M.fit=subset(SSL.fit,Focal.Sex=="M"&Chromosome!="T4")
M.fit$rankedFitness=NA

#separating into groups so that I can rank them
HM.fit=subset(SSL.fit,Focal.Sex=="M"&Treatment=="HIGH"&Chromosome!="T4")
BM.fit=subset(SSL.fit,Focal.Sex=="M"&Treatment=="BENIGN"&Chromosome!="T4")

HF.fit=subset(SSL.fit,Focal.Sex=="F"&Treatment=="HIGH"&Chromosome!="T4")
BF.fit=subset(SSL.fit,Focal.Sex=="F"&Treatment=="BENIGN"&Chromosome!="T4")

#ordering the groups for ranking by Relative.Fit2 for B/H split
HM.fit=HM.fit[order(HM.fit$Relative.Fit2,decreasing=TRUE),]
BM.fit=BM.fit[order(BM.fit$Relative.Fit2,decreasing=TRUE),]
HF.fit=HF.fit[order(HF.fit$Relative.Fit2,decreasing=TRUE),]
BF.fit=BF.fit[order(BF.fit$Relative.Fit2,decreasing=TRUE),]

#For males go through the ranking and assign the ranked point
counter=0
for(line in HM.fit$Rep.ID){
  counter=counter+1
  M.fit[M.fit$Rep.ID==line,"rankedFitness"]=counter}
counter=0
for(line in BM.fit$Rep.ID){
  counter=counter+1
  M.fit[M.fit$Rep.ID==line,"rankedFitness"]=counter}
counter=0
#For females go through the ranking and assign the ranked point
for(line in HF.fit$Rep.ID){
  counter=counter+1
  F.fit[F.fit$Rep.ID==line,"rankedFitness"]=counter}
counter=0
for(line in BF.fit$Rep.ID){
  counter=counter+1
  F.fit[F.fit$Rep.ID==line,"rankedFitness"]=counter}

FM.ranked.fit=merge(F.fit[,c("Line.ID","Chromosome","Population","Treatment","Relative.Fit2","rankedFitness")],M.fit[,c("Line.ID","Relative.Fit2","rankedFitness")],by="Line.ID")
#col5 relative fit 2 (female), col6 rankedFit femalles
colnames(FM.ranked.fit)[5]="Relative.Fit.F"
colnames(FM.ranked.fit)[6]="rankedFitness.F"
colnames(FM.ranked.fit)[7]="Relative.Fit.M"
colnames(FM.ranked.fit)[8]="rankedFitness.M"

ggplot(data=subset(FM.ranked.fit),aes(y=rankedFitness.M,x=rankedFitness.F,colour=Chromosome))+
  geom_point()+theme_classic()+ggtitle("ranked rmf fit")+
  geom_smooth(fullrange=TRUE,span=2,aes(color=Chromosome,linetype=Treatment),na.rm=TRUE,se=FALSE,method="lm")

ggplot(data=subset(FM.ranked.fit,Treatment=="BENIGN"),aes(y=Relative.Fit.M,x=Relative.Fit.F,colour=Chromosome))+
  geom_point()+theme_classic()+ggtitle("relative rmf fit Benign")+
  geom_smooth(fullrange=TRUE,span=2,aes(color=Chromosome,linetype=Treatment),na.rm=TRUE,se=FALSE,method="lm")


###RANKED FITNESS
#get the rmf + significance
cor.test(FM.ranked.fit[FM.ranked.fit$Treatment=="BENIGN"&FM.ranked.fit$Chromosome=="ML","rankedFitness.M"],FM.ranked.fit[FM.ranked.fit$Treatment=="BENIGN"&FM.ranked.fit$Chromosome=="ML","rankedFitness.F"],method="spearman")
cor.test(FM.ranked.fit[FM.ranked.fit$Treatment=="HIGH"&FM.ranked.fit$Chromosome=="ML","rankedFitness.M"],FM.ranked.fit[FM.ranked.fit$Treatment=="HIGH"&FM.ranked.fit$Chromosome=="ML","rankedFitness.F"],method="spearman")

#FL
cor.test(FM.ranked.fit[FM.ranked.fit$Treatment=="BENIGN"&FM.ranked.fit$Chromosome=="FL","rankedFitness.M"],FM.ranked.fit[FM.ranked.fit$Treatment=="BENIGN"&FM.ranked.fit$Chromosome=="FL","rankedFitness.F"],method="spearman")
cor.test(FM.ranked.fit[FM.ranked.fit$Treatment=="HIGH"&FM.ranked.fit$Chromosome=="FL","rankedFitness.M"],FM.ranked.fit[FM.ranked.fit$Treatment=="HIGH"&FM.ranked.fit$Chromosome=="FL","rankedFitness.F"],method="spearman")

###RELATIVE FITNESS
#get the rmf + significance
#ALL
cor.test(FM.ranked.fit[FM.ranked.fit$Treatment=="BENIGN","Relative.Fit.M"],FM.ranked.fit[FM.ranked.fit$Treatment=="BENIGN","Relative.Fit.F"],method="spearman")
cor.test(FM.ranked.fit[FM.ranked.fit$Treatment=="HIGH","Relative.Fit.M"],FM.ranked.fit[FM.ranked.fit$Treatment=="HIGH","Relative.Fit.F"],method="spearman")

#FL
cor.test(FM.ranked.fit[FM.ranked.fit$Treatment=="BENIGN"&FM.ranked.fit$Chromosome=="FL","Relative.Fit.M"],FM.ranked.fit[FM.ranked.fit$Treatment=="BENIGN"&FM.ranked.fit$Chromosome=="FL","Relative.Fit.F"],method="spearman",exact=FALSE)
cor.test(FM.ranked.fit[FM.ranked.fit$Treatment=="HIGH"&FM.ranked.fit$Chromosome=="FL","Relative.Fit.M"],FM.ranked.fit[FM.ranked.fit$Treatment=="HIGH"&FM.ranked.fit$Chromosome=="FL","Relative.Fit.F"],method="spearman",exact=FALSE)
#ML
cor.test(FM.ranked.fit[FM.ranked.fit$Treatment=="BENIGN"&FM.ranked.fit$Chromosome=="ML","Relative.Fit.M"],FM.ranked.fit[FM.ranked.fit$Treatment=="BENIGN"&FM.ranked.fit$Chromosome=="ML","Relative.Fit.F"],method="spearman")
cor.test(FM.ranked.fit[FM.ranked.fit$Treatment=="HIGH"&FM.ranked.fit$Chromosome=="ML","Relative.Fit.M"],FM.ranked.fit[FM.ranked.fit$Treatment=="HIGH"&FM.ranked.fit$Chromosome=="ML","Relative.Fit.F"],method="spearman")

# ####===========================--------center by by CHRxTreatmentxSex (2x2x2)-----####
# 
# 
# ###---RANK BY MEAN SUBTRACTION THAN SPEARMAN RANK
# 
# #have to create the new column before calling the function
# SSL.fit$meanGroupNorm.Fit=0
# #the function won't modify the original data.frame so save new
# SSL.fit=GroupNorm_func(SSL.fit,"meanGroupNorm.Fit","Relative.Fit2",c("F","M"))
# 
# 
# #splitting them to combine them per line
# F.fit=subset(SSL.fit,Focal.Sex=="F"&Chromosome!="T4")
# M.fit=subset(SSL.fit,Focal.Sex=="M"&Chromosome!="T4")
# #merging the fitness to organize it per line
# Per.Line.Fit=merge(F.fit[,c("Line.ID","Chromosome","Population","Treatment","Relative.Fit2","meanGroupNorm.Fit")],M.fit[,c("Line.ID","Relative.Fit2","meanGroupNorm.Fit","Fake.Total.Offspring")],by="Line.ID")
# #renaming the colums to distinguish Female and Male fitness
# colnames(Per.Line.Fit)[5]="Relative.Fit.F"
# colnames(Per.Line.Fit)[6]="GroupNorm.Fit.F"
# 
# colnames(Per.Line.Fit)[7]="Relative.Fit.M"
# colnames(Per.Line.Fit)[8]="GroupNorm.Fit.M"
# 
# #basic visual inspection
# ggplot(data=subset(Per.Line.Fit),aes(y=GroupNorm.Fit.M,x=GroupNorm.Fit.F,colour=Chromosome,shape=Treatment))+
#   geom_point()+theme_classic()+ggtitle("meanNormalized rMF fitness")+ylim(-0.5,0.5)+xlim(-0.5,0.5)+
#   #scale_shape_manual(values=c(0,0,16,16))+
#   geom_smooth(fullrange=TRUE,span=2,aes(color=Chromosome,linetype=Treatment),na.rm=TRUE,se=FALSE,method="lm")
# 
# #examine if increase mortality confers success
# ggplot(data=subset(Per.Line.Fit),aes(y=Competitor.Female.Death,x=Lay.Ancestor.Death,colour=Chromosome,shape=Treatment))+
#   geom_point()+theme_classic()+ggtitle("Line Death Correlations")+#ylim(-0.5,0.5)+#xlim(-0.5,0.5)+
#   #scale_shape_manual(values=c(0,0,16,16))+
#   geom_smooth(fullrange=TRUE,span=2,aes(color=Chromosome,linetype=Treatment),na.rm=TRUE,se=FALSE,method="lm")
# 
# 
# #
# 
# 
# ##Do actual correlation testing split by FL and ML
# #FL
# cor.test(Per.Line.Fit[Per.Line.Fit$Treatment=="HIGH"&Per.Line.Fit$Chromosome=="FL","GroupNorm.Fit.M"],
#          Per.Line.Fit[Per.Line.Fit$Treatment=="HIGH"&Per.Line.Fit$Chromosome=="FL","GroupNorm.Fit.F"],method="spearman")
# cor.test(Per.Line.Fit[Per.Line.Fit$Treatment=="BENIGN"&Per.Line.Fit$Chromosome=="FL","GroupNorm.Fit.M"],
#          Per.Line.Fit[Per.Line.Fit$Treatment=="BENIGN"&Per.Line.Fit$Chromosome=="FL","GroupNorm.Fit.F"],method="spearman")
#          
# #ML
# cor.test(Per.Line.Fit[Per.Line.Fit$Treatment=="HIGH"&Per.Line.Fit$Chromosome=="ML","GroupNorm.Fit.M"],
#          Per.Line.Fit[Per.Line.Fit$Treatment=="HIGH"&Per.Line.Fit$Chromosome=="ML","GroupNorm.Fit.F"],method="spearman")
# cor.test(Per.Line.Fit[Per.Line.Fit$Treatment=="BENIGN"&Per.Line.Fit$Chromosome=="ML","GroupNorm.Fit.M"],
#          Per.Line.Fit[Per.Line.Fit$Treatment=="BENIGN"&Per.Line.Fit$Chromosome=="ML","GroupNorm.Fit.F"],method="spearman")
# 
# ####DOING THE ABOVE PER POPULATION####
# ####MEANS 2TREATx2SEXx2CHRx3REP POPS
# #separating into groups to get mean. (Not re-calculating realtive fitness for each bootstrap though)
# #I *could*?
# #ML
# #Male
# group.means=data_summary(SSL.fit,"Relative.Fit2",c("Chromosome","Treatment","Focal.Sex","Population"))
# 
# groupNorm_func <-function(data, chr, sex, means){
#   Group.Norm.Fitness=subset(data,Focal.Sex==sex&Chromosome==chr)
#   #the normalized fitness column
#   Group.Norm.Fitness$meanGroupNorm.Fit=NA
#   for (pop in means$Population){
#     Group.Norm.Fitness[Group.Norm.Fitness$Population==pop,"meanGroupNorm.Fit"]=Group.Norm.Fitness[Group.Norm.Fitness$Population==pop,"Relative.Fit2"]-
#       means[means$Focal.Sex==sex&means$Population==pop&means$Chromosome==chr,"Relative.Fit2"]
#   }
#   return(Group.Norm.Fitness)  
# }
# #do fitness fo each sex separately to make merging easier
# M.Fit=rbind(groupNorm_func(SSL.fit,"ML","M",group.means),groupNorm_func(SSL.fit,"FL","M",group.means))
# F.Fit=rbind(groupNorm_func(SSL.fit,"ML","F",group.means),groupNorm_func(SSL.fit,"FL","F",group.means))
# #merging the fitness to organize it per line
# Per.Line.Fit=merge(F.Fit[,c("Line.ID","Chromosome","Population","Treatment","Relative.Fit2","meanGroupNorm.Fit")],M.Fit[,c("Line.ID","Fit2","meanGroupNorm.Fit","Offspring.Focal")],by="Line.ID")
# 
# #renaming the colums to distinguish Female and Male fitness
# colnames(Per.Line.Fit)[5]="Fit.F"
# colnames(Per.Line.Fit)[6]="GroupNorm.Fit.F"
# colnames(Per.Line.Fit)[7]="Relative.Fit.M"
# colnames(Per.Line.Fit)[8]="Fit.M"
# 
# #basic visual inspection
# ggplot(data=subset(Per.Line.Fit),aes(y=GroupNorm.Fit.M,x=GroupNorm.Fit.F,colour=Chromosome,shape=Treatment))+
#   geom_point()+theme_classic()+ggtitle("meanNormalized rMF fitness (sep by pop)")+ylim(-0.5,0.5)+xlim(-0.5,0.5)+
#   #scale_shape_manual(values=c(0,0,16,16))+
#   geom_smooth(fullrange=TRUE,span=2,aes(color=Chromosome,linetype=Treatment),na.rm=TRUE,se=FALSE,method="lm")
# 
# ggplot(data=subset(Per.Line.Fit,Population=="H1"),aes(y=GroupNorm.Fit.M,x=GroupNorm.Fit.F,colour=Chromosome))+
#   geom_point()+theme_classic()+ggtitle("meanNormalized rMF fitness (H1)")+#ylim(-0.5,0.5)+xlim(-0.5,0.5)+
#   geom_path(aes(x=GroupNorm.Fit.F, y=GroupNorm.Fit.M,colour=Chromosome), size=1, linetype=2)
#   #scale_shape_manual(values=c(0,0,16,16))+
#   geom_smooth(fullrange=TRUE,span=2,aes(color=Chromosome,linetype=Treatment),na.rm=TRUE,se=FALSE,method="lm")
# 
# #qplot(data=df, x=x,y=y,colour=group)+stat_ellipse(level=0.95)+scale_colour_manual(values=colorpal)
# qplot(data=subset(Per.Line.Fit,Population=="H1"),y=GroupNorm.Fit.M,x=GroupNorm.Fit.F,colour=Chromosome)+
#   ggtitle("meanNormalized rMF fitness (H1)")+stat_ellipse(level=0.95)+theme_classic2()
# qplot(data=subset(Per.Line.Fit,Population=="H2"),y=GroupNorm.Fit.M,x=GroupNorm.Fit.F,colour=Chromosome)+
#   ggtitle("meanNormalized rMF fitness (H2)")+stat_ellipse(level=0.95)+theme_classic2()
# qplot(data=subset(Per.Line.Fit,Population=="H3"),y=GroupNorm.Fit.M,x=GroupNorm.Fit.F,colour=Chromosome)+
#   ggtitle("meanNormalized rMF fitness (H3)")+stat_ellipse(level=0.95)+theme_classic2()
# qplot(data=subset(Per.Line.Fit,Population=="B1"),y=GroupNorm.Fit.M,x=GroupNorm.Fit.F,colour=Chromosome)+
#   ggtitle("meanNormalized rMF fitness (B1)")+stat_ellipse(level=0.95)+theme_classic2()
# qplot(data=subset(Per.Line.Fit,Population=="B2"),y=GroupNorm.Fit.M,x=GroupNorm.Fit.F,colour=Chromosome)+
#   ggtitle("meanNormalized rMF fitness (B2)")+stat_ellipse(level=0.95)+theme_classic2()
# qplot(data=subset(Per.Line.Fit,Population=="B3"),y=GroupNorm.Fit.M,x=GroupNorm.Fit.F,colour=Chromosome)+
#   ggtitle("meanNormalized rMF fitness (B3)")+stat_ellipse(level=0.95)+theme_classic2()
# 
#   ##Do actual correlation testing split by FL and ML
# #FL
# cor.test(Per.Line.Fit[Per.Line.Fit$Treatment=="HIGH"&Per.Line.Fit$Chromosome=="FL","GroupNorm.Fit.M"],
#          Per.Line.Fit[Per.Line.Fit$Treatment=="HIGH"&Per.Line.Fit$Chromosome=="FL","GroupNorm.Fit.F"],method="spearman")
# cor.test(Per.Line.Fit[Per.Line.Fit$Treatment=="BENIGN"&Per.Line.Fit$Chromosome=="FL","GroupNorm.Fit.M"],
#          Per.Line.Fit[Per.Line.Fit$Treatment=="BENIGN"&Per.Line.Fit$Chromosome=="FL","GroupNorm.Fit.F"],method="spearman")
# 
# #ML
# cor.test(Per.Line.Fit[Per.Line.Fit$Treatment=="HIGH"&Per.Line.Fit$Chromosome=="ML","GroupNorm.Fit.M"],
#          Per.Line.Fit[Per.Line.Fit$Treatment=="HIGH"&Per.Line.Fit$Chromosome=="ML","GroupNorm.Fit.F"],method="spearman")
# cor.test(Per.Line.Fit[Per.Line.Fit$Treatment=="BENIGN"&Per.Line.Fit$Chromosome=="ML","GroupNorm.Fit.M"],
#          Per.Line.Fit[Per.Line.Fit$Treatment=="BENIGN"&Per.Line.Fit$Chromosome=="ML","GroupNorm.Fit.F"],method="spearman")
# 
# #ALL (T4 will have NA value for GroupNorm.Fit.M)
# cor.test(Per.Line.Fit[Per.Line.Fit$Treatment=="HIGH","GroupNorm.Fit.M"],
#          Per.Line.Fit[Per.Line.Fit$Treatment=="HIGH","GroupNorm.Fit.F"],method="spearman")
# cor.test(Per.Line.Fit[Per.Line.Fit$Treatment=="BENIGN","GroupNorm.Fit.M"],
#          Per.Line.Fit[Per.Line.Fit$Treatment=="BENIGN","GroupNorm.Fit.F"],method="spearman")
# 
# #Examine non-ranked for fun (don't use)
# cor.test(Per.Line.Fit[Per.Line.Fit$Treatment=="HIGH"&Per.Line.Fit$Chromosome=="FL","GroupNorm.Fit.M"],
#          Per.Line.Fit[Per.Line.Fit$Treatment=="HIGH"&Per.Line.Fit$Chromosome=="FL","GroupNorm.Fit.F"])
# cor.test(Per.Line.Fit[Per.Line.Fit$Treatment=="BENIGN"&Per.Line.Fit$Chromosome=="FL","GroupNorm.Fit.M"],
#          Per.Line.Fit[Per.Line.Fit$Treatment=="BENIGN"&Per.Line.Fit$Chromosome=="FL","GroupNorm.Fit.F"])
# 
# #ML
# cor.test(Per.Line.Fit[Per.Line.Fit$Treatment=="HIGH"&Per.Line.Fit$Chromosome=="ML","GroupNorm.Fit.M"],
#          Per.Line.Fit[Per.Line.Fit$Treatment=="HIGH"&Per.Line.Fit$Chromosome=="ML","GroupNorm.Fit.F"])
# cor.test(Per.Line.Fit[Per.Line.Fit$Treatment=="BENIGN"&Per.Line.Fit$Chromosome=="ML","GroupNorm.Fit.M"],
#          Per.Line.Fit[Per.Line.Fit$Treatment=="BENIGN"&Per.Line.Fit$Chromosome=="ML","GroupNorm.Fit.F"])
# 

####------do rmf correlation test by Bootstrapping---####

#figure out number of lines in each treatment
# reps.B1FL=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="B1","Line.ID"])#32
# reps.B2FL=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="B2","Line.ID"])#30
# reps.B3FL=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="B3","Line.ID"])#29
# reps.B1ML=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="B1","Line.ID"])#27
# reps.B2ML=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="B2","Line.ID"])#27
# reps.B3ML=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="B3","Line.ID"])#33

####changing this section

####here getting groupnormalized fitness
#splitting them to combine them per line
SSL.fit$meanGroupNorm.Fit=0
SSL.fit=GroupNorm_func(SSL.fit,"meanGroupNorm.Fit","Relative.Fit2",c("M","F"))
F.fit=subset(SSL.fit,Focal.Sex=="F"&Chromosome!="T4")
M.fit=subset(SSL.fit,Focal.Sex=="M"&Chromosome!="T4")
#merging the fitness to organize it per line
Per.Line.Fit=merge(F.fit[,c("Line.ID","Chromosome","Population","Treatment","Relative.Fit2","meanGroupNorm.Fit")],M.fit[,c("Line.ID","Relative.Fit2","meanGroupNorm.Fit")],by="Line.ID")
#renaming the colums to distinguish Female and Male fitness
colnames(Per.Line.Fit)[5]="Relative.Fit.F"
colnames(Per.Line.Fit)[6]="GroupNorm.Fit.F"

colnames(Per.Line.Fit)[7]="Relative.Fit.M"
colnames(Per.Line.Fit)[8]="GroupNorm.Fit.M"

#Pearson now
#getting a boostraped rho/corr of the rmf
BFL.boot=replicate(1000,{
  #pull out replicates per population. The default size is n. So it might be uncessary to calculate length?
  samp.B1=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="B1","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="B1","Line.ID"]),replace=TRUE)#32
  samp.B2=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="B2","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="B2","Line.ID"]),replace=TRUE)#30
  samp.B3=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="B3","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="B3","Line.ID"]),replace=TRUE)#29
  
  samp.BFL=append(samp.B1,samp.B2)
  samp.BFL=append(samp.BFL,samp.B3)
  
 BFL.rmf=cor.test(Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.BFL,"GroupNorm.Fit.M"],Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.BFL,"GroupNorm.Fit.F"])
 #rho is value four
 BFL.rmf[[4]]#this is just value
})

BML.boot=replicate(1000,{
  #pull out replicates (all)
  
  samp.B1=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="B1","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="B1","Line.ID"]),replace=TRUE)#32
  samp.B2=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="B2","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="B2","Line.ID"]),replace=TRUE)#30
  samp.B3=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="B3","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="B3","Line.ID"]),replace=TRUE)#29
  
  samp.BML=append(samp.B1,samp.B2)
  samp.BML=append(samp.BML,samp.B3)
  
  #rho is value four
  BML.rmf=cor.test(Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.BML,"GroupNorm.Fit.M"],Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.BML,"GroupNorm.Fit.F"])
  #rho is value four
  #BFL.rmf[4]#this will keep that it's rho
  BML.rmf[[4]]#this is just value
})

HFL.boot=replicate(1000,{
  #pull out replicates (all)
  
  samp.H1=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="H1","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="H1","Line.ID"]),replace=TRUE)#32
  samp.H2=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="H2","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="H2","Line.ID"]),replace=TRUE)#30
  samp.H3=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="H3","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="FL"&Per.Line.Fit$Population=="H3","Line.ID"]),replace=TRUE)#29
  
  samp.HFL=append(samp.H1,samp.H2)
  samp.HFL=append(samp.HFL,samp.H3)
  
  #rho is value four
  HFL.rmf=cor.test(Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.HFL,"GroupNorm.Fit.M"],Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.HFL,"GroupNorm.Fit.F"])
  #rho is value four
  #BFL.rmf[4]#this will keep that it's rho
  HFL.rmf[[4]]#this is just value
})

HML.boot=replicate(1000,{
  #pull out replicates (all)
  
  samp.H1=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="H1","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="H1","Line.ID"]),replace=TRUE)#32
  samp.H2=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="H2","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="H2","Line.ID"]),replace=TRUE)#30
  samp.H3=sample(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="H3","Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome=="ML"&Per.Line.Fit$Population=="H3","Line.ID"]),replace=TRUE)#29
  
  samp.HML=append(samp.H1,samp.H2)
  samp.HML=append(samp.HML,samp.H3)
  
  #rho is value four
  HML.rmf=cor.test(Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.HML,"GroupNorm.Fit.M"],Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.HML,"GroupNorm.Fit.F"])
  #rho is value four
  #BFL.rmf[4]#this will keep that it's rho
  HML.rmf[[4]]#this is just value
})


##new and easier way
boot_corr<-function(Per.Line.Fit,chr,pop,col1,col2){
  BFL.Adboot=replicate(1000,{
    #pull out replicates per population. The default size is n. So it might be uncessary to calculate length?
    samp.B1=sample(Per.Line.Fit[Per.Line.Fit$Chromosome==chr&Per.Line.Fit$Population==pop[1],"Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome==chr&Per.Line.Fit$Population==pop[1],"Line.ID"]),replace=TRUE)#32
    samp.B2=sample(Per.Line.Fit[Per.Line.Fit$Chromosome==chr&Per.Line.Fit$Population==pop[2],"Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome==chr&Per.Line.Fit$Population==pop[2],"Line.ID"]),replace=TRUE)#30
    samp.B3=sample(Per.Line.Fit[Per.Line.Fit$Chromosome==chr&Per.Line.Fit$Population==pop[3],"Line.ID"],size=length(Per.Line.Fit[Per.Line.Fit$Chromosome==chr&Per.Line.Fit$Population==pop[3],"Line.ID"]),replace=TRUE)#29
    
    samp.BFL=append(samp.B1,samp.B2)
    samp.BFL=append(samp.BFL,samp.B3)
    
    #rho is value four
    BFL.rmf=cor.test(Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.BFL,col1],Per.Line.Fit[Per.Line.Fit$Line.ID%in%samp.BFL,col2])
    #rho is value four
    #BFL.rmf[4]#this will keep that it's rho
    BFL.rmf[[4]]#this is just value
  })
  
  
  return(BFL.Adboot)
}
diff.way=TRUE
if(diff.way){
  
  SSL.fit$GN.Fit2=0
  SSL.fit=GroupNorm_func(SSL.fit,"GN.Fit2","Relative.Fit2",c("M","F"))
  M.Fit=subset(SSL.fit,Focal.Sex=="M")
  F.Fit=subset(SSL.fit,Focal.Sex=="F")
  
  Per.Line.Fit2=merge(F.Fit[,c("Line.ID","Chromosome","Population","Treatment","Fit2","GN.Fit2")],M.Fit[,c("Line.ID","Fit2","GN.Fit2")],by="Line.ID")
  #renaming the colums to distinguish Female and Male fitness
  colnames(Per.Line.Fit2)[5]="Fit2.F"
  colnames(Per.Line.Fit2)[6]="GroupNorm.Fit2.F"
  colnames(Per.Line.Fit2)[7]="Fit2.M"
  colnames(Per.Line.Fit2)[8]="GroupNorm.Fit2.M"
  
  BFL.boot=boot_corr(Per.Line.Fit2,"FL",c("B1","B2","B3"),"GroupNorm.Fit2.F","GroupNorm.Fit2.M")
  BML.boot=boot_corr(Per.Line.Fit2,"ML",c("B1","B2","B3"),"GroupNorm.Fit2.F","GroupNorm.Fit2.M")
  HFL.boot=boot_corr(Per.Line.Fit2,"FL",c("H1","H2","H3"),"GroupNorm.Fit2.F","GroupNorm.Fit2.M")
  HML.boot=boot_corr(Per.Line.Fit2,"ML",c("H1","H2","H3"),"GroupNorm.Fit2.F","GroupNorm.Fit2.M")
  
  
}


#compare the distributions
#can do a corr of corrs? or cor1-cor2?
##start with ML. ##does the pairing matter?????
##if the corrs are the same the difference will be 0. If there is no 0 diff, then SIG DIFF
diff.ML=HML.boot-BML.boot
diff.ML.CI=quantile(diff.ML,c(0.975,0.025))# interval(0.5441168 <-> 0.1124657) SIG
diff.FL=HFL.boot-BFL.boot
diff.FL.CI=quantile(diff.FL,c(0.975,0.025))#interval (0.414297915 <-> -0.003945277) NOT SIG

#now get the average per treatment
Ball.boot=(BFL.boot+BML.boot)/2
Hall.boot=(HFL.boot+HML.boot)/2

diff.all=Hall.boot-Ball.boot
diff.all.CI=quantile(diff.all,c(0.975,0.025))

#---GRAPH RMF-----"Fit2"---(use)#####
plotting_values=data.frame(mean=c(mean(BFL.boot),mean(BML.boot),mean(HFL.boot),mean(HML.boot),mean(Ball.boot),mean(Hall.boot)),Chromosome=c("BenignFL","BenignML","NovelFL","NovelML","B.avg","N.avg"),
                           upper_CI=c(quantile(BFL.boot,0.975),quantile(BML.boot,0.975),quantile(HFL.boot,0.975),quantile(HML.boot,0.975),quantile(Ball.boot,0.975),quantile(Hall.boot,0.975)),
                           lower_CI=c(quantile(BFL.boot,0.025),quantile(BML.boot,0.025),quantile(HFL.boot,0.025),quantile(HML.boot,0.025),quantile(Ball.boot,0.025),quantile(Hall.boot,0.025)))
#this is to group generally if needed
plotting_values$Chr=c("FL","ML","FL","ML","Average","Average")
plotting_values$Treatment=c("Benign","Benign","Novel","Novel","Benign","Novel")
Fig2=ggplot(plotting_values,aes(x=Chr,y=as.numeric(mean),fill=Chromosome,colour=Chromosome))+
  geom_point(size=3,alpha=.9,shape=21)+ 
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI,colour=Chromosome), width=.7,position=position_dodge(.9))+
  scale_fill_manual(values=c("black","#E41A1C","#3e86ff","black","#E41A1C","#3e86ff"))+
  scale_colour_manual(values=c("black","#E41A1C","#3e86ff","black","#E41A1C","#3e86ff"))+#"#f6564b
  geom_hline(yintercept=0,colour="grey",linewidth=0.5,linetype="dashed")+
  facet_wrap(~Treatment,scales="free_x")+
  #ggtitle("groupNorm Relative.Fit2 rmf (corr) w/ \nstructure CI Pearson")+
  ylab("Intersexual Correlation for Fitness")+xlab("Chromosome Pool")+
  theme_classic() +theme(text = element_text(size = 14,color="black"),legend.position="none",panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7))+
  theme(axis.text.x = element_text(size = 12,color="black"), axis.text.y = element_text(size = 12,color="black"))






#####Overall Fecundity per Individual and Surviving Female####

#first just total death and per F fecundity
F.fit=subset(SSL.fit,Chromosome!="T4")
death.GroupNorm=data_summary(F.fit,"Surviving.Females",c("Population","Chromosome","Treatment"))
fec.GroupNorm=data_summary(F.fit,"Avg.Offspring.per.F",c("Population","Chromosome","Treatment"))
populations=unique(death.GroupNorm$Population)
F.fit$GroupNorm.Death=0
F.fit$GroupNorm.Avg.F=0
counter=0

for(pop in populations){
  F.fit[F.fit$Population==pop&F.fit$Chromosome=="ML","GroupNorm.Death"]=
    F.fit[F.fit$Population==pop&F.fit$Chromosome=="ML","Surviving.Females"]-death.GroupNorm$Surviving.Females[[counter+1]]
  F.fit[F.fit$Population==pop&F.fit$Chromosome=="FL","GroupNorm.Death"]=
    F.fit[F.fit$Population==pop&F.fit$Chromosome=="FL","Surviving.Females"]-death.GroupNorm$Surviving.Females[[counter+2]]
  print(pop)
  counter=counter+2
}
counter=0
for(pop in populations){
  F.fit[F.fit$Population==pop&F.fit$Chromosome=="ML","GroupNorm.Avg.F"]=
    F.fit[F.fit$Population==pop&F.fit$Chromosome=="ML","Avg.Offspring.per.F"]-fec.GroupNorm$Avg.Offspring.per.F[[counter+1]]
  F.fit[F.fit$Population==pop&F.fit$Chromosome=="FL","GroupNorm.Avg.F"]=
    F.fit[F.fit$Population==pop&F.fit$Chromosome=="FL","Avg.Offspring.per.F"]-fec.GroupNorm$Avg.Offspring.per.F[[counter+2]]
  counter=counter+2
  print(pop)
}





BFL.fecoot=replicate(1000,{
  #pull out replicates per population. The default size is n. So it might be uncessary to calculate length?
  samp.B1=sample(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="B1","Line.ID"],size=length(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="B1","Line.ID"]),replace=TRUE)#32
  samp.B2=sample(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="B2","Line.ID"],size=length(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="B2","Line.ID"]),replace=TRUE)#30
  samp.B3=sample(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="B3","Line.ID"],size=length(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="B3","Line.ID"]),replace=TRUE)#29
  
  samp.BFL=append(samp.B1,samp.B2)
  samp.BFL=append(samp.BFL,samp.B3)
  
  #rho is value four
  BFL.rmf=cor.test(F.fit[F.fit$Line.ID%in%samp.BFL,"GroupNorm.Death"],F.fit[F.fit$Line.ID%in%samp.BFL,"GroupNorm.Avg.F"])
  #rho is value four
  #BFL.rmf[4]#this will keep that it's rho
  BFL.rmf[[4]]#this is just value
})
BML.fecoot=replicate(1000,{
  #pull out replicates per population. The default size is n. So it might be uncessary to calculate length?
  samp.B1=sample(F.fit[F.fit$Chromosome=="ML"&F.fit$Population=="B1","Line.ID"],size=length(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="B1","Line.ID"]),replace=TRUE)#32
  samp.B2=sample(F.fit[F.fit$Chromosome=="ML"&F.fit$Population=="B2","Line.ID"],size=length(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="B2","Line.ID"]),replace=TRUE)#30
  samp.B3=sample(F.fit[F.fit$Chromosome=="ML"&F.fit$Population=="B3","Line.ID"],size=length(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="B3","Line.ID"]),replace=TRUE)#29
  
  samp.BML=append(samp.B1,samp.B2)
  samp.BML=append(samp.BML,samp.B3)
  
  #rho is value four
  BML.rmf=cor.test(F.fit[F.fit$Line.ID%in%samp.BML,"GroupNorm.Death"],F.fit[F.fit$Line.ID%in%samp.BML,"GroupNorm.Avg.F"])
  #rho is value four
  #BFL.rmf[4]#this will keep that it's rho
  BML.rmf[[4]]#this is just value
})

HFL.fecoot=replicate(1000,{
  #pull out replicates per population. The default size is n. So it might be uncessary to calculate length?
  samp.B1=sample(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="H1","Line.ID"],size=length(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="H1","Line.ID"]),replace=TRUE)#32
  samp.B2=sample(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="H2","Line.ID"],size=length(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="H2","Line.ID"]),replace=TRUE)#30
  samp.B3=sample(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="H3","Line.ID"],size=length(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="H3","Line.ID"]),replace=TRUE)#29
  
  samp.BFL=append(samp.B1,samp.B2)
  samp.BFL=append(samp.BFL,samp.B3)
  
  #rho is value four
  BFL.rmf=cor.test(F.fit[F.fit$Line.ID%in%samp.BFL,"GroupNorm.Death"],F.fit[F.fit$Line.ID%in%samp.BFL,"GroupNorm.Avg.F"])
  #rho is value four
  #BFL.rmf[4]#this will keep that it's rho
  BFL.rmf[[4]]#this is just value
})
#got lazy and didn't probably change all variable names
HML.fecoot=replicate(1000,{
  #pull out replicates per population. The default size is n. So it might be uncessary to calculate length?
  samp.B1=sample(F.fit[F.fit$Chromosome=="ML"&F.fit$Population=="H1","Line.ID"],size=length(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="H1","Line.ID"]),replace=TRUE)#32
  samp.B2=sample(F.fit[F.fit$Chromosome=="ML"&F.fit$Population=="H2","Line.ID"],size=length(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="H2","Line.ID"]),replace=TRUE)#30
  samp.B3=sample(F.fit[F.fit$Chromosome=="ML"&F.fit$Population=="H3","Line.ID"],size=length(F.fit[F.fit$Chromosome=="FL"&F.fit$Population=="H3","Line.ID"]),replace=TRUE)#29
  
  samp.BML=append(samp.B1,samp.B2)
  samp.BML=append(samp.BML,samp.B3)
  
  #rho is value four
  BML.rmf=cor.test(F.fit[F.fit$Line.ID%in%samp.BML,"GroupNorm.Death"],F.fit[F.fit$Line.ID%in%samp.BML,"GroupNorm.Avg.F"])
  #rho is value four
  #BFL.rmf[4]#this will keep that it's rho
  BML.rmf[[4]]#this is just value
})

diff.ML=HML.fecoot-BML.fecoot
diff.ML.CI=quantile(diff.ML,c(0.975,0.025))# interval(0.5441168 <-> 0.1124657) SIG
diff.FL=HFL.fecoot-BFL.fecoot
diff.FL.CI=quantile(diff.FL,c(0.975,0.025))#interval (0.414297915 <-> -0.003945277) NOT SIG

#now get the average per treatment
Ball.fecoot=(BFL.fecoot+BML.fecoot)/2
Hall.fecoot=(HFL.fecoot+HML.fecoot)/2

diff.all=Hall.fecoot-Ball.fecoot
diff.all.CI=quantile(diff.all,c(0.975,0.025))

plotting_values=data.frame(mean=c(mean(BFL.fecoot),mean(BML.fecoot,na.rm=TRUE),mean(HFL.fecoot),mean(HML.fecoot),mean(Hall.fecoot),mean(Ball.fecoot)),
                           Chromosome=c("BFL","BML","HFL","HML","Havg","Bavg"),
                           upper_CI=c(quantile(BFL.fecoot,0.975),quantile(BML.fecoot,0.975,na.rm=TRUE),quantile(HFL.fecoot,0.975),quantile(HML.fecoot,0.975),quantile(Hall.fecoot,0.975),quantile(Ball.fecoot,0.975)),
                           lower_CI=c(quantile(BFL.fecoot,0.025),quantile(BML.fecoot,0.025,na.rm=TRUE),quantile(HFL.fecoot,0.025),quantile(HML.fecoot,0.025),quantile(Hall.fecoot,0.025),quantile(Ball.fecoot,0.025)),
                           median=c(quantile(BFL.fecoot,0.5),quantile(BML.fecoot,0.5,na.rm=TRUE),quantile(HFL.fecoot,0.5),quantile(HML.fecoot,0.5),quantile(Hall.fecoot,0.5),quantile(Ball.fecoot,0.5)))
#overall=c(BFL.corr[[4]],BML.corr[[4]],HML.corr[[4]],HFL.corr[[4]]))
ggplot(plotting_values,aes(x=Chromosome,y=median))+geom_point()+ggtitle("per line corr (GroupNorm) Total Female Survival (F.fit) \n Fecundity (perF)")+
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), width=.2,position=position_dodge(.9), colour="black")+theme_classic()





####==========****************--------------------------->Local Adaptation <----***********####


#take the centered per line fitness
SSL.fit$meanGroupNorm.Fit=0
SSL.fit=GroupNorm_func(SSL.fit,"meanGroupNorm.Fit","Relative.Fit2",c("M","F"))

#splitting them to combine them per line
F.fit=subset(SSL.fit,Focal.Sex=="F"&Chromosome!="T4")
M.fit=subset(SSL.fit,Focal.Sex=="M"&Chromosome!="T4")
#merging the fitness to organize it per line
Per.Line.Fit=merge(F.fit[,c("Line.ID","Chromosome","Population","Treatment","Relative.Fit2","meanGroupNorm.Fit")],
                   M.fit[,c("Line.ID","Relative.Fit2","meanGroupNorm.Fit","Lay.Ancestor.Death","Offspring.Total","Avg.Offspring.per.F")],by="Line.ID")
#renaming the colums to distinguish Female and Male fitness
colnames(Per.Line.Fit)[5]="Relative.Fit.F"
colnames(Per.Line.Fit)[6]="GroupNorm.Fit.F"

colnames(Per.Line.Fit)[7]="Relative.Fit.M"
colnames(Per.Line.Fit)[8]="GroupNorm.Fit.M"


####not maintaining structure---->because using the centered values
#Calculate "Adaptation" improvement of ML chromosomes in males
boot_diff.max<-function(data,chr,pop,col){
  #this substracts two bootstrapped means. (Note the first chr always subtracted from second ignore FL/ML notation)
  #divides by maximum
  Benefit.boot=replicate(1000,{
    #pull out replicates per population. The default size is n. So it might be uncessary to calculate length?
    samp.B1=sample(data[data$Chromosome==chr[1]&data$Population==pop[1],"Line.ID"],size=length(data[data$Chromosome==chr[1]&data$Population==pop[1],"Line.ID"]),replace=TRUE)#32
    samp.B2=sample(data[data$Chromosome==chr[1]&data$Population==pop[2],"Line.ID"],size=length(data[data$Chromosome==chr[1]&data$Population==pop[2],"Line.ID"]),replace=TRUE)#30
    samp.B3=sample(data[data$Chromosome==chr[1]&data$Population==pop[3],"Line.ID"],size=length(data[data$Chromosome==chr[1]&data$Population==pop[3],"Line.ID"]),replace=TRUE)#29
    samp.BFL=append(samp.B1,samp.B2)
    samp.BFL=append(samp.BFL,samp.B3)
    chr1.mean=mean(data[data$Line.ID%in%samp.BFL,col])
 
    #pull out replicates per population. The default size is n. So it might be uncessary to calculate length?
    samp.B1a=sample(data[data$Chromosome==chr[2]&data$Population==pop[1],"Line.ID"],size=length(data[data$Chromosome==chr[2]&data$Population==pop[1],"Line.ID"]),replace=TRUE)#32
    samp.B2a=sample(data[data$Chromosome==chr[2]&data$Population==pop[2],"Line.ID"],size=length(data[data$Chromosome==chr[2]&data$Population==pop[2],"Line.ID"]),replace=TRUE)#30
    samp.B3a=sample(data[data$Chromosome==chr[2]&data$Population==pop[3],"Line.ID"],size=length(data[data$Chromosome==chr[2]&data$Population==pop[3],"Line.ID"]),replace=TRUE)#29
    samp.BML=append(samp.B1a,samp.B2a)
    samp.BML=append(samp.BML,samp.B3a)
    chr2.mean=mean(data[data$Line.ID%in%samp.BML,col])
    
    Benefit=chr1.mean-chr2.mean/(max(chr1.mean,chr2.mean))
    Benefit
  } )
  return(Benefit.boot)
}

LocalAd.BM=boot_diff.max(Per.Line.Fit,c("ML","FL"),c("B1","B2","B3"),"GroupNorm.Fit.M")
LocalAd.HM=boot_diff.max(Per.Line.Fit,c("ML","FL"),c("H1","H2","H3"),"GroupNorm.Fit.M")

LocalAd.BF=boot_diff.max(Per.Line.Fit,c("FL","ML"),c("B1","B2","B3"),"GroupNorm.Fit.F")
LocalAd.HF=boot_diff.max(Per.Line.Fit,c("FL","ML"),c("H1","H2","H3"),"GroupNorm.Fit.F")

LocalAd.M=(LocalAd.BM+LocalAd.HM)/2
LocalAd.F=(LocalAd.BF+LocalAd.HF)/2
LocalAd.B=(LocalAd.BM+LocalAd.BF)/2
LocalAd.H=(LocalAd.HM+LocalAd.HF)/2

plotting_values=data.frame(mean=c(mean(LocalAd.BM),mean(LocalAd.HM,na.rm=TRUE),mean(LocalAd.BF),mean(LocalAd.HF),mean(LocalAd.M),mean(LocalAd.F)),
                           Treatment=c("BM","HM","BF","HFL","Mavg","Favg"),
                           upper_CI=c(quantile(LocalAd.BM,0.975),quantile(LocalAd.HM,0.975,na.rm=TRUE),quantile(LocalAd.BF,0.975),quantile(LocalAd.HF,0.975),quantile(LocalAd.M,0.975),quantile(LocalAd.F,0.975)),
                           lower_CI=c(quantile(LocalAd.BM,0.025),quantile(LocalAd.HM,0.025,na.rm=TRUE),quantile(LocalAd.BF,0.025),quantile(LocalAd.HF,0.025),quantile(LocalAd.M,0.025),quantile(LocalAd.F,0.025)),
                           median=c(quantile(LocalAd.BM,0.5),quantile(LocalAd.HM,0.5,na.rm=TRUE),quantile(LocalAd.BF,0.5),quantile(LocalAd.HF,0.5),quantile(LocalAd.M,0.5),quantile(LocalAd.F,0.5)))
#GRAPH###
ggplot(plotting_values,aes(x=Treatment,y=as.numeric(mean),fill=Treatment,colour=Treatment))+
  geom_point(size=3,alpha=.9,shape=21)+ 
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI,colour=Treatment), width=.4,position=position_dodge(.9))+
  scale_fill_manual(values=c("red","#3e86ff","black","red","#3e86ff","black"))+
  scale_colour_manual(values=c("red","#3e86ff","black","red","#3e86ff","black"))+#"#f6564b
  geom_vline(xintercept=3.5,colour="grey",linewidth=0.5)+
  ggtitle("GroupNorm Relative.Fit2\n WS1L-WS2L/max(WS1L,WS2L)")+ylab("WS1L-WS2L/max(WS1L,WS2L) ")+xlab("Treatment Pool")+
  theme_classic() +theme(text = element_text(size = 14,color="black"),legend.position="none",panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7))+
  theme(axis.text.x = element_text(size = 12,color="black"), axis.text.y = element_text(size = 12,color="black"))



boot_diff.max<-function(data,chr,pop,col){
  #this substracts two bootstrapped means. (Note the first chr always subtracted from second ignore FL/ML notation)
  #divides by maximum
  Benefit.boot=replicate(1000,{
    #pull out replicates per population. The default size is n. So it might be uncessary to calculate length?
    samp.B1=sample(data[data$Chromosome==chr[1]&data$Population==pop[1],"Line.ID"],size=length(data[data$Chromosome==chr[1]&data$Population==pop[1],"Line.ID"]),replace=TRUE)#32
    #samp.B2=sample(data[data$Chromosome==chr[1]&data$Population==pop[2],"Line.ID"],size=length(data[data$Chromosome==chr[1]&data$Population==pop[2],"Line.ID"]),replace=TRUE)#30
    #samp.B3=sample(data[data$Chromosome==chr[1]&data$Population==pop[3],"Line.ID"],size=length(data[data$Chromosome==chr[1]&data$Population==pop[3],"Line.ID"]),replace=TRUE)#29
    #samp.BFL=append(samp.B1,samp.B2)
    #samp.BFL=append(samp.BFL,samp.B3)
    samp.BFL=samp.B1
    chr1.mean=mean(data[data$Line.ID%in%samp.BFL,col])
    
    #pull out replicates per population. The default size is n. So it might be uncessary to calculate length?
    samp.B1a=sample(data[data$Chromosome==chr[2]&data$Population==pop[1],"Line.ID"],size=length(data[data$Chromosome==chr[2]&data$Population==pop[1],"Line.ID"]),replace=TRUE)#32
    #samp.B2a=sample(data[data$Chromosome==chr[2]&data$Population==pop[2],"Line.ID"],size=length(data[data$Chromosome==chr[2]&data$Population==pop[2],"Line.ID"]),replace=TRUE)#30
    #samp.B3a=sample(data[data$Chromosome==chr[2]&data$Population==pop[3],"Line.ID"],size=length(data[data$Chromosome==chr[2]&data$Population==pop[3],"Line.ID"]),replace=TRUE)#29
    #samp.BML=append(samp.B1a,samp.B2a)
    #samp.BML=append(samp.BML,samp.B3a)
    samp.BML=samp.B1a
    chr2.mean=mean(data[data$Line.ID%in%samp.BML,col])
    
    Benefit=chr1.mean-chr2.mean/(max(chr1.mean,chr2.mean))
    Benefit
  } )
  return(Benefit.boot)
}
####--------****************----- (use) Local Adaptation-RAW Fit2 Per Pop w/ Structure (use)----**********####
fitness_col="Relative.Fit2"
#bootstrapped Fit2 without structure
# Fit2.boot=replicate(1000,{
#   samp.ID=sample(SSL.fit$Line.ID,replace=TRUE)
#   data_summary(SSL.fit[SSL.fit$Line.ID%in%samp.ID,],fitness_col,c("Focal.Sex","Chromosome","Population"))$Relative.Fit2
# } )
# Fit2.data=data_summary(SSL.fit,fitness_col,c("Focal.Sex","Chromosome","Population"))
# Fit2.data$upper_CI=apply(Fit2.boot,1,quantile,0.975)
# Fit2.data$lower_CI=apply(Fit2.boot,1,quantile,0.025)
# Fit2.data$Sex.Chr=paste(Fit2.data$Chromosome,Fit2.data$Focal.Sex)

#bootstrp Fit2 WITH structure---
#fitness_col="Fit2"
Fit2.structure=replicate(1000,{
  #sample equally FL and ML from all corresponding blocks
  #the default # will be proportional to the total length/# of replicates done in each block per pool
  samp.ID.bl1=append(sample(SSL.fit[SSL.fit$Block==1&SSL.fit$Chromosome=="FL","Line.ID"],replace=TRUE),
                     sample(SSL.fit[SSL.fit$Block==1&SSL.fit$Chromosome=="ML","Line.ID"],replace=TRUE))
  samp.ID.bl2=append(sample(SSL.fit[SSL.fit$Block==2&SSL.fit$Chromosome=="FL","Line.ID"],replace=TRUE),
                     sample(SSL.fit[SSL.fit$Block==2&SSL.fit$Chromosome=="ML","Line.ID"],replace=TRUE))
  samp.ID.bl3=append(sample(SSL.fit[SSL.fit$Block==3&SSL.fit$Chromosome=="FL","Line.ID"],replace=TRUE),
                     sample(SSL.fit[SSL.fit$Block==3&SSL.fit$Chromosome=="ML","Line.ID"],replace=TRUE))
  samp.ID.bl4=append(sample(SSL.fit[SSL.fit$Block==4&SSL.fit$Chromosome=="FL","Line.ID"],replace=TRUE),
                     sample(SSL.fit[SSL.fit$Block==4&SSL.fit$Chromosome=="ML","Line.ID"],replace=TRUE))
  samp.ID.bl5=append(sample(SSL.fit[SSL.fit$Block==5&SSL.fit$Chromosome=="FL","Line.ID"],replace=TRUE),
                     sample(SSL.fit[SSL.fit$Block==5&SSL.fit$Chromosome=="ML","Line.ID"],replace=TRUE))
  
  samp.ID=append(append(append(samp.ID.bl1,samp.ID.bl2),append(samp.ID.bl3,samp.ID.bl4)),samp.ID.bl5)
  
  data_summary(SSL.fit[SSL.fit$Line.ID%in%samp.ID,],fitness_col,c("Focal.Sex","Chromosome","Population"))$Relative.Fit2
} )


Fit2.data=data_summary(SSL.fit[SSL.fit$Chromosome!="T4",],fitness_col,c("Focal.Sex","Chromosome","Population"))
Fit2.data$upper_CI=apply(Fit2.structure,1,quantile,0.975)
Fit2.data$lower_CI=apply(Fit2.structure,1,quantile,0.025)
Fit2.data$Sex.Chr=paste(Fit2.data$Chromosome,Fit2.data$Focal.Sex)
Fit2.data$median=apply(Fit2.structure,1,quantile,0.5)


###trying to get population Relative Fitness measures

R.Fit2.structure=replicate(1000,{
  #Noting the the default sample is N, this should proportionally sample each subgroup
  #Benign Female Block 1
  B1.FFL1=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  B2.FFL1=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  B3.FFL1=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  B1.FML1=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  B2.FML1=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  B3.FML1=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  
  B.F1=append(append(append(B1.FFL1,B2.FFL1),append(B3.FFL1,B1.FML1)),append(B2.FML1,B3.FML1))
  
  #Beningn Female Block 2
  B1.FFL2=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  B2.FFL2=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  B3.FFL2=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  B1.FML2=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  B2.FML2=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  B3.FML2=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  
  B.F2=append(append(append(B1.FFL2,B2.FFL2),append(B3.FFL2,B1.FML2)),append(B2.FML2,B3.FML2))
  
  #BENIGN MALE BLOCK 3
  B1.MFL3=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  B2.MFL3=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  B3.MFL3=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  B1.MML3=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  B2.MML3=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  B3.MML3=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  
  B.M3=append(append(append(B1.MFL3,B2.MFL3),append(B3.MFL3,B1.MML3)),append(B2.MML3,B3.MML3))
  
  #BENIGN MALE BLOCK 3
  B1.MFL4=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  B2.MFL4=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  B3.MFL4=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  B1.MML4=sample(SSL.fit[SSL.fit$Population=="B1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  B2.MML4=sample(SSL.fit[SSL.fit$Population=="B2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  B3.MML4=sample(SSL.fit[SSL.fit$Population=="B3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  
  B.M4=append(append(append(B1.MFL4,B2.MFL4),append(B3.MFL4,B1.MML4)),append(B2.MML4,B3.MML4))
  
  #HIGH Female Block 1
  H1.FFL1=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  H2.FFL1=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  H3.FFL1=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  H1.FML1=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  H2.FML1=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  H3.FML1=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="1","Rep.ID"],replace=TRUE)
  
  H.F1=append(append(append(H1.FFL1,H2.FFL1),append(H3.FFL1,H1.FML1)),append(H2.FML1,H3.FML1))
  
  #HIGH Female block 2
  H1.FFL2=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  H2.FFL2=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  H3.FFL2=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  H1.FML2=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  H2.FML2=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  H3.FML2=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="F"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="2","Rep.ID"],replace=TRUE)
  
  H.F2=append(append(append(H1.FFL2,H2.FFL2),append(H3.FFL2,H1.FML2)),append(H2.FML2,H3.FML2))
  
  #HIGH MALE BLOCK 3
  H1.MFL3=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  H2.MFL3=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  H3.MFL3=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  H1.MML3=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  H2.MML3=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  H3.MML3=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="3","Rep.ID"],replace=TRUE)
  
  H.M3=append(append(append(H1.MFL3,H2.MFL3),append(H3.MFL3,H1.MML3)),append(H2.MML3,H3.MML3))
  
  #HIGH MALE BLOCK 4
  H1.MFL4=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  H2.MFL4=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  H3.MFL4=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  H1.MML4=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  H2.MML4=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  H3.MML4=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="4","Rep.ID"],replace=TRUE)
  
  H.M4=append(append(append(H1.MFL4,H2.MFL4),append(H3.MFL4,H1.MML4)),append(H2.MML4,H3.MML4))
  
  #HIGH MALE BLOCK 5
  H1.MFL5=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="5","Rep.ID"],replace=TRUE)
  H2.MFL5=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="5","Rep.ID"],replace=TRUE)
  H3.MFL5=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="FL"&SSL.fit$Block=="5","Rep.ID"],replace=TRUE)
  H1.MML5=sample(SSL.fit[SSL.fit$Population=="H1"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="5","Rep.ID"],replace=TRUE)
  H2.MML5=sample(SSL.fit[SSL.fit$Population=="H2"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="5","Rep.ID"],replace=TRUE)
  H3.MML5=sample(SSL.fit[SSL.fit$Population=="H3"&SSL.fit$Focal.Sex=="M"&SSL.fit$Chromosome=="ML"&SSL.fit$Block=="5","Rep.ID"],replace=TRUE)
  
  H.M5=append(append(append(H1.MFL5,H2.MFL5),append(H3.MFL5,H1.MML5)),append(H2.MML5,H3.MML5))
  
  #saving these separate so easier to check if correct
  samp.B=append(append(B.F1,B.F2),append(B.M3,B.M4))
  samp.H=append(append(append(H.F1,H.F2),append(H.M3,H.M4)),H.M5)
  samp.BH=append(samp.B,samp.H)
  #put together the full sample in it's correct proportions
  #samp=SSL.fit[SSL.fit$ep.ID%in%samp.BH,]
  
 
  # 
  # ##For each sample, make sure the appropriate Fit2 is turned into Relative.Fitness
  # samp$Relative.Fitness=0
  # samp[samp$Rep.ID%in%B.F1,"Relative.Fitness"]= samp[samp$Rep.ID%in%B.F1,fitness_col]
  # samp[samp$Rep.ID%in%B.F2,"Relative.Fitness"]= samp[samp$Rep.ID%in%B.F2,fitness_col]
  # samp[samp$Rep.ID%in%B.M3,"Relative.Fitness"]= samp[samp$Rep.ID%in%B.M3,fitness_col]
  # samp[samp$Rep.ID%in%B.M4,"Relative.Fitness"]= samp[samp$Rep.ID%in%B.M4,fitness_col]
  # samp[samp$Rep.ID%in%H.F1,"Relative.Fitness"]= samp[samp$Rep.ID%in%H.F1,fitness_col]
  # samp[samp$Rep.ID%in%H.F2,"Relative.Fitness"]= samp[samp$Rep.ID%in%H.F2,fitness_col]
  # samp[samp$Rep.ID%in%H.M3,"Relative.Fitness"]= samp[samp$Rep.ID%in%H.M3,fitness_col]
  # samp[samp$Rep.ID%in%H.M4,"Relative.Fitness"]= samp[samp$Rep.ID%in%H.M4,fitness_col]
  # samp[samp$Rep.ID%in%H.M5,"Relative.Fitness"]= samp[samp$Rep.ID%in%H.M5,fitness_col] 
  #this will output results? The order doesn't matter here because necessary changes have already been made
  #All of the raw Fit2 has now been strandardized by appropriate control, so can summarize the MEAN
  #data_summary(samp,"Relative.Fitness",c("Chromosome","Treatment","Focal.Sex","Population"))$Relative.Fitness
  data_summary(SSL.fit[SSL.fit$Rep.ID%in%samp.BH,],"Relative.Fit2",c("Focal.Sex","Chromosome","Population"))$Relative.Fit2
})

#note that even though I standardize by pop I don't take the mean by pop


#do this to get the mean Relative.Fit2 for each group which the above should be equivalent to. Make sure ORDER is same as above boot
Fit2.data=data_summary(SSL.fit[SSL.fit$Chromosome!="T4",],fitness_col,c("Focal.Sex","Chromosome","Population"))
Fit2.data$upper_CI=apply(R.Fit2.structure,1,quantile,0.975)
Fit2.data$lower_CI=apply(R.Fit2.structure,1,quantile,0.025)
Fit2.data$median=apply(R.Fit2.structure,1,quantile,0.5)

Fit2.data$Focal.Sex<-revalue(Fit2.data$Focal.Sex,c("F"="Females","M"="Males"))
Fit2.data$Sex.Chr=paste(Fit2.data$Focal.Sex,Fit2.data$Chromosome)
Fit2.data$Population<-revalue(Fit2.data$Population,c("H1"="N1","H2"="N2","H3"="N3"))

title="Relative.Fit2 w/ ALL structure 95 CI"
####Graph
ggplot(subset(Fit2.data,Chromosome!="T4"),aes(x=Population,y=as.numeric(Relative.Fit2),fill=Sex.Chr))+
  geom_bar(position=position_dodge(),stat = "summary", fun = "mean")+ 
  geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), colour="black",width=.4,position=position_dodge(.9))+
  #scale_fill_manual(values=c("darkred","pink","darkblue","#3e86ff"))+
  scale_fill_manual(values=c("#f6564b","#3e86ff","#f6564b","#3e86ff"))+
  scale_colour_manual(values=c("black","black"))+#"#f6564b
  geom_vline(xintercept=3.5,colour="grey",linewidth=0.5)+
  facet_wrap(~Focal.Sex)+
  ggtitle(title)+ylab("Mean Relative Fitness ")+xlab("Population")+
  theme(text = element_text(size = 14,color="black"),legend.position="top",panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7))+
  theme(axis.text.x = element_text(size = 12,color="black"), axis.text.y = element_text(size = 12,color="black"))

###now WSL1-WSL2/Max(WSL1,WSL2)
#from above F data is 1:12, and M data 13:24
populations=unique(Fit2.data$Population)
LocalAd.M=0
LocalAd.F=0

counterF=1
counterM=13
  for(pops in populations){
    #FL-ML
    F.diff=R.Fit2.structure[counterF,]-R.Fit2.structure[counterF+6,]
    max.F=apply(rbind(R.Fit2.structure[counterF,],R.Fit2.structure[counterF+6,]),2,max)
    LocalAd.F=rbind(LocalAd.F,F.diff/max.F)
    #F.diff.boot=rbind(F.diff.boot,F.diff)
    #max.F.boot=rbind(max.F.boot,max.F)
    
    counterF=counterF+1
    
    #ML-FL (FL comes first)
    M.diff=R.Fit2.structure[counterM+6,]-R.Fit2.structure[counterM,]
    max.M=apply(rbind(R.Fit2.structure[counterM,],R.Fit2.structure[counterM+6,]),2,max)
    #M.diff.boot=rbind(M.diff.boot,M.diff)
    #max.M.boot=rbind(max.M.boot,max.M)
    LocalAd.M=rbind(LocalAd.M,M.diff/max.M)
    counterM=counterM+1
  }
#remove blank row
  LocalAd.M= LocalAd.M[-1,]
  LocalAd.F= LocalAd.F[-1,]
  #Getting the "ACTUAL" population values
  counterF=1
  counterM=13
  mean.LocalAd.M=0
  mean.LocalAd.F=0
  for(pops in populations){
    #FL-ML
    F.diff=Fit2.data[counterF,"Relative.Fit2"]-Fit2.data[counterF+6,"Relative.Fit2"]
    max.F=apply(rbind(Fit2.data[counterF,"Relative.Fit2"],Fit2.data[counterF+6,"Relative.Fit2"]),2,max)
    mean.LocalAd.F=rbind(mean.LocalAd.F,F.diff/max.F)
    #F.diff.boot=rbind(F.diff.boot,F.diff)
    #max.F.boot=rbind(max.F.boot,max.F)
    
    counterF=counterF+1
    
    #ML-FL (FL comes first)
    M.diff=Fit2.data[counterM+6,"Relative.Fit2"]-Fit2.data[counterM,"Relative.Fit2"]
    max.M=apply(rbind(Fit2.data[counterM,"Relative.Fit2"],Fit2.data[counterM+6,"Relative.Fit2"]),2,max)
    #M.diff.boot=rbind(M.diff.boot,M.diff)
    #max.M.boot=rbind(max.M.boot,max.M)
    mean.LocalAd.M=rbind(mean.LocalAd.M,M.diff/max.M)
    counterM=counterM+1
  }
  
  #get actual mean WSL1-WSL2/Max(WSL1,WSL2)--Fit2 data will have mean
  #for F and THEN for M calculates mean produces 12 values
  mean.LocalAd=apply(rbind(LocalAd.F,LocalAd.M),1,mean)
  centre.LocalAd=append(mean.LocalAd.F[2:7],mean.LocalAd.M[2:7])
  #try to calculate from actual not bootstrap values

  LocalAd.data=data.frame(Population=rep(populations,2),Sex=c("F","F","F","F","F","F","M","M","M","M","M","M"))
  LocalAd.data$Local.Ad=centre.LocalAd
  LocalAd.data$upper_CI=append(apply(LocalAd.F,1,quantile,0.975),apply(LocalAd.M,1,quantile,0.975))  
  LocalAd.data$lower_CI=append(apply(LocalAd.F,1,quantile,0.025),apply(LocalAd.M,1,quantile,0.025))  
  LocalAd.data$median=append(apply(LocalAd.F,1,quantile,0.5),apply(LocalAd.M,1,quantile,0.5))  
  
  #LocalAd.data$lower_CI=append(apply(LocalAd.F,1,quantile,0.0025),apply(LocalAd.M,1,quantile,0.0025))  
  BF.LA=apply(LocalAd.F[1:3,],2,mean)
  NF.LA=apply(LocalAd.F[4:6,],2,mean)
  BM.LA=apply(LocalAd.M[1:3,],2,mean)
  NM.LA=apply(LocalAd.M[4:6,],2,mean)
  
    plotting_values2=data.frame(Population=c("B.Avg","N.Avg","B.Avg","N.Avg"),Sex=c("F","F","M","M"),
               Local.Ad=c(mean(BF.LA),mean(NF.LA),mean(BM.LA),mean(NM.LA)),
                              median=c(quantile(BF.LA,0.5),quantile(NF.LA,0.5),quantile(BM.LA,0.5),quantile(NM.LA,0.5)),
                              upper_CI=c(quantile(BF.LA,0.975),quantile(NF.LA,0.975),quantile(BM.LA,0.975),quantile(NM.LA,0.975)),
                              lower_CI=c(quantile(BF.LA,0.025),quantile(NF.LA,0.025),quantile(BM.LA,0.025),quantile(NM.LA,0.025)))
  
  plotting_values=rbind(LocalAd.data,plotting_values2)
  plotting_values$Sex<-revalue(plotting_values$Sex,c("F"="Females","M"="Males"))
  plotting_values$Treatment=c("BENIGN","BENIGN","BENIGN","NOVEL","NOVEL","NOVEL","BENIGN","BENIGN","BENIGN","NOVEL","NOVEL","NOVEL","BENIGN","NOVEL","BENIGN","NOVEL")
  
  Fig1C=ggplot(plotting_values,aes(x=Population,y=as.numeric(Local.Ad),fill=Sex,colour=Sex))+
    geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), colour="black",width=.6,stat="identity",position=position_dodge(.4),size=0.5)+
    geom_point(position=position_dodge(.4),stat = "identity",shape=21,size=3)+ 
    scale_fill_manual(values=c("darkred","#3e86ff"))+
    scale_colour_manual(values=c("black","black"))+#"#f6564b
    #geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), colour="black",width=.4,stat="identity",position=position_dodge(.9))+
    geom_vline(xintercept=6.5,colour="black",linewidth=0.5)+
    facet_wrap(~Treatment,scales="free_x")+
    geom_hline(yintercept=0,colour="grey",linewidth=0.5,linetype="dashed")+
    ggtitle("Local Adaptation-(median) Relative.Fit2 w/structure\nWSL1-WSL2/Max(WSL1,WSL2)")+ylab("Local Adaptation")+xlab("Population")+
    theme_classic() +theme(text = element_text(size = 14,color="black"),legend.position="top",panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7))+
    theme(axis.text.x = element_text(size = 12,color="black"), axis.text.y = element_text(size = 12,color="black"))
  
  
  AltFig1C= ggplot(plotting_values,aes(x=Population,y=as.numeric(Local.Ad),fill=Sex,colour=Sex))+
    geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), colour="black",width=.5,stat="identity",position=position_dodge(.5),size=0.5)+
    geom_point(position=position_dodge(.5),stat = "identity",shape=21,size=3)+ 
    scale_fill_manual(values=c("darkred","#3e86ff"))+
    scale_colour_manual(values=c("black","black"))+#"#f6564b
   geom_vline(xintercept=6.5,colour="black",linewidth=0.5)+
    geom_hline(yintercept=0,colour="grey",linewidth=0.5,linetype="dashed")+
    ggtitle("Local Adaptation-(median) Relative.Fit2 w/structure\nWSL1-WSL2/Max(WSL1,WSL2)")+ylab("Local Adaptation")+xlab("Population")+
    theme_classic() +theme(text = element_text(size = 14,color="black"),legend.position="top",panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7))+
    theme(axis.text.x = element_text(size = 12,color="black"), axis.text.y = element_text(size = 12,color="black"))
  
  
  #note we only actually care about the sex diff in LA we don't need the weighted diff (the LA is already weighted)
  LA.diff=LocalAd.M-LocalAd.F
  LA.centre.diff=mean.LocalAd.M[2:7]-mean.LocalAd.F[2:7]
  
  #again, this value does not have the "real" middle average
  plotting_values=data.frame(Population=populations,
                  sex.diff=LA.centre.diff, median=apply(LA.diff,1,quantile,0.5),
                  upper_CI=apply(LA.diff,1,quantile,0.975),lower_CI= apply(LA.diff,1,quantile,0.025))
                  
  ##try to get an averaged B vs N 
  #remember LA.diff has 6 rows for each population. It takes the PAIRED diff (which is appropriate)
  B.LA1=append(append(LA.diff[1,],LA.diff[2,]),LA.diff[3,])#each population contributes equally, NOT weighted
  N.LA1=append(append(LA.diff[4,],LA.diff[5,]),LA.diff[6,])
  
  #here averaging and then comparing
  B.LA2=apply(LA.diff[1:3,],2,mean)
  N.LA2=apply(LA.diff[4:6,],2,mean)
  
  #chooese one
  B.LA=B.LA2
  N.LA=N.LA2
  #Now calculate the means and such
  
  plotting_values1=data.frame(Population=c("ZB.Avg","ZN.Avg"), sex.diff=c(mean(B.LA),mean(N.LA)),
                             median=c(quantile(B.LA,0.5),quantile(N.LA,0.5)), upper_CI=c(quantile(B.LA,0.975),quantile(N.LA,0.975)),
                             lower_CI=c(quantile(B.LA,0.025),quantile(N.LA,0.025)))
  plotting_values2=data.frame(Population=c("B1","B2","B3","N1","N2","N3"),
                             sex.diff=LA.centre.diff,median=apply(LA.diff,1,quantile,0.5),
                             upper_CI=apply(LA.diff,1,quantile,0.975),lower_CI= apply(LA.diff,1,quantile,0.025))
  plotting_values=rbind(plotting_values1,plotting_values2)
  plotting_values$Treatment=c("BENIGN","NOVEL","BENIGN","BENIGN","BENIGN","NOVEL","NOVEL","NOVEL")
  Fig1D=ggplot(plotting_values,aes(x=Population,y=as.numeric(sex.diff)))+
    geom_point(size=3,alpha=.9,shape=19)+ facet_wrap(~Treatment,scales="free_x")+
    geom_errorbar(aes(ymin=lower_CI, ymax=upper_CI), colour="black",width=.4,position=position_dodge(.9))+
    #scale_colour_manual(values=c("darkred","red","pink","#3e86ff","blue","darkblue"))+
   #geom_vline(xintercept=4.5,colour="black",linewidth=0.5)+ 
    geom_hline(yintercept=0,colour="grey",linewidth=0.5,linetype="dashed")+
    ggtitle("Local Adaptation-(mean)Relative.Fit2 w/structure\nLA.M-LA.F avg bt pops")+ylab("Local Adaptation (male - female)")+xlab("Population")+
    theme_classic() +theme(text = element_text(size = 14,color="black"),legend.position="none",panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7))+
    theme(axis.text.x = element_text(size = 12,color="black"), axis.text.y = element_text(size = 12,color="black"))
  
  
  
  
  
  
################---Saving Graphs---(use all)----

##Figure 1
#composed of Fig1A,B,C,D
#first remove title
A=Fig1A+theme(plot.title = element_blank(),axis.title.x=element_blank(),axis.title.y=element_text(size = 13,color="black"),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7))+
  scale_fill_manual(values=c("#E41A1C","#3e86ff","red"))
B=Fig1B+theme(plot.title = element_blank(),axis.title.x=element_blank(),axis.title.y=element_text(size = 13,color="black"),
              legend.margin=margin(-3,0,-3,0),legend.title=element_text(size = 12,color="black"),
              panel.border = element_rect(colour = "black", fill=NA, linewidth=0.7))+
  scale_fill_manual(values=c("#E41A1C","#3e86ff","red"))
C=Fig1C+theme(plot.title = element_blank(),axis.text.x = element_text(size = 9,color="black"),
              axis.title=element_text(size = 12,color="black"),
              legend.position="none")
D=Fig1D+theme(plot.title = element_blank(),axis.text.x = element_text(size = 10,color="black"),
              axis.title=element_text(size = 12,color="black"),axis.text.y = element_text(size = 10,color="black"))
# E=AltFig1C+theme(plot.title = element_blank(),axis.text.x = element_text(size = 9,color="black"),
#               axis.title=element_text(size = 12,color="black"),axis.text.y = element_text(size = 10,color="black"),
#               legend.margin=margin(-3,0,-3,0),legend.title=element_text(size = 12,color="black"))
#  # scale_fill_manual(values=c("#E41A1C","#3e86ff","red"))



Relative.Fitness=ggarrange(A,B,ncol=2,labels=c("A","B"),
          widths=c(0.6,1),common.legend=TRUE)
Local.Adaptation=ggarrange(C,D,ncol=2,labels=c("C","D"),
          widths=c(1,0.9))
Fig1=ggarrange(Relative.Fitness,Local.Adaptation,nrow=2)
Fig1
#save the png
ggsave(
  "SSL_Fit_Fig2r2.png",
  Fig1,
  width     = 9,
  height    = 7.25,
  dpi       = 1200)
###-rmf Fig2--####
##here is where I would make any adjustments to title size etc
ggsave(
  "SSL_Fit_Fig2.png",
  Fig2,
  width     = 4.5,
  height    = 4.25,
  dpi       = 1000)
####Mate Harm 3#####

F3A=mate.harm.fec.plot+ylab("Offspring per Female Mating Partner")+#facet_wrap(~Treatment,scale="free_x")+
  theme(axis.text.x = element_text(size = 9,color="black"), axis.title.y=element_text(size=9,colour="black"),
         axis.text.y = element_text(size = 9,color="black"),legend.margin=margin(-1,0,0,-1))+#legend.key.size = unit(1, "lines"))
  scale_fill_manual(values=c("#E41A1C","#3e86ff","red"))+scale_colour_manual(values=c("#E41A1C","#3e86ff","red"))

F3B=mate.harm.surv.plot+ylab("Female Mating Partner Survival")+#facet_wrap(~Treatment,scale="free_x")+ 
scale_fill_manual(values=c("#E41A1C","#3e86ff","red"))+scale_colour_manual(values=c("#E41A1C","#3e86ff","red"))+
  theme(axis.text.x = element_text(size = 9,color="black"), 
        axis.text.y = element_text(size = 9,color="black"),axis.title.y=element_text(size=9,colour="black"))

F3C=AltFig3C+ ylim(-0.32,0.18)+geom_point(inherit.aes=FALSE,data=plotting_values[1:4,],aes(y=real.means,x=Chromosome,colour=Chromosome,fill=Chromosome),size=5,shape=4)+
  theme(axis.text.x = element_text(size = 9,color="black"), 
        axis.text.y = element_text(size = 9,color="black"),axis.title.y=element_text(size=9,colour="black"))

  


Fig3top=ggarrange(F3B,F3A,labels=c("A","B"),ncol=2,widths=c(1,1),common.legend=TRUE)
#Fig3top

Fig3=ggarrange(Fig3top,F3C,labels=c("","C"),nrow=2,heights=c(1,1.1))
Fig3
#save the png
ggsave(
  "SSL_Fit_Fig3v2.png",
 Fig3,
  width     = 5.30,
  height    = 7.06,
  dpi       = 1200)