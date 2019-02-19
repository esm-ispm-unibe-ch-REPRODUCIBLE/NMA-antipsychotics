
pdf("Meta-regression and sensitivity analyses graphs.pdf")
sink("Meta-regression and sensitivity analyses output.txt")

#*********************************************************************************
#         Unadjusted analysis for comparison               
#*********************************************************************************

cat("\n \n \n UNADJUSTED analysis \n \n \n")

SchizoEFF=SchizoEFF[order(SchizoEFF$Study_No),]

#transform the data into a list suitable for JAGS analysis
NMAdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFF,type="cont",reference = "Placebo")
#run Jags and create a jags object
NMAinJAGS<- jags.parallel(data = NMAdataContinuous, inits = NULL,
                 parameters.to.save = c("SMD","SMD.ref","tau", "SUCRA"), n.chains = 2, n.iter = 100000,
                 n.burnin = 5000,DIC=F,n.thin=10,
                 model.file = modelNMAContinuous)
print(NMAinJAGS)
save(NMAinJAGS,file="NMAinJAGSmorethan100.RData",envir = .GlobalEnv)

#check chain mixing
#traceplot(NMAinJAGS,varname="tau" )
#traceplot(NMAinJAGS,varname="SMD.ref" )

#forestplot against placebo
y=NMAinJAGS$BUGSoutput$mean$SMD.ref
ES=y
seES=NMAinJAGS$BUGSoutput$sd$SMD.ref
#then get the league table and save it
leaguetable=out.jagsNMA.results(NMAinJAGS,"SMD",F, treatnames=sort(unique(SchizoEFF$Drug)))
leaguetableEFF100=leaguetable$leaguetable
write.csv(leaguetableEFF100,file="leaguetableEFFmorethan100.csv")


##################################################################################
#### ANALYSIS PRimary outcome efficacy meta-regressions USING JAGS   
##################################################################################



#*********************************************************************************
#         Year               
#*********************************************************************************
cat("\n \n \n YEAR \n \n \n")
#exclude studies that don't report the year
SchizoEFFyear=SchizoEFF[!is.na(SchizoEFF$Year),] 
round(tapply(SchizoEFF$Year,SchizoEFF$Drug,mean))
range(SchizoEFFyear$Year)
table(SchizoEFFyear$Drug)
tapply(SchizoEFFyear$OverallEfficN,SchizoEFFyear$Drug,sum)
range(table(SchizoEFFyear$Study_No))
#transform the data into a list suitable for JAGS analysis
NMRdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFFyear,type="cont",reference = "Placebo",othervar=2016-SchizoEFFyear$Year)

#run Jags and create a jags object
NMRinJAGSyear<- jags.parallel(data = NMRdataContinuous, inits = NULL,
                 parameters.to.save = c("SMD.ref","tau", "SMD","SUCRA","B","tauB"), n.chains = 2, n.iter = 100000,
                 n.burnin = 5000,DIC=F,n.thin=10,
                 model.file = modelNMRContinuousCommonB)
save(NMRinJAGSyear,file="NMRinJAGSyear.RData",envir = .GlobalEnv)


#check chain mixing
#traceplot(NMRinJAGSyear,varname="tau" )
#traceplot(NMRinJAGSyear,varname="tauB" )
#traceplot(NMRinJAGSyear,varname="B" )
#traceplot(NMRinJAGSyear,varname="SMD.ref" )


#Present results

y=NMRinJAGSyear$BUGSoutput$mean$SMD.ref
ES=y
seES=NMRinJAGSyear$BUGSoutput$sd$SMD.ref
drugs=sort(unique(SchizoEFFyear$Drug))
forest(metagen(ES[order(y)],seES[order(y)],studlab=drugs[order(y)],comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,layout="JAMA",xlab="SMD meta-regressions for year=2016")

NMRinJAGSyear$BUGSoutput$summary[1,c(3,5,7)]

leaguetable=out.jagsNMA.results(NMRinJAGSyear,"SMD",treatnames=sort(unique(SchizoEFFyear$Drug)))
leaguetableEFFyear=leaguetable$leaguetable
write.csv(leaguetableEFFyear,file="leaguetableEFFyear.csv")



#*********************************************************************************
#        age               
#*********************************************************************************

#exclude studies that don't report the year
SchizoEFFage=SchizoEFF[!is.na(SchizoEFF$MeanAge),] 
round(tapply(SchizoEFFage$MeanAge,SchizoEFFage$Drug,mean))
range(SchizoEFFage$MeanAge)
mean(SchizoEFFage$MeanAge)
table(SchizoEFFage$Drug)
tapply(SchizoEFFage$OverallEfficN,SchizoEFFage$Drug,sum)
range(table(SchizoEFFage$Study_No))
#transform the data into a list suitable for JAGS analysis
NMRdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFFage,type="cont",reference = "Placebo",othervar=40-SchizoEFFage$MeanAge)

#run Jags and create a jags object
NMRinJAGSage<- jags.parallel(data = NMRdataContinuous, inits = NULL,
                     parameters.to.save = c("SMD.ref","tau", "SMD","SUCRA","b"), n.chains = 2, n.iter = 100000,
                     n.burnin = 5000,DIC=F,n.thin=10,
                     model.file = modelNMRContinuous)
save(NMRinJAGSage,file="NMRinJAGSage.RData",envir = .GlobalEnv)




#Present results

y=NMRinJAGSage$BUGSoutput$mean$SMD.ref
ES=y
seES=NMRinJAGSage$BUGSoutput$sd$SMD.ref
drugs=sort(unique(SchizoEFFage$Drug))
forest(metagen(ES[order(y)],seES[order(y)],studlab=drugs[order(y)],comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,layout="JAMA",xlab="SMD meta-regressions for age=40")

coef=NMRinJAGSage$BUGSoutput$mean$b[order(y)]
secoef=NMRinJAGSage$BUGSoutput$sd$b[order(y)]
forest(metagen(coef,secoef,studlab=drugs[order(y)],comb.fixed = F, comb.random = F), xlim=c(-0.1,0.1), digits=3,layout="JAMA",xlab="SMD meta-regressions for age=40")

leaguetable=out.jagsNMA.results(NMRinJAGSage,"SMD",treatnames=sort(unique(SchizoEFFyear$Drug)))
leaguetableEFFage=leaguetable$leaguetable
write.csv(leaguetableEFFage,file="leaguetableEFFage.csv")


#*********************************************************************************
#        Baseline Severity              
#*********************************************************************************

#exclude studies that don't report baseline

out=unique(SchizoEFF$Study_No[is.na(SchizoEFF$BaselineSeverity)])
keep=is.na(match(SchizoEFF$Study_No,out))
SchizoEFFbaseline=SchizoEFF[keep,] 
range(table(SchizoEFFpermale$Study_No))


round(tapply(SchizoEFF$BaselineSeverity,SchizoEFF$Drug,mean))
range(SchizoEFFbaseline$BaselineSeverity)
mean(SchizoEFFbaseline$BaselineSeverity)
table(SchizoEFFbaseline$Drug)
tapply(SchizoEFFbaseline$OverallEfficN,SchizoEFFbaseline$Drug,sum)
range(table(SchizoEFFbaseline$Study_ID))


#transform the data into a list suitable for JAGS analysis
NMRdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFFbaseline,type="cont",reference = "Placebo",othervar=94-SchizoEFFbaseline$BaselineSeverity)

#run Jags and create a jags object
NMRinJAGSbaseline<- jags.parallel(data = NMRdataContinuous, inits = NULL,
                    parameters.to.save = c("SMD.ref","tau", "SMD","SUCRA","b"), n.chains = 2, n.iter = 100000,
                    n.burnin = 5000,DIC=F,n.thin=10,
                    model.file = modelNMRContinuous)
save(NMRinJAGSbaseline,file="NMRinJAGSbaseline.RData",envir = .GlobalEnv)




#Present results

y=NMRinJAGSbaseline$BUGSoutput$mean$SMD.ref
ES=y
seES=NMRinJAGSbaseline$BUGSoutput$sd$SMD.ref
drugs=sort(unique(SchizoEFFbaseline$Drug))
forest(metagen(ES[order(y)],seES[order(y)],studlab=drugs[order(y)],comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,layout="JAMA",xlab="SMD meta-regressions for baseline severity=94")

coef=NMRinJAGSbaseline$BUGSoutput$mean$b[order(y)]
secoef=NMRinJAGSbaseline$BUGSoutput$sd$b[order(y)]
forest(metagen(coef,secoef,studlab=drugs[order(y)],comb.fixed = F, comb.random = F), xlim=c(-0.1,0.1), digits=3,layout="JAMA",xlab="SMD meta-regressions for baseline severity=94")

leaguetable=out.jagsNMA.results(NMRinJAGSbaseline,"SMD",treatnames=sort(unique(SchizoEFFbaseline$Drug)))
leaguetableEFFbaseline=leaguetable$leaguetable
write.csv(leaguetableEFFbaseline,file="leaguetableEFFbaseline.csv")


#*********************************************************************************
#        % males             
#*********************************************************************************

#exclude studies that don't report % males
SchizoEFF$permale=SchizoEFF$No_male/SchizoEFF$No_randomised
out=unique(SchizoEFF$Study_No[is.na(SchizoEFF$permale)])
keep=is.na(match(SchizoEFF$Study_No,out))
SchizoEFFpermale=SchizoEFF[keep,] 
range(table(SchizoEFFpermale$Study_No))

round(tapply(SchizoEFF$,SchizoEFF$Drug,mean))
range(SchizoEFFpermale$permale)
mean(SchizoEFFpermale$permale)
table(SchizoEFFpermale$Drug)
tapply(SchizoEFFpermale$OverallEfficN,SchizoEFFpermale$Drug,sum)
tapply(SchizoEFFpermale$permale,SchizoEFFpermale$Drug,range)


#transform the data into a list suitable for JAGS analysis
NMRdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFFpermale,type="cont",reference = "Placebo",othervar=0.5-permale)

#run Jags and create a jags object
NMRinJAGSpermale<- jags.parallel(data = NMRdataContinuous, inits = NULL,
                         parameters.to.save = c("SMD.ref","tau", "SMD","SUCRA","b"), n.chains = 2, n.iter = 100000,
                         n.burnin = 5000,DIC=F,n.thin=10,
                         model.file = modelNMRContinuous)
save(NMRinJAGSpermale,file="NMRinJAGSpermale.RData",envir = .GlobalEnv)


y=NMRinJAGSpermale$BUGSoutput$mean$SMD.ref
ES=y
seES=NMRinJAGSpermale$BUGSoutput$sd$SMD.ref
drugs=sort(unique(SchizoEFFpermale$Drug))
forest(metagen(ES,seES,studlab=drugs,comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,layout="JAMA",xlab="SMD meta-regressions for % male=50%")


#Present results

leaguetable=out.jagsNMA.results(NMRinJAGSpermale,"SMD",treatnames=sort(unique(SchizoEFFpermale$Drug)))
leaguetableEFFpermale=leaguetable$leaguetable
write.csv(leaguetableEFFpermale,file="leaguetableEFFpermale.csv")



#*********************************************************************************
#        Variance            
#*********************************************************************************

NMRdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFF,type="cont",reference = "Placebo",othervar=OverallEfficN)
variance=NMRdataContinuous$pooled.sd*NMRdataContinuous$pooled.sd/NMRdataContinuous$variab
NMRdataContinuous$variab=variance

#run Jags and create a jags object
NMRinJAGSvariance<- jags.parallel(data = NMRdataContinuous, inits = NULL,
                         parameters.to.save = c("SMD.ref","tau", "SMD","SUCRA","b"), n.chains = 2, n.iter = 100000,
                         n.burnin = 5000,DIC=F,n.thin=10,
                         model.file = modelNMRContinuous)
save(NMRinJAGSvariance,file="NMRinJAGSvariance.RData",envir = .GlobalEnv)


#check chain mixing
#traceplot(NMRinJAGSvariance,varname="tau" )
#traceplot(NMRinJAGSvariance,varname="SMD.ref" )
#traceplot(NMRinJAGSvariance,varname="b" )



#Present results

y=NMRinJAGSvariance$BUGSoutput$mean$SMD.ref
ES=y
seES=NMRinJAGSvariance$BUGSoutput$sd$SMD.ref
drugs=sort(unique(SchizoEFFvariance$Drug))
forest(metagen(ES[order(y)],seES[order(y)],studlab=drugs[order(y)],comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,layout="JAMA",xlab="SMD meta-regressions for an infinitely large study")

coef=NMRinJAGSvariance$BUGSoutput$mean$b[order(y)]
secoef=NMRinJAGSvariance$BUGSoutput$sd$b[order(y)]
forest(metagen(coef,secoef,studlab=drugs[order(y)],comb.fixed = F, comb.random = F), xlim=c(-0.1,0.1), digits=3,layout="JAMA",xlab="SMD meta-regression coefficients for variance")

leaguetable=out.jagsNMA.results(NMRinJAGSvariance,"SMD",treatnames=sort(unique(SchizoEFFSponsor$Drug)))
leaguetableEFFvariance=leaguetable$leaguetable
write.csv(leaguetableEFFvariance,file="leaguetableEFFvariance.csv")

#*********************************************************************************
#        Sponsoring             
#*********************************************************************************
#nr of studies with NA in sponsoring 
sum(table(is.na(SchizoEFF$Sponsor),SchizoEFF$Study_No)[2,]>0)
#IDs of studies with NA in sponsoring 
out=unique(SchizoEFF$Study_No[is.na(SchizoEFF$Sponsor)])
keep=is.na(match(SchizoEFF$Study_No,out))
SchizoEFFSponsor=SchizoEFF[keep,] 
range(table(SchizoEFFSponsor$Study_No))
#nr of times a drug is sponsored over the total number of times it appears

cat("\n nr of times a drug is sponsored over the total number of times it appears: \n " )
(tapply(SchizoEFFSponsor$Sponsor,SchizoEFFSponsor$Drug,sum))/table(SchizoEFFSponsor$Drug)


NMRdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFFSponsor,type="cont",reference = "Placebo")
NMRdataContinuousFAKE=make.jagsNMA.data(studyid=Study_No,t=Drug,y=Sponsor,sd=Sponsor,n=OverallEfficN,data=SchizoEFFSponsor,type="cont",reference = "Placebo")
NMRdataContinuous$Sponsor=NMRdataContinuousFAKE$y
rm(NMRdataContinuousFAKE)


#run Jags and create a jags object

NMRinJAGSSponsor<- jags.parallel(data = NMRdataContinuous, inits = NULL,
                         parameters.to.save = c("SMD.ref","tau", "SMD","SUCRA","Beta"), n.chains = 2, n.iter = 100000,
                         n.burnin = 5000,DIC=F,n.thin=10,
                         model.file = modelNMRContinuousSponsor)
save(NMRinJAGSSponsor,file="NMRinJAGSSponsor.RData",envir = .GlobalEnv)



#check chain mixing
#traceplot(NMRinJAGSSponsor,varname="tau" )
#traceplot(NMRinJAGSSponsor,varname="SMD.ref" )
#traceplot(NMRinJAGSSponsor,varname="Beta" )


#Present results

y=NMRinJAGSSponsor$BUGSoutput$mean$SMD.ref
ES=y
seES=NMRinJAGSSponsor$BUGSoutput$sd$SMD.ref
drugs=sort(unique(SchizoEFFSponsor$Drug))
forest(metagen(ES[order(y)],seES[order(y)],studlab=drugs[order(y)],comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,layout="JAMA",xlab="SMD meta-regressions equal sponsoring interests")

leaguetable=out.jagsNMA.results(NMRinJAGSSponsor,"SMD",treatnames=sort(unique(SchizoEFFSponsor$Drug)))
leaguetableEFFSponsor=leaguetable$leaguetable
write.csv(leaguetableEFFSponsor,file="leaguetableEFFSponsor.csv")
cat("\n Common sponsoring coefficient: \n " )
NMRinJAGSSponsor$BUGSoutput$summary[1,c(1,3,7)]

############################################################################################################################################
############################################################################################################################################
# SUBGROUP ANALYSES *** SUBGROUP ANALYSES *** SUBGROUP ANALYSES *** SUBGROUP ANALYSES *** SUBGROUP ANALYSES *** SUBGROUP ANALYSES *** SUBGRO
############################################################################################################################################
############################################################################################################################################


#*********************************************************************************
#        Length of follow-up             
#*********************************************************************************

#exclude studies that don't report floow-up

out=unique(SchizoEFF$Study_No[is.na(SchizoEFF$Duration)])
keep=is.na(match(SchizoEFF$Study_No,out))
SchizoEFFDuration=SchizoEFF[keep,] 

table(SchizoEFFDuration$Duration,SchizoEFFDuration$Drug)
mean(SchizoEFFDuration$Duration)

#####
### keep those with duration 4 to 8  weeks!!!!!!!!!!!!!!!!!!

out=unique(SchizoEFFDuration$Study_No[SchizoEFFDuration$Duration<4|SchizoEFFDuration$Duration>8])
keep=is.na(match(SchizoEFFDuration$Study_No,out))
SchizoEFFDuration=SchizoEFFDuration[keep,] 
table(SchizoEFFDuration$Duration,SchizoEFFDuration$Drug)
tapply(SchizoEFFDuration$Duration,SchizoEFFDuration$Drug,mean)

#transform the data into a list suitable for JAGS analysis
NMRdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFFDuration,type="cont",reference = "Placebo")

#run Jags and create a jags object
NMRinJAGSDuration<- jags.parallel(data = NMRdataContinuous, inits = NULL,
                         parameters.to.save = c("SMD.ref","tau", "SMD","SUCRA"), n.chains = 2, n.iter = 100000,
                         n.burnin = 5000,DIC=F,n.thin=10,
                         model.file = modelNMAContinuous)
save(NMRinJAGSDuration,file="NMRinJAGSDuration.RData",envir = .GlobalEnv)

#traceplot(NMRinJAGSDuration,varname="tau" )
#traceplot(NMRinJAGSDuration,varname="SMD.ref" )


#Present results

y=NMRinJAGSDuration$BUGSoutput$mean$SMD.ref
ES=y
seES=NMRinJAGSDuration$BUGSoutput$sd$SMD.ref
drugs=sort(unique(SchizoEFFDuration$Drug))
forest(metagen(ES[order(y)],seES[order(y)],studlab=drugs[order(y)],comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,layout="JAMA",xlab="SMD meta-regressions for Duration>=6 weeks")

leaguetable=out.jagsNMA.results(NMRinJAGSDuration,"SMD",treatnames=sort(unique(SchizoEFFDuration$Drug)))
leaguetableEFFDuration=leaguetable$leaguetable
write.csv(leaguetableEFFDuration,file="leaguetableEFFDuration.csv")

#*********************************************************************************
#        high risk of bias        
#*********************************************************************************
#nr of studies with high risk of bias
length(unique(SchizoEFF$Study_ID[SchizoEFF$HighRoB==0]))
length(unique(SchizoEFF$Study_ID[SchizoEFF$HighRoB==1]))

SchizoEFFRoB=SchizoEFF[SchizoEFF$HighRoB<1,] 
range(table(SchizoEFFRoB$Study_ID))

#transform the data into a list suitable for JAGS analysis
NMRdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFFRoB,type="cont",reference = "Placebo")

#run Jags and create a jags object
NMRinJAGSRoB<- jags.parallel(data = NMRdataContinuous, inits = NULL,
                         parameters.to.save = c("SMD.ref","tau", "SMD","SUCRA"), n.chains = 2, n.iter = 100000,
                         n.burnin = 5000,DIC=F,n.thin=10,
                         model.file = modelNMAContinuous)
save(NMRinJAGSRoB,file="NMRinJAGSRoB.RData",envir = .GlobalEnv)

#traceplot(NMRinJAGSRoB,varname="tau" )
#traceplot(NMRinJAGSRoB,varname="SMD.ref" )


#Present results

y=NMRinJAGSRoB$BUGSoutput$mean$SMD.ref
ES=y
seES=NMRinJAGSRoB$BUGSoutput$sd$SMD.ref
drugs=sort(unique(SchizoEFFRoB$Drug))
forest(metagen(ES[order(y)],seES[order(y)],studlab=drugs[order(y)],comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,layout="JAMA",xlab="SMD meta-regressions excluding studies with High RoB")

leaguetable=out.jagsNMA.results(NMRinJAGSRoB,"SMD",treatnames=sort(unique(SchizoEFFRoB$Drug)))
leaguetableEFFRoB=leaguetable$leaguetable
write.csv(leaguetableEFFRoB,file="leaguetableEFFRoB.csv")
#*********************************************************************************
#        imputed standard deviations       
#*********************************************************************************
#nr of studies with imputed SD
out=unique(SchizoEFF$Study_No[is.na(SchizoEFF$CalculatedSD)|SchizoEFF$CalculatedSD>2])
keep=is.na(match(SchizoEFF$Study_No,out))
SchizoEFFcalcSD=SchizoEFF[keep,] 

length(table(SchizoEFFcalcSD$Study_No))
min(table(SchizoEFFcalcSD$Study_No))


#transform the data into a list suitable for JAGS analysis
NMRdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFFcalcSD,type="cont",reference = "Placebo")


#run Jags and create a jags object
NMRinJAGScalcSD<- jags.parallel(data = NMRdataContinuous, inits = NULL,
                    parameters.to.save = c("SMD.ref","tau", "SMD","SUCRA"), n.chains = 2, n.iter = 100000,
                    n.burnin = 5000,DIC=F,n.thin=10,
                    model.file = modelNMAContinuous)
save(NMRinJAGScalcSD,file="NMRinJAGScalcSD.RData",envir = .GlobalEnv)

#traceplot(NMRinJAGScalcSD,varname="tau" )
#traceplot(NMRinJAGScalcSD,varname="SMD.ref" )


#Present results

y=NMRinJAGScalcSD$BUGSoutput$mean$SMD.ref
ES=y
seES=NMRinJAGScalcSD$BUGSoutput$sd$SMD.ref
drugs=sort(unique(SchizoEFFcalcSD$Drug))
forest(metagen(ES[order(y)],seES[order(y)],studlab=drugs[order(y)],comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,layout="JAMA",xlab="SMD meta-regressions excluding studies with imputed SD")

leaguetable=out.jagsNMA.results(NMRinJAGScalcSD,"SMD",treatnames=as.character(sort(unique(SchizoEFFcalcSD$treat))))
leaguetableEFFcalcSD=leaguetable$leaguetable
write.csv(leaguetableEFFcalcSD,file="leaguetableEFFcalcSD.csv")

#*********************************************************************************
#        exclude placebo-controlled studies     
#*********************************************************************************
#nhead to head studies
length(unique(SchizoEFF$Study_ID[SchizoEFF$Placebo==0]))
SchizoEFFHtH=SchizoEFF[SchizoEFF$Placebo==0,] 
range(table(SchizoEFFHtH$Study_No))
(table(SchizoEFFHtH$Drug))
#transform the data into a list suitable for JAGS analysis
#check for disconnections
SchizoEFFpairsHtH=pairwise(treat=Drug,mean=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN, data=SchizoEFFHtH, studlab = Study_No, sm= "SMD")
#network plot
EFFHtH=netmeta(SchizoEFFpairsHtH, tol.multiarm=0.01)
netgraph(EFFHtH, plastic=F, thickness="number.of.studies", multiarm = F, points=T, cex.points=table(SchizoEFFHtH$Drug)/4, col=1)

NMRdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFFHtH,type="cont", reference = "Haloperidol")


#run Jags and create a jags object
NMRinJAGSHtH<- jags.parallel(data = NMRdataContinuous, inits = NULL,
                       parameters.to.save = c("SMD.ref","tau", "SMD","SUCRA"), n.chains = 2, n.iter = 100000,
                       n.burnin = 5000,DIC=F,n.thin=10,
                       model.file = modelNMAContinuous)
save(NMRinJAGSHtH,file="NMRinJAGSHtH.RData",envir = .GlobalEnv)

#traceplot(NMRinJAGSHtH,varname="tau" )
#traceplot(NMRinJAGSHtH,varname="SMD.ref" )

#Present results
y=NMRinJAGSHtH$BUGSoutput$mean$SMD.ref
ES=y
seES=NMRinJAGSHtH$BUGSoutput$sd$SMD.ref
drugs=sort(unique(SchizoEFFHtH$Drug))
forest(metagen(ES[order(y)],seES[order(y)],studlab=drugs[order(y)],comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,
layout="JAMA",xlab="SMD from head-to-head studies only")

leaguetable=out.jagsNMA.results(NMRinJAGSHtH,"SMD",treatnames=(sort(unique(SchizoEFFHtH$Drug))))
leaguetableEFFHtH=leaguetable$leaguetable
write.csv(leaguetableEFFHtH,file="leaguetableEFFHtH.csv")

#*********************************************************************************
#        exclude placebo arms  
#*********************************************************************************

SchizoEFFnoPlac=SchizoEFF[SchizoEFF$Drug!="Placebo",] 
length(table(SchizoEFFnoPlac$Study_No))
range(table(SchizoEFFnoPlac$Study_No))
#Exclude the single arm studies
SchizoEFFnoPlac=pool.arms(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFFnoPlac,type="cont")
length(table(SchizoEFFnoPlac$id))

#transform the data into a list suitable for JAGS analysis


NMRdataContinuous=make.jagsNMA.data(studyid=id,t=treat,y=y,sd=sd,n=n,data=SchizoEFFnoPlac,type="cont",reference = "Haloperidol")


#run Jags and create a jags object
NMRinJAGSnoPlac<- jags.parallel(data = NMRdataContinuous, inits = NULL,
                       parameters.to.save = c("SMD.ref","tau", "SMD","SUCRA"), n.chains = 2, n.iter = 100000,
                       n.burnin = 5000,DIC=F,n.thin=10,
                       model.file = modelNMAContinuous)
save(NMRinJAGSnoPlac,file="NMRinJAGSnoPlac.RData",envir = .GlobalEnv)

#traceplot(NMRinJAGSnoPlac,varname="tau" )
#traceplot(NMRinJAGSnoPlac,varname="SMD.ref" )

leaguetable=out.jagsNMA.results(NMRinJAGSnoPlac,"SMD",treatnames=as.character(sort(unique(SchizoEFFnoPlac$treat))))
leaguetableNoPlac=leaguetable$leaguetable
write.csv(leaguetableNoPlac,file="leaguetableNoPlac.csv")

#*********************************************************************************
#        Completers        
#*********************************************************************************
#nr of studies with ITT

length(unique(SchizoEFF$Study_ID[SchizoEFF$ComplITT==1]))
out=unique(SchizoEFF$Study_No[SchizoEFF$ComplITT>1 | is.na(SchizoEFF$ComplITT) ])
keep=is.na(match(SchizoEFF$Study_No,out)) 
SchizoEFFITT=SchizoEFF[keep,] 
length(table(SchizoEFFITT$Study_No))



#transform the data into a list suitable for JAGS analysis
NMRdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFFITT,type="cont",reference = "Placebo")

#run Jags and create a jags object
NMRinJAGSITT<- jags.parallel(data = NMRdataContinuous, inits = NULL,
                    parameters.to.save = c("SMD.ref","tau", "SMD","SUCRA"), n.chains = 2, n.iter = 100000,
                    n.burnin = 5000,DIC=F,n.thin=10,
                    model.file = modelNMAContinuous)
save(NMRinJAGSITT,file="NMRinJAGSITT.RData",envir = .GlobalEnv)

#traceplot(NMRinJAGSITT,varname="tau" )
#traceplot(NMRinJAGSITT,varname="SMD.ref" )


leaguetable=out.jagsNMA.results(NMRinJAGSITT,"SMD",treatnames=sort(unique(SchizoEFFRoB$Drug)))
leaguetableEFFITT=leaguetable$leaguetable
write.csv(leaguetableEFFITT,file="leaguetableEFFITT.csv")

#*********************************************************************************
#        Exclude studies with unfair dose comparisons       
#*********************************************************************************

unfair.dose.fun=function(dosevector){
  dosevector=as.vector(dosevector)
  dosevector=dosevector[!is.na(dosevector)]
  if(length(dosevector)<2) out=T
  if(length(dosevector[dosevector>0])<2) {
   # print("not at least two dosages available")
    out=T}
  else{
  minimum=min(dosevector[dosevector>0])
  difference=(max(dosevector)-minimum)/max(dosevector)
  if(difference>0.5) out=T
  else out=F
  }
  out
}



doseout=tapply(SchizoEFF$DoseOlaEquivalent,SchizoEFF$Study_No,unfair.dose.fun)
out=unique(as.numeric(names(doseout)))
out=out[doseout]
tapply(SchizoEFF$DoseOlaEquivalent,SchizoEFF$Study_No,print)

keep=is.na(match(SchizoEFF$Study_No,out)) 
SchizoEFFdose=SchizoEFF[keep,] 
length(table(SchizoEFFdose$Study_No))
range(table(SchizoEFFdose$Study_No))
SchizoEFFdose=SchizoEFFdose[order(SchizoEFFdose$Study_No),]

#network plot
SchizoEFFpairsdose=pairwise(treat=Drug,mean=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN, data=SchizoEFFdose, studlab = Study_No, sm= "SMD")
EFFdose=netmeta(SchizoEFFpairsdose, tol.multiarm=0.01)
#netrank(EFFdose)
#netleague(EFFdose)
#netgraph(EFFdose, plastic=F, thickness="number.of.studies", multiarm = F, points=T, cex.points=table(SchizoEFFdose$Drug)/4, col=1)

#transform the data into a list suitable for JAGS analysis
NMRdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFFdose,type="cont",reference = "Placebo")

#run Jags and create a jags object
NMRinJAGSdose<- jags.parallel(data = NMRdataContinuous, inits = NULL,
                    parameters.to.save = c("SMD.ref","tau", "SMD","SUCRA"), n.chains = 2, n.iter = 100000,
                    n.burnin = 5000,DIC=F,n.thin=10,
                    model.file = modelNMAContinuous)
save(NMRinJAGSdose,file="NMRinJAGSdose.RData",envir = .GlobalEnv)

#traceplot(NMRinJAGSdose,varname="tau" )
#traceplot(NMRinJAGSdose,varname="SMD.ref" )


leaguetable=out.jagsNMA.results(NMRinJAGSdose,"SMD",treatnames=sort(unique(SchizoEFFdose$Drug)))
leaguetableEFFdose=leaguetable$leaguetable
write.csv(leaguetableEFFdose,file="leaguetableEFFdose.csv")

#Present results
y=NMRinJAGSdose$BUGSoutput$mean$SMD.ref
ES=y
seES=NMRinJAGSdose$BUGSoutput$sd$SMD.ref
drugs=sort(unique(SchizoEFFdose$Drug))
forest(metagen(ES[order(y)],seES[order(y)],studlab=drugs[order(y)],comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,
       layout="JAMA",xlab="SMD from head-to-head studies only")

#SchizoEFFpairsdose[SchizoEFFpairsdose$treat2=="Ziprasidone",]
#SchizoEFFpairsdose[SchizoEFFpairsdose$treat2=="Zuclopenthixol",]
#SchizoEFFpairsdose[SchizoEFFpairsdose$treat2=="Brexpiprazole",]
#SchizoEFFpairsdose[SchizoEFFpairsdose$studlab==115,]

#*********************************************************************************
#        exclude studies before 1990  
#*********************************************************************************

SchizoEFF1990=SchizoEFF[SchizoEFF$Year>=1990,] 
length(table(SchizoEFF1990$Study_No))
range(table(SchizoEFF1990$Study_No))

#transform the data into a list suitable for JAGS analysis
NMRdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFF1990,type="cont",reference = "Placebo")

#run Jags and create a jags object
NMRinJAGS1990<- jags.parallel(data = NMRdataContinuous, inits = NULL,
                              parameters.to.save = c("SMD.ref","tau", "SMD","SUCRA"), n.chains = 2, n.iter = 100000,
                              n.burnin = 5000,DIC=F,n.thin=10,
                              model.file = modelNMAContinuous)
save(NMRinJAGS1990,file="NMRinJAGS1990.RData",envir = .GlobalEnv)

#traceplot(NMRinJAGS1990,varname="tau" )
#traceplot(NMRinJAGS1990,varname="SMD.ref" )


leaguetable=out.jagsNMA.results(NMRinJAGS1990,"SMD",treatnames=sort(unique(SchizoEFF1990$Drug)))
leaguetableEFF1990=leaguetable$leaguetable
write.csv(leaguetableEFF1990,file="leaguetableEFF1990.csv")

#Present results
y=NMRinJAGS1990$BUGSoutput$mean$SMD.ref
ES=y
seES=NMRinJAGS1990$BUGSoutput$sd$SMD.ref
drugs=sort(unique(SchizoEFF1990$Drug))
forest(metagen(ES[order(y)],seES[order(y)],studlab=drugs[order(y)],comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,
       layout="JAMA",xlab="SMD from studies in 1990 or later")

#*********************************************************************************
#        exclude failed studies  
#*********************************************************************************

SchizoEFFnofail=SchizoEFF[SchizoEFF$FailedStudy<1,] 
length(table(SchizoEFFnofail$Study_No))
range(table(SchizoEFFnofail$Study_No))

#transform the data into a list suitable for JAGS analysis
NMRdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFFnofail,type="cont",reference = "Placebo")

#run Jags and create a jags object
NMRinJAGSnofail<- jags.parallel(data = NMRdataContinuous, inits = NULL,
                              parameters.to.save = c("SMD.ref","tau", "SMD","SUCRA"), n.chains = 2, n.iter = 100000,
                              n.burnin = 5000,DIC=F,n.thin=10,
                              model.file = modelNMAContinuous)
save(NMRinJAGSnofail,file="NMRinJAGSnofail.RData",envir = .GlobalEnv)

#traceplot(NMRinJAGSnofail,varname="tau" )
#traceplot(NMRinJAGSnofail,varname="SMD.ref" )


leaguetable=out.jagsNMA.results(NMRinJAGSnofail,"SMD",treatnames=sort(unique(SchizoEFFnofail$Drug)))
leaguetableEFFnofail=leaguetable$leaguetable
write.csv(leaguetableEFFnofail,file="leaguetableEFFnofail.csv")

#Present results
y=NMRinJAGSnofail$BUGSoutput$mean$SMD.ref
ES=y
seES=NMRinJAGSnofail$BUGSoutput$sd$SMD.ref
drugs=sort(unique(SchizoEFFnofail$Drug))
forest(metagen(ES[order(y)],seES[order(y)],studlab=drugs[order(y)],comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,
       layout="JAMA",xlab="SMD excluding failed studies")
dev.off()
sink()
