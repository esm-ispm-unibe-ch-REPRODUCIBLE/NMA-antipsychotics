
pdf("Meta-regression on placebo primary outcome graphs.pdf")
sink("Meta-regression on placebo primary outcome output.txt")

#*********************************************************************************
#                 Placebo response              
#*********************************************************************************

# Select those that report the primary outcome
SchizoEFF=SchizoEFF[order(SchizoEFF$Study_No),]


######exclude drugs with less than 100 patients randomised!
##################################################################################
while(min(tapply(SchizoEFF$OverallEfficN,SchizoEFF$Drug,sum))<=100){
  a=tapply(SchizoEFF$OverallEfficN,SchizoEFF$Drug,sum)
  indrugs=names(a)[a>100]
  SchizoEFF=subset(SchizoEFF,SchizoEFF$Drug%in%indrugs)
  instudy=names(table(SchizoEFF$Study_ID))[table(SchizoEFF$Study_ID)>1]#Take care to exclude studies left with one arm!
  SchizoEFF=subset(SchizoEFF,SchizoEFF$Study_ID%in%instudy)}

write.table(SchizoEFF,"~/Google Drive/_mydrive/Schizophrenia/SchizoR/SchizoWITHOUTopenlabelREPORTINGefficMOREthan100.csv",  sep="\t")

#-------------------------------------
##get and see the data  - and fix them
#-------------------------------------
# start with SchizoEFF
SchizoEFFyear=SchizoEFF[!is.na(SchizoEFF$Year),] 
SchizoEFFyear$Year1=2016-SchizoEFF$Year


cat(" \n MEAN YEAR:")
mean(SchizoEFFyear$PlaceboChange,na.rm=T)

#substitute 2 missing values in the placebo arm with the mean
SchizoEFFyear$PlaceboChange[is.na(SchizoEFF$PlaceboChange) & (SchizoEFF$Drug=="Placebo")]=mean(SchizoEFF$PlaceboChange,na.rm=T)
a=lm(SchizoEFFyear$PlaceboChange~c(SchizoEFFyear$Year))

plot(SchizoEFFyear$Year,SchizoEFFyear$PlaceboChange,xlab="Year of publication",ylab="Change in PANSS for the placebo arm")
abline(a,col="red",lwd=1.5)

#####---------------------
##    Make the data
#####------------------

#put placebo trials first
SchizoEFFyear=SchizoEFFyear[order(-SchizoEFFyear$Placebo),] 
a=SchizoEFFyear$Study_No
id=rep(1:length(table(a)),table(match(a,unique(a))))
SchizoEFFyear$Study_No=id

#create the PlaceboChange at the study level and the data
pchange=tapply(SchizoEFFyear$PlaceboChange,SchizoEFFyear$Study_No,mean,na.rm=T)
year=tapply(SchizoEFFyear$Year,SchizoEFFyear$Study_No,mean,na.rm=T)
pchange[is.na(pchange)]=-0.0001
ns=length(pchange)
nsP=sum(pchange!=-0.0001)

###############################################


NMRdataPlacChange=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFFyear,type="cont",reference = "Placebo",othervar=2016-Year)
NMRdataPlacChange$nsP=nsP
NMRdataPlacChange$pchange=pchange[1:nsP]
NMRdataPlacChange$year=2016-year
NMRdataPlacChange$variab=NULL
NMRdataPlacChange$refvaluecov=-6
names(NMRdataPlacChange)


NMRPlacChange<- jags.parallel(data =NMRdataPlacChange, inits = NULL,
                            parameters.to.save = c("a","b", "SMD", "SMD.ref","tau","breg","SUCRA"), n.chains = 2, n.iter = 100000,
                            n.burnin = 5000,DIC=F,n.thin=10,
                            model.file = modelNMRPlacChange)
save(NMRPlacChange,file="NMRPlacChangereference6unitspchange.RData",envir = .GlobalEnv)

traceplot(NMRPlacChange,varname="tau" )
#traceplot(NMRPlacChange,varname="SMD.ref" )
#traceplot(NMRPlacChange,varname="b" )

SUCRA2=NMRPlacChange$BUGSoutput$mean$SUCRA

#forest plots
y=NMRPlacChange$BUGSoutput$mean$SMD.ref
ES=y[order(y)]
seES=NMRPlacChange$BUGSoutput$sd$SMD.ref[order(y)]
drugs=sort(unique(SchizoEFF$Drug))[order(y)]
forest(metagen(ES,seTE=seES,studlab=drugs,comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,studlab = drugs,layout="JAMA",xlab="SMD meta-regressions for 6 change in PANSS for placebo")

#then get the league table and save it
leaguetable=out.jagsNMA.results(NMRPlacChange,"SMD",F, treatnames=sort(unique(SchizoEFF$Drug)))
leaguetableEFFPlacChange6=leaguetable$leaguetable
write.csv(leaguetableEFFPlacChange6,file="leaguetableEFFPlacChange6.csv")

coef=NMRPlacChange$BUGSoutput$mean$b[order(y)]
secoeff=NMRPlacChange$BUGSoutput$sd$b[order(y)]

forest(metagen(coef,secoeff,comb.fixed = F),digits=3,studlab = drugs,layout="JAMA", xlab="regression coefficients for SMD ~ PANSS for placebo")

#fit a model at placebo change 0
NMRdataPlacChange$refvaluecov=0
names(NMRdataPlacChange)

NMRPlacChange<- jags.parallel(data =NMRdataPlacChange, inits = NULL,
                     parameters.to.save = c("a","b", "SMD", "SMD.ref","tau","breg","SUCRA"), n.chains = 2, n.iter = 100000,
                     n.burnin = 5000,DIC=F,n.thin=10,
                     model.file = modelNMRPlacChange)
save(NMRPlacChange,file="NMRPlacChangereference0unitspchange.RData",envir = .GlobalEnv)
leaguetable=out.jagsNMA.results(NMRPlacChange,"SMD",F, treatnames=sort(unique(SchizoEFF$Drug)))
leaguetableEFFPlacChange0=leaguetable$leaguetable
write.csv(leaguetableEFFPlacChange0,file="leaguetableEFFPlacChange0.csv")

#####################################################################
###  FIT PLACEBO RESPONSE MODEL WITH EXCHANGEABLE COEFFICIENTS
#####################################################################

NMRPlacChangeExchB<- jags.parallel(data =NMRdataPlacChange, inits = NULL,
                     parameters.to.save = c("B","tauB", "SMD", "SMD.ref","tau","SUCRA"), n.chains = 2, n.iter = 100000,
                     n.burnin = 5000,DIC=F,n.thin=10,
                     model.file = modelNMRPlacChangeExchB)
save(NMRPlacChangeExchB,file="NMRPlacChangereference6unitspchangeExchB.RData",envir = .GlobalEnv)


traceplot(NMRPlacChangeExchB,varname="tau" )
#traceplot(NMRPlacChangeExchB,varname="tauB" )
traceplot(NMRPlacChangeExchB,varname="B" )



#forest plots
y=NMRPlacChangeExchB$BUGSoutput$mean$SMD.ref
ES=y[order(y)]
seES=NMRPlacChangeExchB$BUGSoutput$sd$SMD.ref[order(y)]
drugs=sort(unique(SchizoEFF$Drug))[order(y)]
forest(metagen(ES,seTE=seES,studlab=drugs,comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,studlab = drugs,layout="JAMA",xlab="SMD meta-regressions for 6 change in PANSS for placebo with exchangeable coefficients")

#then get the league table and save it
leaguetable=out.jagsNMA.results(NMRPlacChangeExchB,"SMD",F, treatnames=sort(unique(SchizoEFF$Drug)))
leaguetableEFFPlacChange6ExchB=leaguetable$leaguetable
write.csv(leaguetableEFFPlacChange6ExchB,file="leaguetableEFFPlacChange6ExchB.csv")

cat("\n \n \n MEAN COEFFICIENT BETA IN EXCHANGEABLE MODEL: \n \n \n")
MeanBeta=NMRPlacChangeExchB$BUGSoutput$summary[1,c(1,3,7)]
print(MeanBeta)

dev.off()
sink()
