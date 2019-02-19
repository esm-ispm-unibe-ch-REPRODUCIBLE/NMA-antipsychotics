
pdf("Analysis of primary outcome graphs.pdf")
sink("Analysis of primary outcome output.txt")

##################################################################################
#### ANALYSIS Primary outcome efficacy
##################################################################################



##get and see the data and exclude drugs and studies
#---------------------------------------------------


# Select those that report the primary outcome
EFFmissing=(is.na(SchizoOBJ$OverallEfficM)|is.na(SchizoOBJ$OverallEfficSD)|is.na(SchizoOBJ$OverallEfficN))
SchizoEFF=SchizoOBJ[EFFmissing==F,]

write.table(SchizoEFF,"SchizoWITHOUTopenlabelREPORTINGeffic.csv", sep="\t")
SchizoEFFCinema=cbind.data.frame(study=SchizoEFF$Study_No,id=SchizoEFF$Study_No,t=SchizoEFF$Drug, y=SchizoEFF$OverallEfficM,sd=SchizoEFF$OverallEfficSD,n=SchizoEFF$OverallEfficN, rob=SchizoEFF$Allocation+1)
write.table(SchizoEFFCinema,"SchizoEFFcinema.csv",  sep=";",row.names=F,quote = F)


#---------------------------------------------------
#describe the data

cat("\n \n FREQUENCY OF DRUGS \n \n")
print(tapply(SchizoEFF$OverallEfficN,SchizoEFF$Drug,sum))
cat("\n \n TOTAL NUMBER OF DRUGS \n \n")
print(length(unique(SchizoEFF$Drug)))

#------------------
##Describe the data 
#------------------
cat("\n \n TOTAL NUMBER OF STUDIES \n \n")
print(length(table(SchizoEFF$Study_No)))#nr of studies
cat("\n \n FREQUENCY OF MULTIARM STUDIES \n \n")
print(table(table(SchizoEFF$Study_No)))#type of studies with multiarms



cat("\n \n \n ***************************** \n \n \n")
cat("\n \n \n RESULTS USING NETMETA  \n \n \n")
cat("\n \n \n ***************************** \n \n \n")

#data transformation
SchizoEFFpairs=pairwise(treat=Drug,mean=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN, data=SchizoEFF, studlab = Study_No, sm= "SMD")
head(SchizoEFFpairs)

#network plot
EFF=netmeta(SchizoEFFpairs, tol.multiarm=0.01)
netgraph(EFF, plastic=F, thickness="number.of.studies", multiarm = F, points=T, cex.points=table(SchizoEFF$Drug)/12, col=1)


#------------------
##Perform NMA
#------------------
# create an object with the NMA data and output
EFF<-netmeta(TE,seTE,treat1,treat2,studlab,data=SchizoEFFpairs,  sm="SMD",r="Placebo",comb.fixed =F, comb.random = T, tol.multiarm=0.01)
print(EFF, digits=2, ref="Placebo")

#------------------
## Present Results
#------------------
#League table
netleague(EFF, digits=2)


#forest plots
#forest(EFF, ref="Placebo", sortvar = EFF$TE.random[,23],xlab="SMD")

#treatment ranking
cat("\n  P-SCORES  \n")
netrank(EFF)
league0 <- netleague(EFF, digits = 2, bracket = "(", separator = " to ")


#------------------
#Testing for inconsistency
#------------------
cat("\n  INCONSISTENCY  \n")
split=netsplit(EFF) 
names(split)
SIDEp=split$compare.random$p
SIDEp=SIDEp[!is.na(SIDEp)]
#proportion of inconsistent loops
length(SIDEp[SIDEp<=0.1])/length(SIDEp)
print(split,showall = F,digits=2)
print(decomp.design(EFF))



#*********************************************************************************
#             NMA USING JAGS              
#*********************************************************************************


cat("\n \n \n JAGS output \n \n \n")


SchizoEFF=SchizoEFF[order(SchizoEFF$Study_No),]

#transform the data into a list suitable for JAGS analysis
NMAdataContinuous=make.jagsNMA.data(studyid=Study_No,t=Drug,y=OverallEfficM,sd=OverallEfficSD,n=OverallEfficN,data=SchizoEFF,type="cont",reference = "Placebo")
#run Jags and create a jags object
NMAinJAGS<- jags.parallel(data = NMAdataContinuous, inits = NULL,
                 parameters.to.save = c("SMD","SMD.ref","tau", "SUCRA"), n.chains = 2, n.iter = 100000,
                 n.burnin = 5000,DIC=F,n.thin=10,
                 model.file = modelNMAContinuous)
print(NMAinJAGS)
save(NMAinJAGS,file="NMAinJAGSall.RData",envir = .GlobalEnv)

#check chain mixing
traceplot(NMAinJAGS,varname="tau" )
#traceplot(NMAinJAGS,varname="SMD.ref" )


#forestplot against placebo
y=NMAinJAGS$BUGSoutput$mean$SMD.ref
ES=y
seES=NMAinJAGS$BUGSoutput$sd$SMD.ref
drugs=sort(unique(SchizoEFF$Drug))[order(y)]
forest(metagen(ES[order(y)],seTE=seES[order(y)],studlab=drugs,comb.fixed = F, comb.random = F), xlim=c(-2,1), digits=3,studlab = drugs,layout="JAMA",xlab="SMD meta-regressions for 6 change in PANSS for placebo")


#then get the league table and save it
leaguetable=out.jagsNMA.results(NMAinJAGS,"SMD",F, treatnames=sort(unique(SchizoEFF$Drug)))
leaguetableEFF=leaguetable$leaguetable
write.csv(leaguetableEFF,file="leaguetableEFF.csv")

###Compare JAGS with netmeta
NMAinJAGS
cbind(EFF$TE.random[,23],NMAinJAGS$BUGSoutput$mean$SMD.ref)
plot(EFF$TE.random[,23],NMAinJAGS$BUGSoutput$mean$SMD.ref,xlab="SMD vs Placebo netmeta",ylab="SMD vs Placebo NMAJags")

cat("\n  SUCRAS  \n")
print(NMAinJAGS$BUGSoutput$mean$SUCRA)
cat("\n  COMMON HETEROGENEITY  \n")
print(NMAinJAGS$BUGSoutput$mean$tau)



dev.off()
sink()
