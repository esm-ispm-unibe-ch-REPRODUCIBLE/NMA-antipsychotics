modelNMRPlacChange=function()
{
  
  ####################################################################################
  #### Regression model to predict placebo change in head to head studies using the year
  
  #Data: r (length as ns with a 0 at head to head studies), n (length as ns with a 0 at head to head studies )
  # ns total number of studies
  # nsP number of studies with Placebo arms 
  # nsHtH total number of HtH studies 
  #Pstudies a vector of lenght nsP that has the ids for Placebo-controlled studies
  #HtHstudies a vector of lenght nsHtH that has the ids for HtH studies
  #year i the year centralised (2016-year)
  #pchange the placebo change
  #nsP the number of placebo-controlled studies
  
  for (i in 1:nsP)#likelihood for the placebo-arms
  {
    pchange[i]~dnorm(PCHANGE[i],precpchange)
    PCHANGE[i]<-a+breg*year[i]
    npchange[i]<-pchange[i]
  }
  
  for (i in (nsP+1):ns) #predictions for the head-to-head studies (but includes all studies)
  {
    npchange[i]~dnorm(nPCHANGE[i],precpchange)
    nPCHANGE[i]<-a+breg*year[i]
  }
  
  #Priors
  a~dnorm(0,0.001)
  breg~dnorm(0,0.001)
  precpchange~dgamma(0.001,0.001)
  varpchange<-1/precpchange
  taupchange<-sqrt(varpchange)
  
  
  ###### NMR model for the real and predicted placebo change npchange ####
  
  for (i in 1:ns) {
    
    w[i, 1] <- 0
    delta[i, t[i, 1]] <- 0
    u[i] ~ dnorm(0, 1e-04)
    for (k in 1:na[i]) {
      y[i, t[i, k]] ~ dnorm(phi[i, t[i, k]], prec[i, t[i, 
                                                       k]])
      phi[i, t[i, k]] <- (u[i] + delta1[i, t[i, k]]) * 
        pooled.sd[i]
      delta1[i, t[i, k]] <- delta[i, t[i, k]] + beta[t[i, 
                                                       1], t[i, k]] * (npchange[i]-refvaluecov)
    }
    for (k in 2:na[i]) {
      delta[i, t[i, k]] ~ dnorm(md[i, t[i, k]], taud[i, 
                                                     t[i, k]])
      md[i, t[i, k]] <- d[t[i, k]] - d[t[i, 1]] + sw[i, 
                                                     k]
      taud[i, t[i, k]] <- PREC * 2 * (k - 1)/k
      w[i, k] <- (delta[i, t[i, k]] - d[t[i, k]] + d[t[i, 
                                                       1]])
      sw[i, k] <- sum(w[i, 1:(k - 1)])/(k - 1)
    }
  }
  d[ref] <- 0
  for (k in 1:(ref - 1)) {
    d[k] ~ dnorm(0, 1e-04)
  }
  for (k in (ref + 1):nt) {
    d[k] ~ dnorm(0, 1e-04)
  }
  tau ~ dunif(0, 5)
  PREC <- 1/pow(tau, 2)
  for (c in 1:(nt - 1)) {
    for (k in (c + 1):nt) {
      SMD[c, k] <- d[c] - d[k]
    }
  }
  for (c in 1:nt) {
    SMD.ref[c] <- d[c] - d[ref]
  }
  for (c in 1:(ref - 1)) {
    X[c] <- d[c] - d[ref]
    predSMD.ref[c] ~ dnorm(X[c], PREC)
  }
  for (c in (ref + 1):nt) {
    X[c] <- d[c] - d[ref]
    predSMD.ref[c] ~ dnorm(X[c], PREC)
  }
  for (c in 1:(nt - 1)) {
    for (k in (c + 1):nt) {
      predSMD[c, k] ~ dnorm(SMD[c, k], PREC)
    }
  }
  order[1:nt] <- rank(d[1:nt])
  for (k in 1:nt) {
    most.effective[k] <- equals(order[k], 1)
    for (j in 1:nt) {
      effectiveness[k, j] <- equals(order[k], j)
    }
  }
  for (k in 1:nt) {
    for (j in 1:nt) {
      cumeffectiveness[k, j] <- sum(effectiveness[k, 1:j])
    }
  }
  for (k in 1:nt) {
    SUCRA[k] <- sum(cumeffectiveness[k, 1:(nt - 1)])/(nt - 
                                                        1)
  }
  
  
  b[ref] <- 0
  for (k in 1:(ref - 1)) {
    b[k] ~ dnorm(0, 1e-04)
  }
  for (k in (ref + 1):nt) {
    b[k] ~ dnorm(0, 1e-04)
  }
  
  for (i in 1:nt) {
    for (j in 1:nt) {
      beta[i, j] <- b[j] - b[i]
    }
  }
  
  
  
  
}
