#-----PARAMETERS-----------

param.s <- c(sc = 0.01,     # S to C
             sd = 0.01,     # Death
             ss = 0.02)     # S to S
param.c <- c(cp = 0.1,     # C to P
             cc = 0.05,     # C to C
             h = 0.1)      # Controls Deaths in C comp
param.p <- c(pp = 0.1,     # P to P
             pd = 0.05,     # Death
             tau.p = 5)  # time after which P differentiates into R
param.r <- c(tau.r = 100,  # time after which R dies
             rh = 1e-3)     # Factor for converting R to H
tf <- 500
PARAMETERS <- c(param.s,param.c,param.p,param.r,tf)

#----STATE VARIABLES----------

yini <- c(S = 1e5,  # Stem cell/BFU compartment
          C = 0,  # CFU compartment
          P = 0,  # Precursor compartment
          R = 0,  # RBC compartment
          newP = 0,
          newR = 0)
#------WORK HORSE FUNCTION-----------

foo <- function(t, y, param,...) {
  with(as.list(c(y, param)), {
    
    #if (t >= 61) browser()
    dS <- (ss - sc - sd) * S
    cd <- rh * h * R    # Determinant of death in C compartment.
    dC <- (sc * S) + (cc - cp - cd) * C
    dnewP <- (cp * C) - newP  # new P entrant at time t
    #if (t > 0) browser()
    tlagP <- t - as.integer(seq(tau.p,to=0))
    lagnewP <- rep(0,length.out=length(tlagP))
    for (i in 1:length(tlagP)) {
      if (tlagP[i] <= 0) {
        lagnewP[i] <- 0
      } else {
        lagnewP[i] <- lagvalue(tlagP[i],5)
      }
    }
    sum.int <- sum(lagnewP*exp(pp*(tau.p:0)))
    sum.int <- ifelse(sum.int == 0, 1, sum.int)
    pr <- (lagnewP[1]*exp(pp*tau.p))/sum.int
    dP <- newP + (pp - pd - pr) * P
    tlagR <- t - as.integer(tau.r)
    dnewR <- (pr * P) - newR
    #if (lagvalue(tlagR,4) > 0) browser()
    if (tlagR > 0) {
      lagnewR <- lagvalue(tlagR,6)
    } else {
      lagnewR <- 0
    }
    #if (lagR > 0) browser()
    dR <- newR - lagnewR

    return(list(c(dS, dC, dP, dR, dnewP, dnewR), Hb = rh*R))
  })
}

#-------EVENT DATA----------------
eventdat <- data.frame(var = c('S'), time = c(200),value = c(0.5), method = c("mult"))
eventdat2 <- data.frame(var = c('S','S'), time = c(200,300),value = c(0.5,2), method = c("mult","mult"))


#-------DIFFERENTIAL EQUATION--------------

TIME <- seq(from  = 0, to = tf, by = 100)
system.time(sol <- dede(y = yini, times = TIME, func = foo, parms = PARAMETERS,events=list(data=eventdat)))

system.time(sol1 <- dede(y = yini, times = TIME, func = foo, parms = PARAMETERS))

system.time(sol2 <- dede(y = yini, times = TIME, func = foo, parms = PARAMETERS,events=list(data=eventdat2)))

#------DIAGNOSTICS------------
diagnostics(sol)
#------PLOTTINGS--------------
plot(sol1)
plot(sol,sol1,sol2)
par(mfrow=c(1,1))
plot(sol[,1],sol[,'S'],type='l')
matplot(sol[,1],sol[,-1],type='l')

matplot(sol[,1],cbind(sol[,'R'],sol1[,'R'],sol2[,'R']),type='l',col=c('black','red','green'),lty=1:3)