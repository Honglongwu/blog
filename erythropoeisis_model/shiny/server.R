library(shiny)
library(deSolve)
# #-------EVENT DATA----------------
# #eventdat <- data.frame(var = c('S','S'), time = c(200,300),value = c(0.5,2), method = c("mult","mult"))
#print('enetered server.R')
# Workhorse function for calculating differential equations
foo <- function(t, y, param,...) {
  with(as.list(c(y, param)), {
    #if (t >= 500) browser()
    #if (t >= 61) browser()
    dS <- (ss - sc - sd) * S
    #dcd <- (rh * h * R * (hmax - (rh*R))) - cd   # Determinant of death in C compartment.
    cd <- rh * h * R
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
    
    return(list(c(dS, dC, dP, dR, dnewP, dnewR), Hb = rh*R))  #, dcd
  })
}

shinyServer(function(input,output) {
  
  PARAMETERS <- reactive({
    sc <- input$sc
    sd <- input$sd
    ss <- input$ss
    cp <- input$cp
    cc <- input$cc
    h <- input$h
    pp <- input$pp
    pd <- input$pd
    tau.p <- input$tau.p
    tau.r <- input$tau.r
    rh <- input$rh
    p <- c(sc=sc,sd=sd,ss=ss,cp=cp,cc=cc,h=h,pp=pp,pd=pd,tau.p=tau.p,tau.r=tau.r,rh=rh)
    return(p)
  })
  
  tf <- reactive({
    tf <- input$tf
    return(c(tf=tf))
  })
  
  dt <- reactive({
    dt <- input$dt
    return(c(dt=dt))
  })
  
  yini <- reactive({
    S <- input$S
    C <- input$C
    P <- input$P
    R <- input$R
    newP <- input$newP
    newR <- input$newR
    #cd <- input$cd
    y <- c(S=S,C=C,P=P,R=R,newP=newP,newR=newR) #,cd=cd)
    return(y)
  })

  t.print <- reactive({
    t.vec <- c(0,as.numeric(unlist(strsplit(input$dPrint,','))),tf())
    return(t.vec)
  })
  
  df.events <- reactive({
    var.vec <- str_trim(as.character(unlist(strsplit(input$vars,','))))
    time.vec <- as.numeric(unlist(strsplit(input$time,',')))
    val.vec <- as.numeric(unlist(strsplit(input$val,',')))
    method.vec <- str_trim(as.character(unlist(strsplit(input$method,','))))
    return(data.frame(var=var.vec,time=time.vec,value=val.vec,method=method.vec))
  })
  
  sol <- reactive({
    
    TIME <- seq(from  = 0, to = tf(), by = dt())
    if (nrow(df.events()) == 0) {
      sol <- dede(y = yini(), times = TIME, func = foo, parms = PARAMETERS())
      return(sol)  
    } else {
      sol <- dede(y = yini(), times = TIME, func = foo, parms = PARAMETERS(),events=list(data=df.events()))
      return(sol)
    }
    
  })
  
  return.lst <- reactive({
    state.var <- sol()[sol()[,1] %in% t.print(),]
    ret.val <- list(state.var=state.var,initParam=PARAMETERS())
    return(ret.val)
  })
     
  output$mPlot <- renderPlot({
    plot(sol())
  })
  
  output$initParam_stateVar <- renderPrint(return.lst())
  
  output$downloadParam <- downloadHandler(filename=function() {
    paste('data-',Sys.time(),'.RData',sep='')
  }, content=function(file) {
    save(return.lst,file=file)
  })

})

#print('exited server.R')