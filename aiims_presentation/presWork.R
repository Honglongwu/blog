library(tidyr, quietly = TRUE)
library(magrittr, quietly = TRUE)
library(purrr, quietly = TRUE)
library(metafor, quietly = TRUE)

set.seed(seed)

simBin <- function(pa, pb, nSamp, nSim) {
  
  sampA <- replicate(nSim, sample(c(0, 1), replace = TRUE, size = nSamp, prob = c(1 - pa, pa)), simplify = FALSE)
  sampB <- replicate(nSim, sample(c(0, 1), replace = TRUE, size = nSamp, prob = c(1 - pb, pb)), simplify = FALSE)
  #names(sampA) <- paste0("sim", 1:nSim)
  #names(sampB) <- paste0("sim", 1:nSim)
  
  prop <- function(x) sum(x)/length(x)
  propA <- sampA %>% map_v(prop)
  propB <- sampB %>% map_v(prop)
  
  succ <- function(x) sum(x)
  succA <- sampA %>% map_v(succ)
  succB <- sampB %>% map_v(succ)
  
  ci <- function(y, x) {
    r <- prop.test(c(sum(x), sum(y)), c(length(x), length(y)))
    return(data_frame(lcl = r$conf.int[1], ucl = r$conf.int[2], pVal = r$p.value))
  }
  
  dfCI <- map2(sampA, sampB, ci) %>% bind_rows()
  
  df <- data_frame(a = propA, b = propB, succA = succA, succB = succB, diff = b - a, lclDiff = dfCI$lcl, uclDiff = dfCI$ucl, pValDiff = dfCI$pVal, incl0 = inclX(0, lclDiff, uclDiff), inclpDiff = inclX((pb - pa), lclDiff, uclDiff))
  
  pow <-  sum(1*!(df$incl0))/nrow(df)
  alphaErr <- sum(1*!(df$inclpDiff))/nrow(df)
  
  res <- structure(list(data = df, input = list(pa = pa, pb = pb, pDiff = (pb - pa), power = pow, alphaErr = alphaErr, nSamp = nSamp)), class = "simBin")
  return(res)
  
}

plot.simBin <- function(x) {
#-------------------------------------------  
  plot1 <- function(n = 50) {
    stopifnot(nrow(x$data) > n)
    dat <- x$data
  
    cl <- dat %>% summarise(lcl = quantile(diff, 0.025), ucl = quantile(diff, 0.975))
  
    dat1 <- dat %>% select(diff) %>% slice(1:n) %>% mutate(simNo = factor(1:n), z = factor(1*(diff <= 0), c(0,1), c("B > A", "A > B")))
  
    p1 <- ggplot(data = dat1) + geom_point(aes(x = simNo, y = diff, color = z), size = 4) + geom_hline(yintercept = c(x$input$pDiff, cl$lcl, cl$ucl), linetype = c("dashed", "solid", "solid"), colour = c("red", "blue", "blue")) + scale_y_continuous(limits = c(-0.4, 0.8), breaks = seq(-0.4, 0.8, 0.1)) + scale_color_manual(values = c("black", "red"), name = "") + coord_flip() + theme_bw()
  
    dat2 <- dat %>% select(diff)
  
    p2 <- ggplot(data = dat2) + geom_histogram(aes(x = diff, y = ..density..)) + geom_vline(xintercept = c(0, x$input$pDiff, cl$lcl, cl$ucl), linetype = c("dashed", "dashed", "solid", "solid"), color = c("red", "black", "blue", "blue")) + scale_x_continuous(limits = c(-0.4, 0.8), breaks = seq(-0.4, 0.8, 0.1)) + theme_bw()
  
    grid.arrange(p1, p2, nrow = 1)
    
  }
#-----------------------------------------------------------
  plot2 <- function(n = 50) {
    
    stopifnot(nrow(x$data) > n)
    dat <- x$data
    
    dat1 <- dat %>% select(diff, lclDiff, uclDiff, incl0, inclpDiff) %>% slice(1:n) %>% mutate(simNo = factor(1:n), incl0 = factor(1*incl0), inclpDiff = factor(1*inclpDiff))
    
    ggplot(data = dat1) + geom_point(aes(x = simNo, y = diff, shape = inclpDiff), size = 4) + scale_shape_manual(values = c(1, 19), limits = c(0, 1), breaks = c(0, 1), name = "pDiff included", labels = c("No", "Yes")) + geom_errorbar(aes(x = simNo, ymin = lclDiff, ymax = uclDiff, color = incl0)) + scale_color_manual(values = c("red", "black"), limits = c(0, 1), breaks = c(0, 1), name = "Zero included", labels = c("No", "Yes")) + geom_hline(yintercept = c(0, x$input$pDiff), linetype = "dashed", colour = c("red", "blue")) + scale_y_continuous(limits = c(-0.8, 0.8), breaks = seq(-0.8, 0.8, 0.1)) + coord_flip() + theme_bw()
  }
#------------------------------------------------------------
  function(choice = c("withoutCI", "withCI")) {
    
    ch <- match.arg(choice)
    
    if (ch == "withoutCI") {
      plot1
    } else if (ch == "withCI") {
      plot2
    }
  }
  
}

simMA <- function(vSampSize, nTrials, pa, pb) {
  
  vs <- vSampSize
  succA <- succB <- integer(nTrials)
  nA <- nB <- integer(nTrials)
  
  sSize <- sample(vs, size = nTrials, replace = TRUE)
  
  for (i in 1:nTrials) {
    ss <- simBin(pa, pb, sSize[i], 1)
    succA[i] <- ss$data$succA
    succB[i] <- ss$data$succB
    nA[i] <- sSize[i]
    nB[i] <- sSize[i]
  }
  
  return(data_frame(succA = succA, succB = succB, nA = nA, nB = nB))
  
}
