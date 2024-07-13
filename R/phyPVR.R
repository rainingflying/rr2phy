
phyPVR<- function(x, phy = NULL, trait = NULL, envVar = NULL, method = "Moran.I",
                  weights = NULL, scaled = FALSE, sig = TRUE, sig.t = 0.05, lamdba.t = 0.05, psr.t = 0.01, accvalue.t = 0.9, ...)
{
  
  if(method == "Moran.I" ){
    
    pvr <- x@Eigen$vectors
    Naxis <- ncol(pvr)
    if(is.null(weights)){
      
      W <- max(x@phyDist) - x@phyDist
    } else {
      W <- weights
    }
    diag(W) <- 0
    
    N <- .poss(pvr)
    c1 <- 1
    c2 <- 1
    
    tmpRes <- matrix(nrow = Naxis, ncol = 2)
    selected <- matrix(nrow = nrow(pvr), ncol = Naxis)
    reg <- data.frame(Vectors = 0, Moran.I = 0, p = 0)
    tmpTrait <- trait
    
    for(i in 1:N){
      
      pvr <- as.matrix(pvr)
      tmpLM <- lm(tmpTrait ~ pvr[ ,c1])
      #MI <-lambdaTest1(tmpLM$residuals,x@Eigen$vectors)
      MI <- Moran.I(tmpLM$residuals, weight = W, scaled = scaled, alternative = "two.sided")
      #tmpRes[c1, 1] <- MI$p[1]
      #tmpRes[c1, 2] <- MI$Lambda
      tmpRes[c1, 1] <- round(MI$p.value, 5)
      tmpRes[c1, 2] <- MI$observed
      
      if(c1 == ncol(pvr)){
        
        selAxis <- which(tmpRes[ ,2] == max(tmpRes[ ,2]))
        tmpLM <- lm(tmpTrait ~ pvr[ ,selAxis[1]])
        tmpTrait <- tmpLM$residuals
        selected[ ,c2] <- pvr[ ,selAxis[1]]
        reg[c2,1] <- colnames(pvr)[selAxis[1]]
        reg[c2,2] <- tmpRes[selAxis[1], 2]
        reg[c2,3] <- tmpRes[selAxis[1], 1]
        tmpColNames <- colnames(pvr)
        pvr <- pvr[ ,-selAxis[1]]
        tmpColNames <- tmpColNames[-selAxis[1]]
        pvr <- as.matrix(pvr)
        colnames(pvr) <- tmpColNames
        if(sig){
          
          if(tmpRes[selAxis[1], 1] > sig.t | c2 == Naxis){
            
            break
          }
        } else {
          
          if(tmpRes[selAxis[1], 2] > MI.t | c2 == Naxis){
            
            break
          }
        }
        
        
        c2 <- c2 + 1
        c1 <- 0
        
        tmpRes <- matrix(nrow = ncol(pvr), ncol = 2)
      }
      
      c1 <- c1 + 1
    }
    
    selected <- as.matrix(selected[ ,1:nrow(reg)])
    colnames(selected) <- reg$Vectors
    
    selection <- list(Method = "Moran.I", Vectors = selected, Id = reg)
    x@Selection <- selection
    print(x@Selection )
  }
  if(!is.null(envVar)){
    traits=data.frame(trait)
    iv <- list(
      envVar=data.frame(envVar),
      phy=data.frame(selection$Vectors))
    hp <- rdacca.hp(traits,iv, method = 'RDA', type = 'R2',var.part = TRUE)
    print(hp)
    plot(hp)
    
  } else{
    
    pvrOLS <- lm(trait ~ selection$Vectors + envVar)
    x@PVR <- list(R2 = summary(pvrOLS)$r.squared, Residuals = pvrOLS$residuals, p = (summary(pvrOLS))$coefficient[2, 4])
  }
}

.poss <- function(x){
  
  N <- 0
  for(i in 1:ncol(x)){
    
    N <- N + i
  }
  return(N)
}

#lambdaTest1 <- function(x, vcv){
 # lambda.max <- max(vcv)/max(vcv[lower.tri(vcv)])
 # opt <- suppressWarnings(
 #   optimize(pagelLogLik, interval=c(0,0.1, lambda.max), xr = x, vcvr = vcv, maximum = TRUE)
#  )
#  Lambda <- opt$maximum
#  logL0 <- pagelLogLik(0, x, vcv)
#  pvalue <- as.numeric(pchisq(2 * (opt$objective - logL0), df = 1, lower.tail = FALSE))
  
#  return(list(Lambda = Lambda, pvalue = pvalue))
#}


