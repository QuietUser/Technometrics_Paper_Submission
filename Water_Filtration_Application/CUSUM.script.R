load("../Data/headloss_matrix.rda")
load("../Data/headloss_cycles_info.rda")
library(fields)
library(qcc)


cusum.with.reset <- function(data,center=0,decision.interval=4,se.shift=.5){
  data <- data - center
  x.n <- length(data)
  # sd.dd <- sd(data)
  sd.dd <- sd.xbar.one(data)
  # sd.dd <- cusum(detrended.data,decision.interval = 4,se.shift=.5)$std.dev
  slack <- se.shift
  limit <- decision.interval
  cs.pos <- 0
  cs.neg <- 0
  cs.pos.vec <- rep(0,x.n)
  cs.neg.vec <- rep(0,x.n)
  for ( i in 1:x.n){
    cs.pos <- max(cs.pos + data[i]/sd.dd - slack,0)
    cs.neg <- min(cs.neg + data[i]/sd.dd + slack,0)
    cs.pos.vec[i] <- cs.pos
    cs.neg.vec[i] <- cs.neg
    
    if (cs.pos > limit | cs.neg < -limit){
      cs.pos <- 0
      cs.neg <- 0
    }
    
    violations <- sort(c(which(cs.pos.vec>limit),which(cs.neg.vec< -limit)))
  }
  violations
}

all.analyzed.cycles.list.cusum <- vector(mode = "list", length = 3305)
all.cycles.anomaly.table.cusum <- data.frame(cycle.number=1:3305,cycle.length=headloss_cycles_info$cycle_length-6,Anomalous=NA,Identified=NA,s=NA,max.gamma=NA,run.time=NA)

for (cn in 1:3305){
  if (headloss_cycles_info$cycle_length[cn] > 60){
    n <- headloss_cycles_info$cycle_length[cn]+1
    x <- 1:(n-6)
    y <- headloss_matrix[cn,4:(n-3)]
    n <- length(x)
    
    init.time <- proc.time()
    ### Detrend
    detrended.data <- Tps(x,y)$residuals
    ### CUSUM
    violations <- cusum.with.reset(detrended.data,center=0,decision.interval=4,se.shift=.5)
    time <- time <- (proc.time()-init.time)[3]
    
    all.analyzed.cycles.list.cusum[[cn]] <- violations
    names(all.analyzed.cycles.list.cusum)[cn] <- paste0("model.",cn)
    all.cycles.anomaly.table.cusum[cn,4] <- ifelse(length(violations)!=0,"YES","NO")
    # all.cycles.anomaly.table.cusum[cn,5] <- model$s
    # all.cycles.anomaly.table.cusum[cn,6] <- ifelse(length(violations)!=0, max(violations),0)
    all.cycles.anomaly.table.cusum[cn,7] <- time
  }
  
  if (headloss_cycles_info$cycle_length[cn] < 60){
    all.analyzed.cycles.list.cusum[[cn]] <- "Short Cycle"
  }
  
  if(cn%%20==0) cat(cn," ")
  
}

save(all.cycles.anomaly.table.cusum,file = "UWTF.CUSUM.results.Rda")