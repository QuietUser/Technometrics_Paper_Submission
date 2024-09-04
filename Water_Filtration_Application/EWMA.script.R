load("../Data/headloss_matrix.rda")
load("../Data/headloss_cycles_info.rda")
library(fields)
library(qcc)

all.analyzed.cycles.list.ewma <- vector(mode = "list", length = 3305)
all.cycles.anomaly.table.ewma <- data.frame(cycle.number=1:3305,cycle.length=headloss_cycles_info$cycle_length-6,Anomalous=NA,Identified=NA,s=NA,max.gamma=NA,run.time=NA)

for (cn in 1:3305){
  if (headloss_cycles_info$cycle_length[cn] > 60){
    n <- headloss_cycles_info$cycle_length[cn]+1
    x <- 1:(n-6)
    y <- headloss_matrix[cn,4:(n-3)]
    n <- length(x)
    
    init.time <- proc.time()
    ### Detrend
    detrended.data <- Tps(x,y)$residuals
    ### EWMA
    model <- ewma(detrended.data,center=0,plot=FALSE,lambda=.3)
    violations <- unique(c(unique(model$violations), unique(model$violations) + 1, unique(model$violations) - 1)) 
    time <- time <- (proc.time()-init.time)[3]
    
    all.analyzed.cycles.list.ewma[[cn]] <- violations
    names(all.analyzed.cycles.list.ewma)[cn] <- paste0("model.",cn)
    all.cycles.anomaly.table.ewma[cn,4] <- ifelse(length(violations)!=0,"YES","NO")
    # all.cycles.anomaly.table.ewma[cn,5] <- model$s
    # all.cycles.anomaly.table.ewma[cn,6] <- ifelse(length(violations)!=0, max(violations),0)
    all.cycles.anomaly.table.ewma[cn,7] <- time
  }
  
  if (headloss_cycles_info$cycle_length[cn] < 60){
    all.analyzed.cycles.list.ewma[[cn]] <- "Short Cycle"
  }
  
  if(cn%%20==0) cat(cn," ")
  
}

save(all.cycles.anomaly.table.ewma,file = "UWTF.EWMA.results.Rda")