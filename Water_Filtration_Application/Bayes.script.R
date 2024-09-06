load("../Data/headloss_matrix.rda")
load("../Data/headloss_cycles_info.rda")
source("../Functions/hybrid.smoother.freq.R")
source("../Functions/hybrid.smoother.bayes.R")

all.analyzed.cycles.list.bayes <- vector(mode = "list", length = 3305)
all.cycles.anomaly.table.bayes <- data.frame(cycle.number=1:3305,cycle.length=headloss_cycles_info$cycle_length-6,Anomalous=NA,Identified=NA,s=NA,max.gamma=NA,run.time=NA)

for (cn in 1:3305){
  if (headloss_cycles_info$cycle_length[cn] > 60){
    n <- headloss_cycles_info$cycle_length[cn]+1
    x <- 1:(n-6)
    y <- headloss_matrix[cn,4:(n-3)]
    n <- length(x)
    
    init.time <- proc.time()
    model <- hybrid.smoother.bayes(x,y,n.chains = 1,chain.length=1600, burn.in.percent = 1/2)
    time <- time <- (proc.time()-init.time)[3]
    
    all.analyzed.cycles.list.bayes[[cn]] <- model
    names(all.analyzed.cycles.list.bayes)[cn] <- paste0("model.",cn)
    all.cycles.anomaly.table.bayes[cn,4] <- ifelse(length(model$change.points)!=0,"YES","NO")
    all.cycles.anomaly.table.bayes[cn,5] <- model$s
    all.cycles.anomaly.table.bayes[cn,6] <- ifelse(length(model$change.points)!=0, max(abs(model$gammas)),0)
    all.cycles.anomaly.table.bayes[cn,7] <- time
  }
  
  if (headloss_cycles_info$cycle_length[cn] < 60){
    all.analyzed.cycles.list.bayes[[cn]] <- "Short Cycle"
  }
  
  if(cn%%20==0) cat(cn," ")
  
}

save(all.cycles.anomaly.table.bayes,file = "UWTF.Bayes.results.Rda")