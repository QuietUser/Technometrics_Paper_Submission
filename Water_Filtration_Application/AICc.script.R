load("../Data/headloss_matrix.rda")
load("../Data/headloss_cycles_info.rda")
source("../Functions/hybrid.smoother.freq.R")
source("../Functions/hybrid.smoother.bayes.R")

all.analyzed.cycles.list.freq.AIC <- vector(mode = "list", length = 3305)
all.cycles.anomaly.table.freq.AIC <- data.frame(cycle.number=1:3305,cycle.length=headloss_cycles_info$cycle_length-6,Anomalous=NA,Identified=NA,RMSE=NA,Max.CP=NA,run.time=NA)

for (cn in 1:3305){
  if (headloss_cycles_info$cycle_length[cn] > 60){
    n <- headloss_cycles_info$cycle_length[cn]+1
    x <- 1:(n-6)
    y <- headloss_matrix[cn,4:(n-3)]
    n <- length(x)
    
    lambda.grid <- 2^seq(2,-4,-1)
    omega.grid <- 2^seq(4,16,1)
    
    init.time <- proc.time()
    model <- hybrid.smoother.freq(x,y,lambda.grid,omega.grid,tolerance = 1e-4,method="AICc")
    time <- time <- (proc.time()-init.time)[3]
    all.analyzed.cycles.list.freq.AIC[[cn]] <- model
    names(all.analyzed.cycles.list.freq.AIC)[cn] <- paste0("model.",cn)
    all.cycles.anomaly.table.freq.AIC[cn,4] <- ifelse(length(model$change.points)!=0,"YES","NO")
    all.cycles.anomaly.table.freq.AIC[cn,5] <- sqrt(sum(model$residuals^2)/n)
    all.cycles.anomaly.table.freq.AIC[cn,6] <- ifelse(length(model$change.points)!=0, max(abs(model$gamma.coefs)[model$change.points-1]),NA)
    all.cycles.anomaly.table.freq.AIC[cn,7] <- time
  }
  
  if (headloss_cycles_info$cycle_length[cn] < 60){
    all.analyzed.cycles.list.freq.AIC[[cn]] <- "Short Cycle"
  }
  
  if(cn%%50==0) cat(cn," ")
  
}

save(all.cycles.anomaly.table.freq.AIC,file = "UWTF.AICc.results.Rda")