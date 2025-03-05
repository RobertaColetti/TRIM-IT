### Step 1 of TRIM-IT: MOMIP APPLICATION TO GBM DATA 

##...............Function definition ....................
#Run the functions below before executing the main code:
# Data normalization: 
log_norm <-  function(x) log2(x+1)


# Optimization function for solving MOMIP
MOMIP.fun <- function(w) {
  #w: non-negative square matrix with a zero diagonal  
  n <- dim(w)[1]
  # ......................................................
  w.hat <- -3/2*w
  diag(w.hat) <- apply(w,1,sum)  
  # MOMIP solves a QUBO optimization problem.
  # Construction of the coefficients for quadratic terms
  # x_u^2 e  x_u*x_v
  args.list <- list()  
  args.list$A <- matrix(0,nrow=1,ncol=n)
  args.list$Q <- w.hat #Q.aux
  args.list$vtype <- rep('B', n)
  args.list$modelsense="max"
  return(args.list)
} #end of MOMIP.fun

## .................Main Code.....................
# Library section
library(gurobi)  #for MOMIP resolution
# ......................................................

# Load data
load("~/Data/RNAseq-protein_informed-2021WHO.RData")
# Data normalization
DS=apply(GBM_RNA, 2, log_norm)
# Prepare the data for the application of MOMIP (1st)
w.gene=cor(DS)^2  
diag(w.gene) <- 0

# Fix a time limit for optimization
params <- list(OutputFlag = 1, TimeLimit = 10800) 
# Call MOMIP function
CPU.time <- proc.time()
RES.gene1 <- gurobi(MOMIP.fun(w.gene), params) 
CPU.time <- proc.time()-CPU.time
# Result after the 1st step of variable selection
var_sel.1=colnames(DS)[RES.gene1$x==1]

# Prepare the data for the application of MOMIP (2nd)
DS2=DS[,which(colnames(DS) %in% var_sel.1)]
w.gene2=cor(DS2)^2  
diag(w.gene2) <- 0
CPU.time <- proc.time()
RES.gene2 <- gurobi(MOMIP.fun(w.gene2), params) 
CPU.time <- proc.time()-CPU.time
# Result after the 2nd step of variable selection
var_sel.2=colnames(DS2)[RES.gene2$x==1]

#save(var_sel.1,var_sel.2,RES.gene1,RES.gene2,file="MOMIP-GBM-results.RData") #uncomment this line to save your results
