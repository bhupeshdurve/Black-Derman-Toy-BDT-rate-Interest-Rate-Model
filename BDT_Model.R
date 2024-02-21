library(tidyverse)
#————————————————
# construct short rate tree (r0)
#————————————————
generate_BDT_tree = function(v.rd, v.rv, nmat) {
  
  df.rt  = matrix(0, nmat, nmat) 
  df.rt[1,1] = v.rd[1] # r0
  for(t in 2:nmat) {
    df.rt[t,t] = v.rd[t]
    for(s in (t-1):1) { # cross-sectional iteration
      df.rt[s,t] = df.rt[s+1,t]*exp(2*v.rv[t])
    }
  }
  return(df.rt)
}

#————————————————
# discount cash flow at maturity using BDT tree
#————————————————
pricing_BDT_tree <- function(df.BDT_tree, v.cf, nmat) {
  
  df.rt <- df.BDT_tree
  t     <- nmat
  df.pt <- matrix(0, t+1, t+1)
  df.pt[,t+1] <- v.cf # cash flow at maturity
  
  for(b in t:1) { # backward induction
    for(s in 1:b) { # cross-sectional iteration
      df.pt[s,b] <- 0.5*
        (df.pt[s,b+1]+df.pt[s+1,b+1])/(1+df.rt[s,b])
    }
  }
  return(df.pt)
}

#————————————————
# objective function
#————————————————
# v.unknown : rd, rdd, …, rv2, rv3, …
# df.mkt    : data.frame for mat, y, yv
#————————————————
objfunc_BDT_tree <- function(v.unknown, df.mkt) {
  
  v.unknown <- v.unknown^2 # non-negativity 
  nmat <- length(df.mkt$mat)
  
  # short rate and its volatility 
  v.rd <- c(df.mkt$y[1], v.unknown[1:(nmat-1)])
  v.rv <- c(df.mkt$yv[1],v.unknown[nmat:(2*(nmat-1))])
  
  # construct rate tree(rt)
  df.rt <- generate_BDT_tree(v.rd, v.rv, nmat)
  
  #————————————————
  # make a price tree for each maturity
  # calculate discount factor difference
  # calculate yield volatility difference
  #————————————————
  v.df_diff <- v.yv_diff <- rep(0,nmat-1)
  
  for(t in 2:nmat) {
    df.pt <- pricing_BDT_tree(df.rt, 1, t) 
    
    # difference between model and 
    # market discount factor at time 0
    v.df_diff[t-1] <- df.pt[1,1] - 1/(1+df.mkt$y[t])^t
    
    # difference between yield volatilities 
    # from up & down discount factors (model) and
    # market discount factors
    yu <- (1/df.pt[1,2])^(1/(t-1))-1
    yd <- (1/df.pt[2,2])^(1/(t-1))-1
    yv <- log(yu/yd)/2
    v.yv_diff[t-1] <- yv - df.mkt$yv[t]
  }
  return(sqrt(sum(v.df_diff^2) + sum(v.yv_diff^2)))
}

#————————————————--------------------------
#                   Main
#————————————————-----------------------------

# Input data : maturity (1-30), yield (continuously compounded, give D(T)), yield volatility (given sigmas) 
df.mkt <- data.frame(
  mat = c(1:30),            
  y   = c(0.0566, 0.0578, 0.0591, 0.0603, 0.0614, 0.0622, 0.0629, 0.0635, 0.0639, 0.0643, 0.0647, 0.0651, 0.0654, 0.0658, 0.0661, 0.0664, 0.0667, 0.067, 0.0672, 0.0674, 0.0676, 0.0678, 0.068, 0.0681, 0.0682, 0.0683, 0.0684, 0.0685, 0.0686, 0.0687),
  yv    = c(0,0.1,0.12,0.135,0.15,0.16,0.162,0.164,0.162,0.16,0.157,0.154,0.151,0.148,0.145,0.142,0.139,0.136,0.133,0.13,0.127,0.125,0.123,0.12,0.118,0.116,0.114,0.112,0.11,0.108)
)

# initial guess for short rate and volatility
# short rate : rd, rdd, rddd, rdddd
# short rate volatilities : rv2, rv3, rv4, rv5 
v.initial_guess <- c( 0.05, 0.05,  0.05,  0.05,0.05,0.05 ,0.05 ,0.05 ,0.05 ,0.05 ,0.05 ,0.05 ,0.05 ,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05, # rd
                      0.0211, 0.0591, 0.0603, 0.0614, 0.0622, 0.0629, 0.0635, 0.0639, 0.0643, 0.0647, 0.0651, 0.0654, 0.0658, 0.0661, 0.0664, 0.0667, 0.067, 0.0672, 0.0674, 0.0676, 0.0678, 0.068, 0.0681, 0.0682, 0.0683, 0.0684, 0.0685, 0.0686, 0.0687) # rv
v.initial_guess <- sqrt(v.initial_guess)

# calibration of BDT tree
m<-optim(v.initial_guess, objfunc_BDT_tree,
         control = list(maxit=5000, trace=2, reltol = 1e-16),
         method=c("BFGS"),df.mkt = df.mkt)

# transform and split parameters
v.param_trans <- m$par^2
nmat <- 30

# final calibrated parameters
# short rate : r0, rd, rdd, rddd, rdddd
(v.rd <- c(df.mkt$y[1],  v.param_trans[1:(nmat-1)]))
# short rate volatilities : for 29 periods
(v.rv <- c(df.mkt$yv[1], v.param_trans[nmat:(2*(nmat-1))]))



# final BDT short rate tree
(df.rt <- generate_BDT_tree(v.rd, v.rv, nmat))

write.csv(df.rt, file = "BDT_short_rate_tree.csv")
BDT_rates <- colSums(df.rt)/c(1:30)

write.csv(BDT_rates, file = "short_rates.csv")
# final BDT price tree (30-periods)
write.csv(Dt <- pricing_BDT_tree(df.rt, 1, nmat), file = "BDT_discount_Tree.csv")

#  list of 6 month forward rates calculated in excel

forward_rates <- c(0.0566, 0.0591, 0.0616, 0.0639, 0.0656, 0.0666, 0.067, 0.0672, 0.0675, 0.0679, 0.0686, 0.0692, 0.0698, 0.0703, 0.0707, 0.0711, 0.0713, 0.0715, 0.0716, 0.0716, 0.0716, 0.0714, 0.0713, 0.0712, 0.071, 0.0709, 0.0709, 0.0709, 0.0709, 0.0711)
rates <- sapply(seq(0.05, 0.19, by = 0.0047),head,30)
length(rates)
library(ggplot2)
time <- seq(0.5, 15, by=0.5)
data <- data.frame(time, BDT_rates, forward_rates,rates)
ggplot(data, aes(x=time)) + 
  geom_line(aes(y = BDT_rates), color = "tomato") + 
  geom_line(aes(y = forward_rates), color="steelblue", linetype="twodash") 
