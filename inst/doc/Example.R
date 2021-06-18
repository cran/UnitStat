## ----UnitStat, fig.width=8, message=FALSE, warning=FALSE----------------------
library(UnitStat)
# Example 1 : A simple demostration on uniformly distributed data
# y = runif(50,1,49)
# UnitStat(y)
# UnitStat(y,View_results = "T") #To view results at all lags

#Example 2: How the code handles non-stationarity problem
set.seed(1)
a = 0.5
mu= 40
d = 0.5
n = 100
dt = 1
s = matrix(ncol = 1, nrow = n)
s[1] = 40
ep = NULL
pr = NULL
for (i in 1:n){
  pr[i] = runif(1,0.001,0.999)
  ep[i] = qnorm(pr[i],0,1)
  s[i+1] = (s[i]+ (a*(mu-s[i])*dt)+(d*sqrt(dt)*ep[i]))
}

t = 1:99
s1 = s[1:99]
s2 = s1+t

UnitStat(s2, View_results = "T")


