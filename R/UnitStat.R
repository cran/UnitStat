#' @title Performs Unit Root Test Statistics

#' @description 'A test to understand the stability of the underlying stochastic data.Helps the user understand whether the random variable under consideration is stationary or non-stationary without any manual interpretation of the results.It further ensures to check all the prerequisites and assumptions which are underlying the unit root test statistics and if the underlying data is found to be non-stationary in all the 4 lags the function diagnoses the input data and returns with an optimised.
#' @param y Univariate time series or vector to be tested
#' @param lag Numeric.Default is 0.Select Lags to view results at different lags. Maximum number of lags is 4
#' @param View_results Boolean.Default is False. If True is selected the function returns results for all the 4 lags.

#' @return An object with class UnitStat().
#' Retunrs with a statement explaining the type of input data and its stability.
#' lag - Displays results for the lag number selected
#' View_results - Shows all lag results
#' @examples
#' y = runif(50,1,49)
#' UnitStat(y)
#' UnitStat(y,View_results = "T") #To view results at all lags
#' @author Ankita Sharma
#' @references Dickey, D. A. and Fuller, W. A. (1981), Likelihood Ratio Statistics for Autoregressive Time Series with a Unit Root, Econometrica, 49, 1057--1072. Hamilton (1994), Time Series Analysis, Princeton University Press


#' @importFrom lmtest bgtest
#' @importFrom stats Box.test
#' @importFrom graphics par
#' @importFrom stats anova lm resid
#' @export



UnitStat = function(y, lag = 0 ,View_results = "True"){
  Stationary = function(y){
    if (any(is.na(y)))
      stop("\nNAs in y.\n")

    y = as.vector(y)
    l = length(y)
    if (l < 25){stop("\nSample size too small.\n")}


    tau1 = c(-1.95,-1.95,-1.95,-1.95)
    tau2 = c(-3.00,-2.93,-2.89,-2.88)
    phi1 = c(5.18, 4.86, 4.71,4.63)
    tau3 = c(-3.6,-3.5,-3.45,-3.43)
    phi2 = c(5.68,5.13,4.88,4.75)
    phi3 = c(7.24,6.73,6.49,6.49)
    c_vals = rbind(tau1,tau2,phi1,tau3,phi2,phi3)
    colnames(c_vals) = c(25,50,100,250)


    n50 = c_vals[,2]
    n100 = c_vals[,3]
    n250 = c_vals[,4]


    if (25 <= l & l < 50) {n = n50}
    if (50 <= l & l < 100){n = n100}
    if (l >= 100 ){n = n250}


    lags = 4

    for (i in 1:lags){
      dy = diff(y)
      tt = 1:length(y)

      if (i == 1) {
        dy1 = dy [2:length(dy)]
        yt_1 = y [2: (length(y)-1)]
        lag1 = dy[1:(length(dy)-1)]
        tt   = tt[3: length(tt)]

        reg1 = lm(dy1~yt_1 -1 + lag1)
        reg2 = lm(dy1~ yt_1 +1 + lag1)
        reg3 = lm(dy1~ yt_1 +1 +tt +lag1)

        t_r1 = round(summary(reg1)$coefficients[1,3],digits = 2)
        t_r2 = round(summary(reg2)$coefficients[2,3],digits = 2)
        t_r3 = round(summary(reg3)$coefficients[2,3],digits = 2)

        if (reg1$coefficients[1] < 0) {wn1 = Box.test(resid(reg1), type = "Ljung-Box", lag = 1)} else {wn1 = "0"}
        if (wn1[1:1] != "0" & wn1[3:3] > 0.05) {wn1 = "1"} else {wn1 = "0"}
        if (reg2$coefficients[2] < 0) {wn2 = Box.test(resid(reg2), type = "Ljung-Box", lag = 1)} else {wn2 = "0"}
        if (wn2[1:1] != "0" & wn2[3:3] > 0.05) {wn2 = "1"} else {wn2 = "0"}
        if (reg3$coefficients[2] < 0) {wn3 = Box.test(resid(reg3), type = "Ljung-Box", lag = 1)} else {wn3 = "0"}
        if (wn3[1:1] != "0" & wn3[3:3] > 0.05) {wn3 = "1"} else {wn3 = "0"}

        if (reg1$coefficients[1] < 0) {ar1 = bgtest(reg1,type = "Chisq", order = 1)} else {ar1 = "0"}
        if (ar1[1:1] != "0" & ar1[4:4] > 0.05) {ar1 = "1"} else {ar1 = "0"}
        if (reg2$coefficients[2] < 0) {ar2 = bgtest(reg1,type = "Chisq", order = 1)} else {ar2 = "0"}
        if (ar2[1:1] != "0" & ar2[4:4] > 0.05) {ar2 = "1"} else {ar2 = "0"}
        if (reg3$coefficients[2] < 0) {ar3 = bgtest(reg1,type = "Chisq", order = 1)} else {ar3 = "0"}
        if (ar3[1:1] != "0" & ar3[4:4] > 0.05) {ar3 = "1"} else {ar3 = "0"}


        phi11 = lm(dy1~ -1 + lag1)
        phi1 = round(anova(phi11, reg2)$F[2], digits = 2)


        phi21 =  lm(dy1~ -1 + lag1)
        phi31 = lm(dy1~ 1 + lag1)
        phi2 = round(anova(phi21,reg3)$F[2], digits = 2)
        phi3 = round(anova(phi31,reg3)$F[2], digits = 2)


        Rho =rbind(round((reg1$coefficients[1]), digits = 4),
                   round((reg2$coefficients[2]), digits = 4),
                   round((reg3$coefficients[2]), digits = 4),
                   round((reg3$coefficients[2]), digits = 4))


        rr = rbind(ar1,ar2,ar3,ar3)
        rw = rbind(wn1,wn2,wn3,wn3)

        Test_name = rbind((paste("wn1","-","ar1")),(paste("wn2","-","ar2")),(paste("wn3","-","ar3")),(paste("wn3","-","ar3")))
        Equation =rbind("None","Drift equation","Trend equation","Trend equation")
        t = rbind(t_r1,t_r2,t_r3,t_r3)
        f = rbind("na", phi1, phi2,phi3)
        r1 = data.frame(Equation,t,f)
        r1[is.na(r1)] = 0

        r2 = ifelse((rr == "1")|(rw == "1"), "Assumption passed", "Assumption failed")

        t1 = ifelse(r2 == "Assumption passed" & t_r1 < n[1], "stationary", "non-stationary")
        t2 = ifelse(r2 == "Assumption passed" & t_r2 < n[2], "stationary", "non-stationary")
        t3 = ifelse(r2 == "Assumption passed" & t_r3 < n[4], "stationary", "non-stationary")
        r3 = rbind(t1[1],t2[1],t3[1],t3[1])
        r4 = ifelse(r2 == "Assumption passed" & phi1 > n[3], "drift present", "insignificant drift")
        r5 = ifelse(r2 == "Assumption passed" & phi2 > n[5], "drift & trend present", "insignificant drift & trend")
        r6 = ifelse(r2 == "Assumption passed" & phi3 > n[6], "trend present", "insignificant trend")

        f1 = rbind("No deterministic drift & trend", r4[2],r5[2],r6[3])
        f_0 = paste("Lag",i,sep=" ")
        f0 = rbind(f_0,f_0,f_0,f_0)

        result1 = data.frame(f0,r1,r2,r3,f1)
        rownames(result1) = c("1","2","3","4")
        colnames(result1) = c("Lags","Equation type","T-test",
                              "phi1/phi2/phi3","Assumption Outcome" ,"Tau Result","Phi Result")


      }


      if (i == 2) {
        dy1 = dy[3:length(dy)]
        yt_1 = y[3:(length(y)-1)]
        lag1 = dy[2:(length(dy)-1)]
        lag2 = dy[1:(length(dy)-2)]
        tt   = tt[4 :length(tt)]


        reg1 = lm(dy1~yt_1 -1 + lag1 + lag2)
        reg2 = lm(dy1~ yt_1 +1 + lag1 + lag2)
        reg3 = lm(dy1~ yt_1 +1 +tt +lag1+ lag2)

        t_r1 = round(summary(reg1)$coefficients[1,3],digits = 2)
        t_r2 = round(summary(reg2)$coefficients[2,3],digits = 2)
        t_r3 = round(summary(reg3)$coefficients[2,3],digits = 2)

        if (reg1$coefficients[1] < 0) {wn1 = Box.test(resid(reg1), type = "Ljung-Box", lag = 1)} else {wn1 = "0"}
        if (wn1[1:1] != "0" & wn1[3:3] > 0.05) {wn1 = "1"} else {wn1 = "0"}
        if (reg2$coefficients[2] < 0) {wn2 = Box.test(resid(reg2), type = "Ljung-Box", lag = 1)} else {wn2 = "0"}
        if (wn2[1:1] != "0" & wn2[3:3] > 0.05) {wn2 = "1"} else {wn2 = "0"}
        if (reg3$coefficients[2] < 0) {wn3 = Box.test(resid(reg3), type = "Ljung-Box", lag = 1)} else {wn3 = "0"}
        if (wn3[1:1] != "0" & wn3[3:3] > 0.05) {wn3 = "1"} else {wn3 = "0"}

        if (reg1$coefficients[1] < 0) {ar1 = bgtest(reg1,type = "Chisq", order = 1)} else {ar1 = "0"}
        if (ar1[1:1] != "0" & ar1[4:4] > 0.05) {ar1 = "1"} else {ar1 = "0"}
        if (reg2$coefficients[2] < 0) {ar2 = bgtest(reg1,type = "Chisq", order = 1)} else {ar2 = "0"}
        if (ar2[1:1] != "0" & ar2[4:4] > 0.05) {ar2 = "1"} else {ar2 = "0"}
        if (reg3$coefficients[2] < 0) {ar3 = bgtest(reg1,type = "Chisq", order = 1)} else {ar3 = "0"}
        if (ar3[1:1] != "0" & ar3[4:4] > 0.05) {ar3 = "1"} else {ar3 = "0"}



        phi11 = lm(dy1~ -1 + lag1)
        phi1 = round(anova(phi11, reg2)$F[2], digits = 2)


        phi21 =  lm(dy1~ -1 + lag1)
        phi31 = lm(dy1~ 1 + lag1)
        phi2 = round(anova(phi21,reg3)$F[2], digits = 2)
        phi3 = round(anova(phi31,reg3)$F[2], digits = 2)


        Rho =rbind(round((reg1$coefficients[1]), digits = 4),
                   round((reg2$coefficients[2]), digits = 4),
                   round((reg3$coefficients[2]), digits = 4),
                   round((reg3$coefficients[2]), digits = 4))


        rr = rbind(ar1,ar2,ar3,ar3)
        rw = rbind(wn1,wn2,wn3,wn3)

        Test_name = rbind((paste("wn1","-","ar1")),(paste("wn2","-","ar2")),(paste("wn3","-","ar3")),(paste("wn3","-","ar3")))
        Equation =rbind("None","Drift equation","Trend equation","Trend equation")
        t = rbind(t_r1,t_r2,t_r3,t_r3)
        f = rbind("na", phi1, phi2,phi3)
        r1 = data.frame(Equation,t,f)
        r1[is.na(r1)] = 0

        r2 = ifelse((rr == "1")|(rw == "1"), "Assumption passed", "Assumption failed")

        t1 = ifelse(r2 == "Assumption passed" & t_r1 < n[1], "stationary", "non-stationary")
        t2 = ifelse(r2 == "Assumption passed" & t_r2 < n[2], "stationary", "non-stationary")
        t3 = ifelse(r2 == "Assumption passed" & t_r3 < n[4], "stationary", "non-stationary")
        r3 = rbind(t1[1],t2[1],t3[1],t3[1])
        r4 = ifelse(r2 == "Assumption passed" & phi1 > n[3], "drift present", "insignificant drift")
        r5 = ifelse(r2 == "Assumption passed" & phi2 > n[5], "drift & trend present", "insignificant drift & trend")
        r6 = ifelse(r2 == "Assumption passed" & phi3 > n[6], "trend present", "insignificant trend")

        f1 = rbind("No deterministic drift & trend", r4[2],r5[2],r6[3])
        f_0 = paste("Lag",i,sep=" ")
        f0 = rbind(f_0,f_0,f_0,f_0)

        result2 = data.frame(f0,r1,r2,r3,f1)
        rownames(result2) = c("1","2","3","4")
        colnames(result2) = c("Lags","Equation type","T-test",
                              "phi1/phi2/phi3","Assumption Outcome" ,"Tau Result","Phi Result")

      }

      if (i == 3) {
        dy1 = dy[4:length(dy)]
        yt_1 = y[4: (length(y)-1)]
        lag1 = dy[3:(length(dy)-1)]
        lag2 = dy[2:(length(dy)-2)]
        lag3 = dy[1: (length(dy)-3)]
        tt  = tt[5 : length(tt)]

        reg1 = lm(dy1~yt_1 -1 + lag1 + lag2 + lag3)
        reg2 = lm(dy1~ yt_1 +1 + lag1 + lag2 + lag3)
        reg3 = lm(dy1~ yt_1 +1 +tt +lag1+ lag2 + lag3)

        t_r1 = round(summary(reg1)$coefficients[1,3],digits = 2)
        t_r2 = round(summary(reg2)$coefficients[2,3],digits = 2)
        t_r3 = round(summary(reg3)$coefficients[2,3],digits = 2)

        if (reg1$coefficients[1] < 0) {wn1 = Box.test(resid(reg1), type = "Ljung-Box", lag = 1)} else {wn1 = "0"}
        if (wn1[1:1] != "0" & wn1[3:3] > 0.05) {wn1 = "1"} else {wn1 = "0"}
        if (reg2$coefficients[2] < 0) {wn2 = Box.test(resid(reg2), type = "Ljung-Box", lag = 1)} else {wn2 = "0"}
        if (wn2[1:1] != "0" & wn2[3:3] > 0.05) {wn2 = "1"} else {wn2 = "0"}
        if (reg3$coefficients[2] < 0) {wn3 = Box.test(resid(reg3), type = "Ljung-Box", lag = 1)} else {wn3 = "0"}
        if (wn3[1:1] != "0" & wn3[3:3] > 0.05) {wn3 = "1"} else {wn3 = "0"}

        if (reg1$coefficients[1] < 0) {ar1 = bgtest(reg1,type = "Chisq", order = 1)} else {ar1 = "0"}
        if (ar1[1:1] != "0" & ar1[4:4] > 0.05) {ar1 = "1"} else {ar1 = "0"}
        if (reg2$coefficients[2] < 0) {ar2 = bgtest(reg1,type = "Chisq", order = 1)} else {ar2 = "0"}
        if (ar2[1:1] != "0" & ar2[4:4] > 0.05) {ar2 = "1"} else {ar2 = "0"}
        if (reg3$coefficients[2] < 0) {ar3 = bgtest(reg1,type = "Chisq", order = 1)} else {ar3 = "0"}
        if (ar3[1:1] != "0" & ar3[4:4] > 0.05) {ar3 = "1"} else {ar3 = "0"}



        phi11 = lm(dy1~ -1 + lag1)
        phi1 = round(anova(phi11, reg2)$F[2], digits = 2)


        phi21 =  lm(dy1~ -1 + lag1)
        phi31 = lm(dy1~ 1 + lag1)
        phi2 = round(anova(phi21,reg3)$F[2], digits = 2)
        phi3 = round(anova(phi31,reg3)$F[2], digits = 2)


        Rho =rbind(round((reg1$coefficients[1]), digits = 4),
                   round((reg2$coefficients[2]), digits = 4),
                   round((reg3$coefficients[2]), digits = 4),
                   round((reg3$coefficients[2]), digits = 4))


        rr = rbind(ar1,ar2,ar3,ar3)
        rw = rbind(wn1,wn2,wn3,wn3)

        Test_name = rbind((paste("wn1","-","ar1")),(paste("wn2","-","ar2")),(paste("wn3","-","ar3")),(paste("wn3","-","ar3")))
        Equation =rbind("None","Drift equation","Trend equation","Trend equation")
        t = rbind(t_r1,t_r2,t_r3,t_r3)
        f = rbind("na", phi1, phi2,phi3)
        r1 = data.frame(Equation,t,f)
        r1[is.na(r1)] = 0

        r2 = ifelse((rr == "1")|(rw == "1"), "Assumption passed", "Assumption failed")

        t1 = ifelse(r2 == "Assumption passed" & t_r1 < n[1], "stationary", "non-stationary")
        t2 = ifelse(r2 == "Assumption passed" & t_r2 < n[2], "stationary", "non-stationary")
        t3 = ifelse(r2 == "Assumption passed" & t_r3 < n[4], "stationary", "non-stationary")
        r3 = rbind(t1[1],t2[1],t3[1],t3[1])
        r4 = ifelse(r2 == "Assumption passed" & phi1 > n[3], "drift present", "insignificant drift")
        r5 = ifelse(r2 == "Assumption passed" & phi2 > n[5], "drift & trend present", "insignificant drift & trend")
        r6 = ifelse(r2 == "Assumption passed" & phi3 > n[6], "trend present", "insignificant trend")

        f1 = rbind("No deterministic drift & trend", r4[2],r5[2],r6[3])
        f_0 = paste("Lag",i,sep=" ")
        f0 = rbind(f_0,f_0,f_0,f_0)

        result3 = data.frame(f0,r1,r2,r3,f1)
        rownames(result3) = c("1","2","3","4")
        colnames(result3) = c("Lags","Equation type","T-test",
                              "phi1/phi2/phi3","Assumption Outcome" ,"Tau Result","Phi Result")

      }


      if (i == 4) {
        dy1 =  dy[5:length(dy)]
        yt_1 =  y[5: (length(y)-1)]
        lag1 = dy[4:(length(dy)-1)]
        lag2 = dy[3:(length(dy)-2)]
        lag3 = dy[2: (length(dy)-3)]
        lag4 = dy[1:(length(dy)-4)]
        tt  =  tt[6 : length(tt)]

        reg1 = lm(dy1~yt_1 -1 + lag1 + lag2 + lag3 + lag4)
        reg2 = lm(dy1~ yt_1 +1 + lag1 + lag2 + lag3 + lag4)
        reg3 = lm(dy1~ yt_1 +1 +tt +lag1+ lag2 + lag3 + lag4)

        t_r1 = round(summary(reg1)$coefficients[1,3],digits = 2)
        t_r2 = round(summary(reg2)$coefficients[2,3],digits = 2)
        t_r3 = round(summary(reg3)$coefficients[2,3],digits = 2)

        if (reg1$coefficients[1] < 0) {wn1 = Box.test(resid(reg1), type = "Ljung-Box", lag = 1)} else {wn1 = "0"}
        if (wn1[1:1] != "0" & wn1[3:3] > 0.05) {wn1 = "1"} else {wn1 = "0"}
        if (reg2$coefficients[2] < 0) {wn2 = Box.test(resid(reg2), type = "Ljung-Box", lag = 1)} else {wn2 = "0"}
        if (wn2[1:1] != "0" & wn2[3:3] > 0.05) {wn2 = "1"} else {wn2 = "0"}
        if (reg3$coefficients[2] < 0) {wn3 = Box.test(resid(reg3), type = "Ljung-Box", lag = 1)} else {wn3 = "0"}
        if (wn3[1:1] != "0" & wn3[3:3] > 0.05) {wn3 = "1"} else {wn3 = "0"}

        if (reg1$coefficients[1] < 0) {ar1 = bgtest(reg1,type = "Chisq", order = 1)} else {ar1 = "0"}
        if (ar1[1:1] != "0" & ar1[4:4] > 0.05) {ar1 = "1"} else {ar1 = "0"}
        if (reg2$coefficients[2] < 0) {ar2 = bgtest(reg1,type = "Chisq", order = 1)} else {ar2 = "0"}
        if (ar2[1:1] != "0" & ar2[4:4] > 0.05) {ar2 = "1"} else {ar2 = "0"}
        if (reg3$coefficients[2] < 0) {ar3 = bgtest(reg1,type = "Chisq", order = 1)} else {ar3 = "0"}
        if (ar3[1:1] != "0" & ar3[4:4] > 0.05) {ar3 = "1"} else {ar3 = "0"}



        phi11 = lm(dy1~ -1 + lag1)
        phi1 = round(anova(phi11, reg2)$F[2], digits = 2)


        phi21 =  lm(dy1~ -1 + lag1)
        phi31 = lm(dy1~ 1 + lag1)
        phi2 = round(anova(phi21,reg3)$F[2], digits = 2)
        phi3 = round(anova(phi31,reg3)$F[2], digits = 2)


        Rho =rbind(round((reg1$coefficients[1]), digits = 4),
                   round((reg2$coefficients[2]), digits = 4),
                   round((reg3$coefficients[2]), digits = 4),
                   round((reg3$coefficients[2]), digits = 4))


        rr = rbind(ar1,ar2,ar3,ar3)
        rw = rbind(wn1,wn2,wn3,wn3)

        Test_name = rbind((paste("wn1","-","ar1")),(paste("wn2","-","ar2")),(paste("wn3","-","ar3")),(paste("wn3","-","ar3")))
        Equation =rbind("None","Drift equation","Trend equation","Trend equation")
        t = rbind(t_r1,t_r2,t_r3,t_r3)
        f = rbind("na", phi1, phi2,phi3)
        r1 = data.frame(Equation,t,f)
        r1[is.na(r1)] = 0

        r2 = ifelse((rr == "1")|(rw == "1"), "Assumption passed", "Assumption failed")

        t1 = ifelse(r2 == "Assumption passed" & t_r1 < n[1], "stationary", "non-stationary")
        t2 = ifelse(r2 == "Assumption passed" & t_r2 < n[2], "stationary", "non-stationary")
        t3 = ifelse(r2 == "Assumption passed" & t_r3 < n[4], "stationary", "non-stationary")
        r3 = rbind(t1[1],t2[1],t3[1],t3[1])
        r4 = ifelse(r2 == "Assumption passed" & phi1 > n[3], "drift present", "insignificant drift")
        r5 = ifelse(r2 == "Assumption passed" & phi2 > n[5], "drift & trend present", "insignificant drift & trend")
        r6 = ifelse(r2 == "Assumption passed" & phi3 > n[6], "trend present", "insignificant trend")

        f1 = rbind("No deterministic drift & trend", r4[2],r5[2],r6[3])
        f_0 = paste("Lag",i,sep=" ")
        f0 = rbind(f_0,f_0,f_0,f_0)

        result4 = data.frame(f0,r1,r2,r3,f1)
        rownames(result4) = c("1","2","3","4")
        colnames(result4) = c("Lags","Equation type","T-test",
                              "phi1/phi2/phi3","Assumption Outcome" ,"Tau Result","Phi Result")

      }

    }

    z = rbind(result1, result2,result3,result4)

    r_ns = rbind(paste(result1$`Tau Result`[1]),paste(result2$`Tau Result`[1]),
                 paste(result3$`Tau Result`[1]),paste(result4$`Tau Result`[1]))


    r_ds = rbind(paste(result1$`Tau Result`[2]),paste(result2$`Tau Result`[2]),
                 paste(result3$`Tau Result`[2]),paste(result4$`Tau Result`[2]))

    r_d = rbind(paste(result1$`Phi Result`[2]),paste(result2$`Phi Result`[2]),
                paste(result3$`Phi Result`[2]),paste(result4$`Phi Result`[2]))

    r_dts = rbind(paste(result1$`Tau Result`[3]),paste(result2$`Tau Result`[3]),
                  paste(result3$`Tau Result`[3]),paste(result4$`Tau Result`[3]))

    r_dt = rbind(paste(result1$`Phi Result`[3]),paste(result2$`Phi Result`[3]),
                 paste(result3$`Phi Result`[3]),paste(result4$`Phi Result`[3]))


    r_ts = rbind(paste(result1$`Tau Result`[4]),paste(result2$`Tau Result`[4]),
                 paste(result3$`Tau Result`[4]),paste(result4$`Tau Result`[4]))

    r_t = rbind(paste(result1$`Phi Result`[4]),paste(result2$`Phi Result`[4]),
                paste(result3$`Phi Result`[4]),paste(result4$`Phi Result`[4]))



    s_0 = "stationary with No deterministic drift & trend"
    s_1 = "stationary with drift present"
    s_2 = "stationary with drift & trend present"
    s_3 = "stationary with trend present"

    p = z$`Phi Result`
    l = z$Lags
    t = z$`Tau Result`
    xx = cbind(paste(t,"with",p))
    res1 = 0
    res2 = 0
    res = 0
    c2 = 0
    for (i in 1: length(xx)){
      if ((xx[i] == s_0)||(xx[i] == s_1)||(xx[i] == s_2)||(xx[i] == s_3))
      {
        c2 = c2+1
      }
      if (c2!= 0)
      {
        res1 = xx[i]
        res2 = l[i]
        res = paste("Time-series is",res1,"at",res2)
        break
      }
    }

    cd  = 0
    cdt = 0
    ct  = 0

    for (i in 1:4){

      if(res!= 0){
        out = res
      }

      if(res == 0){

        if(r_ns[i] == "stationary")
        {
          out =(paste("Time-series is stationary at Lag",i))
          break
        }
        if(r_ns[i] == "non-stationary")
        {

          if(r_ds[i] == "stationary")
          {
            if(r_d[i] == "drift present")
            {
              out = (paste("Time-series is stationary with drift present at Lag",i))

              break
            }
          }

          if(r_ds[i] == "stationary")
          {
            if(r_d[i] == "insignificant drift")
            {
              cd = cd + 1
              if(cd != 0){out = paste("Time-series is stationary with insignificant drift at Lag", i)}
              break
            }
          }

          if(r_ds[i] == "non-stationary")
          {
            if(r_d[i] == "drift present")
            {
              out =(paste("Time-series is non-stationary with drift present at Lag",i))

              break
            }
          }
          if(r_ds[i] == "non-stationary")
          {
            if(r_d[i]== "insignificant drift")
            {

              if(r_dts[i] == "stationary")
              {
                if(r_dt[i] == "drift & trend present")
                {
                  out =(paste("Time-series is stationary with drift and trend present at Lag",i))
                  break
                }
              }
              if(r_dts[i] == "stationary")
              {
                if(r_dt[i] == "insignificant drift & trend")
                {

                  cdt = cdt + 1

                  if(cdt != 0){out = paste("Time-series is stationary with insignificant drift and trend at Lag", i)}
                  break
                }
              }

              if(r_dts[i] == "non-stationary")
              {
                if(r_dt[i] == "drift & trend present")
                {
                  out =(paste("Time-series is non-stationary with drift & trend present at Lag",i))

                  break
                }
              }
              if(r_dts[i] == "non-stationary"){
                {
                  if(r_dt[i] == "insignificant drift & trend")
                  {

                    if(r_ts[i] == "stationary")
                    {
                      if(r_t[i] == "trend present")
                      {
                        out =(paste("Time-series is stationary with trend present at Lag",i))
                        break
                      }
                    }

                    if(r_ts[i] == "stationary")
                    {
                      if(r_t[i] == "insignificant trend")
                      {

                        ct = ct + 1
                        if (ct != 0){out = paste("Time-series is stationary with insignificant trend at Lag", i)}
                        break
                      }
                    }

                    if(r_ts[i] == "non-stationary")
                    {
                      if(r_t[i] == "trend present")
                      {
                        out =(paste("Time-series is non-stationary with trend present at Lag",i))
                        break
                      }
                    }
                    if(r_ts[i] == "non-stationary")

                    {
                      if(r_t[i] == "insignificant trend")
                      {
                        out =(paste("Time-series is non-stationary at Lag",i))


                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    if (out == "Time-series is non-stationary at Lag 4"){
      out = "Time-series is non-stationary at all 4 Lags"
    }


    rownames(z) = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
    z1 = list(All_results = z, Outcome = out, critical_values = n)
    return(z1)

  }

  par(mfrow = c(1,2))
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  p1 = plot(y,type ="l", main = "Random variable_Input", xlab = "No. of Data points",ylab = "Input variable")

  result1 = Stationary(y)
  result_f = result1[[1]]
  result = result1[[2]]

  ns_d1 = "Time-series is non-stationary with drift present at Lag 1"
  ns_d2 = "Time-series is non-stationary with drift present at Lag 2"
  ns_d3 = "Time-series is non-stationary with drift present at Lag 3"
  ns_d4 = "Time-series is non-stationary with drift present at Lag 4"

  ns_dt1 = "Time-series is non-stationary with drift and trend present at Lag 1"
  ns_dt2 = "Time-series is non-stationary with drift and trend present at Lag 2"
  ns_dt3 = "Time-series is non-stationary with drift and trend present at Lag 3"
  ns_dt4 = "Time-series is non-stationary with drift and trend present at Lag 4"

  ns_t1 = "Time-series is non-stationary with trend present at Lag 1"
  ns_t2 = "Time-series is non-stationary with trend present at Lag 2"
  ns_t3 = "Time-series is non-stationary with trend present at Lag 3"
  ns_t4 = "Time-series is non-stationary with trend present at Lag 4"


  s1 = "Time-series is stationary at lag 1"
  s2 = "Time-series is stationary at lag 2"
  s3 = "Time-series is stationary at lag 3"
  s4 = "Time-series is stationary at lag 4"



  NS = "Time-series is non-stationary at all 4 Lags"



  trend = 1:length(y)
  ny1 = 0
  ny2 = 0
  ny3 = 0
  for (j in 1:4){
    if((result == s1)||(result == s2)||(result == s3)||(result == s3)){
      break
    }
    if ((result == ns_d1)||(result == ns_d2)||(result == ns_d3)||(result == ns_d4)){
      ny1 = diff(y)
    }
    if((result == ns_dt1)||(result == ns_dt2)||(result == ns_dt3)||(result == ns_dt4)){
      ny2 = resid(lm(y~trend+1))
    }
    if((result == ns_t1)||(result == ns_t2)||(result == ns_t3)||(result == ns_t4)){
      ny3 = resid(lm(y~trend+1))
    }
    if(result == NS){
      break
    }
  }



  z = 0
  if (ny1[1:1] != 0){
    p2 = plot(ny1,type ="l", main = "Random variable_Transformed")
    z = ny1
    q = print(paste("The input", result,"however, after taking the First Difference"))
  }
  if (ny2[1:1] != 0){
    p2 = plot(ny2,type ="l", main = "Random variable_Transformed")
    z = ny2
    q = print(paste("The input", result,"however,after De-trending/De-drifting"))
  }
  if (ny3[1:1] != 0){
    p2 = plot(ny3,type ="l", main = "Random variable_Transformed")
    z = ny3
    q = print(paste("The input", result,"however,after De-trending"))
  }

  if (z[1:1] != 0){Stationary(y)}


  Out1 = Stationary(y)
  Out2 = Out1[[1]]
  Out3 = Out1[[2]]
  res1 = Out2[1:4,]
  res2 = Out2[5:8,]
  res3 = Out2[9:12,]
  res4 = Out2[13:16,]
  cval = Out1[[3]]


  if (View_results == "True"){
    print(list(All_Results = Out2, Critical_values = cval))
  } else {
    print(Out3)
  }

  lg = ""
  if (lag == 1){
    lg = list(Results = res1, Critical_Values = cval)
    print(lg)
  }

  if (lag == 2){
    lg = list(Results = res2, Critical_Values = cval)
    print(lg)
  }
  if (lag == 3){
    lg = list(Results = res3, Critical_Values = cval)
    print(lg)
  }
  if (lag == 4){
    lg = list(Results = res4, Critical_Values = cval)
    print(lg)
  }


}



























