#################################################################################
### Code to test the sensitivity of the Schoolfield model to Tref ###############

#####
Length <- seq(5,20,by=1)
Noise <- seq(0.001,0.5,by=0.025)
TestData <- data.frame(merge(Length,Noise))

ETest <- cbind()
for (i in 1:nrow(TestData)){
    ETest <- cbind(ETest,TestTref(TestData$x[i],TestData$y[i]))}

TestTref <- function(Length,Noise){

temp <- seq(0,30,length.out=Length)+273.15
k <- 8.62e-5
B0 <- 0.5
E_D <- 1.2
E <- 0.8
T_ref <- 273.15 + 5
T_h <- 273.15+18
data <- B0 - E/k * (1/temp - 1/T_ref) - log(1 + exp((E_D/k) * (1/T_h - 1/temp)))

set.seed(123)

NewData <- data + rnorm(length(data), Noise)

NewData <- exp(NewData)

MyData <- data.frame(temp,NewData)
# Estimate T.h as being approximately T.peak.
T.h.st  <- GetTpk(tmp=temp, rate=NewData)
E.st    <- GetE(tmp=temp, rate=NewData, T.p=T.h.st)
B.st <- GetB0(tmp=temp, rate=NewData, T.ref=273.15 + 0)

schoolfield_nls0 <- NA
try( 
    schoolfield_nls0 <- nlsLM(
        log(NewData) ~ Schoolfield(B0, E, E_D, T_h, temp = temp,T_ref=273.15+0), 
        start=c(B0 = B.st, E = E.st, E_D = 4*E.st, T_h=T.h.st),
        lower=c(B0=-Inf,   E=0,    E.D=0, Th=0),
        upper=c(B0=Inf,    E=Inf,  E.D=Inf, Th=Inf),
        data=MyData, control=list(minFactor=1 / 2^16, maxiter=1024)),
    silent=FALSE)

T.h.st  <- GetTpk(tmp=temp, rate=NewData)
E.st    <- GetE(tmp=temp, rate=NewData, T.p=T.h.st)
B.st <- GetB0(tmp=temp, rate=NewData, T.ref=273.15 + 5)

schoolfield_nls5 <- NA
try( 
    schoolfield_nls5 <- nlsLM(
        log(NewData) ~ Schoolfield(B0, E, E_D, T_h, temp = temp,T_ref=273.15+5), 
        start=c(B0 = B.st, E = E.st, E_D = 4*E.st, T_h=T.h.st),
        lower=c(B0=-Inf,   E=0,    E.D=0, Th=0),
        upper=c(B0=Inf,    E=Inf,  E.D=Inf, Th=Inf),
        data=MyData, control=list(minFactor=1 / 2^16, maxiter=1024)),
    silent=TRUE)


T.h.st  <- GetTpk(tmp=temp, rate=NewData)
E.st    <- GetE(tmp=temp, rate=NewData, T.p=T.h.st)
B.st <- GetB0(tmp=temp, rate=NewData, T.ref=273.15 + 10)

schoolfield_nls10 <- NA
try( 
    schoolfield_nls10 <- nlsLM(
        log(NewData) ~ Schoolfield(B0, E, E_D, T_h, temp = temp,T_ref=273.15+10), 
        start=c(B0 = B.st, E = E.st, E_D = 4*E.st, T_h=T.h.st),
        lower=c(B0=-Inf,   E=0,    E.D=0, Th=0),
        upper=c(B0=Inf,    E=Inf,  E.D=Inf, Th=Inf),
        data=MyData, control=list(minFactor=1 / 2^16, maxiter=1024)),
    silent=TRUE)

T.h.st  <- GetTpk(tmp=temp, rate=NewData)
E.st    <- GetE(tmp=temp, rate=NewData, T.p=T.h.st)
B.st <- GetB0(tmp=temp, rate=NewData, T.ref=273.15 + 20)

schoolfield_nls20 <- NA
try( 
    schoolfield_nls20 <- nlsLM(
        log(NewData) ~ Schoolfield(B0, E, E_D, T_h, temp = temp,T_ref=273.15+20), 
        start=c(B0 = B.st, E = E.st, E_D = 4*E.st, T_h=T.h.st),
        lower=c(B0=-Inf,   E=0,    E.D=0, Th=0),
        upper=c(B0=Inf,    E=Inf,  E.D=Inf, Th=Inf),
        data=MyData, control=list(minFactor=1 / 2^16, maxiter=1024)),
    silent=TRUE)

T.h.st  <- GetTpk(tmp=temp, rate=NewData)
E.st    <- GetE(tmp=temp, rate=NewData, T.p=T.h.st)
B.st <- GetB0(tmp=temp, rate=NewData, T.ref=273.15 + 30)

schoolfield_nls30 <- NA
try( 
    schoolfield_nls30 <- nlsLM(
        log(NewData) ~ Schoolfield(B0, E, E_D, T_h, temp = temp,T_ref=273.15+30), 
        start=c(B0 = B.st, E = E.st, E_D = 4*E.st, T_h=T.h.st),
        lower=c(B0=-Inf,   E=0,    E.D=0, Th=0),
        upper=c(B0=Inf,    E=Inf,  E.D=Inf, Th=Inf),
        data=MyData, control=list(minFactor=1 / 2^16, maxiter=1024)),
    silent=TRUE)

if (is.na(schoolfield_nls0[1]))  E0 <- NA else E0 <- coef(schoolfield_nls0)["E"]
if (is.na(schoolfield_nls5[1]))  E5 <- NA else E5 <- coef(schoolfield_nls5)["E"]
if (is.na(schoolfield_nls10[1]))  E10 <- NA else E10 <- coef(schoolfield_nls10)["E"]
if (is.na(schoolfield_nls20[1]))  E20 <- NA else E20 <- coef(schoolfield_nls20)["E"]
if (is.na(schoolfield_nls30[1]))  E30 <- NA else E30 <- coef(schoolfield_nls30)["E"]
  
return(c(E0,E5,E10,E20,E30))
}


