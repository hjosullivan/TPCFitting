##### Self-start function test


model_function<-function(Tpred,a,TO,Tl){a*Tpred*(Tpred-TO)*sqrt(Tl-Tpred)}
#function   for the model

model_initial<-function(mCall,LHS,data){
xy<-sortedXyData(mCall[['Tpred']],LHS,data)
fit<-lm(xy[,'y']~xy[,'x'])
coefs<-coef(fit)
a<-coefs[1]    
TO<-coefs[2]
Tl<-coefs[3]
value<-c(a,TO,Tl)

names(value)<-mCall[c('a','TO','Tl')]
value
}

Schoolfield <- function(x,B0,E, E_D, T_h) {
    T_ref <- 273.15 + 10
    
    # Boltzmann's constant. Units imply that E and E_D are in eV.
    k <- 8.62e-5

    log(B0 * exp(-E * ((1/(k * x)) - (1/(k * T_ref)))) / 
        (1 + exp(E_D/k * ((1/T_h) - (1/x)))))

}

School.Init <- function(mcall,LHS,data){
    xy <- sortedXyData(mCall[["x"]],LHS,data)
    T_h <- max(xy[["x"]][which.max(xy[["y"]])])
    tmp.w <- which(xy[["x"]] <= T_h)
    if (length(tmp.w)==1) {tmp.w <- 1:2}
    m <- lm(log(rate[tmp.w]) ~ I(1 / (k * (xy[["x"]][tmp.w]))))
    E <- (abs(summary(m)$coefficients[2, 1]))
    T.ref <- 283.15
    if (min(xy[["x"]]) > T.ref){
    B0 <- min(xy[["y"]][1]) } else {
    B0 <- max(xy[["y"]][which(xy[["x"]] <= T.ref)])}

    E_D <- 4*E
    structure(list(B0, E, E_D, T_h),
                          names=as.character(mCall[c("B0",
                                "E", "E_D", "T_h")]))
}




Self_starter<-selfStart(Schoolfield,School.Init,c('B0','E','E_D','T_h'))

bsrData <- data.frame(
    x = current_dataset$K,
    y = current_dataset$OriginalTraitValue)

schoolfield_nls <- NA
    try( 
        schoolfield_nls<- nlsLM(
            log(y) ~ Self_starter(x, B0, E, E_D, T_h), 
            
            data = bsrData
           
            )
	)

