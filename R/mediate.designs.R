
#source("C:/Users/dtingley/Documents/mediation/software/setupBoundsData.R")


#WRAPPER FUNCTIONS

#single experiment design
mediate.sed<-function(outcome,mediator,treatment,encouragement=NULL,SI=TRUE) {

    if(SI) {
    out<-mediate.np()#this will be luke's function
    } else {
    out<-mechanism.bounds(outcome,mediator,treatment,encouragement=NULL,design=SED)
    }

#I'm not sure how to return things correctly so things play nicely with the summary functions. right now we return "out" from mechanism.bounds function. I'm forgetting how to make this sort of wrapper thing work
return(out)
}



#parallel design
mediate.pd<-function(outcome,mediator,treatment,encouragement,NINT=TRUE,conf.level=.95) {

    if(NINT) {     
     out<-boot.pd(outcome,mediator,treatment,encouragement,NINT,conf.level)
    } else {
    out<-mechanism.bounds(outcome,mediator,treatment,encouragement,design="PD")
    }
    
return(out)
}


#parallel encouragement design
mediate.ped<-function(outcome,mediator,treatment,encouragement=NULL) {

    out<-mechanism.bounds(outcome,mediator,treatment,encouragement=NULL,design="PED")
    return(out)
    
}




#WORKHORSE FUNCTIONS

#parallel design under no interaction assumption
    boot.pd<-function(outcome,mediator, treatment,encouragement,sims=100, conf.level=.95) {
    
    n.o <- length(outcome)
    n.t <- length(treatment)
    n.m <- length(mediator)
    n.z<-length(encouragement)
   

        if(n.o != n.t | n.t != n.m |n.m !=n.z){
            stop("Error: Number of observations not the same in treatment, outcome, mediator, and encouragement data")
        }
    orig.length<-length(outcome)
    data<-matrix(,nrow=length(outcome),ncol=3)

    data[,1]<-outcome
    data[,2]<-treatment
    data[,3]<-mediator
    data[,4]<-encouragement
    data<-as.data.frame(data)
    names(data)<-c("Y","T","M","D")
    
    
    acme<- matrix(NA, sims, 1)
        # Bootstrap function.
                for(b in 1:sims){
                    index <- sample(1:n, n, replace = TRUE)
                    d<-data[index,]
            
            d0.temp<-mean(d$Y[D$T==1 & d$D==0) - mean(d$Y[d$T==0 & d$D==0)
            
            weight.m1<-sum(d$M==1 & d$D==1 )/sum(d$D==1)
            weight.m0<-sum(d$M==0 & d$D==1 )/sum(d$D==1)
            m1<-mean(d$Y[d$T==1 & d$M==1 & d$D==1) - mean(d$Y[d$T==0 & d$M==1 & d$D==1)
            m0<-mean(d$Y[d$T==1 & d$M==0 & d$D==1) - mean(d$Y[d$T==0 & d$M==0 & d$D==1)
            acme[b]<-weight.m1*m1 + weight.m0*m0
            }

        acme[acme==-Inf]<-NA
        acme[acme==Inf]<-NA
        
        low <- (1 - conf.level)/2
        high <- 1 - low

        acme.mu<-median(acme, na.rm=TRUE)
        acme.ci <- quantile(acme.mu, c(low,high), na.rm=TRUE)
        
        out<-list(acme.mu=acme.mu, acme.ci=acme.ci,acme=acme,conf.level=conf.level,sims=sims)
        out

    }
    
    
    

#Bounds Function for parallel, parallel encouragement, and single experiment designs
mechanism.bounds<-function(outcome,mediator,treatment,encouragement=NULL,design=NULL) {

    n.o <- length(outcome)
    n.t <- length(treatment)
    n.m <- length(mediator)

    if(is.null(design)){
        stop("Error: Design not indicated.")
    }

    if(design=="SED") {
        if(n.o != n.t | n.t != n.m){
            stop("Error: Number of non-missing observations not the same in treatment, outcome, and mediator data")
        }
    }

    if(design=="PED"|design=="PD") {
        n.z<-length(encouragement)
        if(n.o != n.t | n.t != n.m |n.m !=n.z){
            stop("Error: Number of observations not the same in treatment, outcome, and mediator data")
        }
        if(is.null(encouragement)){
            stop("Error: Design with mediator manipulation chosen but no indicator for manipulation status chosen.")
        }
    }


    orig.length<-length(outcome)
    d<-matrix(,nrow=length(outcome),ncol=3)

    d[,1]<-outcome
    d[,2]<-treatment
    d[,3]<-mediator
    d<-as.data.frame(d)
    names(d)<-c("Y","T","M")
    #encouragement indicates the variable that assigns manipulation direction. Z in the JRSSA paper. This is left to null for the single experiment design.
    manip<-!is.null(encouragement)


        #type of design used
        design<-design

            if(manip){
                d$Z<-d$D<-encouragement
                d<-na.omit(d)
                nobs<-nrow(d)
                if(nobs<orig.length) {
                    warning("NA's in data. Observations with missing data removed.")
                    }
                }

                #Single Experiment Design
                if(design=="SED") {
                d$Z<-0
                d<-na.omit(d)
                nobs<-nrow(d)
                if(nobs<orig.length) {
                    warning("NA's in data. Observations with missing data removed.")
                }
                P111<-sum(d$Y==1 & d$M==1 & d$Z==0 & d$T==1)/sum(d$T==1 & d$Z==0)
                P101<-sum(d$Y==1 & d$M==0 & d$Z==0 & d$T==1)/sum(d$T==1 & d$Z==0)
                P001<-sum(d$Y==0 & d$M==0 & d$Z==0 & d$T==1)/sum(d$T==1 & d$Z==0)
                P011<-sum(d$Y==0 & d$M==1 & d$Z==0 & d$T==1)/sum(d$T==1 & d$Z==0)

                P110<-sum(d$Y==1 & d$M==1 & d$Z==0 & d$T==0)/sum(d$T==0 & d$Z==0)
                P100<-sum(d$Y==1 & d$M==0 & d$Z==0 & d$T==0)/sum(d$T==0 & d$Z==0)
                P000<-sum(d$Y==0 & d$M==0 & d$Z==0 & d$T==0)/sum(d$T==0 & d$Z==0)
                P010<-sum(d$Y==0 & d$M==1 & d$Z==0 & d$T==0)/sum(d$T==0 & d$Z==0)
                #Calculate Bounds Expressions
                #Single Experiment
                l.delta.t <- max(-P001-P011, -P000-P001-P100, -P011-P010-P110)
                u.delta.t <- min(P101+P111, P000+P100+P101, P010+P110+P111)
                l.delta.c <- max(-P100-P110, -P001-P101-P100, -P011-P111-P110)
                u.delta.c <- min(P000+P010, P011+P111+P010, P000+P001+P101)
                S.t <- cbind(l.delta.t,u.delta.t)
                S.c <- cbind(l.delta.c,u.delta.c)
                }


            #Parallel Design
            if(design=="PD") {
                #D indicates if there was manipulation
                #Calculate Quantities necessary for computation of bounds
                Z111 <- sum(d$Y==1 & d$T==1 & d$M==1 & d$D==1)/sum(d$T==1 & d$M==1 & d$D==1)
                Z011 <- sum(d$Y==0 & d$T==1 & d$M==1 & d$D==1)/sum(d$T==1 & d$M==1 & d$D==1)
                Z001 <- sum(d$Y==0 & d$T==1 & d$M==0 & d$D==1)/sum(d$T==1 & d$M==0 & d$D==1)
                Z110 <- sum(d$Y==1 & d$T==0 & d$M==1 & d$D==1)/sum(d$T==0 & d$M==1 & d$D==1)
                Z010 <- sum(d$Y==0 & d$T==0 & d$M==1 & d$D==1)/sum(d$T==0 & d$M==1 & d$D==1)
                Z000 <- sum(d$Y==0 & d$T==0 & d$M==0 & d$D==1)/sum(d$T==0 & d$M==0 & d$D==1)
                Z101 <- sum(d$Y==1 & d$T==1 & d$M==0 & d$D==1)/sum(d$T==1 & d$M==0 & d$D==1)
                Z100 <- sum(d$Y==1 & d$T==0 & d$M==0 & d$D==1)/sum(d$T==0 & d$M==0 & d$D==1)

                P111 <- sum(d$Y==1 & d$M==1 & d$T==1 & d$D==0)/sum(d$T==1 & d$D==0)
                P101 <- sum(d$Y==1 & d$M==0 & d$T==1 & d$D==0)/sum(d$T==1 & d$D==0)
                P001 <- sum(d$Y==0 & d$M==0 & d$T==1 & d$D==0)/sum(d$T==1 & d$D==0)
                P011 <- sum(d$Y==0 & d$M==1 & d$T==1 & d$D==0)/sum(d$T==1 & d$D==0)
                P110 <- sum(d$Y==1 & d$M==1 & d$T==0 & d$D==0)/sum(d$T==0 & d$D==0)
                P100 <- sum(d$Y==1 & d$M==0 & d$T==0 & d$D==0)/sum(d$T==0 & d$D==0)
                P000 <- sum(d$Y==0 & d$M==0 & d$T==0 & d$D==0)/sum(d$T==0 & d$D==0)
                P010 <- sum(d$Y==0 & d$M==1 & d$T==0 & d$D==0)/sum(d$T==0 & d$D==0)

                #Parallel Design
                P.l.delta.t <- max(-P001-P011, -P011-P010-P110-P001+Z001, -P000-P001-P100-P011+Z011, -P001-P011+Z001-Z111, -P001+P101-Z101, -P011+P111-Z111)
                P.u.delta.t <- min(P101+P111, P010+P110+P101+P111-Z101, P000+P100+P101+P111-Z111, P101+P111+Z001-Z111, P111-P011+Z011, P101-P001+Z001)
                P.l.delta.c <- max(-P100-P110, -P011-P111-P110-P100+Z000, -P001-P101-P100-P110+Z110,-P100-P110+Z100-Z010, -P110+P010-Z010, -P100+P000-Z100)
                P.u.delta.c <- min(P000+P010, P011+P111+P010+P000-Z100, P000+P001+P101+P010-Z010, P000+P010+Z100-Z010, P010-P110+Z110, P000-P100+Z100)
                P.t <- cbind(P.l.delta.t,P.u.delta.t)
                P.c <- cbind(P.l.delta.c,P.u.delta.c)
            }


            #Parallel Encouragement design
            if(design=="PED") {
                    #Calculate Quantities necessary for computation of bounds
                    dP000 <- sum(d$Y==0 & d$M==0 & d$T==0 & d$Z==-1)/sum(d$T==0 & d$Z==-1)
                    dP001 <- sum(d$Y==0 & d$M==0 & d$T==1 & d$Z==-1)/sum(d$T==1 & d$Z==-1)
                    dP010 <- sum(d$Y==0 & d$M==1 & d$T==0 & d$Z==-1)/sum(d$T==0 & d$Z==-1)
                    dP011 <- sum(d$Y==0 & d$M==1 & d$T==1 & d$Z==-1)/sum(d$T==1 & d$Z==-1)
                    dP100 <- sum(d$Y==1 & d$M==0 & d$T==0 & d$Z==-1)/sum(d$T==0 & d$Z==-1)
                    dP101 <- sum(d$Y==1 & d$M==0 & d$T==1 & d$Z==-1)/sum(d$T==1 & d$Z==-1)
                    dP110 <- sum(d$Y==1 & d$M==1 & d$T==0 & d$Z==-1)/sum(d$T==0 & d$Z==-1)
                    dP111 <- sum(d$Y==1 & d$M==1 & d$T==1 & d$Z==-1)/sum(d$T==1 & d$Z==-1)

                   tP000 <- sum(d$Y==0 & d$M==0 & d$T==0 & d$Z==0)/sum(d$T==0 & d$Z==0)
                   tP001 <- sum(d$Y==0 & d$M==0 & d$T==1 & d$Z==0)/sum(d$T==1 & d$Z==0)
                   tP010 <- sum(d$Y==0 & d$M==1 & d$T==0 & d$Z==0)/sum(d$T==0 & d$Z==0)
                   tP011 <- sum(d$Y==0 & d$M==1 & d$T==1 & d$Z==0)/sum(d$T==1 & d$Z==0)
                   tP100 <- sum(d$Y==1 & d$M==0 & d$T==0 & d$Z==0)/sum(d$T==0 & d$Z==0)
                   tP101 <- sum(d$Y==1 & d$M==0 & d$T==1 & d$Z==0)/sum(d$T==1 & d$Z==0)
                   tP110 <- sum(d$Y==1 & d$M==1 & d$T==0 & d$Z==0)/sum(d$T==0 & d$Z==0)
                   tP111 <- sum(d$Y==1 & d$M==1 & d$T==1 & d$Z==0)/sum(d$T==1 & d$Z==0)

                   sP000 <- sum(d$Y==0 & d$M==0 & d$T==0 & d$Z==1)/sum(d$T==0 & d$Z==1)
                   sP001 <- sum(d$Y==0 & d$M==0 & d$T==1 & d$Z==1)/sum(d$T==1 & d$Z==1)
                   sP010 <- sum(d$Y==0 & d$M==1 & d$T==0 & d$Z==1)/sum(d$T==0 & d$Z==1)
                   sP011 <- sum(d$Y==0 & d$M==1 & d$T==1 & d$Z==1)/sum(d$T==1 & d$Z==1)
                   sP100 <- sum(d$Y==1 & d$M==0 & d$T==0 & d$Z==1)/sum(d$T==0 & d$Z==1)
                   sP101 <- sum(d$Y==1 & d$M==0 & d$T==1 & d$Z==1)/sum(d$T==1 & d$Z==1)
                   sP110 <- sum(d$Y==1 & d$M==1 & d$T==0 & d$Z==1)/sum(d$T==0 & d$Z==1)
                   sP111 <- sum(d$Y==1 & d$M==1 & d$T==1 & d$Z==1)/sum(d$T==1 & d$Z==1)


                    P0 <- c(dP000,dP010,dP100,dP110,tP000,tP010,tP100,tP110,sP000,sP010,sP100,sP110)
                    P1 <- c(dP001,dP011,dP101,dP111,tP001,tP011,tP101,tP111,sP001,sP011,sP101,sP111)
                    Q0 <- c(dP010+dP110,tP010+tP110,sP010+sP110)
                    Q1 <- c(dP011+dP111,tP011+tP111,sP011+sP111)


                #Linear programming functions for PED
                bounds.bidir.cmpl <- function(P,Q,dir){
                    f.obj <- c(rep(0,16), 0,0,1,0, 0,0,1,0, 0,-1,0,0, 0,-1,0,0,
                            0,0,-1,0, 0,0,-1,0, 0,1,0,0, 0,1,0,0, rep(0,16))
                    f.con <- matrix(c(
                            rep(c(1,1,1,0),8), rep(0,32), #dP00t
                            rep(c(rep(c(0,0,0,1),4), rep(0,16)),2), #dP01t
                            rep(0,32), rep(c(1,1,1,0),8), #dP10t
                            rep(c(rep(0,16), rep(c(0,0,0,1),4)),2), #dP11t
                            rep(c(1,1,0,0),8), rep(0,32), #tP00t
                            rep(c(rep(c(0,0,1,1),4), rep(0,16)),2), #tP01t
                            rep(0,32), rep(c(1,1,0,0),8), #tP10t
                            rep(c(rep(0,16), rep(c(0,0,1,1),4)),2), #tP11t
                            rep(c(1,0,0,0),8), rep(0,32), #sP00t
                            rep(c(rep(c(0,1,1,1),4), rep(0,16)),2), #sP01t
                            rep(0,32), rep(c(1,0,0,0),8), #sP10t
                            rep(c(rep(0,16), rep(c(0,1,1,1),4)),2), #sP11t
                            rep(c(rep(0,12), rep(1,4)),4), #dQ
                            rep(c(rep(0,8), rep(1,8)),4), #tQ
                            rep(c(rep(0,4), rep(1,12)),4), #sQ
                            rep(1,64) #sum=1
                            ), nrow=16, byrow=T)
                    f.dir <- rep("=", 16)
                    f.rhs <- c(P,Q,1)
                    r <- lp(dir, f.obj, f.con, f.dir, f.rhs)
                    r$objval
                }

                bounds.bidir.popl <- function(P,Q,dir){
                    f.obj <- c(rep(0,16), 0,0,1,1, 0,0,1,1, -1,-1,0,0, -1,-1,0,0,
                            0,0,-1,-1, 0,0,-1,-1, 1,1,0,0, 1,1,0,0, rep(0,16))
                    f.con <- matrix(c(
                            rep(c(1,1,1,0),8), rep(0,32), #dP00t
                            rep(c(rep(c(0,0,0,1),4), rep(0,16)),2), #dP01t
                            rep(0,32), rep(c(1,1,1,0),8), #dP10t
                            rep(c(rep(0,16), rep(c(0,0,0,1),4)),2), #dP11t
                            rep(c(1,1,0,0),8), rep(0,32), #tP00t
                            rep(c(rep(c(0,0,1,1),4), rep(0,16)),2), #tP01t
                            rep(0,32), rep(c(1,1,0,0),8), #tP10t
                            rep(c(rep(0,16), rep(c(0,0,1,1),4)),2), #tP11t
                            rep(c(1,0,0,0),8), rep(0,32), #sP00t
                            rep(c(rep(c(0,1,1,1),4), rep(0,16)),2), #sP01t
                            rep(0,32), rep(c(1,0,0,0),8), #sP10t
                            rep(c(rep(0,16), rep(c(0,1,1,1),4)),2), #sP11t
                            rep(c(rep(0,12), rep(1,4)),4), #dQ
                            rep(c(rep(0,8), rep(1,8)),4), #tQ
                            rep(c(rep(0,4), rep(1,12)),4), #sQ
                            rep(1,64) #sum=1
                            ), nrow=16, byrow=T)
                    f.dir <- rep("=", 16)
                    f.rhs <- c(P,Q,1)
                    r <- lp(dir, f.obj, f.con, f.dir, f.rhs)
                    r$objval
                }

                BE.t <- c(bounds.bidir.popl(P1, Q0, "min"), bounds.bidir.popl(P1, Q0, "max"))
                BE.c <- c(-bounds.bidir.popl(P0 ,Q1, "max"), -bounds.bidir.popl(P0, Q1, "min"))
                num.BET.t.lo <- bounds.bidir.cmpl(P1, Q0, "min")
                num.BET.t.up <- bounds.bidir.cmpl(P1, Q0, "max")
                num.BET.c.lo <- -bounds.bidir.cmpl(P0, Q1, "max")
                num.BET.c.up <- -bounds.bidir.cmpl(P0, Q1, "min")
                denom.BET.t <- P1[1] + P1[3] - P1[9] - P1[11]
                denom.BET.c <- P0[1] + P0[3] - P0[9] - P0[11]
                BET.t <- c(num.BET.t.lo/denom.BET.t, num.BET.t.up/denom.BET.t)
                BET.c <- c(num.BET.c.lo/denom.BET.c, num.BET.c.up/denom.BET.c)

            }#end computation section

                #function output
                if(design=="SED"){
                out<-list(Single.t = S.t, Single.c = S.c, manip=manip,design=design, nobs=nobs)
                }
                if(design=="PD"){
                out<-list(Parallel.t = P.t, Parallel.c = P.c, manip=manip,design=design, nobs=nobs)
                }
                if(design=="PED") {
                out<-list(BE.t=BE.t, BE.c=BE.c,  BET.t=BET.t, BET.c=BET.c, manip=manip, design=design, nobs=nobs)
                }


    rm(d)
    class(out) <- "mechanism.bounds"
    out
}


summary.mechanism.bounds <- function(object, ...){
    structure(object, class = c("summary.mechanism.bounds", class(object)))
}



print.summary.mechanism.bounds <- function(x, ...){
    cat("\n Bounds for Causal Mediation Effects \n\n")
    cat("Design / Treatment Condition / Lower and Upper Bounds \n\n")

        # Print values
        if(x$design=="SED"){
        cat("Single Experiment Design, Treatment ", format(x$Single.t, digits=4), "\n")
        cat("Single Experiment Design, Control ", format(x$Single.c, digits=4), "\n")
            cat("\n\n")
    cat("Sample Size Used:", x$nobs,"\n\n")
        }

        if(x$design=="PD"){
            cat("Parallel Design, Treatment ", format(x$Parallel.t, digits=4), "\n")
            cat("Parallel Design, Control ", format(x$Parallel.c, digits=4), "\n")
            cat("\n\n")
            cat("Sample Size Used:", x$nobs,"\n\n")
            }

        if(x$design=="PED"){
            cat("Parallel Encouragment Design, Treatment, Population ", format(x$BE.t, digits=4), "\n")
            cat("Parallel Encouragment Design, Treatment, Complier ", format(x$BET.t, digits=4), "\n")
            cat("Parallel Encouragment Design, Control, Population ", format(x$BE.c, digits=4), "\n")
            cat("Parallel Encouragment Design, Control, Complier ", format(x$BET.c, digits=4), "\n")
            cat("\n\n")
            cat("Sample Size Used:", x$nobs,"\n\n")
            }


    invisible(x)
}
