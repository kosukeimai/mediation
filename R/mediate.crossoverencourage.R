
crossoverencourage<-function(outcome=outcome,mediator1=mediator1,mediator2=mediator2,treatment=treatment,encouragement=encouragement,sims=100, conf.level=.95) {

        data<-matrix(,nrow=length(outcome),ncol=5)

        data[,1]<-outcome
        data[,2]<-treatment
        data[,3]<-mediator1
        data[,4]<-mediator2
        data[,5]<-encouragement
        data<-as.data.frame(data)
        names(data)<-c("Y2","T","M1","M2","V")
        d<-data
        d<-na.omit(d)
        n <- length(d$Y2)

        # Storage
            d.p.c<-d.p.t <- matrix(NA, sims, 1)
        # Bootstrap function.
        for(b in 1:sims){
            index <- sample(1:n, n, replace = TRUE)
            d<-data[index,]

                #Atmv
                A111=sum(d$M1==1 & d$T==1 & d$M2==1 & d$V==1)/sum(d$T==1 & d$M2==1 & d$V==1)
                A011=sum(d$M1==1 & d$T==0 & d$M2==1 & d$V==1)/sum(d$T==0 & d$M2==1 & d$V==1)
                A101=sum(d$M1==1 & d$T==1 & d$M2==0 & d$V==1)/sum(d$T==1 & d$M2==0 & d$V==1)
                A001=sum(d$M1==1 & d$T==0 & d$M2==0 & d$V==1)/sum(d$T==0 & d$M2==0 & d$V==1)

                A110=sum(d$M1==1 & d$T==1 & d$M2==1 & d$V==0)/sum(d$T==1 & d$M2==1 & d$V==0)
                A010=sum(d$M1==1 & d$T==0 & d$M2==1 & d$V==0)/sum(d$T==0 & d$M2==1 & d$V==0)
                A100=sum(d$M1==1 & d$T==1 & d$M2==0 & d$V==0)/sum(d$T==1 & d$M2==0 & d$V==0)
                A000=sum(d$M1==1 & d$T==0 & d$M2==0 & d$V==0)/sum(d$T==0 & d$M2==0 & d$V==0)

                #Gtm1vm2
                G1111=mean(d$Y2[d$T==1 & d$M1==1 & d$V==1 & d$M2==1])
                G0111=mean(d$Y2[d$T==0 & d$M1==1 & d$V==1 & d$M2==1])
                G1011=mean(d$Y2[d$T==1 & d$M1==0 & d$V==1 & d$M2==1])
                G0011=mean(d$Y2[d$T==0 & d$M1==0 & d$V==1 & d$M2==1])

                G1101=mean(d$Y2[d$T==1 & d$M1==1 & d$V==0 & d$M2==1])
                G0101=mean(d$Y2[d$T==0 & d$M1==1 & d$V==0 & d$M2==1])
                G1001=mean(d$Y2[d$T==1 & d$M1==0 & d$V==0 & d$M2==1])
                G0001=mean(d$Y2[d$T==0 & d$M1==0 & d$V==0 & d$M2==1])

                G1110=mean(d$Y2[d$T==1 & d$M1==1 & d$V==1 & d$M2==0])
                G0110=mean(d$Y2[d$T==0 & d$M1==1 & d$V==1 & d$M2==0])
                G1010=mean(d$Y2[d$T==1 & d$M1==0 & d$V==1 & d$M2==0])
                G0010=mean(d$Y2[d$T==0 & d$M1==0 & d$V==1 & d$M2==0])

                G1100=mean(d$Y2[d$T==1 & d$M1==1 & d$V==0 & d$M2==0])
                G0100=mean(d$Y2[d$T==0 & d$M1==1 & d$V==0 & d$M2==0])
                G1000=mean(d$Y2[d$T==1 & d$M1==0 & d$V==0 & d$M2==0])
                G0000=mean(d$Y2[d$T==0 & d$M1==0 & d$V==0 & d$M2==0])

                temp11<-sum(d$T==1 & d$M1==1)/sum(d$T==1|d$T==0)
                temp01<-sum(d$T==0 & d$M1==1)/sum(d$T==1|d$T==0)
                temp10<-sum(d$T==1 & d$M1==0)/sum(d$T==1|d$T==0)
                temp00<-sum(d$T==0 & d$M1==0)/sum(d$T==1|d$T==0)


                PH1<-((A111-A110)*temp11)/((A111-A110)*temp11+(1-2*A101)*temp10)
                PH0<-((A011-A010)*temp01)/((A011-A010)*temp01+(1-2*A001)*temp00)

                d.p.c[b]<-(PH1/(A111-A110))*(G1111+G1110-G1100-G1110*A111-G1101*A110)+((1-PH1)/(A100-A101))*(G1010-G1001-G1000+G1000*A100+G1011*A101)
                d.p.t[b]<-(PH0/(A010-A011))*(G0111+G0110-G0100-G0110*A011-G0101*A010)+((1-PH0)/(A001-A000))*(G0010-G0001-G0000+G0000*A000+G0011*A001)

        }#bootstraploop
        #replace any instances of infinity
        d.p.c[d.p.c==-Inf]<-NA
        d.p.t[d.p.t==-Inf]<-NA
        d.p.c[d.p.c==Inf]<-NA
        d.p.t[d.p.t==Inf]<-NA

        low <- (1 - conf.level)/2
        high <- 1 - low
        #We probably want to take the median rather than the mean so outlier draws have less influence?
        d.p.c.mu<-mean(d.p.c, na.rm=TRUE)
        d.p.t.mu<-mean(d.p.t, na.rm=TRUE)
        d.p.c.ci <- quantile(d.p.c, c(low,high), na.rm=TRUE)
        d.p.t.ci <- quantile(d.p.t, c(low,high), na.rm=TRUE)

        out<-list(d.p.c.mu = d.p.c.mu, d.p.t.mu = d.p.t.mu,d.p.c.ci = d.p.c.ci, d.p.t.ci = d.p.t.ci,d.p.c=d.p.c,d.p.t=d.p.t,conf.level=conf.level)

        class(out) <- "crossoverencourage"
        out
    }



summary.crossoverencourage <- function(object, ...){
    structure(object, class = c("summary.crossoverencourage", class(object)))
}


print.summary.crossoverencourage <- function(x, ...){
    clp <- 100 * x$conf.level
    cat("\n Crossover Encouragment Design \n\n")
                cat("Mediation Effect_0: ", format(x$d.p.c.mu, digits=4), clp, "% CI ",
                format(x$d.p.c.ci, digits=4), "\n")
                cat("Mediation Effect_1: ", format(x$d.p.t.mu, digits=4), clp, "% CI ",
                format(x$d.p.t.ci, digits=4), "\n")
    invisible(x)
}
