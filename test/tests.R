
# gee support @(#) tests.q 1.2 96/09/27

library(gee)

source("testgee.R")

myr <- gee(gsind~x2+x3+x4,id=id,data=testgee)
myr
summary(myr)
summary(gee(gs~x2+x4,id=id,data=testgee,corstr="exchangeable"))

summary(gee(gsind~sin(x2)+log(x3)+x4,id=id,data=testgee))
#summary(gee(gsind~gam(x2~x1)$fitted.values+x4,id=id,data=testgee,subset=id<20))

summary(gee(lgind~x2+x3+x4,id=id,data=testgee,family=binomial))
glm(lgind~x2+x3+x4,data=testgee, start=rep(0,4), # needed in R
    family=quasi(link=logit,variance="mu(1-mu)"))

summary(gee(lgind~x2+x3+x4,id=id,data=testgee,family=binomial,scale.fix=T))
#glm(lgind~x2+x3+x4,data=testgee,family=binomial)

summary(gee(lg~x2+x3+x4,id=id,data=testgee,family=binomial,corstr="exchangeable"))
summary(gee(lg~x2+x3+x4,id=id,data=testgee,family=binomial,corstr="stat_M_dep",M=1))

summary(gee(lg~x2+x3+x4,id=id,data=testgee,family=binomial(link=probit)))
# next line does not work in R's glm: needs starting values for quasi
#summary(gee(lg~x2+x3+x4,id=id,data=testgee,family=quasi(link=probit,variance="mu(1-mu)")))
summary(gee(lg~x2+x3+x4, id=id, data=testgee, b=rep(0,4), family=quasi(link=probit,variance="mu(1-mu)")))

logmat <- cbind(testgee[,"lg1"],testgee[,"nn"])
xmat   <- cbind(testgee[,"x2"], testgee[,"x3"])
summary(gee(logmat~xmat,id = testgee[,"id"],family=binomial(link=probit)))

sub<-rep(c(F,F,F,F,T),100)
summary(gee(lg~x2+x4,id=id,data=testgee,family=binomial,subset=sub))

R <- c(1.0,0.2,0.3,0.1,
       0.2,1.0,0.2,0.3,
       0.3,0.2,1.0,0.2,
       0.1,0.3,0.2,1.0)
R <- matrix(R,nrow=4)
summary(gee(lg1~x2+x4,id=id,data=testgee,family=poisson(link=log),
             subset=sub,corstr="fixed",R=R))
summary(gee(lg1~x2+x4,id=id,data=testgee,family=Gamma(link=log),
             subset=sub,corstr="fixed",R=R))
