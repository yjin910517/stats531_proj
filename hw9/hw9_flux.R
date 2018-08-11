require(pomp)
require(doParallel)
require(foreach)


bsflu_data<-read.table("https://ionides.github.io/531w18/12/bsflu_data.txt")

bsflu_statenames <- c("S","I","R1","R2")
bsflu_paramnames <- c("Beta","mu_I","rho","mu_R1","mu_R2")
(bsflu_obsnames <- colnames(bsflu_data)[1:2])


bsflu_dmeasure <- "
  lik = dpois(B,rho*R1+1e-6,give_log);
"

bsflu_rmeasure <- "
  B = rpois(rho*R1+1e-6);
  C = rpois(rho*R2);
"

bsflu_rprocess <- "
  double t1 = rbinom(S,1-exp(-Beta*I*dt));
  double t2 = rbinom(I,1-exp(-dt*mu_I));
  double t3 = rbinom(R1,1-exp(-dt*mu_R1));
  double t4 = rbinom(R2,1-exp(-dt*mu_R2));
  S -= t1;
  I += t1 - t2;
  R1 += t2 - t3;
  R2 += t3 - t4;
"

bsflu_fromEstimationScale <- "
 TBeta = exp(Beta);
 Tmu_I = exp(mu_I);
 Trho = expit(rho);
"

bsflu_toEstimationScale <- "
 TBeta = log(Beta);
 Tmu_I = log(mu_I);
 Trho = logit(rho);
"

bsflu_initializer <- "
 S=762;
 I=1;
 R1=0;
 R2=0;
"

bsflu2 <- pomp(
  data=bsflu_data,
  times="day",
  t0=0,
  rprocess=euler.sim(
    step.fun=Csnippet(bsflu_rprocess),
    delta.t=1/12
  ),
  rmeasure=Csnippet(bsflu_rmeasure),
  dmeasure=Csnippet(bsflu_dmeasure),
  fromEstimationScale=Csnippet(bsflu_fromEstimationScale),
  toEstimationScale=Csnippet(bsflu_toEstimationScale),
  obsnames = bsflu_obsnames,
  statenames=bsflu_statenames,
  paramnames=bsflu_paramnames,
  initializer=Csnippet(bsflu_initializer)
)


run_level <- 3
switch(run_level,
       {bsflu_Np=100; bsflu_Nmif=10; bsflu_Neval=10; bsflu_Nglobal=10; bsflu_Nlocal=10}, 
       {bsflu_Np=20000; bsflu_Nmif=100; bsflu_Neval=10; bsflu_Nglobal=10; bsflu_Nlocal=10}, 
       {bsflu_Np=60000; bsflu_Nmif=300; bsflu_Neval=10; bsflu_Nglobal=100; bsflu_Nlocal=20},
       {bsflu_Np=50000; bsflu_Nmif=220; bsflu_Neval=10; bsflu_Nglobal=100; bsflu_Nlocal=20}
)



bsflu_params <- data.matrix(read.table("https://ionides.github.io/531w18/12/mif_bsflu_params.csv",row.names=NULL,header=TRUE))
bsflu_mle <- bsflu_params[which.max(bsflu_params[,"logLik"]),][bsflu_paramnames]
bsflu_fixed_params <- c(mu_R1=1/(sum(bsflu_data$B)/512),mu_R2=1/(sum(bsflu_data$C)/512))


cores <- 20  # The number of cores on this machine 
registerDoParallel(cores)
mcopts <- list(set.seed=TRUE)

set.seed(396658101,kind="L'Ecuyer")


### Local Search
bsflu_rw.sd <- 0.02
bsflu_cooling.fraction.50 <- 0.5

stew(file=sprintf("local_search-%d.rda",run_level),{
  
  t_local <- system.time({
    mifs_local <- foreach(i=1:bsflu_Nlocal,.packages='pomp', .export = ls(globalenv()), .combine=c, .options.multicore=mcopts) %dopar%  {
      mif2(
        bsflu2,
        start=bsflu_mle,
        Np=bsflu_Np,
        Nmif=bsflu_Nmif,
        cooling.type="geometric",
        cooling.fraction.50=bsflu_cooling.fraction.50,
        transform=TRUE,
        rw.sd=rw.sd(
          Beta=bsflu_rw.sd,
          mu_I=bsflu_rw.sd,
          rho=bsflu_rw.sd
        )
      )
      
    }
  })
  
},seed=900242057,kind="L'Ecuyer")

stew(file=sprintf("lik_local-%d.rda",run_level),{
    t_local_eval <- system.time({
    liks_local <- foreach(i=1:bsflu_Nlocal,.packages='pomp',.export = ls(globalenv()),.combine=rbind) %dopar% {
      evals <- replicate(bsflu_Neval, logLik(pfilter(bsflu2,params=coef(mifs_local[[i]]),Np=bsflu_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=900242057,kind="L'Ecuyer")


results_local <- data.frame(logLik=liks_local[,1],logLik_se=liks_local[,2],t(sapply(mifs_local,coef)))
summary(results_local$logLik,digits=5)


### Global Search
bsflu_box <- rbind(
  Beta=c(0.001,0.01),
  mu_I=c(0.5,2),
  rho = c(0.5,1)
)

stew(file=sprintf("box_eval-%d.rda",run_level),{
  
  t_global <- system.time({
    mifs_global <- foreach(i=1:bsflu_Nglobal,.packages='pomp', .export = ls(globalenv()),.combine=c, .options.multicore=mcopts) %dopar%  mif2(
      mifs_local[[1]],
      start=c(apply(bsflu_box,1,function(x)runif(1,x[1],x[2])),bsflu_fixed_params)
    )
  })
},seed=1270401374,kind="L'Ecuyer")


stew(file=sprintf("lik_global_eval-%d.rda",run_level),{
  t_global_eval <- system.time({
    liks_global <- foreach(i=1:bsflu_Nglobal,.packages='pomp', .export = ls(globalenv()),.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(bsflu_Neval, logLik(pfilter(bsflu2,params=coef(mifs_global[[i]]),Np=bsflu_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")

results_global <- data.frame(logLik=liks_global[,1],logLik_se=liks_global[,2],t(sapply(mifs_global,coef)))


#### Q9.1
run_level <- 4


stew(file=sprintf("Mif-9.1-%d.rda",run_level),{
  
  t_global.2 <- system.time({
    mifs_global.2 <- foreach(i=1:bsflu_Nglobal,.packages='pomp', .export = ls(globalenv()),.combine=c, .options.multicore=mcopts) %dopar% 
      mif2(
        mifs_global[[1]],
        start=c(apply(bsflu_box,1,function(x)runif(1,x[1],x[2])),bsflu_fixed_params),
        Np=bsflu_Np,
        Nmif=bsflu_Nmif,
        cooling.fraction.50=bsflu_cooling.fraction.50        
      )
  })
},seed=1270401374,kind="L'Ecuyer")


stew(file=sprintf("lik-9.1-%d.rda",run_level),{
  t_global_eval.2 <- system.time({
    liks_global.2 <- foreach(i=1:bsflu_Nglobal,.packages='pomp',.export = ls(globalenv()),.combine=rbind, .options.multicore=mcopts) %dopar% {
      evals <- replicate(bsflu_Neval, logLik(pfilter(bsflu2,params=coef(mifs_global.2[[i]]),Np=bsflu_Np)))
      logmeanexp(evals, se=TRUE)
    }
  })
},seed=442141592,kind="L'Ecuyer")

results_global.2 <- data.frame(logLik=liks_global.2[,1],logLik_se=liks_global.2[,2],t(sapply(mifs_global.2,coef)))



#### Q9.3
stew(file=sprintf("profile beta-%d.rda",It),{
  
  t_global.4 <- system.time({
      prof.llh<- foreach(i=1:(It*nprof),.packages='pomp', .export = ls(globalenv()),.combine=rbind, .options.multicore=mcopts) %dopar%{
        # Find MLE
        mif2(
          mifs_global[[1]],
          start=c(unlist(profile.box[i,]),bsflu_fixed_params),
          Np=5000,Nmif=100,
          rw.sd=rw.sd(
            mu_I=bsflu_rw.sd,
            rho=bsflu_rw.sd
          )
        )->mifs_global.4
        # evaluate llh
        evals = replicate(10, logLik(pfilter(mifs_global.4,Np=10000)))
        ll=logmeanexp(evals, se=TRUE)        
        
        data.frame(as.list(coef(mifs_global.4)),
                   loglik = ll[1],
                   loglik.se = ll[2])
      }
  })
},seed=931129,kind="L'Ecuyer")



prof.llh %>% 
  ddply(~Beta,subset,rank(-loglik)<=10) %>%
  subset(select=bsflu_paramnames) -> pars


## mif2 again
stew(file=sprintf("profile beta-2-%d.rda",It),{
  
  t_global.5 <- system.time({
    prof.llh<- foreach(i=1:(nrow(pars)),.packages='pomp', .export = ls(globalenv()),.combine=rbind, .options.multicore=mcopts) %dopar%{
      # Find MLE
      mif2(
        mifs_global[[1]],
        start=unlist(pars[i,]),
        Np=5000,Nmif=50,
        rw.sd=rw.sd(
          mu_I=bsflu_rw.sd,
          rho=bsflu_rw.sd
        )
      )->mifs_global.5
      # evaluate llh 
      pf= replicate(10,pfilter(mifs_global.5,Np=5000))
      evals=sapply(pf,logLik)
      ll=logmeanexp(evals, se=TRUE)  
      nfail=sapply(pf,getElement,"nfail")
      
      data.frame(as.list(coef(mifs_global.5)),
                 loglik = ll[1],
                 loglik.se = ll[2],
                 nfail.max=max(nfail))
    }
  })
},seed=931129,kind="L'Ecuyer")