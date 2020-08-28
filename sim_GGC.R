# simulation of G|G|c queue (with loss/no-loss)
# provides simple time-averaged steady state results
# version april 2020
# rich wood (richard dot wood16 at nhs dot net)

rm(list=ls())

require(parallel)
require(ggplot2)

####################################################################################################################

WUP<-100
DUR<-365
nruns<-10

arr_dist<-"pois"
arr_params<-c(11)
srv_dist<-"exp"
srv_params<-c(1)
cap<-15

loss<-TRUE

###################################################################################################################

start.time<-Sys.time()

print(paste("traffic intensity:",round(mean(do.call(paste0("r",arr_dist),c(list(n=1000000),arr_params)))/(cap*1/mean(do.call(paste0("r",srv_dist),c(list(n=1000000),srv_params)))),2)),quote=FALSE)

simfn<-function(runs) {
  set.seed(runs)
  #get num arrivals by day
  day_arr_times<-round(do.call(paste0("r",arr_dist),c(list(n=WUP+DUR),arr_params)))
  arr_neg<-sum(day_arr_times<0)/length(day_arr_times)
  day_arr_times[which(day_arr_times<0)]<-0
  arr_times<-unlist(sapply(1:length(day_arr_times), function(x) {
    sort(runif(day_arr_times[x],0,1)+x-1)
  }))
  cal<-data.frame(id=1:length(arr_times),time=arr_times,event="arrival")
  tx<-0
  res<-data.frame(time=1:(WUP+DUR),occ=NA,niq=NA,arr_admit=0,arr_no_admit=0)
  niq<-0  #number in queue
  occ<-0  #occupancy (number in unit)
  while (tx<=(WUP+DUR) & nrow(cal)>0) {
    ind1<-which(cal$time>tx & cal$event %in% c("arrival","endsrv"))
    ind<-ind1[which.min(cal$time[ind1])]
    niq_old<-niq
    occ_old<-occ
    tx_old<-tx
    tx<-cal$time[ind]
    if (tx>(WUP+DUR) | nrow(cal)==0) break
    tx_day<-ceiling(tx)
    if (cal$event[ind]=="arrival") {
      if (occ<cap) {
        res$arr_admit[tx_day]<-res$arr_admit[tx_day]+1
        #admit patient
        cal<-rbind(cal,data.frame(id=cal$id[ind],time=tx,event="startsrv"))
        los<-do.call(paste0("r",srv_dist),c(list(n=1),srv_params))
        cal<-rbind(cal,data.frame(id=cal$id[ind],time=tx+los,event="endsrv"))
        occ<-occ+1
      } else {
        res$arr_no_admit[tx_day]<-res$arr_no_admit[tx_day]+1
        if (loss==FALSE) {
          #patient wait in queue
          niq<-niq+1
        } else {
          cal<-cal[-which(cal$id==cal$id[ind]),] 
        }
      }
    } else if (cal$event[ind]=="endsrv") {
      cal<-cal[-which(cal$id==cal$id[ind]),] 
      if (niq==0) {
        occ<-occ-1
      } else {
        #admit patient (backfill bed)
        los<-do.call(paste0("r",srv_dist),c(list(n=1),srv_params))
        #select patient who's been waiting longest (and has not started/finished service)
        poss_ids<-setdiff(unique(cal$id),cal$id[which(cal$event=="startsrv")])
        waits<-data.frame(id=poss_ids,waits=cal$time[which(cal$id %in% poss_ids)]-tx)
        admit_id<-waits$id[which.min(waits$waits)]
        cal<-rbind(cal,data.frame(id=admit_id,time=tx,event="startsrv"))
        cal<-rbind(cal,data.frame(id=admit_id,time=tx+los,event="endsrv"))
        niq<-niq-1
      }
    }
    cal<-cal[order(cal$time),]
    #save results, extract performance measures
    wt_new<-(tx-tx_old)/tx
    res$niq[tx_day]<-ifelse(is.na(res$niq[tx_day]),(tx-floor(tx))*niq_old+(ceiling(tx)-tx)*niq,wt_new*niq+(1-wt_new)*res$niq[tx_day])
    res$occ[tx_day]<-ifelse(is.na(res$occ[tx_day]),(tx-floor(tx))*occ_old+(ceiling(tx)-tx)*occ,wt_new*occ+(1-wt_new)*res$occ[tx_day])
  }
  res<-res %>%
    mutate(niq=ifelse(time==1 & is.na(niq),0,niq)) %>%
    mutate(occ=ifelse(time==1 & is.na(occ),0,occ)) %>%
    fill(niq) %>%
    fill(occ)
  return(list(res,arr_neg))
}

cl<-makeCluster(detectCores()-1)
clusterExport(cl=cl,varlist=c("WUP","DUR","arr_dist","arr_params","srv_dist","srv_params","cap","loss"))
clusterEvalQ(cl=cl,c(library(tidyr),library(dplyr)))
RES<-parLapply(cl,1:nruns,simfn)
stopCluster(cl)

###################################################################################################################

#number of negative arrivals that needed to be bounded
avg_arr_neg<-mean(unlist(lapply(RES,function(x) x<-x[[2]])))
print(paste0("average number of negative arrival days: ",round(avg_arr_neg*100,2),"%"),quote=FALSE)

#graphical output for single run
run<-1
print(ggplot(RES[[run]][[1]],aes(x=time)) +
  geom_line(aes(y=occ),colour="red") +
  geom_line(aes(y=niq),colour="blue") +
  geom_vline(xintercept=WUP,linetype="dashed",colour="darkgrey") +
  xlab("Time") +
  ylab("Number") +
  ggtitle("Occupancy (red) and queuing (blue)"))

#get performance measures for each run
agg_RES<-lapply(RES,function(x) {
  x<-x[[1]]
  x<-x[x$time>WUP,]
  avg_occ<-mean(x$occ)
  avg_occ_prop<-avg_occ/cap
  avg_niq<-mean(x$niq)
  avg_wait<-avg_niq/mean(do.call(paste0("r",arr_dist),c(list(n=1000000),arr_params)))
  arr_no_admit<-sum(x$arr_no_admit)
  arr_no_admit_prop<-arr_no_admit/(arr_no_admit+sum(x$arr_admit))
  crisis_days<-length(which(x$arr_no_admit>0))
  crisis_days_prop<-crisis_days/nrow(x)
  return(data.frame(avg_occ=avg_occ,avg_occ_prop=avg_occ_prop,avg_niq=avg_niq,avg_wait=avg_wait,arr_no_admit=arr_no_admit,arr_no_admit_prop=arr_no_admit_prop,crisis_days=crisis_days,crisis_days_prop=crisis_days_prop))
})
agg_RES<-do.call("rbind",agg_RES)

#average over runs
sum_RES<-signif(colMeans(agg_RES),4)
print(sum_RES)
print(paste("note: a crisis day is a day in which at least one arriving patient cannot be immediately admitted"),quote=FALSE)

#processing time
print(difftime(Sys.time(),start.time),quote=FALSE)












