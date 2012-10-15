homol.search <-
function(  peaklist,iso_list,minmz=3,maxmz=60,minrt=0.05,maxrt=2,ppm=TRUE,mztol=3.5,
                  rttol=0.4,minlen=5,plotit=TRUE,restr=FALSE){
  
      ##########################################################################
      # (0) Issue warnings: check arguments ####################################
      if(ppm==TRUE & mztol>10){cat("Such a big mztol?")};
      if(length(plotit)>1){stop("plotit either TRUE or FALSE")};
      if(minlen<3){stop("invalid minlen argument. set minlen >= 3")};
      if(mztol<=0){warning("mztol should be >0!")};
      # (1) Function parameters ################################################
      this<-(c(iso_list[,3]-round(iso_list[,3],digits=0))/(iso_list[,3]+(1/(iso_list[,5])*12)));
      maxup<-max(this);
      maxdown<-min(this);
      peak2<-peaklist[order(peaklist[,1],decreasing=FALSE),];
      len<-length(peak2[,1]);
      dif<-c();
      ret<-c();
      from<-c();
      to<-c();
      getit<-seq(1,len,1);
      if(plotit==TRUE){
        #plot.new();
        sc<-close.screen();if(sc[1]!=FALSE){for(m in 1:length(sc)){close.screen(sc[m])}};
        def.par <- par(no.readonly = TRUE) # save default, for resetting...
        split.screen(c(2,1));
        make<-seq(1,len,20);
        screen(1);
        par(mar=c(4.5,4,2,0.3));
        plot.window(xlim=c(0,(length(make)+1)),ylim=c(-0.2,0.1))
        text(2,0.11,labels="Progress",pos=4,cex=1.2);
        text(2,0.03,labels="(1) Screen peak list for relevant m/z increments...",pos=4,cex=0.9);
        points(seq(1,length(make),1),rep(-0.02,length(make)),cex=1.2,col="darkgrey",pch=22,bg="darkgrey");
      }
      ##########################################################################
      # (2) screen for differences ############################################# 
      if(ppm==TRUE){
        delt<-(mztol*max(peak2[,1])/1e6)*2;
      }else{
        delt<-mztol;
      }
      neg2<-c(0);  # how many rejected in outer loop?
      neg1<-c(0);  # how many rejected in inner loop?
      within<-rep(FALSE,len);
      for(i in 1:len){
        if(plotit==TRUE){
          if(any(make==i)){
           points(seq(1,length(make),1)[make==i],-0.02,cex=1.2,col="green",pch=22,bg="green");
          }
        }
        # (2.1) Find all peaks within range ####################################
        til<-getit[
            (peak2[,1]-peak2[i,1])<=maxmz   &
            (peak2[,1]-peak2[i,1])>=minmz   &
            (peak2[,3]-peak2[i,3])>=minrt   &
            (peak2[,3]-peak2[i,3])<=maxrt   
        ];
        if(length(til)>0){
        # (2.2) check for feasible mass defects ################################
            that1<-c(peak2[til,1]-peak2[i,1]);
            m<-c(that1*maxup);   # upper limit
            n<-c(that1*maxdown); # lower limit
            this<-c(that1-round(that1,digits=0));   
            til<-til[this>=n & this<=m];
            it1<-length(that1)
            that1<-that1[this<=m & this>=n];
            neg2<-c(neg2+(it1-length(that1)));
            # (2.3) check for second peak occurrence within mass range #########
            if(length(til)>0){
                getthat<-c();
                within[]<-FALSE # store once and use over all j ################
                within[
                  peak2[,3]<=(((max(peak2[til,3])-peak2[i,3])*2)+rttol)+peak2[i,3] & 
                  peak2[,3]>=(((min(peak2[til,3])-peak2[i,3])*2)-rttol)+peak2[i,3]     
                ]<-TRUE;
                for(j in 1:length(til)){
                    if(
                      any ( # check for available masses within RT window ######
                            # use broad (not ppm) tolerance ####################
                          peak2[,1]<=((peak2[til[j],1]+(that1[j]))+2*delt) &
                          peak2[,1]>=((peak2[til[j],1]+(that1[j]))-2*delt) & 
                          within
                          )
                      ){
                          getthat<-c(getthat,TRUE);
                      }else{
                          getthat<-c(getthat,FALSE);
                      }   # if 1
                }       # for j
                neg1<-c(neg1+length(getthat[getthat==FALSE])); 
                til<-til[getthat];
                if(length(til)>0){
                  ret<-c(ret,peak2[til,3]-peak2[i,3]);
                  from<-c(from,rep(i,length(til)));
                  to<-c(to,til);
                  dif<-c(dif,that1[getthat]);
                } # if length til.3
            } # if length til.2
        } # if length til.1
      } 
      len<-length(dif);                               
      if(plotit==TRUE){
        text(length(make),0.03,labels="     done.",pos=2);
        text(2,-0.08,labels=paste(len," m/z increments included & ",neg1+neg2," m/z increments rejected.",sep=""),pos=4);
        make<-seq(1,len,20);
        plot.window(xlim=c(0,(length(make)+1)),ylim=c(-0.2,0.1))
        text(2,-0.15,labels="(2) Combine increments to homologue series...",pos=4);
        points(seq(1,length(make),1),rep(-0.2,length(make)),cex=1.2,col="darkgrey",pch=22,bg="darkgrey");
        points(seq(1,length(make),1),rep(-0.23,length(make)),cex=0.3,col="lightgrey",pch=22,bg="lightgrey");
      }
      ##########################################################################
      # (3) screen differences for homologues series of differences ############
      ##########################################################################
      ret<-ret[order(dif,decreasing=FALSE)];
      from<-from[order(dif,decreasing=FALSE)];
      to<-to[order(dif,decreasing=FALSE)];
      dif<-dif[order(dif,decreasing=FALSE)];
      done<-rep(FALSE,len);
      often1<-c(0); # how often gone into the first loop?
      often2<-c(0); # how often gone into the second loop?
      often3<-c(0); # how often gone into the third loop?
      peakID<-list(0);
      mzincr<-c();
      retincr<-c();
      retmin<-c();
      retmax<-c();
      last<-c(1);
      jup=1;
      jlow=1;
      for(i in 1:(length(dif))){
        if(plotit==TRUE){
          if(any(make==i)){
            points(seq(1,length(make),1)[make==i],-0.2,cex=1.2,col="green",pch=22,bg="green");
            last<-seq(1,length(make),1)[make==i];
          }
        }
        ########################################################################
        if(done[i]==FALSE){ # do not redo the same #############################
          # (3.1) within m/z window of increments ? (not ppm-precise) ##########
          while(  ((dif[jup]-dif[i]) <= 3*delt) & jup<len   ){jup=jup+1};
          if( jup!=len  ){jup=jup-1};
          while(  ((dif[i]-dif[jlow]) > 3*delt) & jlow<len  ){jlow=jlow+1};
          get1<-c(jlow:jup);
          get1<-get1[get1!=i];
          if((length(get1))>=(minlen)){
            # get1:   all feasible peaks, mz-tol = del
            # get12:  peaks related to dif[x], mz-tol = ppm
            # (3.2) are peaks related at all ? #################################
            if(any(from[get1]==to[i])){  # i.e. at least minlen >=3 
              # (3.3) first series step within RT-time window ? (seems weak) ###
              if(any(abs(ret[i]-ret[get1])<=rttol)){ 
                  often1<-often1+1;                           
                  # (3.4) within tolerances ? (ppm-precise) ####################
                  if(ppm==TRUE){
                    get12<-get1[
                      (abs(dif[get1]-dif[i])<=(mztol*2*peak2[from[get1],1]/1e6))  &
                      (abs(ret[i]-ret[get1])<=rttol)                              
                    ];
                  }else{
                    get12<-get1[
                      (abs(dif[get1]-dif[i])<=(mztol*2))                          &
                      (abs(ret[i]-ret[get1])<=rttol)                              
                    ];            
                  }
                  if(length(get12)>=(minlen)){
                  # (3.5) are peaks still related ? ############################
                  if(any(from[get12]==to[i])){
                        ########################################################
                        # build serie ##########################################
                        often2<-often2+1;
                        # (3.6) extent forward: get anything, as list ##########
                        # not yet referenced to peak IDs! ######################
                        listit1<-as.list(i);
                        listit1[[2]]<-get12[from[get12]==to[i]];
                        it<-c(2);
                        doit<-TRUE;
                        while(doit==TRUE){
                            get2<-c();
                            for(b in 1:length(listit1[[it]])){
                              # reset ppm-window to concerned peak: increasing mass ! 
                              if(ppm==TRUE){
                                get12<-get1[
                                  (abs(dif[get1]-dif[i])<=(mztol*2*peak2[from[listit1[[it]][b]],1]/1e6)) &
                                  (abs(ret[listit1[[it]][b]]-ret[get1])<=rttol)
                                  ]
                              }else{
                                get12<-get1[
                                  (abs(dif[get1]-dif[i])<=(mztol*2)) &
                                  (abs(ret[listit1[[it]][b]]-ret[get1])<=rttol)
                                  ]
                              }
                              get2<-c(get2,get12[from[get12]==to[listit1[[it]]][b]])
                            }
                            if(length(get2)>0){
                              it<-it+1;
                              listit1[[it]]<-get2;
                            }else{
                              doit<-FALSE;
                            } 
                        }
                        if(it>2){
                            # (3.7) extent backward: get anything as second list 
                            # from it=2 backward, must be part of get1 #########
                            # (3.7.1) 3->2 #####################################
                            get2<-c();
                            for(k in 1:length(listit1[[3]])){
                              if(ppm==TRUE){
                                get12<-get1[
                                  (abs(dif[get1]-dif[i])<=(mztol*2*peak2[from[listit1[[3]][k]],1]/1e6)) &
                                  (abs(ret[listit1[[3]][k]]-ret[get1])<=rttol)
                                  ]
                              }else{
                                get12<-get1[
                                  (abs(dif[get1]-dif[i])<=(mztol*2)) &
                                  (abs(ret[listit1[[3]][k]]-ret[get1])<=rttol)
                                  ]
                              }
                              get2<-c(get2,get12[to[get12]==from[listit1[[3]]][k]])                      
                            }
                            listit1[[2]]<-as.numeric(levels(as.factor(c(listit1[[2]],get2))));
                            # (3.7.2) 2->1 #####################################
                            get2<-c();
                            for(k in 1:length(listit1[[2]])){
                              if(ppm==TRUE){
                                get12<-get1[
                                  (abs(dif[get1]-dif[i])<=(mztol*2*peak2[from[listit1[[2]][k]],1]/1e6)) &
                                  (abs(ret[listit1[[2]][k]]-ret[get1])<=rttol)
                                  ]
                              }else{
                                get12<-get1[
                                  (abs(dif[get1]-dif[i])<=(mztol*2)) &
                                  (abs(ret[listit1[[2]][k]]-ret[get1])<=rttol)
                                  ]
                              }
                              get2<-c(get2,get12[to[get12]==from[listit1[[2]]][k]])                      
                            }
                            listit1[[1]]<-as.numeric(levels(as.factor(c(listit1[[1]],get2))));
                            # (3.7.3) 1->X #####################################
                            listit2<-as.list(listit1[[1]][1]);
                            listit2[[1]]<-listit1[[1]];
                            it2<-c(1);
                            doit<-TRUE; 
                            while(doit==TRUE){
                                get2<-c();
                                for(b in 1:length(listit2[[it2]])){
                                  # reset ppm-window to concerned peak: increasing mass ! 
                                  if(ppm==TRUE){
                                    get12<-get1[
                                      (abs(dif[get1]-dif[i])<=(mztol*2*peak2[from[listit2[[it2]][b]],1]/1e6)) &
                                      (abs(ret[listit2[[it2]][b]]-ret[get1])<=rttol)
                                      ]
                                  }else{
                                    get12<-get1[
                                      (abs(dif[get1]-dif[i])<=(mztol*2)) &
                                      (abs(ret[listit2[[it2]][b]]-ret[get1])<=rttol)
                                      ]
                                  }
                                  get2<-c(get2,get12[to[get12]==from[listit2[[it2]]][b]])
                                }
                                if(length(get2)>0){
                                  it2<-it2+1;
                                  listit2[[it2]]<-get2;
                                  listit2[[it2]]<-as.numeric(levels(as.factor(listit2[[it2]])));
                                }else{
                                  doit<-FALSE;
                                } 
                            }
                            # (3.8) merge (3.6) forward with (3.7) backward ####
                            if(length(listit2)>1){
                                listit<-as.list(1); # start list
                                it<-c(1);
                                for(b in 1:length(listit2)){
                                  listit[[it]]<-listit2[[(length(listit2)+1-it)]];
                                  it<-c(it+1);
                                }
                                for(b in 2:length(listit1)){
                                  listit[[it]]<-listit1[[b]];
                                  it<-c(it+1);
                                
                                }
                                it<-c(it-1);
                            }else{
                                listit<-listit1
                            }
                            ####################################################
                            if(length(listit)>=(minlen-1)){
                                # (3.9) insert last missing link, if available #
                                # warning: if available: NOT part of get1 ######
                                get3<-from[listit[[it]]]
                                get4<-c();
                                for(b in 1:length(get3)){
                                    if(ppm==TRUE){
                                        get4<-c(get4,
                                            seq(1,length(peak2[,1]),1)[
                                              peak2[,1]<=(peak2[get3[b],1]+(mztol*2*peak2[get3[b],1]/1e6)+dif[i])   &
                                              peak2[,1]>=(peak2[get3[b],1]-(mztol*2*peak2[get3[b],1]/1e6)+dif[i])   &
                                              peak2[,3]<=(peak2[get3[b],3]+rttol+(max(ret[listit[[it]]])))    &
                                              peak2[,3]>=(peak2[get3[b],3]-rttol+(max(ret[listit[[it]]]))) 
                                            ]                  
                                        )
                                    }else{
                                        get4<-c(get4,
                                            seq(1,length(peak2[,1]),1)[
                                              peak2[,1]<=(peak2[get3[b],1]+(mztol*2)+dif[i])   &
                                              peak2[,1]>=(peak2[get3[b],1]-(mztol*2)+dif[i])   &
                                              peak2[,3]<=(peak2[get3[b],3]+rttol+(max(ret[listit[[it]]])))    &
                                              peak2[,3]>=(peak2[get3[b],3]-rttol+(max(ret[listit[[it]]]))) 
                                            ]                  
                                        )
                                    }  
                                }
                                get4<-as.numeric(levels(as.factor(get4)));
                                if(length(get4)>0){
                                  for(b in 1:length(get4)){
                                    dif3<-c(peak2[get4[b],1]-peak2[from[listit[[it]]],1]);
                                    ret3<-c(peak2[get4[b],3]-peak2[from[listit[[it]]],3]);
                                  }
                                  it<-it+1;
                                }
                                # (3.10) is minimum series length reached ? ####
                                if((it)>=minlen){
                                    # stop("yes!")
                                    often3<-often3+1;   
                                    if(plotit==TRUE){                   
                                        points(last,-0.23,cex=0.3,col="red",pch=22,bg="red");
                                    }
                                    ############################################
                                    # (3.11) mark difs if included in a series # 
                                    dif2<-c();
                                    ret2<-c();
                                    for(k in 1:length(listit)){
                                      done[listit[[k]]]<-TRUE
                                      dif2<-c(dif2,dif[listit[[k]]]);
                                      ret2<-c(ret2,ret[listit[[k]]])
                                    };
                                    if(length(get4)>0){
                                      dif2<-c(dif2,dif3);
                                      ret2<-c(ret2,ret3);
                                    }
                                    ############################################
                                    # (3.12) retrieve mean(dif), mean (RT) & min.max(RT)
                                    mzincr<-c(mzincr,round(mean(dif2),digits=5));
                                    retincr<-c(retincr,round(mean(ret2),digits=3));
                                    retmin<-c(retmin,round(min(ret2),digits=3));
                                    retmax<-c(retmax,round(max(ret2),digits=3));
                                    ############################################
                                    # (3.13) store list of IDs of concerned peaks in list peakID
                                    listit3<-as.list(0);
                                    listit3[[1]]<-as.numeric(levels(as.factor(from[listit[[1]]])));
                                    for(k in 2:length(listit)){
                                      listit3[[length(listit3)+1]]<-as.numeric(levels(as.factor(from[listit[[k]]])))
                                    }
                                    if(length(get4)>0){                               
                                      listit3[[length(listit3)+1]]<-as.numeric(levels(as.factor(get4)))               
                                    }                                                 
                                    peakID[[length(peakID)+1]]<-listit3
                                } # if it >= minlen
                            } # length(listit)>minlen
                        } # it it > 2
                ################################################################   
                } # if two entries linked after set ppm
                } # if get 1 > minlen       
              } # if two entries linked
            } # if get1 > minlen
          } # if all > minlen
        } # if done
        ########################################################################
      } # while
      # plot(mzincr,retincr,cex=0.3)
      # remove first dummy entry in peakID list ################################
      peakID<-peakID[-1];
      ##########################################################################
      
      ##########################################################################
      # Generate outputs #######################################################
      ##########################################################################
      # (4) resort to dataset with increasing RT ###############################
      len<-length(peak2[,1]);
      convert<-seq(1,len,1);
      convert<-convert[order(peak2[,3],decreasing=FALSE)];
      getit<-seq(1,length(peak2[,3]),1)
      peak3<-peak2[order(peak2[,3],decreasing=FALSE),]
      # reassign peak IDs
      for(i in 1:length(mzincr)){
        for(j in 1:length(peakID[[i]])){
          for(k in 1:length(peakID[[i]][[j]])){
            peakID[[i]][[j]][k]<-getit[convert==peakID[[i]][[j]][k]];
          }
        }
      }
      ##########################################################################
      # (5) generate peaklist with links & #####################################
      # (6) generate component list with relevant m/z & RT increments ##########
      group1<-c();  # group ID
      group2<-c();  # peak IDs
      group3<-c();  # mzincr
      group4<-c();  # retincr
      group5<-c();  # retmin
      group6<-c();  # retmax
      getit1<-rep("0",len);      # (1) which group?
      getit2<-rep("0",len);      # (2) to which peak?
      getit3<-rep("none",len);   # (3) which mzincr?
      getit4<-rep("none",len);   # (4) which retincr?
      getit5<-rep("0",len);      # (5) level in series
      if(plotit==TRUE){
          text(length(make),-0.15,labels="     done.",pos=2);
          text(2,-0.26,labels=paste(length(peakID)," homologue series detected",sep=""),pos=4,col="red");
          screen(2);
          par(mar=c(4.5,4,0,0.3));
          plot(peak2[,3],peak2[,1],pch=19,col="lightgrey",cex=0.3,xlab="Retention time",ylab="m/z")#,xlim=c(15,20),ylim=c(450,500))
      }
      count<-c(1); # groupID
      getID<-rep(FALSE,length(peakID)); 
      for(i in 1:length(mzincr)){
         if(restr[1]==FALSE){
            # on (5) ###########################################################
            for(j in 1:(length(peakID[[i]])-1)){
                # generate links
                get1<-c("0");
                get2<-c("0");
                get3<-c("0");
                get4<-c("0");
                for(k in 1:length(peakID[[i]][[j+1]])){
                   get1<-paste(get1,"/",peakID[[i]][[j+1]][k],sep="");
                   get2<-paste(get2,"/",mzincr[i],sep="");
                   get3<-paste(get3,"/",retincr[i],sep="");
                   get4<-paste(get4,"/",j,sep="");
                }
                get1<-substr(get1,2,nchar(get1));
                get2<-substr(get2,2,nchar(get2));
                get3<-substr(get3,2,nchar(get3));
                get4<-substr(get4,2,nchar(get4));
                # fill in links
                for(k in 1:length(peakID[[i]][[j]])){
                     getit1[peakID[[i]][[j]][k]]<-paste(getit1[peakID[[i]][[j]][k]],"/",count,sep="");
                     getit2[peakID[[i]][[j]][k]]<-paste(getit2[peakID[[i]][[j]][k]],get1,sep="");
                     getit3[peakID[[i]][[j]][k]]<-paste(getit3[peakID[[i]][[j]][k]],get2,sep="");
                     getit4[peakID[[i]][[j]][k]]<-paste(getit4[peakID[[i]][[j]][k]],get3,sep="");
                     getit5[peakID[[i]][[j]][k]]<-paste(getit5[peakID[[i]][[j]][k]],get4,sep="");
                }
            }
            # last link without "to ID"
            for(k in 1:length(peakID[[i]][[j+1]])){
              getit1[peakID[[i]][[j+1]][k]]<-paste(getit1[peakID[[i]][[j+1]][k]],"/",count,sep="");
              getit3[peakID[[i]][[j+1]][k]]<-paste(getit3[peakID[[i]][[j+1]][k]],"/",mzincr[i],sep="");
              getit4[peakID[[i]][[j+1]][k]]<-paste(getit4[peakID[[i]][[j+1]][k]],"/",retincr[i],sep="");
              getit5[peakID[[i]][[j+1]][k]]<-paste(getit5[peakID[[i]][[j+1]][k]],"/",j+1,sep="");
            }      
            # on (6) ###########################################################
            getpeak<-c()
            for(j in 1:length(peakID[[i]])){
              for(k in 1:length(peakID[[i]][[j]])){
                  getpeak<-c(getpeak,peakID[[i]][[j]][k])
              }
            }
            if(plotit==TRUE){
                    points(peak3[getpeak,3],peak3[getpeak,1],pch=19,col="red",cex=0.45)
                    lines(peak3[getpeak,3],peak3[getpeak,1],col="red")
            }
            get1<-c(as.character(getpeak[1]))
            for(j in 2:length(getpeak)){
             get1<-paste(get1,",",getpeak[j],sep="")
            }
            group1<-c(group1,count);
            group2<-c(group2,get1);
            group3<-c(group3,mzincr[i]);
            group4<-c(group4,retincr[i]);
            group5<-c(group5,retmin[i]);
            group6<-c(group6,retmax[i]);
            ####################################################################
            count<-c(count+1);
         }else{
              if(ppm==TRUE){
                del<-(peakID[[i]][[1]][1]*mztol*2/1e6)
              }else{
               del<-mztol
              }
             if(any(
             (mzincr[i]-del)<=restr &
             (mzincr[i]+del)>=restr 
             )){
                  getID[i]<-TRUE;
                  # on (5) #####################################################
                  for(j in 1:(length(peakID[[i]])-1)){
                      # generate links
                      get1<-c("0");
                      get2<-c("0");
                      get3<-c("0");
                      get4<-c("0");
                      for(k in 1:length(peakID[[i]][[j+1]])){
                         get1<-paste(get1,"/",peakID[[i]][[j+1]][k],sep="");
                         get2<-paste(get2,"/",mzincr[i],sep="");
                         get3<-paste(get3,"/",retincr[i],sep="");
                         get4<-paste(get4,"/",j,sep="");
                      }
                      get1<-substr(get1,2,nchar(get1));
                      get2<-substr(get2,2,nchar(get2));
                      get3<-substr(get3,2,nchar(get3));
                      get4<-substr(get4,2,nchar(get4));
                      # fill in links
                      for(k in 1:length(peakID[[i]][[j]])){
                           getit1[peakID[[i]][[j]][k]]<-paste(getit1[peakID[[i]][[j]][k]],"/",count,sep="");
                           getit2[peakID[[i]][[j]][k]]<-paste(getit2[peakID[[i]][[j]][k]],get1,sep="");
                           getit3[peakID[[i]][[j]][k]]<-paste(getit3[peakID[[i]][[j]][k]],get2,sep="");
                           getit4[peakID[[i]][[j]][k]]<-paste(getit4[peakID[[i]][[j]][k]],get3,sep="");
                           getit5[peakID[[i]][[j]][k]]<-paste(getit5[peakID[[i]][[j]][k]],get4,sep="");
                      }
                  }
                  # last link without "to ID"
                  for(k in 1:length(peakID[[i]][[j+1]])){
                    getit1[peakID[[i]][[j+1]][k]]<-paste(getit1[peakID[[i]][[j+1]][k]],"/",count,sep="");
                    getit3[peakID[[i]][[j+1]][k]]<-paste(getit3[peakID[[i]][[j+1]][k]],"/",mzincr[i],sep="");
                    getit4[peakID[[i]][[j+1]][k]]<-paste(getit4[peakID[[i]][[j+1]][k]],"/",retincr[i],sep="");
                    getit5[peakID[[i]][[j+1]][k]]<-paste(getit5[peakID[[i]][[j+1]][k]],"/",j+1,sep="");
                  }
                  # on (6) #####################################################
                  getpeak<-c()
                  for(j in 1:length(peakID[[i]])){
                    for(k in 1:length(peakID[[i]][[j]])){
                        getpeak<-c(getpeak,peakID[[i]][[j]][k])
                    }
                  }
                  if(plotit==TRUE){
                          points(peak3[getpeak,3],peak3[getpeak,1],pch=19,col="red",cex=0.45)
                          lines(peak3[getpeak,3],peak3[getpeak,1],col="red")
                  }
                  get1<-c(as.character(getpeak[1]))
                  for(j in 2:length(getpeak)){
                   get1<-paste(get1,",",getpeak[j],sep="")
                  }
                  group1<-c(group1,count);
                  group2<-c(group2,get1);
                  group3<-c(group3,mzincr[i]);
                  group4<<-c(group4,retincr[i]);
                  group5<-c(group5,retmin[i]);
                  group6<-c(group6,retmax[i]);
                  count<-c(count+1);
            }
        }
      }
      if(plotit==TRUE){
        sc<-close.screen();if(sc[1]!=FALSE){for(m in 1:length(sc)){close.screen(sc[m])}};
        par<-def.par#- reset to default 
      }
      if(restr[1]!=FALSE){peakID<-peakID[getID];}
      grouping<-data.frame(group1,group2,group3,group4,group5,group6,group6-group5);
      names(grouping)<-c("group ID","peak IDs","m/z increment","RT increment","min. RT in series","max. RT in series","max.-min. RT");
      for(i in 1:length(getit1)){
        if(getit1[i]!="0"){
          getit1[i]<-substr(getit1[i],3,nchar(getit1[i]))
          getit3[i]<-substr(getit3[i],6,nchar(getit3[i]))
          getit4[i]<-substr(getit4[i],6,nchar(getit4[i]))
          getit5[i]<-substr(getit5[i],3,nchar(getit5[i]))
        }
        if(getit2[i]!="0"){
          getit2[i]<-substr(getit2[i],3,nchar(getit2[i]))
        }
      }
      grouped_samples<-data.frame(peak3,seq(1,len,1),getit1,getit5,getit2,getit3,getit4);
      names(grouped_samples)<-c("mz","intensity","RT","peak ID","group ID","series level","to ID","m/z increment","RT increment")
      # (4) store parameters used ##############################################
      parameters<-c(rttol,mztol,ppm,minmz,maxmz,minrt,maxrt,minlen);
      names(parameters)<-c("rttol","mztol","ppm","minmz","maxmz","minrt","maxrt","minlen");
      homol<-list(grouped_samples,parameters,grouping,restr,peakID);
      names(homol)<-c("Homologue Series","Parameters","Peaks in homologue series","m/z restrictions used","Peaks per level")
      ##########################################################################
      return(homol); ###########################################################
      ##########################################################################
}
