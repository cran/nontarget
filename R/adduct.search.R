adduct.search <-
function(
                    peaklist,adduct_list,
                    rttol=0,mztol=2,massfrac=0.1,ppm=TRUE,
                    adducts=c("M+H","M+K","M+Na"),ion_mode="positive",
                    entry=20
                    ){

    ############################################################################
    # (0) check inputs #########################################################
    if(massfrac>1 || massfrac<=0){ stop("massfrac must be >0 and <=1") };
    if(ion_mode!="positive" & ion_mode!="negative"){stop("ion mode: positive or negative?")}
    for(i in 1:length(adducts)){if(any(adduct_list[,1]==adducts[i])!=TRUE &  any(adduct_list[adduct_list[,1]==adducts[i],6]==ion_mode)){stop(paste("Adduct ",adducts[i]," not in adduct_list!",sep=""))}};
    for(i in 1:length(adducts)){if((adduct_list[adduct_list[,1]==adducts[i],6]!=ion_mode)){stop(paste(adducts[i]," not in ion mode ",ion_mode,sep=""))}};
    if(length(peaklist)>3){stop("peaklist with > 3 columns not allowed")}
    ############################################################################
    cat("\n (1) Assemble lists...");
    # (1.1) sort peaklist etc
    getback<-order(peaklist[,3],peaklist[,1],decreasing=FALSE);
    samples<-peaklist[getback,];
    alls<-length(peaklist[,1]);
    # (1.2) storage ...
    ID<-seq(1:alls);
    getit1<-rep("none",alls);     # (1) which adduct removed <-> added?
    getit2<-rep("0",alls);        # (2) from which peak?
    getit4<-rep("0",alls);        # (3) to which peak?
    getit5<-rep("0",alls);        # (4) within [1] large or [2] small mass tolerance?
    # (1.3) retrieve selected subset of adduct_list
    adducts<-as.character(levels(as.factor(adducts)))
    these<-match(adducts,adduct_list[,1]);
    these<-these[is.na(these)==FALSE];
    add<-adduct_list[these,];
    add<-add[add[,6]==ion_mode,];
    if(length(add[,1])<1){stop("No selected adducts in list!")};
    add2<-data.frame(0,0,0,0,0);
    add3<-c();
    names(add2)<-c("mult1","mass1","mult2","mass2","count")
    that<-c(2);
    for(i in 1:(length(add[,1]))){
      for(j in 1:length(add[,1])){
          if(i!=j){
          add2<-rbind(add2,rep(0,5));
          add2[that,1]<-add[i,4];
          add2[that,2]<-add[i,5];
          add2[that,3]<-add[j,4];
          add2[that,4]<-add[j,5];
          that<-c(that+1);
          add3<-c(add3,paste(add[i,1],add[j,1],sep="<->"));
          }
       }; #j
    }; #i
    add2<-add2[-1,];
    #data.frame(add3,add2);
    cat("done.");
    ############################################################################

    ############################################################################
    # (2) run search ###########################################################
    cat("\n (2) Screen for mass increments...");
    if(ppm==TRUE){ppm2=1}else{ppm2=2};
    getit1a<-rep(0,alls*entry);
    getit2a<-rep(0,alls*entry);
    getit4a<-rep(0,alls*entry);
    getit5a<-rep(0,alls*entry);
    maxmass<-max(peaklist[,1]);
    #dyn.load("C:\\Program Files\\R\\R-2.13.1\\bin\\i386\\adductCpp.dll");
    #dyn.load(paste(.libPaths(),"/nontarget/temp/adductCpp.dll",sep=""));
    result<-.C("adduct",
      as.double(samples[,1]),as.double(samples[,3]),as.integer(length(samples[,1])),  # 3
      as.double(mztol*2),as.double(massfrac*2),as.double(rttol),    # 6
      as.integer(length(add2[,1])),                                                               # 7
      as.double(add2[,1]),as.double(add2[,2]),as.double(add2[,3]),as.double(add2[,4]),as.integer(add2[,5]),    # 12
      as.integer(entry),as.integer(ppm2),                                       # 14
      as.integer(getit1a),as.integer(getit2a),as.integer(getit4a),as.integer(getit5a) # 18
      ,PACKAGE="nontarget"
    );
    # (1) which adduct?
    for(i in 1:(alls-1)){for(j in 1:entry){if(result[15][[1]][(i-1)*entry+j]!=0){getit1[i]<-paste(getit1[i],result[15][[1]][(i-1)*entry+j],sep="/")}}};
    # (2) from which peak?
    for(i in 1:(alls-1)){for(j in 1:entry){if(result[16][[1]][(i-1)*entry+j]!=0){getit2[i]<-paste(getit2[i],result[16][[1]][(i-1)*entry+j],sep="/")}}};
    # (3) to which peak?
    for(i in 1:(alls-1)){for(j in 1:entry){if(result[17][[1]][(i-1)*entry+j]!=0){getit4[i]<-paste(getit4[i],result[17][[1]][(i-1)*entry+j],sep="/")}}};
    # (4) tolerance: small or large?
    for(i in 1:(alls-1)){for(j in 1:entry){
      if(result[18][[1]][(i-1)*entry+j]==1){getit5[i]<-paste(getit5[i],"small",sep="/")};
      if(result[18][[1]][(i-1)*entry+j]==2){getit5[i]<-paste(getit5[i],"large",sep="/")};
    }};
    if(result[13][[1]]!=entry){cat("WARNING: entry overflow -> links missing!")};
    #data.frame(ID,getit4,getit2,getit1,getit5)
    rm(result);
    #dyn.unload(paste(.libPaths(),"/nontarget/temp/adductCpp.dll",sep=""));
    ############################################################################
    
    ############################################################################
    # correct outputs for missing adduct combis (only submatrix searched!)
    for(i in 1:alls){
      if(getit1[i]!="none"){
        this12<-as.numeric(strsplit(as.character(getit1[i]),"/")[[1]][-1]);
        if(length(this12)==1){
          getit1[i]<-paste("none//",add3[this12],"//",sep="");
        }else{
          getit1[i]<-paste("none//",add3[this12[1]],"//",sep="");
          for(j in 2:length(this12)){getit1[i]<-paste(getit1[i],add3[this12[j]],"//",sep="");}
        };
      };
    };
    for(i in 1:alls){
      if(getit4[i]!="0"){
        this1<-strsplit(as.character(getit4[i]),"/")[[1]];
        this2a<-strsplit(as.character(getit1[i]),"//")[[1]];
        this5<-strsplit(as.character(getit5[i]),"/")[[1]];
        this2a<-this2a[this1!="0"];
        this5<-this5[this1!="0"];
        this1<-this1[this1!="0"];
        this2<-c();this3<-c();
        for(j in 1:length(this2a)){
              this2<-c(this2,strsplit(as.character(this2a[j]),"<->")[[1]][1]);
              this3<-c(this3,strsplit(as.character(this2a[j]),"<->")[[1]][2]);
              };
        for(j in 1:length(this1)){
          this10<-strsplit(as.character(getit4[as.numeric(this1[j])]),"/")[[1]]
          this11<-strsplit(as.character(getit1[as.numeric(this1[j])]),"//")[[1]]
          if((any(this10==as.character(i)) & any(this11==paste(this3[j],"<->",this2[j],sep=""))) == FALSE ){
               getit4[as.numeric(this1[j])]<-paste(getit4[as.numeric(this1[j])],i,sep="/");
               if(getit1[as.numeric(this1[j])]=="none"){
                getit1[as.numeric(this1[j])]="none//"
                getit1[as.numeric(this1[j])]<-paste(getit1[as.numeric(this1[j])],paste(this3[j],"<->",this2[j],"//",sep=""),sep="");
               }else{
                getit1[as.numeric(this1[j])]<-paste(getit1[as.numeric(this1[j])],paste(this3[j],"<->",this2[j],"//",sep=""),sep="");
               }
               getit5[as.numeric(this1[j])]<-paste(getit5[as.numeric(this1[j])],this5[j],sep="/");
          };
        };
      };
    };
    #data.frame(ID,getit4,getit1);
    cat("done.");
    ############################################################################

    ############################################################################
    # (4) group ################################################################
    cat("\n (3) Group peaks...");
    ############################################################################
    group1<-c(); # groupnumber?
    group2<-c(); # which peaks?
    group3<-c(); # which adducts?
    group4<-rep(0,alls); # groupnumber? 1-alls
    groupnumber<-c(1);
    getit1b<-getit1;
    getit4b<-getit4;
    for(i in 1:alls){
    if(getit4b[i]!="0"){
        this1<-as.numeric(strsplit(as.character(getit4b[i]),"/")[[1]][-1]);
        this2a<-strsplit(as.character(getit1b[i]),"//")[[1]][-1];
        this2<-c();this3<-c();
        for(j in 1:length(this2a)){
              this2<-c(this2,strsplit(as.character(this2a[j]),"<->")[[1]][1]);
              this3<-c(this3,strsplit(as.character(this2a[j]),"<->")[[1]][2]);
              };
        this4<-levels(as.factor(this2));
        for(j in 1:length(this4)){
          # assemble group information  ########################################
          group1<-c(group1,groupnumber);
          group2b<-(as.character(i)); ##########################################
          for(k in 1:length(this1[this2==this4[j]])){group2b<-paste(group2b,"/",this1[this2==this4[j]][k],sep="");};
          group2<-c(group2,group2b);
          group3b<-(as.character(this4[j])); ###################################
          for(k in 1:length(this2[this2==this4[j]])){group3b<-paste(group3b,"/",this3[this2==this4[j]][k],sep="");};
          group3<-c(group3,group3b);
          this5<-as.numeric(strsplit(as.character(group2[groupnumber]),"/")[[1]]); ###########
          for(k in 1:length(this5)){group4[this5[k]]<-paste(group4[this5[k]],"/",groupnumber,sep="")};
          ######################################################################
          # excise these entries in ALL group members ##########################
          for(k in 1:length(this1[this2==this4[j]])){
              g4<-strsplit(as.character(getit4b[this1[this2==this4[j]][k]]),"/")[[1]][-1];
              g1<-strsplit(as.character(getit1b[this1[this2==this4[j]][k]]),"//")[[1]][-1];
              getit<-rep(TRUE,length(g4));
              getit[(g4==as.character(i) & g1==paste(this3[this2==this4[j]][k],"<->",this4[j],sep=""))]<-FALSE
              for(m in 1:length(this1[this2==this4[j]])){
                getit[(g4==as.character(this1[this2==this4[j]][m]) & g1==paste(this3[this2==this4[j]][k],"<->",this3[this2==this4[j]][m],sep=""))]<-FALSE;
              }
              getit1b[this1[this2==this4[j]][k]]<-"none";
              getit4b[this1[this2==this4[j]][k]]<-"0";
              g1<-g1[getit];
              g4<-g4[getit];
              if(length(g1)>0){
                for(n in 1:length(g1)){
                    getit1b[this1[this2==this4[j]][k]]<-paste(getit1b[this1[this2==this4[j]][k]],"//",g1[n],sep="");
                    getit4b[this1[this2==this4[j]][k]]<-paste(getit4b[this1[this2==this4[j]][k]],"/",g4[n],sep="");
                    };
                getit1b[this1[this2==this4[j]][k]]<-paste(getit1b[this1[this2==this4[j]][k]],"//",sep="");
              };              
          };
          ######################################################################
          groupnumber<-c(groupnumber+1);
        } # for j
    } # if
    } # for i
    #data.frame(group1,group2,group3);
    cat("done.");
    ############################################################################

    ############################################################################
    # (4) apply rules ##########################################################
    #cat("\n (4) Check plausibility...");

    #cat("done.");
    ############################################################################

    ############################################################################
    # (4) generate output ######################################################
    cat("\n (4) Generate output...");
    ############################################################################
    parameters<-data.frame(rttol,mztol,massfrac,ppm,ion_mode);
    ############################################################################
    # correct entries:
    for(i in 1:alls){
      if(getit2[i]!="0"){getit2[i]<-substr(getit2[i],3,nchar(getit2[i]))};
      if(getit4[i]!="0"){getit4[i]<-substr(getit4[i],3,nchar(getit4[i]))};
      if(getit5[i]!="0"){getit5[i]<-sub("0/","",getit5[i])};
      if(getit1[i]!="none"){getit1[i]<-sub("none//","",getit1[i])};
      if(group4[i]!="0"){group4[i]<-substr(group4[i],3,nchar(group4[i]))};
    };
    for(i in 1:length(group2)){
      leng<-length(strsplit(as.character(group2[i]),"/")[[1]])
      for(k in 1:leng){group2[i]<-sub("/",",",group2[i],fixed=TRUE);};
    };
    ############################################################################
    # count hits !
    hits<-data.frame(adducts,rep(0,length(adducts)));
    names(hits)<-c("names","counts");
    for(i in 1:length(group1)){
      this1<-strsplit(as.character(group3[i]),"/")[[1]];
      for(j in 1:length(this1)){hits[hits[,1]==this1[j],2]<-hits[hits[,1]==this1[j],2]+1;};
    };
    ############################################################################
    # overlaps!
    overlaps<-data.frame(seq(1:100),rep(0,100));
    names(overlaps)<-c("number of groups in overlap","counts");
    for(i in 1:alls){if(group4[i]!="0"){overlaps[length(strsplit(as.character(group4[i]),"/")[[1]]),2]<-overlaps[length(strsplit(as.character(group4[i]),"/")[[1]]),2]+1;}};
    overlaps<-overlaps[overlaps[,2]!=0,];
    ############################################################################
    # groups!
    for(k in 1:length(group1)){
        group1[k]<-paste("/",group1[k],"/",sep="")
    }
    grouping<-data.frame(group1,group2,group3);
    names(grouping)<-c("group ID","peak IDs","adducts");
    ############################################################################
    adducts<-data.frame(samples[,1:3],ID,group4,getit4,getit1,getit5);
    names(adducts)<-c(names(samples)[1:3],"peak ID","group ID","to ID","adduct(s)","mass tolerance");
    adduct<-list(adducts,parameters,grouping,hits,overlaps);
    names(adduct)<-c("Adducts","Parameters","Peaks in adduct groups","Number of adducts","Number of peaks with grouped adducts overlapping");
    cat("done.\n\n");
    ############################################################################

    ############################################################################
    return(adduct);
    ############################################################################

}
