make.isos <-
function(iso_list,elements=c("C","N","Cl","S","Br","K"),charges=c(1:2),rm.isotopes=c("33S","36S","36S","40K","40K","15N"),rm.charges=c(2,1,2,1,2,2)){

    ############################################################################
    if(rm.isotopes[1]!=FALSE || rm.charges[1]!=FALSE){if(length(rm.isotopes)!=length(rm.charges)){stop("more rm.isotopes than rm.charges")}};
    if(length(rm.charges)>0 & rm.charges[1]!=FALSE){for(i in 1:length(rm.charges)){if(any(rm.charges[i]==charges)==FALSE){stop(paste("rm.charges ",rm.charges[i]," not in charges!",sep=""))};};};
    ############################################################################
    if(length(rm.charges)>0 & rm.charges[1]!=FALSE){rm.charges<-abs(rm.charges);}
    get_iso<-rep(FALSE,length(iso_list[,1]));
    for(i in 1:length(elements)){get_iso[iso_list[,1]==elements[i]]<-TRUE};
    isos<-iso_list[get_iso,];
    isos<-isos[isos[,4]!=0,];
    manyisos<-c(0);
    for(i in 1:length(elements)){
      that<-seq(1:length(isos[,1]))[isos[,1]==elements[i]];
      manyisos<-c(manyisos+length(that)-1);
      }
    isomat<-data.frame(matrix(ncol=6,nrow=manyisos,0));
    isomat[,5]<-1;
    colnames(isomat)<-c("name","dmass","dabund","how_often","#atoms","/C");
    counter<-c(1);
    for(i in 1:length(elements)){
      if(length(rep(1:length(isos[,1]))[isos[,1]==elements[i]])>1){
      isosub<-isos[,][isos[,1]==elements[i],];
      for(j in 2:(length(rep(1:length(isos[,1]))[isos[,1]==elements[i]]))){
            isomat[counter,1]<-c(paste(isosub[j,2]));
            isomat[counter,2]<-c(abs(isosub[1,3]-isosub[j,3]));
            isomat[counter,3]<-c((1/isosub[1,4])*isosub[j,4]);
            isomat[counter,6]<-c(isosub[j,5]);
            counter<-c(counter+1);
            }
      }else{stop(paste("not >1 isotope in list for ",elements[i]))};
      }
     isomat<-cbind(isomat,rep(abs(charges[1]),length(isomat[,1])));
     names(isomat)[7]<-c("z");
     isomat[,2]<-c(isomat[,2]/abs(charges[1]));
     if(length(charges)>1){
      isomat2<-isomat;
      for(i in 2:length(charges)){
          isomat3<-isomat2;
          isomat3[,2]<-c(isomat3[,2]/abs(charges[i]));
          isomat3[,7]<-abs(charges[i]);
          isomat<-rbind(isomat,isomat3);
      }
    }
    isomat<-isomat[order(isomat[,2]),];
    if(length(rm.isotopes)>0 & rm.isotopes[1]!=FALSE){
    for(i in 1:length(rm.isotopes)){
           isomat<-isomat[(isomat[,1]==rm.isotopes[i] & isomat[,7]==rm.charges[i])==FALSE,];
    }}
    manyisos<-length(isomat[,2]);
    iso<-list(isos,isomat,charges,manyisos,elements);
    names(iso)<-c("list of isotopes","list of isotope masses","charges","number of isotope m/z","elements");
    return(iso)
}
