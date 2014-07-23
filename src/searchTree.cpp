#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <vector>
#include <algorithm>
#include <deque>

#define RMATRIX(m,i,j) (REAL(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])
#define RMATRIX2(m,i,j) (INTEGER(m)[ INTEGER(GET_DIM(m))[0]*(j)+(i) ])
#define RRow(m) (INTEGER(GET_DIM(m))[0])
#define RCol(m) (INTEGER(GET_DIM(m))[1])


inline int inbound (SEXP data,int n,SEXP bounds){
    int a=1,i,ncol;
    ncol=RCol(data);
    for(i=0;i<ncol;i++){
        if(( RMATRIX(data,n,i)<RMATRIX(bounds,i,0) )||( RMATRIX(data,n,i)>RMATRIX(bounds,i,1) )){
            a=0;
            break;
        }
    }
    return a;
}

inline int boundbound (SEXP data,int n,SEXP bounds){
    int a=1,i,nrow;
    nrow=RRow(bounds);
    for(i=0;i<nrow;i++){
        if(( RMATRIX(data,n,(i*2))>RMATRIX(bounds,i,1) )||( RMATRIX(data,n,((i*2)+1))<RMATRIX(bounds,i,0) )){
            a=0;
            break;
        }
    }
    return a;
}

void search_tree_sub(SEXP data, SEXP tree, SEXP bounds, std::deque<int> &found){

    int n,nrow;
    double value;
    nrow=RRow(data);
    std::deque<int> at_nodes;

    /* find root node & initialize */
    for(n=0;n<nrow;n++){
        if(RMATRIX(tree,n,2)==1){
            break;
        }
    }
    at_nodes.push_back(n);
    while(at_nodes.size()>0){
        /* first node within search bounds? */
        if( inbound(data,at_nodes.front(),bounds)==1 ){
                found.push_back(at_nodes.front());
        }
        /* ... on LOSON */
        if(RMATRIX(tree,at_nodes.front(),0)!=0){
            /* still within bounds? */
            value=RMATRIX(data,at_nodes.front(),int((RMATRIX(tree,at_nodes.front(),3))-1));
            if(RMATRIX(bounds,int(RMATRIX(tree,at_nodes.front(),3)-1),0)<=value){
                at_nodes.push_back(int(RMATRIX(tree,at_nodes.front(),0)-1));
            }
        }
        /* ... on HISON */
        if(RMATRIX(tree,at_nodes.front(),1)!=0){
            /* still within bounds? */
            value=RMATRIX(data,at_nodes.front(),int(RMATRIX(tree,at_nodes.front(),3)-1));
            if(RMATRIX(bounds,int(RMATRIX(tree,at_nodes.front(),3)-1),1)>=value){
                /* include HISON & update its bounds */
                at_nodes.push_back(int(RMATRIX(tree,at_nodes.front(),1)-1));
            }
        }
        /* remove current = first node */
        at_nodes.pop_front();
    }

}

double *qua_a2;
int qua_b2;
bool smaller2 (int i,int j) {return(qua_a2[qua_b2+i]<qua_a2[qua_b2+j]);}

extern "C"{


/******************************************************************************/
/* build groups from peak-peak relations **************************************/
/******************************************************************************/

      SEXP metagroup(
        SEXP proffrom, /* must be sorted */
        SEXP profto
    ){

            PROTECT(proffrom = AS_INTEGER(proffrom));
            PROTECT(profto = AS_INTEGER(profto));
            int *proffrom2;
            proffrom2 = INTEGER_POINTER(proffrom);
            int *profto2;
            profto2 = INTEGER_POINTER(profto);
            int leng = LENGTH(proffrom);
            int n,m,g,k,atwhat,atto_n,atfrom_n,stay;

            SEXP group;
            PROTECT(group = NEW_INTEGER(leng));
            int *grouped;
            grouped = INTEGER_POINTER(group);
            for(n=0;n<leng;n++){*(grouped+n) = 0;}
            SEXP from;
            PROTECT(from = NEW_INTEGER(leng));
            int *atfrom;    /* stores indices */
            atfrom = INTEGER_POINTER(from);
            for(n=0;n<leng;n++){*(atfrom+n) = 0;}
            SEXP to;        /* stores profile IDs */
            PROTECT(to = NEW_INTEGER(leng));
            int *atto;
            atto = INTEGER_POINTER(to);
            for(n=0;n<leng;n++){*(atto+n) = 0;}

            g=1;
            for(n=0;n<leng;n++){
                if(*(grouped+n)==0){
                    *(grouped+n)=g;
                    /* initialize  */
                    *(atto)=*(proffrom2+n); /* cause it could be there repeatedly */
                    *(atto+1)=*(profto2+n);
                    atto_n=2;
                    stay=1;
                    while(stay>0){
                        /* to -> from */
                        /* given peak IDs in to, search indices in from having these IDs */
                        atfrom_n=0;
                        atwhat=n;
                        for(m=0;m<atto_n;m++){
                           if(*(atto+m)>=*(proffrom2+atwhat)){
                               for(k=atwhat;k<leng;k++){
                                    if( *(proffrom2+k) > *(atto+m) ){
                                        atwhat=k;
                                        break;
                                    }
                                    if( *(proffrom2+k) == *(atto+m) ){
                                        if(*(grouped+k)==0){
                                            *(grouped+k)=g;
                                            *(atfrom+atfrom_n)=k;
                                            atfrom_n++;
                                            atwhat=k;
                                        }
                                    }
                                }
                            }else{
                               for(k=atwhat;k>n;k--){
                                    if( *(proffrom2+k) < *(atto+m) ){
                                        atwhat=k;
                                        break;
                                    }
                                    if( *(proffrom2+k) == *(atto+m) ){
                                        if(*(grouped+k)==0){
                                            *(grouped+k)=g;
                                            *(atfrom+atfrom_n)=k;
                                            atfrom_n++;
                                            atwhat=k;
                                        }
                                    }
                                }
                            }
                        }
                        /* from -> to */
                        /* write with indices in from peak IDs into to */
                        atto_n=0;
                        if(atfrom_n>0){
                            for(k=0;k<atfrom_n;k++){
                                *(atto+atto_n)=*(profto2+*(atfrom+k));
                                atto_n++;
                            }
                        }else{
                            stay=0;
                        }
                    }
                    g++;
                }
            }

            UNPROTECT(5);
            return group;

}

/******************************************************************************/
/* peak -> peak search - ******************************************************/
/******************************************************************************/

    SEXP peak_search(
        SEXP peaklist,
        SEXP peakTree,
        SEXP mass_slots,
        SEXP int_slots,
        SEXP prec_mass,
        SEXP prec_ppm,
        SEXP prec_intens,
        SEXP RT_tol,
        SEXP pBar
    ){

            PROTECT(peaklist = AS_NUMERIC(peaklist));
            PROTECT(peakTree = AS_NUMERIC(peakTree));
            PROTECT(mass_slots = AS_NUMERIC(mass_slots));
            PROTECT(int_slots = AS_NUMERIC(int_slots));
            PROTECT(prec_mass = AS_NUMERIC(prec_mass));
            double prec_mass2 = NUMERIC_VALUE(prec_mass);
            PROTECT(prec_ppm = AS_NUMERIC(prec_ppm));
            int prec_ppm2 = NUMERIC_VALUE(prec_ppm);
            PROTECT(prec_intens = AS_NUMERIC(prec_intens));
            double prec_intens2 = NUMERIC_VALUE(prec_intens);
            PROTECT(RT_tol = AS_NUMERIC(RT_tol));
            double RT_tol2 = NUMERIC_VALUE(RT_tol);

            SEXP utilsPackage; /* definitions for the progres bar */
            PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
            SEXP percentComplete;
            PROTECT(percentComplete = NEW_NUMERIC(1));
            double *rPercentComplete;
            rPercentComplete = NUMERIC_POINTER(percentComplete);

            int n,m,s1,s2;
            int leng_peaks = RRow(peaklist);
            int leng_slots = RRow(mass_slots);
            double prec_dm;
            std::deque<int> from_peak (0);
            std::deque<int> to_peak (0);
            SEXP bounds_peaks;
            PROTECT(bounds_peaks = allocMatrix(REALSXP, 3, 2));

            for(n=0;n<leng_peaks;n++){
                *rPercentComplete = n;
                eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
                if(prec_ppm2==1){
                    prec_dm=(2*(RMATRIX(peaklist,n,0)*prec_mass2/1E6));
                }else{
                    prec_dm=(2*prec_mass2);
                }
                RMATRIX(bounds_peaks,2,0)=(RMATRIX(peaklist,n,2)-RT_tol2);
                RMATRIX(bounds_peaks,2,1)=(RMATRIX(peaklist,n,2)+RT_tol2);
                s1=to_peak.size();
                for(m=0;m<leng_slots;m++){
                    RMATRIX(bounds_peaks,1,0)=((RMATRIX(peaklist,n,0)-(prec_intens2*RMATRIX(peaklist,n,1)))/RMATRIX(int_slots,m,1));
                    RMATRIX(bounds_peaks,1,1)=((RMATRIX(peaklist,n,0)+(prec_intens2*RMATRIX(peaklist,n,1)))/RMATRIX(int_slots,m,0));
                    RMATRIX(bounds_peaks,0,0)=(RMATRIX(peaklist,n,0)+RMATRIX(mass_slots,m,0)-prec_dm);
                    RMATRIX(bounds_peaks,0,1)=(RMATRIX(peaklist,n,0)+RMATRIX(mass_slots,m,1)+prec_dm);
                    search_tree_sub(peaklist,peakTree,bounds_peaks,to_peak);
                } /* over mass_slots */
                s2=to_peak.size();
                for(m=s1;m<s2;m++){
                    from_peak.push_back(n);
                }
            } /* over peaks */

            s1=from_peak.size();
            SEXP results;
            PROTECT(results = allocMatrix(REALSXP, s1, 2));
            double *results2;
            results2 = REAL(results);
            for(n=0;n<s1;n++){
                results2[n]=(from_peak.front()+1);
                from_peak.pop_front();
                results2[s1+n]=(to_peak.front()+1);
                to_peak.pop_front();
            }
            UNPROTECT(12);
            return(results);
    }

/******************************************************************************/
/* Build random box tree ******************************************************/
/******************************************************************************/

    SEXP boxtree(
        SEXP data
    ){

            PROTECT(data = AS_NUMERIC(data));
            int n,nrow,ncol,n_disc,level,disc,doit,parent;
            ncol=RCol(data);
            n_disc=(ncol/2);
            nrow=RRow(data);
            SEXP results;
            PROTECT(results = allocMatrix(REALSXP, nrow, 5));
            double *results2;
            results2 = REAL(results);
            for(n=0;n<nrow;n++){
                results2[n]=0;          /* LOSON */
                results2[(1*nrow)+n]=0; /* MIDSON */
                results2[(2*nrow)+n]=0; /* HISON */
                results2[(3*nrow)+n]=0; /* level */
                results2[(4*nrow)+n]=0; /* disc */
            }

            /* first entry = root node */
            n=0;
            level=1;
            disc=0;
            results2[(3*nrow)+n]=level;
            results2[(4*nrow)+n]=(disc+1);

            for(n=1;n<nrow;n++){
                doit=2;
                parent=0;
                level=2;
                disc=0;
                while(doit==2){

                    /* LOSON *****************************************************************/
                    if( RMATRIX(data,n,((disc*2)+1))<=RMATRIX(data,parent,((disc*2)+1)) ){
                        if(results2[parent]==0){
                            results2[parent]=n+1;
                            results2[(3*nrow)+n]=level;
                            results2[(4*nrow)+parent]=(disc+1);
                            doit=1;
                        }else{
                            parent=(results2[parent]-1);
                            disc++;
                            if(disc>=n_disc){
                                disc=0;
                            }
                            level++;
                        }
                        continue;
                    }
                    /* HISON *****************************************************************/
                    if( RMATRIX(data,n,((disc*2)+0))>=RMATRIX(data,parent,((disc*2)+0)) ){
                        if(results2[(2*nrow)+parent]==0){
                            results2[(2*nrow)+parent]=n+1;
                            results2[(3*nrow)+n]=level;
                            results2[(4*nrow)+parent]=(disc+1);
                            doit=1;
                        }else{
                            parent=(results2[(2*nrow)+parent]-1);
                            disc++;
                            if(disc>=n_disc){
                                disc=0;
                            }
                            level++;
                        }
                        continue;
                    }
                    /* MIDSON *****************************************************************/
                    if(results2[(1*nrow)+parent]==0){
                        results2[(1*nrow)+parent]=n+1;
                        results2[(3*nrow)+n]=level;
                        results2[(4*nrow)+parent]=(disc+1);
                        doit=1;
                    }else{
                        parent=(results2[(1*nrow)+parent]-1);
                        disc++;
                        if(disc>=n_disc){
                            disc=0;
                        }
                        level++;
                    }
                }
            }

            UNPROTECT(2);
            return(results);

    }

/******************************************************************************/
/* Build random kd tree *******************************************************/
/******************************************************************************/

    SEXP kdtree4(
        SEXP data,
        SEXP pBar
    ){

            PROTECT(data = AS_NUMERIC(data));
            int n,m,ncol,nrow,doit,parent,level,disc;
            double mindist;
            ncol=RCol(data);
            nrow=RRow(data);
            SEXP results;
            PROTECT(results = allocMatrix(REALSXP, nrow, 6));
            double *results2;
            results2 = REAL(results);
            for(n=0;n<nrow;n++){
                    results2[n]=0;
                    results2[(1*nrow)+n]=0;
                    results2[(2*nrow)+n]=0;
                    results2[(3*nrow)+n]=0;
                    results2[(4*nrow)+n]=0;
                    results2[(5*nrow)+n]=0;
            }
            SEXP utilsPackage; /* definitions for the progres bar */
            PROTECT(utilsPackage = eval(lang2(install("getNamespace"), ScalarString(mkChar("utils"))), R_GlobalEnv));
            SEXP percentComplete;
            PROTECT(percentComplete = NEW_NUMERIC(1));
            double *rPercentComplete;
            rPercentComplete = NUMERIC_POINTER(percentComplete);

            /* first entry = root node */
            n=0;
            level=1;
            disc=0;
            results2[(2*nrow)+n]=level;
            results2[(3*nrow)+n]=(disc+1);

            for(n=1;n<nrow;n++){
                doit=2;
                parent=0;
                level=2;
                disc=0;
                while(doit==2){
                    if( RMATRIX(data,n,disc)>=RMATRIX(data,parent,disc) ){ /* HISON */
                        if(results2[(1*nrow)+parent]==0){
                            results2[(1*nrow)+parent]=n+1;
                            results2[(2*nrow)+n]=level;
                            results2[(3*nrow)+parent]=(disc+1);
                            mindist=0;
                            for(m=0;m<ncol;m++){
                                mindist=(mindist+fabs(RMATRIX(data,n,m)-RMATRIX(data,parent,m)));
                            }
                            results2[(4*nrow)+n]=mindist;
                            results2[(5*nrow)+n]=(parent+1);
                            if( (results2[(4*nrow)+parent]>mindist) | (results2[(5*nrow)+parent]==0) ){
                                results2[(4*nrow)+parent]=mindist;
                                results2[(5*nrow)+parent]=(n+1);
                            }
                            doit=1;
                        }else{
                            parent=(results2[(1*nrow)+parent]-1);
                            disc++;
                            if(disc>=ncol){
                                disc=0;
                            }
                            level++;
                        }
                    }else{ /* LOSON */
                        if(results2[parent]==0){
                            results2[parent]=n+1;
                            results2[(2*nrow)+n]=level;
                            results2[(3*nrow)+parent]=(disc+1);
                            mindist=0;
                            for(m=0;m<ncol;m++){
                                mindist=(mindist+fabs(RMATRIX(data,n,m)-RMATRIX(data,parent,m)));
                            }
                            results2[(4*nrow)+n]=mindist;
                            results2[(5*nrow)+n]=(parent+1);
                            if( (results2[(4*nrow)+parent]>mindist) | (results2[(5*nrow)+parent]==0) ){
                                results2[(4*nrow)+parent]=mindist;
                                results2[(5*nrow)+parent]=(n+1);
                            }
                            doit=1;
                        }else{
                            parent=(results2[parent]-1);
                            disc++;
                            if(disc>=ncol){
                                disc=0;
                            }
                            level++;
                        }
                    }

                }
                *rPercentComplete = n+1;
                eval(lang4(install("setTxtProgressBar"), pBar, percentComplete, R_NilValue), utilsPackage);
            }

            UNPROTECT(4);
            return(results);

    }


/******************************************************************************/
/* Build kd tree - partial sort on n_th element - only sorting indices ********/
/******************************************************************************/

    SEXP kdtree3(
        SEXP data
    ){

            PROTECT(data = AS_NUMERIC(data));
            qua_a2 = NUMERIC_POINTER(data);
            int n,m,ncol,nrow,splitit,last_k=0;
            ncol=RCol(data);
            nrow=RRow(data);

            SEXP results;
            PROTECT(results = allocMatrix(REALSXP, nrow, 4));
            double *results2;
            results2 = REAL(results);
            for(n=0;n<nrow;n++){
                    results2[n]=0;
                    results2[(1*nrow)+n]=0;
                    results2[(2*nrow)+n]=0;
                    results2[(3*nrow)+n]=0;
            }

            std::vector<int> from (1);
            from[0]=0;
            std::vector<int> to (1);
            to[0]=nrow;
            std::vector<int> level (1);
            level[0]=1;
            std::vector<int> disc (1);
            disc[0]=0;
            std::vector<int> parent (1);
            parent[0]=-1;

            int *here;
            here = new int[nrow];
            for(n=0;n<nrow;n++){here[n]=n;}

            while(from.size()>0){
                void R_CheckUserInterrupt(void);
                if((to[0]-from[0])>1){ /* at least 2 entries */
                    qua_b2=(disc[0]*nrow);
                    splitit=(from[0]+(((to[0]-1)-from[0])/2));
                    std::partial_sort(&here[from[0]],&here[splitit+1],&here[to[0]],smaller2);
                    while(splitit>from[0]){
                        if(RMATRIX(data,here[splitit],disc[0])==RMATRIX(data,here[splitit-1],disc[0])){
                            splitit--;
                        }else{
                            break;
                        }
                    };
                    /* set new parent */
                    results2[(2*nrow)+here[splitit]]=level[0];
                    results2[(3*nrow)+here[splitit]]=(disc[0]+1);
                    /* complete old parent */
                    if(parent[0]!=-1){
                        last_k=(disc[0]-1);
                        if(last_k<0){
                            last_k=(ncol-1);
                        }
                        if(RMATRIX(data,parent[0],last_k)<=RMATRIX(data,here[splitit],last_k)){ /* HISON & = */
                            results2[(1*nrow)+parent[0]]=(here[splitit]+1);
                        }else{ /* LOSON */
                            results2[parent[0]]=(here[splitit]+1);
                        }
                    }
                    /* new HISON */
                    from.push_back(splitit+1);
                    to.push_back(to[0]);
                    level.push_back(level[0]+1);
                    m=(disc[0]+1);
                    if(m>(ncol-1)){
                        disc.push_back(0);
                    }else{
                        disc.push_back(m);
                    }
                    parent.push_back(here[splitit]);
                    /* new LOSON */
                    if(splitit>from[0]){
                        from.push_back(from[0]);
                        to.push_back(splitit);
                        level.push_back(level[0]+1);
                        m=(disc[0]+1);
                        if(m>(ncol-1)){
                            disc.push_back(0);
                        }else{
                            disc.push_back(m);
                        }
                        parent.push_back(here[splitit]);
                    }
                    /* remove old subspace */
                    from.erase (from.begin());
                    to.erase (to.begin());
                    level.erase (level.begin());
                    disc.erase (disc.begin());
                    parent.erase (parent.begin());
                }else{ /* terminal node */
                    if(parent[0]!=-1){
                        last_k=(disc[0]-1);
                        if(last_k<0){
                            last_k=(ncol-1);
                        }
                        if( RMATRIX(data,parent[0],last_k) <= RMATRIX(data,here[from[0]],last_k) ){
                            results2[(1*nrow)+parent[0]]=(here[from[0]]+1);
                        }else{
                            results2[parent[0]]=(here[from[0]]+1);
                        }
                    }
                    results2[(2*nrow)+here[from[0]]]=level[0];
                    /* remove old subspace */
                    from.erase (from.begin());
                    to.erase (to.begin());
                    level.erase (level.begin());
                    disc.erase (disc.begin());
                    parent.erase (parent.begin());
                }
            }

            delete[] here;
            UNPROTECT(2);
            return(results);

    }


/******************************************************************************/
/* Build kd tree - partial sort on n_th element *******************************/
/******************************************************************************/

    SEXP kdtree2(
        SEXP data
    ){

            PROTECT(data = AS_NUMERIC(data));
            int n,m,ncol,nrow,splitit,last_k=0,counter,started,ended,gotnode=0,gotnode2;
            double median;
            ncol=RCol(data);
            nrow=RRow(data);

            SEXP results;
            PROTECT(results = allocMatrix(REALSXP, nrow, 4));
            double *results2;
            results2 = REAL(results);
            for(n=0;n<nrow;n++){
                    results2[n]=0;
                    results2[(1*nrow)+n]=0;
                    results2[(2*nrow)+n]=0;
                    results2[(3*nrow)+n]=0;
            }

            std::vector<int> from (1);
            from[0]=0;
            std::vector<int> to (1);
            to[0]=(nrow-1);
            std::vector<int> level (1);
            level[0]=1;
            std::vector<int> disc (1);
            disc[0]=0;
            std::vector<int> parent (1);
            parent[0]=-1;
            std::vector<int> swapit (1);
            swapit[0]=0;

            int *here;
            here = new int[nrow];
            for(n=0;n<nrow;n++){here[n]=n;}
            int *there;
            there = new int[nrow];
            double *store;
            store = new double[nrow];


            while(from.size()>0){
                void R_CheckUserInterrupt(void);
                if((to[0]-from[0])>0){ /* at least 2 entries */

                    /* get values */
                    counter=0;
                    if(swapit[0]==0){ /* here */
                        for(n=from[0];n<=to[0];n++){
                                store[counter]=RMATRIX(data,here[n],disc[0]);
                                counter++;
                        }
                    }else{ /* there */
                        for(n=from[0];n<=to[0];n++){
                                store[counter]=RMATRIX(data,there[n],disc[0]);
                                counter++;
                        }
                    }
                    splitit=((to[0]-from[0]+1)/2);
                    /*std::nth_element(store.begin(), store.begin()+(splitit-1), store.begin()+(counter));*/
                    std::nth_element(&store[0], &store[splitit-1], &store[counter]);
                    median=store[splitit-1];
                    started=from[0];
                    ended=to[0];
                    if(swapit[0]==0){ /* here -> there */
                        for(n=from[0];n<=to[0];n++){ /* find new parent node */
                            if(RMATRIX(data,here[n],disc[0])==median){
                                gotnode=n;
                                gotnode2=here[n];
                                break;
                            }
                        }
                        for(n=from[0];n<=to[0];n++){
                            if(n!=gotnode){
                                if(RMATRIX(data,here[n],disc[0])<median){
                                    there[started]=here[n];
                                    started++;
                                }else{
                                    there[ended]=here[n];
                                    ended--;
                                }
                            }
                        }
                    }else{ /* there -> here */
                        for(n=from[0];n<=to[0];n++){ /* find new parent node */
                            if(RMATRIX(data,there[n],disc[0])==median){
                                gotnode=n;
                                gotnode2=there[n];
                                break;
                            }
                        }
                        for(n=from[0];n<=to[0];n++){
                            if(n!=gotnode){
                                if(RMATRIX(data,there[n],disc[0])<median){
                                    here[started]=there[n];
                                    started++;
                                }else{
                                    here[ended]=there[n];
                                    ended--;
                                }
                            }
                        }
                    }
                    /* set new parent */
                    results2[(2*nrow)+gotnode2]=level[0];
                    results2[(3*nrow)+gotnode2]=(disc[0]+1);
                    /* complete old parent */
                    if(parent[0]!=-1){
                        last_k=(disc[0]-1);
                        if(last_k<0){
                            last_k=(ncol-1);
                        }
                        if(RMATRIX(data,parent[0],last_k)<=RMATRIX(data,gotnode2,last_k)){ /* HISON & = */
                            results2[(1*nrow)+parent[0]]=(gotnode2+1);
                        }else{ /* LOSON */
                            results2[parent[0]]=(gotnode2+1);
                        }
                    }
                    /* new HISON */
                    from.push_back(ended+1);
                    to.push_back(to[0]);
                    level.push_back(level[0]+1);
                    m=(disc[0]+1);
                    if(m>(ncol-1)){
                        disc.push_back(0);
                    }else{
                        disc.push_back(m);
                    }
                    parent.push_back(gotnode2);
                    if(swapit[0]==0){
                        swapit.push_back(1);
                    }else{
                        swapit.push_back(0);
                    }
                    /* new LOSON */
                    if(started>from[0]){
                        from.push_back(from[0]);
                        to.push_back(started-1);
                        level.push_back(level[0]+1);
                        m=(disc[0]+1);
                        if(m>(ncol-1)){
                            disc.push_back(0);
                        }else{
                            disc.push_back(m);
                        }
                        parent.push_back(gotnode2);
                        if(swapit[0]==0){
                            swapit.push_back(1);
                        }else{
                            swapit.push_back(0);
                        }
                    }
                    /* remove old subspace */
                    from.erase (from.begin());
                    to.erase (to.begin());
                    level.erase (level.begin());
                    disc.erase (disc.begin());
                    parent.erase (parent.begin());
                    swapit.erase (swapit.begin());

                }else{ /* terminal node */

                    if(parent[0]!=-1){
                        last_k=(disc[0]-1);
                        if(last_k<0){
                            last_k=(ncol-1);
                        }
                        if(swapit[0]==0){
                            if( RMATRIX(data,parent[0],last_k) <= RMATRIX(data,here[from[0]],last_k) ){
                                results2[(1*nrow)+parent[0]]=(here[from[0]]+1);
                            }else{
                                results2[parent[0]]=(here[from[0]]+1);
                            }
                        }else{
                           if( RMATRIX(data,parent[0],last_k) <= RMATRIX(data,there[from[0]],last_k) ){
                                results2[(1*nrow)+parent[0]]=(there[from[0]]+1);
                            }else{
                                results2[parent[0]]=(there[from[0]]+1);
                            }
                        }
                    }
                    if(swapit[0]==0){
                        results2[(2*nrow)+here[from[0]]]=level[0];
                        results2[(2*nrow)+here[to[0]]]=level[0];
                    }else{
                        results2[(2*nrow)+there[from[0]]]=level[0];
                        results2[(2*nrow)+there[to[0]]]=level[0];
                    }
                    /* remove old subspace */
                    from.erase (from.begin());
                    to.erase (to.begin());
                    level.erase (level.begin());
                    disc.erase (disc.begin());
                    parent.erase (parent.begin());
                    swapit.erase (swapit.begin());

             }
           }

            delete[] here;
            delete[] there;
            delete[] store;
            UNPROTECT(2);
            return(results);

    }


/******************************************************************************/
/* Build kd tree - full sort **************************************************/
/******************************************************************************/

    SEXP kdtree(
        SEXP data
    ){

            PROTECT(data = AS_NUMERIC(data));
            int n,m,ncol,nrow,splitit,last_k;
            ncol=RCol(data);
            nrow=RRow(data);

            std::vector<int> from (1);
            from[0]=0;
            std::vector<int> to (1);
            to[0]=(nrow-1);
            std::vector<int> level (1);
            level[0]=1;
            std::vector<int> disc (1);
            disc[0]=0;
            std::vector<int> parent (1);
            parent[0]=0;

            SEXP ordit;
            PROTECT(ordit = NEW_INTEGER(nrow));
            int *atordit;
            atordit = INTEGER_POINTER(ordit);
            for(n=0;n<nrow;n++){*(atordit+n) = n;};

            SEXP ordit_med;
            PROTECT(ordit_med = NEW_INTEGER(nrow));
            int *atordit_med;
            atordit_med = INTEGER_POINTER(ordit_med);

            SEXP ordit_med2;
            PROTECT(ordit_med2 = NEW_INTEGER(nrow));
            int *atordit_med2;
            atordit_med2 = INTEGER_POINTER(ordit_med2);

            SEXP intermed;
            PROTECT(intermed = NEW_NUMERIC(nrow));
            double *intermed2;
            intermed2 = NUMERIC_POINTER(intermed);
            for(n=0;n<nrow;n++){*(intermed2+n) = 0;};

            SEXP results;
            PROTECT(results = allocMatrix(REALSXP, nrow, 4));
            double *results2;
            results2 = REAL(results);
            for(n=0;n<nrow;n++){
                    results2[n]=0;
                    results2[(1*nrow)+n]=0;
                    results2[(2*nrow)+n]=0;
                    results2[(3*nrow)+n]=0;
            }

            while(from.size()>0){
                void R_CheckUserInterrupt(void);
                if((to[0]-from[0])>0){ /* at least 2 entries */
                    /* get split indice */
                    m=0;
                    for(n=from[0];n<=to[0];n++){
                        *(intermed2+m)=RMATRIX(data,*(atordit+n),disc[0]);
                        m++;
                    }
                    R_orderVector(atordit_med,m,Rf_lang1(intermed),FALSE,FALSE);
                    for(n=0;n<m;n++){
                        *(atordit_med2+n)=*(atordit+from[0]+*(atordit_med+n));
                    }
                    for(n=0;n<m;n++){
                        *(atordit+from[0]+n)=*(atordit_med2+n);
                    }
                    splitit=(from[0]+((to[0]-from[0])/2));
                    /* extented split indice to lower bound = assign equal values to HISON */
                    while(splitit>from[0]){
                        if(RMATRIX(data,*(atordit+splitit),disc[0])==RMATRIX(data,*(atordit+splitit-1),disc[0])){
                            splitit--;
                        }else{
                            break;
                        }
                    }
                    /* set new parent */
                    results2[(2*nrow)+*(atordit+splitit)]=level[0];
                    results2[(3*nrow)+*(atordit+splitit)]=(disc[0]+1);
                    /* complete old parent */
                    if(parent[0]!=0){
                        last_k=(disc[0]-1);
                        if(last_k<0){
                            last_k=(ncol-1);
                        }
                        if( RMATRIX(data,*(atordit+parent[0]),last_k) <= RMATRIX(data,*(atordit+splitit),last_k) ){ /* HISON & = */
                            results2[(1*nrow)+*(atordit+parent[0])]=(*(atordit+splitit)+1);
                        }else{ /* LOSON */
                            results2[*(atordit+parent[0])]=(*(atordit+splitit)+1);
                        }
                    }
                    /* new HISON */
                    from.push_back(splitit+1);
                    to.push_back(to[0]);
                    level.push_back(level[0]+1);
                    m=(disc[0]+1);
                    if(m>(ncol-1)){
                        disc.push_back(0);
                    }else{
                        disc.push_back(m);
                    }
                    parent.push_back(splitit);
                    /* new LOSON */
                    if(splitit!=from[0]){
                        from.push_back(from[0]);
                        to.push_back(splitit-1);
                        level.push_back(level[0]+1);
                        m=(disc[0]+1);
                        if(m>(ncol-1)){
                            disc.push_back(0);
                        }else{
                            disc.push_back(m);
                        }
                        parent.push_back(splitit);
                    }
                    /* remove old subspace */
                    from.erase (from.begin());
                    to.erase (to.begin());
                    level.erase (level.begin());
                    disc.erase (disc.begin());
                    parent.erase (parent.begin());
                }else{ /* terminal node */
                    if(parent[0]!=0){
                        last_k=(disc[0]-1);
                        if(last_k<0){
                            last_k=(ncol-1);
                        }
                        if( RMATRIX(data,*(atordit+parent[0]),last_k) <= RMATRIX(data,*(atordit+from[0]),last_k) ){ /* HISON & = */
                            results2[(1*nrow)+*(atordit+parent[0])]=(*(atordit+from[0])+1);
                        }else{ /* LOSON */
                            results2[*(atordit+parent[0])]=(*(atordit+from[0])+1);
                        }
                    }
                    results2[(2*nrow)+*(atordit+from[0])]=level[0];
                    results2[(2*nrow)+*(atordit+to[0])]=level[0];
                    /* remove old subspace */
                    from.erase (from.begin());
                    to.erase (to.begin());
                    level.erase (level.begin());
                    disc.erase (disc.begin());
                    parent.erase (parent.begin());
                }
            }

            UNPROTECT(6);
            return(results);

    }

/******************************************************************************/
/* Build random box tree ******************************************************/
/******************************************************************************/

    SEXP search_boxtree(
        SEXP data,
        SEXP tree,
        SEXP bounds,
        SEXP return_all
    ){

            PROTECT(data = AS_NUMERIC(data));
            PROTECT(tree = AS_NUMERIC(tree));
            PROTECT(bounds = AS_NUMERIC(bounds));
            double *bounds2;
            bounds2 = NUMERIC_POINTER(bounds);


            PROTECT(return_all = AS_INTEGER(return_all));
            int n,ret,disc,ncol,a,s;
            ncol=(RCol(data)/2);
            ret=INTEGER_VALUE(return_all);

            std::deque<int> doing;
            std::deque<int> found;

            doing.push_back(0);
            while(!doing.empty()){

                /* check old **************************************************************/
                a=1;
                for(n=0;n<ncol;n++){
                    if(
                       (bounds2[(n*2)+0]>RMATRIX(data,doing.front(),((n*2)+1))) ||
                       (bounds2[(n*2)+1]<RMATRIX(data,doing.front(),((n*2)+0)))
                    ){
                        a=0;
                        break;
                    }
                }
                if(a==1){
                    if(ret==1){
                        found.push_back(doing.front());
                    }else{
                        SEXP results;
                        PROTECT(results = allocMatrix(REALSXP, 1, 1));
                        double *results2;
                        results2 = REAL(results);
                        results2[0]=-2;
                        UNPROTECT(5);
                        return(results);
                    }
                }
                disc=int(RMATRIX(tree,doing.front(),4)-1);
                if(disc==-1){
                    doing.pop_front();
                    continue;
                }

                /* <UB ********************************************************************/
                if(bounds2[(disc*2)+1]<RMATRIX(data,doing.front(),((disc*2)+0))){
                    if(RMATRIX(tree,doing.front(),0)!=0){
                        doing.push_back(int(RMATRIX(tree,doing.front(),0)-1));
                    }
                    if(RMATRIX(tree,doing.front(),1)!=0){
                        doing.push_back(int(RMATRIX(tree,doing.front(),1)-1));
                    }
                    doing.pop_front();
                    s=doing.size();
                    continue;
                }
                /* >LB ********************************************************************/
                if(bounds2[(disc*2)+0]>RMATRIX(data,doing.front(),((disc*2)+1))){
                     if(RMATRIX(tree,doing.front(),1)!=0){
                        doing.push_back(int(RMATRIX(tree,doing.front(),1)-1));
                    }
                    if(RMATRIX(tree,doing.front(),2)!=0){
                        doing.push_back(int(RMATRIX(tree,doing.front(),2)-1));
                    }
                    doing.pop_front();
                    continue;
                }
                /* else ********************************************************************/
                if(RMATRIX(tree,doing.front(),0)!=0){
                    doing.push_back(int(RMATRIX(tree,doing.front(),0)-1));
                }
                if(RMATRIX(tree,doing.front(),1)!=0){
                    doing.push_back(int(RMATRIX(tree,doing.front(),1)-1));
                }
                if(RMATRIX(tree,doing.front(),2)!=0){
                    doing.push_back(int(RMATRIX(tree,doing.front(),2)-1));
                }
                doing.pop_front();
            }

            if(ret==1){ /* return all matches */
                SEXP results;
                s=found.size();
                PROTECT(results = allocMatrix(REALSXP, s, 1));
                double *results2;
                results2 = REAL(results);
                for(n=0;n<s;n++){
                        results2[n]=(found.front()+1);
                        found.pop_front();
                }
                UNPROTECT(5);
                return(results);
            }else{ /* any matches? -1 = none, -2 = at least one match */
                SEXP results;
                PROTECT(results = allocMatrix(REALSXP, 1, 1));
                double *results2;
                results2 = REAL(results);
                results2[0]=-1;
                UNPROTECT(5);
                return(results);
            }
    }


/******************************************************************************/
/* Search kd tree: nearest neighbour, euclidean *******************************/
/******************************************************************************/

    SEXP search_kdtree2(
        SEXP data,
        SEXP tree,
        SEXP ID,
        SEXP scaled
    ){

            PROTECT(data = AS_NUMERIC(data));
            PROTECT(tree = AS_NUMERIC(tree));
            PROTECT(ID = AS_INTEGER(ID));
            int i;
            i = (INTEGER_VALUE(ID)-1);
            PROTECT(scaled = AS_NUMERIC(scaled));
            double *scaled2;
            scaled2 = NUMERIC_POINTER(scaled);
            int n,nearest,ncol,nrow,disc,LOSON,HISON,son,parent,here,sized;
            double dist_global,dist_local,distance;
            dist_global=RMATRIX(tree,i,4);
            nearest=(int(RMATRIX(tree,i,5))-1);
            ncol=RCol(data);
            nrow=RRow(data);

            /* initialize */
            std::vector<int> SON;
            std::vector<int> PARENT;
            std::vector<double> DISTANCE;

            for(n=0;n<nrow;n++){ /* find starting point */
                if(RMATRIX(tree,n,2)==1){
                    break;
                }
            }
            SON.push_back(n);
            PARENT.push_back(-1);
            DISTANCE.push_back(0);

            /* run search */
            while(SON.size()>0){
                /* find subtree with smallest potential distance */
                dist_local= R_PosInf;
                here=0;
                sized=DISTANCE.size();
                for(n=0;n<sized;n++){
                    if(DISTANCE[n]<dist_local){
                        dist_local=DISTANCE[n];
                        here=n;
                    }
                }
                son=SON[here];
                parent=PARENT[here];
                distance=DISTANCE[here];
                SON.erase(SON.begin()+here);
                PARENT.erase(PARENT.begin()+here);
                DISTANCE.erase(DISTANCE.begin()+here);
                if(distance<=dist_global){
                    if(son!=i){
                        dist_local=0;
                        for(n=0;n<ncol;n++){
                            dist_local=(dist_local+(pow(((RMATRIX(data,son,n)-RMATRIX(data,i,n))/scaled2[n]),2)));
                        }
                        dist_local=sqrt(dist_local);
                        if(dist_local<dist_global){
                            dist_global=dist_local;
                            nearest=son;
                        }
                        if(RMATRIX(tree,son,4)>dist_local){
                            RMATRIX(tree,son,4)=dist_local;
                        };
                    }
                    disc=int(RMATRIX(tree,son,3)-1);
                    LOSON=int(RMATRIX(tree,son,0));
                    HISON=int(RMATRIX(tree,son,1));
                    if(LOSON!=0){
                        LOSON=(LOSON-1);
                        if(RMATRIX(data,i,disc)>RMATRIX(data,son,disc)){
                            dist_local=distance;
                            if(parent!=-1){
                                dist_local=(dist_local-pow(((RMATRIX(data,parent,disc)-RMATRIX(data,i,disc))/scaled2[disc]),2));
                            }
                            dist_local=(dist_local+pow(((RMATRIX(data,son,disc)-RMATRIX(data,i,disc))/scaled2[disc]),2));
                            dist_local=sqrt(dist_local);
                            if(dist_local<=dist_global){
                                SON.push_back(LOSON);
                                PARENT.push_back(son);
                                DISTANCE.push_back(dist_local);
                            }
                        }else{
                            SON.push_back(LOSON);
                            PARENT.push_back(son);
                            DISTANCE.push_back(distance);
                        }
                    }
                    if(HISON!=0){
                        HISON=(HISON-1);
                        if(RMATRIX(data,i,disc)<RMATRIX(data,son,disc)){
                            dist_local=distance;
                            if(parent!=-1){
                                dist_local=(dist_local-pow(((RMATRIX(data,parent,disc)-RMATRIX(data,i,disc))/scaled2[disc]),2));
                            }
                            dist_local=(dist_local+pow(((RMATRIX(data,son,disc)-RMATRIX(data,i,disc))/scaled2[disc]),2));
                            dist_local=sqrt(dist_local);
                            if(dist_local<=dist_global){
                                SON.push_back(HISON);
                                PARENT.push_back(son);
                                DISTANCE.push_back(dist_local);
                            }
                        }else{
                            SON.push_back(HISON);
                            PARENT.push_back(son);
                            DISTANCE.push_back(distance);
                        }
                    }
                }
            }

            SEXP results;
            PROTECT(results = allocMatrix(REALSXP, 1, (ncol+1)));
            double *results2;
            results2 = REAL(results);
            results2[0] = double(nearest+1);
            for(n=1;n<=ncol;n++){
                results2[n] = fabs(RMATRIX(data,nearest,n-1)-RMATRIX(data,i,n-1));
            }
            UNPROTECT(5);
            return(results);

    }

/******************************************************************************/
/* Search kd tree: range ******************************************************/
/******************************************************************************/

    SEXP search_kdtree(
        SEXP data,
        SEXP tree,
        SEXP bounds
    ){

            PROTECT(data = AS_NUMERIC(data));
            PROTECT(tree = AS_NUMERIC(tree));
            PROTECT(bounds = AS_NUMERIC(bounds));
            std::deque<int> found;
            int s,n;

            search_tree_sub(data, tree, bounds, found);

            s=found.size();
            SEXP results;
            PROTECT(results = allocMatrix(REALSXP, s, 1));
            double *results2;
            results2 = REAL(results);
            for(n=0;n<s;n++){
                    results2[n]=(found.front()+1);
                    found.pop_front();
            }
            UNPROTECT(4);
            return(results);

    }

}
