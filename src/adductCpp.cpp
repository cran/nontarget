#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
using namespace std;
struct dat{
       vector<double> mass;
       vector<double> retent;
       bool operator() (int i,int j) { return (mass[i]<mass[j]);};
       };
extern "C" {   
void adduct(double *mass, double *retent, int *a, double *masstol, double *massfrac, double *rttol, int *many, double *add1, double *add2, double *add3, 
                   double *add4, int *add5,  int *entry, int *ppm2, int *getit1, int *getit2, int *getit4, int *getit5){
      int i = 0, j = 0, k, l, m, rtlow = 0, rtup = 0, howmany=0, upcount, lowcount, entry2=*entry;
      double uptol, lowtol, thismass, thismass2;
      vector<int> index;
      vector<int> index2;
      vector<int> index3;
      vector<int> getit1b; getit1b.assign ((*a+1),0);
      vector<int> getit2b; getit2b.assign ((*a+1),0);
      vector<int> getit4b; getit4b.assign ((*a+1),0);
      vector<int> getit5b; getit5b.assign ((*a+1),0);
      dat dat1;
      dat dat2;
      for(i=0;i<*a;i++){dat1.mass.push_back(mass[i]);};
      for(i=0;i<*a;i++){dat1.retent.push_back(retent[i]);};
      for(i=0;i<(*a);i++) {                                                                                                                   
       if((dat1.retent[i]!=dat1.retent[i-1]) || (i==0)){
            uptol=dat1.retent[i]+*rttol;
            lowtol=dat1.retent[i]-*rttol;
            while((dat1.retent[rtup]<=uptol) & (rtup<(*a-1))){rtup=(rtup+1);};
            rtup=(rtup-1);
            while((dat1.retent[rtlow]<lowtol) & (rtlow<(*a-1))){rtlow=(rtlow+1);};
            index.erase (index.begin(),index.end());
            for(j=rtlow; (j<=(rtup)) & (j<(*a)); j++){index.push_back(j);};
            sort(index.begin(),index.end(),dat1);                                                                   
           };
       howmany = index.size();
      if(howmany>0){       
           dat2.mass.erase (dat2.mass.begin(),dat2.mass.end());                               
           for(j=0;j<*many;j++){dat2.mass.push_back(((((dat1.mass[i]-add2[j])*(1.0/add1[j]))*add3[j])+add4[j]));};
           index3.erase (index3.begin(),index3.end());
           for(j=0;j<*many;j++){index3.push_back(j);};      
           sort(index3.begin(),index3.end(),dat2);
           if(*ppm2==1){
                           thismass=(dat1.mass[i]**masstol/1E6);
                           thismass2=(((dat1.mass[i]**masstol/1E6))**massfrac);
                           }else{
                           thismass=(*masstol);
                           thismass2=(*masstol**massfrac);
                           };
           upcount=0;
           lowcount=0;
          for(k = 0; k<*many; k++){
                 if(upcount<=(howmany-1)){     
                   while( (dat1.mass[index[upcount]]<=(dat2.mass[index3[k]]+thismass)) & (upcount<(howmany-1)) ){upcount=(upcount+1);};               
                   while( (dat1.mass[index[lowcount]]<(dat2.mass[index3[k]]-thismass)) & (lowcount<(howmany-1)) ){lowcount=(lowcount+1);};
                if(lowcount==(howmany-1)){
                   if( (dat1.mass[index[lowcount]]<=(dat2.mass[index3[k]] + thismass)) & (dat1.mass[index[lowcount]]>=(dat2.mass[index3[k]] - thismass)) ){
                          upcount=upcount+1;
                       }else{
                          upcount=lowcount;
                       };
               };
               if((upcount-lowcount)>0.0){
                                           index2.erase (index2.begin(),index2.end());
                                           for(l=lowcount;l<upcount;l++){index2.push_back(index[l]);};
                                           for(m=0;m<((int)index2.size());m++){                       
                                            if(getit2b[index2[m]]<(*entry+1)){ 
                                            getit2[(index2[m]**entry)+(getit2b[index2[m]])]=(i+1);                                                        
                                            getit2b[index2[m]] = getit2b[index2[m]]+1;                                                  
                                            };
                                           if(getit4b[i]<(*entry+1)){         
                                            getit4[i**entry+getit4b[i]]=(index2[m]+1);            
                                            getit4b[i] = getit4b[i]+1;
                                            };
                                            
                                            
                                            if(getit1b[i]<(*entry+1)){      
                                            getit1[i**entry+getit1b[i]]=(index3[k]+1);                                   
                                            getit1b[i] = getit1b[i]+1;  
                                            };                       
                                            
                                            
                                           if((dat1.mass[index2[m]]<=(dat2.mass[index3[k]] + thismass2)) & (dat1.mass[index2[m]]>=(dat2.mass[index3[k]] - thismass2))){                                                                                           
                                                 getit5[i**entry+getit5b[i]]=1;                       
                                                 getit5b[i] = getit5b[i]+1;                                                                                
                                            }else{                                               
                                                 getit5[i**entry+getit5b[i]]=2;                     
                                                 getit5b[i] = getit5b[i]+1;                                                                                                                                                                                                                                                                                   
                                            };                                                                                                                            
                                          add5[index3[k]]=add5[index3[k]]+1;                               
                                           } 
                                          if(upcount>howmany){upcount=lowcount;};
               };               
          };                
        };                                                                                                                                                                 
      };                                                        
      }; 
      for(j=0;j<(*a-1);j++){if(getit1b[j]>entry2){*entry=getit1b[j];};};
      for(j=0;j<(*a-1);j++){if(getit2b[j]>entry2){*entry=getit2b[j];};};
      for(j=0;j<(*a-1);j++){if(getit4b[j]>entry2){*entry=getit4b[j];};};
      for(j=0;j<(*a-1);j++){if(getit5b[j]>entry2){*entry=getit5b[j];};};
     }; 
};   
