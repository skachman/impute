#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <list>
#include <random>
#include <unordered_map>
#include <algorithm>
//#include <matvec/mim.h>
//#include <matvec/parmMap.h>
//#include <matvec/session.h>
//#include <matvec/statdist.h>
#include <time.h>
#include <dirent.h>
#include <cmath>

#include <Eigen/Dense>

#include "Configuration.h"


enum variable_t { VAR_ID,VAR_DEP_VAR,VAR_RINVERSE,VAR_COVARIATE,VAR_CLASS,VAR_SKIP};

using namespace std;

using namespace Eigen;




class regionAnim{
 public:
  vector<double> gHat;
  vector<vector<double> > stateExp;
  regionAnim(int nr,int ns){vector <double> sp(ns,0.);gHat.resize(nr,0.),stateExp.resize(nr,sp);};
};

class regionQTL{
 public:
  vector<int> qtlLoci;
  int midpoint;
};

class locusMap{
 public:
  string name;
  int chrom,isSNP,isQTL;
  long pos,start,stop;
  locusMap(string s,int c, long p,int SNP,int QTL){name=s;chrom=c;pos=p;isSNP=SNP;isQTL=QTL;};
  locusMap(string s,int c, long p){name=s;chrom=c;pos=p;isSNP=1;isQTL=-1;};
  locusMap(){};
  
};

class haploMapLocus{
 public:
  vector<int> start;
  vector<int> a;
  haploMapLocus(int ns=0,int np=0){start.resize(ns+1,0);a.resize(2*np);};
};

bool locusMapCompare(locusMap &A,locusMap &B);

class qtlLocus{
 public:
  double b;
  vector<int> delta;
  int active;
  list<long> *activePos;
  void init(int nc,double sig2b,double pi);
  void init(int nc);
  void updateSum(const qtlLocus &A);
};

class qtlResultLocus{
 public:
  double b;
  vector<double> delta;
  double modelFreq;
  void init(int ns){delta.resize(ns);};

};

class qtl2Locus{
 public:
  Vector2d b;
  vector<int> delta;
  int active;
  vector<int> activeTrait;
  list<long> *activePos;
  void init(int nc,Matrix2d bL,double pi,double rho);
  void init(int nc);
  void updateSum(const qtl2Locus &A);
};

class haploLocus{
public:
  vector<double> p;
  haploLocus(int nc){p.assign(nc,0);};
};

class hmmLoci{
public:
  vector<double> e,f,b,E,pState,piVec;
  int newChrom,pos;
  hmmLoci() {hmmLoci(0,0);};
  hmmLoci(int ns,int nc){
    newChrom=0;
    e.resize(ns);
    f.resize(nc);
    b.resize(nc);
    E.resize(ns);
    pState.resize(ns);
    piVec.resize(nc);
  };
   
};


class hmm{
public:
  int nStates,nComb;
  long nLoci;
  double c,lambda;
  
  vector<vector<double> > p,pComb; 
  vector<double> pi,piComb;
  vector<int> stateI,stateJ;
  vector<hmmLoci> loci;
  vector<locusMap> *lociMapPt;

  void resize(const long nL){
    hmmLoci locus(nStates,nComb);
    nLoci=nL;
    loci.resize(nLoci,locus);
  };
  void initPComb(){
    pComb.resize(nComb);
    piComb.assign(nComb,0.0);
    for(int c=0;c<nComb;c++){
      pComb[c].assign(nComb,0.0);
    }

    for(int i1=0;i1<nStates;i1++){
      for(int j1=0;j1<nStates;j1++){
	int c1=i1*(i1+1)/2+j1;
	if(j1>i1) c1=j1*(j1+1)/2+i1;
	piComb[c1]+=pi[i1]*pi[j1];
	if(j1<=i1){
	  for(int i2=0;i2<nStates;i2++){
	    for(int j2=0;j2<nStates;j2++){
	      int c2=i2*(i2+1)/2+j2;
	      if(j2>i2) c2=j2*(j2+1)/2+i2;
	      pComb[c1][c2]+=p[i1][i2]*p[j1][j2];
	    } 
	  }
	}
      }
    }
  };

  void resetP(){
    for(int i=0;i<nStates;i++){
      for(int j=0;j<nStates;j++){
	p[i][j]=((1-c)*pi[j]/(1-pi[i]));
      }
      p[i][i]=c; 
    }
  };
  hmm(const long nL,const int nS,const double alpha,const double beta){
    nLoci=nL;
    nStates=nS;
    initComb();
    initLoci(alpha,beta);
    pi.resize(nStates);
    piComb.resize(nComb);
    p.resize(nStates);
    for(int k=0;k<nStates;k++) p[k].resize(nStates,0.0);
  };

  void initComb(){
    nComb=nStates*(nStates+1)/2;
    stateI.resize(nComb);
    stateJ.resize(nComb);
    for(int i=0;i<nStates;i++) {
      for(int j=0;j<=i;j++){
	int c=i*(i+1)/2+j;
	stateI[c]=i;
	stateJ[c]=j;
      }
    }
  };
  void initBW(const double priorCount){
    for(long i=0;i<nLoci;i++){
      loci[i].E.assign(nStates,priorCount);
      loci[i].pState.assign(nStates,2.*priorCount);
      loci[i].piVec.assign(nComb,priorCount);
    }
  };
  void initSeq(){
    for(long i=0;i<nLoci;i++){
      loci[i].f.assign(nComb,0.0);
      loci[i].b.assign(nComb,0.0);
    }
  };
  
  void initLoci(double alpha,double beta) {
    random_device rd;
    mt19937 gen(rd());
    gamma_distribution<double> gamA(alpha,1.),gamB(beta,1.);
    
    hmmLoci locus(nStates,nComb);
    loci.resize(nLoci,locus);
    for(long i=0;i<nLoci;i++) {
      loci[i].e.assign(nStates,0.0);
      for(int l=0;l<nStates;l++){
	double val;
	double g1=gamA(gen);
	double g2=gamB(gen);
	val=g1/(g1+g2);
	loci[i].e[l]=val;
      }
    }
  };
  
  double emit(const long i, const int l, const int x){
    int I=stateI[l];
    int J=stateJ[l];
    double val=1;
    if(x==0) val=(1.-loci[i].e[I])*(1.-loci[i].e[J]);
    else if(x==1) val=loci[i].e[I]*(1.-loci[i].e[J])+(1.-loci[i].e[I])*loci[i].e[J];
    else if(x==2) val=loci[i].e[I]*loci[i].e[J];
    return(val);
  } 
  int write(const string &filename);
  int read(const string &filename);
};
  
typedef unordered_map<string,int> idmap;

void  forward(const long start, const long end,hmm &HMM, const vector<int> &X,const vector<double> &picomb);
void backward(const long start, const long end,hmm &HMM, const vector<int> &X,const vector<double> &picomb);
void  forwardVec(const long start, const long end,hmm &HMM, const vector<int> &X,const vector<double> &picomb,vector< vector<double> > &fVec);
void backwardVec(const long start, const long end,hmm &HMM, const vector<int> &X,const vector<double> &picomb,vector< vector<double> > &bVec);
void calcPComb(int nComb,hmmLoci &HMMlocus,vector<double> &P);
void calcDeltaProposal(const double sig2e,const double sig2b,const vector<int> &delta,
			 vector<int> &dVec,vector<double> &rhsV,vector<double> &lhsV,vector<double> &lhsVs,const double pi,const int nDeltaStates,const int nStates,
			 double &pInactive,double &maxAoverI,vector<double> &AoverIVec,double &psum);

  int writeX(const string &filename);
  int readX(const string &filename);
