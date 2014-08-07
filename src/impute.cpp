#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <matvec/mim.h>
#include <matvec/parmMap.h>
#include <matvec/session.h>
#include <matvec/statdist.h>
#include <time.h>
#include <dirent.h>
#include <cmath>



using namespace std;



class qtlLocus{
public:
  double b;
  vector<int> delta;
  int active;
  void init(int nc,double sig2b,double pi);
};

class haploLocus{
public:
  vector<double> p;
  haploLocus(int nc){p.assign(nc,0);};
};

class hmmLoci{

public:
  vector<double> e,f,b,E,pState;
  hmmLoci() {hmmLoci(0,0);};
  hmmLoci(int ns,int nc){
    e.resize(ns);
    f.resize(nc);
    b.resize(nc);
    E.resize(ns);
    pState.resize(ns);
  };
   
};


class hmm{
public:
  int nStates,nComb;
  long nLoci;
  double c;
  
  vector<vector<double> > p,pComb; 
  vector<double> pi,piComb;
  vector<int> stateI,stateJ;
  vector<hmmLoci> loci;

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
	for(int i2=0;i2<nStates;i2++){
	  for(int j2=0;j2<nStates;j2++){
	    int c2=i2*(i2+1)/2+j2;
	    if(j2>i2) c2=j2*(j2+1)/2+i2;
	    pComb[c1][c2]+=p[i1][i2]*p[j1][j2];
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
    }
  };
  void initSeq(){
    for(long i=0;i<nLoci;i++){
      loci[i].f.assign(nComb,0.0);
      loci[i].b.assign(nComb,0.0);
    }
  };
  
  void initLoci(double alpha,double beta) {
    matvec::BetaDist betaDist(alpha,beta);
    hmmLoci locus(nStates,nComb);
    loci.resize(nLoci,locus);
    for(long i=0;i<nLoci;i++) {
      loci[i].e.assign(nStates,0.0);
      for(int l=0;l<nStates;l++){
	double val;
	val=betaDist.sample();
	loci[i].e[l]=val;
      }
    }
  };
  
  double emit(long i, int l, int x){
    int I=stateI[l];
    int J=stateJ[l];
    double val=1;
    if(x==0) val=(1.-loci[i].e[I])*(1.-loci[i].e[J]);
    if(x==1) val=loci[i].e[I]*(1.-loci[i].e[J])+(1.-loci[i].e[I])*loci[i].e[J];
    if(x==2) val=loci[i].e[I]*loci[i].e[J];
    return(val);
  } 
  
};
  
typedef unordered_map<string,int> idmap;

void forward(const long start, const long end,hmm &HMM, const vector<int> &X,const vector<double> picomb);
void backward(const long start, const long end,hmm &HMM, const vector<int> &X);
void calcPComb(long i,int nComb,hmmLoci *HMMlocus,vector<double> P);

int main(int argc,char **argv){

  idmap seqMap,phenMap;
  idmap::iterator seqMapIt,phenMapIt;

  ifstream Geno,Pheno;
  string filename;
  string line,id;
  vector<vector<int> > X;
  vector<vector<haploLocus> > haploProb;
  //vector<vector<double> > y;
  vector<string > ID;
  double val;
  int iVal;
  unsigned seed=3434241;
  seed=(unsigned int) time(NULL);
  DIR *dirPtr = opendir("matvec_trash");
        if (!dirPtr){
                umask(0);
                cout << "making matvec_trash directory" << endl;
                mkdir("./matvec_trash",0777);
        } 
  matvec::SESSION.initialize("matvec_trash");
  filename="Chrom1.dat";
  Geno.open(filename);
  //getline(Geno,line);
  int seq=0;
  while(getline(Geno,line)){
    stringstream linestr(line);
    linestr >> id;
    ID.push_back(id);
    seqMap[id]=seq++;
    vector<int> row;
    X.push_back(row);
    while(linestr >> iVal){
      iVal+=10;
      iVal/=10;
      X.back().push_back(iVal);
    }
    
  }
  
  cout <<  X.size() << endl;
  cout << X[0].size() << endl;
  cout << X[X.size()-1].size() << endl;
  
  for(int i=0;i<10;i++){
    cout << ID[i] ;
    for(int j=0;j<5;j++){
      cout << " " << X[i][j];
    }
    cout << endl;
  }

  //
  // Read Pheno
  //
  filename="BWTdrgEPDACC_ra.dat";
  Pheno.open(filename);
  getline(Pheno,line);
  vector<double> y,mu,rinverse;
  vector<int> phenSeq;
  while(getline(Pheno,line)){
    stringstream linestream(line);
    linestream >> id;
    seqMapIt=seqMap.find(id);
    if(seqMapIt != seqMap.end()){
      phenSeq.push_back(seqMapIt->second);
      linestream >> val;
      y.push_back(val);
      linestream >> val;
      mu.push_back(val);
      linestream >> val;
      rinverse.push_back(val);
      
    }
  }
  int nPheno=y.size();

  cout << "Matched Phenotypes "<< nPheno  << endl;
  for(int i=0;i< 10;i++){
    cout << phenSeq[i] << " " << ID[phenSeq[i]] << " " << y[i] << " " << mu[i] << " " << rinverse[i] << endl;
  }

  int nStates=6;
  long nLoci=X[0].size();
  hmmLoci locus;
  hmm HMM(nLoci,nStates,0.5,0.5);
  double priorCount=5;
  vector<double> doubleRow;
  

  double sum=0,pival=0.30,delta=-.05;
  for(int i=0;i<nStates;i++){
    HMM.pi[i]=pival;
    sum+=pival;
    pival+=delta;
    //HMM.pi[i]=1./((double) nStates);
  }
  for(int i=0;i<nStates;i++) HMM.pi[i]/=sum;

  double c=.95;
  HMM.c=c;
  HMM.resetP();
  HMM.initPComb();
  
  
  

  //Estimation
  int nComb=HMM.nComb;
  int nIter=10;
  for(int iter=0;iter<nIter;iter++){
    cout << "Iteration " << iter <<endl;
    HMM.initBW(priorCount);
    if(0)cout << HMM.piComb[0] << " " << HMM.emit(0,0,1) << " " << X[0][0]<< endl;
    
    for(int seq=0;seq< X.size();seq++){
      HMM.initSeq();

      forward(0, nLoci,HMM,X[seq],HMM.piComb);
      backward(0, nLoci,HMM,X[seq]);

      
      //Compute Probabilities and estimates
      for(long i=0;i<nLoci;i++){
	double Psum=0;
	double Pval;
	for(int l=0;l<nComb;l++){
	  Psum+=HMM.loci[i].f[l]*HMM.loci[i].b[l];
	}
	for(int l=0;l<nComb;l++){
	  int I=HMM.stateI[l];
	  int J=HMM.stateJ[l];
	  Pval=HMM.loci[i].f[l]*HMM.loci[i].b[l]/Psum;
	 
	  HMM.loci[i].pState[I]+=Pval;
	  HMM.loci[i].pState[J]+=Pval;
	  if(X[seq][i]==2){
	    HMM.loci[i].E[I]+=Pval;
	    HMM.loci[i].E[J]+=Pval;
	  }
	  if(X[seq][i]==1){
	    double p1,p2,wt1,wt2;
	    p1=HMM.loci[i].e[I];
	    p2=HMM.loci[i].e[J];
	    wt1=p1*(1.-p2);
	    wt2=p2*(1.-p1);
	    if((wt1+wt2)>0){
	      wt1=wt1/(wt1+wt2);
	      wt2=wt2/(wt1+wt2);
	      HMM.loci[i].E[I]+=Pval*wt1;
	      HMM.loci[i].E[J]+=Pval*wt2;
	    }
	  }
	}
      }
      
    }
    double nSeq;
    nSeq=(double) X.size();
    for(long i=0;i<nLoci;i++){
      for(int l=0;l<nStates;l++){
	if(HMM.loci[i].pState[l] > 0){
	  HMM.loci[i].e[l]=HMM.loci[i].E[l]/HMM.loci[i].pState[l];
	  if(isnan(HMM.loci[i].e[l])) cout << i << " " << l << " "<< HMM.loci[i].E[l]<< " " <<HMM.loci[i].pState[l]<< endl;
	}
	else{
	  HMM.loci[i].e[l]=0;
	}
      }
    }
    
    for(int l=0;l<nStates;l++){
      HMM.pi[l]=HMM.loci[0].pState[l]/(2.*(nSeq+priorCount));
    }
    HMM.resetP();
    HMM.initPComb();


    cout << std::fixed << std::setprecision(2);
    for(int l=0;l<nStates;l++){
      cout << "e  "<< l;
      for(long i=0;i<10;i++){
	cout << " " << HMM.loci[i].e[l];
      }
      cout << "    ";
      for(long i=nLoci-10;i<nLoci;i++){	
	cout << " " << HMM.loci[i].e[l];
      }
      cout << endl;
    }
    cout << endl;

cout << std::fixed << std::setprecision(2);
      for(int l=0;l<nStates;l++){
	cout << "P  "<< l;
	for(long i=0;i<10;i++){
	  cout << " " << HMM.loci[i].pState[l]/(2.*(nSeq+priorCount));
	}
	cout << "    ";
	for(long i=nLoci-10;i<nLoci;i++){
	  cout << " " << HMM.loci[i].pState[l]/(2.*(nSeq+priorCount));
	}
	cout << endl;
      }
      cout <<endl;



      cout << std::fixed << std::setprecision(2);
      vector<double> pSumV;
      pSumV.assign(nLoci,0.0);
      for(int l=0;l<nComb;l++){
	for(long i=0;i<10;i++) pSumV[i]+=HMM.loci[i].f[l]*HMM.loci[i].b[l];	
	for(long i=nLoci-10;i<nLoci;i++)pSumV[i]+=HMM.loci[i].f[l]*HMM.loci[i].b[l];
      }
      
      
      for(int l=0;l<nComb;l++){
	cout << "pc "<< HMM.stateI[l] << "/" << HMM.stateJ[l] << ":";
	for(long i=0;i<10;i++){
	  cout << " " << HMM.loci[i].f[l]*HMM.loci[i].b[l]/pSumV[i];
	}
	cout << "    ";
	for(long i=nLoci-10;i<nLoci;i++){
	  cout << " " << HMM.loci[i].f[l]*HMM.loci[i].b[l]/pSumV[i];
	}
	cout << endl;
      }
      cout <<endl;



      if(0){
      cout << std::fixed << std::setprecision(2);
      for(int l=0;l<nStates;l++){
	cout << "E "<< l;
	for(long i=0;i<10;i++){
	  cout << " " << HMM.loci[i].E[l];
	}
	cout << endl;
      }


      cout << std::fixed << std::setprecision(2);
      for(int l=0;l<nStates;l++){
	cout << "P "<< l;
	for(long i=0;i<10;i++){
	  cout << " " << HMM.loci[i].pState[l]/nSeq;
	}
	cout << endl;
      }
      
      cout << std::fixed << std::setprecision(2);
      for(int l=0;l<nComb;l++){
	cout << "f "<< l;
	for(long i=0;i<10;i++){
	  cout << " " << HMM.loci[i].f[l];
	}
	cout << endl;
      }
      cout << std::fixed << std::setprecision(2);
      for(int l=0;l<nComb;l++){
	cout << "b " << l;
	for(long i=0;i<10;i++){
	  cout << " " << HMM.loci[i].b[l];
	}
	cout << endl;
      }
    }
      
  }


  //
  // MCMC
  //

  matvec::UniformDist u;
  double sig2b=10,sig2e=15;
  double pi=.95;
  vector<qtlLocus> qtlVec;
  qtlVec.resize(nLoci);
  for(long i=0;i<nLoci;i++) qtlVec[i].init(nComb,sig2b,pi);
  
  
  
}

  

void qtlLocus::init(int nc,double sig2b,double pi){
  matvec::UniformDist u;
  matvec::NormalDist x(0,sig2b);
  b=x.sample();
  delta.assign(nc,0);
  for(int i=0;i<nc;i++) {
    if(u.sample() < 0.5) delta[i]=1;
  }
  active=0;
  if(u.sample()>pi) active=1;
}


void forward(const long start, const long end,hmm &HMM, const vector<int> &X,const vector<double> picomb){
  int nComb=HMM.nComb,nStates=HMM.nStates;
  double sum=0,val;
  for(int l=0;l<nComb;l++){
    val=picomb[l]*HMM.emit(start,l,X[start]);
      HMM.loci[start].f[l]=val;
      sum+=val;
   
  }
  for(int l=0;l<nComb;l++) HMM.loci[start].f[l]/=sum;//scale
  for(int i=start+1;i<end;i++){
    sum=0;
    for(int l=0;l<nComb;l++) {
      val=0;
      for(int k=0;k<nComb;k++){
	val+=HMM.loci[i-1].f[k]*HMM.pComb[k][l];
      }
      val*=HMM.emit(i,l,X[i]);
      HMM.loci[i].f[l]=val;
      sum+=val;
    }
    for(int l=0;l<nComb;l++) HMM.loci[i].f[l]/=sum;//scale
  }
}

void backward(const long start, const long end,hmm &HMM, const vector<int> &X){
  int nComb=HMM.nComb,nStates=HMM.nStates;
  double sum,val;
  for(int l=0;l<nComb;l++) HMM.loci[end-1].b[l]=1;
  for(long i=end-2;i>=start;i--){
    sum=0;
    for(int k=0;k<nComb;k++){
      val=0;
      for(int l=0;l<nComb;l++){
	val+=HMM.pComb[k][l]*HMM.emit(i+1,l,X[i+1])*HMM.loci[i+1].b[l]; 
      }
      HMM.loci[i].b[k]=val;
      sum+=val;
      
    }
    for(int k=0;k<nComb;k++) HMM.loci[i].b[k]/=sum;
  }
  
}


void calcPComb(long i,int nComb,hmmLoci *HMMlocus,vector<double> P){
  P.resize(nComb);
  double Pval,Psum=0;
  for(int l=0;l<nComb;l++){
    Pval=HMMlocus->f[l]*HMMlocus->b[l];
    Psum+=Pval;
    P[l]=Pval;
  }
  for(int l=0;l<nComb;l++) P[l]/=Psum;
}

