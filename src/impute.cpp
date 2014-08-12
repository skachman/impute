#include "impute.h"

int main(int argc,char **argv){

  matvec::UniformDist u;
  matvec::NormalDist Z;
  idmap seqMap,phenMap;
  idmap::iterator seqMapIt,phenMapIt;

  ifstream Geno,Pheno;
  string filename;
  string line,id;
  vector<vector<int> > X,XHaplo;
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
  vector<double> y,Xmu,rinverse;
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
      Xmu.push_back(val);
      linestream >> val;
      rinverse.push_back(val);
      
    }
  }
  int nPheno=y.size();

  cout << "Matched Phenotypes "<< nPheno  << endl;
  for(int i=0;i< 10;i++){
    cout << phenSeq[i] << " " << ID[phenSeq[i]] << " " << y[i] << " " << Xmu[i] << " " << rinverse[i] << endl;
  }

  int nStates=4;
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
  
  XHaplo.resize(nLoci);
  for(long i=0;i<nLoci;i++) XHaplo[i].resize(nPheno); 
  

  //Estimation
  int nComb=HMM.nComb;
  int nIter=10;
  vector<double> pVec(nComb);
  discrete_distribution<int> classDist;
  for(int iter=0;iter<nIter;iter++){
    cout << "Iteration " << iter <<endl;
    HMM.initBW(priorCount);
    if(0)cout << HMM.piComb[0] << " " << HMM.emit(0,0,1) << " " << X[0][0]<< endl;
    
    for(int seq=0;seq< X.size();seq++){
      HMM.initSeq();

      forward(0, nLoci,HMM,X[seq],HMM.piComb);
      backward(0, nLoci,HMM,X[seq]);

      vector<double> Pvec(nComb);
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
	  Pvec[l]=Pval;
	  HMM.loci[i].piVec[l]+=Pval;
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
      for(int l=0;l<nComb;l++) HMM.loci[i].piVec[l]/=(nSeq+priorCount);
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
  // Fill in HaploVec
  //
  for(int a=0;a<nPheno;a++){
    int seq=phenSeq[a];
    int ranClass;
    forward(0, nLoci,HMM,X[seq],HMM.piComb);
    //backward(0, nLoci,HMM,X[seq]);
    vector<double> Pvec(nComb);
    {
      double Psum=0;
      double Pval;
      for(int l=0;l<nComb;l++){
	Psum+=HMM.loci[nLoci-1].f[l];
      }
      for(int l=0;l<nComb;l++){
	Pval=HMM.loci[nLoci-1].f[l]/Psum;
	Pvec[l]=Pval;
      }
      double uSmp=u.sample();
      double Thresh=0;
      for(int l=0;uSmp>Thresh && l< nComb;l++){
	ranClass=l;
	Thresh+=Pvec[l];
      }
      XHaplo[nLoci-1][a]=ranClass;
    }
    for(long i=nLoci-2;i>=0;i--){
      double Psum=0;
      double Pval;
      for(int l=0;l<nComb;l++){
	Psum+=HMM.loci[i].f[l]*HMM.pComb[l][ranClass];
      }
      for(int l=0;l<nComb;l++){
	Pval=HMM.loci[i].f[l]*HMM.pComb[l][ranClass]/Psum;
	Pvec[l]=Pval;
      }
      double uSmp=u.sample();
      double Thresh=0;

      for(int l=0;uSmp>Thresh && l< nComb;l++){
	ranClass=l;
	Thresh+=Pvec[l];
	//if(i <5 && a < 5) cout << "Pv l "<< Pvec[l] << " "<< l << endl; 
      }
      //if(i <5 && a < 5) cout << i << " " << a << " " << uSmp << " " << Thresh  <<endl << Pvec[ranClass]<< endl << ranClass<<  endl; 
      XHaplo[i][a]=ranClass;  
    }
  }

  cout << "XHalpo " << XHaplo.size() <<"x" << XHaplo[0].size() << endl;

  for(int i=0;i<10;i++) 
    {
      cout << i ;
      for(int a=0;a<20; a++){ 
	cout << " " << setw(2) << XHaplo[i][a];
      }
      cout <<endl;
    }
  cout <<endl;


  //
  // MCMC
  //

  double sig2bPrior=.2,sig2ePrior=21;

  double pi=.95,mu=0;
  int nSamples=10;
  int nBurnIn=0;

  double nusig2e=10,nusig2b=4;
  double sig2e=sig2ePrior,sig2b=sig2bPrior;
  double scaleRes=sig2ePrior*(nusig2e-2.);
  double scaleB=sig2bPrior*(nusig2b-2.);
  double nuTilde,uSmp;
  double flipDelta=2./((double) nStates);
  vector<qtlLocus> qtlVec;
  vector<long> activeLoci;
  long activePos=0;
  uSmp=u.sample();
  qtlVec.resize(nLoci);
  activeLoci.resize(0);
  for(long i=0;i<nLoci;i++) {
    qtlVec[i].init(nStates,sig2b,pi);
    if(qtlVec[i].active){
      activeLoci.push_back(i);
    }
  }

  
  vector<double> yDev(nPheno),xVec(nPheno),probClass(nComb);
  

  for(int s=0;s<nSamples;s++){
    double ssb=scaleB*nusig2b,sse=scaleRes*nusig2e;
    cout <<"ssb sse "<< ssb <<  " " << sse << " " << scaleB << " " <<nusig2b<< endl;
    long nQTL=activeLoci.size();
    yDev=y;
    for(int a=0;a<nPheno;a++){
      int seq=phenSeq[a];
      yDev[a]-=mu;
      for(long i=0;i<nQTL;i++){
	long locus=activeLoci[i];
	int seqClass=XHaplo[locus][a];
	int I=HMM.stateI[seqClass];
	int J=HMM.stateJ[seqClass];
	yDev[a]-=(qtlVec[locus].delta[I]+qtlVec[locus].delta[J])*qtlVec[locus].b;
      }
    }

    activeLoci.resize(0);
    for(long i=0;i<nLoci;i++){
      int active=qtlVec[i].active;
      int nowActive=0;
      double logAoverI=log((1.-pi)/pi);
      double b=qtlVec[i].b;
      for(int a=0;a<nPheno;a++){
	int seq=phenSeq[a];
	int seqClass=XHaplo[i][a];
	int I=HMM.stateI[seqClass];
	int J=HMM.stateJ[seqClass];
  	xVec[a]=(qtlVec[i].delta[I]+qtlVec[i].delta[J]);
	double xb=xVec[a]*b;
	if(active)  yDev[a]+=xb;
	logAoverI+=(2.*yDev[a]-xb)*xb*rinverse[a]/(2.*sig2e);
	
      }
      uSmp=u.sample();
      if(log(uSmp/(1.-uSmp))<logAoverI) nowActive=1; 
      double bDelta=0;
      qtlVec[i].active=nowActive;
      if(nowActive) activeLoci.push_back(i); 
      double sumXY,sumXX;
      sumXY=0;
      sumXX=sig2e/sig2b;
      if(nowActive){
	for(int a=0;a<nPheno;a++){
	  double x=xVec[a];
	  sumXY+=x*yDev[a]*rinverse[a];
	  sumXX+=rinverse[a]*x*x;
	}
      }
      b=sumXY/sumXX+Z.sample()*sqrt(sig2e/sumXX); //sample b
     
      qtlVec[i].b=b;
      ssb+=b*b;
      //
      // Now sample delta
      //
      vector<double> logPdelta(nComb,0.0);
      //
      // cout << nComb << " " << logPdelta.size() << endl;
      if(nowActive) {	
	for(int a=0;a<nPheno;a++){
	  int seq=phenSeq[a];
	  int seqClass=XHaplo[i][a];
	  int I=HMM.stateI[seqClass];
	  int J=HMM.stateJ[seqClass];
	  logPdelta[I*(I+1)/2+J]-=b*b*rinverse[a]/(2.*sig2e);
	  logPdelta[I*(I+1)/2+I]+=(2.*yDev[a]*b-b*b)*rinverse[a]/(2.*sig2e);
	  logPdelta[J*(J+1)/2+J]+=(2.*yDev[a]*b-b*b)*rinverse[a]/(2.*sig2e);
	}
      }

      vector<int> deltaNew(nStates,0);
      vector<int> deltaOld=qtlVec[i].delta;
      for(int l=0;l<nStates;l++){
	if(u.sample()<flipDelta) deltaNew[l]=1-deltaOld[l];  
      }
      double logSwitch=0;
      for(int I=0;I<nStates;I++){
	for(int J=0;J<=I;J++){
	  int l=I*(I+1)/2+J;
	  logSwitch+=logPdelta[l]*(deltaNew[I]*deltaNew[J]-deltaOld[I]*deltaOld[J]);
	}
      }
      uSmp=u.sample();
      if(log(uSmp/(1.-uSmp))<logSwitch)qtlVec[i].delta=deltaNew;
      

      if(nowActive){
	for(int a=0;a<nPheno;a++){
	  int seq=phenSeq[a];
	  int seqClass=XHaplo[i][a];
	  int I=HMM.stateI[seqClass];
	  int J=HMM.stateJ[seqClass];
	  yDev[a]-=(qtlVec[i].delta[I]+qtlVec[i].delta[J])*b;
	}
      }
    }
    //update mu;
    double sumXY=0,sumXX=0;
    for(int a=0;a<nPheno;a++){
      sumXY+=rinverse[a]*(yDev[a]+mu);
      sumXX+=rinverse[a];
    }
    double newMu=sumXY/sumXX+Z.sample()*sqrt(sig2e/sumXX);
    for(int a=0;a<nPheno;a++) yDev[a]-=newMu-mu;
    mu=newMu;
    //update sig2e;
    nuTilde=((double) nPheno)+nusig2e;
    cout << "sse " << sse << "  -> ";
    for(int a=0;a<nPheno;a++) sse+=yDev[a]*yDev[a]*rinverse[a];
    cout << sse <<endl;
    double X2;
    X2=2.*matvec::sgamma(nuTilde/2.0); //Chi-squared;
    sig2e=sse/X2;

    //update sig2b;
    nuTilde=((double) nLoci)+nusig2b;
    X2=2.*matvec::sgamma(nuTilde/2.0); //Chi-squared;
    sig2b=ssb/X2;
    cout << endl;
    cout << "Sample: " << s << endl;
    cout << setprecision(8) << "Sig2b: " << sig2b << " Sig2e: " << sig2e << endl;
    cout << "mu " << mu << endl;
    cout << "number Active: " <<  activeLoci.size() << endl;
    for(int i=0; i< activeLoci.size()  ; i++ ) {
      cout <<setw(5) << " " << activeLoci[i];
      if((i %10)==9) cout <<endl;
    }
    cout << endl <<endl;
  }
	
      
}


  

void qtlLocus::init(int ns,double sig2b,double pi){
  matvec::UniformDist u;
  matvec::NormalDist x(0,sig2b);
  b=x.sample();
  delta.assign(ns,0);
  for(int i=0;i<ns;i++) {
    if(u.sample() < .5) delta[i]=1;
  }
  active=0;
  if(u.sample()>pi) active=1;
}


void forward(const long start, const long end,hmm &HMM, const vector<int> &X,const vector<double> &picomb){
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


void calcPComb(int nComb,hmmLoci &HMMlocus,vector<double> &P){
  P.resize(nComb);
  double Pval,Psum=0;
  for(int l=0;l<nComb;l++){
    Pval=HMMlocus.f[l]*HMMlocus.b[l];
    Psum+=Pval;
    P[l]=Pval;
  }
  for(int l=0;l<nComb;l++) P[l]/=Psum;
}






