#include "impute.h"
#include <limits.h>

int main(int argc,char **argv){

  

  //cout << "Maximum value for int:  " << numeric_limits<int>::max() << '\n';
  //cout << "Maximum value for long: " << numeric_limits<long>::max() << '\n';

  matvec::UniformDist u;
  matvec::NormalDist Z;
  idmap seqMap,phenMap;
  idmap::iterator seqMapIt,phenMapIt;

  ifstream Geno,Pheno,MapFile;
  ofstream MCMCSamples,QTLResults,gHatResults;
  string filename;
  string line,id;
  vector<vector<int> > X,XHaplo;
  vector<vector<haploLocus> > haploProb;
  vector<locusMap> lociMap;
  locusMap aMapLocus;

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
  filename="Map.txt";
  MapFile.open(filename);
  getline(MapFile,line);
  while(getline(MapFile,line)){
    stringstream linestr(line);
    linestr >> aMapLocus.name;
    linestr >> aMapLocus.chrom;
    linestr >> aMapLocus.pos;
    lociMap.push_back(aMapLocus);
  }
  cout << "read in Map " << endl;
  sort(lociMap.begin(),lociMap.end(),locusMapCompare);
  map<string,long> mapOrder;
  for(long i=0;i<lociMap.size();i++) {
    mapOrder[lociMap[i].name]=i;
  }

  vector<int> chromStart;
  int lastChrom=lociMap[0].chrom;
  chromStart.push_back(0);
  for(long i=1;i<lociMap.size();i++){
    if(lastChrom!=lociMap[i].chrom){
      chromStart.push_back(i);
      lastChrom=lociMap[i].chrom;
    }
  }
  chromStart.push_back(lociMap.size());
  int nChrom=chromStart.size()-1;
  for(int chr=0;chr<nChrom;chr++) cout << chromStart[chr] << " ";
  cout << endl;

  vector<long> posVector(lociMap.size(),-1);

  filename="geno.dat";
  Geno.open(filename);
  getline(Geno,line);
  stringstream linestr(line);
  linestr >> id;
  long ipos=0;
  while(linestr >> id){
    posVector[ipos++]=mapOrder[id];
  }

  int seq=0;
  while(getline(Geno,line)){
    stringstream linestr(line);
    linestr >> id;
    ID.push_back(id);
    seqMap[id]=seq++;
    vector<int> row(mapOrder.size());
    X.push_back(row);
    long i=0;
    if((seq%200)==0) cout << seq << endl;
    while(linestr >> iVal){
      iVal+=10;
      iVal/=10;
      X.back()[posVector[i++]]=iVal;
    }
    
  }
  cout << seq <<endl;
  cout <<  X.size() << endl;
  cout << X[0].size() << endl;
  cout << X[X.size()-1].size() << endl;
  double nSeq;
  nSeq=(double) X.size();
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
  filename="BWT.dat";
  Pheno.open(filename);
  getline(Pheno,line);
  vector<double> y,Xmu,rinverse;
  vector<int> phenSeq;
  vector<string> phenID;
  int notFound=0;
  while(getline(Pheno,line)){
    stringstream linestream(line);
    linestream >> id;
    seqMapIt=seqMap.find(id);
    if(seqMapIt != seqMap.end()){
      phenID.push_back(id);
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
    }
    for(int i=0;i<nStates;i++) HMM.pi[i]/=sum;

    for(int chr=0;chr<nChrom;chr++) HMM.loci[chromStart[chr]].newChrom=1;

  double c=.95;
 
    HMM.c=c;
    HMM.resetP();

    string hmmFileName;
    int failCode;
    hmmFileName="HMM" + to_string(HMM.nStates) + "_" + to_string(HMM.nLoci) +".bin";
    if((failCode=HMM.read(hmmFileName))!=0) cout << hmmFileName << " failed to be read  with code " << failCode << endl;

    HMM.initPComb();
  
  XHaplo.resize(nLoci);
  for(long i=0;i<nLoci;i++) XHaplo[i].resize(nPheno); 
  

  //Estimation
  int nComb=HMM.nComb;
  int nIter=0;
  vector<double> pVec(nComb);

  for(int iter=0;iter<nIter;iter++){
    cout << "Iteration " << iter <<endl;
    HMM.initBW(priorCount);
    
    for(int seq=0;seq< X.size();seq++){
      HMM.initSeq();

      forward(0,nLoci,HMM,X[seq],HMM.piComb);
      backward(0, nLoci,HMM,X[seq],HMM.piComb);

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
      HMM.pi[l]=0;
      for(int chr=0;chr<nChrom;chr++){
	HMM.pi[l]+=HMM.loci[0].pState[l]/(2.*(nSeq+priorCount));
      }
      HMM.pi[l]/=((double) nChrom);
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
    
    
   
  }

  if(nIter) HMM.write(hmmFileName);
  
  //
  // MCMC
  //


  MCMCSamples.open("MCMCSamples.txt");
  QTLResults.open("QTLResults.txt");
  gHatResults.open("gHatResults.txt");

  //double sig2bPrior=.5,sig2ePrior=20;
  double sig2bPrior=.05,sig2ePrior=5;
  double pi=.99,mu=0;
  int nSamples=4100;
  int nBurnIn=100;
  int FreqToSampleHaplo=100; 
  int outputFreq=4;

  double nusig2e=10,nusig2b=4;
  double sig2e=sig2ePrior,sig2b=sig2bPrior,sig2g;
  double scaleRes=sig2ePrior*(nusig2e-2.);
  double scaleB=sig2bPrior*(nusig2b-2.);
  double nuTilde,uSmp;
  double flipDelta=1./((double) nStates);
  vector<qtlLocus> qtlVec,qtlSumVec;
  vector<long> activeLoci;
  long activePos=0;
  uSmp=u.sample();

  qtlVec.resize(nLoci);
  qtlSumVec.resize(nLoci);
  activeLoci.resize(0);
  for(long i=0;i<nLoci;i++) {
    qtlVec[i].init(nStates,sig2b,pi);
    qtlSumVec[i].init(nStates);
    if(qtlVec[i].active){
      activeLoci.push_back(i);
    }
  }
  long nQTL=activeLoci.size();

  
  vector<double> yDev(nPheno),gHat(nPheno),gHatSum(nPheno,0.0),xVec(nPheno),probClass(nComb);
  
  MCMCSamples << "Sample\tmu\tsig2b\tsig2e\tsig2g" << endl;
  QTLResults << "Loci\tName\tChrom\tPos\tmodelFreq\tb";
  for(int l=0;l<nStates;l++) QTLResults << "\tDelta" << l ;
  for(int l=0;l<nStates;l++) QTLResults << "\thaploFreq" << l ;
  for(int l=0;l<nStates;l++) QTLResults << "\thaploTemplate" << l ;
  QTLResults << endl;
  gHatResults << "ID\tgHat"<< endl;
  vector<double> logPHaplo(nPheno);

  vector<int> dVec(nStates,0),deltaZero(nStates,0);
  vector<vector<int> > deltaStates;
  int s;
  do{
    for(s=0;s<nStates && dVec[s]==1;s++){
      dVec[s]=0;
    }
    if(s<nStates) {
      dVec[s]=1;
      for(int i=0;i<nStates;i++) cout << dVec[i] << " ";
      cout << endl;
      deltaStates.push_back(dVec);
    }
  }while(s<nStates);
  deltaStates.pop_back(); // Drop all out and all in
  int computeLhsV=1;
  vector<double> rhsV,lhsV,lhsVs;
  rhsV.assign(nStates,0.0);
  lhsV.assign(nComb,0.0);
  lhsVs.assign(nStates,0.0);
  vector<vector<double> > lhsVArray(nLoci,lhsV),lhsVsArray(nLoci,lhsVs);
  
  for(int s=0;s<nSamples;s++){

    //
    // Fill in HaploVec
    //
    if((s%FreqToSampleHaplo)==0){
      computeLhsV=1;
      int nFlipped=0;
      vector<int> XHaploNew(nLoci);
      for(int a=0;a<nPheno;a++){
	int seq=phenSeq[a];
	int ranClass,ranClassOld;
	double logPNewvsOld;
	forward(0, nLoci,HMM,X[seq],HMM.piComb);
	vector<double> Pvec(nComb);
	{
	  double Psum=0;
	  double Pval;
	  for(int l=0;l<nComb;l++){
	    Pval=HMM.loci[nLoci-1].f[l];
	    Psum+=Pval;
	    Pvec[l]=Pval;
	  }
	  for(int l=0;l<nComb;l++){
	    Pvec[l]/=Psum;
	  }
	  double uSmp=u.sample();
	  double Thresh=0;
	  for(int l=0;uSmp>Thresh && l< nComb;l++){
	    ranClass=l;
	    Thresh+=Pvec[l];
	  }
	  XHaploNew[nLoci-1]=ranClass;
	  logPNewvsOld=log(Pvec[ranClass]);
	}
	for(long i=nLoci-2;i>=0;i--){
	  double Psum=0;
	  double Pval;
	  for(int l=0;l<nComb;l++){
	    Pval=HMM.loci[i].f[l]*HMM.pComb[l][ranClass];
	    Pvec[l]=Pval;
	    Psum+=Pval;
	  }
	  for(int l=0;l<nComb;l++){
	    Pvec[l]/=Psum;
	  }
	  double uSmp=u.sample();
	  double Thresh=0;
	  
	  for(int l=0;uSmp>Thresh && l< nComb;l++){
	    ranClass=l;
	    Thresh+=Pvec[l];
	  }
	  XHaploNew[i]=ranClass;
	  logPNewvsOld+=log(Pvec[ranClass]);
	}
	double yDevNew=y[a]-mu;
	for(long i=0;i<nQTL;i++){
	  long locus=activeLoci[i];
	  int l=XHaploNew[locus];
	  int I=HMM.stateI[l]; 
	  int J=HMM.stateJ[l];
	  yDevNew-=(qtlVec[locus].delta[I]+qtlVec[locus].delta[J])*qtlVec[locus].b;
	}
	uSmp=u.sample();
	if(s==0 || log(uSmp)<-.5*(yDevNew*yDevNew-yDev[a]*yDev[a])*rinverse[a]/sig2e){
	  for(long i=0;i<nLoci;i++) XHaplo[i][a]=XHaploNew[i];
	  logPHaplo[a]=logPNewvsOld;
	  nFlipped++;
	}
      }
      cout << nFlipped << " of " << nPheno << endl;

      cout << "XHalpo " << XHaplo.size() <<"x" << XHaplo[0].size() << endl;
      
      for(int i=0;i<10;i++) 
	{
	  
	  cout << i ;
	  for(int a=0;a<20; a++){ 
	    cout << " " << setw(1) << HMM.stateI[XHaplo[i][a]] << "/" << HMM.stateJ[XHaplo[i][a]];
	  }
	  cout <<endl;
	}
      cout <<endl;
    }

    //
    // Samplers
    //



    double ssb=scaleB,sse=scaleRes,ssg=0;
    
    nQTL=activeLoci.size();
    yDev=y;
    gHat.assign(nPheno,0.);
    for(int a=0;a<nPheno;a++){
      int seq=phenSeq[a];
      yDev[a]-=mu;
      gHat[a]=0;
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
      
      rhsV.assign(nStates,0.0);
      if(computeLhsV){
	lhsV.assign(nComb,0.0);
	lhsVs.assign(nStates,0.0);
      }
      else{
	lhsV=lhsVArray[i];
	lhsVs=lhsVsArray[i];
      }
      for(int a=0;a<nPheno;a++){
	//double ydev=yDev[a];
	//int seq=phenSeq[a];
	int seqClass=XHaplo[i][a];
	int I=HMM.stateI[seqClass];
	int J=HMM.stateJ[seqClass];
	if(active){
	  double xb=(qtlVec[i].delta[I]+qtlVec[i].delta[J])*b;
	  yDev[a]+=xb;
	}
	rhsV[I]+=yDev[a]*rinverse[a];
	rhsV[J]+=yDev[a]*rinverse[a];
	if(computeLhsV){
	  lhsVs[I]+=rinverse[a];
	  lhsVs[J]+=rinverse[a];
	  lhsV[seqClass]+=2.0*rinverse[a];
	}
      }
      if(computeLhsV){
	lhsVArray[i]=lhsV;
	lhsVsArray[i]=lhsVs;
      }

      int nDeltaStates=deltaStates.size();
      double psum=1.;
      vector<double> AoverIVec(nDeltaStates,log((1.-pi)/pi)-.5*log(sig2b)-log((double) nDeltaStates));
      for(int sd=0;sd<nDeltaStates;sd++){
	double lhs=sig2e/sig2b;
	double rhs=0;
	dVec=deltaStates[sd];
	int IJ=0;
	for(int I=0;I<nStates;I++){
	  if(dVec[I]){
	    rhs+=rhsV[I];
	    lhs+=lhsVs[I];
	  }
	  for(int J=0;J<=I;J++,IJ++){
	    if(dVec[I] && dVec[J]) lhs+=lhsV[IJ];      
	  }
	}
	rhs/=sig2e;
	lhs/=sig2e;
	AoverIVec[sd]+=0.5*(rhs*rhs/lhs-log(lhs));
	AoverIVec[sd]=exp(AoverIVec[sd]);
	psum+=AoverIVec[sd];
      }

      

      uSmp=u.sample();
      uSmp*=psum;


      if(uSmp<1.){
	qtlVec[i].active=0;
	qtlVec[i].b=0;
	qtlVec[i].delta=deltaZero;
      }
      else{
	qtlVec[i].active=1;
	activeLoci.push_back(i);
	uSmp-=1.;

	//nowActive=1;
	
	
	
	//
	// Sample Delta
	//
	int sd;
	for(sd=0; uSmp>AoverIVec[sd]   && sd<nDeltaStates ;sd++){
	  uSmp-=AoverIVec[sd];
	}
	if(sd==nDeltaStates) {
	  cout <<"Delta Sampler Failed!!!"<< uSmp <<   " " << AoverIVec[sd-1] << endl;
	  exit(1);
	}
	dVec=deltaStates[sd];
	qtlVec[i].delta=dVec;
	double sumXY,sumXX;
	sumXY=0;
	sumXX=sig2e/sig2b;
	int IJ=0;
	for(int I=0;I<nStates;I++){
	  if(dVec[I]) {
	    sumXY+=rhsV[I];
	    sumXX+=lhsVs[I];
	  }
	  for(int J=0;J<=I;J++,IJ++){
	    if(dVec[I]&& dVec[J]) sumXX+=lhsV[IJ];
	  }
	}
	
	b=sumXY/sumXX+Z.sample()*sqrt(sig2e/sumXX); //sample b
	
	qtlVec[i].b=b;
	ssb+=b*b;
	for(int a=0;a<nPheno;a++){
	  //int seq=phenSeq[a];
	  int seqClass=XHaplo[i][a];
	  int I=HMM.stateI[seqClass];
	  int J=HMM.stateJ[seqClass];
	  yDev[a]-=(dVec[I]+dVec[J])*b;
	  gHat[a]+=(dVec[I]+dVec[J])*b;
	}
	if(s>=nBurnIn)qtlSumVec[i].updateSum(qtlVec[i]);
      }
      //if(qtlVec[i].active && !nowActive) cout << "Switch locus "<< i << " " << qtlVec[i].active << "=>" << nowActive << endl;

    
     
    
    
   
    }
    computeLhsV=0;
    //update mu;
    double sumXY=0,sumXX=0;
    for(int a=0;a<nPheno;a++){
      sumXY+=rinverse[a]*yDev[a];
      sumXX+=rinverse[a];
    }
    double deltaMu=sumXY/sumXX+Z.sample()*sqrt(sig2e/sumXX);
    for(int a=0;a<nPheno;a++) yDev[a]-=deltaMu;
    mu+=deltaMu;

    //Calc sig2g
    ssg=0;
    for(int a=0;a<nPheno;a++) ssg+=gHat[a]*gHat[a];
    sig2g=ssg/((double) nPheno);

    //update gHatSum
    if(s>=nBurnIn){
      for(int a=0;a<nPheno;a++) {
	gHatSum[a]+=gHat[a];
      }
    }

    //update sig2e;
    nuTilde=((double) nPheno)+nusig2e;

    for(int a=0;a<nPheno;a++) sse+=yDev[a]*yDev[a]*rinverse[a];

    double X2;
    X2=2.*matvec::sgamma(nuTilde/2.0); //Chi-square;
    sig2e=sse/X2;

    //update sig2b;
    nQTL=activeLoci.size();
    nuTilde=((double) nQTL)+nusig2b;
    X2=2.*matvec::sgamma(nuTilde/2.0); //Chi-square;
    sig2b=ssb/X2;
    cout << endl;
    cout << "Sample: " << s << endl;
    cout << setprecision(8) << "Sig2b: " << sig2b << " Sig2e: " << sig2e << " Sig2g: " << sig2g << endl;
    cout << "mu " << mu << endl;
    cout << "number Active: " <<  nQTL << endl;
    for(int i=0; i< 100  ; i++ ) {
      cout <<setw(5) << " " << activeLoci[i];
      if((i %10)==9) cout <<endl;
    }
    cout << endl <<endl;
    if(s%outputFreq==0) MCMCSamples << s << "\t" << mu <<"\t" << sig2b << "\t" <<sig2e <<"\t" << sig2g << endl; 
  }
  double numSampled;
  numSampled=(double)(nSamples-nBurnIn);
  for(long i=0;i<nLoci;i++){
    double numActive=1;
    if(qtlSumVec[i].active) numActive=(double) qtlSumVec[i].active;
    
    QTLResults << i << "\t" << lociMap[i].name << "\t" << lociMap[i].chrom << "\t"<< lociMap[i].pos<<"\t"<< ((double) qtlSumVec[i].active)/numSampled << "\t" << qtlSumVec[i].b/numSampled;
    for(int l=0;l<nStates;l++) QTLResults << "\t" <<  2.*((double) qtlSumVec[i].delta[l])/numActive-1.;
    for(int l=0;l<nStates;l++) QTLResults << "\t" <<HMM.loci[i].pState[l]/(2.*(nSeq+priorCount));
    for(int l=0;l<nStates;l++) QTLResults << "\t" <<HMM.loci[i].e[l];
    QTLResults << endl;
  }
  
  for(int a=0;a<nPheno;a++){
    gHatResults << phenID[a] << "\t" << gHatSum[a]/numSampled << endl;
  }
}




void qtlLocus::updateSum(const qtlLocus &A){
  
  if(A.active){
    active++;
    if(A.b >0){
      b+=A.b;
      for(int i=0;i<A.delta.size();i++){
	delta[i]+=A.delta[i];
      }
    }
    else{
      b-=A.b;
      for(int i=0;i<A.delta.size();i++){
	delta[i]+=1-A.delta[i];
      }
    }
  }
}

void qtlLocus::init(int ns){
  b=0;
  delta.assign(ns,0);
  active=0;

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
    for(int l=0;l<nComb;l++){
      if(HMM.loci[i].newChrom){
	val=picomb[l]*HMM.emit(i,l,X[start]);
      }
      else{
	  val=0;
	  for(int k=0;k<nComb;k++){
	    val+=HMM.loci[i-1].f[k]*HMM.pComb[k][l];
	  }
	  val*=HMM.emit(i,l,X[i]);
      }
      HMM.loci[i].f[l]=val;
      sum+=val;
    }
    for(int l=0;l<nComb;l++) HMM.loci[i].f[l]/=sum;//scale
  }
}

void backward(const long start, const long end,hmm &HMM, const vector<int> &X,const vector<double> &picomb){
  int nComb=HMM.nComb,nStates=HMM.nStates;
  double sum,val;
  for(int l=0;l<nComb;l++) HMM.loci[end-1].b[l]=1;
  for(long i=end-2;i>=start;i--){
    sum=0;
    for(int k=0;k<nComb;k++){
      val=0;
      for(int l=0;l<nComb;l++){
	if(HMM.loci[i+1].newChrom){
	  val+=picomb[l]*HMM.emit(i+1,l,X[i+1])*HMM.loci[i+1].b[l];
	}
	else{
	  val+=HMM.pComb[k][l]*HMM.emit(i+1,l,X[i+1])*HMM.loci[i+1].b[l]; 
	}
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




bool locusMapCompare(locusMap &A,locusMap &B){
  return((A.chrom < B.chrom) || ((A.chrom == B.chrom) && (A.pos < B.pos)));
};

int hmm::write(const string &filename){
  ofstream hmmStream;
  hmmStream.open(filename,ifstream::binary);
  if(hmmStream.fail()) return(1);
  hmmStream.write(reinterpret_cast<char*>(&nStates),sizeof(nStates));
  hmmStream.write(reinterpret_cast<char*>(&nLoci),sizeof(nLoci));
  hmmStream.write(reinterpret_cast<char*>(&pi[0]),nStates*sizeof(pi[0]));
  for(int l=0;l<nStates;l++){
      hmmStream.write(reinterpret_cast<char*>(&p[l][0]),nStates*sizeof(p[l][0]));
  }
  for(long i=0;i<nLoci;i++){
    hmmStream.write(reinterpret_cast<char*>(&loci[i].e[0]),nStates*sizeof(loci[i].e[0]));
  }
  return(0);
}


int hmm::read(const string &filename){
  int ns;
  long nl;
  ifstream hmmStream;
  hmmStream.open(filename,ofstream::binary);
  if(hmmStream.fail()) return(1);
  hmmStream.read(reinterpret_cast<char*>(&ns),sizeof(ns));
  hmmStream.read(reinterpret_cast<char*>(&nl),sizeof(nl));
  if(ns!=nStates || nl != nLoci) return(2); 
  if(pi.size() !=nStates || p.size() != nStates) return(3); 
  hmmStream.read(reinterpret_cast<char*>(&pi[0]),nStates*sizeof(pi[0]));
  for(int l=0;l<nStates;l++){
    if(p[l].size() != nStates) return(4);
    hmmStream.read(reinterpret_cast<char*>(&p[l][0]),nStates*sizeof(p[l][0]));
  }
  if(loci.size() != nLoci) return(5);
  for(long i=0;i<nLoci;i++){
    if(loci[i].e.size() != nStates) return(6);
    hmmStream.read(reinterpret_cast<char*>(&loci[i].e[0]),nStates*sizeof(loci[i].e[0]));
  }
  return(0);
}


 

