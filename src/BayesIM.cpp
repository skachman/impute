#include "impute.h"
#include <limits.h>

int main(int argc,char **argv){

  string genoName,phenoName,mapName;
  mapName="Map.txt";
  genoName="geno.dat";
  phenoName="BWT.dat";
  int freqQTLKb=25; 
  int nStates=4;
  double lambdaKb=500.;
  
  double c=.95; // Not Used
  int nIter=0; // Number of iteration to build HMM
  string baseName,MCMCName,QTLName,gHatName;
  baseName="IM";
  int nSamples=81000;
  int nBurnIn=2000;
  int FreqToSampleHaplo=200; 
  int printFreq=10;
  int outputFreq=8;
  int windowSize=10;  //windowWidth=2*windowSize+1
  double sig2bPrior=.05,sig2ePrior=5;
  double pi=.99,mu=0,inactiveProposal=0.9;
  double nusig2e=10,nusig2b=4;

  
MCMCName=baseName+"_MCMCSamples.txt";
  QTLName=baseName+"_QTLResults.txt";
  gHatName=baseName+"_gHatResults.txt";


  Configuration config;

  if(argc>1){
    if(!config.Load(argv[1])) exit(1);
    config.Get("genoName",genoName);
    config.Get("phenoName",phenoName);
    config.Get("mapName",mapName);
    config.Get("freqQTLKb",freqQTLKb);
    config.Get("nStates",nStates);
    config.Get("lambdaKb",lambdaKb);
    config.Get("nIter",nIter);
    config.Get("baseName",baseName);
    MCMCName=baseName+"_MCMCSamples.txt";
    QTLName=baseName+"_QTLResults.txt";
    gHatName=baseName+"_gHatResults.txt";
    config.Get("MCMCName",MCMCName);
    config.Get("QTLName",QTLName);
    config.Get("gHatName",gHatName);
    config.Get("nSamples",nSamples);
    config.Get("nBurnIn",nBurnIn);
    config.Get("FreqToSampleHaplo",FreqToSampleHaplo);
    config.Get("printFreq",printFreq);
    config.Get("outputFreq",outputFreq);
    config.Get("c",c); //Not Used
    config.Get("windowSize",windowSize); //Not Used
    config.Get("nusig2e",nusig2e);
    config.Get("sig2ePrior",sig2ePrior);
    config.Get("nusig2b",nusig2b);
    config.Get("sig2bPrior",sig2bPrior);
    config.Get("pi",pi);
    config.Get("mu",mu);
    config.Get("inactiveProposal",inactiveProposal);

  }
  double lambda=lambdaKb*1000.;
  int freqQTL=freqQTLKb*1000;


  cout << "Input Parameters" << endl <<endl;
  cout << setw(22) << "genoName = " << " " << genoName << endl;
  cout << setw(22) << "phenoName = " << " " << phenoName << endl;
  cout << setw(22) << "mapName = " << " " << mapName << endl;
  cout << setw(22) << "freqQTLKb = " << " " << freqQTLKb << endl;
  cout << setw(22) << "nStates = " << " " << nStates << endl;
  cout << setw(22) << "lambdaKb = " << " " << lambdaKb << endl;
  cout << setw(22) << "nIter = " << " " << nIter << endl;
  cout << setw(22) << "baseName = " << " " << baseName << endl;
  cout << setw(22) << "MCMCName = " << " " << MCMCName <<endl;
  cout << setw(22) << "QTLName = " << " " << QTLName <<endl;
  cout << setw(22) << "gHatName = " << " " << gHatName <<endl;
  cout << setw(22) << "nSamples = " << " " << nSamples << endl;
  cout << setw(22) << "nBurnIn = " << " " << nBurnIn << endl;
  cout << setw(22) << "FreqToSampleHaplo = " << " " << FreqToSampleHaplo << endl;
  cout << setw(22) << "printFreq = " << " " << printFreq << endl;
  cout << setw(22) << "outputFreq = " << " " << outputFreq << endl;
  //cout << setw(22) << "c = " << " " << c << endl;
  //cout << setw(22) << "windowSize = " << " " << windowSize << endl;
  cout << setw(22) << "nusig2e = " << " " << nusig2e << endl;
  cout << setw(22) << "sig2ePrior = " << " " << sig2ePrior << endl;
  cout << setw(22) << "nusig2b = " << " " << nusig2b << endl;
  cout << setw(22) << "sig2bPrior = " << " " << sig2bPrior << endl;
  cout << setw(22) << "pi = " << " " << pi << endl;
  cout << setw(22) << "mu = " << " " << mu << endl;
  cout << setw(22) << "inactiveProposal = " << " " << inactiveProposal << endl;
  cout << endl << endl;
  
  //  cout << "Maximum value for int:  " << numeric_limits<int>::max() << '\n';
  //  cout << "Maximum value for long: " << numeric_limits<long>::max() << '\n';

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
  string sVal;
  unsigned seed=3434241;
  if(0) seed =(unsigned)(time(NULL));
  DIR *dirPtr = opendir("matvec_trash");
        if (!dirPtr){
                umask(0);
                cout << "making matvec_trash directory" << endl;
                mkdir("./matvec_trash",0777);
        } 
  matvec::SESSION.initialize("matvec_trash");
  //filename=mapName;
  MapFile.open(mapName);
  getline(MapFile,line);
  while(getline(MapFile,line)){
    stringstream linestr(line);
    linestr >> aMapLocus.name;
    linestr >> aMapLocus.chrom;
    linestr >> aMapLocus.pos;
    lociMap.push_back(aMapLocus);
  }
  cout << "Read in map file" << endl << endl;;
  sort(lociMap.begin(),lociMap.end(),locusMapCompare);
  map<string,long> mapOrder;
  for(long i=0;i<lociMap.size();i++) {
    lociMap[i].isSNP=i;
    lociMap[i].isQTL=-1;
    mapOrder[lociMap[i].name]=i;
  }

  vector<int> chromStart,Chrom;
  int lastChrom=lociMap[0].chrom;
  Chrom.push_back(lastChrom);
  chromStart.push_back(0);
  for(long i=1;i<lociMap.size();i++){
    if(lastChrom!=lociMap[i].chrom){
      chromStart.push_back(i);
      Chrom.push_back(lociMap[i].chrom);
      lastChrom=lociMap[i].chrom;
    }
  }
  chromStart.push_back(lociMap.size());
  int nChrom=chromStart.size()-1;

  if(0){
    for(int chr=0;chr<nChrom;chr++) cout << chromStart[chr] << " ";
    cout << endl;
  }
  vector<long> posVector(lociMap.size(),-1);

  //filename="geno.dat";
  Geno.open(genoName);
  cout << "Read in genotype file" << endl;
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
  //cout << seq <<endl;
  cout <<  "Number of genotyped animals: " << X.size() << endl;
  cout <<  "              Number of SNP: " << X[0].size() << endl<<endl;
  //cout << X[X.size()-1].size() << endl;
  double nSeq;
  nSeq=(double) X.size();
  if(0){
    for(int i=0;i<10;i++){
      cout << ID[i] ;
      for(int j=0;j<5;j++){
	cout << " " << X[i][j];
      }
      cout << endl;
    }
  }


  //
  // Read Pheno
  //
  //filename="BWT.dat";
  cout << "Read in phenotype file" << endl;
  Pheno.open(phenoName);
  getline(Pheno,line);
  vector<double> y,Xmu,rinverse;
  vector<int> phenSeq;
  vector<string> phenID;
  vector<string> phenLabels;
  vector<variable_t> varType;
  int nColumns=0;
  int rinversePos=0;
  stringstream leader(line);
  int next;
  int nCovariate=0;
  int nClass=0;
  while(line.length() && (next=line.find_first_not_of(" \t")) !=string::npos){
    line=line.substr(next);
    string label=line;
    int last=line.find_first_of(" \t");
    if(last != string::npos) label=line.substr(0,last);
    phenLabels.push_back(label);
    switch(nColumns){
    case 0:
      varType.push_back(VAR_ID);
      break;
    case 1:
      varType.push_back(VAR_DEP_VAR);
      break;
    default:
      if(label=="rinverse") {
	varType.push_back(VAR_RINVERSE);
	rinversePos=nColumns;
      }
      else{
	char lastChar=label[label.length()-1];
	switch(lastChar){
	case '$': 
	  varType[nColumns]=VAR_CLASS;
	  nClass++;
	  break;
	case '#':
	  varType[nColumns]=VAR_SKIP;
	  break;
	default:
	  varType.push_back(VAR_COVARIATE);
	  nCovariate++;
	}
      }
    }
    if(last != string::npos){
      line=line.substr(last);
    }
    else line="";
    nColumns++;
  }

 
  nColumns=phenLabels.size();
  if(0)
    {
      for(int col=0;col<nColumns;col++){
	cout << phenLabels[col] <<" " <<varType[col];
	if(col==rinversePos) cout  << " rinverse";
	cout << endl;
      }
    }

  vector<  map<string,int> > classMatrix(nClass);
  map<string,int>::iterator classIt;
  vector< vector<int> > posMatrix(nClass);
  vector< vector<double> > valMatrix(nCovariate);
  vector<int > nLevels(nClass,0);
  
  
  int notFound=0;
  while(getline(Pheno,line)){
    int iCl=0;
    int iCv=0;
    
    stringstream linestream(line);
    linestream >> id;
    seqMapIt=seqMap.find(id);
    if(seqMapIt != seqMap.end()){
      phenID.push_back(id);
      phenSeq.push_back(seqMapIt->second);
      if(!rinversePos) rinverse.push_back(1.);
      for(int col=1;col<nColumns;col++){
	switch(varType[col]){
	case VAR_DEP_VAR:
	  linestream >> val;
	  y.push_back(val);
	  break;
	case VAR_RINVERSE:
	  linestream >> val;
	  rinverse.push_back(val);
	  break;
	case VAR_COVARIATE:
	  linestream >> val;
	  valMatrix[iCv++].push_back(val);
	  break;
	case VAR_CLASS:
	  linestream >> sVal;
	  classIt=classMatrix[iCl].find(sVal);
	  if(classIt==classMatrix[iCl].end())classMatrix[iCl][sVal]=nLevels[iCl]++;
	  posMatrix[iCl++].push_back(classMatrix[iCl][sVal]);
	  break;
	default:
	  linestream >> sVal;
	}
	  
      }
    }
  }
  
  int nPheno=y.size();

  cout << "Matched Phenotypes "<< nPheno  << endl << endl;
  if(1){
    for(int i=0;i< 10;i++){
      cout << setw(5) << phenSeq[i] << " ID " << setw(10) <<  ID[phenSeq[i]] << " y " << setw(10) 
	   << setprecision(5) << y[i] ;
      if(rinversePos) cout << " Rinv " << setw(10) << setprecision(5) << rinverse[i] ;
      for(int iCl=0;iCl<nClass;iCl++) cout << " Pos " << setw(2) << posMatrix[iCl][i];
      for(int iCv=0;iCv<nCovariate;iCv++) cout << " Val " << setw(6) << setprecision(4) << valMatrix[iCv][i];
      cout << endl;
    }
    cout << endl;
  }
  
  long nLoci=X[0].size();
  hmmLoci locus;
  hmm HMM(nLoci,nStates,0.5,0.5);
  HMM.lambda=lambda;
  HMM.lociMapPt=&lociMap;
  double priorCount=5;
  vector<double> doubleRow;
  

    double sum=0,pival=0.30,delta=-.05;
    /*   for(int i=0;i<nStates;i++){
      HMM.pi[i]=pival;
      sum+=pival;
      pival+=delta;
      }*/
    for(int i=0;i<nStates;i++) HMM.pi[i]=1./((double) nStates);

    for(int chr=0;chr<nChrom;chr++)  {
      HMM.loci[chromStart[chr]].newChrom=1;
      for(int i=chromStart[chr];i<chromStart[chr+1];i++) HMM.loci[i].pos=lociMap[i].pos;
    }

  
 
    HMM.c=c;
    HMM.resetP();

    string hmmFileName;
    int failCode;
    hmmFileName="HMMIM" + to_string(HMM.nStates) + "_" + to_string(HMM.nLoci) +".bin";
    if((failCode=HMM.read(hmmFileName))!=0) {
      cout << hmmFileName << " failed to be read  with code " << failCode << endl;
      string hmmFileNameAlt;
      hmmFileNameAlt="HMM" + to_string(HMM.nStates) + "_" + to_string(HMM.nLoci) +".bin";
      if((failCode=HMM.read(hmmFileNameAlt))!=0) {
	cout << hmmFileNameAlt << " also failed to be read  with code " << failCode << endl;
      }
    }

    for(int k=0;k<nStates;k++) HMM.pi[k]=1./((double) nStates);
    HMM.initPComb();
  
  
  

  //Estimation
  int nComb=HMM.nComb;
  
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
    /*
    for(int l=0;l<nStates;l++){
      HMM.pi[l]=0;
      for(int chr=0;chr<nChrom;chr++){
	HMM.pi[l]+=HMM.loci[0].pState[l]/(2.*(nSeq+priorCount));
      }
      HMM.pi[l]/=((double) nChrom);
    }
    */
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



  // Add QTL loci;
  int q=0;
  int nQTLLoci;
  int nSNPLoci=nLoci;


  for(int chr=0;chr<nChrom;chr++)  HMM.loci[chromStart[chr]].newChrom=0;
  for(int chr=0;chr<nChrom;chr++){
    aMapLocus.chrom=Chrom[chr];
    //cout << "qpos "<< chromStart[chr] << " "  << lociMap[chromStart[chr]].pos+1 << " " << lociMap[chromStart[chr+1]-1].pos << endl;
    for(long qPos=lociMap[chromStart[chr]].pos+1;qPos<lociMap[chromStart[chr+1]-1].pos;qPos+=freqQTL,q++){
      aMapLocus.name="QTL_" + to_string(chr) + "_" + to_string(qPos);
      aMapLocus.pos=qPos;
      aMapLocus.isSNP=-1;
      aMapLocus.isQTL=q;
      nQTLLoci=q+1;
      lociMap.push_back(aMapLocus);
      //if(q< 100) cout << q <<" " << aMapLocus.name << endl;
    }
  }
  
  sort(lociMap.begin(),lociMap.end(),locusMapCompare);
  chromStart.resize(0);
  lastChrom=lociMap[0].chrom;
  chromStart.push_back(0);
  for(long i=1;i<lociMap.size();i++){
    if(lastChrom!=lociMap[i].chrom){
      chromStart.push_back(i);
      lastChrom=lociMap[i].chrom;
    }
  }
  chromStart.push_back(lociMap.size());
  
  nChrom=chromStart.size()-1;
  nLoci=lociMap.size();
  cout << "  Number of QTL loci: "<< nQTLLoci<< endl;
  cout << "  Number of SNP loci: " << nSNPLoci << endl;
  cout << "Total Number of loci: " <<  nLoci << endl << endl;

  HMM.resize(nLoci);
  for(int chr=0;chr<nChrom;chr++) {
    HMM.loci[chromStart[chr]].newChrom=1;
      for(int i=chromStart[chr];i<chromStart[chr+1];i++) HMM.loci[i].pos=lociMap[i].pos;
  }

  for(int chr=0;chr<nChrom;chr++){
    for(long i=chromStart[chr];i<chromStart[chr+1];i++){
      int q=lociMap[i].isQTL;
      if(q>=0){
	for(int j=i;j>=chromStart[chr] && lociMap[j].pos>lociMap[i].pos-500*1000 ; j--){
	  int qj=lociMap[j].isQTL;
	  if(qj>=0) lociMap[q].start=qj;
	}
	
	for(int j=i;j<chromStart[chr+1] && lociMap[j].pos<lociMap[i].pos+500*1000 ; j++){
	  int qj=lociMap[j].isQTL;
	  if(qj>=0) 
	  lociMap[q].stop=qj;
	}
      }
    }
  }
  
  if(0){
    for(int i=0;i<nQTLLoci;i+=1000){
      cout<< i << " " << lociMap[i].start << " " << lociMap[i].stop << endl;
    }
  }

  XHaplo.resize(nQTLLoci);
  for(long i=0;i<nQTLLoci;i++) XHaplo[i].resize(nPheno); 
  
  //
  // MCMC
  //


  MCMCSamples.open(MCMCName);
  QTLResults.open(QTLName);
  gHatResults.open(gHatName);

  //double sig2bPrior=.5,sig2ePrior=20;
  
  
  //int nSamples=40;
  //int nBurnIn=2;
  
  double sig2e=sig2ePrior,sig2b=sig2bPrior,sig2g;
  double scaleRes=sig2ePrior*(nusig2e-2.);
  double scaleB=sig2bPrior*(nusig2b-2.);
  double nuTilde,uSmp;
  double flipDelta=1./((double) nStates);
  vector<qtlLocus> qtlVec,qtlSumVec;
  vector<int> windowSumVec,windowVec;
  vector<long> activeLoci;
  long activePos=0;
  uSmp=u.sample();

  qtlVec.resize(nQTLLoci);
  qtlSumVec.resize(nQTLLoci);
  windowSumVec.assign(nQTLLoci,0);
  windowVec.resize(nQTLLoci);
  activeLoci.resize(0);
  for(long i=0;i<nQTLLoci;i++) {
    qtlVec[i].init(nStates,sig2b,pi);
    qtlSumVec[i].init(nStates);
    if(qtlVec[i].active){
      activeLoci.push_back(i);
    }
  }
  long nQTL=activeLoci.size();

  
  vector<double> yDev(nPheno),gHat(nPheno),gHatSum(nPheno,0.0),xVec(nPheno),probClass(nComb);
  
  MCMCSamples << "Sample\tmu\tsig2b\tsig2e\tsig2g" << endl;
  QTLResults << "Loci\tName\tChrom\tPos\tmodelFreq\tb\twindowFreq";
  for(int l=0;l<nStates;l++) QTLResults << "\tDelta" << l ;
  //for(int l=0;l<nStates;l++) QTLResults << "\thaploFreq" << l ;
  //for(int l=0;l<nStates;l++) QTLResults << "\thaploTemplate" << l ;
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
      if(0){
	for(int i=0;i<nStates;i++) cout << dVec[i] << " ";
	cout << endl;
      }
      deltaStates.push_back(dVec);
    }
  }while(s<nStates);
  deltaStates.pop_back(); // Drop all out and all in
  int computeLhsV=1;
  vector<double> rhsV,lhsV,lhsVs;
  rhsV.assign(nStates,0.0);
  lhsV.assign(nComb,0.0);
  lhsVs.assign(nStates,0.0);
  vector<vector<double> > lhsVArray(nQTLLoci,lhsV),lhsVsArray(nQTLLoci,lhsVs);
  


  for(int s=0;s<nSamples;s++){

    //
    // Fill in HaploVec
    //
    if((s%FreqToSampleHaplo)==0){
      cout << endl << "Sampling haplotypes" << endl;
      computeLhsV=1;
      int nFlipped=0;
      vector<int> XHaploNew(nQTLLoci);
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
	  int q=lociMap[nLoci-1].isQTL;
	  if(q>=0)XHaploNew[q]=ranClass;
	  logPNewvsOld=log(Pvec[ranClass]);
	}
	for(long i=nLoci-2;i>=0;i--){
	  double Psum=0;
	  double Pval;
	  double S=(double) nStates,lambda=HMM.lambda;
	  lambda*=(S-1.)/S;
	  vector<double> transSingle(2); 
	  double eMult=0.;
	  if(! HMM.loci[i+1].newChrom) eMult=exp(-((double)(HMM.loci[i+1].pos-HMM.loci[i].pos))/lambda);
	  transSingle[1]=(1.+(S-1.)*eMult)/S;
	  transSingle[0]=(1.-eMult)/S;

	  for(int l=0;l<nComb;l++){

	  int I1=HMM.stateI[l]; 
	  int J1=HMM.stateJ[l];
	  int I2=HMM.stateI[ranClass]; 
	  int J2=HMM.stateJ[ranClass];
	  Pval=transSingle[I1==I2]*transSingle[J1==J2];
	  if(I2!=J2) Pval+=transSingle[I1==J2]*transSingle[J1==I2];
	  Pval*=HMM.loci[i].f[l];
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
	  int q=lociMap[i].isQTL;
	  if(q>=0)XHaploNew[q]=ranClass;
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
	  for(long i=0;i<nQTLLoci;i++) XHaplo[i][a]=XHaploNew[i];
	  logPHaplo[a]=logPNewvsOld;
	  nFlipped++;
	}
      }
      cout << "Accepted "<< nFlipped << " of " << nPheno << " proposed haplotypes." << endl;

      //cout << "XHalpo " << XHaplo.size() <<"x" << XHaplo[0].size() << endl;
      cout << endl << "Haplotypes " <<endl <<endl;;
      for(int a=0;a<10;a++) 
	{
	  
	  cout << a ;
	  for(int i=0;i<20; i++){ 
	    cout << " " << setw(1) << HMM.stateI[XHaplo[i][a]] << "/" << HMM.stateJ[XHaplo[i][a]];
	  }
	  cout << "   ";
	  for(int i=0;i<20; i++){ 
	    cout << " " << setw(1) << HMM.stateI[XHaplo[nQTLLoci-30+i][a]] << "/" << HMM.stateJ[XHaplo[nQTLLoci-30+i][a]];
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
    int nRejectI2A=0;
    int nRejectA2I=0;
    for(long i=0;i<nQTLLoci;i++){
      int active=qtlVec[i].active;
      int nowActive=0;
      int proposeActive=0;
      
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
      uSmp=u.sample();
      if(uSmp>inactiveProposal) proposeActive=1;
      if(proposeActive || active || computeLhsV){
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
      }
      if(active || proposeActive){
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
     
      

	int accept=0;
	if(proposeActive && active) accept=1;
	else{
	  if(proposeActive && u.sample() < (1.-1./psum)/(1-inactiveProposal)){
	    accept=1;
	  }
	  if(active && u.sample() < (1.-inactiveProposal)/(1.-1./psum)){
	    accept=1;
	  }
	}
	if(!accept) {
	  if(active) nRejectA2I++;
	  else nRejectI2A++;
	}
	if(accept || active){	
	  
	  
	  if(accept && !proposeActive){
	    qtlVec[i].active=0;
	    qtlVec[i].b=0;
	    qtlVec[i].delta=deltaZero;
	  }
	  else{
	    
	    uSmp=u.sample();
	    uSmp*=psum-1.;
	    qtlVec[i].active=1;
	    activeLoci.push_back(i);
	    //uSmp-=1.;
	    
	    nowActive=1;
	    
	    
	    
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
	    
	  }
	}
      }
      
      //if(qtlVec[i].active && !nowActive) cout << "Switch locus "<< i << " " << qtlVec[i].active << "=>" << nowActive << endl;
      
      if(qtlVec[i].active){
	
      }
      
      if(s>=nBurnIn)qtlSumVec[i].updateSum(qtlVec[i]);
      
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
    double sumg=0;
    for(int a=0;a<nPheno;a++) {
      ssg+=gHat[a]*gHat[a];
      sumg+=gHat[a];
    }
    sumg/=(double) nPheno;
    sig2g=ssg/((double) nPheno)-(sumg*sumg);

    //update gHatSum
    if(s>=nBurnIn){
      for(int a=0;a<nPheno;a++) {
	gHatSum[a]+=gHat[a];
      }
    }

    //update sig2e;
    nuTilde=((double) nPheno)+nusig2e;
    //cout << "sse " << sse << endl; 
    for(int a=0;a<nPheno;a++)  sse+=yDev[a]*yDev[a]*rinverse[a];
     
    

    double X2;
    X2=2.*matvec::sgamma(nuTilde/2.0); //Chi-square;
    sig2e=sse/X2;
    //cout << "Sig2E " << sse << " " << nuTilde << " " << X2 << " " << sig2e << endl; 

    //update sig2b;
    nQTL=activeLoci.size();
    nuTilde=((double) nQTL)+nusig2b;
    X2=2.*matvec::sgamma(nuTilde/2.0); //Chi-square;
    sig2b=ssb/X2;
    if(s>=nBurnIn){
      windowVec.assign(nQTLLoci,0);
      for(int i=0;i<nQTLLoci;i++){
	if(qtlVec[i].active){
	  int start=lociMap[i].start;
	  int end=lociMap[i].stop;
	  for(int j=start;j<=end;j++) windowVec[j]=1;
	}
      }
      for(int i=0;i<nQTLLoci;i++) {
	if(windowVec[i]) windowSumVec[i]++;
      }
    }
  

    if((s % printFreq)==0){
      cout << endl;
      cout << "Sample: " << s << endl;
      cout << setprecision(8) << "Sig2b: " << sig2b << " Sig2e: " << sig2e << " Sig2g: " << sig2g << endl;
      cout << "mu " << mu << endl;
      cout << "number Active: " <<  nQTL << " Number Rejected A2I: " << nRejectA2I<< " Number Rejected I2A: " << nRejectI2A << endl;
      for(int i=0; i< 100  ; i++ ) {
	cout  << " " << setw(6) << activeLoci[i];
	if((i %10)==9) cout <<endl;
      }
      cout << endl;
      for(int i=0; i< 100  ; i++ ) {
	cout  << " " << setw(6) << activeLoci[nQTL-100+i];
	if((i %10)==9) cout <<endl;
      }
      cout << endl <<endl;
      if(s%outputFreq==0) MCMCSamples << s << "\t" << mu+sumg <<"\t" << sig2b << "\t" <<sig2e <<"\t" << sig2g << endl; 
    }
  }
  double numSampled;
  numSampled=(double)(nSamples-nBurnIn);
  for(long i=0;i<nLoci;i++){
    int q=lociMap[i].isQTL;
    if(q>=0){
      //cout << i << " " << q << endl;
      double numActive=1;
      if(qtlSumVec[q].active) numActive=(double) qtlSumVec[q].active;
      
      QTLResults << q << "\t" << lociMap[i].name << "\t" << lociMap[i].chrom << "\t"<< lociMap[i].pos<<"\t"<< ((double) qtlSumVec[q].active)/numSampled << "\t" << qtlSumVec[q].b/numSampled << "\t" << ((double) windowSumVec[q])/numSampled;
      for(int l=0;l<nStates;l++) QTLResults << "\t" <<  2.*((double) qtlSumVec[q].delta[l])/numActive-1.;
      //for(int l=0;l<nStates;l++) QTLResults << "\t" <<HMM.loci[i].pState[l]/(2.*(nSeq+priorCount));
      //for(int l=0;l<nStates;l++) QTLResults << "\t" <<HMM.loci[i].e[l];
      QTLResults << endl;
    }
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
  vector<locusMap> *lociMapPt=HMM.lociMapPt;
  double S=(double) nStates,lambda=HMM.lambda;
  lambda*=(S-1.)/S;
  double sum=0,val;
  vector<double> transSingle(2); 
  for(int l=0;l<nComb;l++){
    val=picomb[l];
    int snpLoc=(*lociMapPt)[start].isSNP;
    if(snpLoc>=0) val*=HMM.emit(snpLoc,l,X[snpLoc]);
    HMM.loci[start].f[l]=val;
    sum+=val;

  }
  for(int l=0;l<nComb;l++) HMM.loci[start].f[l]/=sum;//scale

  for(int i=start+1;i<end;i++){
    
    double eMult=exp(-((double)(HMM.loci[i].pos-HMM.loci[i-1].pos))/lambda);
    transSingle[1]=(1.+(S-1.)*eMult)/S;
    transSingle[0]=(1.-eMult)/S;
    if(0 ) cout << i << " " << (HMM.loci[i].pos-HMM.loci[i-1].pos)/lambda <<  " " << eMult << " " << transSingle[0] << " " << transSingle[1] << endl;
    sum=0;
    int l=0;
    for(int I1=0;I1<nStates;I1++){
      for(int J1=0;J1<=I1;J1++,l++){
	if(HMM.loci[i].newChrom){
	  val=picomb[l];

	  int snpLoc=(*lociMapPt)[i].isSNP;
	  if(snpLoc>=0) val*=HMM.emit(snpLoc,l,X[snpLoc]);
	}
	else{
	  val=0;
	  int k=0;
	  for(int I2=0;I2<nStates;I2++){
	    for(int J2=0;J2<=I2;J2++,k++){
	      val+=HMM.loci[i-1].f[k]*transSingle[I1==I2]*transSingle[J1==J2];
	      if(I2!=J2) val+=HMM.loci[i-1].f[k]*transSingle[I1==J2]*transSingle[J1==I2];
	    }
	  }
	  int snpLoc=(*lociMapPt)[i].isSNP;
	  if(snpLoc>=0) val*=HMM.emit(snpLoc,l,X[snpLoc]);
	}
	HMM.loci[i].f[l]=val;
	sum+=val;
      }
    }
    for(int l=0;l<nComb;l++) HMM.loci[i].f[l]/=sum;//scale
  }
}



void backward(const long start, const long end,hmm &HMM, const vector<int> &X,const vector<double> &picomb){
  int nComb=HMM.nComb,nStates=HMM.nStates;
  vector<locusMap> *lociMapPt=HMM.lociMapPt;
  double S=(double) nStates,lambda=HMM.lambda;
  lambda*=(S-1.)/S;
  vector<double> transSingle(2); 
  double sum,val;
  for(int l=0;l<nComb;l++) HMM.loci[end-1].b[l]=1;
  for(long i=end-2;i>=start;i--){
    double eMult=exp(-((double)(HMM.loci[i+1].pos-HMM.loci[i].pos))/lambda);
    transSingle[1]=(1.+(S-1.)*eMult)/S;
    transSingle[0]=(1.-eMult)/S;
    sum=0;
    int k=0;
    for(int I1=0;I1<nStates;I1++){
      for(int J1=0;J1<=I1;J1++,k++){
	val=0;
	int l=0;
	for(int I2=0;I2<nStates;I2++){
	  for(int J2=0;J2<=I2;J2++,l++){
	    int snpLoc=(*lociMapPt)[i+1].isSNP;
	    double eProb=1.;
	    if(snpLoc>=0) eProb=HMM.emit(snpLoc,l,X[snpLoc]);
	    if(HMM.loci[i+1].newChrom){
	      val+=picomb[l]*HMM.loci[i+1].b[l]*eProb;	      
	    }
	    else{
	      val+=transSingle[I1==I2]*transSingle[J1==J2]*HMM.loci[i+1].b[l]*eProb; 
	      if(I2 !=J2)val+=transSingle[I1==J2]*transSingle[J1==I2]*HMM.loci[i+1].b[l]*eProb; 
	    }
	  }
	}
	HMM.loci[i].b[k]=val;
	sum+=val;
      }
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


 


