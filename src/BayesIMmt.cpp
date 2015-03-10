#include "impute.h"
#include <limits.h>


uniform_real_distribution<double> u(0.,1.);
normal_distribution<double> Z(0.,1.);
mt19937 gen;

int main(int argc,char **argv){

  string genoName,phenoName,mapName,SNPName,QTLMapName;
  string randomString;
  bool QTLMap;
  int nRandom=0;
  int nTraits=2;
  vector<Matrix2d> sig2rPrior,sig2r;
  vector <double> nusig2r;
  map<string,int> sig2rVar;
  int failCode;
  mapName="Map.txt";
  genoName="geno.dat";
  int freqQTLKb=25; 
  int nStates=4;
  double lambdaKb=500.;
  double rho=.5;
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
  Matrix2d sig2bPrior,sig2ePrior;
  double pi=.99,inactiveProposal=0.9;
  Vector2d mu;
  double nusig2e=10,nusig2b=4;
  int deltaSampler=2;
  
  int enableSwapActive=0;
  
  double piPrior=0.99,piPriorCount=0;
  MCMCName=baseName+"_MCMCSamples.txt";

  QTLName=baseName+"_QTLResults.txt";
  gHatName=baseName+"_gHatResults.txt";


  int nQTLClasses=0;
  vector<int> QTLClassVec,QTLClassCount;
  map<string,int> QTLClasses;
  string QTLClassStr;
  vector<Matrix2d> sig2bPriorVec,sig2bVec;
  vector<double> nusig2bVec,piPriorVec,piPriorCountVec,piClassVec;

  
  Configuration config;
  string QTLString;
  if(argc>1){
    if(!config.Load(argv[1])) exit(1);
    
    config.Get("nTraits",nTraits);
    config.Get("nusig2b",nusig2b);
    string varianceString;
    config.Get("sig2bPrior",varianceString);
    sig2bPrior.resize(nTraits,nTraits);
    {
      double val;
      stringstream line(varianceString);
      for(int i=0;i<nTraits;i++){
	for(int j=0;j<i;j++){
	  line >> val;
	  sig2bPrior(i,j)=val;
	  sig2bPrior(j,i)=val;
	}
	  line >> val;
	  sig2bPrior(i,i)=val;
      }
    }
    bool setPi=config.Get("pi",pi);
    if(!config.Get("piPrior",piPrior)) piPrior=pi;
    if(!setPi)pi=piPrior;
    config.Get("piPriorCount",piPriorCount);
    config.Get("rho",rho);
    config.Get("genoName",genoName);
    config.Get("phenoName",phenoName);
    config.Get("mapName",mapName);
    QTLMap=config.Get("QTLMapName",QTLMapName);
    if(config.Get("QTLClasses",QTLString)){
      int next;
      Matrix2d sigtmp(nTraits,nTraits);
      while(QTLString.length() && (next=QTLString.find_first_not_of(" \t")) !=string::npos){
	double val;
	QTLString=QTLString.substr(next);
	string label=QTLString,varLabel;
	int last=QTLString.find_first_of(" \t");
	if(last != string::npos) label=QTLString.substr(0,last);
	QTLClasses[label]=nQTLClasses;
	QTLClassCount.push_back(0);
	if(!config.Get("QTLnu_"+label,val)) val=nusig2b;
	nusig2bVec.push_back(val);
	if(!config.Get("sig2bPrior_"+label,varianceString)) {
	  sig2bPriorVec.push_back(sig2bPrior);
	}
	else{
	  stringstream line(varianceString);
	  for(int i=0;i<nTraits;i++){
	    for(int j=0;j<i;j++){
	      line >> val;
	      sigtmp(i,j)=val;
	      sigtmp(j,i)=val;
	    }
	    line >> val;
	    sigtmp(i,i)=val;
	  }
	  sig2bPriorVec.push_back(sigtmp);
	}
	sig2bVec.push_back(sigtmp);
	if(!config.Get("piPrior_"+label,val)) val=piPrior;
	piPriorVec.push_back(val);
	if(!config.Get("piPriorCount_"+label,val)) val=piPriorCount;
	piPriorCountVec.push_back(val);
	piClassVec.push_back(pi);
	nQTLClasses++;

    if(last != string::npos){
      QTLString=QTLString.substr(last);
    }
    else QTLString="";
      }
    }
    if(!QTLMap){
      QTLClasses[""]=nQTLClasses;
	QTLClassCount.push_back(0);
	nusig2bVec.push_back(nusig2b);
	sig2bPriorVec.push_back(sig2bPrior);
	sig2bVec.push_back(sig2bPrior);
	piPriorVec.push_back(piPrior);
	piPriorCountVec.push_back(piPriorCount);
	piClassVec.push_back(pi);
	nQTLClasses++;
    }
    
    config.Get("freqQTLKb",freqQTLKb);
    config.Get("nStates",nStates);
    config.Get("lambdaKb",lambdaKb);
    config.Get("nIter",nIter);
    config.Get("baseName",baseName);
    MCMCName=baseName+"_MCMCSamples.txt";
    QTLName=baseName+"_QTLResults.txt";
    gHatName=baseName+"_gHatResults.txt";
    SNPName=baseName+"_SNPResults.txt";
    config.Get("MCMCName",MCMCName);
    config.Get("QTLName",QTLName);
    config.Get("gHatName",gHatName);
    config.Get("nSamples",nSamples);
    config.Get("nBurnIn",nBurnIn);
    config.Get("deltaSampler",deltaSampler);
    config.Get("enableJiggle",enableSwapActive);
    config.Get("enableSwapActive",enableSwapActive);
    if(enableSwapActive && (deltaSampler==1 || deltaSampler==2)) deltaSampler+=2;
    enableSwapActive=0;
    config.Get("FreqToSampleHaplo",FreqToSampleHaplo);
    config.Get("printFreq",printFreq);
    config.Get("outputFreq",outputFreq);
    config.Get("c",c); //Not Used
    config.Get("windowSize",windowSize); //Not Used
    config.Get("nusig2e",nusig2e);

    config.Get("sig2ePrior",varianceString);
    sig2ePrior.resize(nTraits,nTraits);
    {
      stringstream line(varianceString);
      double val;
      for(int i=0;i<nTraits;i++){
	for(int j=0;j<i;j++){
	  line >> val;
	  sig2ePrior(i,j)=val;
	  sig2ePrior(j,i)=val;
	}
	line >> val;
	sig2ePrior(i,i)=val;
      }
    }
    string muString;
    mu=Vector2d::Zero();
    if(config.Get("mu",muString)){
      double val;
      stringstream line(muString);
      for(int t=0;t<nTraits;t++) {
	line >>mu[t];
      }
    }
    
    config.Get("inactiveProposal",inactiveProposal);
    if(config.Get("randomEffects",randomString)){
      int next;
      
      while(randomString.length() && (next=randomString.find_first_not_of(" \t")) !=string::npos){
	double val;
	randomString=randomString.substr(next);
	string label=randomString,varLabel;
	int last=randomString.find_first_of(" \t");
	if(last != string::npos) label=randomString.substr(0,last);
	sig2rVar[label]=nRandom;
	if(!config.Get("nu_"+label,val)) val=4;
	nusig2r.push_back(val);
	if(!config.Get("sig2Prior_"+label,varianceString)) {
	  sig2rPrior.push_back(sig2ePrior);
	  sig2r.push_back(sig2ePrior);
	}
	else{
	  Matrix2d sigtmp(nTraits,nTraits);
	  stringstream line(varianceString);
	  double val;
	  for(int i=0;i<nTraits;i++){
	    for(int j=0;j<i;j++){
	      line >> val;
	      sigtmp(i,j)=val;
	      sigtmp(j,i)=val;
	    }
	    line >> val;
	    sigtmp(i,i)=val;
	  }
	  sig2rPrior.push_back(sigtmp);
	  sig2r.push_back(sigtmp);
	}
	nRandom++;

    if(last != string::npos){
      randomString=randomString.substr(last);
    }
    else randomString="";
      }
    }

    }
  double lambda=lambdaKb*1000.;
  int freqQTL=freqQTLKb*1000;

  cout << argv[0] << ": Version 2.0" << " Mar. 10, 2015" << endl << endl;
  
  cout << "Input Parameters" << endl <<endl;
  
  cout << setw(22) << "nTraits = " << " " << nTraits << endl;
  cout << setw(22) << "genoName = " << " " << genoName << endl;
  cout << setw(22) << "phenoName = " << " " << phenoName << endl;
  cout << setw(22) << "mapName = " << " " << mapName << endl;
  if(QTLMap) {
    cout << setw(22) << "QTLMapName = " << " " << QTLMapName << endl;
    cout << setw(22) << "QTLClasses =" ;
    for (auto it=QTLClasses.begin(); it!=QTLClasses.end(); ++it){ 
      cout <<" " << it->first; 
    }
    cout << endl;
    for (auto it=QTLClasses.begin(); it!=QTLClasses.end(); ++it){
      string label=it->first;
      int c=it->second;
      cout << setw(22) << "QTLnu_"+label+" = "  << " " << nusig2bVec[c] << endl; 
      cout << setw(22) << "sig2bPrior_"+label+" = "  ;
      for(int i=0;i<nTraits;i++){
	for(int j=0;j<=i;j++){
	  cout << " " << sig2bPriorVec[c](i,j);
	}
      }
      cout << endl; 
      cout << setw(22) << "piPrior_"+label+" = "  << " " << piPriorVec[c] << endl; 
      cout << setw(22) << "piPriorCount_"+label+" = "  << " " << piPriorVec[c] << endl; 
    }
  }
  
  cout << setw(22) << "freqQTLKb = " << " " << freqQTLKb << endl;
  cout << setw(22) << "nStates = " << " " << nStates << endl;
  cout << setw(22) << "lambdaKb = " << " " << lambdaKb << endl;
  cout << setw(22) << "nIter = " << " " << nIter << endl;
  cout << setw(22) << "baseName = " << " " << baseName << endl;
  cout << setw(22) << "MCMCName = " << " " << MCMCName <<endl;
  cout << setw(22) << "QTLName = " << " " << QTLName <<endl;
  cout << setw(22) << "gHatName = " << " " << gHatName <<endl;
  cout << setw(22) << "SNPName = " << " " << SNPName <<endl;
  cout << setw(22) << "nSamples = " << " " << nSamples << endl;
  cout << setw(22) << "nBurnIn = " << " " << nBurnIn << endl;
  cout << setw(22) << "deltaSampler = " << " " << deltaSampler << " # 1=Use full sampler, 2=Use nState sampler"<<endl;
  cout <<  "                         # 3=Use full sampler with locus swap, 2=Use nState sampler with locus swap"<<endl;
  // cout << setw(22) << "enableSwapActive = " << " " << enableSwapActive << " # 0=Don't swap active locus, 1=Enable swap active locus sampler"<<endl;
  cout << setw(22) << "FreqToSampleHaplo = " << " " << FreqToSampleHaplo << endl;
  cout << setw(22) << "printFreq = " << " " << printFreq << endl;
  cout << setw(22) << "outputFreq = " << " " << outputFreq << endl;
  //cout << setw(22) << "c = " << " " << c << endl;
  //cout << setw(22) << "windowSize = " << " " << windowSize << endl;
  cout << setw(22) << "nusig2e = " << " " << nusig2e << endl;
  cout << setw(22) << "sig2ePrior = ";
  for(int i=0;i<nTraits;i++){
    for(int j=0;j<=i;j++){
      cout << " " << sig2ePrior(i,j);
    }
  }
  cout << endl; 
  cout << setw(22) << "nusig2b = " << " " << nusig2b << endl;
  cout << setw(22) << "sig2bPrior = " ;
  for(int i=0;i<nTraits;i++){
    for(int j=0;j<=i;j++){
      cout << " " << sig2bPrior(i,j);
    }
  }
  cout << endl; 
  
  if(nRandom){
    cout << setw(22) << "randomEffects =" ;
    for (std::map<string,int>::iterator it=sig2rVar.begin(); it!=sig2rVar.end(); ++it){ 
      cout <<" " << it->first; 
    }
    cout << endl;
    for (std::map<string,int>::iterator it=sig2rVar.begin(); it!=sig2rVar.end(); ++it){
      string label=it->first;
      int rEffect=it->second;
      cout << setw(22) << "nu_"+label+" = "  << " " << nusig2r[rEffect] << endl; 
      cout << setw(22) << "sig2Prior_"+label+" = ";
      for(int i=0;i<nTraits;i++){
	for(int j=0;j<=i;j++){
	  cout << " " << sig2rPrior[rEffect](i,j);
	}
      }
      cout << endl; 
    }
  }
  cout << setw(22) << "pi = " << " " << pi << endl;
  cout << setw(22) << "piPrior = " << " " << piPrior << endl;
  cout << setw(22) << "piPriorCount = " << " " << piPriorCount << " # 0=Fixed"<< endl;
  
  cout << setw(22) << "rho = " << " " << rho << " # an active QTL is active for a specific trait" << endl;
  cout << setw(22) << "mu = " ;
  printVector(cout,mu," ");
  cout << endl;
  cout << setw(22) << "inactiveProposal = " << " " << inactiveProposal << endl;
  cout << endl << endl;

  if(nTraits!=2) {
    cout << "Ntraits is not equal to 2!" << endl;
    exit(835);
  }
  
  //cout << "Maximum value for int:  " << numeric_limits<int>::max() << '\n';
  //cout << "Maximum value for long: " << numeric_limits<long>::max() << '\n';

  //matvec::UniformDist u;
  //uniform_real_distribution<double> u(0.,1.);
  //matvec::NormalDist Z;
  //normal_distribution<double> Z(0.,1.);
  idmap seqMap,phenMap;
  idmap::iterator seqMapIt,phenMapIt;

  ifstream Geno,Pheno,MapFile,QTLMapFile;
  ofstream MCMCSamples,QTLResults,gHatResults,snpResults;
  string filename;
  string line,id;
  vector<vector<int> > X;
  vector<vector<int> > XHaplo;
  vector<vector<haploLocus> > haploProb;
  vector<locusMap> lociMap;
  locusMap aMapLocus;

  vector<string > ID;
  double val;
  int iVal;
  string sVal;
  
  Vector2d zVec;
  unsigned seed=3434241;
  random_device rd;
  //mt19937 gen;
  gen.seed(rd());
  //if(0) gen(seed);
  /*
    DIR *dirPtr = opendir("matvec_trash");
        if (!dirPtr){
                umask(0);
                cout << "making matvec_trash directory" << endl;
                mkdir("./matvec_trash",0777);
        } 
  matvec::SESSION.initialize("matvec_trash");
 */
  //filename=mapName;
  MapFile.open(mapName);
  if(MapFile.fail()) {
    cout << "Failed to open map file " + mapName << endl;
    exit(100);
  }
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
    //if(i<10) cout << lociMap[i].name << " " << i << endl;
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
  vector<long> posVector;

  //filename="geno.dat";
  Geno.open(genoName);
  if(Geno.fail()) {
    cout << "Failed to open genotype file " + genoName << endl;
    exit(101);
  }
  cout << "Read in genotype file" << endl;
  getline(Geno,line);
  stringstream linestr(line);
  linestr >> id;
  long ipos=0;
  while(linestr >> id){
    if(mapOrder.find(id)!=mapOrder.end()){
      posVector.push_back(mapOrder[id]);
      iVal=mapOrder[id];
    }
    else{
      iVal=-1;
      posVector.push_back(-1);
    }
    //if(ipos<10) cout << ipos << " "  << id << " " << iVal << endl;
    //if(iVal >-1 && iVal < 10) cout << ipos << " "  << id <<  " " << iVal << " " << posVector[ipos] <<"*" <<endl;
    ipos++;
  }
  int nPos=posVector.size();
  int seq=0;
  string binaryGenoFileName;
  binaryGenoFileName= genoName + ".bin";
  if(readX(binaryGenoFileName, X,ID,seqMap)!=0){
    while(getline(Geno,line)){
      stringstream linestr(line);
      linestr >> id;
      ID.push_back(id);
      seqMap[id]=seq++;
      vector<int> row(mapOrder.size(),-2);
      //    X.push_back(row);
      long i=0;
      if((seq%200)==0) cout << seq << endl;
      while(linestr >> val){
	if(i >=nPos) {
	  cout <<"Genotype record is longer than expected, check record with ID: " << id << " " << i << endl;
	  exit(101);
	}
	iVal=-1;
	if(val==10 || val==0 || val==-10){
	  iVal=val+10;
	  iVal/=10;
	}
	
	//  if(seq < 1 && i<5) cout << i << " " << val << " " << iVal << " " <<posVector[i] << endl;
	//if(seq < 5 && posVector[i]<5 && posVector[i]> -1) cout << i << " " << val << " " << iVal << " " <<posVector[i] <<"*"<< endl;
	
	if(posVector[i]>-1) row[posVector[i]]=iVal;
	//X.back()[posVector[i]]=iVal;
	i++;
      }
      X.push_back(row);
      
    }
    writeX(binaryGenoFileName,X,ID);
  }
  
  cout << endl;
  for(int i=0;i<10;i++){
    cout << setw(10) << ID[i];
    for(int j=0;j<10;j++){
      cout << " " <<setw(2) << X[i][j];
    } 
    cout <<"    ";
    for(int j=10;j>0;j--){
      
      cout << " " <<setw(3) << X[i][mapOrder.size()-j];
    }
    cout << endl;
  }
  cout << endl;
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
  if(Pheno.fail()){
    cout << "Failed to open phenotype file " + phenoName << endl;
    exit(102);
  }
  getline(Pheno,line);
  vector<Vector2d> y;
  Vector2d ytmp(nTraits);
  Vector2d rinvtmp(nTraits);
  for(int t=0;t<nTraits;t++) rinvtmp[t]=1.;
  vector<double> Xmu;
  vector<Vector2d> rinverse;
  vector<int> phenSeq;
  vector<string> phenID;
  vector<string> phenLabels;
  vector<string> className,covName;
  vector<int > classRandom;
  vector<variable_t> varType;
  int nColumns=0;
  vector<int> rinversePos;
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
    if(nColumns <= nTraits){
      if(nColumns==0){
	varType.push_back(VAR_ID);
      }
      else{
	varType.push_back(VAR_DEP_VAR);
      }
    }
    else {
      if(label=="rinverse") {
	varType.push_back(VAR_RINVERSE);
	rinversePos.push_back(nColumns);
      }
      else{
	char lastChar=label[label.length()-1];
	switch(lastChar){
	case '$': 
	  varType.push_back(VAR_CLASS);
	  iVal=-1;
	  if(sig2rVar.find(label)!=sig2rVar.end()) {
	    iVal=sig2rVar[label];
	  }
	  classRandom.push_back(iVal);
	  className.push_back(label);
	  nClass++;
	  break;
	case '#':
	  varType.push_back(VAR_SKIP);
	  break;
	default:
	  varType.push_back(VAR_COVARIATE);
	  covName.push_back(label);
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

  int nrinv=rinversePos.size();
  if(nrinv!=0 && nrinv !=nTraits){
    cout << "The number of rinverse columns " << nrinv << " must be equal to 0 or the number of traits." << endl;
    exit(3431);
  }
  nColumns=phenLabels.size();
  

  vector<  map<string,int> > classMatrix(nClass);
  vector<vector<string> > classValues(nClass);
  
  map<string,int>::iterator classIt;
  vector< vector<int> > posMatrix(nClass);
  vector< vector<double> > valMatrix(nCovariate);
  vector<vector<Vector2d> > betaClass(nClass);
  vector<Vector2d> betaCov(nCovariate,Vector2d::Zero());
  vector<int > nLevels(nClass,0);
  

  
  int notFound=0;
  while(getline(Pheno,line)){
    int iCl=0;
    int iCv=0;
    int iRinv=0;
    stringstream linestream(line);
    linestream >> id;
    seqMapIt=seqMap.find(id);
    if(seqMapIt != seqMap.end()){
      phenID.push_back(id);
      phenSeq.push_back(seqMapIt->second);
      if(!nrinv) rinverse.push_back(rinvtmp);
      int t=0;
      for(int col=1;col<nColumns;col++){
	switch(varType[col]){
	case VAR_DEP_VAR:
	  linestream >> val;
	  ytmp[t++]=val;
	  if(t==nTraits) y.push_back(ytmp);
	  break;
	case VAR_RINVERSE:
	  linestream >> val;
	  rinvtmp[iRinv++]=sqrt(val);
	  if(iRinv==nTraits)rinverse.push_back(rinvtmp);
	  break;
	case VAR_COVARIATE:
	  linestream >> val;
	  valMatrix[iCv++].push_back(val);
	  break;
	case VAR_CLASS:
	  linestream >> sVal;
	  classIt=classMatrix[iCl].find(sVal);
	  if(classIt==classMatrix[iCl].end()){
	    classValues[iCl].push_back(sVal);
	    classMatrix[iCl][sVal]=nLevels[iCl]++;
	    betaClass[iCl].push_back(Vector2d::Zero());
	  }
	  iVal=classMatrix[iCl][sVal];
	  posMatrix[iCl++].push_back(iVal);
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
      cout << setw(5) << phenSeq[i] << " ID " << setw(10) <<  ID[phenSeq[i]];
      for(int t=0;t<nTraits;t++){
	cout << " y " << setw(10) << setprecision(5) << y[i][t] ;
      }
      for(int iRinv=0;iRinv<nrinv;iRinv++) cout << " Rinv " << setw(10) << setprecision(5) << rinverse[i](iRinv) ;
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
    hmmFileName="HMMIM" + to_string(HMM.nStates) + "_" + to_string(HMM.nLoci) +".bin";
    if((failCode=HMM.read(hmmFileName))!=0) {
      cout << hmmFileName << " failed to be read  with code " << failCode << endl;
      string hmmFileNameAlt;
      hmmFileNameAlt="HMM" + to_string(HMM.nStates) + "_" + to_string(HMM.nLoci) +".bin";
      if((failCode=HMM.read(hmmFileNameAlt))==0) {
	cout << hmmFileNameAlt << " was found however." << endl;
      }
    }

    for(int k=0;k<nStates;k++) HMM.pi[k]=1./((double) nStates);
    HMM.initPComb();
  
  
  

  //Estimation
  int nComb=HMM.nComb;

  vector<double> locusf(nComb),locusE(nStates);
  vector< vector<double> > f(nLoci,locusf),b(nLoci,locusf);
  vector<vector<double> > E(nLoci,locusE),pState(nLoci,locusE),piVec(nLoci,locusf);
  vector<double> recombFrac(nLoci,0.),rf(nLoci,0.);
  vector<double> pVec(nComb);
  for(int iter=0;iter<nIter;iter++){
    cout << "Iteration " << iter <<endl;
    HMM.initBW(priorCount);
     for(long i=0;i<nLoci;i++){
      E[i].assign(nStates,0.);
      pState[i].assign(nStates,0.);
      piVec[i].assign(nComb,0.);
    }
    HMM.initSeq();

    if(0){
      forwardVec(0,nLoci,HMM,X[0],HMM.piComb,f);
      forward(0,nLoci,HMM,X[0],HMM.piComb);
      for(int l=0;l<nComb;l++){
	cout << "f " << HMM.stateI[l] << "/" << HMM.stateJ[l] ;
	for(long i=0;i<4;i++){
	  
	  cout << " " << fixed << setprecision(4) << f[i][l] << " : " << HMM.loci[i].f[l];
	}
	cout << endl;
      }
      backwardVec(0, nLoci,HMM,X[0],HMM.piComb,b);
      backward(0, nLoci,HMM,X[0],HMM.piComb);
      for(int l=0;l<nComb;l++){
	cout << "b " << HMM.stateI[l] << "/" << HMM.stateJ[l];
	for(long i=nLoci-1;i>nLoci-5;i--){
	  
	  cout << " " << fixed << setprecision(4) << b[i][l] << " : " << HMM.loci[i].b[l];
	}
	cout << endl;
      }

    }
#pragma omp parallel firstprivate(f,b,E,pState,piVec,rf)
    {
      #pragma omp for schedule(static)
      for(int seq=0;seq< X.size();seq++){
	
	forwardVec(0,nLoci,HMM,X[seq],HMM.piComb,f);
	backwardVec(0, nLoci,HMM,X[seq],HMM.piComb,b);
	
	vector<double> Pvec(nComb);
	
	//Compute Probabilities and estimates
	if(iter==(nIter-1)){
	  
	}
	for(long i=0;i<nLoci;i++){
	  double Psum=0;
	  double Pval;
	  for(int l=0;l<nComb;l++){
	    Psum+=f[i][l]*b[i][l];
	  }
	  for(int l=0;l<nComb;l++){
	    int I=HMM.stateI[l];
	    int J=HMM.stateJ[l];
	    Pval=f[i][l]*b[i][l]/Psum;
	    Pvec[l]=Pval;
	    piVec[i][l]+=Pval;
	    pState[i][I]+=Pval;
	    pState[i][J]+=Pval;
	    if(X[seq][i]<0){
	      E[i][I]+=Pval*HMM.loci[i].e[I];
	      E[i][J]+=Pval*HMM.loci[i].e[J];
	    }
	    if(X[seq][i]==2){
	      E[i][I]+=Pval;
	      E[i][J]+=Pval;
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
		E[i][I]+=Pval*wt1;
		E[i][J]+=Pval*wt2;
	      }
	    }
	  }
	}
	
      }//end seq
      #pragma omp critical(updateHMMdata)
      {
	for(long i=0;i<nLoci;i++){
	  for(int l=0;l<nComb;l++){
	    HMM.loci[i].piVec[l]+=piVec[i][l];
	    HMM.loci[i].f[l]+=f[i][l];
	    HMM.loci[i].b[l]+=b[i][l];
	  }
	  for(int I=0;I<nStates;I++){ 
	    HMM.loci[i].E[I]+=E[i][I];
	    HMM.loci[i].pState[I]+=pState[i][I];
	  }
	}
      }
      
    }//end parallel
    
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


  if(nIter) {
    snpResults.open(SNPName);
    snpResults << "Loci\tName\tChrom\tPos";
    if(nIter) for(int l=0;l<nStates;l++) snpResults << "\tframeTemplate" << l ;
    for(int l=0;l<nStates;l++) snpResults<< "\tframeFreq" << l ;
    snpResults << endl;
    
    for(long i=0;i<nLoci;i++){
      snpResults << i << "\t" << lociMap[i].name << "\t" << lociMap[i].chrom << "\t"<< lociMap[i].pos;
      for(int l=0;l<nStates;l++) snpResults << "\t" <<HMM.loci[i].e[l];
      if(nIter) for(int l=0;l<nStates;l++) snpResults << "\t" <<HMM.loci[i].pState[l]/(2.*(nSeq+priorCount));
      snpResults << endl;
    }
    snpResults.close();
  }
    // Add QTL loci;
    int q=0;
  int nQTLLoci=0;
  int nSNPLoci=nLoci;
  
  if(QTLMap){
    QTLMapFile.open(QTLMapName);
    if(QTLMapFile.fail()) {
      cout << "Failed to open QTL map file " + QTLMapName << endl;
      exit(101);
    }
    getline(QTLMapFile,line);
    while(getline(QTLMapFile,line)){
      stringstream linestr(line);
      linestr >> aMapLocus.name;
      linestr >> aMapLocus.chrom;
      linestr >> aMapLocus.pos;
      linestr >> QTLClassStr;
      
      aMapLocus.isSNP=-1;
      aMapLocus.isQTL=nQTLLoci++;
      if(QTLClasses.find(QTLClassStr)!=QTLClasses.end()){
	iVal=QTLClasses[QTLClassStr];
      }
      else{
	cout << "Couldn't find " + QTLClassStr + " in QTLClasses = " << endl;
	exit(312);
      }
      lociMap.push_back(aMapLocus);
      QTLClassVec.push_back(iVal);
      QTLClassCount[iVal]++;
    }
    cout << "Read in QTL map file" << endl << endl;;
  }
  else{
    nQTLClasses=1;
    QTLClasses[""]=0;
    QTLClassCount.push_back(0);
    for(int chr=0;chr<nChrom;chr++)  HMM.loci[chromStart[chr]].newChrom=0;
    for(int chr=0;chr<nChrom;chr++){
      aMapLocus.chrom=Chrom[chr];
      //cout << "qpos "<< chromStart[chr] << " "  << lociMap[chromStart[chr]].pos+1 << " " << lociMap[chromStart[chr+1]-1].pos << endl;
      for(long qPos=lociMap[chromStart[chr]].pos+1;qPos<lociMap[chromStart[chr+1]-1].pos;qPos+=freqQTL,q++){
	aMapLocus.name="QTL_" + to_string(Chrom[chr]) + "_" + to_string(qPos);
	aMapLocus.pos=qPos;
	aMapLocus.isSNP=-1;
	aMapLocus.isQTL=q;
	nQTLLoci=q+1;
	lociMap.push_back(aMapLocus);
	QTLClassVec.push_back(0);
	QTLClassCount[0]++;
	//if(q< 100) cout << q <<" " << aMapLocus.name << endl;
      }
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
    for(int i=chromStart[chr];i<chromStart[chr+1];i++) {
      HMM.loci[i].pos=lociMap[i].pos;
      HMM.loci[i].newChrom=0;
    }
    HMM.loci[chromStart[chr]].newChrom=1;
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
  //vector<vector<double> > XHaploT(nPheno);
  //for(int a=0;a<nPheno;a++) XHaploT[a].resize(nLoci);
  //haploMapLocus XHaploMapLocus(nStates,nPheno);
  //vector<haplMapLocus > XHaploMap(nQTLLoci,XHaploMapLocus);
  //vector<int> tmpState(2*nPheno);
  //vector<int> tmpNLoci(nStates);
  //vector<vector<int> > tmpLocus(nStates,tmpState)
 
  
  
  //
  // MCMC
  //

  
  MCMCSamples.open(MCMCName);


  //double sig2bPrior=.5,sig2ePrior=20;
  
  
  //int nSamples=40;
  //int nBurnIn=2;
  
  vector<double> AoverIVecP,AoverIVec;
  Matrix2d sig2e=sig2ePrior,sig2b=sig2bPrior,sig2g;
  Matrix2d scaleRes=sig2ePrior*(nusig2e-((double) nTraits)-1.);
  Matrix2d scaleB=sig2bPrior*(nusig2b-((double) nTraits)-1.);
  double nuTilde,uSmp;
  double flipDelta=1./((double) nStates);
  vector<qtlXLocus> qtlVec;
  vector<qtlXLocusSum> qtlSumVec;
  vector<int> windowVec,windowLocus(nTraits,0);
  vector<vector<int> > windowSumVec;
  vector<long> activeLoci;
  long activePos=0;

  //vector<double> piVec(nQTLClasses,pi);
  //vector<double> sig2bVec(nQTLClassses,sig2b);
  vector<int> nQTLVec(nQTLClasses);
  
  uSmp=u(gen);

  qtlVec.resize(nQTLLoci);
  qtlSumVec.resize(nQTLLoci);
  windowSumVec.assign(nQTLLoci,windowLocus);
  windowVec.resize(nQTLLoci);
  activeLoci.resize(0);
  LLT<Matrix2d> lltOfsig2b(sig2b);
  Matrix2d Bl=lltOfsig2b.matrixL();
  for(long i=0;i<nQTLLoci;i++) {
    qtlVec[i].init(nStates,nTraits,Bl,pi,rho);
    qtlSumVec[i].init(nStates,nTraits);
    if(qtlVec[i].active){
      activeLoci.push_back(i);
    }
  }
  long nQTL=activeLoci.size();

  vector<Vector2d> yDev(nPheno,Vector2d::Zero()),gHat(nPheno,Vector2d::Zero()),
    gHatSum(nPheno,Vector2d::Zero());
  vector<Matrix2d> gHatSumSq(nPheno,Matrix2d::Zero());
  vector<double> xVec(nPheno),probClass(nComb);
  Vector2d *yDevpt,*gHatpt;
  //MCMCSamples << "Sample\tmu\tsig2b\tsig2e\tsig2g\tpi";
  
  MCMCSamples << "Sample";
  for(int i=0;i<nTraits;i++) MCMCSamples << "\tmu_"+phenLabels[i+1];
  if(QTLMap){
    for(auto it=QTLClasses.begin();it !=QTLClasses.end();it++){
      int c=it->second;
      MCMCSamples <<"\t"<<  it->first << "_Pi" ;
      for(int i=0;i<nTraits;i++)for(int j=0;j<=i;j++)
				  MCMCSamples <<"\t"<<  it->first << "_sig2b["<< phenLabels[i+1] << "][" << phenLabels[j+1] << "]";
    }
  }
  else{
    MCMCSamples <<"\tPi" ;
    for(int i=0;i<nTraits;i++)for(int j=0;j<=i;j++)
				MCMCSamples <<"\tsig2b["<< phenLabels[i+1] << "][" << phenLabels[j+1] << "]";
  }
  
  for(int i=0;i<nTraits;i++)for(int j=0;j<=i;j++)
			      MCMCSamples << "\tsig2e["<<  phenLabels[i+1] << "][" << phenLabels[j+1] << "]";
  for(int i=0;i<nTraits;i++)for(int j=0;j<=i;j++)
			      MCMCSamples << "\tsig2g["<< phenLabels[i+1] << "][" << phenLabels[j+1]  << "]";
  MCMCSamples << "\trho";
  for(int iCl=0;iCl<nClass;iCl++) {
    if(classRandom[iCl] > -1){
       for(int i=0;i<nTraits;i++)for(int j=0;j<=i;j++)
			      MCMCSamples << "\tsig2_"+className[iCl] +"["<< phenLabels[i+1] << "][" << phenLabels[j+1] << "]";
    }
  }
  for(int iCl=0;iCl<nClass;iCl++) 
    for(int l=0;l<nLevels[iCl];l++)
      for(int t=0;t<nTraits;t++)
	MCMCSamples << "\t" << className[iCl] + "_" + classValues[iCl][l] +"_" + phenLabels[t+1];  
  
  for(int iCv=0;iCv<nCovariate;iCv++)
    for(int t=0;t<nTraits;t++) MCMCSamples << "\t" << covName[iCv] +"_" + phenLabels[t+1];
  MCMCSamples << endl;

  
 

  
  vector<double> logPHaplo(nPheno);

  vector<int> dVec(nStates,0),dVec1(nStates,0),deltaZero(nStates,0),deltaOne(nStates,1),dsum(nStates,0),dsum2(nStates,0),deltaBase;
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

  
  int nDeltaStates=deltaStates.size();
  uniform_int_distribution<int> uDelta(0,nDeltaStates+1);

  deltaStates.push_back(deltaOne);
  deltaStates.push_back(deltaZero);
  int computeLhsV=1;

  vector<Vector2d> rhsV;
  vector<Matrix2d> lhsV,lhsVs;
  rhsV.assign(nStates,Vector2d::Zero());
  lhsV.assign(nComb,Matrix2d::Zero());
  lhsVs.assign(nStates,Matrix2d::Zero());
  vector<Vector2d> rhsVtmp=rhsV;
  vector<Matrix2d> lhsVtmp=lhsV,lhsVstmp=lhsVs;
  vector<Vector2d> yDevtmp=y;
  vector<vector<Matrix2d> > lhsVArray(nQTLLoci,lhsV),lhsVsArray(nQTLLoci,lhsVs);
  Vector2d *rhsVpt;
  Matrix2d *lhsVpt,*lhsVspt;
  

  f.resize(nLoci,locusf);
  b.resize(nLoci,locusf);

  for(int s=0;s<nSamples;s++){
    int printCnt=0;
    vector<Matrix2d> Binv(nQTLClasses,Matrix2d::Zero());
    vector<double> BinvLogDet(nQTLClasses);
    double RinvLogDet;
    Matrix2d Rinv=sig2e.inverse();
    RinvLogDet=log(Rinv.determinant());
    for(int c=0;c<nQTLClasses;c++) {
      Binv[c]=sig2bVec[c].inverse();
      BinvLogDet[c]=log(Binv[c].determinant());
    }
    //
    // Fill in HaploVec ---------------------------------------------
    //
    if((s%FreqToSampleHaplo)==0){


      /*for(long i=0;i<nQTLLoci;i++) {
	for(int l=0;l<nStates;l++){
	  XHaploMap[i][l].resize(0);
	}
	}*/
      cout << endl << "Sampling haplotypes" << endl;
      computeLhsV=1;
      int nFlipped=0;

      vector<int> vecFlipped(nPheno,0);
      vector <int>   XHaploNew(nQTLLoci);
#pragma omp parallel for firstprivate(f,XHaploNew)
      for(int a=0;a<nPheno;a++){
	int seq=phenSeq[a];
	int ranClass,ranClassOld;
	double logPNewvsOld;
	//forward(0, nLoci,HMM,X[seq],HMM.piComb);
	forwardVec(0, nLoci,HMM,X[seq],HMM.piComb,f);
	vector<double> Pvec(nComb);
	{
	  double Psum=0;
	  double Pval;
	  for(int l=0;l<nComb;l++){
	    //Pval=HMM.loci[nLoci-1].f[l];
	    Pval=f[nLoci-1][l];
	    Psum+=Pval;
	    Pvec[l]=Pval;
	  }
	  for(int l=0;l<nComb;l++){
	    Pvec[l]/=Psum;
	  }
	  double uSmp=u(gen);
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
	    //Pval*=HMM.loci[i].f[l];
	    Pval*=f[i][l];
	    Pvec[l]=Pval;
	    Psum+=Pval;
	  }
	  for(int l=0;l<nComb;l++){
	    Pvec[l]/=Psum;
	  }
	  double uSmp=u(gen);
	  double Thresh=0;
	  
	  for(int l=0;uSmp>Thresh && l< nComb;l++){
	    ranClass=l;
	    Thresh+=Pvec[l];
	  }
	  int q=lociMap[i].isQTL;
	  if(q>=0)XHaploNew[q]=ranClass;
	  logPNewvsOld+=log(Pvec[ranClass]);
	}
	
	Vector2d yDevNew=y[a]-mu;
	for(int iCl=0;iCl<nClass;iCl++){
	  yDevNew-=betaClass[iCl][posMatrix[iCl][a]];
	}
	for(int iCv=0;iCv<nCovariate;iCv++) yDevNew-=valMatrix[iCv][a]*betaCov[iCv]; 
	
	for(long i=0;i<nQTL;i++){
	  long locus=activeLoci[i];
	  int l=XHaploNew[locus];
	  int I=HMM.stateI[l]; 
	  int J=HMM.stateJ[l];
	  for(int t=0;t<nTraits;t++){
	    if(qtlVec[locus].activeTrait[t])yDevNew(t)-=(qtlVec[locus].delta[I]+qtlVec[locus].delta[J])*qtlVec[locus].b(t);
	  }
	}
	uSmp=u(gen);
	yDevNew=yDevNew.cwiseProduct(rinverse[a]);
	yDev[a]=yDev[a].cwiseProduct(rinverse[a]);
	if(s==0 || log(uSmp)<-.5*(yDevNew.dot(Rinv*yDevNew)-yDev[a].dot(Rinv*yDev[a]))){
	  for(long i=0;i<nQTLLoci;i++) {
	    XHaplo[i][a]=XHaploNew[i];
	    //XHaploT[a][i]=XHaploNew[i];
	  }
	  logPHaplo[a]=logPNewvsOld;
	  vecFlipped[a]=1;
	}
      }
      nFlipped=0;
      for(int a=0;a<nPheno;a++) {
	if(vecFlipped[a]==1) nFlipped++;
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



    //double ssb=scaleB;
    Matrix2d sse=scaleRes,ssg=Matrix2d::Zero();
    vector<Matrix2d> ssbVec(nQTLClasses);
    for(int c=0;c<nQTLClasses;c++){
      ssbVec[c]=sig2bPriorVec[c]*(nusig2bVec[c]-((double) nTraits)-1.);
    }
    
    nQTL=activeLoci.size();
    yDev=y;
    gHat.assign(nPheno,Vector2d::Zero());
#pragma omp parallel for schedule(static) 
    for(int a=0;a<nPheno;a++){
      Vector2d ysum=Vector2d::Zero();
      int seq=phenSeq[a];
      //#pragma omp parallel for reduction (+:ysum)
      for(long i=0;i<nQTL;i++){
	long locus=activeLoci[i];
	int seqClass=XHaplo[locus][a];
	int I=HMM.stateI[seqClass];
	int J=HMM.stateJ[seqClass];
	for(int t=0;t<nTraits;t++) {
	  if(qtlVec[locus].activeTrait[t])ysum(t)-=(qtlVec[locus].delta[I]+qtlVec[locus].delta[J])*qtlVec[locus].b(t);
	}
      }

      ysum-=mu;
      for(int iCl=0;iCl<nClass;iCl++){
	ysum-=betaClass[iCl][posMatrix[iCl][a]];
      }
      for(int iCv=0;iCv<nCovariate;iCv++) ysum-=valMatrix[iCv][a]*betaCov[iCv]; 
      

   
      yDev[a]+=ysum;
    }
    

    activeLoci.resize(0);
    int nRejectI2A=0;
    int nRejectA2I=0;
    nQTLVec.assign(nQTLClasses,0);
    for(long i=0;i<nQTLLoci;i++){

      
      if(deltaSampler==4 || deltaSampler==3){
	if(i<(nQTLLoci-1)){
	  if(qtlVec[i].active || qtlVec[i+1].active){
	    int a0=qtlVec[i].active;
	    Vector2d b=qtlVec[i].b;
	    dVec=qtlVec[i].delta;
	    int a1=qtlVec[i+1].active;
	    Vector2d b1=qtlVec[i+1].b;
	    dVec1=qtlVec[i+1].delta;
	    double ssSw;
	    for(int a=0;a<nPheno;a++){
	      	int seqClass=XHaplo[i][a];
	      	int seq1Class=XHaplo[i+1][a];
		if(seqClass != seq1Class){
		  int I=HMM.stateI[seqClass];
		  int J=HMM.stateJ[seqClass];
		  int I1=HMM.stateI[seq1Class];
		  int J1=HMM.stateJ[seq1Class];
		  Vector2d yd=yDev[a];
		  if(a0){
		    for(int t=0;t<nTraits;t++){
		      if(qtlVec[i].activeTrait[t])yd(t)-=((dVec[I1]+dVec[I1])-(dVec[I]+dVec[I]))*b(t);
		    }
		  }
		  if(a1){
		    for(int t=0;t<nTraits;t++){
		      if(qtlVec[i+1].activeTrait[t])yd(t)-=((dVec1[I]+dVec1[I])-(dVec1[I1]+dVec1[I1]))*b1(t);
		    }
		  }
		  yd=yd.cwiseProduct(rinverse[a]);
		  ssSw+=yd.dot(Rinv*yd);
		  yd=yDev[a].cwiseProduct(rinverse[a]);
		  ssSw-=yd.dot(Rinv*yd);
		}
	    }
	    double probSw=exp(-.5*ssSw);
	    if(u(gen)*(1.+probSw)<probSw){
	      for(int a=0;a<nPheno;a++){
	      	int seqClass=XHaplo[i][a];
	      	int seq1Class=XHaplo[i+1][a];
		if(seqClass != seq1Class){
		  int I=HMM.stateI[seqClass];
		  int J=HMM.stateJ[seqClass];
		  int I1=HMM.stateI[seq1Class];
		  int J1=HMM.stateJ[seq1Class];
		  if(a0){
		    for(int t=0;t<nTraits;t++){
		      if(qtlVec[i].activeTrait[t])yDev[a](t)-=((dVec[I1]+dVec[I1])-(dVec[I]+dVec[I]))*b(t);
		    }
		  }
		  if(a1){
		    for(int t=0;t<nTraits;t++){
		      if(qtlVec[i+1].activeTrait[t])yDev[a](t)-=((dVec1[I]+dVec1[I])-(dVec1[I1]+dVec1[I1]))*b1(t);
		    }
		  }
		}
	      }
	      vector<int>  activeTrait=qtlVec[i].activeTrait;
	      int T=qtlVec[i].t;
	      qtlVec[i].active=a1;
	      qtlVec[i].b=b1;
	      qtlVec[i].activeTrait=qtlVec[i+1].activeTrait;
	      qtlVec[i].t=qtlVec[i+1].t;
	      qtlVec[i].delta=dVec1;
	      qtlVec[i+1].active=a0;
	      qtlVec[i+1].b=b;
	      qtlVec[i+1].activeTrait=activeTrait;
	      qtlVec[i+1].delta=dVec;
	      qtlVec[i+1].t=T;
	    }
	    
	  }
	}
      }

      
      int active=qtlVec[i].active;
      int nowActive=0;
      int proposeActive=0;
      
      //double logAoverI=log((1.-pi)/pi);
      Vector2d b=qtlVec[i].b;
      
      rhsV.assign(nStates,Vector2d::Zero());
      if(computeLhsV){
	lhsV.assign(nComb,Matrix2d::Zero());
	lhsVs.assign(nStates,Matrix2d::Zero());
      }
      else{
	lhsV=lhsVArray[i];
	lhsVs=lhsVsArray[i];
      }
      uSmp=u(gen);
      
      if(!active && uSmp>inactiveProposal) proposeActive=1;
      computeLhsV=1;
      if(proposeActive || active /*|| computeLhsV*/){
	dVec=qtlVec[i].delta;
	
	rhsVtmp.assign(nStates,Vector2d::Zero());
	lhsVtmp.assign(nComb,Matrix2d::Zero());
	lhsVstmp.assign(nStates,Matrix2d::Zero());
	//yDevtmp.assign(nPheno,0.0);
	if(active){
	  for(int t=0;t<nTraits;t++){
	    if(qtlVec[i].activeTrait[t]){
	      for(int a=0;a<nPheno;a++){
		int seqClass=XHaplo[i][a];
		int I=HMM.stateI[seqClass];
		int J=HMM.stateJ[seqClass];
		yDev[a](t)+=(dVec[I]+dVec[J])*b(t);
	      }
	    }
	  }
	}

	//#pragma omp parallel firstprivate(rhsVtmp,lhsVtmp,lhsVstmp)
	{
	  rhsVpt=rhsVtmp.data();
	  lhsVpt=lhsVtmp.data();
	  lhsVspt=lhsVstmp.data();
	  //#pragma omp for schedule(static)
	  
	  Matrix2d RinvScaled=Matrix2d::Zero();
	  for(int a=0;a<nPheno;a++){
	    Vector2d ydev(nTraits);
	    //ydev=(Rinv*yDev[a].cwiseProduct(rinverse[a])).cwiseProduct(rinverse[a]);
	    for(int ti=0;ti<nTraits;ti++){
	      for(int tj=0;tj<nTraits;tj++){
		RinvScaled(ti,tj)=rinverse[a][ti]*Rinv(ti,tj)*rinverse[a][tj];
	      }
	    }
	    ydev=RinvScaled*yDev[a];
	    int seqClass=XHaplo[i][a];
	    int I=HMM.stateI[seqClass];
	    int J=HMM.stateJ[seqClass];
	    *(rhsVpt+I)+=ydev;
	    *(rhsVpt+J)+=ydev;
	    if(computeLhsV){
	      *(lhsVspt+I)+=RinvScaled;
	      *(lhsVspt+J)+=RinvScaled;
	      *(lhsVpt+seqClass)+=2.0*RinvScaled;
	    }
	  }

	  //#pragma omp critical(updateRhsLhs)
	  {
	    for(int I=0;I<nStates;I++){
	      rhsV[I]+=rhsVtmp[I];
	      // cout << "rhsV " << I << "=" << rhsV[I] <<endl;
	    } 
	    if(computeLhsV){
	      for(int I=0;I<nStates;I++){
		lhsVs[I]+=lhsVstmp[I];
		
		//  cout << "lhsVs " << I << "=" << lhsVs[I] <<endl;
	      }
	      for(int l=0;l<nComb;l++){
		lhsV[l]+=lhsVtmp[l];
		
		// cout << "lhsV " << l << "=" << lhsV[l] <<endl;
	      }
	    }
	  }
	}
	  // matches omp parallel
	  if(computeLhsV){
	    lhsVArray[i]=lhsV;
	    lhsVsArray[i]=lhsVs;
	  }
	}


      if(active || proposeActive){
	
	double psum=0,pInactive=1,ptot=0;
	double maxAoverI=0.;
	double psumP=0,pInactiveP=1,ptotP=0;
	double maxAoverIP=0.;
	int SD=0,T=0,SDT=0;
	if(!active) qtlVec[i].delta=deltaZero;
	int QTLClass=QTLClassVec[i];
	switch(deltaSampler){
	case 3:
	case 1:
	  
	  calcDeltaProposalFullX(sig2e,sig2bVec[QTLClass],Rinv,Binv[QTLClass],qtlVec[i].delta,deltaStates,
				dVec,rhsV,lhsV,lhsVs,piClassVec[QTLClass],rho,nDeltaStates,nStates,
				pInactive,maxAoverI,AoverIVec,psum);
	  
	  ptot=psum;
	  if(active)ptot+=pInactive;
	  if(u(gen) < psum/ptot) proposeActive=1;
	  deltaBase=deltaZero;
	  if(proposeActive){
	    //Propose Delta
	    uSmp=u(gen);
	    uSmp*=psum;
	    int sdt,sd,T;
	    for(sdt=0; uSmp>AoverIVec[sdt]   && sdt<nDeltaStates*(nTraits+1) ;sdt++){
	      uSmp-=AoverIVec[sdt];
	    }
	    if(sdt< nDeltaStates*(nTraits+1)){
	      SD=sdt/(nTraits+1);
	      T=sdt%(nTraits+1);
	      deltaBase=deltaStates[sd];
	    }
	    else{
	      cout << endl << "Full Delta Sampler Failed." << endl;
	      exit(314);
	    }
	    SD=sd;
	  }
	  psumP=psum;
	  pInactiveP=pInactive;
	  maxAoverIP=maxAoverI;
	  ptotP=psumP;
	  if(proposeActive)ptotP+=pInactiveP;
	  break;
	case 4:
	case 2:
	  
	  
	  if(!active) {
	    qtlVec[i].delta=deltaStates[uDelta(gen)];
	  }
	  if(0){  
	    cout << QTLClass << endl;
	    cout << sig2bVec[QTLClass] << endl;
	    cout << piClassVec[QTLClass] << endl;
	  }
	  
	  //cout << QTLClass << " " << Binv[0](0,0) << endl;
	  calcDeltaProposalX(sig2e,sig2bVec[QTLClass],Rinv,Binv[QTLClass],qtlVec[i].delta,
			    dVec,rhsV,lhsV,lhsVs,piClassVec[QTLClass],rho,nDeltaStates,nStates,
			    pInactive,maxAoverI,AoverIVec,psum,RinvLogDet,BinvLogDet[QTLClass]);
	  
	 
	  ptot=psum;
	  if(active)ptot+=pInactive;
	  
	  
	  if(u(gen) < psum/ptot) proposeActive=1;
	  
	  deltaBase=deltaZero;
	  if(proposeActive){
	    //Propose Delta
	    uSmp=u(gen);
	    uSmp*=psum;
	    int sdt,sd,t;
	    for(sdt=0; uSmp>AoverIVec[sdt]   && sdt<(nStates+1)*(nTraits+1) ;sdt++){
	      uSmp-=AoverIVec[sdt];
	    }
	    deltaBase=qtlVec[i].delta;
	    
	    sd=sdt/(nTraits+1);
	    t=sdt%(nTraits+1);
	    if(sd != nStates){
	      deltaBase[sd]=1-qtlVec[i].delta[sd];
	    }
	    SD=sd;
	    T=t;
	  }
	  
	  //cout << QTLClass << " " << Binv[0](0,0) << endl;
	  calcDeltaProposalX(sig2e,sig2bVec[QTLClass],Rinv,Binv[QTLClass],deltaBase,
			     dVec,rhsV,lhsV,lhsVs,piClassVec[QTLClass],rho,nDeltaStates,nStates,
			     pInactiveP,maxAoverIP,AoverIVecP,psumP,RinvLogDet,BinvLogDet[QTLClass]);
	  ptotP=psumP;
	  if(proposeActive)ptotP+=pInactiveP;
	  break;
	default:
	  cout << endl << "Invalid deltaSampler=" << deltaSampler << "specified." <<endl;
	  exit(313);
	}
	

	
	
	int accept=0;
	double inactiveMult=1.,inactiveMultP=1.;
	if(!active ) inactiveMult=1.-inactiveProposal;
	if(!proposeActive) inactiveMultP=1.-inactiveProposal;
	if(u(gen)*ptotP*inactiveMult < inactiveMultP*ptot*exp(maxAoverI-maxAoverIP)) accept=1;
	
	
	
	if(!accept) {
	  if(active) nRejectA2I++;
	  else {
	    nRejectI2A++;
	    qtlVec[i].delta=deltaZero;
	    qtlVec[i].activeTrait.assign(nTraits,0);
	    qtlVec[i].b.setZero();
	  }
	}
	if(accept || active){	
	  if(accept && !proposeActive){
	    qtlVec[i].active=0;
	    qtlVec[i].b.setZero();
	    qtlVec[i].delta=deltaZero;
	    qtlVec[i].activeTrait.assign(nTraits,0);
	  }
	  else{ //
	    qtlVec[i].active=1;
	    activeLoci.push_back(i);
	    nowActive=1;
	    if(accept) {
	      qtlVec[i].delta=deltaBase;
	      if(T< nTraits){
		qtlVec[i].activeTrait.assign(nTraits,0);
		qtlVec[i].activeTrait[T]=1;
	      }
	      else{
		qtlVec[i].activeTrait.assign(nTraits,1);
	      }
	      qtlVec[i].t=T;
	    }
	    dVec=qtlVec[i].delta;
	    T=qtlVec[i].t;
	    Vector2d RHS=Vector2d::Zero();
	    Matrix2d LHS=Matrix2d::Zero();
	    double sumXY,sumXX;
	    sumXY=0;
	    sumXX=0;
	    int IJ=0;
	    for(int I=0;I<nStates;I++){
	      if(dVec[I]) {
		if(T< nTraits){
		  LHS(T,T)+=lhsVs[I](T,T);
		  RHS(T)+=rhsV[I](T);
		}
		else{
		  LHS+=lhsVs[I];
		  RHS+=rhsV[I];
		}
	      }
	      for(int J=0;J<=I;J++,IJ++){
		if(dVec[I]&& dVec[J]) {
		  if(T< nTraits){
		    LHS(T,T)+=lhsV[IJ](T,T);
		  }
		  else{
		    LHS+=lhsV[IJ];
		  }
		}
	      }
	    }
	   
	    Matrix2d LHSinv=(LHS+Binv[QTLClass]).inverse();
	    
	    LLT<Matrix2d> llt(LHSinv);
	    Matrix2d L=llt.matrixL(); 
	    //b=LHSinv*RHS;
	    //cout << "b " << LHS(0,0) << " " << RHS(0) << " " << b(0) << endl;
	    
	    Zvec(zVec);
	    b=LHSinv*RHS+L*zVec;
	    /*
	    for(int tj=0;tj<nTraits;tj++) {
	      val=Z(gen);
	      for(int ti=tj;ti<nTraits;ti++){
		b[ti]+=L(ti,tj)*val;
	      }
	    }
	    */
	    
	  
	    //sample b
	    
	    qtlVec[i].b=b;
	    ssbVec[QTLClass]+=b*b.transpose();
	    
	    //if(printCnt++<5) cout<<  "b["<<i<<"] " <<  b << "\n LHS \n" << LHS << "\n Binv \n" << Binv[QTLClass] << "\n L \n" << L << "\n Z \n" << zVec << "\n RHS \n" << RHS << "\n ssbVec\n"<<  ssbVec[QTLClass] <<endl;
	    nQTLVec[QTLClass]++;
	    yDevpt=yDev.data();
            gHatpt=gHat.data();
	    //#pragma omp parallel for schedule(static)
	    for(int t=0;t<nTraits;t++){
	      if(T==t || T==nTraits){
		for(int a=0;a<nPheno;a++){
		  //int seq=phenSeq[a];
		  int seqClass=XHaplo[i][a];
		  int I=HMM.stateI[seqClass];
		  int J=HMM.stateJ[seqClass];
		  yDev[a](t)-=(dVec[I]+dVec[J])*b(t);
		  gHat[a](t)+=(dVec[I]+dVec[J])*b(t);
		}
	      }
	    }
	    
	}
	}
    }
      
      //if(qtlVec[i].active && !nowActive) cout << "Switch locus "<< i << " " << qtlVec[i].active << "=>" << nowActive << endl;
      
      
      
      if(s>=nBurnIn)qtlSumVec[i].updateSum(qtlVec[i],nTraits);
      
    }

    //Calc sig2g
    ssg.setZero();
    Vector2d sumg(Vector2d::Zero());
    //#pragma omp parallel for schedule(static) reduction(+:ssg,sumg)
    for(int a=0;a<nPheno;a++) {
      ssg+=gHat[a]*gHat[a].transpose();
      sumg+=gHat[a];
    }
    sumg/=(double) nPheno;
    sig2g=ssg/((double) nPheno)-(sumg*sumg.transpose());
    
    //update gHatSum
    if(s>=nBurnIn){
      //#pragma omp parallel for schedule(static)
      for(int a=0;a<nPheno;a++) {
	gHat[a]-=sumg;
	gHatSum[a]+=gHat[a];
	gHatSumSq[a]+=gHat[a]*gHat[a].transpose();
      }
    }

    
 	
    // Update pi
    for(int c=0;c<nQTLClasses;c++){
      if(piPriorCountVec[c] ){
	gamma_distribution<double> gammaA(piPriorVec[c]*piPriorCountVec[c]+((double)(nQTLLoci-nQTLVec[c])),1.0);
	gamma_distribution<double> gammaB((1.-piPriorVec[c])*piPriorCountVec[c]+((double) nQTLVec[c]),1.0);
	double numPi=gammaA(gen);
	double denPi=numPi+gammaB(gen);
	piClassVec[c]=numPi/denPi;
      }
    }
    
    
    computeLhsV=0;
    //update mu;
    Vector2d sumXY=Vector2d::Zero();
    Matrix2d sumXX=Rinv*0.0001;
    //#pragma omp parallel for reduction (+:sumXY,sumXX) schedule(static)
    for(int a=0;a<nPheno;a++){
      yDev[a]+=mu;
      for(int ti=0;ti<nTraits;ti++){
	for(int tj=0;tj<nTraits;tj++){
	  sumXY(ti)+=rinverse[a](ti)*Rinv(ti,tj)*rinverse[a](tj)*yDev[a](tj);
	  sumXX(ti,tj)+=rinverse[a](ti)*Rinv(ti,tj)*rinverse[a](tj);
	}
      }
    }
    Matrix2d sumXXinv=sumXX.inverse();
    LLT<Matrix2d> L(sumXXinv);
    Zvec(zVec);
    Vector2d deltaMu=sumXXinv*sumXY+L.matrixL()*zVec;


    yDevpt=yDev.data();
    //#pragma omp parallel for schedule(static)
    for(int a=0;a<nPheno;a++) *(yDevpt+a)-=deltaMu;
    // cout << "mu " << mu << " " << deltaMu<< " " << sumXY << " " << sumXX << " " << zVec  << endl;
    mu=deltaMu;
    
    //update Class effects;
    
    Matrix2d Dinv=Matrix2d::Zero();
    for(int iCl=0;iCl<nClass;iCl++){
      if(classRandom[iCl]>-1) {
	Dinv=sig2r[classRandom[iCl]].inverse();
      }
      else{
	Dinv=Rinv*0.0001;
      }
      vector<Matrix2d> sumXX(nLevels[iCl],Dinv);
      vector<Vector2d> sumXY(nLevels[iCl],Vector2d::Zero());
      //cout << "iCl " << iCl << " lambdaC " << lambdaC<< endl;` 
      //vector<Vector2d> deltaBeta(nLevels[iCl]);
      for(int a=0;a<nPheno;a++){
	int l=posMatrix[iCl][a];
	yDev[a]+=betaClass[iCl][l];
	for(int ti=0;ti<nTraits;ti++){
	  for(int tj=0;tj<nTraits;tj++){
	    sumXY[l](ti)+=rinverse[a](ti)*Rinv(ti,tj)*rinverse[a](tj)*yDev[a](tj);
	    sumXX[l](ti,tj)+=rinverse[a](ti)*Rinv(ti,tj)*rinverse[a](tj);
	  }
	}
	
      }
      for(int l=0;l<nLevels[iCl];l++){
	Matrix2d sumXXinv=sumXX[l].inverse();
	LLT<Matrix2d> L(sumXXinv);
	Zvec(zVec);
	betaClass[iCl][l]=sumXXinv*sumXY[l]+L.matrixL()*zVec;
      }
      //#pragma omp parallel for schedule(static)
      for(int a=0;a<nPheno;a++)  yDev[a]-=betaClass[iCl][posMatrix[iCl][a]];
    }
    //update Covariate effects;
    for(int iCv=0;iCv<nCovariate;iCv++){
      Vector2d sumXY=Vector2d::Zero();
      Matrix2d sumXX=Matrix2d::Zero();
      //#pragma omp parallel for schedule(static) reduction(+:sumXX,sumXY)
      for(int a=0;a<nPheno;a++){
	yDev[a]+=valMatrix[iCv][a]*betaCov[iCv];
	for(int ti=0;ti<nTraits;ti++){
	  for(int tj=0;tj<nTraits;tj++){
	    sumXY(ti)+=valMatrix[iCv][a]*rinverse[a](ti)*Rinv(ti,tj)*rinverse[a](tj)*yDev[a](tj);
	    sumXX(ti,tj)+=valMatrix[iCv][a]*rinverse[a](ti)*Rinv(ti,tj)*rinverse[a](tj)*valMatrix[iCv][a];
	  }
	}
      }
      
      //#pragma omp paralllel for schedule(static)
      Matrix2d sumXXinv=sumXX.inverse(); 
      LLT<Matrix2d> L(sumXXinv);
      Zvec(zVec);
      Vector2d deltaBeta=sumXXinv*sumXY+L.matrixL()*zVec;
      for(int a=0;a<nPheno;a++) yDev[a]-=deltaBeta*valMatrix[iCv][a];
      betaCov[iCv]=deltaBeta;
    }

  
    // update sig2r;
    for(int iCl=0;iCl<nClass;iCl++){
      int iR=classRandom[iCl];
      if(iR>-1){
	nuTilde=((double)betaClass[iCl].size())+nusig2r[iR];
	Matrix2d ssr=sig2rPrior[iR]*(nusig2r[iR]-1.-((double) nTraits));
	for(int l=0;l<nLevels[iCl];l++)  {
	  ssr+=betaClass[iCl][l]*betaClass[iCl][l].transpose();
	}
	Matrix2d wMat;
	rWishartX(ssr.inverse(),nuTilde,wMat);
	sig2r[iR]=wMat.inverse();
	//	cout << "ssr " << ssr << " iR " <<  iR << " Scale " << sig2rPrior[iR]*(nusig2r[iR]-2.) << " nuTilde " << nuTilde << endl;
      }
    }

    //update sig2e;
    nuTilde=((double) nPheno)+nusig2e;
    //cout << "sse " << sse << endl; 
    Matrix2d ssePart=Matrix2d::Zero();
    //#pragma omp parallel for schedule(static) reduction(+:ssePart)
    double ssMax=0,ydMax,rinvMax;
    int aMax;
    Matrix2d ss;
    Vector2d yd;
    for(int a=0;a<nPheno;a++)  {
      yd=yDev[a].cwiseProduct(rinverse[a]);
      ss=yd*yd.transpose();
      ssePart+=ss;
      if(ss(0,0) > ssMax){
	ssMax=ss(0,0);
	aMax=a;
      }
      // if(a<2 || (a<1150 and a>1147)) cout << a <<" sse Calc " << yDev[a](0) << " * " << rinverse[a](0) << " = " << yd << " ss=" << ss << " ssePart=" << ssePart(0,0) << endl;
    }
    //cout << "aMax " << aMax << " " << yDev[aMax](0) << " " << rinverse[aMax](0) << " " << ssMax << endl;
    sse+=ssePart;
    rWishartX(sse.inverse(),nuTilde,Rinv);
    sig2e=Rinv.inverse(); 
    //cout << "Sig2E " << sse(0,0) << " " << nuTilde << " " << Rinv(0,0) << " " << sig2e(0,0) << endl; 

    //update sig2b;
    for(int c=0;c<nQTLClasses;c++){
      
      nuTilde=((double) nQTLVec[c])+nusig2bVec[c];
      Matrix2d wMat;
      rWishartX(ssbVec[c].inverse(),nuTilde,wMat);
      sig2bVec[c]=wMat.inverse();
      // cout << c << " nuTidle=" << nuTilde << "=" << nQTLVec[c] << "+" << nusig2bVec[c] << "\n ssbVec=" << ssbVec[c] << "\n sig2b=" << sig2bVec[c] << endl;
    }
    if(s>=nBurnIn){
      for(int t=0;t<nTraits;t++){
	windowVec.assign(nQTLLoci,0);
	for(int i=0;i<nQTLLoci;i++){
	  if(qtlVec[i].active && qtlVec[i].activeTrait[t]){
	    int start=lociMap[i].start;
	    int end=lociMap[i].stop;
	    for(int j=start;j<=end;j++) windowVec[j]=1;
	  }
	}
	for(int i=0;i<nQTLLoci;i++) {
	  if(windowVec[i]) windowSumVec[i][t]++;
	}
      }
    }


    if((s % printFreq)==0){
      cout << endl;
      cout << "Sample: " << s << endl;
      cout << " Rho: " << rho <<endl;
      if(QTLMap){
	for(auto it=QTLClasses.begin();it !=QTLClasses.end();it++){
	  int c=it->second;
	  cout << setw(10) << it->first << " ";
	  cout << it->first << " " << setprecision(8) << " Pi: " << piClassVec[c] <<"\n Sig2b: \n" << sig2bVec[c] << endl;
	}
	cout << setprecision(8) << "\n Sig2e: \n" << sig2e << "\n Sig2g: \n" << sig2g;
      }
      else{
	cout << setprecision(8) << " Pi: " << piClassVec[0] <<"\n Sig2b: \n" << sig2bVec[0] << "\n Sig2e: \n" << sig2e << "\n Sig2g: \n" << sig2g;
      }
      for(int iCl=0;iCl<nClass;iCl++) {
	if(classRandom[iCl] > -1){
	  cout << "\n Sig2_"+className[iCl]+": \n" << sig2r[classRandom[iCl]];
	}
      }
      cout << endl;
      cout << " mu \n" << mu << endl;
      cout << "number Active: " <<  nQTL << " Number Rejected A2I: " << nRejectA2I<< " Number Rejected I2A: " << nRejectI2A << endl;
      int nQTLPrint=100;
      if(nQTLPrint>nQTL) nQTLPrint=nQTL;
      for(int i=0; i< nQTLPrint  ; i++ ) {
	cout  << " " << setw(6) << activeLoci[i];
	if((i %10)==9) cout <<endl;
      }
      cout << endl;
      if(nQTLPrint>nQTL-100)nQTLPrint=nQTL-100;
      for(int i=0; i< nQTLPrint  ; i++ ) {
	cout  << " " << setw(6) << activeLoci[nQTL-nQTLPrint+i];
	if((i %10)==9) cout <<endl;
      }
      
      cout << endl <<endl;
      if(s%outputFreq==0) {
	MCMCSamples << s ;
	printVector(MCMCSamples,mu+sumg,"\t") ;
	if(QTLMap){
	  for(auto it=QTLClasses.begin();it !=QTLClasses.end();it++){
	    int c=it->second;
	    MCMCSamples << "\t" << piClassVec[c]  ;
	    printLower(MCMCSamples,sig2bVec[c],"\t");
	  }
	  cout << "\t" <<sig2e <<"\t" << sig2g;
	}
	else{
	  MCMCSamples << "\t" << piClassVec[0]  ;
	  printLower(MCMCSamples,sig2bVec[0],"\t");	  
	}

	
	printLower(MCMCSamples,sig2e,"\t");
	printLower(MCMCSamples,sig2g,"\t");
	MCMCSamples << "\t" << rho;
	for(int iCl=0;iCl<nClass;iCl++) {
	  if(classRandom[iCl] > -1){
	    printLower(MCMCSamples,sig2r[classRandom[iCl]],"\t");
	  }
	}
	for(int iCl=0;iCl<nClass;iCl++) {
	  for(int l=0;l<nLevels[iCl];l++)
	    for(int t=0;t<nTraits;t++)
	      MCMCSamples << "\t" << betaClass[iCl][l][t];  
	}
	for(int iCv=0;iCv<nCovariate;iCv++) printVector(MCMCSamples,betaCov[iCv],"\t");
	MCMCSamples << endl; 
      }
    }
  }
  
  QTLResults.open(QTLName);
  QTLResults << "Loci\tName\tChrom\tPos\tmodelFreq\tAll";
  for(int t=0;t<nTraits;t++) {
    QTLResults << "\twindowFreq" << t;
    QTLResults << "\tb" << t;
    QTLResults << "\tActive" << t ;
    for(int l=0;l<nStates;l++) QTLResults << "\tDelta" << t << "_" << l ;
  }
  //for(int l=0;l<nStates;l++) QTLResults << "\thaploFreq" << l ;
  //for(int l=0;l<nStates;l++) QTLResults << "\thaploTemplate" << l ;
  QTLResults << endl;
  
  double numSampled;
  numSampled=(double)(nSamples-nBurnIn);
  for(long i=0;i<nLoci;i++){
    int q=lociMap[i].isQTL;
    if(q>=0){
      //cout << i << " " << q << endl;
      double numActive=1;
      if(qtlSumVec[q].active) numActive=(double) qtlSumVec[q].active;
      
      QTLResults << q << "\t" << lociMap[i].name << "\t" << lociMap[i].chrom << "\t"<< lociMap[i].pos<<"\t"<< ((double) qtlSumVec[q].active)/numSampled<<"\t"<< ((double) qtlSumVec[q].all)/numSampled  ;
      for(int t=0;t<nTraits;t++)  {
	QTLResults << "\t" << ((double) windowSumVec[q][t])/numSampled;
	QTLResults << "\t" << qtlSumVec[q].b[t]/numSampled << "\t" <<  qtlSumVec[q].activeTrait[t]/numActive;
	for(int l=0;l<nStates;l++) QTLResults << "\t" <<  2.*((double) qtlSumVec[q].delta[t][l])/((double) qtlSumVec[q].activeTrait[t])-1.;
      }
      //for(int l=0;l<nStates;l++) QTLResults << "\t" <<HMM.loci[i].pState[l]/(2.*(nSeq+priorCount));
      //for(int l=0;l<nStates;l++) QTLResults << "\t" <<HMM.loci[i].e[l];
      QTLResults << endl;
    }
  }
  
  gHatResults.open(gHatName);
  gHatResults << "ID";
  for(int t=0;t<nTraits;t++) gHatResults << "\tgHat_" + phenLabels[t+1] +"\tPEV_" + phenLabels[t+1];
  gHatResults << endl;
  for(int a=0;a<nPheno;a++){
    gHatResults << phenID[a];
    for(int t=0;t<nTraits;t++) gHatResults << "\t" << gHatSum[a](t)/numSampled << "\t" <<  (gHatSumSq[a](t)-gHatSum[a](t)*gHatSum[a](t)/numSampled)/numSampled << endl;
  }
}

