#include "impute.h"
int main(int argc,char **argv){
  
  string genoName,mapName;
  int nStates=4;
  double lambdaKb=500.;
  string regionString,line;
  int nRegions=0;
  string id;
  int failCode;
  vector<string> ID;
  idmap regions,seqMap;
  
  idmap::iterator seqMapIt,phenMapIt;
  
  vector<long> regionStart,regionEnd;
  vector<double> regionStartMb,regionEndMb;
  vector<int> regionChrom;
  //vector<string> regions;
  string baseName,QTLName,regionName;

  vector<vector<int> > X;


  Configuration config;
  int iVal;
  double dVal,val;
  string sVal;
 if(argc>1){
    if(!config.Load(argv[1])) exit(1);
    config.Get("genoName",genoName);
    config.Get("mapName",mapName);
    config.Get("nStates",nStates);
    config.Get("lambdaKb",lambdaKb);
    config.Get("baseName",baseName);
    QTLName=baseName+"_QTLResults.txt";
    regionName=baseName+"_regionResults.txt";
    config.Get("QTLName",QTLName);
    //config.Get("gHatName",gHatName);
    config.Get("regionName",regionName);
    if(config.Get("regions",regionString)){
      int next;
      
      while(regionString.length() && (next=regionString.find_first_not_of(" \t")) !=string::npos){
	double val;
	regionString=regionString.substr(next);
	string label=regionString,varLabel;
	int last=regionString.find_first_of(" \t");
	if(last != string::npos) label=regionString.substr(0,last);
	regions[label]=nRegions;
	if(!config.Get("chrom_"+label,iVal)) {
	  cout << "chrom_"+label+"= is missing from config file" << endl;
	  exit(201);
	}
	regionChrom.push_back(iVal);
	if(!config.Get("startMb_"+label,val)) {
	  cout << "startMb_"+label+"= is missing from config file" << endl;
	  exit(201);
	}
	regionStart.push_back((long) (val*1000.*1000.));
	regionStartMb.push_back(val);
	if(!config.Get("endMb_"+label,val)) {
	  cout << "endMb_"+label+"= is missing from config file" << endl;
	  exit(202);
	}
	regionEnd.push_back((long) (val*1000.*1000.));
	regionEndMb.push_back(val);




	nRegions++;
	
	if(last != string::npos){
	  regionString=regionString.substr(last);
	}
	else regionString="";
      }
    }


 }
 double lambda=lambdaKb*1000.;

 
 cout << "Input Parameters" << endl <<endl;
 cout << setw(22) << "genoName = " << " " << genoName << endl;
 cout << setw(22) << "mapName = " << " " << mapName << endl;
 cout << setw(22) << "nStates = " << " " << nStates << endl;
 cout << setw(22) << "lambdaKb = " << " " << lambdaKb << endl;
 cout << setw(22) << "baseName = " << " " << baseName << endl;
 cout << setw(22) << "QTLName = " << " " << QTLName <<endl;
 cout << setw(22) << "regionName = " << " " << regionName <<endl;
 cout << setw(22) << "regions = " ;
 for (auto it=regions.begin(); it!=regions.end(); ++it){ 
   cout <<" " << it->first; 
 }
 cout << endl;
 for (auto it=regions.begin(); it!=regions.end(); ++it){
   string label=it->first;
   int rEffect=it->second;
   cout << setw(22) << "chrom_"+label+" = "  << " " << regionChrom[rEffect] << endl; 
   cout << setw(22) << "startMb_"+label+" = "  << " " << regionStartMb[rEffect] << endl; 
   cout << setw(22) << "endMb_"+label+" = "  << " " << regionEndMb[rEffect] << endl; 
 }
 cout << endl << endl;

 ifstream MapFile,Geno,QTLResults;
 ofstream regionFile;
 
 
 
 vector<locusMap> lociMap;
 locusMap aMapLocus;
 
 MapFile.open(mapName);
 if(MapFile.fail()) {
   cout << "Failed to open map file " + mapName << endl;
   exit(100);
 }
 getline(MapFile,line);
 int nQTLLoci=0;
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
	
	if(posVector[i]>-1) row[posVector[i]]=iVal;
	
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
  cout <<  "Number of genotyped animals: " << X.size() << endl;
  cout <<  "              Number of SNP: " << X[0].size() << endl<<endl;
  double nSeq;
  nSeq=(double) X.size();
  int nGeno=X.size();
  long nLoci=X[0].size();
  long nSNPLoci=nLoci;
  
  hmmLoci locus;
  hmm HMM(nLoci,nStates,0.5,0.5);
  HMM.lambda=lambda;
  HMM.lociMapPt=&lociMap;
  double priorCount=5;
  vector<double> doubleRow;
  

    double sum=0,pival=0.30,delta=-.05;
 
    for(int i=0;i<nStates;i++) HMM.pi[i]=1./((double) nStates);

 
  
 
    
    HMM.resetP();




    string hmmFileName;
    hmmFileName="HMMIM" + to_string(HMM.nStates) + "_" + to_string(HMM.nLoci) +".bin";
    if((failCode=HMM.read(hmmFileName))!=0) {
      cout << hmmFileName << " failed to be read  with code " << failCode << endl;
      exit(207);
    }
    
    for(int k=0;k<nStates;k++) HMM.pi[k]=1./((double) nStates);
    HMM.initPComb();

    QTLResults.open(QTLName);
  if(QTLResults.fail()) {
    cout << "Failed to open QTL results file " + QTLName << endl;
    exit(108);
  }



  cout << "Read in QTL results file" << endl;
  qtlResultLocus qLocus;
  qLocus.init(nStates);
  vector<qtlResultLocus> qtlVec;
  vector<regionQTL> regionQTLoci(nRegions);
  getline(QTLResults,line);
  while(getline(QTLResults,line)){
    stringstream linestr(line);

    linestr >> iVal;
    linestr >> aMapLocus.name;
    linestr >> aMapLocus.chrom;
    linestr >> aMapLocus.pos;
    aMapLocus.isSNP=-1;
    linestr >> qLocus.modelFreq;
    linestr >> qLocus.b;
    for(int l=0;l<nStates;l++) linestr >> qLocus.delta[l];
    int active=0;
    for (auto it=regions.begin(); it!=regions.end(); ++it){
      int rEffect=it->second;
      if(aMapLocus.chrom==regionChrom[rEffect]
	 && aMapLocus.pos >= regionStart[rEffect] 
	 &&  aMapLocus.pos <= regionEnd[rEffect]) {
	active=1;
	regionQTLoci[rEffect].qtlLoci.push_back(nQTLLoci);
      }
    }
    if(active) {
      aMapLocus.isQTL=nQTLLoci;
      qtlVec.push_back(qLocus);
      //cout << aMapLocus.isQTL <<   " " << aMapLocus.chrom << " " << aMapLocus.pos << endl;
      lociMap.push_back(aMapLocus);
      nQTLLoci++;
    }
  }

  //Midpoints

  for (auto it=regions.begin(); it!=regions.end(); ++it){
    string label=it->first;
    int rEffect=it->second;
    aMapLocus.isQTL=nQTLLoci;
    aMapLocus.name=label+"_MidPoint";
    aMapLocus.chrom=regionChrom[rEffect];
    aMapLocus.pos=(regionStart[rEffect]+regionEnd[rEffect])/2;
    aMapLocus.isSNP=-1;

    // cout << aMapLocus.isQTL <<   " " << aMapLocus.chrom << " " << aMapLocus.pos << " mp"<< endl;
    lociMap.push_back(aMapLocus);
    regionQTLoci[rEffect].midpoint=nQTLLoci;
    nQTLLoci++;
  }
    
  // Sort the Map

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
  HMM.resize(nLoci);
  for(int chr=0;chr<nChrom;chr++)  {
    HMM.loci[chromStart[chr]].newChrom=1;
    //cout << chr << " " << chromStart[chr] << " " <<chromStart[chr+1]<< " " << HMM.loci.size() << endl;
    for(int i=chromStart[chr];i<chromStart[chr+1];i++) HMM.loci[i].pos=lociMap[i].pos;
    }


  
  cout << " Number of Region QTL and midpoint loci: "<< nQTLLoci<< endl;
  cout << "                     Number of SNP loci: " << nSNPLoci << endl;
  cout << "                   Total Number of loci: " <<  nLoci << endl << endl;
  
    int nComb=HMM.nComb;
    vector<double> locusfb(nComb);
    vector<double> locusp(nStates,0);
    vector<vector<double> > f(nLoci,locusfb),b(nLoci,locusfb),p(nQTLLoci,locusfb);
    regionAnim regionAnimBase(nRegions,nStates);
    vector<regionAnim> regionAnimVec(nGeno,regionAnimBase);

   

    #pragma omp parallel for firstprivate(f,b,p)
    for(int seq=0;seq<nGeno;seq++){
      
      //cout << f.size() << " " << X[seq].size() << endl; 
	forwardVec(0, nLoci,HMM,X[seq],HMM.piComb,f);
	backwardVec(0, nLoci,HMM,X[seq],HMM.piComb,b);
	regionAnim *regionAnimpt=regionAnimVec.data();
	regionAnimpt+=seq;
	for(long i=0;i<nLoci;i++) {
	  int q=lociMap[i].isQTL;
	  if(q>-1){
	    double pSum=0,pVal;
	    int l=0;
	    p[q].assign(nStates,0.);
	    for(int I=0;I<nStates;I++) {
	      for(int J=0;J<=I;J++,l++) {
		//cout << i << " " << l<< " " << f[i].size() <<endl;
		//cout << i << " " << l<< " " << b[i].size() <<endl;
		pVal=f[i][l]*b[i][l];
		pSum+=pVal;
		p[q][I]+=pVal;
		p[q][J]+=pVal;
	      }
	    }
	    for(int I=0;I<nStates;I++) {
	      p[q][I]/=pSum;
	    }
	  }
	}// nLoci loop
	
	for (int rEffect=0;rEffect<nRegions;rEffect++){
	  regionAnimpt->gHat[rEffect]=0;
	    for(auto it=regionQTLoci[rEffect].qtlLoci.begin();it!=regionQTLoci[rEffect].qtlLoci.end();it++){
	      for(int I=0;I<nStates;I++){
		regionAnimpt->gHat[rEffect]+=p[*it][I]*(qtlVec[*it].b)*(qtlVec[*it].delta[I]);
	      }
	    }
	  for(int I=0;I<nStates;I++){
	    regionAnimpt->stateExp[rEffect][I]=p[regionQTLoci[rEffect].midpoint][I];
	  }
	}// region Loop

    }

    regionFile.open(regionName);
    if(regionFile.fail()) {
      cout << "Failed to open region file " + regionName + " for output"<< endl;
      exit(201);
    }
    regionFile << "Animal\tRegion\tChrom\tStart\tEnd\tgHat";
  for(int I=0;I<nStates;I++) regionFile << "\texpState" << I;
  regionFile << endl;
  for(int seq=0;seq<nGeno;seq++){
    for (auto it=regions.begin(); it!=regions.end(); ++it){
      string label=it->first;
      int rEffect=it->second;
      regionFile  << ID[seq] << "\t" << label << "\t"
		  << regionChrom[rEffect] << "\t"
		  << regionStart[rEffect] << "\t"
		  << regionEnd[rEffect]  << "\t"
		  << regionAnimVec[seq].gHat[rEffect];
      for(int I=0;I<nStates;I++) regionFile << "\t" << regionAnimVec[seq].stateExp[rEffect][I];
      regionFile << endl;
    }
    
  }
}

