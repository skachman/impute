#include "impute.h"

#define EIGEN_DONT_PARALLELIZE

void calcDeltaProposalFullX(const Matrix2d &sig2e,const Matrix2d &sig2b,const Matrix2d &Rinv,const Matrix2d &Binv,const vector<int> &delta,vector<vector<int> > &deltaStates,
			    vector<int> &dVec,vector<Vector2d> &rhsV,vector<Matrix2d> &lhsV,vector<Matrix2d> &lhsVs,const double pi, const double rho,const int nDeltaStates,const int nStates,
			  double &pInactive,double &maxAoverI,vector<double> &AoverIVec,double &psum){
  
  //vector<int> dsum(nDeltaStates,0);
  int nTraits=sig2e.rows();
  double logDetRinv=log(Rinv.determinant());
  Matrix2d LHS(Matrix2d::Zero());
  AoverIVec.assign(nDeltaStates*(nTraits+1),log((1.-pi)/pi)-log((double) nDeltaStates)+.5*log(Binv.determinant()));
  int sdt=0;
  for(int sd=0;sd<nDeltaStates;sd++){
 
  maxAoverI=0.;
  //dsum.assign(nDeltaStates,0);
  sdt=0;
  for(int sd=0;sd<nDeltaStates;sd++){
    Matrix2d lhs(Matrix2d::Zero());
    Vector2d rhs(Vector2d::Zero());
    dVec=deltaStates[sd];
    //if(sd != nStates) dVec[sd]=1-dVec[sd];
    int IJ=0;
    for(int I=0;I<nStates;I++){
      //dsum[sd]+=dVec[I];
      if(dVec[I]){
	rhs+=rhsV[I];
	lhs+=lhsVs[I];
      }
      for(int J=0;J<=I;J++,IJ++){
	if(dVec[I] && dVec[J]) lhs+=lhsV[IJ];      
      }
    }

  
    for(int t=0;t<nTraits;t++){
      LHS.setZero();
      LHS(t,t)=lhs(t,t);
      LHS+=Binv;
      AoverIVec[sdt]+=0.5*(rhs.dot(LHS.llt().solve(rhs))-log(LHS.determinant()))+log((1.-rho)/((double) nTraits));
      if(AoverIVec[sdt]> maxAoverI) {
	maxAoverI=AoverIVec[sdt];
      }
      sdt++;
    }
    LHS=lhs+Binv;
    AoverIVec[sdt]+=0.5*(rhs.dot(LHS.llt().solve(rhs))-log(LHS.determinant()))+log(rho);
    if(AoverIVec[sdt]> maxAoverI) {
      maxAoverI=AoverIVec[sdt];
    }
    sdt++;
    
  }
  
  psum=0;
  for(sdt=0;sdt<nDeltaStates*(nTraits+1);sdt++){
    AoverIVec[sdt]=exp(AoverIVec[sdt]-maxAoverI);
    psum+=AoverIVec[sdt];
  }
  
  pInactive=exp(-maxAoverI);
  }
}

void calcDeltaProposalX(const Matrix2d sig2e,const Matrix2d sig2b,const Matrix2d &Rinv,const Matrix2d &Binv,const vector<int> &delta,
			vector<int> &dVec,vector<Vector2d> &rhsV,vector<Matrix2d> &lhsV,vector<Matrix2d> &lhsVs,const double pi,const double rho,const int nDeltaStates,const int nStates,
			double &pInactive,double &maxAoverI,vector<double> &AoverIVec,double &psum,const double RinvLogDet,const double BinvLogDet){
  
  int nTraits=sig2e.rows();
  double logDetRinv=RinvLogDet;//log(Rinv.determinant());
  //cout << Binv(0,0) << endl;
  Matrix2d LHS(Matrix2d::Zero());
  Vector2d RHS(Vector2d::Zero());
  vector<int> dsum(nStates+1,0);
  AoverIVec.assign((nTraits+1)*(nStates+1),log((1.-pi)/pi)-log((double) nDeltaStates)+0.5*BinvLogDet);
  int sdt=0;
  /*for(int sd=0;sd<=nStates;sd++){
    for(int t=0;t<nTraits;t++){
      AoverIVec[sdt++]-=.5*log(sig2b(t,t));
    }
    AoverIVec[sdt++]-=.5*log(sig2b.determinant());
    }*/
  maxAoverI=0.;
  dsum.assign(nStates+1,0);
  
  sdt=0;
  for(int sd=0;sd<=nStates;sd++){
    Matrix2d lhs(Matrix2d::Zero());
    Vector2d rhs(Vector2d::Zero());
    dVec=delta;
    if(sd != nStates) dVec[sd]=1-dVec[sd];
    int IJ=0;
    for(int I=0;I<nStates;I++){
      dsum[sd]+=dVec[I];
      if(dVec[I]){
	rhs+=rhsV[I];
	lhs+=lhsVs[I];
      }
      for(int J=0;J<=I;J++,IJ++){
	if(dVec[I] && dVec[J]) lhs+=lhsV[IJ];      
      }
     }

      for(int t=0;t<nTraits;t++){
      LHS.setZero();
      RHS.setZero();
      RHS(t)=rhs(t);
      LHS(t,t)=lhs(t,t);
      LHS+=Binv;
      //cout << sdt << " " << sd<< endl;
      AoverIVec[sdt]+=0.5*(RHS.dot(LHS.llt().solve(RHS))-log(LHS.determinant()))+log((1.-rho)/((double) nTraits));
      if(AoverIVec[sdt]> maxAoverI) {
	maxAoverI=AoverIVec[sdt];
      }
      sdt++;
    }
    LHS=lhs+Binv;
    RHS=rhs;
    // cout << sdt << endl;
    AoverIVec[sdt]+=0.5*(RHS.dot(LHS.llt().solve(RHS))-log(LHS.determinant()))+log(rho);
    if(AoverIVec[sdt]> maxAoverI) {
      maxAoverI=AoverIVec[sdt];
    }
    sdt++;
  
  }
  
  psum=0;
  sdt=0;
  for(int sd=0;sd<=nStates;sd++){
    for(int t=0;t<=nTraits;t++){
      AoverIVec[sdt]=exp(AoverIVec[sdt]-maxAoverI);
      if(dsum[sd]==0 || dsum[sd]==nStates) AoverIVec[sdt]=0;
      psum+=AoverIVec[sdt];
      sdt++;
    }
  }
  
  pInactive=0.;
  //if(dsum[nStates]==1) pInactive=exp(-maxAoverI);
  pInactive=exp(-maxAoverI)*((double) (nStates+1))/((double) nDeltaStates+2);
}





void calcDeltaProposalSelective(const double sig2e,const double sig2b,const vector<int> &delta,
				vector<int> &dVec,vector<int> &dVecBack,
				vector<double> &rhsV,vector<double> &lhsV,vector<double> &lhsVs,
				const double pi,const int nDeltaStates,const int nStates,
				const double &inactiveProposal,double &logAlphaAA,double &logAlphaIA){
  
  
  vector<int> dsum(2,0);
  double p0,p1;
  double logPSum01,logPSum10;
  double logP01=0,logP10=0;
  double logMult=log((1.-pi)/pi)-log((double) nDeltaStates-.5*log(sig2b));
  dVec=delta;
  double pCurrent=0,logP=0;
  for(int I=0;I<nStates;I++) dsum[0]=delta[I];
  for(int sd=-1;sd<=nStates;sd++){
    double lhs=sig2e/sig2b;
    double rhs=0;
    if(sd != -1) dVec[sd]=1-dVec[sd];
    int IJ=0;
    for(int I=0;I<nStates;I++){
      if(sd==-1) dsum[0]+=dVec[I];
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
    double p=exp(0.5*(rhs*rhs/lhs-log(lhs)));
    
    if(sd==-1 || u(gen)*(pCurrent+p) < p){
      if(sd!=-1) logP01+=log(p/(pCurrent+p));
      pCurrent=p;
    }
    else{
      dVec[sd]=1-dVec[sd]; // Then dVec goes back to what it was
      logP01+=log(pCurrent/(pCurrent+p));
    }
  }
  for(int I=0;I<nStates;I++) dsum[1]+=dVec[I];
  p1=pCurrent;
  logPSum01=log(p1)+logMult-logP01;
  double logPCurrentBack=0;
  dVecBack=dVec;
  for(int sd=-1;sd<=nStates;sd++){
    double lhs=sig2e/sig2b;
    double rhs=0,p;
    if(sd != -1) dVecBack[sd]=1-dVecBack[sd];
    int IJ=0;
    for(int I=0;I<nStates;I++){
      //if(sd==-1) dsum[0]+=dVec[I];
      if(dVecBack[I]){
	rhs+=rhsV[I];
	lhs+=lhsVs[I];
      }
      for(int J=0;J<=I;J++,IJ++){
	if(dVecBack[I] && dVecBack[J]) lhs+=lhsV[IJ];      
      }
    }
    rhs/=sig2e;
    lhs/=sig2e;
    p=exp(0.5*(rhs*rhs/lhs-log(lhs)));
    
    if(sd==-1 || delta[sd]==dVecBack[sd]){
      if(sd!=-1) logP10+=log(p/(pCurrent+p));
      pCurrent=p;
    }
    else{
      dVecBack[sd]=1-dVecBack[sd]; // Then dVecBack goes back to what it was
      logP10+=log(pCurrent/(pCurrent+p));
    }
  }
  p0=pCurrent;
  logPSum10=log(p0)+logMult-logP10;
  double maxLogP=0;
  if(logPSum10>maxLogP) maxLogP=logPSum10;
  if(logPSum01>maxLogP) maxLogP=logPSum01;
  logAlphaAA=log(p1)-log(p0)+logPSum10-logPSum01+log((exp(-maxLogP)+exp(logPSum01-maxLogP))/(exp(-maxLogP)+exp(logPSum10-maxLogP)));
  logAlphaIA=log(p1)+((double) nStates)*log(2.)+maxLogP-log(exp(-maxLogP)+exp(logPSum10-maxLogP))-log(1.-inactiveProposal)-logP01;
  
  
  
 
}




void calcDeltaProposalFull(const double sig2e,const double sig2b,const vector<int> &delta,vector<vector<int> > &deltaStates,
			 vector<int> &dVec,vector<double> &rhsV,vector<double> &lhsV,vector<double> &lhsVs,const double pi,const int nDeltaStates,const int nStates,
			  double &pInactive,double &maxAoverI,vector<double> &AoverIVec,double &psum){
  
  //vector<int> dsum(nDeltaStates,0);
  AoverIVec.assign(nDeltaStates,log((1.-pi)/pi)-log((double) nDeltaStates)-.5*log(sig2b));
  maxAoverI=0.;
  //dsum.assign(nDeltaStates,0);
  for(int sd=0;sd<nDeltaStates;sd++){
    double lhs=sig2e/sig2b;
    double rhs=0;
    dVec=deltaStates[sd];
    //if(sd != nStates) dVec[sd]=1-dVec[sd];
    int IJ=0;
    for(int I=0;I<nStates;I++){
      //dsum[sd]+=dVec[I];
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
    if(AoverIVec[sd]> maxAoverI) {
       maxAoverI=AoverIVec[sd];
    }
  }
  
  psum=0;
  for(int sd=0;sd<nDeltaStates;sd++){
    AoverIVec[sd]=exp(AoverIVec[sd]-maxAoverI);
    //if(dsum[sd]==0 || dsum[sd]==nStates) AoverIVec[sd]=0;
    psum+=AoverIVec[sd];
  }
  
  pInactive=exp(-maxAoverI);
}

void calcDeltaProposal(const double sig2e,const double sig2b,const vector<int> &delta,
			 vector<int> &dVec,vector<double> &rhsV,vector<double> &lhsV,vector<double> &lhsVs,const double pi,const int nDeltaStates,const int nStates,
			 double &pInactive,double &maxAoverI,vector<double> &AoverIVec,double &psum){
  
  vector<int> dsum(nStates+1,0);
  AoverIVec.assign(nStates+1,log((1.-pi)/pi)-log((double) nDeltaStates)-.5*log(sig2b));
  maxAoverI=0.;
  dsum.assign(nStates+1,0);
  for(int sd=0;sd<=nStates;sd++){
    double lhs=sig2e/sig2b;
    double rhs=0;
    dVec=delta;
    if(sd != nStates) dVec[sd]=1-dVec[sd];
    int IJ=0;
    for(int I=0;I<nStates;I++){
      dsum[sd]+=dVec[I];
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
    if(AoverIVec[sd]> maxAoverI) {
       maxAoverI=AoverIVec[sd];
    }
  }
  
  psum=0;
  for(int sd=0;sd<=nStates;sd++){
    AoverIVec[sd]=exp(AoverIVec[sd]-maxAoverI);
    if(dsum[sd]==0 || dsum[sd]==nStates) AoverIVec[sd]=0;
    psum+=AoverIVec[sd];
  }
  
  pInactive=0.;
  //if(dsum[nStates]==1) pInactive=exp(-maxAoverI);
  pInactive=exp(-maxAoverI)*((double) (nStates+1))/((double) nDeltaStates+2);
}


void qtlLocus::updateSum(const qtlLocus &A){
  
  if(A.active){
    active++;
    gVar+=A.gVar;
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
  gVar=0;
  delta.assign(ns,0);
  active=0;

}
  

void qtlLocus::init(int ns,double sig2b,double pi){
  default_random_engine gen;
  uniform_real_distribution<double> u(0.,1.);
  normal_distribution<double> z(0.,1.);
 
  b=0;
  gVar=0;
  delta.assign(ns,0);
  active=0;
  if(u(gen)>pi) {
    active=1;
    b=sqrt(sig2b)*z(gen);
    for(int i=0;i<ns;i++) {
      if(u(gen) < .5) delta[i]=1;
    }
  }
}


void qtlXLocus::init(int ns,int nt){
  activeTrait.assign(nt,0);
  b.resize(nt);
  b.setZero();
  delta.assign(ns,0);
  t=nt;
  active=0;

}
void qtlXLocusSum::init(int ns,int nt){
  activeTrait.assign(nt,0);
  b.resize(nt);
  b.setZero();
  all=0;
  delta.resize(nt);
  for(int l=0;l<nt;l++) delta[l].assign(ns,0);
  t=nt;
  active=0;

}


void qtlXLocusSum::updateSum(const qtlXLocus &A,const int nt){
  int nact=0;
  if(A.active){
    active++;
    for(int l=0;l<nt;l++){
      if(A.activeTrait[l]){
	nact++;
	activeTrait[l]++;
	if(A.b[l] >0){
	  b[l]+=A.b[l];
	  for(int i=0;i<A.delta.size();i++){
	    delta[l][i]+=A.delta[i];
	  }
	}
	else{
	  b[l]-=A.b[l];
	  for(int i=0;i<A.delta.size();i++){
	    delta[l][i]+=1-A.delta[i];
	  }
	}
      }
    }
    if(nact==nt) all++;
  }
}

void qtlXLocus::init(int ns,int nt,Matrix2d &bL,double pi,double rho){
  
  activeTrait.assign(nt,0);
  b.resize(nt);
  b.setZero();
  delta.assign(ns,0);
  active=0;
  t=0;
  if(u(gen)>pi) {
    active=1;
    double val;
    for(int j=0;j<nt;j++){
      val=Z(gen);
      for(int i=j;i<nt;i++){
	b[i]+=bL(i,j)*val;
      }
    }
    for(int i=0;i<ns;i++) {
      if(u(gen) < .5) delta[i]=1;
    }
    if(u(gen) < rho) {
      t=nt;
      for(int i=0;i<nt;i++) activeTrait[i]=1;
    }
    else{
      t=nt*u(gen);
      activeTrait[t]=1;
    }
    
    
  }
}


void forward(const long start, const long end,hmm &HMM, const vector<int> &X,const vector<double> &picomb){
  int nComb=HMM.nComb,nStates=HMM.nStates;
  vector<locusMap> *lociMapPt=HMM.lociMapPt;
  double S=(double) nStates,lambda=HMM.lambda;
  lambda*=(S-1.)/S;
  double sum=0,val;
  vector<double> transSingle(2),f(nComb),fnew(nComb); 
  for(int l=0;l<nComb;l++){
    val=picomb[l];
    int snpLoc=(*lociMapPt)[start].isSNP;
    if(snpLoc>=0) val*=HMM.emit(snpLoc,l,X[snpLoc]);
    f[l]=val;
    sum+=val;

  }
  for(int l=0;l<nComb;l++) {
    f[l]/=sum;
    HMM.loci[start].f[l]=f[l];//scale
  }
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
	    
	    for(int J2=0;J2<I2;J2++,k++){
	      val+=f[k]*(transSingle[I1==I2]*transSingle[J1==J2]+transSingle[I1==J2]*transSingle[J1==I2]);
	    }
	    val+=f[k]*transSingle[I1==I2]*transSingle[J1==I2];
	    k++;
	  }
	  int snpLoc=(*lociMapPt)[i].isSNP;
	  if(snpLoc>=0) val*=HMM.emit(snpLoc,l,X[snpLoc]);
	}
	fnew[l]=val;
	sum+=val;
      }
    }
    for(int l=0;l<nComb;l++) {
      f[l]=fnew[l]/sum;
      HMM.loci[i].f[l]=f[l];//scale
    }
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
    //if(HMM.loci[i+1].newChrom) eMult=0.;
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

int writeX(const string &filename,vector<vector<int> > &X,vector<string > &ID){
  ofstream XStream;
  long nLoci,nSeq;
  
  vector<int> row;
  nSeq=X.size();
  nLoci=X[0].size();
  XStream.open(filename,ofstream::binary);
  if(XStream.fail()) return(1);
  XStream.write(reinterpret_cast<char*>(&nSeq),sizeof(nSeq));
  XStream.write(reinterpret_cast<char*>(&nLoci),sizeof(nLoci));
  for(long i=0;i<nSeq;i++) {
    int sl=ID[i].length();
      XStream.write(reinterpret_cast<char*>(&sl),sizeof(sl));
      XStream.write(reinterpret_cast<const char*>(ID[i].c_str()),sl*sizeof(char));
      XStream.write(reinterpret_cast<char*>(X[i].data()),nLoci*sizeof(int));
  }
  
  return(0);
}


int readX(const string &filename,vector<vector<int> > &X,vector<string > &ID,idmap &seqmap){
  ifstream XStream;
  long nLoci,nSeq;
  vector<int> row;
  
  XStream.open(filename,ifstream::binary);
  if(XStream.fail()) return(1);
  XStream.read(reinterpret_cast<char*>(&nSeq),sizeof(nSeq));
  XStream.read(reinterpret_cast<char*>(&nLoci),sizeof(nLoci));
  row.resize(nLoci);
  X.resize(nSeq,row);
  ID.resize(nSeq);
  for(long i=0;i<nSeq;i++) {
    int sl;
      XStream.read(reinterpret_cast<char*>(&sl),sizeof(sl));
      char* buffer = new char[sl];
      XStream.read(reinterpret_cast<char*>(buffer),sl*sizeof(char));
      ID[i].assign(buffer,sl);
      delete[] buffer;
      seqmap[ID[i]]=i;
      XStream.read(reinterpret_cast<char*>(X[i].data()),nLoci*sizeof(int));
  }
  
  return(0);
}

int hmm::write(const string &filename){
  ofstream hmmStream;
  hmmStream.open(filename,ofstream::binary);
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
  hmmStream.open(filename,ifstream::binary);
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


 

void recombVec(hmm &HMM, const vector<int> &X,const vector<double> &picomb,const vector< vector<double> > &fVec,const vector< vector<double> > &bVec,vector<double> &rf){
  int nComb=HMM.nComb,nStates=HMM.nStates;
  vector<locusMap> *lociMapPt=HMM.lociMapPt;
  double S=(double) nStates,lambda=HMM.lambda;
  lambda*=(S-1.)/S;
  
  vector<double> transSingle(2),transDouble(3),f(nComb),fnew(nComb),valVecI(nStates),valVecJ(nStates); 
}

 void forwardVec(const long start, const long end,hmm &HMM, const vector<int> &X,const vector<double> &picomb,vector< vector<double> > &fVec){
  int nComb=HMM.nComb,nStates=HMM.nStates;
  vector<locusMap> *lociMapPt=HMM.lociMapPt;
  double S=(double) nStates,lambda=HMM.lambda;
  lambda*=(S-1.)/S;
  double sum=0,val;
  
  vector<double> transSingle(2),transDouble(3),f(nComb),fnew(nComb),valVecI(nStates),valVecJ(nStates); 
  for(int l=0;l<nComb;l++){
    val=picomb[l];
    int snpLoc=(*lociMapPt)[start].isSNP;
    if(snpLoc>=0) val*=HMM.emit(snpLoc,l,X[snpLoc]);
    f[l]=val;
    sum+=val;

  }
  for(int l=0;l<nComb;l++) {
    f[l]/=sum;
    fVec[start][l]=f[l];//scale
  }
  for(int i=start+1;i<end;i++){
    double eMult=exp(-((double)(HMM.loci[i].pos-HMM.loci[i-1].pos))/lambda);

    transSingle[1]=(1.+(S-1.)*eMult)/S;
    transSingle[0]=(1.-eMult)/S;
    transDouble[0]= transSingle[0]* transSingle[0];
    transDouble[1]= transSingle[0]* transSingle[1] - transDouble[0];
    transDouble[2]= transSingle[1]* transSingle[1]-2.*transDouble[1]-transDouble[0];
    if(0 ) cout << i << " " << (HMM.loci[i].pos-HMM.loci[i-1].pos)/lambda <<  " " << eMult << " " << transSingle[0] << " " << transSingle[1] << endl;
    sum=0;
    val=0;
    int l;
    if(HMM.loci[i].newChrom){
      for(int l=0;l<nComb;l++) fnew[l]=picomb[l];
    }
    else{
      l=0;
      val=0.;
      for(int I2=0;I2<nStates;I2++){
	for(int J2=0;J2<I2;J2++,l++){
	  val+=2.*f[l];
	}
	val+=f[l];
	l++;
      }
      val*=transDouble[0];
      l=0;
      for(int I1=0;I1<nStates;I1++){
	for(int J1=0;J1<I1;J1++,l++){
	  fnew[l]=val;
	}
	fnew[l]=val;
	l++;
      }
      valVecI.assign(nStates,0.);
      valVecJ.assign(nStates,0.);
      l=0;
      for(int I2=0;I2<nStates;I2++){
	for(int J2=0;J2<I2;J2++,l++){
	  valVecI[I2]+=f[l];
	  valVecJ[J2]+=f[l];
	  valVecJ[I2]+=f[l];
	  valVecI[J2]+=f[l];
	}
	  valVecI[I2]+=f[l];
	  valVecJ[I2]+=f[l];
	  l++;
      }
      for(int I=0;I<nStates;I++){
	valVecI[I]*=transDouble[1];
	valVecJ[I]*=transDouble[1];
      }
      l=0;
      for(int I1=0;I1<nStates;I1++){
	for(int J1=0;J1<I1;J1++,l++){
	  fnew[l]+=valVecI[I1];
	  fnew[l]+=valVecJ[J1];
	}
	fnew[l]+=valVecI[I1];
	fnew[l]+=valVecJ[I1];
	l++;
      }
      l=0;
      for(int I1=0;I1<nStates;I1++){
	for(int J1=0;J1<=I1;J1++,l++){
	  fnew[l]+=f[l]*transDouble[2];
	}
      }
      
    }
    sum=0;
    for(int l=0;l<nComb;l++) {
      int snpLoc=(*lociMapPt)[i].isSNP;
      if(snpLoc>=0) fnew[l]*=HMM.emit(snpLoc,l,X[snpLoc]);
      sum+=fnew[l];
    }
    /*   
      
      l=0;
    for(int I1=0;I1<nStates;I1++){
      for(int J1=0;J1<=I1;J1++,l++){
	if(HMM.loci[i].newChrom){
	  val=picomb[l];

	 
	}
	else{
	  val=0;
	  int k=0;
	  for(int I2=0;I2<nStates;I2++){
	    
	    for(int J2=0;J2<I2;J2++,k++){
	      val+=f[k]*(transSingle[I1==I2]*transSingle[J1==J2]+transSingle[I1==J2]*transSingle[J1==I2]);
	      //val+=f[k]*(transDouble[I1==I2+J1==J2]+transDouble[I1==J2+J1==I2]);
	    }
	    val+=f[k]*transSingle[I1==I2]*transSingle[J1==I2];
	    //val+=f[k]*transDouble[I1==I2+J1==I2];
	    k++;
	  }
	}
	int snpLoc=(*lociMapPt)[i].isSNP;
	if(snpLoc>=0) val*=HMM.emit(snpLoc,l,X[snpLoc]);
	fnew[l]=val;
	sum+=val;
      }
    }
    */
    for(int l=0;l<nComb;l++) {
      f[l]=fnew[l]/sum;
      fVec[i][l]=f[l];//scale
    }
  }
}



void backwardVec(const long start, const long end,hmm &HMM, const vector<int> &X,const vector<double> &picomb, vector<vector <double>> &bVec){
  int nComb=HMM.nComb,nStates=HMM.nStates;
  vector<locusMap> *lociMapPt=HMM.lociMapPt;
  double S=(double) nStates,lambda=HMM.lambda;
  vector<double> b(nComb),bNew(nComb);
  lambda*=(S-1.)/S;
   vector<double> transSingle(2),transDouble(3),valVecI(nStates),valVecJ(nStates);
   vector<double> eProbVec(nComb);
  double sum,val;
  sum=(double )nComb;
  for(int l=0;l<nComb;l++) {
    bVec[end-1][l]=1./sum;
    b[l]=1./sum;
  }
  double eMult=0;
  for(long i=end-2;i>=start;i--){
    if(!HMM.loci[i+1].newChrom) eMult=exp(-((double)(HMM.loci[i+1].pos-HMM.loci[i].pos))/lambda);
    transSingle[1]=(1.+(S-1.)*eMult)/S;
    transSingle[0]=(1.-eMult)/S;
    
    transDouble[0]= transSingle[0]* transSingle[0];
    transDouble[1]= transSingle[0]* transSingle[1] - transDouble[0];
    transDouble[2]= transSingle[1]* transSingle[1]-2.*transDouble[1]-transDouble[0];
    int k=0;
    val=0;
    int l=0;
    for(int I2=0;I2<nStates;I2++){
      for(int J2=0;J2<=I2;J2++,l++){
	double eProb=1.;
	int snpLoc=(*lociMapPt)[i+1].isSNP;
	if(snpLoc>=0) eProb=HMM.emit(snpLoc,l,X[snpLoc]);
	eProbVec[l]=eProb;
	if(HMM.loci[i+1].newChrom){
	  val+=picomb[l]*b[l]*eProb;	      
	}
	else{
	  val+=transDouble[0]*b[l]*eProb; 
	  if(I2 !=J2)val+=transDouble[0]*b[l]*eProb; 
	}
      }
    }
    k=0;
    for(int I2=0;I2<nStates;I2++){
      for(int J2=0;J2<I2;J2++,k++){
	bNew[k]=val;
      }
      bNew[k]=val;
      k++;
    }
    if(!HMM.loci[i+1].newChrom){
      valVecI.assign(nStates,0.);
      valVecJ.assign(nStates,0.);
      l=0;
      for(int I1=0;I1<nStates;I1++){
	for(int J1=0;J1<I1;J1++,l++){
	  valVecI[I1]+=b[l]*eProbVec[l];
	  valVecJ[J1]+=b[l]*eProbVec[l];
	  valVecJ[I1]+=b[l]*eProbVec[l];
	  valVecI[J1]+=b[l]*eProbVec[l];
	}
	  valVecI[I1]+=b[l]*eProbVec[l];
	  valVecJ[I1]+=b[l]*eProbVec[l];
	  l++;
      }
      for(int I=0;I<nStates;I++){
	valVecI[I]*=transDouble[1];
	valVecJ[I]*=transDouble[1];
      }
      l=0;
      for(int I1=0;I1<nStates;I1++){
	for(int J1=0;J1<I1;J1++,l++){
	  bNew[l]+=valVecI[I1];
	  bNew[l]+=valVecJ[J1];
	}
	bNew[l]+=valVecI[I1];
	bNew[l]+=valVecJ[I1];
	l++;
      }
      l=0;
      for(int I1=0;I1<nStates;I1++){
	for(int J1=0;J1<=I1;J1++,l++){
	  bNew[l]+=b[l]*transDouble[2]*eProbVec[l];
	}
      }
    }
    /*  
    k=0;
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
	      val+=picomb[l]*b[l]*eProb;	      
	    }
	    else{
	      val+=transSingle[I1==I2]*transSingle[J1==J2]*b[l]*eProb; 
	      if(I2 !=J2)val+=transSingle[I1==J2]*transSingle[J1==I2]*b[l]*eProb; 
	    }
	  }
	}
	bNew[k]=val;
	sum+=val;
      }
    }
    */
    
    sum=0;
    for(k=0;k<nComb;k++){
      sum+=bNew[k];
    }
    for(int k=0;k<nComb;k++) {
      b[k]=bNew[k]/sum;
      bVec[i][k]=b[k];
    }
  }
  
}


void rWishartX(const Matrix2d &Sigma, const double n,Matrix2d &SL,Matrix2d &S){
  int p=Sigma.rows();
  MatrixXd A = Sigma.llt().matrixL();
  MatrixXd L(Matrix2d::Zero());
  double alpha=n/2.;
  for(int i=0;i<p;i++){
    for(int j=0;j<i;j++){
      L(i,j)=Z(gen);
    }
    gamma_distribution<double> X2(alpha,2.);
    L(i,i)=sqrt(X2(gen));
    alpha-=0.5;
  }
  SL=A*L;
  S=SL*SL.transpose();
}

void rWishartX(const Matrix2d &Sigma, const double n,Matrix2d &S){
  int p=Sigma.rows();
  MatrixXd A = Sigma.llt().matrixL();
  MatrixXd L(Matrix2d::Zero());
  double alpha=n/2.;
  for(int i=0;i<p;i++){
    for(int j=0;j<i;j++){
      L(i,j)=Z(gen);
    }
    gamma_distribution<double> X2(alpha,2.);
    L(i,i)=sqrt(X2(gen));
    alpha-=0.5;
  }
  MatrixXd SL=A*L;
  S=SL*SL.transpose();
}


// Sorensen and Gianola (2002) Chapter 14 page 626 see also Korsgaard et al. (1999)
void rInvWishartCond2d(const Matrix2d &Sigma, const double n,Matrix2d &S){

  gamma_distribution<double> X2(n/2.,2.*Sigma(1,1));
  double x1=X2(gen);
  double x2=Sigma(0,1)/Sigma(1,1)+Z(gen)*sqrt((Sigma(0,0)-Sigma(0,1)*Sigma(1,0)/Sigma(1,1))/x1);
  
  S(1,1)=x2*x2+1./x1;
  S(0,1)=-x2;
  S(1,0)=-x2;
  S(0,0)=1.;
}


void Zvec(Vector2d &zVec,const int n){
  /*if(n){
    zVec.resize(n);
    }*/
  int N=zVec.size();
  for(int i=0;i<N;i++) zVec(i)=Z(gen); 
}


void printLower(ostream &out,const MatrixXd & A,const string sep){
  for(int i=0;i<A.rows();i++)
    for(int j=0;j<=i;j++)
      out << sep << A(i,j);
}
void printVector(ostream &out,const VectorXd & A,const string sep){
  for(int i=0;i<A.size();i++)
      out << sep << A[i];
}
