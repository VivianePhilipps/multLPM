#include <math.h>
#include <RcppArmadillo.h>

extern "C"{
  void sadmvn_(int *n, double *lower, double *upper,
	       int *infin, double *correl,
	       int *maxpts, double *abseps, double *releps,
	       double *error, double *value, int *inform );
};


// ordre parametres:
// submod1, .., submodK, covRE, seuil1, gammas1, seuil2, gamma2
// avec prm submodl : nef+ncontr+nvc+ncor+ny+nalea+ntrtot
//      seuil : seuil, seuilX1, ..., seuilXnvarD
//      gammas : gamma1, gamma1X1, ..., gamma1XnvarD, ..., gammaK, gammaKX1, ..., gammaKXnvarD

RcppExport SEXP loglikUACV(SEXP b0, SEXP bfix0, SEXP fix0, SEXP Y0, SEXP X0, SEXP Xd0, SEXP idD0, SEXP D0, SEXP Xseuil0, SEXP nmes0, SEXP nv0, SEXP idx0, SEXP idiag0, SEXP ncor0, SEXP ny0, SEXP nalea0, SEXP ntr0, SEXP link0, SEXP nodes0, SEXP epsY0, SEXP nRE0, SEXP nBM0, SEXP nbevt0, SEXP idcause0, SEXP entreRetard0, SEXP nbTronc0)
{
  using namespace arma;


  Rcpp::NumericVector b00(b0);
  Rcpp::NumericVector bfix00(bfix0);
  Rcpp::IntegerVector fix00(fix0);
  Rcpp::NumericVector Y00(Y0);
  Rcpp::List X(X0);
  Rcpp::List Xd(Xd0);
  Rcpp::IntegerVector idD00(idD0);
  Rcpp::NumericMatrix D00(D0);
  Rcpp::IntegerMatrix nmes00(nmes0);
  Rcpp::IntegerVector nv00(nv0);
  Rcpp::IntegerMatrix idx00(idx0);
  Rcpp::IntegerVector idiag00(idiag0);
  Rcpp::IntegerVector ncor00(ncor0);
  Rcpp::IntegerVector ny00(ny0);
  Rcpp::IntegerVector nalea00(nalea0);
  Rcpp::IntegerVector ntr00(ntr0);
  Rcpp::IntegerVector link00(link0);
  Rcpp::NumericVector nodes00(nodes0);
  Rcpp::NumericVector epsY00(epsY0);
  int nRE = Rcpp::as<int>(nRE0);
  int nBM = Rcpp::as<int>(nBM0);
  Rcpp::NumericMatrix Xseuil00(Xseuil0);
  int nbevt = Rcpp::as<int>(nbevt0);
  Rcpp::IntegerVector idcause00(idcause0);
  Rcpp::IntegerVector entreRetard00(entreRetard0);
  Rcpp::IntegerMatrix nbTronc00(nbTronc0);

  int ns = nmes00.nrow();
  int K = ny00.size();
  int nvtot = idx00.ncol();
  int nvarD = D00.ncol()-1;
  int ndept = idcause00.size()/nbevt - nvarD; 
  int M = nmes00.ncol()-nbevt-1;
  
  vec b(b00.begin(),b00.size(),false);
  vec bfix(bfix00.begin(),bfix00.size());
  ivec fix(fix00.begin(),fix00.size());
  imat nmes(nmes00.begin(),ns,M+nbevt+1);
  ivec nv(nv00.begin(),nv00.size());
  imat idx(idx00.begin(),4,nvtot);
  ivec idiag(idiag00.begin(),idiag00.size());
  ivec ncor(ncor00.begin(),ncor00.size());
  ivec ny(ny00.begin(),ny00.size());
  ivec nalea(nalea00.begin(),nalea00.size());
  vec Y(Y00.begin(),Y00.size());
  ivec link(link00.begin(),link00.size());
  vec nodes(nodes00.begin(),nodes00.size());
  ivec ntr(ntr00.begin(),ntr00.size());
  vec epsY(epsY00.begin(),epsY00.size());
  mat D(D00.begin(),ns,nvarD+1);
  mat Xseuil(Xseuil00.begin(),Xseuil00.nrow(),Xseuil00.ncol());
  ivec idcause(idcause00.begin(),idcause00.size());
  ivec idD(idD00.begin(),idD00.size());
  ivec entreRetard(entreRetard00.begin(),entreRetard00.size());
  imat nbTronc(nbTronc00.begin(),ns,nbevt);
  int npm = b.size();
  int npmfix = sum(fix);
  int npmtot = npm+npmfix;
  vec btot(npmtot);

  int i;
  int iest=0;
  int ifix=0;

  for(i=0; i<npmtot; i++)
    {
      if(fix(i)==0) 
  	{
  	  btot(i) = b(iest);
  	  iest += 1;
  	}
      else
  	{
  	  btot(i) = bfix(ifix);
  	  ifix += 1;
  	}
    }

  // creer nef ncontr nea nvc ntrtot subnpm
  ivec nef(K);
  ivec ncontr(K);
  ivec nvarcontr(K);
  ivec nea(K);
  ivec nvc(K);
  ivec ntrtot(K);
  ivec subnpm(K);

  nef.fill(0);
  ncontr.fill(0);
  nvarcontr.fill(0);
  nea.fill(0);
  nvc.fill(0);
  ntrtot.fill(0);
  subnpm.fill(0);

  int j,k,l,p,q,ll,kk,m,endj,l0,m0,jexcl1,jexcl2;
  int jcurr,jdcurr,iidx,jseuil1curr,jseuil2curr,nmesdcurr,ncortotBM,npmMM;
  int nvarD1,nvarD2,ndept1,ndept2;
  int j1,j2,j3,kbm,kkbm,nu;
  int endj2,tmpntr,pcontr,ibeta,ibcontr,imatB,imatBM,ibtot,iRE,iBM,ispl,tmpcontr;
  int inodes,tmpnmes,ni,nik,ni0,ni01,ni02;

  double sig1,sig2,corr,yy,fct;

  double aatmp,bbtmp,aa,bb,cc,dd,jac;
  double ht,htm,ht2,ht3,h1,hh,h2,h3,h2n,hn,hht,mm,mm1,mm2,im,im1,im2,som;
  double vraistot,vrais_i,vrais_y,vrais_d;
  vec loglik;

  vec splaa,zitr,muK,muD,beta,bcontr,eta1,eta2,eta,eta1_dept,eta2_dept,eta_dept,varD;
  ivec nobs,nmescurr,nmescurry,niparK;
  mat gamma1,gamma2,B,U,BM,Ubm;
  vec Hyi,muiK,tcor,muiD,tdcor,yim,Hyim,ytemp1,ytemp3;
  mat RiK,Xi,Xcontri,Zi,RiD,Xdi,Zdi,RiKD,Xk,Xdk,VK,invVK;
  Rcpp::NumericVector ytemp2;
  int dimPhi,maxpts,inform,nTronc1,nTronc2;
  vec upper,lower,delta,seuil_i,upper1;
  ivec infin;
  uvec js;
  mat Gam1,GG,Gam2,Gam,tGam,GamTronc;
  mat VKD,VD,VD_K,V0,V1,I0,I1,ZCC,P,H;
  vec muiD_K,correl,un,h;
  double abseps,releps,error,value0,value1,valuer,retardi;
  std::vector<double> lower2,upper2,correl2,delta2;
  std::vector<int> infin2;
  rowvec g;

  loglik.zeros(2*ns);

  l0=0;
  m0=0;
  ncortotBM=0;
  for(k=0; k<K; k++)
    {
      nef(k) = -1; // pour ne pas compter l'intercept

      for(l=l0; l<l0+nv(k); l++)
	{
	  if(idx(0,l)==1) 
	    {
	      nef(k) += 1;  
	    }

	  if(idx(1,l)==1) 
	    {
	      nvarcontr(k) += 1;
	      ncontr(k) += ny(k)-1;
	    }

	  if(idx(2,l)==1) 
	    {
	      nea(k) += 1;	 
	    }
	}

      if(idiag(k)==1) 
	{
	  nvc(k)=nea(k)-1;
	}
      else
	{
	  nvc(k) = nea(k)*(nea(k)+1)/2-1;
	}

      for(m=m0; m<m0+ny(k); m++)
	{
	  if(link(m)==0)
	    {
	      ntrtot(k) += 2;
	    }
	  if(link(m)==1)
	    {
	      ntrtot(k) += 4;
	    }
	  if(link(m)==2)
	    {
	      ntrtot(k) += ntr(m);
	    }
	}

      if(ncor(k)==1)
	{
	  ncortotBM += 1;
	}

      l0 += nv(k);
      m0 += ny(k);
      subnpm(k) = nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+ntrtot(k);
    }


  splaa.zeros(max(ntrtot)*K);

  npmMM = sum(subnpm);

  nobs = sum(nmes,0).t(); // nb d'obs por chaque Y

  // nvar par evt
  nvarD1 = 0;
  nvarD2 = 0;
  ndept1 = 0;
  ndept2 = 0; 
  endj=(nvarD+ndept)*nbevt;
  
  for(j=0;j<endj;j++)
    {
      if(j<(nvarD+ndept)) //evt1
	{
	  if(idcause(j)==1)
	    {
	      if(j<nvarD)
		{
		  nvarD1 += 1;
		}
	      else
		{
		  ndept1 += 1;
		}
	    }
	}
      else //evt2
	{
	  if(idcause(j)==1)
	    {
	      if(j<(nvarD+ndept)+nvarD)
		{
		  nvarD2 += 1;
		}
	      else
		{
		  ndept2 += 1;
		}
	    }
	}
    }


  beta.zeros(sum(nef)+K);
  bcontr.zeros(sum(ncontr)+sum(nvarcontr));
  eta1.zeros(nvarD1+1);
  eta2.zeros(nvarD2+1);
  eta.zeros(nbevt+nvarD1+nvarD2);
  eta1_dept.zeros(ndept1);
  eta2_dept.zeros(ndept2);
  eta_dept.zeros(ndept1+ndept2);
  gamma1.zeros(K,nvarD1+1);
  gamma2.zeros(K,nvarD2+1);
  varD.zeros(nvarD1+nvarD2+nbevt);
  B.zeros(sum(nea),sum(nea));
  U.zeros(sum(nea),sum(nea));
  BM.zeros(ncortotBM,ncortotBM);
  Ubm.zeros(ncortotBM,ncortotBM);


  ibeta=0;
  ibcontr=0;
  imatB=0;
  imatBM=0;
  ibtot=0;
  iRE=0;
  iBM=0;
  ispl=0;
  tmpcontr=0;



  // creer vecteurs beta (effets fixes des K MM)
  //       vecteur  eta (K seuils) 
  //       matrice B (mat par blocs, varcov des effets aleatoires des K MM)
  for(k=0; k<K; k++)
    {
      // effets fixes
      beta(ibeta) = 0;
      beta(span(ibeta+1,ibeta+nef(k))) = btot(span(ibtot,ibtot+nef(k)-1));

      // contrastes
      for(l=0; l<nvarcontr(k); l++)
	{
	  bcontr(span(ibcontr,ibcontr+ny(k)-2)) = btot(span(ibtot+nef(k)+tmpcontr,ibtot+nef(k)+tmpcontr+ny(k)-2));
	  bcontr(ibcontr+ny(k)-1) = -sum(btot(span(ibtot+nef(k)+tmpcontr,ibtot+nef(k)+tmpcontr+ny(k)-2)));

	  tmpcontr += ny(k)-1;
	  ibcontr += ny(k);
	}

      // varcov du sous-modele k 
      U(imatB,imatB) = 1; // var(u0)=1
      if(idiag(k)==1)
  	{
  	  for(j=1; j<nea(k); j++)
  	    {
  	      U(imatB+j,imatB+j) = btot(ibtot+nef(k)+ncontr(k)+j-1);
  	    }
  	}
      else
  	{
  	  j=0;
  	  for(j2=0; j2<nea(k); j2++)
  	    {
	      if(j2>0)
		{
		  for(j1=0; j1<=j2; j1++)
		    {
		      U(imatB+j1,imatB+j2) = btot(ibtot+nef(k)+ncontr(k)+j);
		      j += 1;
		    }
		}

	      if((imatB>0) & (nRE>0))
		{
		  for(j3=0; j3<sum(nea(span(0,k-1))); j3++)
		    {
		      U(j3,imatB+j2) = btot(npmMM+iRE);
		      iRE +=1;
		    }
		}
  	    }
  	}
      
      // BM du sous-modele k
      if(ncor(k)==1)
	{
	  Ubm(imatBM,imatBM) = btot(ibtot+nef(k)+ncontr(k)+nvc(k));

	  if((imatBM>0) & (nBM>0))
	    {
	      for(j3=0; j3<imatBM; j3++)
		{
		  Ubm(j3,imatBM)=btot(npmMM+nRE+iBM);
		  iBM += 1;
		}
	    }

	  imatBM += 1;
	}

	
      ibeta += nef(k)+1;
      ibtot += nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+ntrtot(k);
      imatB += nea(k);
    }


  // effets aleatoires
  endj2=U.n_cols;
  for(j2=0;j2<endj2;j2++)
    {
      sig2 = U(j2,j2);

      for(j1=0;j1<=j2;j1++)
	{
	  if(j1==j2)
	    {
	      B(j1,j2)=sig2*sig2; //? ou exp(sig2)?
	    }
	  else
	    {
	      corr=(exp(U(j1,j2))-1)/(1+exp(U(j1,j2)));
	      sig1 = U(j1,j1);
	      B(j1,j2) = corr*fabs(sig1)*fabs(sig2);
	      B(j2,j1)=B(j1,j2);
	    }
	}
    }

  // BM
  endj2=Ubm.n_cols;
  for(j2=0;j2<endj2;j2++)
    {
      sig2 = Ubm(j2,j2);

      for(j1=0;j1<=j2;j1++)
	{
	  if(j1==j2)
	    {
	      BM(j1,j2)=sig2*sig2; //? ou exp(sig2)?
	    }
	  else
	    {
	      corr=(exp(Ubm(j1,j2))-1)/(1+exp(Ubm(j1,j2)));
	      sig1 = Ubm(j1,j1);
	      BM(j1,j2) = corr*fabs(sig1)*fabs(sig2);
	      BM(j2,j1)=BM(j1,j2);

	      if(nBM>0) 
		{
		  btot(npmMM+nRE+iBM) = BM(j1,j2);
		  iBM += 1;
		}
	    }
	}
    }

  eta1 = btot(span(npmMM+nRE+nBM,npmMM+nRE+nBM+nvarD1));
  eta(span(0,nvarD1)) = eta1;
  if(ndept1>0) 
    {
      eta1_dept = btot(span(npmMM+nRE+nBM+1+nvarD1,npmMM+nRE+nBM+1+nvarD1+ndept1-1));
      eta_dept(span(0,ndept1-1)) = eta1_dept;
    }
  if(nbevt==2)
    {
      eta2 = btot(span(npmMM+nRE+nBM+1+nvarD1+ndept1+K*(1+nvarD1),npmMM+nRE+nBM+1+nvarD1+ndept1+K*(1+nvarD1)+nvarD2));
      if(ndept2>0) 
	{
	  eta2_dept = btot(span(npmMM+nRE+nBM+1+nvarD1+ndept1+K*(1+nvarD1)+1+nvarD2,npmMM+nRE+nBM+1+nvarD1+ndept1+K*(1+nvarD1)+1+nvarD2+ndept2-1));
	  eta_dept(span(ndept1,ndept1+ndept2-1)) = eta2_dept;
	}
      eta(span(nvarD1+1,nvarD1+nvarD2+1)) = eta2;
    }

  l=0;
  for(k=0;k<K;k++)
    {
      gamma1(k,span(0,nvarD1)) = trans(btot(span(npmMM+nRE+nBM+1+nvarD1+ndept1+l,npmMM+nRE+nBM+1+nvarD1+ndept1+l+nvarD1)));
      l += nvarD1+1;
    }

  if(nbevt==2)
    {
      l=0;
      for(k=0;k<K;k++)
	{
	  gamma2(k,span(0,nvarD2)) = trans(btot(span(npmMM+nRE+nBM+1+nvarD1+ndept1+K*(1+nvarD1)+1+nvarD2+ndept2+l,npmMM+nRE+nBM+1+nvarD1+ndept1+K*(1+nvarD1)+1+nvarD2+ndept2+l+nvarD2)));
	  l += nvarD2+1;
	}
    }

  
  // nmescurr pour  X
  // nmesdcurr pour Xd
  // nmescurry pour Y
  nmescurr=zeros<ivec>(K);
  nmescurry=zeros<ivec>(K);
  nmescurr(0)=0;
  nmescurry(0)=0;
  m0=0;
  if(K>1)
    {
      for(k=1; k<K; k++)
	{
	  nmescurr(k) = 0;
	  nmescurry(k) = sum(nobs(span(0,m0+ny(k-1)-1)));
	  m0 += ny(k-1);
	}
    } // obs rangees par dimension, sujet, outcome puis temps


  nmesdcurr=0;
  jseuil1curr=0;
  jseuil2curr=nobs[M];



  vraistot=0.0; 

  vrais_i = 0.0;
  vrais_y = 0.0;
  vrais_d = 0.0;
  
  for(i=0; i<ns; i++)
    {
      vrais_i=0.0;
      ni01 = nmes(i,M); // nombre d'evaluation de la demence
      ni02=0;
      if(nbevt==2) ni02 = nmes(i,M+1);
      ni0 = nmes(i,M+nbevt);
      ni = sum(nmes.row(i))-ni0-ni01-ni02;// ni=sum(nikm) 
      niparK.zeros(K);
      
      // vecteurs et matrices a ni lignes //? faire des maxmes plutot que de redefinir?
      Hyi.zeros(ni); 
      muiK.zeros(ni);
      RiK.zeros(ni,ni); // contient cor, alpha et sigma
      Xi.zeros(ni,beta.size());
      Xcontri.zeros(ni,bcontr.size());
      Zi.zeros(ni,sum(nea));
      tcor.zeros(ni);
      // vecteurs et matrices a ni0 lignes
      muiD.zeros(K*ni0);
      RiD.zeros(K*ni0,K*ni0);// contient cor
      Xdi.zeros(K*ni0,beta.size());
      Zdi.zeros(K*ni0,B.n_cols);
      tdcor.zeros(K*ni0);

      RiKD.zeros(RiK.n_rows,RiD.n_cols);// contient cor


      inodes=0;
      ibtot=0;
      jcurr=0;
      jdcurr=0;
      iidx=0;
      p=0;
      pcontr=0;
      q=0;
      jac=0; 
      ispl=0;
      m0=0;
      kbm=0;
      


      // remplir les vecteurs et matrices outcome par outcome
      for(k=0; k<K; k++)
  	{
	  nik = sum(nmes(i,span(m0,m0+ny(k)-1)));
	  
	  tmpnmes=0;
	  tmpntr=0;


	  //transfos
	  for(m=m0; m<m0+ny(k); m++)
	    {
	      if(nmes(i,m)==0) continue; // !ne devrait pas arriver, gerer dans R!

	      niparK(k) += nmes(i,m);

	      // construire H(Y) et calculer le jacobien
	      yim = Y(span(nmescurry(k)+tmpnmes,nmescurry(k)+tmpnmes+nmes(i,m)-1));
	      Hyim.zeros(nmes(i,m));
	      if(link(m)==0)
		{
		  Hyim = (yim - btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr)) / btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr+1);
		  jac += -nmes(i,m)*log(btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr+1));
		  tmpntr += 2;
		  inodes += 2;	  		  
		}
	      else
		{
		  if(link(m)==1)
		    {
		      aatmp = exp(btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr)) / (1+exp(btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr)));
		      bbtmp = exp(btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr+1)) / (1+exp(btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr+1)));
		      bbtmp = aatmp*(1-aatmp)*bbtmp;
		      cc = btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr+2);
		      cc = fabs(cc);
		      dd = btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr+3);
		      dd = fabs(dd);

		      aa = aatmp*aatmp*(1-aatmp)/bbtmp-aatmp;
		      bb = aa*(1-aatmp)/aatmp;
		  
		      ytemp1.zeros(yim.size());
		      ytemp1 = (yim-nodes(inodes)+epsY(m)) / (nodes(inodes+1)-nodes(inodes)+2*epsY(m));

		      ytemp2 = Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(ytemp1));
		      Hyim=Rcpp::pbeta(ytemp2,aa,bb);

		      ytemp3 = Rcpp::dbeta(ytemp2,aa,bb);

		      endj2=yim.size();
		      for(j2=0;j2<endj2;j2++)
			{
			  yy = ytemp3(j2);//Rcpp::dbeta(ytemp2(j2),aa,bb);
			  jac += std::log(std::abs(yy)/dd);
			  ll = std::abs(nodes(inodes+1)-nodes(inodes)+2*epsY(m));
			  jac -= std::log(ll);
			}
		  
		      tmpntr += 4;
		      inodes += 2;
		    }
		  else
		    {
		      if(link(m)==2)
			{
			  /// prm splines dans splaa
			  splaa.zeros(ntr(m));
			  splaa(0) = btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr);
			  for(kk=1; kk<ntr(m); kk++)
			    {
			      splaa(kk) = btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr+kk)*btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr+kk);
			    }

			  // repeter les noeuds dans zitr
			  zitr.zeros(ntr(m)+2);
			  zitr(0)=nodes(inodes);
			  zitr(1)=zitr(0);
			  zitr(span(2,ntr(m)-1))=nodes(span(inodes,inodes+ntr(m)-3));
			  zitr(ntr(m))=zitr(ntr(m)-1);
			  zitr(ntr(m)+1)=zitr(ntr(m)-1);
			  
			  // calcul de H(y)
			  for(j=0;j<nmes(i,m);j++)
			    {
			      ll=0;
			      if(yim(j)==nodes(inodes+ntr(m)-3))
				{
				  ll=ntr(m)-2;
				}
			    
			      som=0.0;
			      for(kk=3;kk<ntr(m);kk++)
				{
				  if((yim(j)>=zitr(kk-1)) && (yim(j)<zitr(kk)))
				    {
				      ll=kk-1;
				    }
				}
			      
			      if((ll<2) || (ll>ntr(m)-2))
				{
				  vraistot=-1E9;
				  goto fin;
				}

			      som=splaa(0);
			      if(ll>2)
				{
				  som=som+sum(splaa(span(1,ll-2)));
				}


			      ht2 = zitr(ll+1)-yim(j);
			      htm= yim(j)-zitr(ll-1);
			      ht = yim(j)-zitr(ll);
			      ht3 = zitr(ll+2)-yim(j);
			      hht = yim(j)-zitr(ll-2);
			      h1 = zitr(ll+1)-zitr(ll);
			      hh= zitr(ll+1)-zitr(ll-1);
			      hn= zitr(ll+1)-zitr(ll-2);
			      h2n=zitr(ll+2)-zitr(ll-1);
			      h2= zitr(ll+2)-zitr(ll);
			      h3= zitr(ll+3)-zitr(ll);
			    
			      if(yim(j)==zitr(ntr(m)-1))
				{
				  mm2=0.0;
				  mm1=0.0;
				  mm=3.0/h1;
				}
			      else
				{
				  mm2=(3.0*ht2*ht2)/(hh*h1*hn);
				  mm1=(3.0*htm*ht2)/(h2n*hh*h1)+(3.0*ht*ht3)/(h2*h1*h2n);
				  mm=(3.0*ht*ht)/(h3*h2*h1);
				}

			      if((mm<0) || (mm1<0) || (mm2<0))
				{
				  vraistot=-1E9;
				  goto fin;
				}

			      im2 = hht*mm2/(3.0) + h2n*mm1/(3.0) + h3*mm/(3.0);
			      im1 = htm*mm1/(3.0) + h3*mm/(3.0);
			      im = ht*mm/(3.0);


			      Hyim(j) = som + splaa(ll-1)*im2 + splaa(ll)*im1 + splaa(ll+1)*im;

			      jac += log(splaa(ll-1)*mm2+splaa(ll)*mm1+splaa(ll+1)*mm);
			    }

			  inodes += ntr(m)-2;
			  tmpntr += ntr(m);
			}
		    }
		}
	      Hyi(span(jcurr+tmpnmes,jcurr+tmpnmes+nmes(i,m)-1)) = Hyim;
	      
	      tmpnmes += nmes(i,m);
	    }// fin boucle m

	  // toutes les variables du k-ieme MM sont dans l'element k de la liste X
	  Xk = Rcpp::as<arma::mat>(X[k]);
	  Xdk = Rcpp::as<arma::mat>(Xd[k]);
	  for(l=0; l<nv(k); l++)
	    {
	      if(idx(0,iidx+l)==1) // effet fixe
		{
		  //cout << "p=" << p << " l=" << l << endl;
		  Xi(span(jcurr,jcurr+nik-1),p) = Xk(span(nmescurr(k),nmescurr(k)+nik-1),l);
		  Xdi(span(jdcurr,jdcurr+ni0-1),p) = Xdk(span(nmesdcurr,nmesdcurr+ni0-1),l);
		  
		  p += 1;
		}

	      if(idx(1,iidx+l)==1) // contraste
		{
		  tmpnmes=0;
		  ll=0;
		  for(m=m0; m<m0+ny(k); m++)
		    {
		      Xcontri(span(jcurr+tmpnmes,jcurr+tmpnmes+nmes(i,m)-1),ny(k)*pcontr+ll) = Xk(span(nmescurr(k)+tmpnmes,nmescurr(k)+tmpnmes+nmes(i,m)-1),l);
		      ll += 1;
		      tmpnmes += nmes(i,m);
		    }
		  
		  pcontr += 1;
		}
	      
	      
	      if(idx(2,iidx+l)==1) // effet aleatoire
		{
		  Zi(span(jcurr,jcurr+nik-1),q) = Xk(span(nmescurr(k),nmescurr(k)+nik-1),l);

		  Zdi(span(jdcurr,jdcurr+ni0-1),q) = Xdk(span(nmesdcurr,nmesdcurr+ni0-1),l);

		  q += 1;
		}
	      
	      if(idx(3,iidx+l)==1) // temps pour BM ou AR
		{
		  tcor(span(jcurr,jcurr+nik-1)) = Xk(span(nmescurr(k),nmescurr(k)+nik-1),l);

		  tdcor(span(jdcurr,jdcurr+ni0-1)) = Xdk(span(nmesdcurr,nmesdcurr+ni0-1),l);
		}
	      
	    }
	  iidx += nv(k);
	  
	  // construire Rk, Rd et Rkd (correlations BM ou AR )
	  if(ncor(k)>0)
	    {
	      if(ncor(k)==1) // BM
		{
		  for(j1=0; j1<nik; j1++)
		    {
		      for(j2=0; j2<nik; j2++)
			{
			  RiK(jcurr+j1,jcurr+j2) = BM(kbm,kbm)*
			    std::min(tcor(jcurr+j1),tcor(jcurr+j2));
			}
		      
		      
		      for(j3=0; j3<ni0; j3++)
			{
			  RiKD(jcurr+j1,jdcurr+j3) = BM(kbm,kbm)*
			    std::min(tcor(jcurr+j1),tdcor(jdcurr+j3));
			}
		    }

		 
		  for(j1=0; j1<ni0; j1++)
		    {
		      for(j2=0; j2<ni0; j2++)
			{
			  RiD(jdcurr+j1,jdcurr+j2) = BM(kbm,kbm)*
			    std::min(tdcor(jdcurr+j1),tdcor(jdcurr+j2));
			}
		    }

		  // correlation des BM : a refaire
		  if((kbm>0) & (nBM>0))
		    {
		      kk=0;
		      tmpnmes=0;
		      kkbm=0;
		      // dans Rik et Rikd:
		      for(j1=0; j1<jcurr; j1++)
			{
			  if(j1==tmpnmes+niparK(kk))
			    {
			      tmpnmes += niparK(kk);
			      kk += 1;
			      if(ncor(kk)==1) kkbm += 1;
			    }

			  for(j2=jcurr; j2<jcurr+nik; j2++)
			    {
			      RiK(j1,j2) = BM(kkbm,kbm)*std::min(tcor(j1),tcor(j2));
			    }

			  for(j3=jdcurr; j3<jdcurr+ni0; j3++)
			    {
			      RiKD(j1,j3) = BM(kkbm,kbm)*std::min(tcor(j1),tdcor(j3));
			    }
			}

		      
		      kk=0;
		      tmpnmes=0;
		      kkbm=0;
		      // dans RiD:
		      for(j1=0; j1<jdcurr; j1++)
			{
			  if(j1==tmpnmes+ni0)
			    {
			      tmpnmes += ni0;
			      kk += 1;
			      if(ncor(kk)==1) kkbm += 1;
			    }

			  // cov entre 2 BM tdem
			  for(j3=jdcurr; j3<jdcurr+ni0; j3++)
			    {
			      RiD(j1,j3) = BM(kkbm,kbm)*std::min(tdcor(j1),tdcor(j3));
			    }
			}
		      
		    }

		  kbm += 1;

		}
	      else // AR
		{
		  for(j1=0; j1<nik; j1++)
		    {
		      for(j2=0; j2<nik; j2++)
			{ 
			  RiK(jcurr+j1,jcurr+j2) = 
			    btot(ibtot+nef(k)+ncontr(k)+nvc(k)+1)*
			    btot(ibtot+nef(k)+ncontr(k)+nvc(k)+1)*
			    exp(-btot(ibtot+nef(k)+ncontr(k)+nvc(k))*
				fabs(tcor(jcurr+j1)-tcor(jcurr+j2)));
			}
		      
		      
		      for(j3=0; j3<ni0; j3++)
			{
			  RiKD(jcurr+j1,jdcurr+j3) = 
			    btot(ibtot+nef(k)+ncontr(k)+nvc(k)+1)*
			    btot(ibtot+nef(k)+ncontr(k)+nvc(k)+1)*
			    exp(-btot(ibtot+nef(k)+ncontr(k)+nvc(k))*
				fabs(tcor(jcurr+j1)-tdcor(jdcurr+j3)));
			}
		    }


		  for(j1=0; j1<ni0; j1++)
		    {
		      for(j2=0; j2<ni0; j2++)
			{
			  RiD(jdcurr+j1,jdcurr+j2) = 
			    btot(ibtot+nef(k)+ncontr(k)+nvc(k)+1)*
			    btot(ibtot+nef(k)+ncontr(k)+nvc(k)+1)*
			    exp(-btot(ibtot+nef(k)+ncontr(k)+nvc(k))*
				fabs(tdcor(jdcurr+j1)-tdcor(jdcurr+j2)));
			}
		    }
		  
		} // fin AR
	      
	    } // fin ncor


	  // alpha (intercept aleatoire specif au marqueur)
	  if(nalea(k)>0)
	    {
	      ll=0;
	      tmpnmes=0;
	      for(m=m0; m<m0+ny(k); m++)
		{
		  for(j1=0; j1<nmes(i,m); j1++)
		    {
		      for(j2=0; j2<nmes(i,m); j2++)
			{
			  RiK(jcurr+tmpnmes+j1,jcurr+tmpnmes+j2) += 
			    btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+ll)*
			    btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+ll);
			}
		    }
		  tmpnmes += nmes(i,m);
		  ll += 1;
		}
	    }

	  // erreur de mesure
	  ll=0;
	  tmpnmes=0;
	  for(m=m0; m<m0+ny(k); m++)
	    {
	      for(j1=0; j1<nmes(i,m); j1++)
		{
		  RiK(jcurr+tmpnmes+j1,jcurr+tmpnmes+j1) += 
		    btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ll)*
		    btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ll);
		}
	      tmpnmes += nmes(i,m);
	      ll += 1;
	    }
	  
	  nmescurr(k) += nik;
	  jcurr += nik;
	  ibtot += nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+ntrtot(k);
	  nmescurry(k) += nik;
	  jdcurr += ni0;
	  m0 += ny(k);	      
	} // fin boucle k


      // definir les quantites qui interviennent dans f(y)
      muiK = Xi*beta + Xcontri*bcontr;
 
      VK = Zi*B*trans(Zi) + RiK ; //valea+sigma contenu dans Rik

      aa = det(VK);
      if(aa<1.0E-10) 
      	{
	  cout << "det(VK)=" <<  det(VK)  << endl;
	  vraistot=-1E9;
      	  goto fin;
	}
      invVK = inv_sympd(VK);


      // calcul de f(y) en log
      vraistot += as_scalar( -ni*log(2*datum::pi) -log(det(VK)) - 
			     trans(Hyi-muiK)*invVK*(Hyi-muiK))/2 + jac;
 
      vrais_i = as_scalar( -ni*log(2*datum::pi) -log(det(VK)) - 
			   trans(Hyi-muiK)*invVK*(Hyi-muiK))/2 + jac;
      //if(i<10)  cout << "vrais_i=" << vrais_i << endl;

      vrais_y += as_scalar( -ni*log(2*datum::pi) -log(det(VK)) - 
			    trans(Hyi-muiK)*invVK*(Hyi-muiK))/2 + jac;


      // debut evts //
      vrais_d = 0.0;
      
      // definir les quantites qui interviennent dans Phi(K*ni0)
      dimPhi = ni01+ni02;
      upper.zeros(dimPhi);
      lower.zeros(dimPhi);
      delta.zeros(dimPhi);
      nu=0;
      infin.zeros(dimPhi);

      varD.fill(0); // variables explicatives non dep du tps sur D
      varD(0) = 1.0;
      if(nbevt==2) varD(nvarD1+1) = 1.0;
      if(nvarD>0)
	{
	  l=1;
	  for(j=0;j<nvarD;j++)
	    {
	      if(idcause(j)==1)
		{
		  varD(l) = D(i,j+1);
		  l += 1;
		}
	    }
	  if(nbevt==2)
	    {
	      l+=1; // car intercept de D2 (en place nvarD1+1) deja rempli
	      for(j=0;j<nvarD;j++)
		{
		  if(idcause(nvarD+ndept+j)==1)
		    {
		      varD(l) = D(i,j+1);
		      l += 1;
		    }
		}
	    }
	}
      
      Gam1.zeros(ni01,K*ni0); // contribution des dimension sur Dem
      GG.zeros(ni01,ni0);

      l=0;
      for(j=0;j<ni0;j++)
	{
	  if(idD(nmesdcurr+j)==1)
	    {
	      GG(l,j)=1;
	      l+=1;
	    }
	}
      Gam1 = kron(trans(gamma1*varD(span(0,nvarD1))),GG);

      Gam2.zeros(ni02,K*ni0); // contribution des dimension sur Dc
      if(nbevt==2)
	{
	  GG.zeros(ni02,ni0);
	  l=0;
	  for(j=0;j<ni0;j++)
	    {
	      if(idD(nobs[M+2]+nmesdcurr+j)==1)
		{
		  GG(l,j)=1;
		  l+=1;
		}
	    }
	  Gam2 = kron(trans(gamma2*varD(span(nvarD1+1,nvarD1+1+nvarD2))),GG);
	}

      Gam = join_vert(Gam1,Gam2);
      tGam=trans(Gam);
      seuil_i.ones(ni01+ni02); // seuil a chaque temps
      seuil_i(span(0,ni01-1)) *= dot(eta1,varD(span(0,nvarD1)));
      if(nbevt==2) seuil_i(span(ni01,ni01+ni02-1)) *= dot(eta2,varD(span(nvarD1+1,nvarD1+nvarD2+1)));
      if(ndept>0)
	{
	  if(ndept1>0)
	    {
	      for(j=0; j<ni01; j++)
		{
		  seuil_i(j) += dot(eta1_dept,trans(Xseuil(jseuil1curr+j,span(0,ndept1-1))));
		}
	    }
	  if(ndept2>0)
	    {
	      for(j=0; j<ni02; j++)
		{
		  seuil_i(ni01+j) += dot(eta2_dept,trans(Xseuil(jseuil2curr+j,span(ndept1,ndept1+ndept2-1))));
		}     
	    }
	}




      // on conditionne sur Lambda pour calculer la vrais //
      
      muiD = Xdi*beta; //moyenne de Lambda
      VKD = Zi*B*trans(Zdi) + RiKD; //covariance de Y et Lambda
      VD = Zdi*B*trans(Zdi) + RiD;	//variance de Lambda 
      
      muiD_K = muiD + trans(VKD)*invVK*(Hyi-muiK); // moyenne de Lambda sachant Y
      VD_K = VD-trans(VKD)*invVK*VKD; // variance de Lambda sachant Y
      
      upper = seuil_i - Gam*muiD_K; // borne de l'integration pour Phi
      
      V0.zeros(ni01+ni02,ni01+ni02);
      I0.eye(ni01+ni02,ni01+ni02);  
      V0 = I0 + Gam*VD_K*tGam; // variance de Phi
      
      
      lower = upper;
      maxpts = 10000*dimPhi;
      abseps = 0.001;
      releps = 0.0001;
      error =0;
      value0 = 0;
      inform = 0;

      // creation du vecteur des correlations a partir de V0
      correl.zeros((ni01+ni02)*(ni01+ni02+1)/2-ni01-ni02);
      endj2 = ni01+ni02; 
      dimPhi = ni01+ni02;
      
      jexcl1 = -1;
      jexcl2 = -1;
      if(nbTronc(i,0)>0) 
	{
	  jexcl1 = 0;
	  dimPhi -= 1;
	}
      if(nbevt==2)
	{
	  if(nbTronc(i,1)>0) 
	    {
	      jexcl2 = ni01;
	      dimPhi -= 1;
	    }
	}
      
      upper1.zeros(dimPhi);
      kk=0;
      ll=0;
      for(j2=0; j2<endj2; j2++)
	{
	  if(j2==jexcl2) continue;
	  
	  if(j2!=jexcl1)
	    {	      
	      upper1(kk) = upper(j2)/pow(V0(j2,j2),0.5);
	      kk += 1;
	    }
	  
	  if(j2>0)
	    {
	      for(j1=0; j1<j2; j1++)
		{
		  if(j1==jexcl1) continue;
		  correl(ll) = V0(j1,j2)/(pow(V0(j1,j1),0.5)*pow(V0(j2,j2),0.5));
		  ll += 1;
		}
	    }
	}
      
      lower2 = conv_to< std::vector<double> >::from(lower);
      upper2 = conv_to< std::vector<double> >::from(upper1);
      infin2 = conv_to< std::vector<int> >::from(infin);
      correl2 = conv_to< std::vector<double> >::from(correl);
      delta2 = conv_to< std::vector<double> >::from(delta);

      // calcul de Phi(K*ni0) 
      sadmvn_(&dimPhi,&lower2[0],&upper2[0],&infin[0],&correl2[0],&maxpts,&abseps,&releps,&error,&value0,&inform); 


      if(inform!=0)  
	{
	  cout << "probleme integration 0 :" << "value=" << value0 << " error=" << error << endl;

	  vraistot=-1E9;
	  goto fin;
	}


      if(D(i,0)==0)
	{
	  vrais_d += log(value0);
	  vraistot += log(value0);
	}

      //Phi(K*(ni0-1)) calcule dans le cas Dem
      if(D(i,0)==1)
      	{ 
	  dimPhi -= 1;
	  if(dimPhi>0)
	    {
	      maxpts=10000*dimPhi;


	      // matrice pr enlever la ligne correspondant au tps de Dem
	      P = join_vert(join_horiz(eye<mat>(ni01-1,ni01-1),
				       zeros<mat>(ni01-1,ni02+1)),
			    join_horiz(zeros<mat>(ni02,ni01),eye<mat>(ni02,ni02)));
	      


	      un.ones(ni0-1);
	      upper = P*(seuil_i - Gam*muiD_K);
	      
	      I1.eye(ni01+ni02-1,ni01+ni02-1);  
	      V1.zeros(ni01+ni02-1,ni01+ni02-1);
	      V1 = I1 + P*Gam*VD_K*tGam*trans(P); 
	      	      
	      //correlations
	      jexcl2 -= 1;

	      correl.zeros();
	      endj2 = V1.n_cols;
	      upper1.zeros(dimPhi);
	      kk=0;
	      ll=0;
	      for(j2=0; j2<endj2; j2++)
		{
		  if(j2==jexcl2) continue;
		  
		  if(j2!=jexcl1)
		    {	      
		      upper1(kk) = upper(j2)/pow(V1(j2,j2),0.5);
		      kk += 1;
		    }
		  
		  if(j2>0)
		    {
		      for(j1=0; j1<j2; j1++)
			{
			  if(j1==jexcl1) continue;
			  correl(ll) = V1(j1,j2)/(pow(V1(j1,j1),0.5)*pow(V1(j2,j2),0.5));
			  ll += 1;
			}
		    }
		}

	      
	      upper2 = conv_to< std::vector<double> >::from(upper1);
	      correl2 = conv_to< std::vector<double> >::from(correl);
	      
	      value1 = 0;

	      // calcul de Phi Dem
	      sadmvn_(&dimPhi,&lower2[0],&upper2[0],&infin[0],&correl2[0],&maxpts,&abseps,&releps,&error,&value1,&inform);

	      
 	      if(inform!=0)      
		{
		  cout << "probleme integration 1 :" << inform << " error=" << error << endl;
		  vraistot=-1E9;
		  goto fin;
		}
	    }
	  else
	    {
	      value1 = 1;
	    }

	  vrais_d += log(value1-value0);
	  vraistot += log(value1-value0);   
	  
	}

      // cas Dc
      if(D(i,0)==2)
	{
	  dimPhi -= 1;
	  if(dimPhi>0)
	    {
	      maxpts=10000*dimPhi;

	      // matrice pr enlever la ligne correspondant au tps de Dc
	      P = join_horiz(eye<mat>(ni01+ni02-1,ni01+ni02-1),zeros<vec>(ni01+ni02-1));
	  
	      upper = P*(seuil_i - Gam*muiD_K);
	      
	      V1.zeros(ni01+ni02-1,ni01+ni02-1);
	      I1.eye(ni01+ni02-1,ni01+ni02-1);  
	      V1 = I1 + P*Gam*VD_K*tGam*trans(P);
	      
	      correl.zeros();
	      endj2 = V1.n_cols;
	      
	      upper1.zeros(dimPhi);
	      kk=0;
	      ll=0;
	      for(j2=0; j2<endj2; j2++)
		{
		  if(j2==jexcl2) continue;
		  
		  if(j2!=jexcl1)
		    {	      
		      upper1(kk) = upper(j2)/pow(V1(j2,j2),0.5);
		      kk += 1;
		    }
		  
		  if(j2>0)
		    {
		      for(j1=0; j1<j2; j1++)
			{
			  if(j1==jexcl1) continue;
			  correl(ll) = V1(j1,j2)/(pow(V1(j1,j1),0.5)*pow(V1(j2,j2),0.5));
			  ll += 1;
			}
		    }
		}

	      upper2 = conv_to< std::vector<double> >::from(upper1);
	      correl2 = conv_to< std::vector<double> >::from(correl);
	      
	      value1 = 0;

	      sadmvn_(&dimPhi,&lower2[0],&upper2[0],&infin[0],&correl2[0],&maxpts,&abseps,&releps,&error,&value1,&inform);
		  
	      if(inform!=0)      
		{
		  cout << "probleme integration 2 :" << inform << " error=" << error << endl;
		  vraistot=-1E9;
		  goto fin;
		}
	      vrais_d += log(value1-value0);
	      vraistot += log(value1-value0); 
	    }
	  else
	    { // si 1 seule mesure de D
	      value1 = 1;
	      vrais_d += log(value1-value0);
	      vraistot += log(value1-value0); 
	    }
	}

      
      loglik(i) = vrais_d;

      // entree retardee
      if(any(entreRetard==1))
	{
	  retardi = 0.0;

	  if(nbevt==1)
	    {
	      dimPhi=nbTronc(i,0);
	      upper.zeros(nbevt);
	      V0.zeros(nbevt,nbevt);
	      
	      upper = seuil_i(0) - Gam1.row(0)*muiD;
	      V0 = eye<mat>(nbevt,nbevt) + Gam1.row(0) * ( Zdi*B*trans(Zdi) + RiD ) * trans(Gam1.row(0));
	      
	      upper(0) = upper(0)/pow(V0(0,0),0.5);
	      correl.zeros();
	    }
	  else
	    { // 2 evts
	      nTronc1 = nbTronc(i,0);
	      nTronc2 = nbTronc(i,1);
	      dimPhi = nTronc1+nTronc2;

	      js.zeros(dimPhi);
	      GamTronc.zeros(dimPhi,K*ni0);

	      if((nTronc1>0) & (nTronc2>0))
		{
		  GamTronc = join_vert(Gam1.row(0), Gam2.row(0));
		}
	      else
		{
		  if(nTronc1>0) GamTronc = Gam1.row(0);
		  if(nTronc2>0) GamTronc = Gam2.row(0);
		}

	      if(nTronc1>0) js(0) = 0;
	      if(nTronc2>0) js(nTronc1) = ni01;
	      
	      upper.zeros(dimPhi);
	      V0.zeros(dimPhi,dimPhi);

	      if(dimPhi>0)
		{
		  upper = seuil_i(js) - GamTronc*muiD;
		  V0 = eye<mat>(dimPhi,dimPhi) + GamTronc * ( Zdi*B*trans(Zdi) + RiD ) * trans(GamTronc);
	      
		  correl.zeros();
		  endj2 = V0.n_cols;
		  for(j2=0; j2<endj2; j2++)
		    {
		      if(j2>0)
			{
			  for(j1=0; j1<j2; j1++)
			    {
			      correl(j1+((j2-1)*j2)/2) = V0(j1,j2)/(pow(V0(j1,j1),0.5)*pow(V0(j2,j2),0.5));
			    }
			}
		      
		      upper(j2) = upper(j2)/pow(V0(j2,j2),0.5);
		    } 
		}
	    }

	  lower=upper;		  
	  
	  upper2 = conv_to< std::vector<double> >::from(upper);
	  lower2 = conv_to< std::vector<double> >::from(lower);
	  correl2 = conv_to< std::vector<double> >::from(correl);

	  if(dimPhi>0)
	    {
	      sadmvn_(&dimPhi,&lower2[0],&upper2[0],&infin[0],&correl2[0],&maxpts,&abseps,&releps,&error,&retardi,&inform);
	      
	      if(inform!=0)      
		{
		  cout << "probleme integration retard :" << inform << " error=" << error << endl;
		  vraistot=-1E9;
		  goto fin;
		}
	      
	      vraistot -= log(retardi);
	      loglik(i) -= log(retardi);
	    }
	  
	} // fin entree retardee

      loglik(ns+i) = loglik(i) + vrais_i;
      
      jseuil1curr += ni01;
      jseuil2curr += ni02;
      nmesdcurr += ni0;
    } // fin boucle i
  
 fin:

  return Rcpp::wrap(loglik);
  

}



