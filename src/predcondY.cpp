#include <math.h>
#include <RcppArmadillo.h>


// a faire ajouter cov BM

RcppExport SEXP predcondY(SEXP b0, SEXP bfix0, SEXP fix0, SEXP Y0, SEXP X0, SEXP nmes0, SEXP nv0, SEXP idx0, SEXP idiag0, SEXP ncor0, SEXP ny0, SEXP nalea0, SEXP ntr0, SEXP link0, SEXP nodes0, SEXP nRE0, SEXP nBM0, SEXP chol0)
{
  using namespace arma;

  Rcpp::NumericVector b00(b0);
  Rcpp::NumericVector bfix00(bfix0);
  Rcpp::IntegerVector fix00(fix0);
  Rcpp::NumericVector Y00(Y0);
  Rcpp::List X(X0);
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
  int nRE = Rcpp::as<int>(nRE0);
  int nBM= Rcpp::as<int>(nBM0);
  int chol = Rcpp::as<int>(chol0);

  int ns = nmes00.nrow();
  int K = ny00.size();
  int nvtot = idx00.ncol();
  int M = nmes00.ncol()-1;
  
  vec b(b00.begin(),b00.size(),false);
  vec bfix(bfix00.begin(),bfix00.size());
  ivec fix(fix00.begin(),fix00.size());
  imat nmes(nmes00.begin(),ns,M+1);
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

  //btot.print("btot");
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

  int j,k,l,p,q,ll,kk,m;
  int jcurr,iidx,iBM;

  int j1,j2,j3;
  int endj2,tmpntr,pcontr;

  double sig1,sig2,corr,aa;

  double ht,htm,ht2,ht3,h,hh,h2,h3,h2n,hn,hht,mm,mm1,mm2,im,im1,im2,som;


  int l0 =0;
  int m0 =0;
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

      l0 += nv(k);
      m0 += ny(k);
      subnpm(k) = nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+ntrtot(k);
    }

  // cout << "btot ok " << btot << endl;

  vec splaa(max(ntrtot)*K,fill::zeros);

  int npmMM = sum(subnpm);


  ivec nobs = sum(nmes,0).t(); // nb d'obs por chaque Y



  vec muK, muD;
  vec beta(sum(nef)+K);
  vec bcontr(sum(ncontr)+sum(nvarcontr));
  mat B = zeros(sum(nea),sum(nea));
  mat U = zeros(sum(nea),sum(nea));

  //  mat Sigma = zeros(sum(ny),sum(ny));
  //  mat valea = zeros(sum(ny*nalea));

  int ibeta=0;
  int ibcontr=0;
  int imatB=0;
  int ibtot=0;
  int iRE=0;
  int ispl=0;
  int tmpcontr=0;
  //  int isigma=0;
  //  int ivalea=0;



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
      
	ibeta += nef(k)+1;
	ibtot += nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+ntrtot(k);
	imatB += nea(k);
	//isigma += ny(k);
	//ivalea += nalea(k);
    }


  if(chol==1)
    {
      B = U.t()*U;
    }
  else
    {
      endj2=U.n_cols;
      for(j2=0;j1<endj2;j2++)
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
    }

  //U.print("U=");
  //B.print("B=");	  


  //cout << "B beta et Sigma ok" << endl;
  // nmescurr pour  X
  // nmesdcurr pour Xd
  // nmescurry pour Y
  ivec nmescurr(K);
  ivec nmescurry(K);
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

  ivec niparK(K);

  mat Rk,V;

  // output
  vec predHY = zeros(Y.size());


  int inodes,tmpnmes,ni,nik,ni0;
  //cout << "avt boucle i" << endl;
  for(i=0; i<ns; i++)
    {
      ni = sum(nmes.row(i)); // ni=sum(nikm)
      niparK.fill(0);


      // vecteurs et matrices a ni lignes //? faire des maxmes plutot que de redefinir?
      vec Hyi(ni); 
      vec muiK = zeros(ni);
      mat RiK = zeros(ni,ni); // contient cor et alpha 
      mat SigmaiK = zeros(ni,ni); // erreurs de mesure
      mat Xi = zeros(ni,beta.size());
      mat Xcontri = zeros(ni,bcontr.size());
      mat Zi = zeros(ni,sum(nea));
      vec tcor = zeros(ni);


      inodes=0;
      ibtot=0;
      jcurr=0;
      iidx=0;
      p=0;
      pcontr=0;
      q=0;
      ispl=0;
      m0=0;
      iBM=0;
      


      // remplir les vecteurs et matrices outcome par outcome
      for(k=0; k<K; k++)
  	{	  
  	  //cout << "i=" << i << "  k=" << k  <<endl;
	  //cout << "nmescurr=" << nmescurr(k)  <<endl;
	  
	  nik = sum(nmes(i,span(m0,m0+ny(k)-1)));
	  //cout << "m0=" << m0 << "  ny=" << ny(k) << "nmes=" << nmes.row(i) <<endl;
	  

	  tmpnmes=0;
	  tmpntr=0;


	  //transfos
	  for(m=m0; m<m0+ny(k); m++)
	    {
	      if(nmes(i,m)==0)
		{
		  if(link(m)==0)
		    {
		      tmpntr += 2;
		      inodes += 2;
		    }
		  
		  if(link(m)==1)
		    {
		      tmpntr += 4;
		      inodes += 2;
		    }
		  
		  if(link(m)==2)
		    {
		      tmpntr += ntr(m);
		      inodes += ntr(m)-2;
		    }
		  
		  continue; 
		}

	      niparK(k) += nmes(i,m);

	      // construire H(Y) et calculer le jacobien
	      //cout << "avt yim" << endl;
	      vec yim = Y(span(nmescurry(k)+tmpnmes,nmescurry(k)+tmpnmes+nmes(i,m)-1));
	      //cout << "nmescurry=" << nmescurry(k) << " tmpnmes=" << tmpnmes << " nmes=" << nmes(i,m) << endl;
	      vec Hyim(nmes(i,m));
	      //cout << "yim " << yim  << endl;
	      //if(i<3) {cout << "yim " << yim  << endl;}
	      //cout << "yim ok" << endl;
	      if(link(m)==0)
		{
		  Hyim = (yim - btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr)) / btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr+1);
		  
		  tmpntr += 2;
		  inodes += 2;	  
		  //cout << "Hyim " << Hyim  << endl;
		  
		}
	      else
		{
		  if(link(m)==1)
		    {
		      // pas fait
		    }
		  else
		    {
		      if(link(m)==2)
			{
			  /// prm splines dans splaa
			  vec splaa(ntr(m));
			  splaa(0) = btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr);
			  for(kk=1; kk<ntr(m); kk++)
			    {
			      splaa(kk) = btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr+kk)*btot(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr+kk);
			    }

			  // repeter les noeuds dans zitr
			  vec zitr(ntr(m)+2);
			  zitr(0)=nodes(inodes);
			  zitr(1)=zitr(0);
			  zitr(span(2,ntr(m)-1))=nodes(span(inodes,inodes+ntr(m)-3));
			  zitr(ntr(m))=zitr(ntr(m)-1);
			  zitr(ntr(m)+1)=zitr(ntr(m)-1);
			  // cout << "zitr=" << zitr << endl;
			  // calcul de H(y)
			  for(j=0;j<nmes(i,m);j++)
			    {
			      ll=0;
			      if(yim(j)==nodes(inodes+ntr(m)-3))
				{
				  ll=ntr(m)-2;
				}
			      // cout << "y=" << yim(j) << endl;
			      som=0.0;
			      for(kk=3;kk<ntr(m);kk++)
				{
				  if((yim(j)>=zitr(kk-1)) && (yim(j)<zitr(kk)))
				    {
				      ll=kk-1;
				    }
				}

			      //cout << "ll=" << ll << endl;

			      if((ll<2) || (ll>ntr(m)-2))
				{
				  //vraistot=-1E9;
				  goto fin;
				}

			      som=splaa(0);
			      if(ll>2)
				{
				  som=som+sum(splaa(span(1,ll-2)));
				}

			    
			      // cout << "splaa=" << splaa << endl; 
			      //cout << "som=" << som << endl;

			      ht2 = zitr(ll+1)-yim(j);
			      htm= yim(j)-zitr(ll-1);
			      ht = yim(j)-zitr(ll);
			      ht3 = zitr(ll+2)-yim(j);
			      hht = yim(j)-zitr(ll-2);
			      h = zitr(ll+1)-zitr(ll);
			      hh= zitr(ll+1)-zitr(ll-1);
			      hn= zitr(ll+1)-zitr(ll-2);
			      h2n=zitr(ll+2)-zitr(ll-1);
			      h2= zitr(ll+2)-zitr(ll);
			      h3= zitr(ll+3)-zitr(ll);
			    
			      if(yim(j)==zitr(ntr(m)-1))
				{
				  mm2=0.0;
				  mm1=0.0;
				  mm=3.0/h;
				}
			      else
				{
				  mm2=(3.0*ht2*ht2)/(hh*h*hn);
				  mm1=(3.0*htm*ht2)/(h2n*hh*h)+(3.0*ht*ht3)/(h2*h*h2n);
				  mm=(3.0*ht*ht)/(h3*h2*h);
				}

			      if((mm<0) || (mm1<0) || (mm2<0))
				{
				  //vraistot=-1E9;
				  goto fin;
				}

			      im2 = hht*mm2/(3.0) + h2n*mm1/(3.0) + h3*mm/(3.0);
			      im1 = htm*mm1/(3.0) + h3*mm/(3.0);
			      im = ht*mm/(3.0);


			      Hyim(j) = som + splaa(ll-1)*im2 + splaa(ll)*im1 + splaa(ll+1)*im;

			      //cout << "y=" << yim(j) << endl;  
			      //cout << "h(y)=" << Hyim(j) << endl;
			      //cout << "jac=" << log(splaa(ispl+ll-1)*mm2+splaa(ispl+ll)*mm1+splaa(ispl+ll+1)*mm) << endl;
			    }

			  inodes += ntr(m)-2;
			  tmpntr += ntr(m);
			  //cout << "jac=" << jac << endl;
			}
		    }
		}
	      //cout << "avt Hyi" << endl;
	      Hyi(span(jcurr+tmpnmes,jcurr+tmpnmes+nmes(i,m)-1)) = Hyim;
	      //cout << "jcurr=" << jcurr << " size Hyi=" << Hyi.size() << endl;
	      
	      tmpnmes += nmes(i,m);
	    }// fin boucle m

	  

	  // toutes les variables du k-ieme MM sont dans l'element k de la liste X
	  mat Xk = X[k];
	  
	  //cout << "nmes(i,)"<< nmes.row(i) << endl;
	  // creer X, Z et tcor a partir de Xk et idx
	  // creer Xd, Zd et tdcor a partir de Xdk et idx
	  //cout << "avt creation Xet Xd" << endl;
	  for(l=0; l<nv(k); l++)
	    {
	      if(idx(0,iidx+l)==1) // effet fixe
		{
		  Xi(span(jcurr,jcurr+nik-1),p) = Xk(span(nmescurr(k),nmescurr(k)+nik-1),l);
		  p += 1;
		  //cout << "p= "<< p << endl;
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
		  //cout << "pcontr= "<< pcontr << endl;
		}
	      
	      
	      if(idx(2,iidx+l)==1) // effet aleatoire
		{
		  Zi(span(jcurr,jcurr+nik-1),q) = Xk(span(nmescurr(k),nmescurr(k)+nik-1),l);

		  q += 1;
		  //cout << "q= "<< q << endl;
		}
	      
	      if(idx(3,iidx+l)==1) // temps pour BM ou AR
		{
		  tcor(span(jcurr,jcurr+nik-1)) = Xk(span(nmescurr(k),nmescurr(k)+nik-1),l);

		  //cout << "tcor ok "<< endl;
		}
	      
	    }
	  iidx += nv(k);
	  // cout << "X et Z ok" << endl;
	  
	  // construire Rk (correlations BM ou AR )
	  if(ncor(k)>0)
	    {
	      //cout << "cor="<< btot(ibtot+nef(k)+nvc(k)+ntrtot(k))<< endl;
	      //tcor.print("tcor=");

	      if(ncor(k)==1) // BM
		{
		  for(j1=0; j1<nik; j1++)
		    {
		      for(j2=0; j2<nik; j2++)
			{
			  RiK(jcurr+j1,jcurr+j2) = 
			    btot(ibtot+nef(k)+ncontr(k)+nvc(k))*
			    btot(ibtot+nef(k)+ncontr(k)+nvc(k))*
			    std::min(tcor(jcurr+j1),tcor(jcurr+j2));
			}
		    }

		  // correlation des BM
		  if((k>0) & (nBM>0))
		    {
		      kk=0;
		      tmpnmes=0;
		      iBM=(k-1)*k/2;

		      for(j1=0; j1<jcurr; j1++)
			{
			  if(j1==tmpnmes+niparK(kk))
			    {
			      tmpnmes += niparK(kk);
			      kk += 1;
			      iBM += 1;
			    }

			  for(j2=jcurr; j2<jcurr+nik; j2++)
			    {
			      RiK(j1,j2) = btot(npmMM+nRE+iBM)*std::min(tcor(j1),tcor(j2));
			    }
			}
		    }
		  // cout << "R ok" << endl;
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
		  SigmaiK(jcurr+tmpnmes+j1,jcurr+tmpnmes+j1) += 
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
	  m0 += ny(k);
	   		      
	} // fin boucle k


      // definir les quantites qui interviennent dans L1(y)=f(y)
      muiK = Xi*beta + Xcontri*bcontr;

 
      mat VK = Zi*B*trans(Zi) + RiK + SigmaiK ; //valea contenu dans Rik


      aa = det(VK);
      if(aa<1.0E-10) 
      	{
	  cout << "i=" << i << " ni=" << ni << endl;
	  cout << "VK=" << VK << endl;
	  Zi.print("Zi=");
	  B.print("B=");
	  //I.print("I=");	  
      	  VK.print("VK=");
      	  //sv.print("pour VK eig=");
	  cout << "det(VK)=" <<  det(VK)  << endl;
	  btot.print("btot=");
	  //vraistot=-1E9;
      	  goto fin;
      	}
      mat invVK = inv_sympd(VK);

      mat cov1 = Zi*B*trans(Zi) + RiK;

      vec predHYi = muiK + cov1*invVK*(Hyi-muiK);

      // dans le vecteur resultat
      jcurr=0;
      for(k=0;k<K;k++)
      	{
	  predHY(span(nmescurry(k)-niparK(k),nmescurry(k)-1)) = predHYi(span(jcurr,jcurr+niparK(k)-1));
	  
	  jcurr += niparK(k);	      
      	}

    } // fin boucle i

 fin:
  return Rcpp::wrap(predHY);
  

}
