#include <math.h>
#include <RcppArmadillo.h>

extern "C"{    // fct permettant de calculer la fct de repartition d'une loi multivariee gaussienne
  void sadmvn_(int *n, double *lower, double *upper,
               int *infin, double *correl,
               int *maxpts, double *abseps, double *releps,
               double *error, double *value, int *inform );
};

RcppExport SEXP probaEC(SEXP Bs0, SEXP Bst0, SEXP b0, SEXP Y0, SEXP X0, SEXP Xd0, SEXP D0, SEXP Xseuil0, SEXP nmes0, SEXP nv0, SEXP idx0, SEXP idiag0, SEXP ncor0, SEXP ny0, SEXP nalea0, SEXP ntr0, SEXP link0, SEXP nodes0, SEXP nRE0, SEXP nBM0) 
{ 
  
  using namespace arma;
 
  
  /////////////////////////////////////////////
  // initialisation des parametres de la fct //
  /////////////////////////////////////////////
  
  // Premiere initialisation : Rcpp
  double Bs = Rcpp::as<double>(Bs0); // identifiant de l'intervalle juste avant s
  double Bst = Rcpp::as<double>(Bst0); // identifiant de l'intervalle juste avant s+t
  Rcpp::NumericVector b00(b0);
  Rcpp::NumericVector Y00(Y0); 
  Rcpp::List X(X0); 
  Rcpp::List Xd(Xd0);
  Rcpp::NumericVector D00(D0); 
  Rcpp::NumericMatrix Xseuil00(Xseuil0); 
  Rcpp::IntegerVector nmes00(nmes0);
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
  int nBM = Rcpp::as<int>(nBM0); // jamais utilise


  // Deuxieme initialisation : Armadillo
  int K = ny00.size(); // nbre de dimensions
  int nvtot = idx00.ncol(); // somme des nv = nbre total de variables explicatives ds ttes les dimensions
  int nvarD = D00.size()-1; // nbre de variables explicatives independantes du temps
  int ndept = Xseuil00.ncol(); // nbre de variables explicatives dependantes du temps
  int M = nmes00.size(); // nbre de marqueurs
  vec b(b00.begin(),b00.size(),false);
  vec Y(Y00.begin(),Y00.size());
  // X et Xd =listes _ pas a reinitialiser
  vec D(D00.begin(),nvarD+1); 
  mat Xseuil(Xseuil00.begin(),Xseuil00.nrow(),Xseuil00.ncol()); 
  ivec nmes(nmes00.begin(),M); // ivec = vecteur avec des entiers = (contenu,taille)
  ivec nv(nv00.begin(),nv00.size()); 
  imat idx(idx00.begin(),4,nvtot);  // imat = matrice avec des entiers (0 ou 1) = (contenu,nb_lignes,nb_colonnes)
  ivec idiag(idiag00.begin(),idiag00.size()); 
  ivec ncor(ncor00.begin(),ncor00.size()); 
  ivec ny(ny00.begin(),ny00.size()); 
  ivec nalea(nalea00.begin(),nalea00.size()); 
  ivec ntr(ntr00.begin(),ntr00.size());
  ivec link(link00.begin(),link00.size()); 
  vec nodes(nodes00.begin(),nodes00.size()); 
  int npm = b.size(); // nbre total de parametres du modele conjoint

  
  ///////////////////////////////////////////////////////////////
  // initialisation de la variable qui contiendra notre sortie //
  ///////////////////////////////////////////////////////////////
  
  double proba;
  proba = 0;
  
  //////////////////////////////////////
  // initialisation de nvls variables //
  //////////////////////////////////////
  
  
  ivec nef(K); // tous des vecteurs d'entiers de taille K (=nbre de dimensions)
  ivec ncontr(K); 
  ivec nvarcontr(K); 
  ivec nea(K); 
  ivec nvc(K); 
  ivec ntrtot(K); 
  ivec subnpm(K); 
  
  nef.fill(0); // initialises avec que des 0
  ncontr.fill(0); 
  nvarcontr.fill(0); 
  nea.fill(0); 
  nvc.fill(0); 
  ntrtot.fill(0); 
  subnpm.fill(0); 
  
  
  int j,k,l,p,q,ll,kk,m,l0,m0,j1,j2,j3,kbm;   // pr les boucles
  int jcurr,jdcurrs,jdcurrst,iidx,tmpntr,pcontr,ibeta,ibcontr,imatB,imatBM,ibtot,iRE,iBM,tmpcontr,tmpnmes,inodes; // reperes
  int ncortotBM,npmMM,ni,nik; // nombres/comptages
  ivec nmescurry;
  // construction de beta, eta, gamma
  vec beta,bcontr,eta,eta_dept,varD; 
  mat gamma;
  // construction de B et BM
  mat B,BM; 
  // construction de H(Y^(s))
  vec splaa,zitr;
  double aa,ht,htm,ht2,ht3,h1,hh,h2,h3,h2n,hn,hht,mm,mm1,mm2,im,im1,im2,som;
  vec Hyi,yim,Hyim;
  // construction de X^(s),~X,~X+,Z^(s),~Z,~Z+
  mat Xk,Xdk,Xi,Xcontri,Xdis,Xdist,Zi,Zdis,Zdist;
  // construction de R+A+E,~R,~R+,~~R,~~R+
  vec tcor,tdcors,tdcorst;
  mat RiK,RiDs,RiDst,RiKDs,RiKDst;
  // construction de moyenne et variance de H(Y^(s))
  vec muiK; 
  mat VK,invVK;
  // construction de seuils, Gamma et Gamma+
  vec seuil,seuilPlus; 
  mat Gam,GG,GamPlus,GGPlus,tGam,tGamPlus; 
  // construction de moyennes, variances et covariance de H(Y^(s)) et Lambda/Lambda+
  vec muiDs,muiDst,muiD_Ks,muiD_Kst; 
  mat VKDs,VKDst,VDs,VDst,VD_Ks,VD_Kst; 
  // calculs des fcts de repartition
  int dimPhi1,dimPhi2,maxpts1,inform1,maxpts2,inform2;
  int endj1,endj2; 
  vec upper1,lower1,upper11,upper2,lower2,upper22;
  ivec infin1,infin2; 
  mat V1,V2,I1,I2;
  vec correl1,correl2;
  double abseps,releps,error1,error2,value1,value2; 
  std::vector<double> lower111,upper111,correl111,lower222,upper222,correl222;  // std:: namespace standard, type d'elements reclame par la fonction sadmvn
  std::vector<int> infin111,infin222; 

  
  ////////////////////////////////////////////////////////
  // calculs des nbres de parametres dans chq dimension //   
  ////////////////////////////////////////////////////////
  

  l0=0; 
  m0=0; 
  ncortotBM=0; 
  
  
  for(k=0; k<K; k++){ // pr chq dimension
    
    
    nef(k) = -1; // pour ne pas compter l'intercept 
    
    for(l=l0; l<l0+nv(k); l++){ // pr chq variable explicative de la dimension k
      
      
      if(idx(0,l)==1){ // si EF = 1
        nef(k) += 1;  // alors nef = nef + 1  // nef(k) = nbre d'EF de la dimension k
      } 
      if(idx(1,l)==1){ // si contrast = 1
        nvarcontr(k) += 1; // alors nvarcontr  = nvarcontr + 1  // nvarcontr(k) = nbre de variables explicatives avec contraste ds la dimension k
        ncontr(k) += ny(k)-1;  // alors ncontr = ncontr + nbre de marqueurs - 1 (pr identifiabilite) // ncontr(k) = nbre total de contrastes de la dimension k
      } 
      if(idx(2,l)==1){ // si EA = 1
        nea(k) += 1;	 // alors nea = nea + 1  // nea(k) = nbre d'EA de la dimension k
      } 
    } 
    
    if(idiag(k)==1){ // si les EA de la dimension k sont independants
      nvc(k)=nea(k)-1;  // alors nvc = nea - 1 
    } 
    else{   // sinon (si les EA de la dimension k sont correles)
      nvc(k) = nea(k)*(nea(k)+1)/2-1;  // nvc(k) = nbre d'elements de la matrice de var-cov des EA 
    } 
    
    for(m=m0; m<m0+ny(k); m++){ // pr chq marqueur de la dimension k
      

      if(link(m)==0){ // si la fct lien est lineaire
        ntrtot(k) += 2; // alors le nbre total de parametres des fonctions liens, des marqueurs de la dimension k, augmente de 2 parametres 
      } 
      if(link(m)==1){ // si la fct lien est beta
        ntrtot(k) += 4; // alors ntrtot(k) augmente de 4 parametres
      } 
      if(link(m)==2){ // si la fct lien est splines
        ntrtot(k) += ntr(m); // alors ntrtot(k) augmente du nbre de parametres de la fct lien du marqueur m
      } 
    }  // ntrtot(k) = nbre total de parametres des fonctions de lien de ts les marqueurs de la dimension k
    

    if(ncor(k)==1){  // si la dimension k contient un mvmt brownien
      ncortotBM += 1;  // alors ncortotBM augmente de 1  // ncortotBM = nbre total de mvmt brownien <= K
    } 
    

    l0 += nv(k);  // l0 = repere des variables explicatives de toutes les dimensions
    m0 += ny(k);  // m0 = repere des marqueurs de toutes les dimensions
    subnpm(k) = nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+ntrtot(k); // nbre de parametres de la dimension k  //le ny(k) represente les erreurs de mesure

  } 
  

  splaa.zeros(max(ntrtot)*K); // splaa = vecteur, de taille = nbre maxi de parametres de fcts de lien d'une dimension, rempli de 0
  
  npmMM = sum(subnpm); // nbre total de parametres (intra-dim), ttes dimensions confondues


  ///////////////////////////////////////////////////////
  // les variables explicatives jouant sur l'evenement //   
  ///////////////////////////////////////////////////////
  
  // nombre de variables explicatives par evenement 
  //nvarD // nbre de variables explicatives independantes du temps pr event
  //ndept // nbre de variables explicatives dependantes du temps pr event
  
  beta.zeros(sum(nef)+K);  // beta = vecteur, de taille = nbre total d'EF + nbre de dimensions (pr intercept), compose que de 0
  bcontr.zeros(sum(ncontr)+sum(nvarcontr)); 
  
  eta.zeros(nvarD+1); // eta, VarExpl indep + 'intercept'
  eta_dept.zeros(ndept); // eta_dept, VarExpl dependantes du temps
  
  gamma.zeros(K,nvarD+1); // matrice
  
  varD.zeros(nvarD+1); 
  B.zeros(sum(nea),sum(nea)); 
  BM.zeros(ncortotBM,ncortotBM); 
  
  ibeta=0; // initialisation
  ibcontr=0; 
  imatB=0; 
  imatBM=0; 
  ibtot=0; 
  iRE=0; 
  iBM=0; 
  tmpcontr=0; 
  
  
  
  //////////////////////////////////////////////////////////////////
  // creation des vecteurs beta et bcontr, ainsi que la matrice B //  
  //////////////////////////////////////////////////////////////////  
  
  // on recupere les parametres du vecteur b (btot auparavant)
  
  for(k=0; k<K; k++){  // pr chq dimension
    // vecteur beta = effets fixes des K dimensions
    beta(ibeta) = 0; // 1er element = 0
    beta(span(ibeta+1,ibeta+nef(k))) = b(span(ibtot,ibtot+nef(k)-1)); // span(first,last)  // beta=(0,k1 1er elements de b,0,...) ou k1=nef(1)
    
    // vecteur bcontr = contrastes  des K dimensions
    for(l=0; l<nvarcontr(k); l++){  // pr chq variable explicative avec contraste
      bcontr(span(ibcontr,ibcontr+ny(k)-2)) = b(span(ibtot+nef(k)+tmpcontr,ibtot+nef(k)+tmpcontr+ny(k)-2)); 
      bcontr(ibcontr+ny(k)-1) = -sum(b(span(ibtot+nef(k)+tmpcontr,ibtot+nef(k)+tmpcontr+ny(k)-2))); // identifiablite : sum(contrastes)=0
      
      tmpcontr += ny(k)-1; 
      ibcontr += ny(k); 
    } 
    
    
    // varcov des EA du sous-modele (= de la dimension) k  
    B(imatB,imatB) = 1; // 1 element de la matrice, en haut a gauche   // matrice U permet de retrouver les correlations des EA
    if(idiag(k)==1){ // si les EA sont independants
      for(j=1; j<nea(k); j++){ 
        B(imatB+j,imatB+j) = b(ibtot+nef(k)+ncontr(k)+j-1); // diagonale de la matrice = parametres de b correspondants
      } 
    } 
    else{ // sinon, si les EA sont correles
      j=0; 
      for(j2=0; j2<nea(k); j2++){ 
        if(j2>0){ 
          for(j1=0; j1<=j2; j1++){ 
            B(imatB+j1,imatB+j2) = b(ibtot+nef(k)+ncontr(k)+j); // on remplit la partie triangulaire superieure (de la matrice) correspondant a la dimension
            j += 1; 
          } 
        } 
        if((imatB>0) & (nRE>0)){ // si des EA sont correles entre les dimensions
          for(j3=0; j3<sum(nea(span(0,k-1))); j3++){ 
            B(j3,imatB+j2) = b(npmMM+iRE); // on remplit la partie triangulaire superieure (de la matrice)
            iRE +=1; 
          } 
        } 
      } 
    } 
    
    // BM du sous-modele k 
    if(ncor(k)==1){ // si la dimension k contient un mvmt brownien
      BM(imatBM,imatBM) = b(ibtot+nef(k)+ncontr(k)+nvc(k)); 
      if((imatBM>0) & (nBM>0)){  // jamais utilise : mvmts browniens correles
        for(j3=0; j3<imatBM; j3++){ 
          BM(j3,imatBM)=b(npmMM+nRE+iBM); 
          iBM += 1; 
        } 
      } 
      imatBM += 1; 
    } 
    
    
    ibeta += nef(k)+1; 
    ibtot += nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+ntrtot(k); 
    imatB += nea(k); 
  } 
  

  // matrice B (mat par blocs, varcov des effets aleatoires des K dimensions)
  endj2=B.n_cols;
  for(j2=0;j2<endj2;j2++){
    for(j1=0;j1<=j2;j1++){
        B(j2,j1)=B(j1,j2);
    }
  }
  
  //B.print("B = ");
  

  // BM
  endj2=BM.n_cols;
  for(j2=0;j2<endj2;j2++){
    for(j1=0;j1<=j2;j1++){
        BM(j2,j1)=BM(j1,j2);
    }
  }


  /////////////////////////////////////////////////////
  // creation des vecteurs eta et des matrices gamma //  
  ///////////////////////////////////////////////////// 
  
  // vecteur  eta (K seuils) 
  
  eta = b(span(npmMM+nRE+nBM,npmMM+nRE+nBM+nvarD)); // recupere dans b les thresholds (intercept + var expl indep)
  if(ndept>0){ // si il y a des variables explicatives dependantes du temps pr l'event
    eta_dept = b(span(npmMM+nRE+nBM+nvarD+1,npmMM+nRE+nBM+nvarD+1+ndept-1)); 
  } 

  
  //gamma
  
  l=0; 
  for(k=0;k<K;k++){ // pr chq dimension
    gamma(k,span(0,nvarD)) = trans(b(span(npmMM+nRE+nBM+1+nvarD+ndept+l,npmMM+nRE+nBM+1+nvarD+ndept+l+nvarD)));  // matrice dimension*contribution
    l += nvarD+1; 
  } 

  //////////////////////////////////////
  // initialisation de pls parametres //  
  ////////////////////////////////////// 
  
  nmescurry=zeros<ivec>(K);  // nmescurry pour Y 
  nmescurry(0)=0; 
  m0=0; 
  if(K>1){  // si il y a plus d'une dimension
    for(k=1; k<K; k++){ // pr chq dimension
      nmescurry(k) = sum(nmes(span(0,m0+ny(k-1)-1))); // nbre d'obs des marqueurs de la dimension k
      m0 += ny(k-1); 
    } 
  } // obs rangees par dimension, outcome puis temps 
  // nmescurry = vecteur de taille K avec nbre d'obs des marqueurs pour chq dimension
  
  ni = sum(nmes); // nbre total de valeurs obs des marqueurs 
  
  // vecteurs et matrices a ni lignes
  Hyi.zeros(ni);  
  muiK.zeros(ni); 
  RiK.zeros(ni,ni);
  Xi.zeros(ni,beta.size()); 
  Xcontri.zeros(ni,bcontr.size()); 
  Zi.zeros(ni,sum(nea));
  tcor.zeros(ni); 
  
  
  /// pr ce qui concerne AVANT le tps s
  muiDs.zeros(K*(Bs+1)); 
  RiDs.zeros(K*(Bs+1),K*(Bs+1));
  RiKDs.zeros(K*(Bs+1),ni);
  Xdis.zeros(K*(Bs+1),beta.size()); 
  Zdis.zeros(K*(Bs+1),B.n_cols);    
  tdcors.zeros(K*(Bs+1)); 
  
  
  /// pr ce qui concerne AVANT le tps s+t
  muiDst.zeros(K*(Bst+1)); 
  RiDst.zeros(K*(Bst+1),K*(Bst+1));
  RiKDst.zeros(K*(Bst+1),ni);
  Xdist.zeros(K*(Bst+1),beta.size()); 
  Zdist.zeros(K*(Bst+1),B.n_cols);    
  tdcorst.zeros(K*(Bst+1)); 

  
  inodes=0; 
  ibtot=0; // repere des parametres de la dimension en cours (ds b)
  jcurr=0; // repere du nbre d'obs des dimensions precedentes (pr X)
  jdcurrs=0; // repere du nbre d'obs des dimensions precedentes (pr Xd)
  jdcurrst=0; // repere du nbre d'obs des dimensions precedentes (pr Xd1)
  iidx=0; // repere des differentes dimensions dans idx
  p=0; 
  pcontr=0; 
  q=0; 
  m0=0; 
  kbm=0; 
  
  
  //////////////////////////
  // Construction de H(Y) //  
  ////////////////////////// 
  
  
  for(k=0; k<K; k++){	  // pr chq dimension  
    
    nik = sum(nmes(span(m0,m0+ny(k)-1))); // nbre de valeurs obs des marqueurs de la dimension k
    
    tmpnmes=0; // repere des obs du marqueur en cours (ds Y et HY)
    tmpntr=0; // repere des parametres de la fct lien du marqueur en cours (ds b)
    
    //transfos : calcul de Y
    for(m=m0; m<m0+ny(k); m++){ // pr chq marqueur
      
      if(nmes(m)==0) continue; // !ne devrait pas arriver, gerer dans R! , ie. si il n'y a pas de valeur obs pr ce marqueur
      
      yim = Y(span(nmescurry(k)+tmpnmes,nmescurry(k)+tmpnmes+nmes(m)-1)); // obs de Y du marqueur m de la dimension k // Y^(s)
      
      Hyim.zeros(nmes(m)); 
      
      if(link(m)==0){ // si la fct de lien du marqueur est lineaire
        Hyim = (yim - b(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr)) / b(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr+1);
        tmpntr += 2; // 2 parametres si lineaire :
        inodes += 2; // min et max	
      } 
      /// omis : fct de lien beta
      else{
        if(link(m)==2){ // si la fct de lien du marqueur est splines
          // pr les splines, il faut regarder ou se situe le Y considere pra noeuds
          // additionner les splines finis avant
          // calculer les valeurs des splines en cours a ce moment : valeurs dans mm, mm1 et mm2 (car 3 splines en cours max)
          // ne pas se preocupper des splines qui n'ont pas encore commences (ceux apres)
          // H(Y) = som + prm1*spline1 + prm2*spline2 + prm3*spline3
          splaa.zeros(ntr(m)); // splaa = vecteur des parametres splines de la fct de lien
          splaa(0) = b(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr);
          for(kk=1; kk<ntr(m); kk++){ // pr chq parametre splines
            splaa(kk) = b(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr+kk)*b(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+tmpntr+kk);
          }
          // repeter les noeuds dans zitr (~splaa en repetant pls fois le min et le max)
          zitr.zeros(ntr(m)+2); 
          zitr(0)=nodes(inodes); // premier noeud-splines du marqueur 
          zitr(1)=zitr(0); // REPET : premier noeud-splines du marqueur
          zitr(span(2,ntr(m)-1))=nodes(span(inodes,inodes+ntr(m)-3));
          zitr(ntr(m))=zitr(ntr(m)-1); 
          zitr(ntr(m)+1)=zitr(ntr(m)-1); 
          // calcul de H(y) 
          for(j=0;j<nmes(m);j++){ // pr chq valeur du marqueur m
            ll=0; // on regarde ou se situe y pra noeuds-splines et on stocke la position dans ll
            if(yim(j)==nodes(inodes+ntr(m)-3)){ // si y = max
              ll=ntr(m)-2;   
            } 
            for(kk=3;kk<ntr(m);kk++){ // pr chq noeud
              if((yim(j)>=zitr(kk-1)) && (yim(j)<zitr(kk))){ // si y compris entre ce noeud et le precedent
                ll=kk-1; 
              } 
            } 
            if((ll<2) || (ll>ntr(m)-2)){ // ERREUR
              proba=-1E9;
              goto fin; 
            } 
                 som=splaa(0); 
                 if(ll>2){ 
                   som=som+sum(splaa(span(1,ll-2))); // on somme les splines finis
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
                 if(yim(j)==zitr(ntr(m)-1)){ 
                   mm2=0.0; 
                   mm1=0.0; 
                   mm=3.0/h1; 
                 } 
                 else{ 
                   mm2=(3.0*ht2*ht2)/(hh*h1*hn); 
                   mm1=(3.0*htm*ht2)/(h2n*hh*h1)+(3.0*ht*ht3)/(h2*h1*h2n); 
                   mm=(3.0*ht*ht)/(h3*h2*h1); 
                 }
                 if((mm<0) || (mm1<0) || (mm2<0)){ 
                   proba=-1E9;
                   goto fin; 
                 } 
                 im2 = hht*mm2/(3.0) + h2n*mm1/(3.0) + h3*mm/(3.0); 
                 im1 = htm*mm1/(3.0) + h3*mm/(3.0); 
                 im = ht*mm/(3.0); 
                 
                 Hyim(j) = som + splaa(ll-1)*im2 + splaa(ll)*im1 + splaa(ll+1)*im; 
                 
                 
          } 
          
          inodes += ntr(m)-2; 
          tmpntr += ntr(m); 
          
        } 
      } 
      
      Hyi(span(jcurr+tmpnmes,jcurr+tmpnmes+nmes(m)-1)) = Hyim;  // vecteur 'general' H(y) _ ts marqueurs, ttes dimensions confondus
      // jcurr = nbre d'obs des dimensions precedentes (pr X)
      // tmpnmes = nbre d'obs des marqueurs precedents de la dimension en cours
      // nmes = nbre d'obs du marqueur en cours de la dimension en cours
      
      tmpnmes += nmes(m); 
      

    }// fin boucle m 
    

    //////////////////////////////////////////////////////
    // Construction de X^(s), ~X, ~X+, Z^(s), ~Z et ~Z+ //  
    //////////////////////////////////////////////////////
    
    
    // toutes les variables de la k-ieme dimension sont dans l'element k de la liste X   /// indicage liste commence a 1 ??
    Xk = Rcpp::as<arma::mat>(X[k]); // k-ieme matrice de la liste X
    Xdk = Rcpp::as<arma::mat>(Xd[k]); // k-ieme matrice de la liste Xd
    
    for(l=0; l<nv(k); l++){ // pr chq variable explicative de la dimension k
      
      if(idx(0,iidx+l)==1){ // si cette variable explicative est un effet fixe  
        // X^(s)
        Xi(span(jcurr,jcurr+nik-1),p) = Xk(span(0,nik-1),l); 
        //~X
        Xdis(span(jdcurrs,jdcurrs+Bs),p) = Xdk(span(0,Bs),l); 
        //~X+
        Xdist(span(jdcurrst,jdcurrst+Bst),p) = Xdk(span(0,Bst),l);
        p += 1; 
      } 
      if(idx(1,iidx+l)==1){ // si cette variable explicative est un contraste 
        tmpnmes=0; // repere du nbre d'obs du marqueur en cours
        ll=0; 
        for(m=m0; m<m0+ny(k); m++){ // pr chq marqueur de la dimension en cours
          Xcontri(span(jcurr+tmpnmes,jcurr+tmpnmes+nmes(m)-1),ny(k)*pcontr+ll) = Xk(span(tmpnmes,tmpnmes+nmes(m)-1),l); 
          ll += 1; 
          tmpnmes += nmes(m); 
        } 
        pcontr += 1; // repere des marqueurs precedents de la dimension en cours // Xcontri rangee en colonne par dimension, var expl ac contraste, marqueur
      } 
      if(idx(2,iidx+l)==1){ // si cette variable explicative est un effet aleatoire 
        // Z^(s)
        Zi(span(jcurr,jcurr+nik-1),q) = Xk(span(0,nik-1),l); 
        // ~Z
        Zdis(span(jdcurrs,jdcurrs+Bs),q) = Xdk(span(0,Bs),l); 
        // ~Z+
        Zdist(span(jdcurrst,jdcurrst+Bst),q) = Xdk(span(0,Bst),l); 
        q += 1; 
      } 
      if(idx(3,iidx+l)==1){ // si cette variable explicative est un temps pour BM ou AR 
        //W^(s)
        tcor(span(jcurr,jcurr+nik-1)) = Xk(span(0,nik-1),l); // vecteurs des valeurs du BM de la dimension, ac des 0 si la dim ne contient pas de BM
        //~W
        tdcors(span(jdcurrs,jdcurrs+Bs)) = Xdk(span(0,Bs),l); 
        //~W+
        tdcorst(span(jdcurrst,jdcurrst+Bst)) = Xdk(span(0,Bst),l); 
      } 
    }
    
    
    iidx += nv(k);  
    

    //////////////////////////////////
    // Construction de R, ~R et ~R+ //   a partir de BM et tcor pr mvmt brownien, et de b pr intercept aleatoire specifique au marqueur (alpha) et erreurs de mesure (epsilon)
    //////////////////////////////////
    
    
    // construire Rk, Rd et Rkd (correlations BM ou AR ) 
    if(ncor(k)>0){ // si la dimension contient un processus d'autocorrelation
      
      if(ncor(k)==1){ // si il s'agit d'un BM 
        for(j1=0; j1<nik; j1++){ 
          for(j2=0; j2<nik; j2++){ 
            // R^(s)
            RiK(jcurr+j1,jcurr+j2) = BM(kbm,kbm)*std::min(tcor(jcurr+j1),tcor(jcurr+j2));  // cov(w(t),w(u))=sigmaW^2 x min(t,u) ou sigmaW^2 stocke dans B, t et u dans tcor
          } 
          for(j3=0; j3<Bs+1; j3++){ 
            // ~~R
            RiKDs(jdcurrs+j3,jcurr+j1) = BM(kbm,kbm)*std::min(tcor(jcurr+j1),tdcors(jdcurrs+j3));
          }
          for(j3=0; j3<Bst+1; j3++){
            // ~~R+
            RiKDst(jdcurrst+j3,jcurr+j1) = BM(kbm,kbm)*std::min(tcor(jcurr+j1),tdcorst(jdcurrst+j3));
          }
        } 
        for(j1=0; j1<Bs+1; j1++){ 
          for(j2=0; j2<Bs+1; j2++){ 
            // ~R
            RiDs(jdcurrs+j1,jdcurrs+j2) = BM(kbm,kbm)*std::min(tdcors(jdcurrs+j1),tdcors(jdcurrs+j2)); 
          }
        }
        for(j1=0; j1<Bst+1; j1++){
          for(j2=0; j2<Bst+1; j2++){
            // ~R+
            RiDst(jdcurrst+j1,jdcurrst+j2) = BM(kbm,kbm)*std::min(tdcorst(jdcurrst+j1),tdcorst(jdcurrst+j2));  
          }
        }
        kbm += 1; // repere des dimensions avec mvmt brownien pr chercher sigmaW^2 ds mat BM 
      } 
      // possibilite : mvmts browniens correles et AR mais omis ici
    } 
    

    // alpha (intercept aleatoire specifique au marqueur) 
    if(nalea(k)>0){ // si la dimension contient des intercepts aleatoires specifiques a ses marqueurs
      ll=0; // repere de l'intercept aleatoire marqueur-specific ds b
      tmpnmes=0; // repere du nbre d'obs du marqueur en cours
      for(m=m0; m<m0+ny(k); m++){ // pr chq marqueur de la dimension
        for(j1=0; j1<nmes(m); j1++){ // pr chq obs du marqueur
          for(j2=0; j2<nmes(m); j2++){ 
            RiK(jcurr+tmpnmes+j1,jcurr+tmpnmes+j2) +=   
              b(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+ll)*b(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+ll); // au carre pr avoir la variance
          } // ecart-type de l'intercept aleatoire marqueur-specific dans b
        } 
        tmpnmes += nmes(m); 
        ll += 1; 
      } 
    } 

    
    // erreur de mesure 
    ll=0; // repere de l'erreur de mesure ds b
    tmpnmes=0; // repere du nbre d'obs du marqueur en cours
    for(m=m0; m<m0+ny(k); m++){ // pr chq marqueur de la dimension
      for(j1=0; j1<nmes(m); j1++){ // pr chq obs du marqueur
        RiK(jcurr+tmpnmes+j1,jcurr+tmpnmes+j1) +=  // rempli la diagonale de la matrice Ri^(s)
          b(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ll)*b(ibtot+nef(k)+ncontr(k)+nvc(k)+ncor(k)+ll); // au carre pr avoir la variance
      } // ecart-type de l'erreur dans b
      tmpnmes += nmes(m); 
      ll += 1; 
    } 
    
    jcurr += nik; // iteration, pr passer a la dimension suivante
    ibtot += nef(k)+ncontr(k)+nvc(k)+ncor(k)+ny(k)+nalea(k)+ntrtot(k); 
    nmescurry(k) += nik; 
    jdcurrs += Bs+1; 
    jdcurrst += Bst+1; 
    m0 += ny(k); 
    
  } // fin boucle k 

  
 
  //////////////////////////////////////////////////////////////////////////////
  // Construction de muiK (moyenne de H(Y)) et de VK (mat de var-cov de H(Y)) //
  //////////////////////////////////////////////////////////////////////////////
  
  
  muiK = Xi*beta + Xcontri*bcontr; 
  
  VK = Zi*B*trans(Zi) + RiK ; // V_H(Y^(s))
  
  //aa = det(VK); // si le determinant de VK est inferieur a 0 -> matrice non-inversible
  // if(aa<1.0E-10){
  //   Zi.print("Zi=");
  //   B.print("B=");
  //   VK.print("VK=");
  //   cout << "det(VK)=" <<  det(VK)  << endl;
  //   b.print("b=");
  //   proba=-1E9;
  //   goto fin;
  // }
  
  
  invVK = inv_sympd(VK); // inverse de la matrice VK // [V_H(Y^(s))]^(-1)
  

  //////////////////////////////////////////////////////////////////////
  // Construction des elements intervenant dans la fct de repartition //
  //////////////////////////////////////////////////////////////////////
  
  //--------------------------------------------------------------------//
  // D = vecteur des valeurs des var expl (survie) independantes du tps //
  //--------------------------------------------------------------------//
  
  
  varD.fill(0);
  varD(0) = 1.0; // pr intercept
  if(nvarD>0){ // si il y a des variables explicatives independantes du temps
    l=1; 
    for(j=0;j<nvarD;j++){ // pr chacune des ces variables
      varD(l) = D(j+1); // alors on recupere sa valeur
      l += 1; 
    } 
  } 

  
  //----------------------------------------------------------------------//
  // Gam = matrice contenant les contributions des dimensions au proc deg //
  //----------------------------------------------------------------------//
  
  Gam.zeros(Bs+1,K*(Bs+1));
  GG.eye(Bs+1,Bs+1);   // matrice identite
  Gam = kron(trans(gamma*varD),GG);
  
  //Gam.print("Gam =");
  
  GamPlus.zeros(Bst+1,K*(Bst+1));
  GGPlus.eye(Bst+1,Bst+1);   // matrice identite
  GamPlus = kron(trans(gamma*varD),GGPlus);
  
  //GamPlus.print("GamPlus =");
    
  
  tGam=trans(Gam); // transposee
  tGamPlus=trans(GamPlus);
  
  //----------------//
  // les seuils eta //
  //----------------//
  
  seuil.ones(Bs+1);
  // pr les var expl independantes du temps
  seuil(span(0,Bs)) *= dot(eta,varD);  // s'assurer que vecteur = 1 nbre revient a vecteur = (nbre,...,nbre)
  // pr les var expl dependantes du temps
  if(ndept>0){
    for(j=0 ; j<Bs+1 ; j++){
      seuil(j) += dot(eta_dept,trans(Xseuil(j,span(0,ndept-1))));
    }
  }
  
  //seuil.print("seuil = ");
  
  
  seuilPlus.ones(Bst+1);
  // pr les var expl independantes du temps
  seuilPlus(span(0,Bst)) *= dot(eta,varD);  // s'assurer que vecteur = 1 nbre revient a vecteur = (nbre,...,nbre)
  // pr les var expl dependantes du temps
  if(ndept>0){
    for(j=0 ; j<Bst+1 ; j++){
      seuilPlus(j) += dot(eta_dept,trans(Xseuil(j,span(0,ndept-1))));
    }
  }
  
  //seuilPlus.print("seuilPlus = ");


  //---------------------------------------------------------------------------//
  // construction des moyennes et des variances intervenant dans la fct de rep //
  //---------------------------------------------------------------------------//
  
  //moyenne de Lambda
  muiDs = Xdis*beta; 
  
  
  //covariance de Y et Lambda
  VKDs = Zdis*B*trans(Zi) + RiKDs;  
  
  
  //variance de Lambda
  VDs = Zdis*B*trans(Zdis) + RiDs;  
  
  
  // moyenne de Lambda sachant Y
  muiD_Ks = muiDs + VKDs*invVK*(Hyi-muiK);  
  
  
  // variance de Lambda sachant Y
  VD_Ks = VDs-VKDs*invVK*trans(VKDs);  

  
  // muiDs.print("muiDs = ");
  // VKDs.print("VKDs = ");
  // VDs.print("VDs = ");
  // muiD_Ks.print("muiD_Ks = ");
  // VD_Ks.print("VD_Ks = ");
  
  
  
  //moyenne de Lambda+
  muiDst = Xdist*beta;
  
  
  //covariance de Y et Lambda+
  VKDst = Zdist*B*trans(Zi) + RiKDst;
  
  
  //variance de Lambda+
  VDst = Zdist*B*trans(Zdist) + RiDst;  
  
  
  // moyennes de Lambda+ sachant Y
  muiD_Kst = muiDst + VKDst*invVK*(Hyi-muiK); 
  
  
  // variances de Lambda+ sachant Y
  VD_Kst = VDst-VKDst*invVK*trans(VKDst);  
  
  
  
  // muiDst.print("muiDst = ");
  // VKDst.print("VKDst = ");
  // VDst.print("VDst = ");
  // muiD_Kst.print("muiD_Kst = ");
  // VD_Kst.print("VD_Kst = ");
  

  /////////////////////////////////////
  // Calculs des fcts de repartition //
  /////////////////////////////////////
  
  
  //~~~~~~~~~~~~~~~~~~~~~~~~~~//
  // FctRep de dimension Bs+1 //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  dimPhi1 = Bs+1;
  upper1.zeros(Bs+1); // borne superieure de la fct de repartition
  lower1.zeros(Bs+1); // borne inferieure de la fct de repartition
  infin1.zeros(Bs+1); // infin1.ones(Bs+1); 
  
  upper1 = seuil - Gam*muiD_Ks; // borne de l'integration  // point sur lequel on calcule la fct de repartition
  
  V1.zeros(Bs+1,Bs+1); 
  I1.eye(Bs+1,Bs+1);   // matrice identite
  V1 = I1 + Gam*VD_Ks*tGam; // variance ds la fct de rep 
  
  //cout << "dimPhi1 =" << dimPhi1 << endl;
  //upper1.print("upper1 = ");
  //V1.print("V1 = ");
  
  lower1 = upper1;  // parametres de la fct sadmvn
  maxpts1 = 10000*dimPhi1; 
  abseps = 0.001; // seuils de tolerance
  releps = 0.0001; 
  error1 = 0; // initialisation
  value1 = 0; 
  inform1 = 0; 
  
  // creation du vecteur des correlations a partir de V1 (pr pouvoir utiliser fct sadmvn)
  correl1.zeros((Bs+1)*(Bs+1+1)/2-Bs-1); // nbre de parametres de correlation a trouver = nbre d'elements de la matrice V1 - la moitie (symetrie) - diago (corr=1)
  
  endj1 = Bs+1;  
  
  upper11.zeros(dimPhi1); 
  
  ll=0; 
  for(j2=0; j2<endj1; j2++){ 
    upper11(j2) = upper1(j2)/pow(V1(j2,j2),0.5); // pow(x,p) = x^p  // seuil <- seuil/variance
    if(j2>0){ 
      for(j1=0; j1<j2; j1++){
        correl1(ll) = V1(j1,j2)/(pow(V1(j1,j1),0.5)*pow(V1(j2,j2),0.5)); // corr(1,2)=cov(1,2)/(|sigma1|.|sigma2|) ou sigma = racine de la variance
        ll += 1; 
      } 
    }
  } 
  

  lower111 = conv_to< std::vector<double> >::from(lower1);  // pr type standard (accepte par sadmvn)
  upper111 = conv_to< std::vector<double> >::from(upper11); 
  infin111 = conv_to< std::vector<int> >::from(infin1); 
  correl111 = conv_to< std::vector<double> >::from(correl1); 


  // calcul de la fct de repartition
  sadmvn_(&dimPhi1,&lower111[0],&upper111[0],&infin111[0],&correl111[0],&maxpts1,&abseps,&releps,&error1,&value1,&inform1); 
  // resultat dans value1
  
  
  if(inform1!=0){
    cout << "prblm integration pr PHI s :" << "value1=" << value1 << " error1=" << error1 << endl;
    //   V1.print("V1=");
    //   upper1.print("upper1=");
    //   correl1.print("correl1=");
    proba=-1E9;
    goto fin;
  }

  
  if(value1==0){
    cout << "prblm value1 nul :" << "value1=" << value1;
    //   V1.print("V1=");
    //   upper1.print("upper1=");
    //   correl1.print("correl1=");
    proba=-1E9;
    goto fin;
  }
  

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  // FctRep de dimension Bst+1 //
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~//
  
  dimPhi2 = Bst+1;
  upper2.zeros(Bst+1); // borne superieure de la fct de repartition
  lower2.zeros(Bst+1); // borne inferieure de la fct de repartition
  infin2.zeros(Bst+1); // infin.ones(Bst+1); 
  
  upper2 = seuilPlus - GamPlus*muiD_Kst; // borne de l'integration  
  
  V2.zeros(Bst+1,Bst+1); 
  I2.eye(Bst+1,Bst+1);   // matrice identite
  V2 = I2 + GamPlus*VD_Kst*tGamPlus; // variance ds la fct de rep 
  
  //cout << "dimPhi2 =" << dimPhi2 << endl;
  //upper2.print("upper2 = ");
  //V2.print("V2 = ");
  
  lower2 = upper2;  // parametres de la fct sadmvn
  maxpts2 = 10000*dimPhi2; 
  abseps = 0.001; // seuils de tolerance
  releps = 0.0001; 
  error2 = 0; // initialisation
  value2 = 0; 
  inform2 = 0;
  
  // creation du vecteur des correlations a partir de V2 (pr pouvoir utiliser fct sadmvn)
  correl2.zeros((Bst+1)*(Bst+1+1)/2-Bst-1); // nbre de parametres de correlation a trouver = nbre d'elements de la matrice V2 - la moitie (symetrie) - diago (corr=1)
  
  endj2 = Bst+1;  
  
  upper22.zeros(dimPhi2); 
  
  ll=0; 
  for(j2=0; j2<endj2; j2++){ 
    upper22(j2) = upper2(j2)/pow(V2(j2,j2),0.5); // pow(x,p) = x^p  // seuil <- seuil/variance
    if(j2>0){ 
      for(j1=0; j1<j2; j1++){
        correl2(ll) = V2(j1,j2)/(pow(V2(j1,j1),0.5)*pow(V2(j2,j2),0.5)); // corr(1,2)=cov(1,2)/(|sigma1|.|sigma2|) ou sigma = racine de la variance
        ll += 1; 
      } 
    } 
  } 
  
  lower222 = conv_to< std::vector<double> >::from(lower2);  // pr type standard (accepte par sadmvn)
  upper222 = conv_to< std::vector<double> >::from(upper22); 
  infin222 = conv_to< std::vector<int> >::from(infin2); 
  correl222 = conv_to< std::vector<double> >::from(correl2); 

  
  // calcul de la fct de repartition
  sadmvn_(&dimPhi2,&lower222[0],&upper222[0],&infin222[0],&correl222[0],&maxpts2,&abseps,&releps,&error2,&value2,&inform2); 
  // resultat dans value2

  if(inform2!=0){
    cout << "prblm integration pr PHI s+t :" << "value2=" << value2 << " error2=" << error2 << endl;
    //   V2.print("V2=");
    //   upper2.print("upper2=");
    //   correl2.print("correl2=");
    
    proba=-1E9;
    goto fin;
  }
  
  
  
  //////////////////////////////
  // Calcul de la probabilite //
  //////////////////////////////
  
  //cout << "value1 = " << value1 << endl;
  //cout << "value2 = " << value2 << endl;

  proba = 1-value2/value1; // proba d avoir l event entre s et s+t
  
  cout << "proba = " << proba << endl;
  

  fin: 
  return Rcpp::wrap(proba); 
  
} // fin RcppExport
