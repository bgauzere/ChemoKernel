/**
 * @file LaplacianKernel.cpp
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#include <algorithm>
#include "LaplacianKernel.h"
#include "GraphEditDistance.h"
#include <cmath>
#include <cfloat>
#include <assert.h>
#include <iomanip>
#include <sys/times.h>
#include <sys/time.h>

#include <unistd.h>

using namespace std;
using namespace cimg_library;
using namespace pandore;


LaplacianKernel::LaplacianKernel(const LaplacianKernel * clone){
  edit = clone->edit;
  sigma = clone->sigma;
  regularization = clone->regularization;
  lambda = clone->lambda;
  K = clone->K;
  W = clone->W;  
  D = clone->D;  
  L = clone->L;
  U_K = clone->U_K;
  _collections = clone->_collections;
}

LaplacianKernel::LaplacianKernel (GraphEditDistance* edit, const Dataset & trainset, const Dataset & testset, double sigma, int regularization, double lambda) : edit(edit), sigma(sigma), regularization(regularization), lambda(lambda)
{	
  unsigned int N = trainset.size();
  unsigned int M = testset.size();
	
  for (unsigned int i=0; i<N; ++i)
    _collections.push_back(trainset.getCollection(i));
	
  for (unsigned int i=0; i<M; ++i)
    _collections.push_back(testset.getCollection(i));
	
  computeKernel();
}

LaplacianKernel::LaplacianKernel (GraphEditDistance* edit, const Dataset & trainset, double sigma, int regularization,  double lambda) : edit(edit), sigma(sigma), regularization(regularization), lambda(lambda)
{	
  unsigned int N = trainset.size();
	
  for (unsigned int i=0; i<N; ++i)
    _collections.push_back(trainset.getCollection(i));
	
  computeKernel();
}

double LaplacianKernel::computeWeight(pandore::Collection* m1, pandore::Collection* m2){
  double dist = (*edit) (m1, m2);
  return exp(-dist/(sigma*sigma)); ///2.0
}

void LaplacianKernel::computeKernel ()
{	
  //Calcul complet
  unsigned int N = _collections.size();
	
  W.resize(N,N); // The weights matrix
  D.resize(N,N); // The degree matrix
  L.resize (N,N); // The graph Laplacian matrix
  U_K.resize(N, N);

  K.resize(N, N);
	
  // Computing the weights and degree matrix

  for (unsigned int i = 0; i<N; ++i)
      for (unsigned int j=0; j<N; ++j)
	{
	  D(i,j) = 0;
	  if (i != j)
	    {
	      W(i,j) = computeWeight(_collections[i], _collections[j]);
	    }
	  else
	    W(i,j)=0;
	}
  
  for (unsigned int i=0; i<N; ++i)
    {
      double d = 0;
      for (unsigned int j=0; j<N; ++j)
	d += W(i,j);
		
      D(i,i) = d;
    }
  // Computing the graph Laplacian matrix
  // //Normalized Version
  // for (unsigned int i=0; i<N; ++i)
  //   for (unsigned int j=0; j<N; j++)
  //     {
  // 	if (i==j)
  // 	  L(i,j)=1;
  // 	else if ((i!=j) && (W(i,j) != 0))
  // 	  L(i,j)=-W(i,j)/sqrt(D(i,i)*D(j,j));
  // 	else 
  // 	  L(i,j)=0;
  //     }

  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; j++)
      {
  	if (i==j)
  	  L(i,j) = D(i,i) - W(i,j);
  	else if ((i!=j) && (W(i,j) != 0))
  	  L(i,j)=-W(i,j);
  	else 
  	  L(i,j)=0;
      }
  //L is the **un** normalized Laplacian

  // Computing the Laplacian Moore-Penrose pseudoinverse
  CImg<double> eigvals;
  CImg<double> eigvects;
  
  //Regularisation
  for (unsigned int i=0; i<N; i++)
    for (unsigned int j=0; j<N; j++)
      {
  	if(i == j)
  	  L(i,j) = 1 + lambda*L(i,j);
  	else
  	  L(i,j) = lambda*L(i,j);
      }
  
  U_K = L;
  
  clock_t start,end;
  struct timeval s_start, s_end;
  start = gettimeofday(&s_start,NULL);
  U_K = U_K.invert();
  end = gettimeofday(&s_end,NULL);//times(&s_end);
  
  long diff = (1000000*s_end.tv_sec + s_end.tv_usec) - 
    (1000000*s_start.tv_sec + s_start.tv_usec);
  cout.precision(DBL_DIG);
  cout << "Standard Inversion computation time (micros) : " << 
    diff << endl;


  // L.symmetric_eigen(eigvals,eigvects);

  // K=(double) 0;
	
  // for (unsigned int i=0; i<N; ++i)
  //   {
  //     //cout << eigvals[i] << endl;
  //     double r;
  //     switch (regularization)
  // 	{
  // 	case 0:
  // 	  r = eigvals[i]; // No regularization
  // 	  break;
  // 	case 1:
  // 	  r = 1 +lambda*eigvals[i]; // Regularized Laplacian
  // 	  break;
  // 	case 2:
  // 	  r = exp(eigvals[i]*lambda*lambda/2); // Diffusion Process
  // 	  break;
  // 	case 3:
  // 	  r = 1/(cos(eigvals[i]*M_PI/4)); // Inverse cosine
  // 	  break;
  // 	case 4:
  // 	  r = 1/(lambda-eigvals[i]); // One-step Random Walk
  // 	  break;
  // 	}
      
  //     assert(fabs(r) > 0.00001);
      
  //     CImg<double> tmp (N,N);
  //     for (unsigned int j=0; j<N; ++j)
  // 	tmp(j,j)= eigvects(i,j)*eigvects(i,j);
  //     for (unsigned int j=0; j<N; ++j)
  // 	for (unsigned int k=j+1; k<N; ++k)
  // 	  {
  // 	    tmp(j,k) = eigvects(i,j)*eigvects(i,k);
  // 	    tmp(k,j)=tmp(j,k);
  // 	  }
  //     K += (1/r)*tmp;
  //     //U_K += (1/r)*tmp;
  //   }

  
  //Normalized Laplacian
  for (unsigned int i=0; i<N; i++)
    for (unsigned int j=0; j<N; j++)
      K(i,j) = sqrt(D(i,i))*U_K(i,j)*sqrt(D(j,j));

  normalizeK();
  
}

void LaplacianKernel::normalizeK(){
  //  Normalization
  unsigned int N  = _collections.size();
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int j=0; j<N; ++j)
      if (i != j)
  	K(i,j) /= sqrt(K(i,i)*K(j,j));  

  for (unsigned int i=0; i<N; ++i)
    K(i,i)=1;
}

//  double LaplacianKernel::operator() (Collection* c1, Collection* c2)
// {	
//   int i=-1;
//   int j=-1;
 	
//   for (unsigned int k=0; k<_collections.size(); ++k)
//     if (_collections[k] == c1)
//       {
// 	i=k;
// 	break;
//       }
//   for (unsigned int k=0; k<_collections.size(); ++k)
//     if (_collections[k] == c2)
//       {
// 	j=k;
// 	break;
//       }
	
//   bool compute = false;

//   if (i == -1)
//     {
//       i=_collections.size();
//       _collections.push_back(c1);
//       cout << "need c1" << endl;
//       compute = true;
//     }
	
//   if (j == -1)
//     {
//       j=_collections.size();
//       _collections.push_back(c2);
//       cout << "need c2" << endl;
//       compute = true;
//     }
	
//   if (compute) //The kernel needs to be recomputed
//     computeKernel();
//   return K(i,j);
// }


double LaplacianKernel::operator() (Collection* c1, Collection* c2)
{	
  int i=-1;
  int j=-1;
	
  for (unsigned int k=0; k<_collections.size(); ++k)
    if (_collections[k] == c1)
      {
	i=k;
	break;
      }
	
  for (unsigned int k=0; k<_collections.size(); ++k)
    if (_collections[k] == c2)
      {
	j=k;
	break;
      }
	
  //bool compute = false;
	
  if (i == -1)
    {
      i=_collections.size();
      //Rajout dans la collection
      clock_t start,end;
      struct timeval s_start, s_end;
      start = gettimeofday(&s_start,NULL);
      fastAdd(c1,0.0001);
      end = gettimeofday(&s_end,NULL);//times(&s_end);
      long diff = (1000000*s_end.tv_sec + s_end.tv_usec) - 
	(1000000*s_start.tv_sec + s_start.tv_usec);
      cout.precision(DBL_DIG);
      cout << "Fast Add (micros) : " << diff << endl;
      
    }
	
  if (j == -1)
    {
      j=_collections.size();
      clock_t start,end;
      struct timeval s_start, s_end;
      start = gettimeofday(&s_start,NULL);
      fastAdd(c2,0.0001);
      end = gettimeofday(&s_end,NULL);//times(&s_end);
      long diff = (1000000*s_end.tv_sec + s_end.tv_usec) - 
	(1000000*s_start.tv_sec + s_start.tv_usec);
      cout.precision(DBL_DIG);
      cout << "Fast Add (micros) : " << diff << endl;
    }
	
  return K(i,j);
}

void LaplacianKernel::printWeightMatrix(){
  unsigned int N = _collections.size();
  cout.precision(DBL_DIG);
  for (unsigned int i = 0; i<N; ++i)
    {
      for (unsigned int j=0; j<N; ++j)
	{
	  double dist = (*edit) (_collections[i], _collections[j]);
	  cout << exp(-dist/(sigma*sigma)) << " ";
	  
	}
      cout << endl;
    }
}

void LaplacianKernel::printMatrix(CImg<double> m){
  if(PRINT_MATRIX)
    {
      for (int i=0; i< m.width(); ++i)
	{
	  for (int j=0; j< m.height(); ++j)
	    {
	      cout.precision(DBL_DIG);
	      cout <<  setw(20) << setiosflags(ios::left) << m(i,j)  << "\t";
	    }
	  cout << endl;
	}
    }
}

void LaplacianKernel::testFastInversion(Collection * c_new){

  //Calcul complet de la matrice de Gram avant ajout de c_new
  unsigned int N = _collections.size();

  // W.resize(N,N); // The weights matrix
  // D.resize(N+1,N+1); // The degree matrix
  // L.resize (N,N); // The graph Laplacian matrix
  // K.resize(N, N); // The normalized gram Matrix
  // U_K.resize(N, N); // The inverse of the unnormalized Laplacian Matrix
	
  // // Computing the weights and degree matrix
  // for (unsigned int i = 0; i<N; ++i)
  //   for (unsigned int j=0; j<N; ++j)
  //     {
  // 	D(i,j) = 0;
  // 	if (i != j)
  // 	  {
  // 	    W(i,j) = computeWeight(_collections[i], _collections[j]);
  // 	  }
  // 	else
  // 	  W(i,j) = 0;
  //     }
	
  // for (unsigned int i=0; i<N; ++i)
  //   {
  //     double d = 0.0;
  //     for (unsigned int j=0; j<N; ++j)
  // 	d += W(i,j);
  //     D(i,i) = d;
  //   }
  // //D and W ok !

  // //l tilde
  // for (unsigned int i=0; i<N; ++i)
  //   for (unsigned int j=0; j<N; j++)
  //     {  	
  // 	if (i==j)
  // 	  L(i,j) = D(i,i) - W(i,j);
  // 	else if ((i!=j) && (W(i,j) != 0))
  // 	  L(i,j)=-W(i,j);
  // 	else 
  // 	  L(i,j)=0;
  //     }

  // if(PRINT_MATRIX)
  //   {
  //     cout << "Delta_n : " << endl;
  //     printMatrix(L);
  //   }
  
  // //Regularisation
  // for (unsigned int i=0; i<N; i++)
  //   for (unsigned int j=0; j<N; j++)
  //     {
  // 	if(i == j)
  // 	  L(i,j) = 1 + lambda*L(i,j);
  // 	else
  // 	  L(i,j) = lambda*L(i,j);
  //     }
  
  // //Verification si inversion de matrice possible
  // U_K = L;
  // double det_L = L.det();
  
  // if(abs(det_L) > 0.0001)
  //   U_K = U_K.invert();
  // else
  //   {
  //     cout << "Matrice presque singuliere (det = " << det_L << "),  abandon..." << endl;
  //     return;
  //   }

  // if(PRINT_MATRIX)
  //   {
  //     cout << "(Delta_n)^-1 : " << endl;
  //     printMatrix(U_K);
  //   }

  // //Normalisation
  // for (unsigned int i=0; i<N; i++)
  //   for (unsigned int j=0; j<N; j++)
  //     K(i,j) = U_K(i,j)*(sqrt(D(i,i))*sqrt(D(j,j)));
  
  computeKernel();

  CImg<double> delta_n(N,N,1,1,0.0);//La Matrice Diagonale ajoutée
  CImg<double> b(N,1,1,1,0.0);//La ligne ajoutée
  double max_w = 0.0;
  double sum_w = 0.0; //D(i,i)
  double sum_b = 0.0;
  for(unsigned int i = 0;i<N;i++)
    {
      double w = computeWeight(_collections[i],c_new);
      sum_w +=  w;
      b(i) = -(lambda * w);
      delta_n(i,i) = b(i);
      sum_b += b(i);
      D(i,i) += w;
      max_w = max(abs(delta_n(i,i)), max_w); //XXX:Attention au lambda
    }
  //colonne = ligne transposée
  CImg<double> b_t = b.transpose();
  double d = 1 - sum_b;

  D.resize(N+1,N+1); // The degree matrix
  D(N,N) = sum_w;

  /* Verifications */
  assert(b(10) == b(0,10));
  assert(max_w <1.0);//XXX:< 1.0
  
  CImg<double> converge = U_K*delta_n;
  double norm = converge.magnitude(2);
  if(norm < 1)
    cout << "Norm ok : ||(Delta_n^-1*delta_n)|| = " << norm << "." << endl;
  else
    cout << "La suite ne converge pas : ||(Delta_n^-1*delta_n)|| = " << norm << "." << endl;
  
  if(PRINT_MATRIX)
    {
      cout << "delta_n : " << endl;
      printMatrix(delta_n);
    }
  
  double epsilon = 0.0001;
  int k_max = (log(epsilon*(1-max_w))/log(max_w));
  //  int k_max = (log(2.0*epsilon) / log(max_w))+1;
  cout << "nb iterations : " << k_max << endl;

  CImg<long double> delta_n_inverse (N,N,1,1,0.0); //GRAND Delta n ou E
  CImg<long double> pow_k (N,N,1,1,0.0);
  CImg<long double> pow_delta_n (N,N,1,1,0.0);
  pow_delta_n = pow_delta_n.identity_matrix();
  pow_k = U_K;
  //k_max << N => Gain en complexité
  for(int k = 0;k < k_max; k++)
    {
      for(unsigned int i = 0;i < N; i++)
	for(unsigned int j = 0;j < N; j++)
	  delta_n_inverse(i,j) += pow_k(i,j)*pow_delta_n(i,i);
      pow_k = pow_k * U_K;
      pow_delta_n = pow_delta_n * delta_n;
    }	
  //Verif Inversion rapide
  verifInversion(delta_n_inverse,L,delta_n);

  //delta_n_inverse(i,j) += pow(-1,k)* pow(U_K(i,j),k+1)*pow(delta_n(i,i),k);

  //XXX:Tests
  // delta_n_inverse = (L+delta_n).invert();
  // cout << "Verif = 0 normalement ..." << endl;
  // verifInversion(delta_n_inverse,L,delta_n);
  //XXX:Fin Test

  //Inversion par bloc
  //http://fr.wikipedia.org/wiki/Matrice_inversible#Inversion_par_bloc
  CImg<double> c_a = CImg<double>(N,1,1,1,0.0);
  for(unsigned int i=0;i<N;i++)
    for(unsigned int k=0;k<N;k++)
      c_a(i) += b(k)*delta_n_inverse(i,k);	//C * A^-1 
  
  CImg<double> a_b = CImg<double>(1,N,1,1,0.0);
  for(unsigned int j=0;j<N;j++)
    for(unsigned int k=0;k<N;k++)
      a_b(0,j) +=  delta_n_inverse(k,j)*b_t(0,k);	// A^-1 * B

  long double schur_complement = 0.0;
  for(unsigned int i=0;i<N;i++)
    schur_complement += c_a(i)*b_t(0,i);
  schur_complement = 1.0/(d - schur_complement);
  cout.precision(DBL_DIG);
  cout << "schur_complement : " << schur_complement << endl ;
  
  CImg<double> Gram = CImg<double>(N+1,N+1,1,1,0.0);
  Gram(N,N) = schur_complement;
  for(unsigned int i=0;i<N;i++)
    {
      Gram(N,i) = - a_b(0,i)*schur_complement;
      Gram(i,N) = - schur_complement*c_a(i);
      for(unsigned int j=0;j<N;j++)
	Gram(i,j) = delta_n_inverse(i,j) + a_b(0,j)*schur_complement*c_a(i);
    }

  //  Normalized Laplacian
  for (unsigned int i=0; i<N+1; ++i)
    for (unsigned int j=0; j<N+1; ++j)
      //if (i != j)
  	Gram(i,j) = sqrt(D(i,i)) * Gram(i,j) * sqrt(D(j,j));

  //  Normalization
  for (unsigned int i=0; i<N+1; ++i)
    for (unsigned int j=0; j<N+1; ++j)
      if (i != j)
  	Gram(i,j) /= sqrt(Gram(i,i)*Gram(j,j));  

  for (unsigned int i=0; i<N+1; ++i)
   Gram(i,i)=1;

  //Done
  // for (unsigned int i=0; i<N+1; ++i)
  //   {
  //     for (unsigned int j=0; j<N+1; ++j)
  // 	{
  // 	  cout.precision(DBL_DIG);
  // 	  cout <<  Gram(i,j)  << " ";
  // 	}
  //     cout << endl;
  //   }
  
  //Verif difference
  _collections.push_back(c_new);
  computeKernel();
  CImg<double> Diff = K - Gram;

  if(PRINT_MATRIX)
    {
      cout << "K : " << endl;
      printMatrix(K);
  
      cout << "Gram : " << endl;
      printMatrix(Gram);
      
      cout << "K - Gram : " << endl;
      printMatrix(Diff);
    }
  
  cout << "Norme 1 de la diff : " << Diff.magnitude(1) << endl;
  cout << "Norme 2 de la diff : " << Diff.magnitude(2) << endl;
  cout << "Norme inf de la diff : " << Diff.magnitude(-1) << endl;
  int i_max,j_max;
  double diff_max  = 0.0;
  double sum_diff = 0.0;
  for(int i = 0;i < Diff.width(); i++)
    for(int j = 0;j < Diff.height(); j++)
      {
	sum_diff += abs(Diff(i,j));
	if(abs(Diff(i,j))>diff_max)
	  {
	    diff_max = abs(Diff(i,j));
	    i_max = i;
	    j_max = j;
	  }
	// if(abs(Diff(i,j)) > 0.000001){
	//   cout << i << "," << j << " : ";
	//   cout << Gram(i,j) << " | ";
	//   cout << K(i,j) << " => " << Diff(i,j) <<  endl;

	// }
	  
	//diff_max = (abs(Diff(i,j))>diff_max)?abs(Diff(i,j)):diff_max;
      }
  
  cout << "Difference Max : " << diff_max << " on (" << i_max <<","<<j_max<<")" << endl;
  cout << "Difference Moyenne : " << sum_diff/(Diff.width()*Diff.height()) << endl;
  cout << Gram(i_max,j_max) << endl;
  cout << K(i_max,j_max) << endl;
  cout << L(i_max,j_max) << endl;

  //TODO : cleanup
	
}


void LaplacianKernel::verifInversion(CImg<double> matrice_inv,
				     CImg<double> Delta_n,
				     CImg<double> delta_n){
  CImg<double> T = Delta_n-delta_n;
  
  double det_T = T.det();
  if(abs(det_T) > 0.0001)
    T = T.invert();
  else
    {
      cout << "Matrice presque singuliere (det = " << det_T << "),  abandon..." << endl;
      return;
    }
  
  CImg<double> Diff = matrice_inv - T;
   
  if(PRINT_MATRIX)
    {
      cout << "(Delta_n+delta_n)^-1: " << endl;
      printMatrix(T);
  
      cout << "Matrice inversee : " << endl;
      printMatrix(matrice_inv);
      
      cout << "Diff : " << endl;
      printMatrix(Diff);
    }

  cout << "Norme 1 de la diff : " << Diff.magnitude(1) << endl;
  cout << "Norme 2 de la diff : " << Diff.magnitude(2) << endl;
  cout << "Norme inf de la diff : " << Diff.magnitude(-1) << endl;
  double diff_max  = 0.0;
  double sum_diff = 0.0;
  for(int i = 0;i < Diff.width(); i++)
    for(int j = 0;j < Diff.height(); j++)
      {
	diff_max = (abs(Diff(i,j))>diff_max)?abs(Diff(i,j)):diff_max;
	sum_diff += abs(Diff(i,j));
      }
  
  cout << "Difference Max : " << diff_max << endl;
  cout << "Difference Moyenne : " << sum_diff/(Diff.width()*Diff.height()) << endl;
  
  
}

//Rajouter lambda pour la regul
void LaplacianKernel::fastAdd(Collection * c_new, double epsilon){
  //The kernel is already computed for current dataset
  //Laplacian (L) regularized with I and Lambda 

  int N = _collections.size();
  CImg<double> delta_n(N,N,1,1,0.0);//La Matrice Diagonale ajoutée
  CImg<double> b(N,1,1,1,0.0);//La ligne ajoutée
  double max_w = 0.0;
  double sum_w = 0.0; //D(i,i)
  double sum_b = 0.0;

  for(int i = 0;i<N;i++)
    {
      double w = computeWeight(c_new,_collections[i]);
      sum_w +=  w;
      b(i) = -(lambda * w);
      delta_n(i,i) = - b(i);
      sum_b += b(i);
      D(i,i) += w;
      max_w = max(abs(delta_n(i,i)), max_w); //XXX:Attention au lambda
    }
  //colonne = ligne transposée
  CImg<double> b_t = b.transpose();
  double d = 1 - sum_b;
  D.resize(N+1,N+1); // The degree matrix
  D(N,N) = sum_w;
  
  //int k_max = (log(epsilon*(1-max_w))/log(max_w));
  int k_max = (log(2.0*epsilon) / log(max_w))+1;
  cout << "k : " << k_max << endl;

  clock_t start,end;
  //struct tms s_start, s_end;
  struct timeval s_start, s_end;
  CImg<double> delta_n_inverse (N,N,1,1,0.0);
  CImg<long double> pow_k[k_max];
  //pre computing of U_K powers
  pow_k[0] = U_K;
  for(int k = 1; k < k_max; k++)
    pow_k[k] = pow_k[k-1] * U_K;

  start = gettimeofday(&s_start,NULL);
  //k_max << N => Gain en complexité
  for(int k = 0;k < k_max; k++)
    {
      for(int i = 0;i < N; i++)
	for(int j = 0;j < N; j++)
	  delta_n_inverse(i,j) += pow_k[k](i,j)*pow(delta_n(i,i),k);
      //delta_n_inverse(i,j) += pow(U_K(i,j),k+1)*pow(delta_n(i,i),k);
    }
  
  CImg<double> c_a(N,1,1,1,0.0);
  for(int i=0;i<N;i++)
    for(int k=0;k<N;k++)
      c_a (i) += b(k)*delta_n_inverse(i,k);	//C * A^-1 
  
  //TODO : doubles boucles a regrouper 
  CImg<double> a_b(1,N,1,1,0.0);
  for(int j=0;j<N;j++)
    for(int k=0;k<N;k++)
      a_b (0,j) +=  delta_n_inverse(k,j)*b_t(0,k);	// A^-1 * B
  
  long double schur_complement = 0.0;
  //(C*A^-1)*B
  for(int i=0;i<N;i++)
    schur_complement += c_a(i)*b_t(0,i);
  schur_complement = 1.0/(d - schur_complement);
  
  _collections.push_back(c_new);
  U_K.resize(_collections.size(), _collections.size(),1,1);
  
  for(int i=0;i<N;i++)
    {
      U_K(N,i) = - a_b(0,i)*schur_complement;
      U_K(i,N) = - schur_complement*c_a(i);
      for(int j=0;j<N;j++)
	U_K(i,j) = delta_n_inverse(i,j) + a_b(0,j)*schur_complement*c_a(i);
        //	K(i,j) = delta_n_inverse(i,j) + s_tmp;
    }
  end = gettimeofday(&s_end,NULL);//times(&s_end);
  //long nb_top =  sysconf(_SC_CLK_TCK);
  long diff = (1000000*s_end.tv_sec + s_end.tv_usec) - 
    (1000000*s_start.tv_sec + s_start.tv_usec);
  cout.precision(DBL_DIG);
  cout << "Fast Inversion computation time (micros) : " << 
    diff << endl;

  // //Normalized Laplacian
  K.resize(_collections.size(), _collections.size(),1,1);
  for (unsigned int i=0; i<N+1; ++i)
    for (unsigned int j=0; j<N+1; ++j)
      K(i,j) = sqrt(D(i,i)) * U_K(i,j) * sqrt(D(j,j));

  normalizeK();
}



