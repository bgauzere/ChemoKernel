/*
 * @file WeightedRegression.cpp
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Wed Jun 15 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 *
 */

#include "WeightedRegression.h"
#include "MoleculeGraph.h"
#include <cassert>
#include <unistd.h>

using namespace std;
using namespace cimg_library;
using namespace pandore;


#define NB_ITERATIONS 10000

bool param_condition(double x)
{
  return (x==x);
}

WeightedRegression::WeightedRegression(MoleculesDataset * trainset)
{
  _trainset = trainset;
  _N = _trainset->size();
  _y = CImg<double> (1,_N);
  for(unsigned int i=0;i<_N;i++)
    _y(0,i) = trainset->getParameter(i);
  #ifdef VISU
  _visu = CImgDisplay(500,1024,"w",0);
  #endif

}

CImg<double> WeightedRegression::getKTreelet(string treelet, int treelet_type)
{
  for(unsigned int i=0;i<_k_to_treelet.size();i++)
    if((_k_to_treelet[i].first.compare(treelet) == 0) && 
       (treelet_type == _k_to_treelet[i].second))
      return _Kks[i];
  cerr << "Treelet G" << treelet_type << "-" << treelet << " not found." << endl;
  return CImg<double> (_N,_N,1,1,0.0);
}

void WeightedRegression::computeKks()
{
  treelet_spectrum * distribution = _trainset->getTreeletDistribution(param_condition);
  _T = 0;
  /*Parcours de tous les motifs de treelets*/
  for(int treelet_type=0;treelet_type < SIZE_SPECTRUM; treelet_type++){
    treelet_spectrum::iterator it = distribution[treelet_type].begin();
    for(;it != distribution[treelet_type].end();it++)
      {
	string treelet = it->first;
	CImg<double> Kk(_N,_N,1,1,0.0);
	//XXX:Symétrique, à optimiser
	for(unsigned int i=0;i<_N;i++)
	  for(unsigned int j=0;j<_N;j++)
	    {
	      //XXX: We assume that the spectrum must be computed
	      pandore::Collection* c1 = _trainset->getCollection(i);
	      pandore::Collection* c2 = _trainset->getCollection(j);
	      treelet_spectrum ** spec_1 = 
		(treelet_spectrum **) c1->GETVALUE("labeled_spectrum",Long);
	      treelet_spectrum ** spec_2 = 
		(treelet_spectrum **) c2->GETVALUE("labeled_spectrum",Long);
	      treelet_spectrum::iterator it1 = spec_1[treelet_type]->find(treelet);
	      treelet_spectrum::iterator it2 = spec_2[treelet_type]->find(treelet);
	      double fi_k, fj_k;
	      //fi_k = 0 si pas de treelets dans la molécule c1
	      (it1 == spec_1[treelet_type]->end())?fi_k=0:fi_k=it1->second; 
	      (it2 == spec_2[treelet_type]->end())?fj_k=0:fj_k=it2->second;

	      if( fi_k !=0 && fj_k !=0)
		{
		  double diff = fi_k - fj_k;
		  double d_rbf = exp(-(diff * diff)/(2*_sigma*_sigma));
		  Kk(i,j) = d_rbf;
		}
	      else//Version K1
		Kk(i,j) = 0.0;
	    }
	if (_normalize) // We normalize the Gram matrix
	  {
	    for (unsigned int i=0; i<_N; ++i)
	      for (unsigned int j=0; j<_N; ++j)
		if (i != j)
		  Kk(i,j) /= sqrt(Kk(i,i)*Kk(j,j));
	    
	    for (unsigned int i=0; i<_N; ++i)
	      Kk (i,i)=1;
	  }
	_Kks.push_back(Kk);
	// //Verification Semi DP
	// CImg<double> eigvals;
	// CImg<double> eigvects;
	// Kk.symmetric_eigen(eigvals,eigvects);
	// cout << "Treelet G"<< treelet_type << "-" << treelet << endl;
	// cout << eigvals[_N-1] << endl;
	_k_to_treelet.push_back(pair<string,int>(treelet,treelet_type));
	_T++;
      }
  }
}

double WeightedRegression::objectiveFunction(CImg<double> w,   CImg<double> alpha){
  assert((unsigned int)w.height() == _T);
  assert((unsigned int)alpha.height() == _N);


  //Calcul de Kw(kgraph,w)
  CImg<double> k_w = Kw(w);
  //Verification Semi DP
  // CImg<double> eigvals;
  // CImg<double> eigvects;
  // k_w.symmetric_eigen(eigvals,eigvects);
  // cout << "KW" << endl;
  // for(int i=0;i<_N;i++)
  //   cout << eigvals[i] << ", ";
  // cout << endl;
  
  CImg<double> k_alpha = k_w*alpha;
  
  //Calcul de reg_1 : 
  double reg_1 = (alpha.get_transpose()*k_alpha)(0,0);
  reg_1 *= _lambda;

  
  // CImg<double> A (_T,_N,1,1,0.0); //A \in R NxT
  // for(int j=0;j<_T;j++) //Vecteur colonne Aj
  //   {
  //     CImg<double>  kk_alpha = _Kks[j]*alpha;
  //     for(int i=0;i<_N;i++)  //Vecteur ligne Ai
  // 	A(j,i) = kk_alpha(0,i);
  //   }
  
  // CImg<double> a_alpha = A.get_transpose()*alpha;
  // double reg_1 = (a_alpha.get_transpose()*w)(0,0);
  // reg_1 *= _lambda;

  
  //Calcul de fit
  double fit=0.0;
  for(int i=0;i<k_alpha.height();i++)
    fit += pow(k_alpha(0,i) - _y(0,i),2);

  // CImg<double> a_w = A*w;  
  // double fit=0.0;
  // for(unsigned int i=0;i<a_w.height();i++)
  //   fit += pow(a_w(0,i) - _y(0,i),2);

  //Calcul de reg_2
  double reg_2 = 0.0;
  for(int i=0;i<w.height();i++)
    reg_2 += abs(w(0,i));
  reg_2 *= _mu;
  
  cerr << "(" << reg_1 << ", " << fit << ", " << reg_2 << ")" ;
  return reg_1+fit+reg_2;
}

CImg<double> WeightedRegression::getOptimalAlpha()
{
  return _optimalAlpha;
}
CImg<double> WeightedRegression::getOptimalWeights()
{
  return _optimalWeights;
}

CImg<double> WeightedRegression::Kw(CImg<double> w)
{
  CImg<double> Kw (_N,_N,1,1,0.0);
  for(unsigned int k=0;k<_T;k++)
    Kw += w(0,k)*_Kks[k];
  if(_normalize)
    {
      for (unsigned int i=0; i<_N; ++i)
	for (unsigned int j=0; j<_N; ++j)
	  if (i != j)
	    Kw(i,j) /= sqrt(Kw(i,i)*Kw(j,j));
	    
      for (unsigned int i=0; i<_N; ++i)
	Kw (i,i)=1;

    }
  return Kw;
}

CImg<double> WeightedRegression::gradientAlpha(CImg<double> w, 
					       CImg<double> alpha)
{
  assert((w!=w) && (alpha!=alpha));
  // CImg<double> K = Kw(w);
  // CImg<double> tmp = (CImg<double>(_N,_N).identity_matrix()*_lambda + K)*alpha.get_transpose();
  // return 2*K*(tmp-_y.get_transpose());
  return CImg<double>::empty();
}

CImg<double> WeightedRegression::descentAlpha(CImg<double> w, 
					      CImg<double> alpha)
{
  CImg<double> tmp = Kw(w) + (CImg<double>(_N,_N).identity_matrix()*_lambda);
  tmp.invert();
  return (tmp*_y);
}


// CImg<double> WeightedRegression::descentAlpha(CImg<double> w, 
// 					      CImg<double> alpha)
// {
//   CImg<double> alpha_t;
//   CImg<double> alpha_t_1 = alpha;
//   CImg<double> gradient;
//   do
//     {
//       alpha_t = alpha_t_1;
//       gradient = gradientAlpha(alpha_t,w);
//       alpha_t_1 = alpha_t - gradient; //XXX: Quel coeff pour le gradient ?
//     }
//   while(gradient.magnitude() > _epsilon_alpha);
  
//   return alpha_t_1;
// }


CImg<double> WeightedRegression::positiveHardTresholding(CImg<double> w, double threshold)
{
  CImg<double> w_st = w;
  for(unsigned int k=0;k<_T;k++)
    if(w(0,k) <= threshold)
      w_st(0,k) = 0;
    else
      w_st(0,k) = w(0,k);
  return w_st;
}

CImg<double> WeightedRegression::hardTresholding(CImg<double> w, double threshold)
{
  CImg<double> w_st = w;
  for(unsigned int k=0;k<_T;k++)
    if(abs(w(0,k)) <= threshold)
      w_st(0,k) = 0;
    else
      w_st(0,k) = w(0,k);
  return w_st;
}


CImg<double> WeightedRegression::positiveSoftTresholding(CImg<double> w, double threshold)
{
  CImg<double> w_st = w;
  for(unsigned int k=0;k<_T;k++)
    if(w(0,k) <= threshold)
      w_st(0,k) = 0;
    else
      w_st(0,k) = w(0,k) - (threshold);
  return w_st;
}


CImg<double> WeightedRegression::softTresholding(CImg<double> w, double threshold)
{
  CImg<double> w_st = w;
  for(unsigned int k=0;k<_T;k++)
    if(abs(w(0,k)) <= threshold) 
      w_st(0,k) = 0;
    else
      {
	int sign = (w(0,k) > 0)?1:-1;
	w_st(0,k) = w(0,k) - (threshold*sign);
      }
  return w_st;
}

CImg<double> WeightedRegression::descentW(CImg<double> w, 
					  CImg<double> alpha)
{
  CImg<double> A (_T,_N,1,1,0.0); //A \in R NxT
  for(unsigned int j=0;j<_T;j++) //Vecteur colonne Aj
    {
      CImg<double>  kk_alpha = _Kks[j]*alpha;
      for(unsigned int i=0;i<_N;i++)  //Vecteur ligne Ai
  	A(j,i) = kk_alpha(0,i);
    }
  

  //Verif
  //A(i=12,j=104) = 
  // double test = 0.0;
  // for(int i=0;i<_N;i++)
  //   test += _Kks[104](i,12)*alpha(0,i);
  
  // cout << "test " << (test == A(104,12)) << endl;

  CImg<double> gradient;
  CImg<double> descente_gradient;
  int test_nb_iterations = 0;
  // double num = 0.0;
  // for(int k=0;k<_T;k++)
  //   {
  //     CImg<double> eigvals;
  //     CImg<double> eigvects;
  //     _Kks[k].symmetric_eigen(eigvals,eigvects);
  //     num += eigvals[0]*eigvals[0];
  //   }
  // double norme_alpha_sq = 0.0;
  // for(int i=0;i<_N;i++)
  //   norme_alpha_sq += alpha(0,i)*alpha(0,i);
  // num *= norme_alpha_sq;
  
  // double tau_t = 1/(num*50);

  CImg<double> eigvals;
  CImg<double> eigvects;
  (A.get_transpose()*A).symmetric_eigen(eigvals,eigvects);
  
  double max = 0.0;
  for(unsigned int i=0;i<eigvals.size();i++)
    {
      max = (abs(eigvals[i]) > max)?abs(eigvals[i]):max;
    }
  // double tau_t = 0.99/(max);
  // double tau_t = 0.5/(max);
  double tau_t = 1.001*(max);
  CImg<double> w_t;
  CImg<double> w_t_1 = w;  
  CImg<double> gradient_comp;


  do
    {
      w_t = w_t_1;
      // CImg<double> p1 = A.get_transpose();
      // CImg<double> p2 = (A*w_t);
      // CImg<double> p3 = _lambda*alpha;
      // CImg<double> p4 = _y-p2;
      // CImg<double> p5 = p4 - p3;
      // gradient = p1*p5;
      // gradient = A.get_transpose()*(_y-(A*w_t)-(_lambda*alpha));
      // descente_gradient = gradient * (2*tau_t);
      
      gradient = A.get_transpose()*((_lambda*alpha)-(2*_y)+(2*A*w_t));
      descente_gradient = -tau_t * gradient;
      
      w_t_1 = hardTresholding(w_t + descente_gradient, _mu*tau_t);
      //w_t_1 = w_t + descente_gradient;
      
      /*VISU*/
#ifdef VISU
      CImg<double> visu_w_t_1 = w_t_1;
      for(int i=0;i<visu_w_t_1.size();i++)
	visu_w_t_1[i] *= 255.0/(w_t_1.max());
      _visu.display(visu_w_t_1);
      /*Fin Visu*/
#endif
      gradient_comp = w_t_1 - w_t;
      cerr <<  "\r" << gradient_comp.magnitude(2) << ", " << objectiveFunction(w_t_1,alpha);
      cerr << "(" << test_nb_iterations << "/" << NB_ITERATIONS << ")"; 
      //usleep(1000000);
      test_nb_iterations++;
    }
  //while(test_nb_iterations < NB_ITERATIONS);
  while((gradient_comp.magnitude() > _epsilon_w) && (test_nb_iterations < NB_ITERATIONS)) ;
  cerr << "\r" << test_nb_iterations << " Itérations" << endl;
  return w_t; //t_1
}
  
  
void WeightedRegression::computeMinimisation(TreeletKernel* kernel,
					     CImg<double> w_0, CImg<double> alpha_0)
{
  cerr << "Initialisation ... " << endl;
  computeKks();
  cerr << "Kks computed ... " << endl;
  
  if(! w_0)
    {
      w_0 = CImg<double> (1,_T,1,1,1.0); 

      // srand (125);
      // w_0 = CImg<double> (1,_T,1,1,0.0);
      // for(int i=0; i<_T;i++)
      // 	w_0(0,i) = 1.0/(rand() % 10 +1);
    }
  if(! alpha_0)
    {
      alpha_0 = descentAlpha(w_0,alpha_0);
    }
  
  assert((unsigned int)w_0.height() == _T);
  assert((unsigned int)alpha_0.height() == _N);
  
  CImg<double> w_t;
  CImg<double> alpha_t;
  
  CImg<double> w_t_1 = w_0;
  CImg<double> alpha_t_1 = alpha_0;
  int nb_iterations =0;
  double cost_t,cost_t1;
  cost_t1= objectiveFunction(w_0, alpha_0);
  cerr << endl <<"Initialisation done : "<< cost_t1 << endl;

  do
    {
      w_t = w_t_1;
      alpha_t = alpha_t_1;
      //Descente de gradient selon alpha
      alpha_t_1 = descentAlpha(w_t, alpha_t);
      
      //Descente de gradient selon w
      w_t_1 = descentW(w_t,alpha_t_1);
      
      cerr << "\r" << "Iteration " << nb_iterations++;
      cost_t = cost_t1;
      cost_t1 = objectiveFunction(w_t_1, alpha_t_1);
      cerr << ", " << cost_t << ","<< cost_t1 << endl;
    }
  while((nb_iterations < (NB_ITERATIONS)) && (cost_t - cost_t1 > _epsilon));
  //while(cost_t - cost_t1 > _epsilon);
  _optimalWeights = w_t_1; //w_t ?
  _optimalAlpha = alpha_t_1; 
  cerr << "done" << endl;
}


double WeightedRegression::objectiveFunctionV2(CImg<double> w)
{
  CImg<double> alpha (1,_N);
  CImg<double> tmp = Kw(w) + (CImg<double>(_N,_N).identity_matrix()*_lambda);
  tmp.invert();
  alpha = (tmp*_y);
  return objectiveFunction(w,alpha);
}

CImg<double> WeightedRegression::descentW_V2(CImg<double> w)
{

  CImg<double> gradient(1,_T);
  CImg<double> descente_gradient;
  int test_nb_iterations = 0;
  
  CImg<double> eigvals;
  CImg<double> eigvects;
  
  CImg<double> w_t;
  CImg<double> w_t_1 = w;  
  CImg<double> gradient_comp;
  double tau_t = 0.0001; //XXX:A determiner autrement
  do {
      w_t = w_t_1;
      
      CImg<double> alpha (1,_N);
      CImg<double> tmp = Kw(w_t) + (CImg<double>(_N,_N).identity_matrix()*_lambda);
      tmp.invert();
      alpha = (tmp*_y);
      
      for(unsigned int k=0;k<_T;k++)
	{
	  gradient[k] = - ((_lambda*alpha.get_transpose())*_Kks[k]*(alpha))(0,0);
	}
      
      descente_gradient = -tau_t * gradient;
      
      w_t_1 = softTresholding(w_t + descente_gradient, _mu*tau_t);
      //w_t_1 = w_t + descente_gradient;
      
      /*VISU*/
#ifdef VISU
      CImg<double> visu_w_t_1 = w_t_1;
      for(int i=0;i<visu_w_t_1.size();i++)
	visu_w_t_1[i] *= 255.0/(w_t_1.max());
      _visu.display(visu_w_t_1);
      /*Fin Visu*/
#endif
      gradient_comp = w_t_1 - w_t;
      cerr <<  "\r" << gradient_comp.magnitude(2) << ", " << objectiveFunctionV2(w_t_1);
      cerr << "(" << test_nb_iterations << "/" << NB_ITERATIONS << ")"; 
      test_nb_iterations++;
    }
  //while(test_nb_iterations < NB_ITERATIONS);
  while(gradient_comp.magnitude() > _epsilon_w);
  cerr << "\r" << test_nb_iterations << " Itérations" << endl;
  return w_t; //t_1
}

 
void WeightedRegression::computeMinimisationV2(TreeletKernel* kernel,
					       CImg<double> w_0)
{
  cerr << "Initialisation ... " << endl;
  computeKks();
  cerr << "Kks computed ... " << endl;
  
  if(! w_0)
    {
      w_0 = CImg<double> (1,_T,1,1,1.0); 
    }
  
  assert((unsigned int)w_0.height() == _T);
  
  CImg<double> w_t  = w_0;
  CImg<double> w_t_optim;
  
  cerr << "Initialisation done : "<< endl;
  w_t_optim = descentW_V2(w_t);
  _optimalWeights = w_t_optim;
  cerr << "done" << endl;
}

 
CImg<double> WeightedRegression::vit_to_weights(vector<string> * vit_list)
{
  computeKks();//Required to initialize
  CImg<double> w(1,_T,1,1,0);
  for(int treelet_type=0;treelet_type<SIZE_SPECTRUM;treelet_type++)
    {
      for(unsigned int i=0;i<vit_list[treelet_type].size();i++)
	{
	  string code =  vit_list[treelet_type][i];
	  for(unsigned int k=0;k<_k_to_treelet.size();k++)
	    {
	      if((_k_to_treelet[k].first.compare(code) == 0) && 
		 _k_to_treelet[k].second == treelet_type)
		{
		  assert(w[k] != 1.0);
		  w[k] = 1.0;
		}
	    }
	}
    }
  return w;
}
