/**
 * @file WeightedRegression.h
 * @author Benoit Gauzere <<benoit.gauzere@greyc.ensicaen.fr>> 
 * @version     0.0.1 - Wed Jun 15 2011
 * 
 * @todo the list of improvements suggested for the file.
 * @bug the list of known bugs.
 *  
 * Description of the program objectives.
 * All necessary references.
 */

#ifndef __WEIGHTEDREGRESSION_H__
#define __WEIGHTEDREGRESSION_H__

#include <vector>
#include <utility>
#include "CImg.h"
#include "MoleculesDataset.h"
#include "TreeletKernel.h"

using namespace std;

class WeightedRegression
{
private:
  cimg_library::CImg<double> _optimalWeights;
  cimg_library::CImg<double> _optimalAlpha;

  //Function to minimize
  double objectiveFunction(cimg_library::CImg<double> w,   
			   cimg_library::CImg<double> alpha);
  double objectiveFunctionV2(cimg_library::CImg<double> w);
  
  //TODO : A modifier pour stocker une fois calculé 
  // Cimg * Kks= new Cimg[T]; Kks[k] = Kk
  void computeKks();
  
  unsigned int _T; //nb treelets
  unsigned int _N; //nb Molecules
  
  //Descente de gradient selon alpha
  cimg_library::CImg<double> gradientAlpha(cimg_library::CImg<double> w, 
					   cimg_library::CImg<double> alpha);
  
  cimg_library::CImg<double> descentAlpha(cimg_library::CImg<double> w, 
					  cimg_library::CImg<double> alpha);
  //Descente de gradient selon W
  cimg_library::CImg<double> gradientW(cimg_library::CImg<double> w, 
				       cimg_library::CImg<double> alpha);
  cimg_library::CImg<double> descentW(cimg_library::CImg<double> w, 
				      cimg_library::CImg<double> alpha);
  cimg_library::CImg<double> descentW_V2(cimg_library::CImg<double> w);
  //Permet de calculer le noyau pondéré
  cimg_library::CImg<double> Kw(cimg_library::CImg<double> w);
  
  //Propriété à prédire
  cimg_library::CImg<double> _y;
  
  //Fonction de Seuillage 
  cimg_library::CImg<double> hardTresholding(cimg_library::CImg<double> w, double threshold);
  cimg_library::CImg<double> positiveHardTresholding(cimg_library::CImg<double> w, double threshold);
  cimg_library::CImg<double> softTresholding(cimg_library::CImg<double> w, double threshold);
  cimg_library::CImg<double> positiveSoftTresholding(cimg_library::CImg<double> w, double threshold);
  cimg_library::CImgDisplay _visu;
public:
  MoleculesDataset * _trainset;

  //Ensembles des noyaux combinés
  vector< cimg_library::CImg<double> > _Kks;
  vector< pair<string,int> > _k_to_treelet; //Association indice/<treelet,type treelet>
  
  // Conditions de convergence 
  double _epsilon; 
  double _epsilon_alpha;
  double _epsilon_w; 
  
  bool _normalize;
  
  //Parametres du noyau
  double _sigma;   // Sigma du noyau gaussien
  double _mu;      // Coeff regul norme 1 de w
  double _lambda;  // Coeff regul alpha 
  
  //Constructeur
  WeightedRegression(MoleculesDataset * trainset);
  
  cimg_library::CImg<double> getKTreelet(string treelet, int treelet_type);
  
  //Récupération des résultats
  cimg_library::CImg<double> getOptimalAlpha();
  cimg_library::CImg<double> getOptimalWeights(); 
  
  //Calcul des poids et alpha optimaux
  void computeMinimisation(TreeletKernel * kernel,
			   cimg_library::CImg<double> w_0 = cimg_library::CImg<double>::empty(), 
			   cimg_library::CImg<double> alpha_0 = cimg_library::CImg<double>::empty() );
  
  void computeMinimisationV2(TreeletKernel * kernel,
			     cimg_library::CImg<double> w_0 = cimg_library::CImg<double>::empty());
  
  cimg_library::CImg<double> vit_to_weights(vector<string> * vit_list);
  
};

#endif // __WEIGHTEDREGRESSION_H__
