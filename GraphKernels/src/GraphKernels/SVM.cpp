/**
 * @file SVM.cpp
 *
 * @author Alice KIJEWSKI <alice.kijewski@ecole.ensicaen.fr>
 * @author David LEMARESQUIER <david.lemaresquier@ecole.ensicaen.fr>
 *
 * @version 1.1.0 (2010-07-21)
 */

#include "SVM.h"
#include "CImg.h"

using namespace cimg_library;
using namespace pandore;

SVM::SVM (const Dataset & trainset, GraphKernel* kgraph, int type, double c, double nu, double epsilon) : kgraph(kgraph)
{
	unsigned int N = trainset.size();
	
	for (unsigned int i=0; i<N; ++i)
		cols.push_back(trainset.getCollection(i));
	
	// Preparation of the SVM parameters

	param.svm_type = type;
	param.kernel_type = PRECOMPUTED;
	param.degree = 3;
	param.gamma = 0;	// 1/num_features
	param.coef0 = 0;
	param.nu = nu;
	param.cache_size = 100;
	param.C = c;
	param.eps = 1e-3;
	param.p = epsilon;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;	
	
	// Construction of the SVM problem
	
	CImg<double> gram = trainset.getGramMatrix(false);
	
	prob.l = N;
	prob.y = new double[N];
	prob.x = new svm_node*[N];
	
	for (unsigned int i=0; i<N; ++i)
	{
		prob.y[i] = trainset.getParameter(i);
		
		prob.x[i] = new svm_node[N+2];
		
		prob.x[i][0].index = 0;
		prob.x[i][0].value = i+1;
		
		for (unsigned int j=0; j<N; ++j)
		{
			prob.x[i][j+1].index = j+1;
			prob.x[i][j+1].value = gram(i,j);
		}
		
		prob.x[i][N+1].index = -1;
	}
	
	// Compute the SVM model
	
	model = svm_train(&prob, &param);
}

SVM::~SVM ()
{
	for (int i=0; i<prob.l; ++i)
		delete[] prob.x[i];
	
	delete[] prob.x;
	delete[] prob.y;
	
	svm_destroy_param(&param);
	svm_destroy_model(model);
}

double SVM::predict (Collection* col) const
{
	int N = cols.size();
	
	svm_node* x = new svm_node[N+2];
	
	x[0].index = 0;
	x[0].value = 0;
	
	for (int j=0; j<N; ++j)
	{
		x[j+1].index = j+1;
		x[j+1].value = (*kgraph)(col, cols[j]);
	}
	
	x[N+1].index = -1;
	
	double pred = svm_predict(model, x);
	
	delete[] x;
	
	return pred;
}
