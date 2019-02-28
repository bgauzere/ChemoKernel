#include <cstring>
#include <stack>
#include <sstream>
#include "TreeletEnumerator.h"
#include "string_utils.h"
#include "utils.h"
#include "separator.h"
#include "MoleculeGraph.h"

using namespace std;
using namespace pandore;
/**************************************************
 *** Adapté de string utils                        
 ***************************************************/
void get_perm_rec(const pair<string,string> * seq, int n, pair<string,string> ** perm_list, int ** perm_pos_init, 
		  int k, int * perm_pos, int start_perm, int nb_mask, int seq_size)
{
  for(int i=start_perm; i<n; i++)
    {
      pair<string,string>* seq_copy = new pair<string,string>[seq_size];
      for(int ii=0;ii<seq_size;ii++){
	seq_copy[ii].first = seq[ii].first;
	seq_copy[ii].second = seq[ii].second;
      }
      seq_copy[i] = pair<string,string>("NULL","NULL"); //XXX: Fuite memoire
      if ((n-nb_mask-1) == k)     //Permutation #perm_pos get !
	{
	  int pos=0;
	  for(int j=0; j < seq_size; j++) 
	    if(! seq_copy[j].first.compare("NULL"))
	      {
		//j : position dans la chaine de caractere en entrée
		perm_pos_init[(*perm_pos)][pos] = j;
		perm_list[(*perm_pos)][pos] = seq_copy[j];
		pos ++;
	      }
	  (*perm_pos) ++; //Permutation computed !
	}
      else if(n-nb_mask > k) //Pruning
	{
	  get_perm_rec(seq_copy,n, perm_list,perm_pos_init,k,perm_pos,i+1,nb_mask+1,seq_size);
	}
      //XXX:
      //delete [] seq_copy;
    }
}


int get_perm(const pair<string,string> * s, int seq_size, 
	     pair<string,string> ** &perm_list, int ** &perm_pos_init, int k)
{
  int nb_perm =  gsl_sf_choose(seq_size, k);
  /* Mem Alloc */
  perm_list = new pair<string,string>*[nb_perm];
  perm_pos_init = new int*[nb_perm];
  for(int i=0;i<nb_perm;i++)
    {
      perm_pos_init[i] = new int[k];
      perm_list[i] = new pair<string,string>[k];
    }
  //Cas particulier n = k
  if (nb_perm ==1)
    {
      for(int i=0;i<k;i++){
	perm_list[0][i] = pair<string,string>(s[i].first,s[i].second);
	perm_pos_init[0][i] = i;
      }
    }
  else
    {
      // int perm_pos = 0;
      // pair<string,string>* seq = new pair<string,string>[seq_size];
      // seq = (pair<string,string>*)(memcpy(seq,s,seq_size*sizeof(pair<string,string>)));
      // get_perm_rec(seq,seq_size,perm_list,perm_pos_init,k,&perm_pos,0,0,seq_size);
    }
  return nb_perm;
}


string pairToChar(const pair< string,string > * pair_list, int nb_pair){
  stringstream res;
  for(int i=0;i<nb_pair;i++)
    {
      if(i!=0)
	res << LBL_SEP;
      res << pair_list[i].first;
      res << LBL_SEP;
      res << pair_list[i].second;
    }
  return res.str();
}

int compPair(const void * l, const void * r)
{
  int c_edge = ((pair<string,string> *)l)->first.compare(((pair<string,string> *)r)->first);
  if(c_edge == 0)
    return ((pair<string,string> *)l)->second.compare(((pair<string,string> *)r)->second);
  else 
    return c_edge;
}

pair<string, string> * pairSort(const pair<string, string> * to_sort, int nb_pair)
{
  pair<string, string> * sort_pair = new pair<string, string>[nb_pair];
  for(int i=0;i<nb_pair;i++){
    sort_pair[i].first = to_sort[i].first;
    sort_pair[i].second = to_sort[i].second;
  }
  qsort (sort_pair,nb_pair,sizeof(pair<string, string>),compPair);
  return sort_pair;  
}

// //XXX: a recoder avec les paires
// char * star_code(char center, char * neighbours, int type)
// {
//   char * result = new char[type + 2];//+2 centre + '\0'
//   result[0] = center;
//   result[1] = '\0';
//   neighbours = string_utils::stringSort(neighbours);
//   result = strcat(result,neighbours);
//   result[type+1] = '\0';
//   return result;
// }

void TreeletEnumerator::countGraphlet(string code, int graphlet,double count_value, treelet_spectrum ** spectre)
{
  assert(code.length() > 0);
  //Les chemins sont stockes sur les SIZE_MAX premieres entrées du spectre
  map<string, double>::iterator it;
  if((it = spectre[graphlet]->find(code)) 
     != spectre[graphlet]->end())
    /*Incrémentation*/
    it->second += count_value;
  else
    /*Creation*/
    spectre[graphlet]->insert(pair<string, double>(code,count_value));
}

/*
 * Enumerating all 3-stars based treelets : G6, G7, G9, G10, G12
 */
void TreeletEnumerator::enumerate3StarsGraphs(GNode<Point3d> * star_center, int nb_neighbours,
					      treelet_spectrum ** spectre, 
					      pandore::Graph3d& g, 
					      vector<string> nodes,
					      vector<string> edges)
{
  
  string center_label = nodes[star_center->Item()];
  //liste des voisins de la 3-star
  //pair<aretes,noeuds>
  pair<string,string> * neighbours = new pair<string,string>[nb_neighbours];

  vector<int> item_neighbours(nb_neighbours);
  int n=0;    
  GEdge * v = star_center->Neighbours();
  while(v != NULL)
    {
      item_neighbours[n] = v->Node();
      // neighbours[n].first = edges[v->Item()];
      // neighbours[n].second = nodes[v->Node()];
      neighbours[n] = pair<string,string>(edges[v->Item()],
					  nodes[v->Node()]);
      n++;
      v = v->Next();
    }
  
  //On recupere toutes les 3 star
  pair<string,string> ** stars;
  int ** stars_items;//XXX:Nom de variable a revoir
  pair<string,string> null_element("","");
  int nb_3star = utils::get_perm(neighbours,nb_neighbours,stars,stars_items,3, null_element);
  for(int i=0; i< nb_3star; i++)  
    {//i-eme 3 star 
      //G6 Get !
      bool * flags_visited = new bool[nodes.size()];
      memset(flags_visited,0,sizeof(bool)*nodes.size());
      stringstream g6_code;
      g6_code << center_label;
      g6_code << LBL_SEP;
      pair<string,string> * sort_pairs = pairSort(stars[i],3);
      g6_code << pairToChar(sort_pairs, 3);
      flags_visited[item_neighbours[stars_items[i][0]]] = 
	flags_visited[item_neighbours[stars_items[i][1]]] = 
      	flags_visited[item_neighbours[stars_items[i][2]]] = 
	flags_visited[star_center->Item()] = true;
      countGraphlet(g6_code.str(),6,1,spectre);
      vector<int> items_g9_leaves(3);
      int nb_g9_leaves = 0;
      //Indexes of current 3-star neighboorhood are in stars_items[i]
      for(int j=0;j<3;j++) //Leaf traversal
      	{
      	  int valence = utils::nbNeighbours(item_neighbours[stars_items[i][j]],g);
      	  if (valence > 1)
      	    {
      	      items_g9_leaves[nb_g9_leaves] = stars_items[i][j];
      	      nb_g9_leaves ++;
      	      //G7 : 3 star + leaf (Factorisable avec G11)
      	      GEdge * v = g[item_neighbours[stars_items[i][j]]]->Neighbours();
      	      while(v != NULL)
      		{
      		  if((g[v->Node()] != star_center) && (flags_visited[v->Node()] ==0))
      		    {//G7 Get ! 1 permutation
      		      stringstream g7_code;
      		      g7_code << nodes[v->Node()];
      		      //current neighbour of 3-star
      		      g7_code << LBL_SEP << edges[v->Item()] << LBL_SEP << stars[i][j].second;
      		      g7_code << LBL_SEP << stars[i][j].first << LBL_SEP << center_label; //leaf	      
      		      //Permutations :
      		      if(compPair((void*)(&stars[i][(j+1)%3]), 
      				  (void*)(&stars[i][(j+2)%3])) < 0 )
      			{
      			  g7_code << LBL_SEP << stars[i][(j+1)%3].first << LBL_SEP << stars[i][(j+1)%3].second;
      			  g7_code << LBL_SEP << stars[i][(j+2)%3].first << LBL_SEP << stars[i][(j+2)%3].second;
      			}
      		      else
      			{
      			  g7_code << LBL_SEP << stars[i][(j+2)%3].first << LBL_SEP << stars[i][(j+2)%3].second;
      			  g7_code << LBL_SEP << stars[i][(j+1)%3].first << LBL_SEP << stars[i][(j+1)%3].second;
      			}
      		      countGraphlet(g7_code.str(),7,1,spectre);
		      
      		      //Now, looking for G10
      		      int valence_v = utils::nbNeighbours(v->Node(),g);
      		      if (valence_v > 1)
      		       	{
      			  //G10 : leaf of a G7 leaf
      		      	  GEdge * v_v = g[v->Node()]->Neighbours();
      		      	  while(v_v != NULL)
      		      	    {
      			      if(flags_visited[v_v->Node()] == 0)
      				{
      				  //G10 get !
      				  stringstream g10_code;
      				  g10_code << nodes[v_v->Node()]; //leaf
      				  g10_code << LBL_SEP;
      				  g10_code << edges[v_v->Item()]; // <<<<<<<<<<<
      				  g10_code << LBL_SEP;
      				  g10_code << g7_code.str();
      				  countGraphlet(g10_code.str(),10,1,spectre);
      				}
      			      v_v = v_v->Next();
      			    }
      			}
      		    }
      		  v = v->Next();
		}//End G7 and G10 search
	      //Looking for G12 treelets, second neighbours traversal(XXX: to factorize)
	      if(valence > 2)
	      	{ //3 star leaf traversal 
	      	  //valence-1 : centre de la 3 star non pris en compte
	      	  pair<string,string> * leaf_neighbours = new pair<string,string>[valence-1]; 
		  int * leaf_item_neighbours = new int[valence];
	      	  n = 0;
	      	  GEdge * v = g[item_neighbours[stars_items[i][j]]]->Neighbours();
	      	  while(v != NULL)
	      	    {
	      	      if((g[v->Node()] != star_center) && (flags_visited[v->Node()] == 0))
	      		{
	      		  leaf_item_neighbours[n] = v->Node();
	      		  leaf_neighbours[n] = pair<string, string>(edges[v->Item()],
								    nodes[v->Node()]); 
	      		  n++;
	      		}
	      	      v = v->Next();
	      	    }
		  
	      	  pair<string,string> ** leaf_perm;
	      	  int ** leaf_perm_item;
		  if(n>=2){
		    int nb_g12 = utils::get_perm(leaf_neighbours,n,leaf_perm, leaf_perm_item,2,null_element);
		    for(int g12_perm=0;g12_perm<nb_g12;g12_perm++)
		      {
			/*    v1        v5
			 *      \      /
			 *       v0--v3
			 *      /      \
			 *     v2       v4
			 *
			 *leaf_perm[g12_perm]                    -> v1 et v2
			 *stars[i][j]                            -> v0, 
			 *star_center                            -> v3, 
			 *stars[i][(j+1)%3] et stars[i][(j+2)%3] -> v4 et v5 */
			stringstream g12_code_1;//left side
			g12_code_1 << stars[i][j].second;
			if(compPair((void*)(&leaf_perm[g12_perm][0]), (void*)(&leaf_perm[g12_perm][1])) < 0)
			  g12_code_1 << LBL_SEP <<
			    leaf_perm[g12_perm][0].first << LBL_SEP << 
			    leaf_perm[g12_perm][0].second << LBL_SEP << 
			    leaf_perm[g12_perm][1].first << LBL_SEP << 
			    leaf_perm[g12_perm][1].second;
			else
			  g12_code_1 << LBL_SEP <<
			    leaf_perm[g12_perm][1].first << LBL_SEP << 
			    leaf_perm[g12_perm][1].second << LBL_SEP << 
			    leaf_perm[g12_perm][0].first << LBL_SEP << 
			    leaf_perm[g12_perm][0].second;
		      
			stringstream g12_code_2;//right side
			g12_code_2 << center_label;
			if(compPair((void*)(&stars[i][(j+1)%3]), (void*)(&stars[i][(j+2)%3])) < 0)
			  g12_code_2 << LBL_SEP <<
			    stars[i][(j+1)%3].first << LBL_SEP << 
			    stars[i][(j+1)%3].second << LBL_SEP <<
			    stars[i][(j+2)%3].first << LBL_SEP <<
			    stars[i][(j+2)%3].second;
			else
			  g12_code_2 << LBL_SEP <<
			    stars[i][(j+2)%3].first << LBL_SEP << 
			    stars[i][(j+2)%3].second << LBL_SEP <<
			    stars[i][(j+1)%3].first << LBL_SEP <<
			    stars[i][(j+1)%3].second;

			stringstream g12_code;
			string left_str = g12_code_1.str();
			string right_str = g12_code_2.str();
			if (left_str.compare(right_str) < 0)
			  g12_code << left_str << LBL_SEP << stars[i][j].first << LBL_SEP << right_str;
			else
			  g12_code << right_str << LBL_SEP << stars[i][j].first << LBL_SEP << left_str;
			countGraphlet(g12_code.str(),12,0.5,spectre);//0.5 => G12 found by two 3-stars
		      
		      }
		    delete [] leaf_neighbours;
		    delete [] leaf_item_neighbours;
		    //XXX: seg fault on aids_riesen
		    // for(int g12_perm=0;g12_perm<nb_g12;g12_perm++){
		    //   delete [] leaf_perm_item[g12_perm];
		    //   delete [] leaf_perm[g12_perm];
		    // }
		    delete [] leaf_perm;
		    delete [] leaf_perm_item;
		    
		  }//End G12 search
		}
	    }
	}//End Leaf traversal
      /*Looking for G9
       *G9 :
       *         v1 
       *          |
       * v3--v2--v0--v4--v5 
       *
       */
      if(nb_g9_leaves > 1)
      	{
      	  for(int j=0;j<3;j++) //Leaf traversal
      	    {//On se positionne sur v1
      	      int index_v2 = item_neighbours[stars_items[i][(j+1) % 3]];
	      int valence_v2 = utils::nbNeighbours(index_v2,g);		
      	      int index_v4 = item_neighbours[stars_items[i][(j+2) % 3]];
      	      int valence_v4 = utils::nbNeighbours(index_v4,g);
      	      if((valence_v2 > 1) && (valence_v4 > 1))
      		{
      		  for(GEdge * v_v2 = g[index_v2]->Neighbours();v_v2!=NULL;v_v2=v_v2->Next())
		    //TODO : test flags_visited[v_v2->Node()] == 0 should be enough : To check
      		    if((g[v_v2->Node()] != star_center) && (flags_visited[v_v2->Node()] == 0))
		      for(GEdge * v_v4 = g[index_v4]->Neighbours();v_v4 != NULL;v_v4 = v_v4->Next())
			if((g[v_v4->Node()] != star_center) && (flags_visited[v_v4->Node()] == 0)
			   && (v_v4->Node() != v_v2->Node())) {
			  //G9 Get !
			  stringstream g9_code;
			  g9_code << center_label << LBL_SEP
				  << stars[i][j].first << LBL_SEP
				  << stars[i][j].second;
			  
			  stringstream g9_left_side;
			  g9_left_side << stars[i][(j+1) % 3].first << LBL_SEP
				       << stars[i][(j+1) % 3].second << LBL_SEP
				       << edges[v_v2->Item()] << LBL_SEP
				       << nodes[v_v2->Node()];
			  
			  stringstream g9_right_side;
			  g9_right_side << stars[i][(j+2) % 3].first << LBL_SEP
					<< stars[i][(j+2) % 3].second << LBL_SEP
					<< edges[v_v4->Item()] << LBL_SEP
					<< nodes[v_v4->Node()];
			  string left_str = g9_left_side.str();
			  string right_str = g9_right_side.str();
			  if(left_str.compare(right_str) < 0)
			    g9_code << LBL_SEP << left_str
				    << LBL_SEP << right_str;
			  else
			    g9_code << LBL_SEP << right_str
				    << LBL_SEP << left_str;
			  countGraphlet(g9_code.str(),9,1.0,spectre);
			}
		}
      	    }
      	}
      //XXX: Invalid delete : 
      //delete [] sort_pairs;
      delete [] flags_visited;
    }//Go for next 3 star
  //Memory free
  // delete [] neighbours;
  // delete [] item_neighbours;
  for (int i=0;i<nb_3star;i++){
    delete [] stars[i];
    delete [] stars_items[i];
  }
  delete [] stars;
  delete [] stars_items; 
}


void TreeletEnumerator::enumerate4StarsGraphs(GNode<Point3d> * star_center, int nb_neighbours,
					  treelet_spectrum ** spectre,
					  Graph3d& g, 
					  vector<std::string> nodes,
					  vector<std::string> edges)
{
  string center_label = nodes[star_center->Item()];
  //liste des voisins de la 4-star
  pair<string,string> * neighbours = new pair<string,string>[nb_neighbours];
  int * item_neighbours= new int[nb_neighbours]; 
  int n=0;
  GEdge * v = star_center->Neighbours();
  while(v != NULL)
    {
      item_neighbours[n] = v->Node();
      neighbours[n] = pair<string,string>(edges[v->Item()],
					  nodes[v->Node()]);
      n++;
      v = v->Next();
    }
  //Next recupere toutes les 4 star
  pair<string,string> ** stars;
  int ** stars_items;
  pair<string,string> null_element("","");
  int nb_4star = utils::get_perm(neighbours,nb_neighbours,stars,stars_items,4,null_element);
  for(int i=0; i<nb_4star; i++)  
    {
      //G8 Get !
      stringstream g8_code;
      bool * flags_visited = new bool[nodes.size()];
      memset(flags_visited,0,sizeof(bool)*nodes.size());
      g8_code << center_label << LBL_SEP;
      pair<string,string> * sort_pairs = pairSort(stars[i],4);
      g8_code << pairToChar(sort_pairs, 4);
      flags_visited[item_neighbours[stars_items[i][0]]] = 
	flags_visited[item_neighbours[stars_items[i][1]]] = 
      	flags_visited[item_neighbours[stars_items[i][2]]] = 
	flags_visited[item_neighbours[stars_items[i][3]]] = 
	flags_visited[star_center->Item()] = true;
      countGraphlet(g8_code.str(),8,1,spectre);
      //Indexes of current 4-star neighboorhood are in stars_items[i]
      for(int j=0;j<4;j++) //Leaf traversal
  	{
  	  int valence = utils::nbNeighbours(item_neighbours[stars_items[i][j]],g);
  	  if (valence > 1)
  	    {//G11 : 4 star + leaf (Factorisable avec G7)
  	      GEdge * v = g[item_neighbours[stars_items[i][j]]]->Neighbours();
  	      while(v != NULL)
  		{
  		  if((g[v->Node()] != star_center) && (flags_visited[v->Node()] == 0))
  		    {
  		      //G11 Get ! 1 permutation
		      stringstream g11_code;
		      g11_code << nodes[v->Node()] << LBL_SEP
			       << edges[v->Item()] << LBL_SEP
			       << stars[i][j].second << LBL_SEP //current neighbour of 4-star
			       << stars[i][j].first << LBL_SEP
			       << center_label; //4-star center ;

  		      //Perms : (other leaves : j+1,j+2,j+3)
  		      pair<string,string> * g11_leaves = new pair<string,string>[3];
  		      g11_leaves[0] = stars[i][(j+1)%4];
  		      g11_leaves[1] = stars[i][(j+2)%4];
  		      g11_leaves[2] = stars[i][(j+3)%4];
		      g11_code << LBL_SEP << pairToChar(pairSort(g11_leaves,3),3);
  		      delete [] g11_leaves;
  		      countGraphlet(g11_code.str(),11,1,spectre);
  		    }
  		  v = v->Next();
  		}
  	    }
  	}
      delete [] flags_visited;
    }
  delete [] neighbours;
  delete [] item_neighbours;
  for (int i=0;i<nb_4star;i++){
    delete [] stars[i];
    delete [] stars_items[i];
  }
  delete [] stars;
  delete [] stars_items; 
}

void TreeletEnumerator::enumerate5Star(GNode<Point3d> * star_center, int nb_neighbours,
				       treelet_spectrum ** spectre,
				       Graph3d& g, 
				       vector<std::string> nodes,
				       vector<std::string> edges)
{
  string center_label = nodes[star_center->Item()];
  //liste des voisins de la 5-star
  pair<string,string> * neighbours = new pair<string,string>[nb_neighbours];
  int * item_neighbours = new int[nb_neighbours]; 
  int n=0;    
  GEdge * v = star_center->Neighbours();
  while(v != NULL)
    {
      item_neighbours[n] = v->Node();      
      neighbours[n] = pair<string,string>(edges[v->Item()],
					  nodes[v->Node()]);
      n++;
      v = v->Next();
    }
  //On recupere toutes les 5 star
  pair<string,string> ** stars;
  int ** stars_items;
  pair<string,string> null_element("","");
  int nb_5star = utils::get_perm(neighbours,nb_neighbours,stars,stars_items,5,null_element);
  for(int i=0; i< nb_5star; i++)  
    {
      //G13 Get!
      stringstream g13_code;
      g13_code << center_label << LBL_SEP;

      g13_code << pairToChar(pairSort(stars[i],5),5);
      countGraphlet(g13_code.str(),13,1,spectre);
    }
  delete [] neighbours;
  delete [] item_neighbours; 
  for (int i=0;i<nb_5star;i++){
    delete [] stars[i];
    delete [] stars_items[i];
  }
  delete [] stars;
  delete [] stars_items; 
}

treelet_spectrum ** TreeletEnumerator::computeTreeletSpectrum(Graph3d& g, 
							      vector<string> nodes,
							      vector<string> edges)
{  
  bool(*fn_pt)(string,string) = string_utils::keyStringComp;
  treelet_spectrum ** spectre = new treelet_spectrum*[SIZE_SPECTRUM];
  for(int i = 0; i < SIZE_SPECTRUM; i++)
    spectre[i] = new treelet_spectrum(fn_pt);
  
  unsigned int n = g.Size();
  //On parcourt tous les noeuds du graphe
  for(unsigned int i=0;i<n;i++)
    {
      //Chemins de taille 6 max a partir du noeud i
      GNode<Point3d> * current_node = g[i];
      if(g[i] !=NULL){//Deleted Node
	bool * flags_visited = new bool[nodes.size()];
	memset(flags_visited,0,sizeof(bool)*nodes.size());
	string init("");
	enumerateKPaths(current_node,SIZE_MAX,init, NULL, spectre, g,nodes,edges,flags_visited);
	delete [] flags_visited;

	int nb_neighbours = utils::nbNeighbours(i,g);
	if (nb_neighbours >= 3) /*3 Star*/
	  {
	    enumerate3StarsGraphs(current_node,nb_neighbours,spectre,g,nodes,edges);
	  }
	if (nb_neighbours >= 4) /*4 Star*/
	  {
	    enumerate4StarsGraphs(current_node, nb_neighbours,spectre,g,nodes,edges);
	  }
	if(nb_neighbours >= 5) /*5 Star*/
	  {
	    enumerate5Star(current_node, nb_neighbours,spectre,g,nodes,edges);
	  }
      }
    }
  
  return spectre;
}

string TreeletEnumerator::reverse_code(std::string code){
  char * tokens = new char[code.length()+1];
  tokens = strcpy(tokens,code.c_str());
  char * token = strtok (tokens,LBL_SEP);
  stack<string> s;

  while (token != NULL)
    {
      s.push(string(token));
      token = strtok (NULL, LBL_SEP);
    }

  string reversed("");
  while(!s.empty()){
    if(!reversed.empty())
      reversed.append(LBL_SEP);
    reversed.append(s.top());
    s.pop();
  }
  delete [] tokens;
  return reversed;
}

void TreeletEnumerator::enumerateKPaths(GNode<Point3d> * current_node, int size, 
					string path_code, GNode<Point3d> * parent_node,
					treelet_spectrum ** spectre,
					pandore::Graph3d& g, 
					std::vector<std::string> nodes,
					std::vector<std::string> edges,
					bool * flags_visited){
  if(size > 0)
    {
      string current_label = nodes[current_node->Item()];
      if(! path_code.empty())
	path_code = path_code.append(LBL_SEP);//"LBL_SEP" separate the labels
      path_code = path_code.append(current_label);
      flags_visited[current_node->Item()] = 1;
      
      if(SIZE_MAX-size == 0)
  	countGraphlet(path_code, 0,1,spectre);// Chemin de longueur 1 traversé 1 fois
      else{//Compute lowest Canonical code
	string path_code_reverse = TreeletEnumerator::reverse_code(path_code);
	string canonical_code = ((path_code.compare(path_code_reverse))< 0? 
				 path_code:path_code_reverse);
	countGraphlet(canonical_code, SIZE_MAX-size,0.5,spectre);
      }
      
      GEdge * v = current_node->Neighbours();
      while(v != NULL)
  	{
  	  if((parent_node == NULL) || ((unsigned long)v->Node() != parent_node->Item()))
  	    {
	      if(!flags_visited[v->Node()])
		{
		  string copy_path = string(path_code);
		  copy_path.append(LBL_SEP);//Ajout de la liaison
		  copy_path.append(edges[v->Item()]);
		  bool * flags_visited_copy = new bool[nodes.size()];
		  memcpy(flags_visited_copy,flags_visited,sizeof(bool)*nodes.size());
		  enumerateKPaths(g[v->Node()],size-1,copy_path,current_node,
				  spectre,g,nodes,edges,flags_visited_copy);
		  delete [] flags_visited_copy;
		}
  	    }
  	  v = v->Next();
  	}
    }
}

// void MoleculeGraph::normalizeLabeledSpectrum()
// {
//   double nb_graphlets = 0.0;
//   for(int i = 0; i < SIZE_SPECTRUM; i++){
//     for(spectrum_map::iterator it = labeled_spectrum[i]->begin(); it != labeled_spectrum[i]->end();it++)
//       nb_graphlets += it->second;
//   }
  
//   for(int i = 0; i < SIZE_SPECTRUM; i++){
//     for(spectrum_map::iterator it = labeled_spectrum[i]->begin(); it != labeled_spectrum[i]->end();it++)
//       {
// 	it->second /= nb_graphlets;
//       }
//   }
// }

// std::map<char *, double, bool (*)(char *, char *)> ** MoleculeGraph::getLabeledTreeletSpectrum(){
//   return labeled_spectrum;
// }

  

int TreeletEnumerator::getNbAtoms(int treelet_type){
  if(treelet_type < 6)//Chemins
    return treelet_type + 1;
  if (treelet_type == 6) 
    return 4;
  if (treelet_type > 8)
    return 6;
  else
    return 5;
}

void createLinearTreelet(int length,string code,Graph3d * g, Collection**  nodes){
  
  for(int i=0;i<length-1;i++)
    g->Link(i, i+1, i, 1.0f,false);
  
}

void createNStar(int n,string code,Graph3d * g, Collection**  nodes){
  for(int i=1;i<(n+1);i++){
    g->Link(0, i, i-1, 1.0f, false);
  }
}

void createG7(int length,string code,Graph3d * g, Collection**  nodes){
  g->Link(0,1,0, 1.0f, false);
  g->Link(1,2,1, 1.0f, false);
  g->Link(2,3,2, 1.0f, false);
  g->Link(2,4,3, 1.0f, false);
}

void createG9(int length,string code,Graph3d * g, Collection**  nodes){
  g->Link(0,1,0, 1.0f, false);
  g->Link(0,2,1, 1.0f, false);
  g->Link(2,3,2, 1.0f, false);
  g->Link(0,4,3, 1.0f, false);
  g->Link(4,5,4, 1.0f, false);
}

void createG10(int length,string code,Graph3d * g, Collection**  nodes){
  g->Link(0,1,0, 1.0f, false);
  g->Link(1,2,1, 1.0f, false);
  g->Link(2,3,2, 1.0f, false);
  g->Link(3,4,3, 1.0f, false);
  g->Link(3,5,4, 1.0f, false);
}

void createG11(int length,string code,Graph3d * g, Collection**  nodes){
  g->Link(0,1,0, 1.0f, false);
  g->Link(1,2,1, 1.0f, false);
  g->Link(2,3,2, 1.0f, false);
  g->Link(2,4,3, 1.0f, false);
  g->Link(2,5,4, 1.0f, false);
}
void createG12(int length,string code,Graph3d * g, Collection**  nodes){
  g->Link(0,1,0, 1.0f, false);
  g->Link(0,2,1, 1.0f, false);
  g->Link(0,3,2, 1.0f, false);
  g->Link(3,4,3, 1.0f, false);
  g->Link(3,5,4, 1.0f, false);
}

Collection * TreeletEnumerator::TreeletToCollection(int treelet_type, string code){
  Graph3d * graph = new Graph3d(false);
  int nb_atoms = getNbAtoms(treelet_type);
  graph->New(nb_atoms,0,0,0);
  pandore::Collection** nodes = new Collection*[nb_atoms];
  pandore::Collection* collection = new Collection;
  //Coloration des noeuds
  for(int i=0;i<nb_atoms;i++){
    nodes[i] = new Collection;
    nodes[i]->SETVALUE("atom", Char, code[4*i]);//C_1_C_1_C ...
    graph->Add(i,i);
  }
  //Coloration des aretes
  int nb_bonds = getNbAtoms(treelet_type)-1 ;
  char * edges = new char[nb_bonds];
  for(int i = 0;i<nb_bonds;i++)
    edges[i] = code[2+(4*i)];
  if(treelet_type < 6){
    createLinearTreelet(treelet_type+1,code,graph,nodes);
  } else {
    switch(treelet_type){
    case 6:
      createNStar(3,code,graph,nodes);
      break;
    case 7:
      createG7(5,code,graph,nodes); 
      break;
    case 8:
      createNStar(4,code,graph,nodes);
      break;
    case 9:
      createG9(6,code,graph,nodes);
      break;
    case 10:
      createG10(6,code,graph,nodes);
      break;
    case 11:
      createG11(6,code,graph,nodes);
      break;
    case 12:
      createG12(6,code,graph,nodes);
      break;
    case 13:
      createNStar(5,code,graph,nodes);
      break;
    }
  }
  collection->SETARRAY("edges", Char, edges, nb_bonds);
  collection->SETPARRAY("nodes", Collection, nodes, nb_atoms);
  collection->SETPOBJECT("graph", Graph3d, graph);
  return collection;

}

