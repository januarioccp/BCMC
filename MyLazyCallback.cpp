#include "MyLazyCallback.h"
#include "NodeInfo.h"
#include "Graph.h"
#include "Util.h"
#include "FloydWarshall.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <limits.h>
#include <stdlib.h>
#include <stack>
#include <algorithm>
#include <exception>

#define P3_EPS 0.000001
//#define P3_EPS 0.001

//construtor da classe
MyLazyCallback::MyLazyCallback(IloEnv env, const IloArray<IloBoolVarArray>& x_ref, const IloBoolVarArray& y_ref, int k, bool use_y) : IloCplex::LazyConstraintCallbackI(env), x(x_ref), y(y_ref), k(k), use_y(use_y), x_vars(env), y_vars(env)
{
   //preenche x_vars e y_vars
   for(int i = 0; i < Graph::getOrder(); i++)
   {
      this->y_vars.add(y[i]);

      for(int j = i+1; j < Graph::getOrder(); j++)
      {
	 this->x_vars.add(this->x[i][j]);
      }
   }
} 

//retorna uma cópia do callback. Este método é uma exigência do CPLEX.
IloCplex::CallbackI* MyLazyCallback::duplicateCallback() const 
{ 
   return new (getEnv()) MyLazyCallback(getEnv(), x, y, k); 
}

//código do callback que é executado pelo cplex.
void MyLazyCallback::main() 
{
   //o cplex chama o LazyCallback quando a solução do GRASP é adicionada. Essa chamada é a primeira do lazy e deve ser descartada,
   //uma vez que a solução do GRASP é viável.
   if(MyLazyCallback::n_calls == 0)
   {
      MyLazyCallback::incNCalls();
      return;
   }

   //dados do node
   NodeInfo *data = dynamic_cast<NodeInfo*>(getNodeData());

   if(!data)
   {
      //se entrou aqui é porque o no atual é a raiz. Cplex não guarda NodeInfo para a raiz, por isso o método acima retorna NULL.
      //O NodeInfo da raiz é um objeto estático da classe NodeInfo. Abaixo, verifica-se se objeto já foi construído e aponta-se o ponteiro data para ele.
      if(NodeInfo::rootData == NULL)
      {
	 NodeInfo::initRootData();
      }
      data = NodeInfo::rootData;
   }

   int depth = data->getDepth();

   //separa P3 e calcula o tempo de separacao
   struct timeval tempoAntes, tempoDepois;
   gettimeofday(&tempoAntes, NULL);

   int n_p3;

   if(Graph::getOrder() <= 100)
   {
      n_p3 = -1;
   }
   else if(Graph::getOrder() <= 150)
   {
      n_p3 = 3;
   }
   else if(Graph::getOrder() <= 200)
   {
      n_p3 = 2;
   }
   else if(Graph::getOrder() <= 250)
   {
      n_p3 = 1;
   }
   else if(Graph::getOrder() <= 300)
   {
      n_p3 = 1;
   }
   else
   {
      n_p3 = 1;
   }

   std::vector<IloConstraint>* cons = separateP3(n_p3);
  
   gettimeofday(&tempoDepois, NULL);
   double p3_sep_time = calculaTempo(&tempoAntes, &tempoDepois);

   MyLazyCallback::addToP3SeparatingTime(p3_sep_time);
   MyLazyCallback::addToNAddedP3(cons->size());
   if(depth == 0)
   {
      MyLazyCallback::addToP3RootSeparatingTime(p3_sep_time);
      MyLazyCallback::addToNAddedP3Root(cons->size());
   }

   //adiciona as violacoes
   for(int i = 0; i < cons->size(); i++)
   {
      add(cons->at(i)).end();
   }
   delete cons;

   MyLazyCallback::incNCalls();
   if(depth == 0)
   {
      MyLazyCallback::incNCallsAtRoot();
   }
}

//separa P3 enumerando todas as inequacoes. tempo O(n^3)
std::vector<IloConstraint>* MyLazyCallback::separateP3(int max_p3_edge)
{
   std::vector<IloConstraint>* cons = new std::vector<IloConstraint>();

   //pede ao cplex de uma vez so os valores das variaveis relaxadas
   IloNumArray x_vals(getEnv(), (Graph::getOrder()*(Graph::getOrder()-1))/2);
   getValues(x_vals, x_vars);

   //pede ao cplex de uma vez so os valores das variaveis relaxadas
   //IloNumArray y_vals(getEnv(), Graph::getOrder());
   //getValues(y_vals, y_vars);

   //armazena os valores das variaveis relaxadas na matriz m.
   std::vector< std::vector<float> > m(Graph::getOrder(), std::vector<float>(Graph::getOrder()));
   int l = 0;
   for(int i = 0; i < Graph::getOrder(); i++)
   {
      for(int j = i+1; j < Graph::getOrder(); j++)
      {
	 m[i][j] = m[j][i] = x_vals[l++];
      }
   }

   //wasAdded[i][j] = true se uma inequacao ja foi adicionada para a violacao da aresta {i,j}, se ela estiver violada
   std::vector< std::vector<int> > wasAdded(Graph::getOrder(), std::vector<int>(Graph::getOrder(), 0));

   for(int i = 0; i < Graph::getOrder(); i++)
   {
      for(int j = i+1; j < Graph::getOrder(); j++)
      {
	 for(int l = j+1; l < Graph::getOrder(); l++)
	 {
	    if(( 1 + m[i][j] - m[i][l] - m[j][l] < -P3_EPS) and (wasAdded[i][j] < max_p3_edge or max_p3_edge == -1))
	    {
	       //if(false)
	       if(this->use_y)
	       {
		  cons->push_back(x[i][j] + 1 >= x[i][l] + x[j][l] + y[l]);
	       }
	       else
	       {
		  cons->push_back(x[i][j] + 1 >= x[i][l] + x[j][l]);
	       }
	       wasAdded[i][j]++;
	    }

	    if((1 + m[i][l] - m[i][j] - m[j][l] < -P3_EPS) and (wasAdded[i][l] < max_p3_edge or max_p3_edge == -1))
	    {
	       cons->push_back(x[i][l] + 1 >= x[i][j] + x[j][l]);
	       wasAdded[i][l]++;
	    }

	    if((1 + m[j][l] - m[i][j] - m[i][l] < -P3_EPS) and (wasAdded[j][l] < max_p3_edge or max_p3_edge == -1))
	    {
	       cons->push_back(x[j][l] + 1 >= x[i][j] + x[i][l]);
	       wasAdded[j][l]++;
	    }
	 }
      }
   }


   return cons;
}

//incrementa o contador de vezes que o lazy foi chamado em toda a arvore. Codigo é uma seção critica, pois n_calls é uma variavel compartilhada
void MyLazyCallback::incNCalls()
{
   lazyMutex.lock();

   MyLazyCallback::n_calls++;

   lazyMutex.unlock();
}

//incrementa o contador de vezes que o lazy foi chamado na raiz. Codigo é uma seção critica, pois n_calls é uma variavel compartilhada
void MyLazyCallback::incNCallsAtRoot()
{
   lazyMutex.lock();

   MyLazyCallback::n_calls_root++;

   lazyMutex.unlock();
}


//adiciona time à varivel p3_separating_time, que acumula o tempo total de separacao de P3 (em toda a arvore)
void MyLazyCallback::addToP3SeparatingTime(double time)
{
   lazyMutex.lock();

   MyLazyCallback::p3_separating_time += time;

   lazyMutex.unlock();
}
 
//adiciona time à varivel p3_root_separating_time, que acumula o tempo total de separacao de P3 na raiz
void MyLazyCallback::addToP3RootSeparatingTime(double time)
{
   lazyMutex.lock();

   MyLazyCallback::p3_root_separating_time += time;

   lazyMutex.unlock();
}
 
//adiciona n à variável n_added_p3, que acumula o número de p3s adicionados em toda a arvore 
void MyLazyCallback::addToNAddedP3(int n)
{
   lazyMutex.lock();

   MyLazyCallback::n_added_p3 += n;

   lazyMutex.unlock();
}
  
//adiciona n à variável n_added_p3_root, que acumula o número de p3s adicionados na raiz
void MyLazyCallback::addToNAddedP3Root(int n)
{
   lazyMutex.lock();

   MyLazyCallback::n_added_p3_root += n;

   lazyMutex.unlock();
}
 
double MyLazyCallback::getP3SeparatingTime()
{
   lazyMutex.lock();

   double time = MyLazyCallback::p3_separating_time;

   lazyMutex.unlock();

   return time;
}

double MyLazyCallback::getP3RootSeparatingTime()
{
   lazyMutex.lock();

   double time = MyLazyCallback::p3_root_separating_time;

   lazyMutex.unlock();

   return time;
}
 
int MyLazyCallback::getNAddedP3()
{
   lazyMutex.lock();

   int n = MyLazyCallback::n_added_p3;

   lazyMutex.unlock();

   return n;
}
 
int MyLazyCallback::getNAddedP3Root()
{
   lazyMutex.lock();

   int n = MyLazyCallback::n_added_p3_root;

   lazyMutex.unlock();

   return n;
}


int MyLazyCallback::getNCalls()
{
   lazyMutex.lock();

   int n = MyLazyCallback::n_calls;

   lazyMutex.unlock();

   return n;
}

int MyLazyCallback::getNCallsAtRoot()
{
   lazyMutex.lock();

   int n = MyLazyCallback::n_calls_root;

   lazyMutex.unlock();

   return n;
}


void MyLazyCallback::addLazyP3(int i, int j, int l, int n)
{
   t_p3 p3;
   p3.i = i; p3.j = j; p3.l = l; p3.n = n;

   MyLazyCallback::notAddedP3.push_front(p3);
}

//determina se o nó atual é uma solução inteira do NACEP
bool MyLazyCallback::isIntegerFeasible()
{
   for(int i = 0; i < Graph::getOrder(); i++)
   {
      for(int j = i+1; j < Graph::getOrder(); j++)
      {
	 if(getFeasibility(x[i][j]) == Infeasible)
	    return false;
      }
   }
   return true;
}

/*
//separa as violações de p3 quando a solução é inteira.
//Para isso, separa as componentes da solução e procura violações dentro das componentes.
std::vector<IloConstraint>* MyLazyCallback::separateIntegerP3()
{
   std::vector<IloConstraint>* p3s = new std::vector<IloConstraint>();
  
   std::vector<std::vector<int>> edges(Graph::getOrder(), std::vector<int>(Graph::getOrder(), Graph::INACTIVE));

   IloNumArray x_vals(getEnv(), (Graph::getOrder()*(Graph::getOrder()-1))/2);
   getValues(x_vals, x_vars);

   //lista de adjacencias dos vertices
   int n = 0;
   std::vector< std::vector<int> > la(Graph::getOrder());
   for(int i = 0; i < Graph::getOrder(); i++)
   {
      for(int j = i+1; j < Graph::getOrder(); j++)
      {
	 if(x_vals[n++] <= P3_EPS)
	 {
	    la[i].push_back(j);
	    la[j].push_back(i);
	 }
      }
   }

   //component_id[i] guarda o id da componente a qual o vértice i pertence. -1 indica que a busca ainda nao o atingiu
   std::vector<int> component_id(Graph::getOrder(), -1); 

   //realiza DFSs iterativas
   int comp_id = 0; //id da componente que sera construida
   for(int v = 0; v < Graph::getOrder(); v++)
   {
      //vertice v ainda nao foi visitado. começa dfs a partir dele.
      if(component_id[v] == -1)
      {
	 stack<int> pilha_dfs;
	 pilha_dfs.push(v);
	 component_id[v] = comp_id;
	 std::vector<int> component;

	 while(not pilha_dfs.empty())
	 {
	    int u = pilha_dfs.top();
	    pilha_dfs.pop();
	    component.push_back(u);

	    //adiciona os vizinhos de u que ainda nao foram visitados a pilha
	    for(int i : la[u])
	    {
	       //vertice i é vizinho e ainda nao foi visitado na busca
	       if(component_id[i] == -1)
	       {
		  component_id[i] = comp_id;
		  pilha_dfs.push(i);
	       }
	    }
	 }

	 //identifica as violacoes de p3 dentro da componente
	 sort(component.begin(), component.end());
	 for(int i = 0; i < component.size(); i++)
	 {
	    for(int j = i+1; j < component.size(); j++)
	    {
	       for(int l = j+1; l < component.size(); l++)
	       {
		  int a = component[i];
		  int b = component[j];
		  int c = component[l];

		  if(getValue(x[a][b]) - getValue(x[a][c]) - getValue(x[b][c]) > P3_EPS)
		     p3s->push_back(x[a][b] <= x[a][c] + x[b][c]);

		  if(getValue(x[a][c]) - getValue(x[a][b]) - getValue(x[b][c]) > P3_EPS)
		     p3s->push_back(x[a][c] <= x[a][b] + x[b][c]);

		  if(getValue(x[b][c]) - getValue(x[a][b]) - getValue(x[a][c]) > P3_EPS)
		     p3s->push_back(x[b][c] <= x[a][b] + x[a][c]);
	       }
	    }
	 }

	 comp_id++;
      }
   }

	 
   return p3s;
}

//separa as violações de p3 quando a solução é inteira.
//Para isso, roda o algoritmo de floydWarshall para determinar os vertices i,j tais que d(i,j)=2.
std::vector<IloConstraint>* MyLazyCallback::separateIntegerP3Floyd()
{
   std::vector<IloConstraint>* p3s = new std::vector<IloConstraint>();

   //pega os valores das variaves x relaxadas
   IloNumArray x_vals(getEnv(), (Graph::getOrder()*(Graph::getOrder()-1))/2);
   getValues(x_vals, x_vars);

   //preenche a matriz edges, que sera passada para o floyd
   std::vector<std::vector<int>> edges(Graph::getOrder(), std::vector<int>(Graph::getOrder(), Graph::INACTIVE));
   int l = 0;
   for( int i = 0; i < Graph::getOrder(); i++)
   {
      for( int j = i+1; j < Graph::getOrder(); j++)
      {
	 if(x_vals[l++] < 0.1)
	 {
	    edges[i][j] = edges[j][i] = Graph::ACTIVE;
	 }
      }
   }

   FloydWarshall fw(Graph::getOrder());
   //retorna todos os paths com 3 vertices
   std::vector< std::vector<int> >* paths_3 = fw.shortestPaths(&edges);

   for(int i = 0; i < paths_3->size(); i++)
   {
      std::vector<int>& p = (*paths_3)[i];

      int a = p[0];
      int b = p[1];
      int c = p[2];

      //////////////////////////
      //adiciona o p3
      if(a < c)
      {
	 if(a < b)
	 {
	    if(c < b)
	       p3s->push_back(x[a][c] <= x[a][b] + x[c][b]);
	    else
	       p3s->push_back(x[a][c] <= x[a][b] + x[b][c]);
	 }
	 else
	    p3s->push_back(x[a][c] <= x[b][a] + x[b][c]);
      }
      else
      {
	 if(a < b)
	 {
	    if(c < b)
	       p3s->push_back(x[c][a] <= x[a][b] + x[c][b]);
	    else
	       p3s->push_back(x[c][a] <= x[a][b] + x[b][c]);
	 }
	 else
	    p3s->push_back(x[c][a] <= x[b][a] + x[c][b]);
      }
      ////////////////////////
   }
   delete paths_3;
   
   return p3s;
}
*/



list<t_p3> MyLazyCallback::notAddedP3;
mutex MyLazyCallback::lazyMutex;
double MyLazyCallback::p3_separating_time = 0;
double MyLazyCallback::p3_root_separating_time = 0;
int MyLazyCallback::n_calls = 0;
int MyLazyCallback::n_calls_root = 0;
int MyLazyCallback::n_added_p3 = 0;
int MyLazyCallback::n_added_p3_root = 0;
