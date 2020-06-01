#include <ilcplex/ilocplex.h>     
#include "NodeInfo.h" 
#include "MyBranchCallback.h"
#include "Graph.h"
#include <iostream>

//construtor da classe BrachCallBack
MyBranchCallback::MyBranchCallback(IloEnv env) : IloCplex::BranchCallbackI(env)  
{
}   

//retorna uma cópia do callback. Este método é uma exigência do cplex.
IloCplex::CallbackI* MyBranchCallback::duplicateCallback() const 
{ 
   return new (getEnv()) MyBranchCallback(getEnv()); 
}   

//codigo do callback que será executado pelo cplex. Esse BranchCallback não muda os branchings que o cplex iria fazer caso ele não existisse.
//O propósito deste callback é apenas colocar o NodeInfo em cada do nó da árvore. Isso é feito no momento do branching.
void MyBranchCallback::main() 
{
   // How many branches would CPLEX create? 
   IloInt const nbranch = getNbranches(); 

   if (nbranch > 0) 
   { 
      // CPLEX would branch. Get the branches CPLEX would create 
      // and create exactly those branches. With each branch store 
      // its NodeInfo. 
      // Note that getNodeData() returns NULL for the root node. 
      NodeInfo *data = dynamic_cast<NodeInfo *>(getNodeData()); 
      unsigned int depth;
   
      if(!data)
      {
	 //se entrou aqui é porque o nó atual é a raiz. Cplex não guarda NodeInfo para a raiz, por isso o método acima retorna NULL.
	 //O NodeInfo da raiz é um objeto estático da classe NodeInfo. Abaixo, verifica-se se objeto já foi construído e aponta-se o ponteiro data para ele.
	 if(NodeInfo::rootData == NULL)
	 {
	    NodeInfo::initRootData();
	 }
	 data = NodeInfo::rootData;

	 //so me interessa guardar o valor da FO na raiz. Por isso o codigo abaixo esta dentro do if(!data)
	 data->setNodeObjValue(getObjValue());

	 //calcula quanto tempo foi gasto na raiz
	 struct timeval tempo_end_root;
	 gettimeofday(&tempo_end_root, NULL);
	 time_at_root = calculaTempo(&tempo_start_root, &tempo_end_root);

	 MyBranchCallback::branched = true;

	 //cout << "root will branch with obj value = " << data->getNodeObjValue() << endl;
	 root_gap = this->getMIPRelativeGap();
	 //abort();
	 //return;
      }

      depth = data->getDepth();

      /*
      //faz o branching na maior variavel y diferente de 1
      if(this->y_branching)
      {
	 double y_value = 0;
	 int y_index = -1;
	 for(int i = 1; i < Graph::getOrder(); i++)
	 {
	    if(getValue(y[i]) < 0.99 and getValue(y[i]) > y_value)
	    {
	       y_value = getValue(y[i]);
	       y_index = i;
	    }
	 }

	 //fazer branching na maior variavel y
	 if(y_index != -1)
	 {
	    makeBranch(y[y_index], 1, IloCplex::BranchUp, this->getObjValue(), new NodeInfo(depth + 1U));
	    makeBranch(y[y_index], 0, IloCplex::BranchDown, this->getObjValue(), new NodeInfo(depth + 1U));
	    return;
	 }
      }
      */

      IloNumVarArray vars(getEnv());
      IloNumArray bounds(getEnv()); 
      IloCplex::BranchDirectionArray dirs(getEnv());   

      //cria as mesmas branches que o cplex criaria, mas coloca o NodeInfo.
      for(IloInt i = 0; i < nbranch; ++i) 
      { 
	 IloNum const est = getBranch(vars, bounds, dirs, i);
	 makeBranch(vars, bounds, dirs, est, new NodeInfo(depth + 1U)); 
      } 

      dirs.end(); 
      bounds.end(); 
      vars.end(); 
   } 

   else 
   { 
      // CPLEX would not create any branch here. Prune the node. 
      prune(); 
   } 
}

double MyBranchCallback::getRootGap()
{
   return root_gap;
}

double MyBranchCallback::getTimeAtRoot()
{
   return time_at_root;
}

bool MyBranchCallback::hasBranched()
{
   return branched;
}

bool MyBranchCallback::branched = false;
double MyBranchCallback::time_at_root = -1;
struct timeval MyBranchCallback::tempo_start_root;
double MyBranchCallback::root_gap;
