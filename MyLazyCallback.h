#ifndef MY_LAZY_CALLBACK_H
#define MY_LAZY_CALLBACK_H

#include <ilcplex/ilocplex.h>
#include <vector>
#include "Graph.h"
#include "Util.h"
#include <list>
#include <mutex>

//guarda os vertices de um p3 que não foi adicionado ao modelo. Esse P3 é x[i][j] <= x[i][l] + x[j][l]. Portanto, i < j < l.
typedef struct
{
      int i,j,l;
      int n;
}t_p3;

/** Lazy callback */ 
class MyLazyCallback : public IloCplex::LazyConstraintCallbackI 
{
   public:
 
      //construtor da classe: recebe o environmet, a referencia da variavel e o numero de clusters
      MyLazyCallback(IloEnv env, const IloArray<IloBoolVarArray>& x_ref, const IloBoolVarArray& y_ref, int k, bool use_y = true); 
      
      //metodo que retorna uma copia do calback. Cplex exige este método.
      IloCplex::CallbackI* duplicateCallback() const;
      
      //insere um p3 a lista de p3s nao adicionados no modelo e que devem ser checados pelos LazyCallbacks. Lista é compatilhada por todos os callbacks. Por
      //isso, é acessada de forma mutuamente exclusiva
      void static addLazyP3(int i, int j, int l, int n);

      //no método main do callback está o código que será executado pelo cplex.
      void main(); 

      static double getP3SeparatingTime();

      static double getP3RootSeparatingTime();

      static int getNAddedP3();
      
      static int getNAddedP3Root();

      static int getNCalls();

      static int getNCallsAtRoot();

   private:

      //quando da construcao do objeto desta classe, x é usado para colocar todas as variaveis no vetor x_vars.
      IloArray<IloBoolVarArray> x;
      //x_vars contem as variaveis x. Com esse vetor, peço todos os valores das variaveis x de uma vez só para o CPLEX. Isso é muito mais rápido que pedir cada valor
      //de uma vez.
      IloNumVarArray x_vars;

      IloBoolVarArray y;
      IloNumVarArray y_vars;

      int k;
      
      //quando o callback é usado no hibrido, nao ha variaveis y. Esse flag diz se o lazy deve considerar ou nao variaveis y.
      bool use_y;

      static int n_calls; //numero de vezes que o lazyCallback foi executado em toda a arvore
      static int n_calls_root;
      static double p3_root_separating_time;
      static double p3_separating_time;
      static int n_added_p3;
      static int n_added_p3_root;

      //lista estatica que contem os p3s que ainda nao foram adicionados ao modelo. Todos os p3s sao adicionados de forma global. Por isso,
      //uma lista é suficiente. Ele deve ser acessada de forma exclusiva pelas threads. Isso é feito com o mutex definido abaixo.
      static list<t_p3> notAddedP3;
      static std::mutex lazyMutex;
      
      bool isIntegerFeasible();

      void incNCalls();

      void incNCallsAtRoot();

      void addToP3SeparatingTime(double);

      void addToP3RootSeparatingTime(double);

      static void addToNAddedP3(int n);
      
      static void addToNAddedP3Root(int n);

      //separa P3 atraves de força bruta: enumera todas as inequacoes. Isso leva tempo O(n^3).
      //O parametro allInequalities diz se todas as violacoes serao retornadas. Se for false, apenas uma inequacao por aresta violada é retornada
      std::vector<IloConstraint>* separateP3(int max_p3_edge = -1);
     
}; 

#endif
