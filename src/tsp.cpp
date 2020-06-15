// TSP classes and global objects
#include "input.h"
#include "solution.h"
#include "construction.h"
#include "neighborhood.h"
#include "perturbation.h"
#include "localsearch.h"
#include "data.h"

// CPLEX classes and global objects
#include <ilcplex/ilocplex.h>
#include "NodeInfo.h"
#include "MyLazyCallback.h"
#include "MyCutCallback.h"
#include "MyBranchCallback.h"
#define MAX_ITER 100
typedef IloArray<IloBoolVarArray> Edges;

// STL classes
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;

IloInt checkTour(IloNumArray2 sol, IloBoolArray seen, IloNum tol);

int toEdge(int x, int y, int n){
	if(x == y)
		return 0;
	if(x>y)
		return toEdge(y,x,n);
	if( x == 0)
		return y-x;
    else
		return (n-1)*x-(x*(x+1)/2)+y;
}

int main(int argc, char **argv)
{
   // Run RVND to identify an UB
   Input in(argc, argv);   
   Solution sol(&in);
	LocalSearch ls(&in);
   int UB = numeric_limits<int>::max();
   if(argc < 3){
      sol = ls.GILSRVND();
      UB = sol.costValueTSP + 1;
   }
   else{
      UB = atoi(argv[2])+1;
   }

   cout<<"Upper bound |- - - - - - - - - - - - - - - -| "<<UB<<endl;

   // Read input to get instance size and initialize MinCutter
   Data input(argc, argv[1]);
   input.readData();
   IloInt n = input.getDimension();

   //Environment
   IloEnv env;
   env.setName("Branch and Cut");

   // create model
   IloModel tsp(env);
   tsp.setName("Traveling Tournament Problem Model");
   try
   {
      // A matrix for storing the input distances
      IloNumArray2 dist(env, n);

      // You are allocating memory. Remmember to free it before leaving
      for (IloInt i = 0; i < n; i++)
         dist[i] = IloNumArray(env, n);

      for (IloInt i = 0; i < n; i++)
         for (IloInt j = 0; j < n; j++)
            dist[i][j] = input.getDistance(i, j);

      // Create variables x[i][j] for all 0 <= i < j < n representing the
      // edge between cities i and j.  A setting of 1 indicates that the edge
      // is part of the tour.
      Edges x(env, n);
      int k = 1;
      for (IloInt i = 0; i < n; i++)
      {
         x[i] = IloBoolVarArray(env, n);
         for (IloInt j = i+1; j < n; j++){
            string name = "x" + to_string(k++);
            x[i][j].setName(name.c_str());
         }
      }

      // Objective is to minimize the sum of edge weights for traveled edges
      IloExpr obj(env);
      for (IloInt i = 0; i < n; i++)
         for (IloInt j = i+1; j < n; j++)
               obj += x[i][j]*dist[i][j];
      tsp.add(IloMinimize(env, obj));

      // degree constraints --- exactly two traveled edges incident to each node
      for (IloInt i = 0; i < n; i++){
         IloExpr expr(env);
         for (IloInt j = 0; j < n; j++){
            if (i < j)
               expr += x[i][j];
            if (j < i)
               expr += x[j][i];
         }
         IloConstraint c(expr == 2);
         tsp.add(c);
      }

      IloCplex cplex(tsp);

      if (argc == 3)
      {
         // Declaring initial solution
         IloNumArray x_start(env, n * n);
         for (int i = 0; i < x_start.getSize(); i++)
            x_start[i] = 0;

         // Build initial solution
         for (int i = 0; i < sol.location.size() - 1; i++)
            x_start[toEdge(sol.location[i] - 1, sol.location[i + 1] - 1, n)] = 1;

         IloNumVarArray y(env);
         for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
               y.add(x[i][j]);

         cplex.addMIPStart(y, x_start);

         // Write initial solution on a file
         cplex.writeMIPStarts("start.mst");
      }

      IloNum tol = cplex.getParam(IloCplex::EpInt);

      // Subtour Elimination Callback
      MyLazyCallback *lazyCbk = new (env) MyLazyCallback(env,x);
      cplex.use(lazyCbk);

      MyCutCallback *cutCbk = new (env) MyCutCallback(env,x);
      cplex.use(cutCbk);

      MyBranchCallback *branchCbk = new (env) MyBranchCallback(env);
      cplex.use(branchCbk);

      cplex.setParam(IloCplex::PreInd, IloFalse);
      cplex.setParam(IloCplex::TiLim, 2 * 60 * 60);
	   cplex.setParam(IloCplex::Threads, 1);
	   cplex.setParam(IloCplex::CutUp, UB + 1);

      // Export the LP model to a txt file to check correctness
      cplex.exportModel("model.lp");
      bool solved = cplex.solve();
      std::cout << "status: " << cplex.getStatus() << std::endl;

      if (solved)
         env.out() << "Optimal tour length "
                   << cplex.getObjValue() << endl;

      return 0;

      IloNumArray2 sol(env, n);
      for (IloInt i = 0; i < n; i++)
      {
         sol[i] = IloNumArray(env, n);
         cplex.getValues(sol[i], x[i]);
      }
      IloBoolArray seen(env);
      IloInt length = checkTour(sol, seen, tol);

      if (length < n)
      {
         IloExpr clique(env);
         for (int i = 0; i < n; i++)
         {
            if (seen[i])
            {
               for (int j = i + 1; j < n; j++)
               {
                  if (seen[j])
                     clique += x[j][i];
               }
            }
         }
         cerr << cplex.getValue(clique) << " <= " << length - 1 << endl;
      }

      // assert (length == n);

#ifdef FULLTEST
      assert(cplex.getImpl()->isConsistent());
      assert(cpxtest.getImpl()->isConsistent());
      assert(cplex.getStatus() == IloAlgorithm::Optimal);
      assert(fabs(cplex.getObjValue() - 11461.0) < 1e-6);
      assert(cutCalled);
      env.out() << "Test completed successfully" << endl;
#endif

      // sec.end();
      for (IloInt i = 0; i < n; i++)
         dist[i].end();
      dist.end();
   }
   catch (const IloException &e)
   {
      cerr << "Exception caught: " << e << endl;
#ifdef FULLTEST
      assert(0);
#endif
   }
   catch (...)
   {
      cerr << "Unknown exception caught!" << endl;
#ifdef FULLTEST
      assert(0);
#endif
   }

   env.end();
   return 0;
}

IloInt checkTour(IloNumArray2 sol, IloBoolArray seen, IloNum tol)
{
   IloInt j, n = sol.getSize();
   IloInt last = -1;
   IloInt length = 0;
   IloInt current = 0;
   seen.clear();
   seen.add(n, 0.0);

   // Search for a subtour if sol[] is integer feasible

   while (seen[current] == 0)
   {
      cout << current << " ";
      length++;
      seen[current] = length;
      for (j = 0; j < current; j++)
      {
         if (j != last && sol[current][j] >= 1.0 - tol)
            break;
      }
      if (j == current)
      {
         for (j = current + 1; j < n; j++)
         {
            if (j != last && sol[j][current] >= 1.0 - tol)
               break;
         }
      }
      if (j == n)
         return (n + 1);
      last = current;
      current = j;
   }
   cout << endl;
   return (length);
}