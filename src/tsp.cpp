// TSP classes
#include "input.h"
#include "solution.h"
#include "construction.h"
#include "neighborhood.h"
#include "perturbation.h"
#include "localsearch.h"

#define MAX_ITER 100
#include "MinCutter.h"
#include <ilcplex/ilocplex.h>
#include "data.h"
#include "NodeInfo.h"
#include <vector>
#include <cmath>
#include <cstdlib>
using namespace std;

#ifdef FULLTEST
#include <assert.h>
IloBool cutCalled = IloFalse;
#endif

ILOSTLBEGIN

typedef IloArray<IloBoolVarArray> Edges;

MinCutter *mc;

IloInt checkTour(IloNumArray2 sol, IloBoolArray seen, IloNum tol);

ILOUSERCUTCALLBACK2(MinCut, Edges, x, IloNum, tol)
{
   if (getCurrentNodeDepth() > 7)
   {
      return;
   }

   IloEnv env = getEnv();

   // Number of vertices
   IloInt n = x.getSize();

   // Retrieve solution information
   vector<vector<double>> w(n);
   for (IloInt i = 0; i < n; i++)
   {
      w[i].resize(n);
      for (IloInt j = 0; j < n; j++)
         w[i][j] = abs(double(getValue(x[i][j])));
   }

   // Creating an adjacency matrix
   for (IloInt i = 0; i < n; i++)
      for (IloInt j = 0; j < i; j++)
         w[j][i] = w[i][j];

   mc->MINIMUMCUTUPDATE(w);

   // Add the constraint
   if (mc->getMinCut() < 2.0 - tol)
   {
      pair<vector<int>, vector<int>> partition = mc->getPartition();

      IloExpr expr1(env);
      for (auto i : partition.first)
         for (auto j : partition.second)
            expr1 += x[i][j] + x[j][i];
      add(expr1 >= 2).end();
      expr1.end();
   }
}

ILOUSERCUTCALLBACK2(MaxBack, Edges, x, IloNum, tol)
{
   if (getCurrentNodeDepth() > 7)
      return;

   IloEnv env = getEnv();
   IloInt n = x.getSize();

   // Retrieve solution information
   vector<vector<double>> sol(n);
   for (IloInt i = 0; i < n; i++)
   {
      sol[i].resize(n);
      for (IloInt j = 0; j < n; j++)
         sol[i][j] = abs(double(getValue(x[i][j])));
   }

   // Creating an adjacency matrix
   for (IloInt i = 0; i < n; i++)
      for (IloInt j = 0; j < i; j++)
         sol[j][i] = sol[i][j];

   vector<int> S;
   vector<int> Smin;
   vector<bool> seen(n);
   vector<double> b(n);
   fill(seen.begin(), seen.end(), false);

   int v = rand() % n;
   S.push_back(v);
   seen[v] = true;
   b[v] = -numeric_limits<double>::infinity();

   double Cutmin;
   double Cutval;
   Cutmin = 0.0;
   for (int i = 0; i < n; i++)
   {
      if (!seen[i] && v != i)
      {
         b[i] = sol[v][i];
         Cutmin += b[i];
      }
   }

   Smin = S;
   Cutval = Cutmin;

   while (S.size() < n)
   {
      // Choose v not in S of maximum max-back value b(v)
      double maxb = -numeric_limits<double>::infinity();
      for (int i = 0; i < n; i++)
         if (!seen[i])
            if (maxb < b[i])
            {
               maxb = b[i];
               v = i;
            }

      S.push_back(v);
      seen[v] = true;
      Cutval = Cutval + 2 - 2 * b[v];

      for (int t = 0; t < n; t++)
      {
         if (!seen[t] && v != t)
            b[t] = b[t] + sol[v][t];
      }
      if (Cutval < Cutmin)
      {
         Cutmin = Cutval;
         Smin = S;
      }
   }

   if (Cutmin < 2.0 - tol)
   {
      fill(seen.begin(), seen.end(), false);
      vector<int> v = Smin;
      for (int i = 0; i < v.size(); i++)
         seen[v[i]] = true;
      vector<int> y;
      for (int i = 0; i < n; i++)
         if (!seen[i])
            y.push_back(i);
      if (v.size() > 0 && y.size() > 0)
      {
         IloExpr expr1(env);
         for (int i = 0; i < v.size(); i++)
            for (int j = 0; j < y.size(); j++)
               expr1 += x[v[i]][y[j]] + x[y[j]][v[i]];
         // cout << "Adding cut " << expr1 << " >= 2 " << endl;
         add(expr1 >= 2).end();
         expr1.end();
      }
   }
   else{
      mc->MINIMUMCUTUPDATE(sol);
      // Add the constraint
      if (mc->getMinCut() < 2.0 - tol)
      {
         pair<vector<int>, vector<int>> partition = mc->getPartition();

         IloExpr expr1(env);
         for (auto i : partition.first)
            for (auto j : partition.second)
               expr1 += x[i][j] + x[j][i];
         add(expr1 >= 2).end();
         expr1.end();
      }
   }
}

ILOLAZYCONSTRAINTCALLBACK2(SubtourEliminationCallback, Edges, x, IloNum, tol)
{

   IloEnv env = getEnv();
   IloInt n = x.getSize();

   IloNumArray2 sol(env, n);
   for (IloInt i = 0; i < n; i++)
   {
      sol[i] = IloNumArray(env, n);
      getValues(sol[i], x[i]);
   }

   for (IloInt i = 0; i < n; i++)
      for (IloInt j = 0; j < i; j++)
         sol[j][i] = sol[i][j];

   // Declares a boolean vector of size n with false
   vector<bool> seen(n, false);

   // An array to store the subtours
   IloNumArray tour(env, n);

   IloInt i, node, len, start;

   // Vector os positions
   vector<pair<int, int>> p;

   for (i = 0; i < n; i++)
      seen[i] = false;

   start = 0;
   node = 0;

   // Start from position 0 in the tour
   while (start < n)
   {
      // Find a node that has not been seen
      for (node = 0; node < n; node++)
         if (!seen[node])
            break;
      // Did you see every node? Time to break
      if (node == n)
         break;

      // Start with lenght 0 and build a subtour
      for (len = 0; len < n; len++)
      {
         // insert the unseen node in the tour
         tour[start + len] = node;

         //Be honest if you saw that node
         seen[node] = true;

         // Check nodes neighbors
         for (i = 0; i < n; i++)
         {
            // Is it connected to someone?
            // First time you see it? ...
            if (sol[node][i] >= 1.0 - tol && !seen[i])
            {
               // ... better catch that guy
               node = i;
               break;
            }
         }

         // Oh man, you could not find a neighbhor? It seens that you closed the loop
         if (i == n)
         {
            // In this case, increase the size of the lenght
            len++;
            pair<int, int> pos(start, len);
            p.push_back(pos);
            //Start a new subtour
            start += len;
            break;
         }
      }
   }

   // Create and add subtour constraint ---
   // No more than 'length-1' edges between members of the subtour
   if (p.size() > 1)
   {
      for (i = 0; i < p.size(); i++)
      {
         IloExpr expr1(env);
         for (int a = p[i].first; a < p[i].first + p[i].second; a++)
            for (int b = a + 1; b < p[i].first + p[i].second; b++)
               expr1 += x[tour[a]][tour[b]] + x[tour[b]][tour[a]];
         // cout << "Adding lazy constraint " << expr1 << " <= " << p[i].second - 1 << endl;
         add(expr1 <= p[i].second - 1).end();
         expr1.end();
      }
   }

   for (IloInt i = 0; i < n; i++)
      for (IloInt j = 0; j < i; j++)
         sol[j][i] = sol[i][j];

   vector<int> S;
   vector<int> Smin;
   vector<double> b(n);
   fill(seen.begin(), seen.end(), false);

   int v = rand() % n;
   S.push_back(v);
   seen[v] = true;
   b[v] = -numeric_limits<double>::infinity();

   double Cutmin;
   double Cutval;
   Cutmin = 0.0;
   for (int i = 0; i < n; i++)
   {
      if (!seen[i] && v != i)
      {
         b[i] = sol[v][i];
         Cutmin += b[i];
      }
   }

   Smin = S;
   Cutval = Cutmin;

   while (S.size() < n)
   {
      // Choose v not in S of maximum max-back value b(v)
      double maxb = -numeric_limits<double>::infinity();
      for (int i = 0; i < n; i++)
         if (!seen[i])
            if (maxb < b[i])
            {
               maxb = b[i];
               v = i;
            }

      S.push_back(v);
      seen[v] = true;
      Cutval = Cutval + 2 - 2 * b[v];

      for (int t = 0; t < n; t++)
      {
         if (!seen[t] && v != t)
            b[t] = b[t] + sol[v][t];
      }
      if (Cutval < Cutmin)
      {
         Cutmin = Cutval;
         Smin = S;
      }
   }

   if (Cutmin < 2 - tol)
   {
      fill(seen.begin(), seen.end(), false);
      vector<int> v = Smin;
      for (int i = 0; i < v.size(); i++)
         seen[v[i]] = true;
      vector<int> y;
      for (int i = 0; i < n; i++)
         if (!seen[i])
            y.push_back(i);
      if (v.size() > 0 && y.size() > 0)
      {
         IloExpr expr1(env);
         for (int i = 0; i < v.size(); i++)
            for (int j = 0; j < y.size(); j++)
               expr1 += x[v[i]][y[j]] + x[y[j]][v[i]];
         // cout << "Adding cut " << expr1 << " >= 2 " << endl;
         add(expr1 >= 2).end();
         expr1.end();
      }
   }
   seen.end();
   tour.end();
   for (IloInt i = 0; i < n; i++)
      sol[i].end();
   sol.end();
}

int main(int argc, char **argv)
{
   Input in(argc, argv);
   Data input(argc, argv[1]);
   input.readData();
   Solution sol(&in);
	LocalSearch ls(&in);
   sol = ls.GILSRVND();

   int UB = sol.costValueTSP + 1;
   cout<<UB<<endl;

   IloInt n = input.getDimension();
   mc = new MinCutter(n);

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

      // Create variables x[i][j] for all 0 <= j < i < n representing the
      // edge between cities i and j.  A setting of 1 indicates that the edge
      // is part of the tour.
      Edges x(env, n);
      for (IloInt i = 0; i < n; i++)
      {
         x[i] = IloBoolVarArray(env, n);
         for (IloInt j = 0; j < n; j++)
         {
            string name = "x_" + to_string(i) + "_" + to_string(j);
            x[i][j].setName(name.c_str());
         }
      }

      // Objective is to minimize the sum of edge weights for traveled edges
      IloExpr obj(env);
      for (IloInt i = 0; i < n; i++)
      {
         obj += IloScalProd(x[i], dist[i]);
      }
      tsp.add(IloMinimize(env, obj));

      // degree constraints --- exactly two traveled edges incident to each node
      for (IloInt i = 0; i < n; i++)
      {
         IloExpr expr(env);
         for (IloInt j = 0; j < n; j++)
         {
            if (i > j)
               expr += x[i][j];
            if (j > i)
               expr += x[j][i];
         }
         IloConstraint c(expr == 2);
         c.setName("Degree");
         tsp.add(c);
      }

      IloCplex cplex(tsp);
      // Export the LP model to a txt file to check correctness
      cplex.exportModel("model.lp");

      IloNum tol = cplex.getParam(IloCplex::EpInt);

      cplex.use(SubtourEliminationCallback(env, x, tol));
      cplex.use(MaxBack(env, x, tol));
      //      cplex.use(MinCut(env, x, tol));
      cplex.setParam(IloCplex::PreInd, IloFalse);
      cplex.setParam(IloCplex::TiLim, 2 * 60 * 60);
	   cplex.setParam(IloCplex::Threads, 1);
	   cplex.setParam(IloCplex::CutUp, UB + 1);
      // cplex.setParam(IloCplex::Param::MIP::Cuts::FlowCovers, -1);
      // cplex.setParam(IloCplex::Param::MIP::Cuts::Gomory, -1);
      // cplex.setParam(IloCplex::Param::MIP::Cuts::ZeroHalfCut, -1);

      bool solved = cplex.solve();

      if (solved)
         env.out() << "Optimal tour length "
                   << cplex.getObjValue() << endl;

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
      delete mc;
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