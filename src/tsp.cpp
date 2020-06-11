// TSP classes and global objects
#include "input.h"
#include "solution.h"
#include "construction.h"
#include "neighborhood.h"
#include "perturbation.h"
#include "localsearch.h"
#include "data.h"
#include "MinCutter.h"
MinCutter *mc;

// CPLEX classes and global objects
#include <ilcplex/ilocplex.h>
#include "NodeInfo.h"
#include "MyLazyCallback.h"
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
		swap(x,y);
	if( x == 0)
		return y-x;
    else
		return (n-1)*x-(x*(x+1)/2)+y;
}

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
      for (IloInt j = 0; j < n; j++){
         if(i < j)
            sol[i][j] = abs(getValue(x[i][j]));
         else
            sol[i][j] = 0;
      }
   }

   for (IloInt i = 0; i < n; i++)
      for (IloInt j = i+1; j < n; j++)
         sol[j][i] = sol[i][j];


   // Declares a boolean vector of size n with false
   vector<bool> seen(n, false);

   // An array to store the subtours
   vector<int> tour(n);

   IloInt i, node, len, start;

   // Vector os positions
   vector<pair<int, int>> position;

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
            position.push_back(pos);
            //Start a new subtour
            start += len;
            break;
         }
      }
   }

   // for(auto city: tour)
   //    cout<<setw(3)<<city;
   // cout<<endl;

   // for(auto p: position)
   //    cout<<setw(3)<<p.first<<setw(3)<<p.second<<endl;
   // cout<<endl;

   // Create and add subtour constraint ---
   // No more than 'length-1' edges between members of the subtour
   if (position.size() > 1)
   {
      for (auto p: position)
      {
         IloExpr expr1(env);
         int i,s;
         for (i = p.first,s=0; s < p.second-1; i++,s++){
            // cout<<"("<<
            //    tour[i]<<","<<
            //    tour[i+1]<<") ="<<
            //    toEdge(tour[i],tour[i+1],n)<<
            //       endl;
               if(tour[i] < tour[i+1])
                    expr1 += x[tour[i]][tour[i+1]];
                else
                    expr1 += x[tour[i+1]][tour[i]];
         }
         // cout<<"("<<
         //       tour[i]<<","<<
         //       tour[p.first]<<") ="<<
         //       toEdge(tour[i],tour[p.first],n)<<
         //          endl;
         if(tour[i] < tour[p.first])
            expr1 += x[tour[i]][tour[p.first]];
         else
            expr1 += x[tour[p.first]][tour[i]];
         // cout << "Adding lazy constraint " << expr1 << " <= " << p.second - 1 << endl;
         add(expr1 <= p.second - 1).end();
         expr1.end();
      }
   }

   // for (IloInt i = 0; i < n; i++)
   //    for (IloInt j = 0; j < i; j++)
   //       sol[j][i] = sol[i][j];

   // vector<int> S;
   // vector<int> Smin;
   // vector<double> b(n);
   // fill(seen.begin(), seen.end(), false);

   // int v = rand() % n;
   // S.push_back(v);
   // seen[v] = true;
   // b[v] = -numeric_limits<double>::infinity();

   // double Cutmin;
   // double Cutval;
   // Cutmin = 0.0;
   // for (int i = 0; i < n; i++)
   // {
   //    if (!seen[i] && v != i)
   //    {
   //       b[i] = sol[v][i];
   //       Cutmin += b[i];
   //    }
   // }

   // Smin = S;
   // Cutval = Cutmin;

   // while (S.size() < n)
   // {
   //    // Choose v not in S of maximum max-back value b(v)
   //    double maxb = -numeric_limits<double>::infinity();
   //    for (int i = 0; i < n; i++)
   //       if (!seen[i])
   //          if (maxb < b[i]){
   //             maxb = b[i];
   //             v = i;
   //          }

   //    S.push_back(v);
   //    seen[v] = true;
   //    Cutval = Cutval + 2 - 2 * b[v];

   //    for (int t = 0; t < n; t++)
   //    {
   //       if (!seen[t] && v != t)
   //          b[t] = b[t] + sol[v][t];
   //    }
   //    if (Cutval < Cutmin)
   //    {
   //       Cutmin = Cutval;
   //       Smin = S;
   //    }
   // }

   // if (Cutmin < 2 - tol)
   // {
   //    fill(seen.begin(), seen.end(), false);
   //    vector<int> v = Smin;
   //    for (int i = 0; i < v.size(); i++)
   //       seen[v[i]] = true;
   //    vector<int> y;
   //    for (int i = 0; i < n; i++)
   //       if (!seen[i])
   //          y.push_back(i);
   //    if (v.size() > 0 && y.size() > 0)
   //    {
   //       IloExpr expr1(env);
   //       for (int i = 0; i < v.size(); i++)
   //          for (int j = 0; j < y.size(); j++)
   //             expr1 += x[v[i]][y[j]] + x[y[j]][v[i]];
   //       // cout << "Adding cut " << expr1 << " >= 2 " << endl;
   //       add(expr1 >= 2).end();
   //       expr1.end();
   //    }
   // }
   seen.end();
   // tour.end();
   for (IloInt i = 0; i < n; i++)
      sol[i].end();
   sol.end();
}

int main(int argc, char **argv)
{
   // Run RVND to identify an UB
   Input in(argc, argv);   
   Solution sol(&in);
	LocalSearch ls(&in);
   sol = ls.GILSRVND();
   int UB = sol.costValueTSP + 1;
   cout<<"Upper bound |- - - - - - - - - - - - - - -| "<<UB<<endl;
   // Read input to get dimmention and initialize MinCutter
   Data input(argc, argv[1]);
   input.readData();
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

      // Declaring initial solution
      IloNumArray x_start(env,(n*(n-1)/2));
      for(int i = 0; i < x_start.getSize(); i++)
         x_start[i]=0;

      // Build initial solution
      for (int i = 0; i < sol.location.size() - 1; i++)
		   x_start[toEdge(sol.location[i],sol.location[i+1],n)] = 1;
      x_start[toEdge(sol.location.back(),sol.location.front(),n)] = 1;

      IloNumVarArray y(env);
	   for(int i = 0; i < n;i++)
         for(int j = i+1; j < n;j++)
		      y.add(x[i][j]);

      cplex.addMIPStart(y, x_start);

      IloNum tol = cplex.getParam(IloCplex::EpInt);

      // Subtour Elimination Callback
      MyLazyCallback *lazyCbk = new (env) MyLazyCallback(env,x);
      cplex.use(lazyCbk);

      // cplex.use(SubtourEliminationCallback(env, x, tol));
      // cplex.use(MaxBack(env, x, tol));
      //      cplex.use(MinCut(env, x, tol));
      cplex.setParam(IloCplex::PreInd, IloFalse);
      cplex.setParam(IloCplex::TiLim, 2 * 60 * 60);
	   cplex.setParam(IloCplex::Threads, 1);
	   cplex.setParam(IloCplex::CutUp, UB + 1);
      // cplex.setParam(IloCplex::Param::MIP::Cuts::FlowCovers, -1);
      // cplex.setParam(IloCplex::Param::MIP::Cuts::Gomory, -1);
      // cplex.setParam(IloCplex::Param::MIP::Cuts::ZeroHalfCut, -1);

      // Export the LP model to a txt file to check correctness
      cplex.exportModel("model.lp");
      bool solved = cplex.solve();

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