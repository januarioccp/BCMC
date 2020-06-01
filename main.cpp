/**
 *  Solving a traveling salesman problem 
 * using lazy constraints and max-back cuts.   
*/

// You might need to change it in order to make it work in your machine
#include "gurobi_c++.h"
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <algorithm>
#include "./src/data.h"
using namespace std;

string itos(int i)
{
  stringstream s;
  s << i;
  return s.str();
}

void findsubtour(int n, double **sol, int *tourlenP, int *tour);
void findsubtour2(int n, double **sol, vector<pair<int, int> > &p, int *tour);
vector<int> maxback(int n, double **sol, double& Cutmin, int v0);

// Subtour elimination callback.  Whenever a feasible solution is found,
// find the smallest subtour, and add a subtour elimination constraint
// if the tour doesn't visit every node.

class subtourelim : public GRBCallback
{
public:
  GRBVar **vars;
  int n;
  subtourelim(GRBVar **xvars, int xn)
  {
    vars = xvars;
    n = xn;
  }

protected:
  void callback()
  {
    try
    {
       if (where == GRB_CB_MIPNODE)
      //  (false)
      {
        if(GRB_CB_MIPNODE_STATUS == GRB_OPTIMAL)
          cout<<"Relaxation if Optimal"<<endl;

        double **x = new double *[n];

        for (int i = 0; i < n; i++)
          x[i] = getSolution(vars[i], n);

        double Cutmin;
        srand(time(0));
        int v0 = rand()%n;
        vector<int> v = maxback(n, x, Cutmin, v0);
        if (Cutmin < 1.9)
        {
          vector<bool> seen(n, false);
          for (int i = 0; i < v.size(); i++)
            seen[v[i]] = true;

          vector<int> y;
          for (int i = 0; i < n; i++)
            if (!seen[i])
              y.push_back(i);

          if(v.size() > 0 && y.size() > 0){
            GRBLinExpr expr = 0;

            for (int i = 0; i < v.size(); i++)
              for (int j = 0; j < y.size(); j++)
                  expr += vars[v[i]][y[j]] + vars[y[j]][v[i]];
                // cout<<expr<<endl;
            addCut(expr,GRB_GREATER_EQUAL,2);
            //  cout<<"cut"<<endl;
          }
        }
          for (int i = 0; i < n; i++)
            delete[] x[i];
          delete[] x;
        
      }
      if (where == GRB_CB_MIPSOL)
      {
        // Found an integer feasible solution - does it visit every node?
        double **x = new double *[n];
        // int *tour2 = new int[n];
        // vector<pair<int, int> > p;
        
        for (int i = 0; i < n; i++)
          x[i] = getSolution(vars[i], n);

        double Cutmin;
        srand(time(0));
        int v0 = rand()%n;
        vector<int> v = maxback(n, x, Cutmin, v0);
        if (Cutmin < 2)
        {
          vector<bool> seen(n, false);
          for (int i = 0; i < v.size(); i++)
            seen[v[i]] = true;

          vector<int> y;
          for (int i = 0; i < n; i++)
            if (!seen[i])
              y.push_back(i);

          if(v.size() > 0 && y.size() > 0){
            GRBLinExpr expr = 0;

            for (int i = 0; i < v.size(); i++)
              for (int j = 0; j < y.size(); j++)
                  expr += vars[v[i]][y[j]];// + vars[y[j]][v[i]];
                // cout<<expr<<endl;
            addLazy(expr,GRB_GREATER_EQUAL,2);
            //  cout<<"cut"<<endl;
          }
        }

        // findsubtour2(n, x, p, tour2);

        // if (p.size() > 1)
        // {
        //   for (i = 0; i < p.size(); i++)
        //   {

        //     GRBLinExpr expr1 = 0;
        //     for (int a = p[i].first; a < p[i].first + p[i].second; a++)
        //       for (int b = a + 1; b < p[i].first + p[i].second; b++)
        //         expr1 += vars[tour2[a]][tour2[b]];
        //     addLazy(expr1 <= p[i].second - 1);

        //     GRBLinExpr expr2 = 0;
        //     for (j = 0; j < p.size(); j++)
        //       if (i != j)
        //         for (int a = p[i].first; a < p[i].first + p[i].second; a++)
        //           for (int b = p[j].first; b < p[j].first + p[j].second; b++)
        //             expr2 += vars[tour2[a]][tour2[b]];

        //     addLazy(expr2 >= 2);
        //   }
        // }

        // delete[] tour2;
        for (int i = 0; i < n; i++)
          delete[] x[i];
        delete[] x;
        
      }
    }
    catch (GRBException e)
    {
      cout << "Error number: " << e.getErrorCode() << endl;
      cout << e.getMessage() << endl;
    }
    catch (...)
    {
      cout << "Error during callback" << endl;
    }
  }
};

// Given an integer-feasible solution 'sol', find the smallest
// sub-tour.  Result is returned in 'tour', and length is
// returned in 'tourlenP'.
void findsubtour(int n,
                 double **sol,
                 int *tourlenP,
                 int *tour)
{
  bool *seen = new bool[n];
  int bestind, bestlen;
  int i, node, len, start;

  for (i = 0; i < n; i++)
    seen[i] = false;

  start = 0;
  bestlen = n + 1;
  bestind = -1;
  node = 0;
  while (start < n)
  {
    for (node = 0; node < n; node++)
      if (!seen[node])
        break;

    if (node == n)
      break;

    for (len = 0; len < n; len++)
    {
      tour[start + len] = node;
      seen[node] = true;
      for (i = 0; i < n; i++)
      {
        if (sol[node][i] > 0.5 && !seen[i])
        {
          node = i;
          break;
        }
      }
      if (i == n)
      {
        len++;
        if (len < bestlen)
        {

          bestlen = len;
          bestind = start;
        }
        start += len;
        break;
      }
    }
  }

  for (i = 0; i < bestlen; i++)
    tour[i] = tour[bestind + i];
  *tourlenP = bestlen;

  delete[] seen;
}

// Given an integer-feasible solution 'sol', find the smallest
// sub-tour.
void findsubtour2(int n,
                  double **sol,
                  vector<pair<int, int> > &p,
                  int *tour)
{
  bool *seen = new bool[n];
  int i, node, len, start;

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

      //Be honest, you saw that node
      seen[node] = true;

      // Check nodes neighbors
      for (i = 0; i < n; i++)
      {
        // Is it connected to someone?
        // First time you see it? ...
        if (sol[node][i] > 0.5 && !seen[i])
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

  delete[] seen;
}

// Given an integer-feasible solution 'sol', find the smallest
// sub-tour.
vector<int> maxback(int n,
                    double **x, double &Cutmin, int v0)
{
  vector<int> S;
  vector<int> Smin;
  vector<bool> seen(n);
  vector<double> b(n);
  fill(seen.begin(), seen.end(), false);

  int v = v0;
  S.push_back(v);
  seen[v] = true;
  b[v] = -numeric_limits<double>::infinity();

  double Cutval;
  Cutmin = 0.0;
  for (int i = 0; i < n; i++)
  {
    if (!seen[i] && v!=i){
      b[i] = x[v][i];
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
        b[t] = b[t] + x[v][t];
    }
    if (Cutval < Cutmin)
    {
      Cutmin = Cutval;
      Smin = S;
    }
  }

  return Smin;
}

int main(int arg1,
         char *arg2[])
{
  if (arg1 < 2)
  {
    cout << "Usage: tsp_c++ instance" << endl;
    return 1;
  }

  Data input(arg1, arg2[1]);
  input.readData();

  int n = input.getDimension();
  double *x = new double[n];
  double *y = new double[n];

  int i;
  for (i = 0; i < n; i++)
  {
    x[i] = i;
    y[i] = i;
  }

  //GRBEnv *env = NULL;
   GRBEnv *env = NULL;
  GRBVar **vars = NULL;

  vars = new GRBVar *[n];
  for (i = 0; i < n; i++)
    vars[i] = new GRBVar[n];

  try
  {
    int j;

    env = new GRBEnv();
    GRBModel model = GRBModel(*env);

    // Algorithm used to solve continuous models

    // Must set LazyConstraints parameter when using lazy constraints
    model.set(GRB_IntParam_LazyConstraints, 1);

    // Create binary decision variables
    for (i = 0; i < n; i++)
      for (j = i; j < n; j++)
      {
        // vars[i][j] = model.addVar(0.0, 1.0, input.getDistance(i,j),GRB_BINARY, "x_" + itos(i) + "_" + itos(j));
        vars[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + itos(i) + "_" + itos(j));
        vars[j][i] = vars[i][j];
      }

    // Degree-2 constraints
    for (i = 0; i < n; i++)
    {
      GRBLinExpr expr = 0;
      for (j = 0; j < n; j++)
        expr += vars[i][j];
      model.addConstr(expr == 2, "deg2_" + itos(i));
    }

    // Set objective function
    GRBLinExpr obj = 0;
    for (i = 0; i < n; i++)
      for (j = i + 1; j < n; j++)
        obj += vars[i][j] * input.getDistance(i, j);
    model.setObjective(obj);

    // Forbid edge from node back to itself
    for (i = 0; i < n; i++)
      vars[i][i].set(GRB_DoubleAttr_UB, 0);

    // Set callback function
    subtourelim cb = subtourelim(vars, n);
    model.setCallback(&cb);

    // Time limit
    // model.set(GRB_DoubleParam_TimeLimit, 100);

    // Remove extra messages on the screen
    // model.set(GRB_IntParam_OutputFlag, false);

    // Disable all heuristics
    // model.set(GRB_DoubleParam_Heuristics, false);

    // Disable all kinds of cuts, including gomory
    // model.set(GRB_IntParam_Cuts, 0);

    // Disable all kinds of presolve
    // model.set(GRB_IntParam_Presolve, false);

    // When adding your own cuts, you must set parameter PreCrush to value 1. 
    // This setting shuts off a few presolve reductions that sometimes prevent 
    // cuts on the original model from being applied to the presolved model.
    model.set(GRB_IntParam_PreCrush, 1);

    // Set the model to minimization
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // Optimize model
    model.optimize();

    // Extract solution

    if (model.get(GRB_IntAttr_SolCount) > 0)
    {
      double **sol = new double *[n];
      for (i = 0; i < n; i++)
        sol[i] = model.get(GRB_DoubleAttr_X, vars[i], n);

      int *tour = new int[n];
      int len;

      findsubtour(n, sol, &len, tour);
      // assert(len == n);

      double distancia = 0;
      for (i = 0; i < len - 1; i++)
      {
        distancia += input.getDistance(tour[i], tour[i + 1]);
      }
      distancia += input.getDistance(tour[len - 1], tour[0]);
      printf("%.*s,", int(strlen(arg2[1])) - 14, arg2[1] + 10);
      cout << distancia << ",";
      cout << model.get(GRB_DoubleAttr_Runtime) << endl;

      for (i = 0; i < n; i++)
        delete[] sol[i];
      delete[] sol;
      delete[] tour;
    }
  }
  catch (GRBException e)
  {
    cout << "Error number: " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
  }
  catch (...)
  {
    cout << "Error during optimization" << endl;
  }

  for (i = 0; i < n; i++)
    delete[] vars[i];
  delete[] vars;
  delete[] x;
  delete[] y;
  delete env;
  return 0;
}