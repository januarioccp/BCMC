#ifndef MYCUTCALLBACK_H
#define MYCUTCALLBACK_H

// User's library
#include "NodeInfo.h"
#include "MinCutter.h"

// CPLEX
#include <ilcplex/ilocplex.h>

// STL
#include <algorithm>
#include <exception>
#include <iostream>
#include <limits.h>
#include <list>
#include <mutex>
#include <stack>
#include <stdlib.h>
#include <vector>
using namespace std;

class MyCutCallback : public IloCplex::UserCutCallbackI
{
public:
  MyCutCallback(IloEnv env, const IloArray<IloBoolVarArray> &par_x);
  ~MyCutCallback();
  IloCplex::CallbackI *duplicateCallback() const;
  void main();

  static std::mutex lazyMutex;

private:
  IloArray<IloBoolVarArray> x;
  MinCutter* myMinCut;
};

#endif