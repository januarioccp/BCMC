#ifndef MYLAZYCALLBACK_H
#define MYLAZYCALLBACK_H

#include <ilcplex/ilocplex.h>
#include <vector>
#include <mutex>
#include <iostream>
using namespace std;

class MyLazyCallback : public IloCplex::LazyConstraintCallbackI
{
public:
  MyLazyCallback(IloEnv env, const IloArray<IloBoolVarArray> &x_ref);
  IloCplex::CallbackI *duplicateCallback() const;
  void main();

  std::vector<IloConstraint> *separate();
  static std::mutex lazyMutex;

private:
  IloArray<IloBoolVarArray> x;
};
#endif