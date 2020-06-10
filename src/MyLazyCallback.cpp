#include "MyLazyCallback.h"
// #include "NodeInfo.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <limits.h>
#include <stdlib.h>
#include <stack>
#include <algorithm>
#include <exception>
#include <list>

MyLazyCallback::MyLazyCallback(IloEnv env, const IloArray<IloBoolVarArray> &par_x) : IloCplex::LazyConstraintCallbackI(env), x(env)
{
    x = par_x;
}

IloCplex::CallbackI *MyLazyCallback::duplicateCallback() const
{
    return new (getEnv()) MyLazyCallback(getEnv(), x);
}

void MyLazyCallback::main()
{
    // You shall not pass... actually, only one thread will pass
    lazyMutex.lock();

    // Use separate function to identify the constraints
    std::vector<IloConstraint> *cons = separate();
    for (int i = 0; i < cons->size(); i++){
        // Add each one of the constraints
        add((*cons)[i]);
        (*cons)[i].end();
    }
    delete cons;
    lazyMutex.unlock();
}

std::vector<IloConstraint> *MyLazyCallback::separate()
{
    double EPSILON = 0.000001;

    std::vector<IloConstraint> *constraints = new std::vector<IloConstraint>();

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

    cout<<"Calling"<<endl;

    return constraints;
}

std::mutex MyLazyCallback::lazyMutex;