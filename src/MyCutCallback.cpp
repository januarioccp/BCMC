#include "MyCutCallback.h"


#define MAX_ITER 100
#define MB_DEPTH 10

MyCutCallback::MyCutCallback(IloEnv env, const IloArray<IloBoolVarArray> &par_x) : IloCplex::UserCutCallbackI(env), x(env)
{
    x = par_x;
}

IloCplex::CallbackI *MyCutCallback::duplicateCallback() const
{
    return new (getEnv()) MyCutCallback(getEnv(), x);
}

void MyCutCallback::main()
{
    double EPSILON = 0.000001;
    lazyMutex.lock();

    NodeInfo *data = dynamic_cast<NodeInfo *>(getNodeData());

    if (!data){
        if (NodeInfo::rootData == NULL)
            NodeInfo::initRootData();
        data = NodeInfo::rootData;
    }

    if (data->getIterations() >= MAX_ITER)
    {
        lazyMutex.unlock();
        return;
    }

    // Insert Code Here!

    data->addIteration();
    lazyMutex.unlock();
}

std::mutex MyCutCallback::lazyMutex;