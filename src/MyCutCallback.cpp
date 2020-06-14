#include "MyCutCallback.h"

#define MAX_ITER 100
#define MB_DEPTH 10

MyCutCallback::MyCutCallback(IloEnv env, const IloArray<IloBoolVarArray> &par_x) : IloCplex::UserCutCallbackI(env), x(env)
{
    x = par_x;
    IloInt n = x.getSize();
    myMinCut = nullptr;
    myMaxBack = nullptr;
}

MyCutCallback::~MyCutCallback()
{
    delete myMinCut;
    delete myMaxBack;
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

    if (!data)
    {
        if (NodeInfo::rootData == NULL)
            NodeInfo::initRootData();
        data = NodeInfo::rootData;
    }

    if (data->getIterations() >= MAX_ITER)
    {
        lazyMutex.unlock();
        return;
    }

    IloEnv env = getEnv();
    IloInt n = x.getSize(); // Number of vertices

    // Retrieve solution information
    vector<vector<double>> w = vector<vector<double>>(n, vector<double>(n, 0));
    for (IloInt i = 0; i < n; i++)
    {
        for (IloInt j = i + 1; j < n; j++)
        {
            w[i][j] = abs(getValue(x[i][j]));
            w[j][i] = w[i][j];
        }
    }

    if (myMaxBack == nullptr)
        myMaxBack = new MaxBacker(w);
     else
         myMaxBack->updateMaxBack(w);

    // // True if maxback heuristic is supposed to be used
    bool useMB = false;

    // // True if mincut algorithm is supposed to be used
    bool useMC = false;

    if (true && data->depth < MB_DEPTH)
    {
        if (myMinCut == nullptr)
            myMinCut = new MinCutter(w);
        else
            myMinCut->updateMinCut(w);

        if (myMinCut->getMinCut() < 2.0 - EPSILON)
            useMC = true;
    }
    else
        useMB = false;

    // Add the cut
    if (useMB || useMC)
    {
        pair<vector<int>, vector<int>> partition;
        if(useMC)
            partition = myMinCut->getPartition();
        else
            partition = myMaxBack->getPartition();
        
        IloExpr expr1(env);
        for (auto i : partition.first)
            for (auto j : partition.second)
            {
                if (i < j)
                    expr1 += x[i][j];
                else
                    expr1 += x[j][i];
            }
        add(expr1 >= 2).end();
        expr1.end();
    }

    data->addIteration();
    lazyMutex.unlock();
}

std::mutex MyCutCallback::lazyMutex;