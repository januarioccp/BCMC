#include "MyCutCallback.h"

#define MAX_ITER 100
#define MB_DEPTH 10

MyCutCallback::MyCutCallback(IloEnv env, const IloArray<IloBoolVarArray> &par_x) : IloCplex::UserCutCallbackI(env), x(env)
{
    x = par_x;
    IloInt n = x.getSize();
    myMinCut = nullptr;
}

MyCutCallback::~MyCutCallback(){
    delete myMinCut;
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

    if (data->getIterations() >= MAX_ITER || data->depth > MB_DEPTH)
    {
        lazyMutex.unlock();
        return;
    }

    IloEnv env = getEnv();
    IloInt n = x.getSize(); // Number of vertices

    // Retrieve solution information
    vector<vector<double>> w = vector<vector<double>>(n,vector<double>(n,0));
    for (IloInt i = 0; i < n; i++){
        for (IloInt j = i+1; j < n; j++){
                w[i][j] = abs(getValue(x[i][j]));
                w[j][i] = w[i][j];
        }
    }       

    if(myMinCut == nullptr)
        myMinCut = new MinCutter(w);
    else
        myMinCut->updateMinCut(w);

    // Add the constraint
    if (myMinCut->getMinCut() < 2.0 - EPSILON)
    {
        pair<vector<int>, vector<int>> partition = myMinCut->getPartition();   
        IloExpr expr1(env);
        for (auto i : partition.first)
            for (auto j : partition.second){
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