#include "MyCutCallback.h"

#define MAX_ITER 100
#define MB_DEPTH 10

MyCutCallback::MyCutCallback(IloEnv env, const IloArray<IloBoolVarArray> &par_x) : IloCplex::UserCutCallbackI(env), x(env)
{
    x = par_x;
    IloInt n = x.getSize();
    
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

    if (data->getIterations() >= MAX_ITER || data->depth > MB_DEPTH)
    {
        lazyMutex.unlock();
        return;
    }

    IloEnv env = getEnv();
    IloInt n = x.getSize(); // Number of vertices

    // Retrieve solution information
    vector<vector<double>> w(n);
    for (IloInt i = 0; i < n; i++)
    {
        w[i].resize(n);
        for (IloInt j = 0; j < n; j++){
            if (i < j){
                w[i][j] = abs(getValue(x[i][j]));
            }
            else
                w[i][j] = 0;
        }
    }

    // Creating an adjacency matrix
    for (IloInt i = 0; i < n; i++)
        for (IloInt j = i+1; j < n; j++)
            w[j][i] = w[i][j];

    MinCutter m(w);
    cout<<m<<endl;

    // Add the constraint
    if (m.getMinCut() < 2.0)
    {
        // cout<<__FILE__<<" "<<__LINE__<<endl;
        pair<vector<int>, vector<int>> partition = m.getPartition();
        // double suma = 0;
        IloExpr expr1(env);
        for (auto i : partition.first)
            for (auto j : partition.second)
            {
                if (i < j)
                    expr1 += x[i][j];
                else
                    expr1 += x[j][i];
                // suma+=w[i][j];
                // cout<<setw(7)<<w[i][j];
            }
            // cout<<endl;
        // cout<<"Cut "<<expr1<<endl;
        //cout<<suma<<endl;
        add(expr1 >= 2).end();
        expr1.end();
    }

    data->addIteration();
    lazyMutex.unlock();
}

std::mutex MyCutCallback::lazyMutex;