#include "MinCutter.h"
#include "DisjSet.h"

ostream &operator<<(ostream &os, const MinCutter &m)
{
    os << m.minCut << endl;
    os << "{";
    for (auto i : m.S1)
    {
        os << i;
        if (i != m.S1.back())
            os << ",";
    }
    os << "}{";
    for (auto i : m.S2)
    {
        os << i;
        if (i != m.S2.back())
            os << ",";
    }
    os << "}";
    return os;
}

MinCutter::MinCutter(vector<vector<double>> &w)
{
    
    n = w.size();

    // You need to find the value of the minimum cut
    this->minCut = numeric_limits<double>::max();

    // Create a vertex set
    G.resize(n);

    this->MINIMUMCUT(w);
}

MinCutter::MinCutter(int size){
    this-> n = size;
    // You need to find the value of the minimum cut
    this->minCut = numeric_limits<double>::max();
    dSet = new DisjSet(n);
    // Create a vertex set
    G.resize(n);
}

void MinCutter::MINIMUMCUTUPDATE( vector<vector<double> > &w){
    // You need to find the value of the minimum cut
    this->minCut = numeric_limits<double>::max();
    dSet->makeSet();
    for (unsigned i = 0; i < G.size(); i++)
        G[i] = i;
    while (G.size() > 1)
        MINIMUMCUTPHASE(w);
}

void MinCutter::MINIMUMCUT(vector<vector<double> > &w)
{
    for (unsigned i = 0; i < G.size(); i++)
        G[i] = i;

    // Do you wand randomize the set of vertices? - use -std=c++2a
    // obtain a time-based seed:
    // unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // shuffle(G.begin(), G.end(),default_random_engine(seed));

    // Use a disjoint set data structure to shrink G later
    dSet = new DisjSet(n);

    // Compute the minimumCut while there is
    // at least 2 vertices -> shrink G
    while (G.size() > 1)
        MINIMUMCUTPHASE(w);
}

void MinCutter::MINIMUMCUTPHASE(vector<vector<double> > &w)
{
    // A container to the vertices in this phase
    vector<int> A;

    // Store a copy of G but in a vector of pair
    vector< pair<double,int> > V;

    // Copying G to V
    for(auto i:G)
        V.push_back(make_pair(0,i));

    // Choose a vector from v=2 to insert in A
    A.push_back(G.front());

    // Remove the vertex inserted in A from G
    V.erase(V.begin());

    // Initialize with the minimum value, because you
    // want to find the maximum cut value
    double cut_of_the_phase;

    // Use this auxiliary variable to help you to find
    // the largest value in each phase
    double cutWeight = 0.0;

    // Store the initial size of G
    int n = V.size() + A.size();

    auto mtcv = V.front();;
    // You need to do until A is as large as the initial size of G
    while (A.size() < n)
    {
        //Update weights
        for(auto &i:V)
            i.first+=w[A.back()][i.second];
        
        // Build heap
        make_heap(V.begin(),V.end());
        
        // Select the most tightly connected vertex - mtcv
        mtcv = V.front();

        // Store current cut_of_the_phase value
        cut_of_the_phase = mtcv.first;
        
        // Add to A the most tightly connected vertex - mtcv
        A.push_back(mtcv.second);

        // Remove mtcv from V
        V.erase(V.begin());
    }

    // Before last
    int s = *(A.end() - 2);
    // Last added
    int last = A.back();

    if (cut_of_the_phase < this->minCut)
    {
        this->minCut = cut_of_the_phase;
        // Clear the partition set every time???
        S1.clear();
        S2.clear();
        for (int i = 0; i < w.size(); i++)
            if (dSet->find(last) == dSet->find(i))
                S1.push_back(i);
            else
                S2.push_back(i);
    }

    // Merge the two last vertex added last
    dSet->Union(s, last);

    // Update weight of merged vertices
    for (auto i : A){
        w[i][s] += w[i][last];
        w[s][i] = w[i][s];
    }

    // Remove mtcv from G
    G.erase(G.begin());
}

MinCutter::~MinCutter()
{
    delete this->dSet;
}