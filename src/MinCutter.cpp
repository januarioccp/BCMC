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

MinCutter::MinCutter(const vector<vector<double>> &wf)
{
    this->w = wf;
    this->MINIMUMCUT();
}

MinCutter::MinCutter(const vector<vector<int>> &wf)
{
    this->w.resize(wf.size());
    for (int i = 0; i < wf.size(); i++)
    {
        this->w[i].resize(wf[i].size());
        for (int j = 0; j < wf[i].size(); j++)
            this->w[i][j] = double(wf[i][j]);
    }
    this->MINIMUMCUT();
}

MinCutter::~MinCutter()
{
    delete this->dSet;
}

void MinCutter::MINIMUMCUT()
{
    int n = this->w.size();
    
    // You need to find the value of the minimum cut
    this->minCut = numeric_limits<double>::max();

    // Create a vertex set
    G.resize(n);
    for (unsigned i = 0; i < G.size(); i++)
        G[i] = i;

    // Randomize the set of vertices - use -std=c++2a
    // obtain a time-based seed:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    shuffle(G.begin(), G.end(),default_random_engine(seed));

    // Use a disjoint set data structure to shrink G later
    dSet = new DisjSet(n);

    // Compute the minimumCut while there is
    // at least 2 vertices
    while (G.size() > 1)
    {
        // Shrink G
        G.erase(remove(G.begin(), G.end(), MINIMUMCUTPHASE()), G.end());
    }
}

int MinCutter::MINIMUMCUTPHASE()
{
    // A container to the vertices in this phase
    vector<int> A;

    // Store a copy of G but in a vector of pair
    vector< pair<double,int> > V;

    // Copying G to V
    for(auto i:G)
        V.push_back(make_pair(0,i));

    // Choose a vector from v=2 to insert in A
    A.push_back(V.front().second);

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

    // You need to do until A is as large as the initial size of G
    while (A.size() < n)
    {
        //Update weights
        for(auto &i:V)
            i.first+=w[A.back()][i.second];
        
        // Build heap
        make_heap(V.begin(),V.end());
        
        // Store current cut_of_the_phase value
        cut_of_the_phase = V.front().first;

        // Select the most tightly connected vertex - mtcv
        pair<double,int> mtcv = V.front();
        
        // Add to A the most tightly connected vertex - mtcv
        A.push_back(mtcv.second);

        // Remove mtcv from V
        V.erase(remove(V.begin(), V.end(), mtcv), V.end());
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

    for (auto i : A)
    {
        w[i][s] += w[i][last];
        w[s][i] = w[i][s];
    }

    return last;
}