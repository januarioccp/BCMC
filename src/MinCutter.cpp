#include "MinCutter.h"

void MinCutter::updateMinCut(vector<vector<double> > &w, vector<int> &S){
    this->minimumCut(w,S);
}

MinCutter::MinCutter(vector<vector<double> > &w, vector<int> &S){
    partition = nullptr;
    this->minimumCut(w,S);
}

MinCutter::~MinCutter()
{
    delete this->partition;
    delete this->best_partition;
}

void MinCutter::minimumCut(vector<vector<double> > &w, vector<int> &S){
    int n = w.size();
    // You need to find the value of the minimum cut
    this->minCut = numeric_limits<double>::max();

    // Create a vertex set
    G.resize(n);
    for (unsigned i = 0; i < G.size(); i++)
        G[i] = i;

    // Use a disjoint set data structure to shrink G later
    if(partition == nullptr){
        partition = new DisjSet(n);
        best_partition = new DisjSet(n);
    }
    else{
        partition->makeSet();
        best_partition->makeSet();
    }

    // Compute the minimumCut while there is
    // at least 2 vertices
    while (G.size() > 1)
    {
        // Shrink G
        G.erase(remove(G.begin(), G.end(), MINIMUMCUTPHASE(w)), G.end());
    }

    /**
     * Retrieving best partition
     * */
    for (int i = 0; i < w.size(); i++)
        if (best_partition->find(best_last) == best_partition->find(i))
            S.push_back(i);
}

int MinCutter::MINIMUMCUTPHASE(vector<vector<double> > &w){
    // A container to the vertices in this phase
    // A.clear();
    A.resize(0);
    

    // Store a copy of G but in a vector of pair
    vector< pair<double,int> > V;

    // Copying G to V
    for(int i:G)
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

    //Update weights based on the last inserted
        for(int i = 0; i < V.size(); i++)
            V[i].first+=w[A.front()][V[i].second];

    // You need to do until A is as large as the initial size of G
    while (A.size() < n)
    {
        // Build heap
        make_heap(V.begin(),V.end());
        
        // Store current cut_of_the_phase value
        cut_of_the_phase = V.front().first;

        // Select the most tightly connected vertex - mtcv
        pair<double,int> mtcv = V.front();
        
        // Add to A the most tightly connected vertex - mtcv
        A.push_back(mtcv.second);

        // Remove mtcv from V
        V.erase(V.begin());

        //Update weights based on the last inserted
        for(int i = 0; i < V.size(); i++)
            V[i].first+=w[mtcv.second][V[i].second];
    }

    // Before last
    int s = *(A.end() - 2);
    // Last added
    int last = A.back();

    if (cut_of_the_phase < this->minCut)
    {
        this->minCut = cut_of_the_phase;
        (*best_partition) = (*partition);
        best_last = last;
    }
    partition->Union(s, last);

    // Merge the two last vertex added last
    for (int i : A)
    {
        w[i][s] += w[i][last];
        w[s][i] = w[i][s];
    }

    return last;
}