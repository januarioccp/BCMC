#include "MaxBacker.h"

void MaxBacker::updateMaxBack(const vector<vector<double> > &wf, vector<int> &S1)
{
    this->maximumBack(wf,S1);
}

MaxBacker::MaxBacker(const vector<vector<double> > &wf, vector<int> &S1)
{
    this->maximumBack(wf,S1);
}

void MaxBacker::maximumBack(const vector<vector<double> > &w, vector<int> &Smin)
{
    /**
     * Number of vertices
     * */
    int n = w.size();

    /**
     * Create a vertex set in a vector of pair
     * */
    vector<pair<double, int> > V;
    for (unsigned i = 0; i < n; i++)
        V.push_back(make_pair(0.0, i));

    /**
     * Initialization
     * */
    Smin.clear();
    vector<int> Set;
    Set.push_back(V.front().second);
    V.erase(V.begin());
    Smin = Set;
    Cutmin = 0.0;
    
    for(int i = 0; i < V.size(); i++)
        Cutmin += w[Set.back()][V[i].second];
    double Cutval = Cutmin;

    /**
     * Update weights based on the last inserted vertex
     * */
    for(int i = 0; i < V.size(); i++)
        V[i].first += w[Set.front()][V[i].second];

    while (Set.size() < n - 1)
    {
        /** 
         * maximum max-back index
         * */
        int mmbi = max_element(V.begin(), V.end()) - V.begin();

        /**
         * Add to S the vertex v that is not in S and v has the maximum max-back value - mmbv
         * */
        Set.push_back(V[mmbi].second);

        /**
         * Compute Cutval
         * */
        Cutval = Cutval + 2.0 - 2.0 * V[mmbi].first;

        /**
         * Update Cutmin, if necessary
         * */
        if (Cutval < Cutmin)
        {
            Cutmin = Cutval;
            Smin = Set;
        }

        /**
        * Update weights based on the last inserted
        * */
        for(int i = 0; i < V.size(); i++)
            V[i].first += w[V[mmbi].second][V[i].second];

        /**
         * Remove mmbv from V
         * */
        V.erase(V.begin() + mmbi);
    }
    // cout<<Cutmin<<endl;
}