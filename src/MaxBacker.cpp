#include "MaxBacker.h"

ostream &operator<<(ostream &os, MaxBacker &m)
{
    pair<vector<int>, vector<int>> S = m.getPartition();
    os << m.Cutmin << endl;
    os << "{";
    for (auto i : S.first)
    {
        os << i + 1;
        if (i != S.first.back())
            os << ",";
    }
    os << "}{";
    for (auto i : S.second)
    {
        os << i + 1;
        if (i != S.second.back())
            os << ",";
    }
    os << "}";
    return os;
}

void MaxBacker::updateMaxBack(const vector<vector<double>> &wf)
{
    this->w = wf;
    this->maximumBack();
}

MaxBacker::MaxBacker(const vector<vector<double>> &wf)
{
    this->w = wf;
    this->maximumBack();
}

void MaxBacker::maximumBack()
{
    // // Number of vertices
//     int a = 1;
//     int N = w.size();
//     std::vector<bool> S1(N);
//     std::vector<bool> S = S1;
//     S[a] = true;
//     int s0_size = 1;
//     std::vector<double> b(N);

//     double cutval = 0;
//     double cutmin = std::numeric_limits<double>::infinity();

//     for (int i = 0; i < N; i++)
//     {
//         if (i == a)
//         {
//             b[i] = -std::numeric_limits<double>::infinity();
//             continue;
//         }
//         b[i] = w[a][i];
//         cutval += b[i];
//     }

//     while (s0_size < S.size() - 1)
//     {
//         int maxB = std::max_element(b.begin(), b.end()) - b.begin();
//         s0_size++;

//         cutval += 2 - 2 * b[maxB];

//         // std::cout << cutval << "\n";

//         S[maxB] = true;
//         b[maxB] = -std::numeric_limits<double>::infinity();

//         for (int i = 0; i < N; i++)
//         {
//             if (S[i] == true)
//             {
//                 continue;
//             }
//             b[i] += w[maxB][i];
//         }

//         if (cutval - cutmin < 0)
//         {
//             cutmin = cutval;
//             S1 = S;
//         }
//     }

//     cout<<"Cutmin "<<cutmin<<endl;


// //////////////////////////////////////////


    // Number of vertices
    int n = w.size();

    // Create a vertex set in a vector of pair
    vector<pair<double, int> > V;
    for (unsigned i = 0; i < n; i++)
        V.push_back(make_pair(0, i));

    S.clear();
    S.push_back(V.front().second);
    V.erase(V.begin());
    Smin = S;

    // You need to find the value of the maxback
    Cutmin = 0.0;
    for (auto i : V)
        Cutmin += w[S.back()][i.second];
    double Cutval = Cutmin;

    while (S.size() < n)
    {
        //Update weights based on the last inserted
        for (auto &i : V)
            i.first += w[S.back()][i.second];

        // Build heap
        make_heap(V.begin(), V.end());

        // Choose v not in S of maximum max-back value - mmbv
        pair<double, int> mmbv = V.front();

        // Add to S v not in S of maximum max-back value - mmbv
        S.push_back(mmbv.second);

        // Remove mmbv from V
        V.erase(V.begin());

        Cutval = Cutval + 2.0 - 2.0 * mmbv.first;

        if (Cutval < Cutmin)
        {
            Cutmin = Cutval;
            Smin = S;
        }
    }
    cout<<"Cutmin"<<Cutmin<<endl;
}

pair<vector<int>, vector<int> > MaxBacker::getPartition(){
    vector<bool> seen(w.size(),false);
    S.clear();
    for(auto i: Smin)
        seen[i] = true;
    for(int i = 0; i < seen.size(); i++)
        if(!seen[i])
            S.push_back(i);
        
    return make_pair(Smin,S);
}