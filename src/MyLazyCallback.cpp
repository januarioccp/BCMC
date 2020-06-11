#include "MyLazyCallback.h"
// #include "NodeInfo.h"

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
      for (IloInt j = 0; j < n; j++){
         if(i < j)
            sol[i][j] = abs(getValue(x[i][j]));
         else
            sol[i][j] = 0;
      }
   }

   for (IloInt i = 0; i < n; i++)
      for (IloInt j = i+1; j < n; j++)
         sol[j][i] = sol[i][j];

    // Declares a boolean vector of size n with false
   vector<bool> seen(n, false);

   // An array to store the subtours
   IloNumArray tour(env, n);

   IloInt i, node, len, start;

   // Vector os positions
   vector<pair<int, int>> p;

   for (i = 0; i < n; i++)
      seen[i] = false;

   start = 0;
   node = 0;

   // Start from position 0 in the tour
   while (start < n)
   {
      // Find a node that has not been seen
      for (node = 0; node < n; node++)
         if (!seen[node])
            break;
      // Did you see every node? Time to break
      if (node == n)
         break;

      // Start with lenght 0 and build a subtour
      for (len = 0; len < n; len++)
      {
         // insert the unseen node in the tour
         tour[start + len] = node;

         //Be honest if you saw that node
         seen[node] = true;

         // Check nodes neighbors
         for (i = 0; i < n; i++)
         {
            // Is it connected to someone?
            // First time you see it? ...
            if (sol[node][i] >= 1.0 - EPSILON && !seen[i])
            {
               // ... better catch that guy
               node = i;
               break;
            }
         }

         // Oh man, you could not find a neighbhor? It seens that you closed the loop
         if (i == n)
         {
            // In this case, increase the size of the lenght
            len++;
            pair<int, int> pos(start, len);
            p.push_back(pos);
            //Start a new subtour
            start += len;
            break;
         }
      }
   }

   // Create and add subtour constraint ---
   // No more than 'length-1' edges between members of the subtour
   if (p.size() > 1)
   {
      for (i = 0; i < p.size(); i++)
      {
         IloExpr expr1(env);
         for (int a = p[i].first; a < p[i].first + p[i].second; a++)
            for (int b = a + 1; b < p[i].first + p[i].second; b++){
                if(tour[a] < tour[b])
                    expr1 += x[tour[a]][tour[b]];
                else
                    expr1 += x[tour[b]][tour[a]];
            }
         cout << "Adding lazy constraint " << expr1 << " <= " << p[i].second - 1 << endl;
         constraints->push_back(expr1 <= p[i].second - 1);
         expr1.end();
      }
   }


   seen.end();
   tour.end();
   for (IloInt i = 0; i < n; i++)
      sol[i].end();
   sol.end();

    return constraints;
}

std::mutex MyLazyCallback::lazyMutex;