// #include<iostream>
// #include<iomanip>
// using namespace std;

// int toEdge(int x, int y, int n){
// 	if(x == y)
// 		return 0;
// 	if(x>y)
// 		swap(x,y);
// 	if( x == 0)
// 		return y-x;
//     else
// 		return (n-1)*x-(x*(x+1)/2)+y;
// }

// int main(int argc, char **argv)
// {
//    int N = 6;
//    for(int i=0; i < N; i++){
//       for(int j=0; j < N; j++)
//          cout<<setw(3)<<toEdge(i,j,N);
//       cout<<endl;
//    }
//    return 0;
// }