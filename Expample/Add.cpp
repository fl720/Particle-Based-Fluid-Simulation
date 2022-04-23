#include <iostream>
#include <vector>

using namespace std;

class LARGE_OBJ
{
    int a;
    double arr[1000];
};

void add_adjacent(vector<LARGE_OBJ> &vec)
{
    if(vec.size() <= 1) return ;
    
    vector<int> tem = vec ;

    for( int i = 1 ; i < vec.size() -1 ; i++ )
    {
        vec[i] = tem[i-1] + tem[i] + tem[i+1] ; 
    }

    vec[0] = tem[0] + tem[1] ;
    vec[vec.size()-1] =  tem[vec.size()-1] + tem[vec.size()-2];
    
    return;
}

int main()
{
    vector<LARGE_OBJ> v = {1,1,1,1,1,1};

    add_adjacent(v); // {2,3,3,3,3,2}

    for(auto i : v)
        cout << v << ' ';

    return 0;
}