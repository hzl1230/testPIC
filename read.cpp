#include "utility.h"

using namespace std;

int main(int argc, char **argv)
{
    vector<Real> vec{1,2,3},info;
    // info.resize(vec.size());
    info.reserve(vec.size());
    transform(vec.begin(),vec.end(),info.begin(),[&](auto& x){return x*2;});
    print(info.begin(), info.end());
    // cout << a << endl;
}