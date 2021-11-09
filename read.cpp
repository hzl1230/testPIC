#include "utility.h"

using namespace std;

int main(int argc, char **argv)
{
    vector<Real> vec{1,2,3},info;
    Real vel, ndens;
    vel = 0.5, ndens = 1e5;
    // info.resize(vec.size());
    info.resize(vec.size());
    transform(std::move(vec).begin(),std::move(vec).end(),info.begin(),[&](auto& x){return x*vel*ndens;});
    cout << vec[0] << endl;
    print(info.begin(), info.end());
    // cout << a << endl;
}