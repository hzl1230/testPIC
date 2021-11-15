#include "cross_section.h"
#include "particle.h"
using namespace std;

int main(int argc, char **argv)
{
    CrossSection* cross_section = new CrossSection("csection.in");
    Reaction* & reaction = cross_section->react_arr[0];
    Real en = 5.8348;
    Real nvt = sqrt(2.*en) * 2414323.5;
    vector<Real> info;
    vector<Particle> particles;
    const Real vth = sqrt(2.*0.01/72820.7);
    auto fnu = [=](auto& x){ return x*nvt; };
    vector<Real> nu;
    info = reaction->en_cs(en);
    transform(info.begin(), info.end(), back_inserter(nu), fnu);

    for(auto in : info)
        cout << in << " ";
    cout << endl;
    for(auto in : nu)
        cout << in << " ";
    cout << endl;

    delete cross_section;
    return 0;
}