#include "particle.h"
#include "utility.h"
#include "cross_section.h"
#include "espic_math.h"
#include "collision.h"
#include "work.h"
#include <iostream>

using std::cout;
using std::endl;
using std::vector;
Bigint ESPIC::Random::seed = (Bigint)time(NULL);
ESPIC::Random ranf;
// const Real ndens_ = 2414323.50534664;
const Real pmass = 1.;

int main(int argc, char** argv)
{
    const size_t npart = 500000;
    vector<Real> x(npart), vx(npart);
    vector<Real> y(npart), vy(npart);
    vector<Real> z(npart), vz(npart);
    Particles particles = Particles(x,y,z,vx,vy,vz);
    CrossSection* cross_section = new CrossSection("csection.in");
    vector<Reaction*>& reaction_arr = cross_section->react_arr; 
    Reaction*& reaction = reaction_arr[0];
    Real nu_max;
    InitParticle(particles, reaction, nu_max);
    
    int coll;
    CollPair newpp;

    
    coll = NullCollisionMethod(particles, nu_max, reaction, pmass, newpp);
    if(!newpp.empty()){
     for(size_t ipair=0; ipair<newpp.size(); ++ipair) {
        particles.append(newpp[ipair].first);
        // particles->append(newpp[ipair].second);0
     } 
    }
    
    
    // if(!newpp.empty()){
    //     for(size_t ipair=0; ipair<newpp.size(); ++ipair) {
    //         particles.append(newpp[ipair].first);
    //         // particles->append(newpp[ipair].second);0
    //     } 
    // }
    
    delete cross_section;
    return 0;
}

