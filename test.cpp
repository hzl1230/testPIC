#include "particle.h"
#include "utility.h"
#include "cross_section.h"
#include "espic_math.h"
#include <iostream>

using std::cout;
using std::endl;
using std::vector;
Bigint ESPIC::Random::seed = (Bigint)time(NULL);

const Real ndens_ = 2414323.50534664;
const Real dt = 0.1;
const Real vth = sqrt(2.);
const Real vthb = sqrt(2.*0.01/72820.7);

int main(int argc, char** argv)
{
    const size_t n = 100000;
    vector<Real> x(n), vx(n);
    vector<Real> y(n), vy(n);
    vector<Real> z(n), vz(n);

    Particles* ps = new Particles(x,y,z,vx,vy,vz);
    CrossSection* arcs = new CrossSection("csection.in");
    Reaction *er(arcs->react_arr[0]);
    // std::ofstream of("vel.dat");
    // std::ofstream ofe("energy.dat");
    
    vector<Real> cs_info, info0, info1;
    Real nu_max = 0, energy1;
    
    for(size_t ip = 0; ip < (*ps).size(); ++ip) {        
        Particle& pt = (*ps)[ip];
        Real vxb, vyb, vzb, energy;
        VelBoltzDistr(vth, pt.vx(), pt.vy(), pt.vz());
        VelBoltzDistr(vthb, vxb, vyb, vzb);
        RelativeVelocity(pt, vxb, vyb, vzb);
        energy = get_energy(pt.vx(), pt.vy(), pt.vz());
        energy1 = get_energy(pt.vxr(), pt.vyr(), pt.vzr());
        if(!cs_info.empty()) cs_info.clear();
        cs_info = er->en_cs(energy1);
    
        Real cstot = vsum(cs_info);
        Real vt = sqrt(2.*energy);

        Real nutot(cstot*vt*ndens_);
        vector<Real>& nu = pt.nu();
        if (!nu.empty()) nu.clear();
        nu.resize(cs_info.size());
        transform(cs_info.begin(), cs_info.end(), 
                  nu.begin(), [&](auto& x){ return x*vt*ndens_; });
        if ( nutot > nu_max ) nu_max = nutot;
    }

    size_t npcoll = static_cast<size_t>(Pcoll(nu_max, dt) * n);
    cout << nu_max << endl;
    cout << Pcoll(nu_max, dt) << " " << npcoll << endl;
    Particles subps = Particles();
    ps->get_sub_particles(npcoll, subps);

    size_t coll = 0;
    ESPIC::Random ranf;
    for(size_t ip = 0; ip < npcoll; ++ip) {
        Particle& particle = subps[ip];
        Real nutot = vsum(particle.nu());
        Real R = ranf();
        if (R < nutot/nu_max) 
            coll+=1;
    }
    cout << coll << endl;
    Real time_aver = 1./((Real)n * dt);
    Real coll_freq = coll*time_aver;
    cout << coll_freq << endl;

    // of.close();
    // ofe.close();
    delete arcs;
    // delete subps;
    return 0;
}

