#include "particle.h"
#include "utility.h"
#include "cross_section.h"
#include <iostream>

using std::cout;
using std::endl;
using std::vector;

const Real ndens_ = 2414323.50534664;
const Real dt = 0.1;

int main(int argc, char** argv)
{
    const size_t n = 100000;
    vector<Real> x, vx;
    vector<Real> y, vy;
    vector<Real> z(n,0), vz(n,0);
    random_sample_Maxwell(n, x, vx);
    random_sample_Maxwell(n, y, vy);
    Particles* ps = new Particles(x,y,z,vx,vy,vz);
    CrossSection* arcs = new CrossSection("csection.in");
    Reaction *er(arcs->react_arr[0]);
    vector<string> name{er->spec_name()};
    // int length(er->size()), size(er->isize());
    vector<Real> eth(er->th()), p_en(ps->get_particles_energy()), cs_info;
    Real nu_max = 0;
    
    for(size_t ip = 0; ip < p_en.size(); ++ip) {
        if(!cs_info.empty()) cs_info.clear();
        cs_info = er->en_cs(p_en[ip]);
        Real cstot = vsum(cs_info);
        Particle& particle = (*ps)[ip];
        Real vt = particle.vel_tot();
        Real nutot(cstot*vt*ndens_);
        vector<Real> nu(cs_info.size());
        transform(std::move(cs_info).begin(),std::move(cs_info).end(), 
                  nu.begin(), [&](auto& x){ return x*vt*ndens_; });
        particle.update_nu(nu);
        if ( nutot > nu_max ) nu_max = nutot;
    }
    size_t npcoll = static_cast<size_t>(Pcoll(nu_max, dt) * n);
    cout << Pcoll(nu_max, dt) << " " << npcoll << endl;
    Particles* subps = new Particles();
    ps->get_sub_particles(npcoll, subps);
    std::ofstream of("out.dat");
    size_t coll = 0;
    for(size_t ip = 0; ip < npcoll; ++ip) {
        Particle& particle = (*subps)[ip];
        Real nutot = vsum(particle.nu());
        Real R = random_G01();
        of << R << " " << nutot/nu_max << endl;
        if (R < nutot/nu_max) 
            coll+=1;
    }
    cout << coll << endl;
    Real time_aver = 1./((Real)n * dt);
    Real coll_freq = coll*time_aver;
    cout << coll_freq << endl;

    of.close();
    delete arcs;
    delete subps;
    return 0;
}

