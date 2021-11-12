#ifndef _WORK_
#define _WORK_
#include "cross_section.h"
#include "utility.h"
#include "particle.h"
#include "collision.h"

using std::cout;
using std::endl;
using std::ofstream;

const Real vth = sqrt(2.);
const Real vthb = sqrt(2.*0.01/72820.7);

const Real dt = 0.1;

void InitParticle(Particles& pts, Reaction*& reaction, Real& nu_max) { 
    nu_max = 0;

    for (size_t ip = 0; ip < pts.size(); ++ip) {
        Particle& pt = pts[ip];
        Real vxb, vyb, vzb, nvt;
        vector<Real> info(reaction->isize()), nu(reaction->isize());
        VelBoltzDistr(vth, pt.vx(), pt.vy(), pt.vz());
        VelBoltzDistr(vthb, vxb, vyb, vzb);
        RelativeVelocity(pt, vxb, vyb, vzb);
        Real en(pt.energy());
        info = reaction->en_cs(en);
        nvt = sqrt(2.0 * en) * ndens_;
        transform(info.begin(), info.end(), nu.begin(), [&](auto& x){return x*nvt;});
        Real nusum = 0;
        for(auto x : nu) { 
            (pt.nu()).emplace_back(x);
            nusum += x;
        }
        if (nusum > nu_max) nu_max = nusum;
    }
}

void ParticleCollision(
    const int type_id, 
    Reaction*& re,
    Collisionpair& cop,
    CollPair& newpp)
{
    Real threshold;
    std::string type = (re->get_types())[type_id];
    if (type_id == 0)  threshold = 0.0;
    else  threshold = (re->th())[type_id-1];

    if ("ela" == type)
        cop.ParticleElasticCollision();
    else if ("exc" == type)
        cop.ParticleExcitatinCollision(threshold);
    else if ("ion" == type) {
        cop.ParticleIonizationCollision(threshold);
        newpp = cop.ion_products();
    }
    else if ("iso" == type) 
        cop.ParticleIsotropicCollision();
    else if ("back" == type)
        cop.ParticleBackwardCollision();
    else {
        std::cout << type << std::endl;
        espic_error("Unknown Collision Type");
    }
}


int NullCollisionMethod(Particles& pts, Real nu_max, 
                         Reaction*& react, Real imass, CollPair& newpp)
{
    int ntype(react->isize());
    size_t ncoll = static_cast<Particles::size_type>(pts.size()*Pcoll(nu_max,dt));
    Particles subpts = Particles();
    pts.get_sub_particles(ncoll, subpts);
    int coll = 0, ela = 0, exc = 0, ion = 0;
    ofstream of("ncoll.dat", ofstream::app);

    for(size_t ip = 0; ip < subpts.size(); ++ip) {
        Particle& ptc = subpts[ip];
        std::vector<Real>& nu = ptc.nu();
        Real rnd = ranf(), nuj = 0;
        int itype = 0;
        while(itype != ntype) {
            nuj += nu[itype];
            if(rnd < nuj / nu_max) {
                Collisionpair cop = Collisionpair(ptc, imass, mass, vth);
                ParticleCollision(itype, imass, react, cop, newpp);
                if(itype == 0)  ela++;
                else if (itype == 1) exc++; 
                else ion++;
                coll++; break;
            }
            itype++;
        }
    }
    of << ela << " " << exc << " " << ion << " " << endl;
    of.close();
    return coll;
}

#endif