#ifndef _SPEC_
#define _SPEC_
#include "particle.h"
#include "utility.h"

const Real bndens = 2414323.50534664;
const Real bmass = 72820.7;
class Species {
public:
    Species(const std::string& nm, Real m, Real q, Real t, Real n)
    : name(nm), mass(m), charge(q), temp(t),
    particles(new Particles()), totenergy(0)
    { }

    Species(const Species& other)
    : name(other.name), mass(other.mass), 
    charge(other.charge), temp(other.temp),
    particles(new Particles(*(other.particles))), 
    totenergy(0)
    { }
    
    ~Species()
    {
        delete particles;
    }

    void reserve_num_particles(Bigint n)
    { 
        particles->reserve(n);
    }

    void update_tot_energy()
    {
        totenergy = 0;
        for (size_t ip = 0; ip < particles->size(); ++ip)
        {
            Particle& pt = (*particles)[ip];
            totenergy += get_energy(pt.vx(), pt.vy(), pt.vz())*mass;
            totenergy += pt.lostenergy();
        }
    }

    std::string name;
    Real mass;
    Real charge;
    Real temp;
    Real ndens;
    class Particles* particles;
    Real totenergy; 
};

#endif