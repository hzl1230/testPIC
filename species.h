#ifndef _SPEC_
#define _SPEC_
#include "particle.h"
#include "utility.h"

class Species {
public:
    Species(const std::string& nm, Real m, Real q, Real t, Real n)
    : name(nm), mass(m), charge(q), temp(t),
    particles(new Particles())
    { }

    Species(const Species& other)
    : name(other.name), mass(other.mass), 
    charge(other), temp(other.temp),
    particles(new Particles(*(other.particles)))
    { }

    void reserve_num_particles(Bigint n)
    { 
        particles->reserve(n);
    }

    std::string name;
    Real mass;
    Real charge;
    Real temp;
    Real ndens;
    class Particles* particles; 
};

#endif