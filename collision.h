#ifndef _COLL_
#define _COLL_
#
#include "cross_section.h"
#include "particle.h"
#include "utility.h"
#include "espic_math.h"

typedef std::vector<std::pair<Particle, Particle>> CollPair;
using namespace ESPIC;

class Collisionpair {
public:    // In class all velocity except for the Update part are relative velocity 
    Collisionpair(Particle& particle, Real mass1, Real mass2, Real vtb)
    : pt(particle), 
    vx(pt.vxr()), vy(pt.vyr()), vz(pt.vzr()),
    vth(vtb),
    F1(mass1/(mass1+mass2)), F2(mass2/(mass1+mass2))
    {
        if (vx == 0) { theta = 0.5 * ESPIC::PI; }
        else { theta = atan2(sqrt(vx*vy+vz*vz), vx); }
        if (vy == 0) { 
            if (vz > 0) phi = 0.5 * ESPIC::PI;
            else phi = -0.5 * ESPIC::PI;
        } else phi = atan2(vz, vy);
        st = sin(theta);
        ct = cos(theta);
        sp = sin(phi);
        cp = cos(phi);
        wx = F1 * vx;
        wy = F1 * vy;
        wz = F1 * vz;
        
        energy = get_energy(vx, vy, vz);
        vel = sqrt(2.0 * energy);
    }

    void ParticleElasticCollision() 
    { 
        chi = acos(1.0 - 2.0*RG01());
        eta = ESPIC::PI2 * RG01();
        UpdateParticleVelInfo();
    }

    void ParticleExcitatinCollision(Real th) 
    {
        energy = fabs(energy - th);
        vel = sqrt(2.0 * energy);
        chi = acos(1.0 - 2.0 * RG01());
        eta = ESPIC::PI2 * RG01();
        UpdateParticleVelInfo();
    }

    void ParticleIonizationCollision(Real th)
    {
        Real en_ej, en_sc;
        Real vel_ej, chi_ej, eta_ej;
        // Real xej, yej, zej;
        // Real vxej, vyej, vzej; 
        // Real sc, cc, se, ce;

        energy = 0.5 * vel * vel;
        energy = fabs(energy - th);
        en_ej = 10.0 * tan(RG01() * atan(energy/20.));
        en_sc = fabs(energy - en_ej);
        vel = sqrt(2.0 * en_sc);
        vel_ej = sqrt(2.0 * en_ej);
        chi = acos(sqrt(en_sc / energy));
        chi_ej = acos(sqrt(en_ej / energy));
        eta = ESPIC::PI2 * RG01();
        eta_ej = eta + ESPIC::PI;

        Particle newelectron = pt;
        Particle newion = pt;
        UpdateParticleVelInfo(chi_ej, eta_ej, vel_ej, newelectron);
        EjectIonReaction(newion);
        newpp.emplace_back(std::make_pair(std::move(newelectron), std::move(newion)));
    }

    void ParticleIsotropicCollision()
    {
        chi = acos(1.0 - 2.0*RG01());
        eta = PI2 * RG01();
        UpdateParticleVelInfo();
    }

    void ParticleBackwardCollision()
    {
        chi = PI;
        eta = PI2 * RG01();
        UpdateParticleVelInfo();
    }

    const CollPair& ion_products() { return newpp; }



private:
    Particle& pt;
    Real vx, vy, vz;
    const Real vth;
    const Real F1, F2;
    Real chi, eta, theta, phi;
    Real wx, wy, wz;
    Real st, ct, sp, cp;
    CollPair newpp;
    Real vel, energy;
    // std::vector<Real> velbuffer;


void UpdateParticleVelInfo(Real chi_, Real eta_, Real vel_, Particle& particle)
{
    Real sc(sin(chi_)), cc(cos(chi_));
    Real se(sin(eta_)), ce(cos(eta_));
    vx = vel_ * (ct * cc - st * sc * ce);
    vy = vel_ * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    vz = vel_ * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
    particle.vx() = wx + F2*vx;
    particle.vy() = wy + F2*vy;
    particle.vz() = wz + F2*vz;
}

void UpdateParticleVelInfo()
{
    Real sc(sin(chi)), cc(cos(chi));
    Real se(sin(eta)), ce(cos(eta));
    vx = vel * (ct * cc - st * sc * ce);
    vy = vel * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    vz = vel * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
    pt.vx() = wx + F2*vx;
    pt.vy() = wy + F2*vy;
    pt.vz() = wz + F2*vz;
}

void EjectIonReaction(Particle& particle)
{
    Real vx, vy, vz;
    VelBoltzDistr(vth, vx, vy, vz);
    particle.vx() = vx; 
    particle.vy() = vy; 
    particle.vz() = vz;
}

};


#endif