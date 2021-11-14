#include "species.h"
#include "utility.h"
#include "cross_section.h"
#include "espic_math.h"
#include "collision.h"
#include <map>
#include <exception>
using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::pair;
typedef size_t size_type;

class Tile {
public:
    Tile(size_type n) 
    : background(new Background(0.01, 72820.7, 2414323.5)),
    npart(n),
    cross_section(new CrossSection("csection.in"))
    {
        InitSpecies();
        InitReaction();
        InitParticle();
    }
    ~Tile() 
    {
        delete cross_section;
        delete background;
        for (size_t ispec = 0; ispec < spec_arr.size(); ++ispec) 
            delete spec_arr[ispec];
        
        spec_arr.clear();
        spec_arr.shrink_to_fit();
    }

    Real UpdateMaxCollFreq(size_type reactid)
    {
        int specid = reaction_arr[reactid].first;
        Species* & species = spec_arr[specid];
        Particles& particles = *(species->particles);
        Reaction* & reaction = reaction_arr[reactid].second;

        Real numax = 0;
    #ifdef DEBUG
        std::ofstream of;
        if (specid == 0) of.open("out/vel_e.dat");
        else of.open("out/vel_p.dat");
    #endif
        // std::ofstream of("ion_E.dat",std::ofstream::app);
        
        for (size_type ip = 0; ip < particles.size(); ++ip)
        {
            Particle& pt = particles[ip];
            pt.gen_relative_vel(background->vth);
            
            vector<Real>& nu = pt.nu();
            if(!nu.empty()) nu.clear();
            Real energy = pt.rel_en(), nvt;
#ifdef DEBUG
            if (specid == 0) {
                of << pt.vx() << " " << pt.vxr() << " "
                << pt.vy() << " " << pt.vyr() << " "
                << pt.vz() << " " << pt.vzr() <<endl;
            }
            else {
                of << pt.vx() << " " << pt.vxr() << " "
                << pt.vy() << " " << pt.vyr() << " "
                << pt.vz() << " " << pt.vzr() <<endl;
            }
            of.close();
#endif
            vector<Real> info;
         
            info = reaction->en_cs(energy);
            
            Real nutot = 0;
            nvt = sqrt(2.0 * energy)*background->ndens;
            auto fnu = [=](auto& x){ return x*nvt; };
            transform(info.begin(), info.end(),back_inserter(nu),fnu);
            for (auto& nuj: nu) nutot += nuj;
            if(nutot > numax) numax = nutot;
        }
        return numax;
    }

    static void ParticleCollision(
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
            newpp.emplace_back(std::move(cop.ion_products()));
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


    class Background 
    {   
     public:
        Background(Real t, Real m, Real n)
        : temp(t), mass(m), ndens(n),
         vth(sqrt(2.0*t/m))
        { }

        const Real temp ;
        const Real mass;
        const Real ndens;
        const Real vth ;
    };
    const Background* background;
    vector<pair<int, Reaction*>> reaction_arr;
    vector<Species*> spec_arr;
    const size_type npart;

    const size_type num_species() const { return spec_arr.size(); }
    const size_type num_coll_species() const { return reaction_arr.size(); }

private:
    CrossSection* cross_section;
    std::map<string, size_type> spec_name_to_indx;
    
     
    void InitSpecies()
    {
        spec_arr.emplace_back(new Species("e", 1, -1, 1, 1));
        spec_arr.emplace_back(new Species("p", bmass, 1, 0.04, 1));
        for (size_t ispec = 0; ispec < spec_arr.size(); ++ispec) {
            Species* & species = spec_arr[ispec];
            string name = species->name;
            spec_name_to_indx[name] = ispec;
        }
    }

    void InitReaction()
    {
        for(int icsp = 0; icsp < cross_section->num_species(); ++icsp){
            Reaction* react = cross_section->react_arr[icsp];
            const string name = react->spec_name();
            int coll_spec_indx = -1;
            try {
                coll_spec_indx = spec_name_to_indx.at(name);
            }
            catch (const std::out_of_range& oor) {
                std::ostringstream oss;
                oss << "Unknown species \"" << name << "\" given to \"cross_section\" command in [csection.in]";
                espic_error(oss.str());
            }
            std::pair<int, Reaction*> temp(std::make_pair(coll_spec_indx, std::move(react)));
            reaction_arr.emplace_back(std::move(temp));
        }
    }

    void InitParticle()
    {
        for (size_t ispec = 0; ispec < spec_arr.size(); ++ispec) 
        {
            Species* & species = spec_arr[ispec];
            vector<Real> x(npart), vx(npart);
            vector<Real> y(npart), vy(npart);
            vector<Real> z(npart), vz(npart);
            Real temp = species->temp, m = species->mass;
            Particles* & particles = species->particles;
            particles = new Particles(x,y,z,vx,vy,vz);
            
            for (size_t ip = 0; ip < particles->size(); ++ip)
            {
                Particle& pt = (*particles)[ip];
                VelBoltzDistr(sqrt(2.0*temp/m), pt.vx(), pt.vy(), pt.vz());
                // Real vxb, vyb, vzb;
                // VelBoltzDistr(background->vth, vxb, vyb, vzb);
                // RelativeVelocity(pt, vxb, vyb, vzb);
            } 
        }
    }

    size_type FindReactionIndx(size_type specid)
    {
        for(size_type ire = 0; ire < reaction_arr.size(); ++ire)
            if (reaction_arr[ire].first == (int)specid)
                return ire;
        espic_error("No corresponding Particle Species");
        return -1;
    }   
};