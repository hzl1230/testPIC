#ifndef _CCSECTION_H
#define _CCSECTION_H
#
#include "utility.h"
#include "parse.h"


typedef std::pair<std::string,class Reaction*> pair_name_react;
extern Real kTe0;
Real kTe0;
class Reaction
{
public:
    Reaction() : infile("e_Ar.dat"), 
                 info_size(0), 
                 name("default"),
                 reaction_id(0),
                 arr_length(3),
                 threshold(0)
    { }

    Reaction (std::string file, int reaction_type_number, std::vector<std::string> type,
                std::string rname, int i) 
    : infile(file),
    info_size(reaction_type_number),
    name(rname),
    types(type),
    reaction_id(i),
    arr_length(3),
    threshold(0)
    {
        info.reserve(info_size);
        FILE* fp = fopen(file.c_str(), "r");
        std::vector<std::string> line;
        if (NULL == fp) {
            std::ostringstream oss;
            oss << "Cannot read file [" << infile << "]";
            espic_error(oss.str());
        }
        while(ParseLine(line, fp)){
            if(!info.empty()) info.clear();
            if (line.empty()) continue; 
            else if("bin"==line.at(0)) {
                arr_length = atoi(line[1].c_str()); 
                info_arr.reserve(arr_length);
                energy.reserve(arr_length);
                if ("de" == line.at(2))
                    de_ = (Real)atof(line[3].c_str());
            }
            else if("threshold"==line.at(0)){
                transform(std::move(line).begin()+1, std::move(line).end(), back_inserter(threshold), toReal);
                resize_threshold();
            }
            else{
                energy.emplace_back((Real)atof(line[0].c_str()));
                transform(std::move(line).begin()+1, std::move(line).begin()+info_size+1, back_inserter(info), toReal);
                info_arr.emplace_back(std::move(info));
            }
        }
    }

    const std::vector<Real>& en_cs(Real en)
    {
        // Get cs_info through 1DLinearInterpoltationMethod
        if(en <= de_) { 
            info ={0.,0.,0.};
            return info;
        } 
        Real deinv(1./de_);
        Real ei(en*deinv-1) ;
        int elo = static_cast<int>(ei);
        if (elo > arr_length) {
            std::string var_name("info_arr_size");
            espic_error(out_bound_info(var_name, infile));
        }
        Real wgt[2], temp;
        // Cal node weight 
        wgt[1] = ei - elo;
        wgt[0] = 1 - wgt[1];
        
        if(!info.empty()) info.clear();
        for(int i = 0; i < info_size; ++i) {
            temp = wgt[0]*info_arr[elo][i] + wgt[1]*info_arr[elo+1][i];
            info.emplace_back(std::move(temp));
        }
        return info;
    }


    const int size() const { return arr_length; }
    const int isize() const { return info_size; }
    const std::vector<Real> th() const { return threshold; }
    const std::string spec_name() const { return name; }
    const int r_index() { return reaction_id; }
    const Real de() { return de_; }
    const std::vector<Real> csection(int i) { return info_arr[i]; }
    const std::vector<std::string> get_types() { return types; }
    
private:
    std::string infile;
    const int info_size;
    const std::string name;
    std::vector<std::string> types;
    int reaction_id, arr_length;
    Real de_;
    std::vector<Real> threshold;
    std::vector<Real> info;
    std::vector<std::vector<Real>> info_arr; // info_arr of every energy
    std::vector<Real> energy;                // energy_arr of total energy bin
    
    void resize_threshold() 
    { 
        if(threshold.size() < static_cast<size_t>(info_size - 1)) {
            std::size_t di = info_size - 1 - threshold.size();
            while(di-- > 0) 
                threshold.emplace_back(0);  
        }
    }
};


class CrossSection
{
    friend class Reaction;

public:
    CrossSection(const std::string &file="csection.in")
    : infile(file) {
        read_input_cross_section();
        react_arr.reserve(num_species());
        get_reaction();
    }

    ~CrossSection()
    {
        for(std::size_t irea = 0; irea < react_arr.size(); ++irea) 
            delete react_arr[irea];
    
        react_arr.clear();
        react_arr.shrink_to_fit();
    }

    std::vector<class Reaction*> react_arr;

    const std::vector<std::string>& name_species() { return name; }
    const std::vector<std::string>& files() { return reaction_file; }
    const std::vector<int>& num_react() const { return reaction_type_number; }
    int num_species() const { return reaction_species_number; }

    /*---------------------begin private method-----------------------*/

private:
    const std::string infile;
    int reaction_species_number;
    std::vector<std::string> name;
    std::vector<std::string> reaction_file;
    std::vector<int> reaction_type_number;
    std::vector<std::string> reaction_types; 

    CrossSection();

    void get_reaction() 
    { 
        std::string file, spec;
        std::vector<std::string> types;
        int type_num;
        for(int i = 0; i < reaction_species_number; ++i) {
            file = reaction_file[i];
            type_num = reaction_type_number[i];
            spec = name[i];
            std::copy(reaction_types.begin(), reaction_types.begin()+type_num, back_inserter(types));
            react_arr.emplace_back(new Reaction(file, type_num, types, spec, i));
            reaction_types.erase(reaction_types.begin(), reaction_types.begin()+type_num);
            types.clear();
        }
    }

    void read_input_cross_section()
    {
        FILE *fp = fopen(infile.c_str(), "r");
        std::string dir("reaction/");
        if (NULL == fp)
        {
            std::ostringstream oss;
            oss << "Cannot read file [" << infile << "]";
            espic_error(oss.str());
        }
        std::vector<std::string> word;
        // std::string cmd;
        while (ParseLine(word, fp))
        {
            if (word.empty())
                continue;
            else if ("total" == word.at(0)) {
                reaction_species_number = static_cast<int>(atoi(word[1].c_str()));
            }
            else if ("pairs" == word.at(0)) {
                name.emplace_back(word[1]);
                word.erase(word.begin(), word.begin() + 2);
                while (word.size() > 0) {
                    if ("reaction" == word.at(0))
                    {
                        int nreact = atoi(std::move(word[1]).c_str());
                        reaction_type_number.emplace_back(nreact);
                        word.erase(word.begin(), word.begin() + 2);
                        for(int i = 0; i < nreact; ++i) 
                            reaction_types.emplace_back(word[i]);
                        word.erase(word.begin(), word.begin()+nreact);
                    }
                    else
                    {
                        std::string cmd(word[0]);
                        espic_error(unknown_cmd_info(cmd, infile));
                    }
                    if ("dir" == word.at(0))
                    {
                        reaction_file.emplace_back(dir+word[1]);
                        word.erase(word.begin(), word.begin() + 2);
                    }
                    else
                    {
                        std::string cmd(word[0]);
                        espic_error(unknown_cmd_info(cmd, infile));
                    }
                }
            }
            else if ("aid_param" == word.at(0)) {
                kTe0 = (Real)atof(word[1].c_str());
            }
            else {
                std::string cmd(word[0]);
                espic_error(unknown_cmd_info(cmd, infile));
            }
        }
    }
};

#endif