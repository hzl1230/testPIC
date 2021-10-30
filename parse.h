#ifndef PARSE_H
#define PARSE_H
#include <cstdio>
#include <iostream>
#include "cstdlib"
#include "str_split.h"

bool ParseLine(std::vector<std::string>& word, FILE* fp)
{
    if (!word.empty()) word.clear();
    
    char delim[16] = " \t\n\r,()";
    std::stringstream ss;
    std::string s;
    
    while (1) {
        char line[MaxLine];
        if (0 == fgets(line, MaxLine, fp)) { return false; }  // file ends

        ss.clear();
        ss.str(std::string(line));
        std::getline(ss, s, '!');
        if (s.empty()) { return true; }   // this is a comment line

        ss.clear();
        ss.str(s);
        std::getline(ss, s, '#');
        if (s.empty()) { return true; }   // this is a comment line

        std::vector<std::string> one = split(s, delim);
        if (one.empty()) { return true; } // this is a blank line

        if (one.back() == "&") {  // catenate next line
            word.insert(word.end(), one.begin(), one.end()-1);
            continue;
        }
        else {
            word.insert(word.end(), one.begin(), one.end());
            break;
        }
    }  
    return true;
}

#endif
