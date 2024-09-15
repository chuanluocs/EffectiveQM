#ifndef QASMREADER_H_INCLUDED
#define QASMREADER_H_INCLUDED

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <stdlib.h>
#include "gate.h"

using namespace std;

class QASMReader{
private:
    string fileName;
    int linenum;
    ifstream fin;
    set<int> vertices;

public:
    QASMReader(){linenum=0;}
    QASMReader(string filename):fin(filename){
        fileName = filename;
        linenum=0;
        if(!fin.is_open()){
             cerr<<fileName<<" does not exist"<<endl;
             exit(-1);
        }

    }
    ~QASMReader(){
        fin.close();
    }

    int getlinenum(){return linenum;}
    int getVerNum(){
        int max_element = -1;
        for(auto i:vertices){
            if(i>max_element) max_element = i;
        }
        return max_element+1;
    }
    /*read reduced quantum circuits containing only CNOT gates*/
    void parse(vector<gate>& gateList);
};

#endif // QASMREADER_H_INCLUDED
