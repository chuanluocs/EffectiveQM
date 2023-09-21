#include <cassert>
#include <queue>
#include <iostream>
#include <utility>
#include <vector>
#include <list>
#include <map>
#include <unordered_map>
#include <ctime>
#include <bitset>
#include <chrono>
//#include <io.h>

#include "gate.h"
#include "QASMReader.h"
#include "CouplingGraph.h"
#include "WeightedGate.h"
#include <algorithm>
#include <random>
#define DEBUG

using namespace std;

// components

// GREEDY
#define GREEDY // ues dynamic greedy tech
//--------------------GREEDY--------------------

#define RANDOM_INIT // ues random init
//--------------------INIT--------------------


// parameters
// attention: parameters defined below will be replaced by the input command, 
// the initial value just show the default value
int IN_LIMIT = 6; // S_f+1 in paper, 6 is the default value
int MAX_LOOK_SIZE = 4; // theta-1 in paper, 4 is the default value
int C_MAX = 500; // I_f in paper, 500 is the default value


// Log file
ofstream logout;

// global random seed, defined in input command
int global_algorithm_seed;
mt19937 global_random_device;

string selected_device_enum;


chrono::system_clock::time_point start_time;
chrono::system_clock::time_point current_time;
// cutoff time with s, not uesd in final version
const double cutoff_time=3600;

int ** interactMat;
int * qubit_gate_cnt;

// debug function
void debugMapping(const vector<int>&qubits){
    logout<<"qubits_mapping:"<<endl;
    for(int i=0;i<qubits.size();i++){
        logout<<i<<"->"<<qubits[i]<<endl;
    }
    logout<<"-------------------------------"<<endl;
}


// mapping function
void getQubitsSequence(vector<int>&qubits,int nqubits,const CouplingGraph &cg){
    #ifdef RANDOM_INIT
        shuffle(qubits.begin(),qubits.end(),global_random_device);
        return ;
    #endif
    logout<<"[ERROR] we should ues random init."<<endl;
    exit(-1);
}

// ues shuffle to perurbate the mapping
void PerurbationByShuffle(
    vector<int> & qubits, 
    vector<int> & locations,
    int nlocs, 
    int nqubits,
    const CouplingGraph &cg)
{
    getQubitsSequence(qubits,nqubits,cg);
    for(int i = 0;i < nlocs; i++)
    {
        if(qubits[i] == -1)
            continue;
        locations[qubits[i]] = i;
    }
}

// ues randomnazation to get initial mapping
void randInitMapping(vector<int> & qubits, vector<int> & locations,int nlocs, int nqubits,const CouplingGraph &cg)
{
    for(int i = 0;i < nlocs; i++)
    {
        if(i < nqubits)
            qubits[i] = i;
        else
            qubits[i] = -1;
    }
    PerurbationByShuffle(qubits,locations,nlocs,nqubits,cg);
}




int ** initInteractionMatrix(int nqubits)
{
    int ** interactionMat = new int*[nqubits];
    
    qubit_gate_cnt = new int[nqubits];
    for(int i=0;i<nqubits;i++){
        qubit_gate_cnt[i]=0;
    }

    for(int i=0;i<nqubits;i++)
        interactionMat[i] = new int[nqubits];
    for(int i=0;i<nqubits;i++)
        for(int j=0;j<nqubits;j++)
            interactionMat[i][j] = 0;
    return interactionMat;
}
void clearInteractionMatrix(int ** interationMat,int nqubits)
{
    for(int i=0;i<nqubits;i++)
        delete [] interationMat[i];
    delete [] qubit_gate_cnt;
    delete [] interationMat;
}


void getInteractionMatrix(const vector<gate> gateList, int ** interactMat, int startIndex = 0)
{
    for(size_t i=startIndex;i<gateList.size();i++)
    {
        gate g1 = gateList[i];
        int c1 = g1.control;
        int t1 = g1.target;
        interactMat[c1][t1]++;
        interactMat[t1][c1]++;
        qubit_gate_cnt[c1]++;
        qubit_gate_cnt[t1]++;
        
    }
    
}


/**
* Partition the logic circuit into n gate arrays,
* lines[i] holds all the gate along the ith qubit.
*/
void genLines(vector<list<gate> > & lines, vector<gate> & gateList)
{
    for(unsigned int i=0; i<gateList.size(); i++)
    {
        gate g1 = gateList[i];
        int c1 = g1.control;
        int t1 = g1.target;
        lines[c1].push_back(g1);
        lines[t1].push_back(g1);
    }
}

// find the active gates in current logic circuit
void findActivegates(vector<list<gate> > & lines,vector<gate> & activeGates,vector<gate> & phyCir,set<int> & activeLines, vector<int> locations, int ** dist)
{
    
    while(!activeLines.empty())
    {
        for(auto it = activeLines.begin();it != activeLines.end(); it++)
        {
            int i = * it;
            int j = 0;
            if(!lines[i].empty())
            {
                gate g1 = lines[i].front();
                int c = g1.control;
                int t = g1.target;
                if(c == i)
                    j = t;
                else
                    j = c;


                gate g2 = lines[j].front();
                if(g1==g2) //active gate
                {
                    int locOfC = locations[c];
                    int locOfT = locations[t];
                    //executable gate
                    if(dist[locOfC][locOfT] == 1)
                    {
                        gate pg(locOfC,locOfT, g1.weight);
                        lines[i].pop_front();
                        lines[j].pop_front();
                        phyCir.push_back(pg);
                        activeLines.insert(j);
                    }
                    else
                    {
                        activeGates.push_back(g1);
                        activeLines.erase(i);
                        activeLines.erase(j);
                    }

                }
                else
                {
                    activeLines.erase(i);
                }
                break;
            }
            else
            {
                activeLines.erase(i);
                break;
            }
        }
    }
}


// find the active gates in current logic circuit
// ues additonal variable to recover the state 
void findActivegates_recover(vector<list<gate> > & lines,vector<gate> & activeGates,vector<gate> & phyCir,set<int> & activeLines, vector<int> locations, int ** dist,vector<list<gate> > & re_line)
{

    while(!activeLines.empty())
    {
        for(auto it = activeLines.begin();it != activeLines.end(); it++)
        {
            int i = * it;
            int j = 0;
            if(!lines[i].empty())
            {
                gate g1 = lines[i].front();
                int c = g1.control;
                int t = g1.target;
                if(c == i)
                    j = t;
                else
                    j = c;
                //始终是i当前，j目标


                gate g2 = lines[j].front();
                if(g1==g2) //active gate
                {
                    int locOfC = locations[c];
                    int locOfT = locations[t];
                    //executable gate
                    if(dist[locOfC][locOfT] == 1)
                    {
                        gate pg(locOfC,locOfT, g1.weight);
                        lines[i].pop_front();
                        lines[j].pop_front();
                        re_line[i].push_back(g1);
                        re_line[j].push_back(g1);
                        phyCir.push_back(pg);
                        activeLines.insert(j);
                    }
                    else
                    {
                        activeGates.push_back(g1);
                        activeLines.erase(i);
                        activeLines.erase(j);
                    }

                }
                else
                {
                    activeLines.erase(i);
                }
                break;
            }
            else
            {
                activeLines.erase(i);
                break;
            }
        }
    }
}

/**
*Get the number of executable gates under the mapping
*lines:gate sequence of each logical qubit/line
*activeGate: array of active gates
*/
int numOfExecuted(vector<list<gate> > & lines, const vector<int> & locations, int ** dist)
{
    int exeNum = 0;
    set<int> activeLines;
    //Save gate deleted from lines for recovery
    vector<gate> executedGates;

    for(size_t i=0 ; i<lines.size(); i++)
        activeLines.insert(i);

    while(!activeLines.empty())
    {
        for(set<int>::iterator it = activeLines.begin();it != activeLines.end(); it++)
        {
            int i = * it;
            int j = 0;
            if(!lines[i].empty()){
                gate g1 = lines[i].front();
                int c = g1.control;
                int t = g1.target;
                if(c == i)
                    j = t;
                else
                    j = c;
                gate g2 = lines[j].front();
                if(g1==g2) //active gate
                {
                    int locOfC = locations[c];
                    int locOfT = locations[t];
                    //executable gate
                    if(dist[locOfC][locOfT] == 1)
                    {
                        exeNum++;
                        lines[i].pop_front();
                        lines[j].pop_front();
                        activeLines.insert(j);
                        executedGates.push_back(g1);
                    }
                    else
                    {
                        activeLines.erase(i);
                        activeLines.erase(j);
                    }

                }
                else
                {
                    // logout<<"!="<<endl;
                    activeLines.erase(i);
                }
                break;
            }
            else
            {
                activeLines.erase(i);
                break;
            }
        }
    }

    //Recover lines
    for(int i=executedGates.size()-1;i>=0;i--)
    {
        gate g1 = executedGates[i];
        int c = g1.control;
        int t = g1.target;
        lines[c].push_front(g1);
        lines[t].push_front(g1);
    }
    return exeNum;
}

set<int>exe_st;
// calculate the number of executed gates and the number of front gates
// and log the executed gates
pair<int,int> numOfExecutedAndFrontGate_3_3_log(vector<list<gate> > & lines, const vector<int> & locations, int ** dist)
{
    int exeNum = 0;
    int gateNum = 0;
    int fixed = 0;
    set<pair<int,int>> frontGate;
    set<int> activeLines;
    //Save gate deleted from lines for recovery
    vector<gate> executedGates;

    for(size_t i=0 ; i<lines.size(); i++)
        activeLines.insert(i);

    while(!activeLines.empty())
    {
        for(set<int>::iterator it = activeLines.begin();it != activeLines.end(); it++)
        {
            int i = * it;
            int j = 0;
            if(!lines[i].empty()){
                gate g1 = lines[i].front();
                int c = g1.control;
                int t = g1.target;
                if(c == i)
                    j = t;
                else
                    j = c;
                gate g2 = lines[j].front();
                if(g1==g2) //active gate
                {
                    int locOfC = locations[c];
                    int locOfT = locations[t];
                    //executable gate
                    if(dist[locOfC][locOfT] == 1)
                    {
                        // logout<<"[INFO] Executed gate: "<<c<<"->"<<t<<" "<<g1.weight<<endl;
                        exe_st.insert(g1.weight);
                        exeNum++;
                        lines[i].pop_front();
                        lines[j].pop_front();
                        activeLines.insert(j);
                        executedGates.push_back(g1);
                    }
                    else
                    {
                        frontGate.insert({min(c,t),max(c,t)});
                        gateNum++;
                        activeLines.erase(i);
                        activeLines.erase(j);
                    }

                }
                else
                {
                    // logout<<"!="<<endl;
                    activeLines.erase(i);
                }
                break;
            }
            else
            {
                activeLines.erase(i);
                break;
            }
        }
    }
    if(gateNum!=frontGate.size()){
        logout<<"[ERROR] gateNum!=frontGate.size()"<<endl;
        exit(-1);
    }

    //Recover lines
    for(int i=executedGates.size()-1;i>=0;i--)
    {
        gate g1 = executedGates[i];
        int c = g1.control;
        int t = g1.target;
        lines[c].push_front(g1);
        lines[t].push_front(g1);
    }
    return {exeNum,gateNum};
}


// remove the executed gates from the logic circuit
void checkActivegates(vector<list<gate> > & lines,vector<gate> & activeGates,vector<gate> & phyCir,set<int> & activeLines, vector<int> locations, int ** dist)
{
    for(vector<gate>::iterator it = activeGates.begin(); it!=activeGates.end(); )
    {
        gate g = * it;
        int c = g.control;
        int t = g.target;
        int locOfC = locations[c];
        int locOfT = locations[t];
        //execute gate
        if(dist[locOfC][locOfT] == 1)
        {
            gate pg(locOfC,locOfT, g.weight);
            lines[c].pop_front();
            lines[t].pop_front();
            phyCir.push_back(pg);
            activeLines.insert(c);
            activeLines.insert(t);
            it = activeGates.erase(it);
        }
        else
        {
            it++;
        }
    }
}

// remove the executed gates from the logic circuit
void checkActivegates_recover(vector<list<gate> > & lines,vector<gate> & activeGates,vector<gate> & phyCir,set<int> & activeLines, vector<int> locations, int ** dist,vector<list<gate> > & re_line)
{
    for(vector<gate>::iterator it = activeGates.begin(); it!=activeGates.end(); )
    {
        gate g = * it;
        int c = g.control;
        int t = g.target;
        int locOfC = locations[c];
        int locOfT = locations[t];
        //execute gate
        if(dist[locOfC][locOfT] == 1)
        {
            gate pg(locOfC,locOfT, g.weight);
            lines[c].pop_front();
            lines[t].pop_front();

            re_line[c].push_back(g);
            re_line[t].push_back(g);

            phyCir.push_back(pg);
            activeLines.insert(c);
            activeLines.insert(t);
            it = activeGates.erase(it);
        }
        else
        {
            it++;
        }
    }
}

// calculate the function first section of SWAP gate
int effOfSWAP(pair<int,int> & swapGate, gate g1, vector<int> locations,  int ** dist)
{
    //two logical qubits of g1
    int i1 = g1.control;
    int j1 = g1.target;

    //two physical qubits of g1 before SWAP gate
    int locOfI1 = locations[i1];
    int locOfJ1 = locations[j1];
    if(locOfI1 > locOfJ1)
        swap(locOfI1, locOfJ1);

    #ifdef DEBUG
    if(locOfI1==-1 || locOfJ1==-1)
        logout<<"[BUG] error in effOfSWAP, active gate must have two physical qubits"<<endl;
    #endif // DEBUG

    // two physical qubits of swapGate
    int s1 = swapGate.first;
    int s2 = swapGate.second;
    if(s1 > s2)
        swap(s1, s2);

    if(s1 == locOfI1 && s2 == locOfJ1)
        return 0;
    else if(s1 == locOfI1)
        // two physical qubits of g1 after SWAP gate
        return dist[locOfI1][locOfJ1] - dist[s2][locOfJ1];
    else if(s1 == locOfJ1)
        return dist[locOfI1][locOfJ1] - dist[locOfI1][s2];
    else if(s2 == locOfI1)
        return dist[locOfI1][locOfJ1] - dist[s1][locOfJ1];
    else if(s2 == locOfJ1)
        return dist[locOfI1][locOfJ1] - dist[locOfI1][s1];
    else
        return 0;
}


/**
* Generate relevant SWAPs or Cnots which have at least one common qubit with some active gate
*/
void genRelevantGates(const CouplingGraph * p_cg,const vector<gate> & activeGates,set<weightedGate> & candiGates,const  vector<int> & locations,int ** dist )
{
    for(unsigned int i=0; i<activeGates.size(); i++)
    {
        gate g1 = activeGates[i];
        int c = g1.control;
        int t = g1.target;
        int cAndT[2] = {locations[c], locations[t]};
        if(dist[cAndT[0]][cAndT[1]] == 2)
        {
            weightedGate cnot1("CNOT",cAndT[0],cAndT[1]);
            candiGates.insert(cnot1);
        }

        //for all possible SWAPs on c and t
        for(int ct=0; ct < 2;ct++)
        {
            int locOfC = cAndT[ct];
            for(unsigned int j=0; j<p_cg->adjList[locOfC].size(); j++)
            {
                int s1 = p_cg->adjList[locOfC][j];
                weightedGate swap1("SWAP");

                if(locOfC < s1)
                {
                    swap1.q1 = locOfC;
                    swap1.q2 = s1;
                }
                else
                {
                    swap1.q1 = s1;
                    swap1.q2 = locOfC;
                }
                candiGates.insert(swap1);
            }
        }
    }
}

/**
* Get the first nth layers of logical circuit
*/
vector<vector<gate> > getFrontLayers(vector<list<gate> > & lines, int layerDepth)
{
    vector<vector<gate> > layers;
    int nqubits = (int) lines.size();
    vector<int> layer(nqubits,-1); //depth of each line during iteration
    vector<list<gate>::iterator> itOfLines(nqubits);
    for(int i=0;i<nqubits;i++)
    {
        itOfLines[i] = lines[i].begin();
    }

    bool isAllLineEmpty = false;
    bool isAllOverWindow = false;
    while(!isAllLineEmpty && !isAllOverWindow)
    {
        isAllLineEmpty = true;
        isAllOverWindow = true;
        for(int i=0;i<nqubits;i++)
        {
            if((itOfLines[i]!=lines[i].end()) && layer[i]<layerDepth)
            {
                isAllLineEmpty = false;
                //isAllOverWindow = false;
                gate g1 = * itOfLines[i];
                int c_g1 = g1.control;
                int t_g1 = g1.target;
                if((*itOfLines[c_g1]) == (*itOfLines[t_g1]))
                {
                    int layerOfG1 = max(layer[c_g1],layer[t_g1])+1;
                    layer[c_g1]=layerOfG1;
                    layer[t_g1]=layerOfG1;
                    if(layerOfG1 <= layerDepth)
                    {
                        if(layers.size()<= (size_t)layerOfG1)
                            layers.push_back(vector<gate>());
                        layers[layerOfG1].push_back(g1);
                        itOfLines[c_g1]++;
                        itOfLines[t_g1]++;
                        isAllOverWindow = false;
                    }

                }
            }
        }
    }
    return layers;
}


/**
* Update the mapping between physical and logical qubits with a SWAP gate
*/
void updateMapping(pair<int,int> & swapGate, vector<int> & qubits, vector<int> & locations )
{
    //Two physical qubits of swapGate
    int pFirst = swapGate.first;
    int pSecond = swapGate.second;
    //Two logical qubits of swapGate
    int lFirst = qubits[pFirst];
    int lSecond = qubits[pSecond];

    //Update qubits array
    int temp = 0;
    temp = qubits[pFirst];
    qubits[pFirst] = qubits[pSecond];
    qubits[pSecond] = temp;

    //update locations array
    if(lFirst != -1) // corresponding to a unoccupied physical qubit
        locations[lFirst] = pSecond;
    if(lSecond != -1)
        locations[lSecond] = pFirst;
    // logout<<"[TEST] ues swap gate: "<<pFirst<<"<->"<<pSecond<<endl;
    // debugMapping(qubits);
}


/**
* Compute cost function for  bridge CNOT gate (dist == 2)
*
*/
int costOfInsertedCnot(vector<list<gate> > & lines, pair<int,int> & cnotGate, const vector<int> & qubits,const vector<int> & locations,  int ** dist)
{
    int totalEffect = 0;
    int c = qubits[cnotGate.first];
    int t = qubits[cnotGate.second];
    gate g1 = lines[c].front();
    lines[c].pop_front();
    gate g2 = lines[t].front();
    lines[t].pop_front();
    if(!(g1 == g2) || !(g1.control==c && g1.target==t))
    {
        logout<<"[BUG] error in 980" <<endl;
        logout<<g1.control<<","<<g1.target<<","<<g1.weight<<endl;
        logout<<g2.control<<","<<g2.target<<","<<g2.weight<<endl;
        exit(-1);
    }
    totalEffect = numOfExecuted(lines,locations,dist);
    lines[c].push_front(g1);
    lines[t].push_front(g1);
    // logout<<"[TEST] function value of cnot: "<<1<<' '<<1+totalEffect<<endl;
    return totalEffect+1+1; // and itself
}

/**
* Compute cost function for each SWAP
* propagate positive effect across layers until an negative occurs
*/
bool has_not_zero;
int costOfInsertedSwap(vector<list<gate> > & lines, pair<int,int> & swapGate,const vector<vector<gate> > & layers, vector<int> qubits,vector<int> locations,  int ** dist)
{
    bool hasExecutable = false;
    bool isNegative = false;
    int funcV1=0,funcV2=0;

    for(size_t i=0;i<layers.size();i++)
    {
        int effOfLayer = 0;
        for(size_t j=0;j<layers[i].size();j++)
        {
            gate g1 = layers[i][j];
            int effectOnG1 = effOfSWAP(swapGate, g1, locations, dist);
            if(i==0) //gate of layer 0
            {
                int loc1_g1 = locations[g1.control];
                int loc2_g1 = locations[g1.target];
                int dist_g1 = dist[loc1_g1][loc2_g1];
                if(dist_g1 == 2 && effectOnG1==1)
                    hasExecutable = true;
            }
            if(effectOnG1 < 0) //Stop propagation at this layer
                isNegative = true;
            effOfLayer += effectOnG1;
        }
        if(isNegative == true)
            break;
        else
            funcV1 += effOfLayer;
        
    }
    if(hasExecutable)
    {
        has_not_zero=true;
        updateMapping(swapGate,qubits,locations);
        int exeGateCounts = numOfExecuted(lines,locations,dist);
        funcV2 = exeGateCounts;
    }
    return funcV1+funcV2;
}



pair<int,int> costOfInsertedCnotSplit(vector<list<gate> > & lines, pair<int,int> & cnotGate, 
                                           const vector<int> & qubits,const vector<int> & locations,  int ** dist)
{
    int totalEffect = 0;
    int c = qubits[cnotGate.first];
    int t = qubits[cnotGate.second];
    gate g1 = lines[c].front();
    lines[c].pop_front();
    gate g2 = lines[t].front();
    lines[t].pop_front();
    if(!(g1 == g2) || !(g1.control==c && g1.target==t))
    {
        logout<<"[BUG] error in 980" <<endl;
        logout<<g1.control<<","<<g1.target<<","<<g1.weight<<endl;
        logout<<g2.control<<","<<g2.target<<","<<g2.weight<<endl;
        exit(-1);
    }
    totalEffect = numOfExecuted(lines,locations,dist);
    lines[c].push_front(g1);
    lines[t].push_front(g1);
    // logout<<"[TEST] function value of cnot: "<<1<<' '<<1+totalEffect<<endl;
    // logout<<"[TEST] qubit of cnot: "<<cnotGate.first<<' '<<cnotGate.second<<endl;
    return {1,totalEffect+1}; // and itself
}


pair<int,int> numOfExecutedAndFrontGate(vector<list<gate> > & lines, const vector<int> & locations, int ** dist)
{
    int exeNum = 0;
    int gateNum = 0;
    set<pair<int,int>> frontGate;
    set<int> activeLines;
    //Save gate deleted from lines for recovery
    vector<gate> executedGates;

    for(size_t i=0 ; i<lines.size(); i++)
        activeLines.insert(i);

    while(!activeLines.empty())
    {
        for(set<int>::iterator it = activeLines.begin();it != activeLines.end(); it++)
        {
            int i = * it;
            int j = 0;
            if(!lines[i].empty()){
                gate g1 = lines[i].front();
                int c = g1.control;
                int t = g1.target;
                if(c == i)
                    j = t;
                else
                    j = c;
                gate g2 = lines[j].front();
                if(g1==g2) //active gate
                {
                    int locOfC = locations[c];
                    int locOfT = locations[t];
                    //executable gate
                    if(dist[locOfC][locOfT] == 1)
                    {
                        // logout<<"[INFO] Executed gate: "<<c<<"->"<<t<<" "<<g1.weight<<endl;
                        exeNum++;
                        lines[i].pop_front();
                        lines[j].pop_front();
                        activeLines.insert(j);
                        executedGates.push_back(g1);
                    }
                    else
                    {
                        frontGate.insert({min(c,t),max(c,t)});
                        gateNum++;
                        activeLines.erase(i);
                        activeLines.erase(j);
                    }

                }
                else
                {
                    // logout<<"!="<<endl;
                    activeLines.erase(i);
                }
                break;
            }
            else
            {
                activeLines.erase(i);
                break;
            }
        }
    }
    if(gateNum!=frontGate.size()){
        logout<<"[ERROR] gateNum!=frontGate.size()"<<endl;
        exit(-1);
    }

    // for(auto i:executedGates){
    //     logout<<i.control<<' '<<i.target<<endl;
    // }
    // logout<<"-----------------"<<endl;
    //Recover lines
    for(int i=executedGates.size()-1;i>=0;i--)
    {
        gate g1 = executedGates[i];
        int c = g1.control;
        int t = g1.target;
        lines[c].push_front(g1);
        lines[t].push_front(g1);
    }
    return {exeNum,gateNum};
}


tuple<int,int,int> costOfInsertedCnot_L0_Nexe_Frontgate(vector<list<gate> > & lines, pair<int,int> & cnotGate, 
                                           const vector<int> & qubits,const vector<int> & locations,  int ** dist)
{
    int c = qubits[cnotGate.first];
    int t = qubits[cnotGate.second];
    gate g1 = lines[c].front();
    lines[c].pop_front();
    gate g2 = lines[t].front();
    lines[t].pop_front();
    if(!(g1 == g2) || !(g1.control==c && g1.target==t))
    {
        logout<<"[BUG] error in 980" <<endl;
        logout<<g1.control<<","<<g1.target<<","<<g1.weight<<endl;
        logout<<g2.control<<","<<g2.target<<","<<g2.weight<<endl;
        exit(-1);
    }
    // logout<<"begin numOfExecuted... "<<endl;
    pair<int,int> gateCost = numOfExecutedAndFrontGate(lines,locations,dist);
    lines[c].push_front(g1);
    lines[t].push_front(g1);
    // logout<<"[TEST] function value of cnot: "<<1<<' '<<1+totalEffect<<endl;
    // logout<<"[TEST] qubit of cnot: "<<cnotGate.first<<' '<<cnotGate.second<<endl;
    return {1,gateCost.first+1,gateCost.second}; // and itself
}

pair<int,int> costOfInsertedSwapSplit(vector<list<gate> > & lines, pair<int,int> & swapGate,const vector<vector<gate> > & layers, 
                                            vector<int> qubits,vector<int> locations,  int ** dist)
{
    bool hasExecutable = false;
    bool isNegative = false;
    int funcV1=0,funcV2=0;
    // logout<<"layers.size() = "<<layers.size()<<endl;

    for(size_t i=0;i<layers.size();i++)
    {
        int effOfLayer = 0;
        for(size_t j=0;j<layers[i].size();j++)
        {
            gate g1 = layers[i][j];
            // logout<<"begin effOfSWAP"<<endl;
            int effectOnG1 = effOfSWAP(swapGate, g1, locations, dist);
            // logout<<"end effOfSWAP"<<endl;
            if(i==0) //gate of layer 0
            {
                int loc1_g1 = locations[g1.control];
                int loc2_g1 = locations[g1.target];
                int dist_g1 = dist[loc1_g1][loc2_g1];
                if(dist_g1 == 2 && effectOnG1==1)
                    hasExecutable = true;
            }
            if(effectOnG1 < 0) //Stop propagation at this layer
                isNegative = true;
            effOfLayer += effectOnG1;
        }
        if(isNegative == true)
            break;
        else
            funcV1 += effOfLayer;
        
    }
    // logout<<"funcV1: "<<funcV1<<endl;
    if(hasExecutable)
    {
        has_not_zero=true;
        updateMapping(swapGate,qubits,locations);
        int exeGateCounts = numOfExecuted(lines,locations,dist);
        funcV2 = exeGateCounts;
    }
    
    return {funcV1,funcV2};
}

// return the cost of inserted swap
tuple<int,int,int> costOfInsertedSwap_L0_Nexe_Frontgate(vector<list<gate> > & lines, pair<int,int> & swapGate,const vector<vector<gate> > & layers, 
                                            vector<int> qubits,vector<int> locations,  int ** dist)
{
    bool hasExecutable = false;
    bool isNegative = false;
    int funcV1=0,funcV2=0,funcV3=0;
    for(size_t i=0;i<layers.size();i++)
    {
        int effOfLayer = 0;
        for(size_t j=0;j<layers[i].size();j++)
        {
            gate g1 = layers[i][j];
            // logout<<"begin effOfSWAP"<<endl;
            int effectOnG1 = effOfSWAP(swapGate, g1, locations, dist);
            // logout<<"end effOfSWAP"<<endl;
            if(i==0) //gate of layer 0
            {
                int loc1_g1 = locations[g1.control];
                int loc2_g1 = locations[g1.target];
                int dist_g1 = dist[loc1_g1][loc2_g1];
                if(dist_g1 == 2 && effectOnG1==1)
                    hasExecutable = true;
            }
            effOfLayer += effectOnG1;
        }
        funcV1 += effOfLayer;
    }
    // logout<<"funcV1: "<<funcV1<<endl;
    if(hasExecutable)
    {
        has_not_zero=true;
        updateMapping(swapGate,qubits,locations);
        pair<int,int> exeGateCounts = numOfExecutedAndFrontGate(lines,locations,dist);
        funcV2 = exeGateCounts.first;
        funcV3 = exeGateCounts.second;
    }
    else{
        funcV3=-1;
    }
    
    return {funcV1,funcV2,funcV3};
}


tuple<int,int,int> costOfInsertedSwap_L0_Nexe_Frontgate_3_3_log(vector<list<gate> > & lines, pair<int,int> & swapGate,const vector<vector<gate> > & layers, 
                                            vector<int> qubits,vector<int> locations,  int ** dist)
{
    bool hasExecutable = false;
    bool isNegative = false;
    int funcV1=0,funcV2=0,funcV3=0;
    for(size_t i=0;i<layers.size();i++)
    {
        int effOfLayer = 0;
        for(size_t j=0;j<layers[i].size();j++)
        {
            gate g1 = layers[i][j];
            // logout<<"begin effOfSWAP"<<endl;
            int effectOnG1 = effOfSWAP(swapGate, g1, locations, dist);
            if(i==0) //gate of layer 0
            {
                int loc1_g1 = locations[g1.control];
                int loc2_g1 = locations[g1.target];
                int dist_g1 = dist[loc1_g1][loc2_g1];
                if(dist_g1 == 2 && effectOnG1==1)
                    hasExecutable = true;
            }
            effOfLayer += effectOnG1;
        }
        funcV1 += effOfLayer;
    }
    // logout<<"funcV1: "<<funcV1<<endl;
    if(hasExecutable)
    {
        has_not_zero=true;
        updateMapping(swapGate,qubits,locations);
        pair<int,int> exeGateCounts = numOfExecutedAndFrontGate_3_3_log(lines,locations,dist);
        funcV2 = exeGateCounts.first;
        funcV3 = exeGateCounts.second;
    }
    else{
        funcV3=-1;
    }

    for(auto layer:layers){
        for(auto g:layer){
            if(exe_st.find(g.weight)!=exe_st.end()){
                int effectOnG1 = effOfSWAP(swapGate, g, locations, dist);
                if(effectOnG1!=-1&&effectOnG1!=0){
                    logout<<"[ERROR]: effectOnG1!=-1 && !=0"<<endl;
                    logout<<"effectOnG1: "<<effectOnG1<<endl;
                    exit(-1);
                }
                funcV1+=effectOnG1;
            }
        }
    }
    exe_st.clear();
    return {funcV1,funcV2,funcV3};
}


/**
*Insert a SWAP gate into the resultant circuit
*
*/
void insertSwap(pair<int,int> & swapGate, vector<gate> & finalList)
{
    int p1 = swapGate.first;
    int p2 = swapGate.second;
    int num = -1;
    finalList.push_back(gate(p1,p2,num));
    finalList.push_back(gate(p2,p1,num));
    finalList.push_back(gate(p1,p2,num));
}

/**
*Insert a bridge cnot gate into the resultant circuit
*/
void inserBridgeCnot(int c, int t, vector<int> & linkPath, vector<gate> & finalList,int number)
{
    int trans_number=-number-2;//0->-2,1->-3,2->-4... ,标注bridge代替的门的编号
    finalList.push_back(gate(linkPath[0],linkPath[1],trans_number));
    finalList.push_back(gate(linkPath[1],linkPath[2],trans_number));
    finalList.push_back(gate(linkPath[0],linkPath[1],trans_number));
    finalList.push_back(gate(linkPath[1],linkPath[2],trans_number));
}



auto getValuedGateList_L0_Nexe_FrontGate(const set<weightedGate > & candiSwap,
                            vector<list<gate> > & lines,
                            int maxLookSize,  const vector<int> & qubits,
                            const vector<int> & locations,  int ** dist)
{
    maxLookSize = MAX_LOOK_SIZE;
    vector<weightedGate > bestGates;
    weightedGate  bestGate;
    vector<vector<gate> > layers = getFrontLayers(lines, maxLookSize);

    using value_gate = pair<weightedGate, tuple<int,int,int> >;
    vector<value_gate>candidate_gate;
    for(set<weightedGate>::iterator it = candiSwap.begin(); it!=candiSwap.end(); it++)
    {
        weightedGate gate1 = * it;
        // logout<<"[TEST] candidate gate: "<<gate1.type<<' '<<gate1.q1<<' '<<gate1.q2<<endl;
        auto effOfGate1 = make_tuple(0,0,0);
        if(gate1.type == "SWAP")
        {
            pair<int,int> equivSwap = make_pair(gate1.q1,gate1.q2);
            // effOfGate1 = costOfInsertedSwapSplit(lines,equivSwap,qubits,locations,dist);
            effOfGate1 = costOfInsertedSwap_L0_Nexe_Frontgate(lines,equivSwap,layers,qubits,locations,dist);

        }
        else if(gate1.type == "CNOT")
        {
            pair<int,int> equivCnot = make_pair(gate1.q1,gate1.q2);
            // effOfGate1 = costOfInsertedCnotSecondOnly(lines,equivCnot,qubits,locations,dist);
            effOfGate1 = costOfInsertedCnot_L0_Nexe_Frontgate(lines,equivCnot,qubits,locations,dist);
        }
        // logout<<"value: "<<effOfGate1.first<<' '<<effOfGate1.second<<endl;
        candidate_gate.push_back({gate1,effOfGate1});
    }
    // logout<<"deal over"<<endl;
    return candidate_gate;
}



// get the valued gate list
// CRG version
auto getValuedGateList_L0_Nexe_FrontGate_3_3_log(const set<weightedGate > & candiSwap,
                            vector<list<gate> > & lines,
                            int maxLookSize,  const vector<int> & qubits,
                            const vector<int> & locations,  int ** dist)
{
    maxLookSize = MAX_LOOK_SIZE;
    vector<weightedGate > bestGates;
    weightedGate  bestGate;
    vector<vector<gate> > layers = getFrontLayers(lines, maxLookSize);

    using value_gate = pair<weightedGate, tuple<int,int,int> >;
    vector<value_gate>candidate_gate;
    for(set<weightedGate>::iterator it = candiSwap.begin(); it!=candiSwap.end(); it++)
    {
        weightedGate gate1 = * it;
        auto effOfGate1 = make_tuple(0,0,0);
        if(gate1.type == "SWAP")
        {
            pair<int,int> equivSwap = make_pair(gate1.q1,gate1.q2);
            effOfGate1 = costOfInsertedSwap_L0_Nexe_Frontgate_3_3_log(lines,equivSwap,layers,qubits,locations,dist);
            

        }
        else if(gate1.type == "CNOT")
        {
            pair<int,int> equivCnot = make_pair(gate1.q1,gate1.q2);
            effOfGate1 = costOfInsertedCnot_L0_Nexe_Frontgate(lines,equivCnot,qubits,locations,dist);
        }
        candidate_gate.push_back({gate1,effOfGate1});
    }
    return candidate_gate;
}

// get the valued gate list
auto getValuedGateList(const set<weightedGate > & candiSwap,
                            vector<list<gate> > & lines,
                            int maxLookSize,  const vector<int> & qubits,
                            const vector<int> & locations,  int ** dist)
{
    vector<weightedGate > bestGates;
    weightedGate  bestGate;
    vector<vector<gate> > layers = getFrontLayers(lines, maxLookSize);
    // logout<<"layer over"<<endl;
    using value_gate = pair<weightedGate, pair<int,int> >;
    vector<value_gate>candidate_gate;
    for(set<weightedGate>::iterator it = candiSwap.begin(); it!=candiSwap.end(); it++)
    {
        weightedGate gate1 = * it;
        // logout<<"[TEST] candidate gate: "<<gate1.type<<' '<<gate1.q1<<' '<<gate1.q2<<endl;
        auto effOfGate1 = make_pair(0,0);
        if(gate1.type == "SWAP")
        {
            pair<int,int> equivSwap = make_pair(gate1.q1,gate1.q2);
            // effOfGate1 = costOfInsertedSwapSplit(lines,equivSwap,qubits,locations,dist);
            effOfGate1 = costOfInsertedSwapSplit(lines,equivSwap,layers,qubits,locations,dist);

        }
        else if(gate1.type == "CNOT")
        {
            pair<int,int> equivCnot = make_pair(gate1.q1,gate1.q2);
            // effOfGate1 = costOfInsertedCnotSecondOnly(lines,equivCnot,qubits,locations,dist);
            effOfGate1 = costOfInsertedCnotSplit(lines,equivCnot,qubits,locations,dist);
        }
        // logout<<"value: "<<effOfGate1.first<<' '<<effOfGate1.second<<endl;
        candidate_gate.push_back({gate1,effOfGate1});
    }
    // logout<<"deal over"<<endl;
    return candidate_gate;
}
/**
* The mapping procedure by using swap and bridge cnot gate
* return the number of ancillary cnot gates
*/

void insert_SWAP_to_phy(weightedGate tar,vector<gate> & phyCir,int &nCnots,vector<int> & qubits,vector<int> & locations){
    pair<int,int> swapg = make_pair(tar.q1,tar.q2);
    insertSwap(swapg, phyCir);
    nCnots += 3;
    updateMapping(swapg, qubits, locations );
}
void insert_CNOT_to_phy(weightedGate tar,vector<gate> & phyCir,int &nCnots,vector<int> & qubits,vector<int> & locations,
                        CouplingGraph * p_cg,vector<list<gate> >  &lines,vector<gate> &activeGates,set<int> & activeLines){
    vector<int> linkedPath = p_cg->getPath(tar.q1,tar.q2);
    
    nCnots +=3;
    int c = qubits[tar.q1];
    int t = qubits[tar.q2];
    gate g1 = lines[c].front();
    lines[c].pop_front();
    gate g2 = lines[t].front();
    lines[t].pop_front();
    if(!(g1 == g2)){
        logout<<"[BUG] error in insert_CNOT_to_phy"<<endl;
        exit(-1);
    }

    for(vector<gate>::iterator it = activeGates.begin(); it!=activeGates.end();it++ )
    {
        gate g2 = *it;
        if(c == g2.control && t == g2.target)
        {
            activeGates.erase(it);
            break;
        }
    }
    activeLines.insert(c);
    activeLines.insert(t);
    inserBridgeCnot(tar.q1,tar.q2,linkedPath, phyCir,g1.weight);
}

// insert CNOT gate, but also record the line situation in re_line, used for later recovery
void insert_CNOT_to_phy_recover(weightedGate tar,vector<gate> & phyCir,int &nCnots,vector<int> & qubits,vector<int> & locations,
                        CouplingGraph * p_cg,vector<list<gate> >  &lines,vector<gate> &activeGates,set<int> & activeLines,
                        vector<list<gate> >  &re_line
                        ){
    
    vector<int> linkedPath = p_cg->getPath(tar.q1,tar.q2);
    
    nCnots +=3;
    int c = qubits[tar.q1];
    int t = qubits[tar.q2];
    gate g1 = lines[c].front();
    lines[c].pop_front();
    re_line[c].push_back(g1);
    gate g2 = lines[t].front();
    lines[t].pop_front();
    re_line[t].push_back(g2);
    if(!(g1 == g2)){
        logout<<"[BUG] error in insert_CNOT_to_phy"<<endl;
        exit(-1);
    }

    for(vector<gate>::iterator it = activeGates.begin(); it!=activeGates.end();it++ )
    {
        gate g2 = *it;
        if(c == g2.control && t == g2.target)
        {
            activeGates.erase(it);
            break;
        }
    }
    activeLines.insert(c);
    activeLines.insert(t);
    inserBridgeCnot(tar.q1,tar.q2,linkedPath, phyCir,g1.weight);
}

int decision_count = 0;

// recover the state after do first insert operation
void doRecoverWork(int &nCnots,int &beforeCnots,
                   set<int> &activeLines,
                   vector<gate> &activeGates,vector<gate> &tmpActiveGates,
                   vector<list<gate> >  &lines,vector<list<gate>> &re_line,
                   int beforePhyCirSize,vector<gate> & phyCir,
                   weightedGate &gateInfo,
                   vector<int> &qubits,vector<int> &locations){
    // recover nCnots
    nCnots=beforeCnots;
    // recover activeLines(set to empty)
    activeLines.clear();
    // recover activeGates
    activeGates=tmpActiveGates;

    // recover lines
    for(int i=0;i<lines.size();i++){
        if(re_line[i].size()){
            auto now=re_line[i].rbegin();
            while(now!=re_line[i].rend()){
                lines[i].push_front(*now);
                now++;
            }
        }
    }

    //recover phyCircle
    while(phyCir.size()>beforePhyCirSize){
        phyCir.pop_back();
    }
    // if SWAP gate is inserted, recover qubits and locations
    if(gateInfo.type == "SWAP"){
        pair<int,int> swap={gateInfo.q1,gateInfo.q2};
        updateMapping(swap,qubits,locations);
    }
}


//final version
bool doTwoInsertOperate_new_func_fixed_abalation( vector<list<gate> >  &lines,
                        vector<gate> &activeGates,set<int> & activeLines,
                        CouplingGraph * p_cg,vector<int> & qubits,vector<int> & locations,
                        int ** dist,vector<gate> & phyCir,int &nCnots){

    //(1) generate and check active gate
    findActivegates(lines,activeGates,phyCir,activeLines, locations,dist);
    if(activeGates.empty()) // all gate have been mapped
        return true;

    //(2) generate candidate swap and bridge gates
    set<weightedGate > candiGates;

    genRelevantGates(p_cg, activeGates, candiGates, locations, dist );

    using value_gate = pair<weightedGate,tuple<int,int,int> >;

    vector<value_gate> gateList=getValuedGateList_L0_Nexe_FrontGate_3_3_log(candiGates,lines,4,qubits,locations,dist);

    pair<weightedGate,weightedGate> bestGate;
    int bestCost = numeric_limits<int>::min();vector<pair<weightedGate,weightedGate> > bestGateList;
    for(auto curGate:gateList){

        auto gateInfo=curGate.first;
        auto gateValue=curGate.second;
        pair<weightedGate,weightedGate> gatePair;
        gatePair.first=gateInfo;
        
        gatePair.first.weight=get<0>(gateValue)+get<1>(gateValue);
        // not consider the third part
        // if(get<2>(gateValue)!=-1)
        //     gatePair.first.weight+=get<2>(gateValue)-((int)activeGates.size());
        if(gatePair.first.weight<=0)  continue; 

        int beforeCnots=nCnots;
        int beforePhyCirSize=phyCir.size();
        vector<list<gate> >re_line(lines.size());

        auto tmpActiveGates=activeGates;

        if(gateInfo.type == "SWAP"){
            insert_SWAP_to_phy(gateInfo,phyCir,nCnots,qubits,locations);
        }
        else{
            insert_CNOT_to_phy_recover(gateInfo,phyCir,nCnots,qubits,locations,p_cg,lines,activeGates,activeLines,re_line);
        }
        checkActivegates_recover(lines,activeGates, phyCir,activeLines, locations, dist,re_line);
            

        findActivegates_recover(lines,activeGates,phyCir,activeLines, locations,dist,re_line);
        if(activeGates.empty()){
            bestCost=gatePair.first.weight;
            gatePair.second=weightedGate();
            gatePair.second.type="NULL";
            bestGateList.clear();
            bestGateList.push_back(gatePair);

            doRecoverWork(nCnots,beforeCnots,activeLines,activeGates,tmpActiveGates,lines,re_line,beforePhyCirSize,phyCir,gateInfo,qubits,locations);
            break;
        }
        set<weightedGate > nextCandiGates;
        genRelevantGates(p_cg, activeGates, nextCandiGates, locations, dist);
        
        auto nextGateList=getValuedGateList_L0_Nexe_FrontGate_3_3_log(nextCandiGates,lines,4,qubits,locations,dist);
        for(auto nextGate:nextGateList){
            auto nextGateInfo=nextGate.first;
            auto nextGateValue=nextGate.second;
            gatePair.second=nextGateInfo;
            

            gatePair.second.weight=get<0>(nextGateValue)+get<1>(nextGateValue);
            // if(get<2>(nextGateValue)!=-1)
            //     gatePair.second.weight+=get<2>(nextGateValue)-((int)activeGates.size());
            
            if(gatePair.second.weight<=0)  continue; 

            if(gatePair.first.weight+gatePair.second.weight>bestCost){
                bestCost=gatePair.first.weight+gatePair.second.weight;
                bestGateList.clear();
                bestGateList.push_back(gatePair);
            }
            else if(gatePair.first.weight+gatePair.second.weight==bestCost){
                bestGateList.push_back(gatePair);
            }
            
        }
        

        doRecoverWork(nCnots,beforeCnots,activeLines,activeGates,tmpActiveGates,lines,re_line,beforePhyCirSize,phyCir,gateInfo,qubits,locations);
    }
   
    if(bestGateList.size()==0){

    }
    else{
        bestGate=bestGateList[rand()%bestGateList.size()];    
    }
        
    if(bestCost<=0){
        int minLenOfFrontCnot = numeric_limits<int>::max();
        int controlOfMinLenGate;
        int targetOfMinlenGate;
        for(size_t i = 0;i < activeGates.size();i++){
            gate frontGate1 = activeGates[i];
            int locOfC = locations[frontGate1.control];
            int locOfT = locations[frontGate1.target];
            if(dist[locOfC][locOfT] < minLenOfFrontCnot)
            {
                minLenOfFrontCnot = dist[locOfC][locOfT];
                controlOfMinLenGate = locOfC;
                targetOfMinlenGate = locOfT;
                
            }
        }
        while(minLenOfFrontCnot!=1){
            vector<int> shortestPath = p_cg->getPath(controlOfMinLenGate,targetOfMinlenGate);
            if(shortestPath.size()<2){
                logout<<"[BUG] shortestPath.size()<2"<<endl;
                exit(-1);
            }

            pair<int,int> swapg1;
            if(global_random_device()%2==0){
                swapg1 = make_pair(shortestPath[0],shortestPath[1]);
                controlOfMinLenGate=shortestPath[1];
            }
                
            else{
                swapg1 = make_pair(shortestPath[shortestPath.size()-1],shortestPath[shortestPath.size()-2]);
                targetOfMinlenGate=shortestPath[shortestPath.size()-2];
            }
            insertSwap(swapg1, phyCir);
            nCnots += 3;
            updateMapping(swapg1, qubits, locations);
            minLenOfFrontCnot--;
        }
    }
    else{
        if(bestGate.first.type == "SWAP"){
            insert_SWAP_to_phy(bestGate.first,phyCir,nCnots,qubits,locations);    
        }
        else{
            insert_CNOT_to_phy(bestGate.first,phyCir,nCnots,qubits,locations,p_cg,lines,activeGates,activeLines);
        }
        checkActivegates(lines,activeGates, phyCir,activeLines, locations, dist);
        findActivegates(lines,activeGates,phyCir,activeLines, locations,dist);

        if(activeGates.empty()) 
            return true;

        if(bestGate.second.type == "SWAP"){
            insert_SWAP_to_phy(bestGate.second,phyCir,nCnots,qubits,locations);    
        }
        else if(bestGate.second.type == "CNOT"){
            insert_CNOT_to_phy(bestGate.second,phyCir,nCnots,qubits,locations,p_cg,lines,activeGates,activeLines);
        }
        else if(bestGate.second.type == "NULL"){

        }
    }
    //(5) check active gate whether some can be inserted into physical circuit
    checkActivegates(lines,activeGates, phyCir,activeLines, locations, dist);
    // logout<<"-----------------"<<endl<<endl;
    return false;
}
int swapAndCnotBasedMapping_advance(int nqubits, vector<gate> & gateList, CouplingGraph * p_cg, 
                                    vector<gate> & phyCir, int ** dist, vector<int> & qubits, vector<int> & locations)
{

    int nCnots = 0;
    set<int> activeLines;
    for(int i=0 ; i<nqubits; i++)
        activeLines.insert(i);
    vector<list<gate> >  lines;
    lines.resize(nqubits);
    genLines(lines, gateList);

    vector<gate> activeGates;
    while(true)
    {
        decision_count++;
        
        //final version
        if(doTwoInsertOperate_new_func_fixed_abalation(lines,activeGates,activeLines,p_cg,qubits,locations,dist,phyCir,nCnots))
            break;
    }
    
    return nCnots;
}



// check whether the physical circuit is valid
void answer_checker(int minCnotCount,vector<gate> &gateList,
                    vector<gate> &phyCirOfBest,const CouplingGraph & cg,
                    int nqubits,vector<int>&qubits,vector<int>&locations){
    for(int i = 0;i < qubits.size(); i++)
    {
        if(qubits[i] == -1)
            continue;
        locations[qubits[i]] = i;
    }

    if(gateList[0].weight>gateList.back().weight){
        reverse(gateList.begin(),gateList.end());
    }
    if(minCnotCount + gateList.size() != phyCirOfBest.size()){
        logout<<"[ERROR] minCnotCount + gateList.size() != phyCirOfBest.size()"<<endl;
        logout<<"minCnotCount: "<<minCnotCount<<endl;
        logout<<"gateList.size(): "<<gateList.size()<<endl;
        logout<<"phyCirOfBest.size(): "<<phyCirOfBest.size()<<endl;
        logout<<"qubitsArrOfBest: "<<endl;
        for (auto i : phyCirOfBest)
        {
            logout << i.control << " " << i.target << endl;
        }
        exit(-1);
    }
    for(auto i:phyCirOfBest)
    {
        // logout<<i.weight<<' '<<i.control<<' '<<i.target<<endl;
        if(i.control == i.target){
            logout<<"[ERROR] control == target"<<endl;
            exit(-1);
        }
        if(i.control<0||i.target<0){
            logout<<"[ERROR] control or target < 0"<<endl;
            exit(-1);
        }
        if(i.control>=cg.getPositions()||i.target>=cg.getPositions()){
            logout<<"[ERROR] control or target >= nlocations"<<endl;
            exit(-1);
        }
        if(cg.dist[i.control][i.target] != 1){
            logout<<"[ERROR] dist != 1"<<endl;
            exit(-1);
        }
    }
    // for(auto i:phyCirOfBest){
    //     logout<<i.weight<<endl;
    // }

    vector<pair<int,int>> graph(gateList.size(),make_pair(-1,-1));
    vector<int>last_gate(nqubits,-1);
    unordered_map<int,bool> is_checked;
    is_checked[-1] = true;
    for(auto i:gateList){
        int target = i.target;
        int control = i.control;
        int number = i.weight;
        graph[number].first = last_gate[target];
        graph[number].second = last_gate[control];
        last_gate[target] = number;
        last_gate[control] = number;
        // logout<<"[INFO] "<<number<<"<-"<<graph[number].first<<"<-"<<graph[number].second<<endl;
    }




    for(int i=0;i<phyCirOfBest.size();i++){
        auto phy_tar = phyCirOfBest[i].target;
        auto phy_con = phyCirOfBest[i].control;
        auto number = phyCirOfBest[i].weight;
        // logout<<"[INFO] "<<number<<" "<<phy_tar<<" "<<phy_con<<endl;
        if(number==-1){
            if(i+2>=phyCirOfBest.size()){
                logout<<"[ERROR] i+2>=phyCirOfBest.size()"<<endl;
                exit(-1);
            }
            swap(qubits[phy_con],qubits[phy_tar]);
            for(int i = 0;i < qubits.size(); i++)
            {
                if(qubits[i] == -1)
                    continue;
                locations[qubits[i]] = i;
            }
            i+=2;
            continue;
        }
        bool is_bridge = false;
        if(number < 0){
            is_bridge=true;
            if(i+3>=phyCirOfBest.size()){
                logout<<"[ERROR] i+3>=phyCirOfBest.size()"<<endl;
                exit(-1);
            }
            number += 2;
            number = -number;
            if(cg.dist[locations[gateList[number].control]][locations[gateList[number].target]] != 2){
                logout<<"[ERROR] dist != 2"<<endl;
                logout<<"number: "<<number<<endl;
                logout<<"control: "<<gateList[number].control<<" target: "<<gateList[number].target<<endl;
                logout<<"location of control: "<<locations[gateList[number].control]<<" location of target: "<<locations[gateList[number].target]<<endl;
                logout<<"cur mapping: "<<endl;
                debugMapping(qubits);
                logout<<"------------"<<endl;


                exit(-1);
            }
            i+=3;
            // continue;
        }
        if(is_checked[number]){
            logout<<"[ERROR] repeated gate in number: "<<number<<endl;
            logout<<"control: "<<phy_con<<" target: "<<phy_tar<<endl;
            exit(-1);
        }
        if(!is_checked[graph[number].first]){
            logout<<"[ERROR] cross error in number: "<<number;
            logout<<"control: "<<phy_con<<" target: "<<phy_tar<<endl;
            logout<<" its constrians is "<<graph[number].first<<" and "<<graph[number].first<<endl;
            logout<<"not checked: "<<graph[number].first<<endl;
            exit(-1);
        }
        if(!is_checked[graph[number].second]){
            logout<<"[ERROR] cross error in number: "<<number;
            logout<<"control: "<<phy_con<<" target: "<<phy_tar<<endl;
            logout<<" its constrians is "<<graph[number].first<<" and "<<graph[number].first<<endl;
            logout<<"not checked: "<<graph[number].first<<endl;
            exit(-1);
        }
        is_checked[number] = true;
        if(is_bridge){
            continue;
        }
        auto origin_gate = gateList[number];
        if(qubits[phy_con] != origin_gate.control){
            logout<<"[ERROR] control error in number: "<<number;
            logout<<"phy control: "<<phy_con<<" phy target: "<<phy_tar<<endl;
            logout<<"logi control: "<<origin_gate.control<<" logi target: "<<origin_gate.target<<endl;
            logout<<"qubits[phy_con]: "<<qubits[phy_con]<<endl;
            exit(-1);
        }
        if(qubits[phy_tar] != origin_gate.target){
            logout<<"[ERROR] target error in number: "<<number;
            logout<<"phy control: "<<phy_con<<" phy target: "<<phy_tar<<endl;
            logout<<"logi control: "<<origin_gate.control<<" logi target: "<<origin_gate.target<<endl;
            logout<<"qubits[phy_tar]: "<<qubits[phy_tar]<<endl;
            exit(-1);
        }
        
        
        
        // logout<<"[INFO] over"<<endl;
    }
    logout<<"[INFO] check over. Circult has no fault."<<endl;
}




// Begin of algorithm
void EffectiveQMFramework(string ciucuitFileName, int outLimit, int inLimit=IN_LIMIT)
{

    // logout<<"[INFO] greedy version"<<endl;
    logout<<"[INFO] parameter: IN_LIMIT: "<<IN_LIMIT<<" MAX_LOOK_SIZE: "<<MAX_LOOK_SIZE<<" C_MAX: "<<C_MAX<<endl;
    vector<gate> gateList;  //logical circuit
    vector<gate> phyCir;   //physical circuit

    vector<list<gate> >  lines;
    vector<gate> activeGates;
    set<int> activeLines;
    CouplingGraph cg;
    int ** dist;

    // Read qasm file into gate list
    ciucuitFileName = "./"+ciucuitFileName;
    QASMReader reader(std::move(ciucuitFileName));
    reader.parse(gateList);
    int nqubits = reader.getVerNum();//number of logical
    
    logout << "Seed: "<<global_algorithm_seed<<endl;
    logout << "The number of logic qubits:" << nqubits << endl;
    logout << "The number of CNOT gates:" << gateList.size() << endl;
    
    //generate coupling map , distance matrix and interaction matrix

    // choose a device 
    // cg.buildQX20();
    // cg.buildQ16();
    // cg.build_IBMQ_GUADALUPE();
    // cg.build_nXn_grid(19);
    // cg.build_Google_Sycamore();
    // cg.build_IBM_Rochester();
    // cg.buildQX5();
    // cg.build_IBM_Johannesburg();
    // cg.build_IBM_Penguin_v1();
    // cg.build_IBM_Penguin_v2();

    if(selected_device_enum == "TOKYO"){
        cg.buildQX20();
    }
    else if(selected_device_enum == "Q16"){
        cg.buildQ16();
    }
    else if(selected_device_enum == "GUADALUPE"){
        cg.build_IBMQ_GUADALUPE();
    }
    else if(selected_device_enum == "ROCHESTER"){
        cg.build_IBM_Rochester();
    }
    else if(selected_device_enum == "19X19"){
        cg.build_nXn_grid(19);
    }
    else if(selected_device_enum == "SYCAMORE"){
        cg.build_Google_Sycamore();
    }
    else{
        logout<<"[ERROR] no such device"<<endl;
        exit(-1);
    }


    if(cg.positions<nqubits){
        logout<<"[INFO] no enough qubits"<<endl;
        exit(0);
    }
    
    int nlocations = cg.getPositions();
    cg.computeDist();
    dist = cg.getDist();
    interactMat = initInteractionMatrix(nqubits);
    getInteractionMatrix(gateList,interactMat);

    //Initialize qubits[] and locations[]
    vector<int> locations(nqubits,-1);
    vector<int> qubits(nlocations,-1);
    int minCnotCount = numeric_limits<int>::max();

    vector<int> qubitsArrOfBest; 
    int dirOfBest = 0;
    vector<gate> phyCirOfBest; 
    int noBetterCounter = 0; 

    randInitMapping(qubits,locations,nlocations,nqubits,cg);
    int iter_count = 0;

    for(int _=0;;_++)
    {


        PerurbationByShuffle(qubits,locations,nlocations,nqubits,cg);
       
        noBetterCounter++;
        bool isNeedNoGate = false;
        int inner_count = 0;
        int epoch_min_cnot_count=numeric_limits<int>::max();

        if(gateList[0].weight>gateList.back().weight){
            reverse(gateList.begin(),gateList.end());
        }

        for(int i=0;;i++)
        {
            
            iter_count++;
            inner_count++;
            vector<int> tmpQubitsArr = qubits;
            int curCnotCount = swapAndCnotBasedMapping_advance(nqubits, gateList, &cg, phyCir, dist, qubits, locations);
            if(curCnotCount < minCnotCount)
            {
                noBetterCounter = 0;
                minCnotCount = curCnotCount;
                phyCirOfBest = phyCir;
                current_time = chrono::system_clock::now();
                auto time_cost=((double)chrono::duration_cast<chrono::milliseconds>(current_time - start_time).count())/1000;
                logout<<"[INFO] better ans: "<<minCnotCount<<" in time: "<<time_cost<<"s."<<endl;
                if(i%2 == 0) //forward traversal
                {
                    dirOfBest = 0;
                    qubitsArrOfBest = tmpQubitsArr;
                }
                else         //reverse traversal
                {
                    dirOfBest = 1;
                    qubitsArrOfBest = qubits;
                }
            }

            reverse(gateList.begin(),gateList.end());
            phyCir.clear();
            if(minCnotCount == 0)
            {
                isNeedNoGate = true;
                break;
            }
            if(epoch_min_cnot_count>curCnotCount){
                epoch_min_cnot_count=curCnotCount;
                inner_count = 0;
            }
            if(inner_count >= inLimit)
            {
                break;
            }
        }
        
        current_time = chrono::system_clock::now();
        auto time_cost=((double)chrono::duration_cast<chrono::milliseconds>(current_time - start_time).count())/1000;
        if(time_cost > cutoff_time){
            break;
        }
        if(noBetterCounter == C_MAX || isNeedNoGate){
            break;
        }
            
    }

    if(dirOfBest == 1){
        reverse(phyCirOfBest.begin(),phyCirOfBest.end());
    }
    answer_checker(minCnotCount,gateList,phyCirOfBest,cg,nqubits,qubitsArrOfBest,locations);

    logout << "Cnot gate count: " << minCnotCount << endl;
    logout << "Size of the best physical circuit:" << phyCirOfBest.size() << endl;
    logout << "Direction of best mapping: " << dirOfBest << endl;//最优方向
    logout << "[INFO] restart count: "<<iter_count<<endl;
    //clear job
    clearInteractionMatrix(interactMat,nqubits);
}

int doExperiments(vector<string>files,int seed,string log_path,int epoch=1){
    

    cout<<"#File:"<<files.size()<<endl;
    auto filenameSplit=[](string &file)->pair<string,string>{
        int index;
        for(int i=file.size()-1;i;i--){
            if(file[i]=='\\'||file[i]=='/'){
                index=i;
                break;
            }
        }
        string path=file.substr(0,index);
        string name=file.substr(index+1,file.size()-index-1);
        return {path,name};
    };

    using namespace chrono;
    for(auto &file:files){
        global_algorithm_seed=seed;
        global_random_device.seed(seed);
        srand(seed);
        auto p=filenameSplit(file);
        decision_count=0;
        cout<<"Benchmark circuit " << p.second << " is in mapping...."<<endl;
        string logfile = "./" +log_path;
        logout.open(logfile,ios::app);
        for(int i=0;i<epoch;i++){
            start_time = system_clock::now();

            #ifdef GREEDY
            EffectiveQMFramework(file,1000);
            #else
                // do nothing
            #endif

            current_time = system_clock::now();

            auto epoch_time=((double)duration_cast<milliseconds>(current_time - start_time).count())/1000;
            logout <<"Time cost :" << epoch_time << " seconds." << endl;
        }
        logout<<endl;
        logout<<endl<<endl;
        
        logout.close();
        cout<<"Benchmark circuit " << p.second << " has been mapped."<<endl;
    }
    return 0;
}


// commmand params:
// 1. seed
// 2. log path
// 3. device enum
// 4. IN_LIMIT
// 5. MAX_LOOK_SIZE
// 6. C_MAX
// 7 ... n. file path
int main(int argc, char* argv[])
{
    if(argc<7){
        cerr<<"no enough args."<<endl;
        exit(-1);
    }


    vector<string>files;
    int seed=atoi(argv[1]);
    string log_path=argv[2];
    selected_device_enum=argv[3];
    IN_LIMIT=atoi(argv[4]);
    MAX_LOOK_SIZE=atoi(argv[5])-1;
    C_MAX=atoi(argv[6]);

    for(int i=7;i<argc;i++){
        files.push_back(argv[i]);
    }
    doExperiments(files,seed,log_path);
    
    return 0;
}
