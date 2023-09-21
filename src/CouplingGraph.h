#ifndef COUPLINGGRAPH_H_INCLUDED
#define COUPLINGGRAPH_H_INCLUDED

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>

using namespace std;

#define MAX_DIST 1000

class CouplingGraph{
    private:
        void initDist();

    public:
        vector<vector<int>> adjList;
        vector<int>deg;
        int positions;  //number of physical qubits
        int ** dist; //distance matrix
        int ** path; //path matrix
    public:
        CouplingGraph() {
            positions = 0;
        }
        CouplingGraph(int pos);
        ~CouplingGraph();

        int ** getDist() {return dist;}

        int getPositions() const {return positions;}

        //build coupling graph for IBM Q Tokyo
        void buildQX20();

        void buildQ9();

        void buildQ10();

        void buildQ16();
        
        void buildQX5();

        void build_IBMQ_GUADALUPE();
        void build_IBM_Rochester();
        void build_Google_Sycamore();
        void build_IBM_Johannesburg();
        void build_IBM_Penguin_v1();
        void build_IBM_Penguin_v2();

        void build_nXn_grid(int n);
        //build testing coupling graph of square
        void buildSquare();

        //get the distance matrix with Floyd-Warshall algorithm
        void computeDist();

        vector<int> getPath(int i,int j);

        void showInfo();

};

#endif
