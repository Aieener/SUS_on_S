#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "square.h"
#include "cells.h"
#include "hardrods.h"
#include "histogram.h"
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <vector>
#include <array>
using namespace std;

#ifndef MC_H
#define MC_H

class MC
{
    private:
    	//data members;
        std::vector<HR> VRodlist; // the list storage the Vertical Rods;
        std::vector<HR> HRodlist; // the list storage the Horizantal Rods;
        std::vector<HR> URodlist; // the list storage the Up Rods;
        int n0,n1,n2; // size of the box
        int length; // length of the rod
        long int step;
        double z; 
        double nh,nv,nu,dh,dv,du,ah,av,au;
        

    public:
        enum INIT {EMPTY,BOX,PLANE};
        MC(long int ST, int LEN, int N0, int N1, int N2, double Z); //(x->n0, y->n1,z->n2)

    	// ********* Getters********//
        const vector<HR>& getVRodlist() const;
        const vector<HR>& getHRodlist() const;
        const vector<HR>& getURodlist() const;

    	// ******** Other Functianality *******//
        // void Add(Cells &s,array<double,100000> &WF,double &prob,double &size,double &S,double &deltaS);
        // void Del(Cells &s,array<double,100000> &WF,double &prob,double &size,double &S,double &deltaS);
        std::tuple<int,int,int,int> Add(Cells &s,double &prob,double &proba);
        void Add_redo(Cells &s,std::tuple<int,int,int,int> Ad);

        std::tuple<int,int,int,int> Del(Cells &s,double &prob,double &probd, double &size);
        void Del_redo(Cells &s,std::tuple<int,int,int,int> De);

    	
        // array<double,2000> MCSUS(); 
        void MCSUS();
        // std::tuple<array<double,2000>,vector<vector<double>>> MCSUS();
};

#endif /* MC_H */