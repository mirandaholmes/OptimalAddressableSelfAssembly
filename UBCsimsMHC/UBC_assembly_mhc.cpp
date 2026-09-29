/* 

Created Summer 2024 by Sushrut Tadwalkar.

Tests assembly into n copies of "UBC" for spanning tree interactions.
*/




#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>             // for timing operations
#include <sstream>

#include "Demo.h"
#include "VMMC.h"
#include "StickySquare.h"

using namespace std;

//typedef std::chrono::high_resolution_clock Clock;

#ifndef M_PI
    #define M_PI 3.1415926535897932384626433832795
#endif


int num(int i, int j, int l0) {
    return l0*j +i;
}


int main(int argc, char** argv)
{
    // Process input
    if(argc <= 3) {    
        cout << "Not enough input arguments! \nUsage: <program name> <ds> <dw> <rep> " << endl;
        return 1;
    }

    // parameters to set in command line
    double ds;
    double dw;
    int rep;

    // extract arguments from command line
    stringstream convert1 {argv[1]};
    stringstream convert2 {argv[2]};
    stringstream convert3 {argv[3]};
    convert1 >> ds;
    convert2 >> dw;
    convert3 >> rep;


    cout << "Checking inputs: " << endl;
    cout << "ds = " << ds << ", dw = " << dw << ", rep = " << rep << endl;


    /* ----------  Parameters  ---------- */

        double dens = 0.05;

        int n0 = 213;   // number of particles in base cell
        int l0 = sqrt(n0);    // side length of the cell (should be sqrt(n0))
        int f0 = n0;   // for reporting; what is completed fragment size?

        string filehead = "Data/ubc_mhc_" + to_string(rep);

        // Simulation parameters.
        int nsteps = 128; //256;    // how many steps of the simulation to run
        int nmc = 2e4; //1e5;       // how many MC sweeps per step (save data after each step)
        int nCopies = 4;     // how many copies of each square there are
        double sval = ds;//20.0; // strong side interaction
        double wval = dw;//2.0; // weak side interaction


        // derived parameters
        double boxLength = round(sqrt(n0*nCopies / dens));  // length of simulation box (must be an integer)
        int nParticles = n0*nCopies;       // number of particles



        /* ----------  Create output files and data structures  ---------- */


        // create filenames
        string statfile, trajfile;
        statfile = filehead + "_stats.txt";
        trajfile = filehead + "_traj.txt";

        // create string describing simulation parameters
        string description;
        stringstream os;
        os << "n0=" << n0 << " nCopies=" << nCopies << " nParticles=" << nParticles << " dens=" << dens
        << " nsteps=" << nsteps << " nmc=" << nmc 
        << " sval=" << sval << " wval=" << wval;
        description = os.str();

        // Display experiment info
        std::cout << "-----------\n  " << description << endl;
        std::cout << statfile << ", " << trajfile << endl;


        vector<double> stats;       // statistics to compute
        vector<int> fragmenthist;   // histogram of fragment sizes
        int nfrag;                  // number of fragments


        // Parameters that generally shouldn't be changed
        bool isLattice = true;              // whether particles must stay on a lattice
        unsigned int dimension = 2;         // dimension of simulation box
        double interactionRange = 1.1;      // interaction range (used in CellList; not used in StickySquare)
        unsigned int maxInteractions = 6;   // maximum number of interactions per particle (orig 15)
        double interactionEnergy = 0;       // (not used; needed to set up StickySphere) pair interaction energy scale (in units of kBT)

        // Initialise random number generator
        MersenneTwister rng;



        /* ----------  Make interactions  ---------- */

        vector<Triple> north0;
        vector<Triple> east0;
        vector<Triple> north;
        vector<Triple> east;
        int num1,num2;  // indices of interaction

        //auto num{ [l0](int i, int j)->int { return (l0*j +i); } };  // calculates square's actual number, from col/row

        // x-axis interactions
        east0 = {{0,1,sval},{1,2,sval},{2,3,sval},{4,5,sval},{5,6,sval},{6,7,sval},{8,9,sval},{10,11,sval},{12,13,sval},{14,15,sval},{16,17,sval},{18,19,sval},{20,21,sval},{22,23,sval},{24,25,sval},{26,27,sval},{28,29,sval},{30,31,sval},{32,33,sval},{34,35,sval},{36,37,sval},{38,39,sval},{40,41,sval},{42,43,sval},{44,45,sval},{46,47,sval},{48,49,sval},{49,50,sval},{51,52,sval},{52,53,sval},{54,55,sval},{55,56,sval},{60,61,sval},{61,62,sval},{64,65,sval},{65,66,sval},{66,67,sval},{67,68,sval},{72,73,sval},{73,74,sval},{74,75,sval},{75,76,sval},{76,77,sval},{77,78,sval},{78,79,sval},{80,81,sval},{81,82,sval},{82,83,sval},{83,84,sval},{84,85,sval},{85,86,sval},{87,88,sval},{89,90,sval},{90,91,sval},{92,93,sval},{94,95,sval},{96,97,sval},{98,99,sval},{100,101,sval},{102,103,sval},{104,105,sval},{106,107,sval},{108,109,sval},{112,113,sval},{113,114,sval},{114,115,sval},{115,116,sval},{116,117,sval},{117,118,sval},{118,119,sval},{119,120,sval},{121,122,sval},{123,124,sval},{125,126,sval},{127,128,sval},{129,130,sval},{131,132,sval},{133,134,sval},{135,136,sval},{137,138,sval},{139,140,sval},{141,142,sval},{143,144,sval},{145,146,sval},{146,147,sval},{147,148,sval},{148,149,sval},{149,150,sval},{150,151,sval},{151,152,sval},{152,153,sval},{79,145,sval},{153,154,sval},{154,155,sval},{155,156,sval},{156,157,sval},{157,158,sval},{158,159,sval},{159,160,sval},{160,161,sval},{161,162,sval},{162,163,sval},{163,207,sval},{164,165,sval},{165,166,sval},{166,167,sval},{167,168,sval},{170,171,sval},{171,172,sval},{173,174,sval},{174,175,sval},{176,177,sval},{183,184,sval},{185,186,sval},{187,188,sval},{189,190,sval},{191,192,sval},{193,194,sval},{195,196,sval},{199,200,sval},{202,203,sval},{203,204,sval},{205,206,sval},{207,208,sval},{208,209,sval},{209,210,sval},{210,211,sval},{211,212,sval},{56,57,wval},{57,58,wval},{58,59,wval},{59,60,wval},{63,64,wval},{68,69,wval},{70,71,wval},{71,72,wval},{110,111,wval},{178,179,wval},{180,181,wval},{197,198,wval}};    

        // y-axis interactions
        north0 = {{9,2,sval},{13,9,sval},{17,13,sval},{21,17,sval},{25,21,sval},{29,25,sval},{33,29,sval},{37,33,sval},{41,37,sval},{45,41,sval},{49,45,sval},{54,49,sval},{63,55,sval},{64,56,sval},{65,57,sval},{70,65,sval},{66,58,sval},{71,66,sval},{67,59,sval},{72,67,sval},{68,60,sval},{69,61,sval},{62,52,sval},{52,46,sval},{46,42,sval},{42,38,sval},{38,34,sval},{34,30,sval},{30,26,sval},{26,22,sval},{22,18,sval},{18,14,sval},{14,10,sval},{10,5,sval},{88,81,sval},{93,88,sval},{97,93,sval},{101,97,sval},{105,101,sval},{109,105,sval},{113,109,sval},{122,113,sval},{126,122,sval},{130,126,sval},{134,130,sval},{138,134,sval},{142,138,sval},{146,142,sval},{89,86,sval},{94,91,sval},{102,98,sval},{106,102,sval},{119,110,sval},{120,111,sval},{111,106,sval},{123,120,sval},{127,124,sval},{131,127,sval},{139,135,sval},{144,139,sval},{153,143,sval},{172,164,sval},{173,168,sval},{175,169,sval},{177,170,sval},{178,174,sval},{179,175,sval},{180,176,sval},{181,177,sval},{182,179,sval},{184,180,sval},{186,184,sval},{188,186,sval},{190,188,sval},{192,190,sval},{194,192,sval},{196,194,sval},{197,196,sval},{199,197,sval},{200,198,sval},{202,200,sval},{206,201,sval},{207,204,sval},{212,205,sval},{8,1,wval},{12,8,wval},{16,12,wval},{20,16,wval},{24,20,wval},{28,24,wval},{32,28,wval},{36,32,wval},{40,36,wval},{44,40,wval},{48,44,wval},{61,51,wval},{53,47,wval},{47,43,wval},{43,39,wval},{39,35,wval},{35,31,wval},{31,27,wval},{27,23,wval},{23,19,wval},{19,15,wval},{15,11,wval},{11,6,wval},{98,94,wval},{135,131,wval},{145,141,wval},{141,137,wval},{137,133,wval},{133,129,wval},{129,125,wval},{125,121,wval},{121,112,wval},{112,108,wval},{108,104,wval},{104,100,wval},{100,96,wval},{96,92,wval},{92,87,wval},{87,80,wval},{195,193,wval},{193,191,wval},{191,189,wval},{189,187,wval},{187,185,wval},{185,183,wval}};
        

        // Done 1st square. Make several copies of each square
        for(int k=0; k<north0.size(); k++) {
            for(int c1=0; c1<nCopies; c1++) {
                for(int c2=0; c2<nCopies; c2++) {
                    north.push_back({north0[k].i+c1*n0, north0[k].j+c2*n0, north0[k].val});            
                }   
            }
        }
        for(int k=0; k<east0.size(); k++) {
            for(int c1=0; c1<nCopies; c1++) {
                for(int c2=0; c2<nCopies; c2++) {
                    east.push_back({east0[k].i+c1*n0, east0[k].j+c2*n0, east0[k].val});
                }
            }   
        }
        

        Interactions interactions(nParticles,north,east);


        /* ----------  Initialise data structures & classes  ---------- */

        // Data structures.
        std::vector<Particle> particles(nParticles);    // particle container
        bool isIsotropic[nParticles];                   // whether the potential of each particle is isotropic
        
        // Create simulation box object.
        std::vector<double> boxSize {boxLength,boxLength};          // simulation box sizes
        Box box(boxSize,isLattice);  // Initialise simulation box object.

        // Initialise cell list.
        CellList cells; 
        cells.setDimension(dimension);
        cells.initialise(box.boxSize, interactionRange);

        // Initialise the sticky square potential model.
        StickySquare StickySquare(box, particles, cells,
            maxInteractions, interactionEnergy, interactionRange,
            interactions);



        /* ----------  Initialise Particles  ---------- */
        // Generate a random particle configuration using Initialise object & MersenneTwister object
        Initialise initialise;
        initialise.random(particles, cells, box, rng, false, isLattice);

        // Initialise data structures needed by the VMMC class.
        double coordinates[dimension*nParticles];
        double orientations[dimension*nParticles];

        // Copy particle coordinates and orientations into C-style arrays.
        for (int i=0;i<nParticles;i++) {
            for (int j=0;j<dimension;j++) {
                coordinates[dimension*i + j] = particles[i].position[j];
                orientations[dimension*i + j] = particles[i].orientation[j];
            }
            // Set all particles as isotropic.
            isIsotropic[i] = true;
        }

        /* ----------  Initialise VMMC functions & object  ---------- */

        // Initialise the VMMC callback functions.
        using namespace std::placeholders;
        vmmc::CallbackFunctions callbacks;
        callbacks.energyCallback =
        std::bind(&StickySquare::computeEnergy, StickySquare, _1, _2, _3);
        callbacks.pairEnergyCallback =
        std::bind(&StickySquare::computePairEnergy, StickySquare, _1, _2, _3, _4, _5, _6);
        callbacks.interactionsCallback =
        std::bind(&StickySquare::computeInteractions, StickySquare, _1, _2, _3, _4);
        callbacks.postMoveCallback =
        std::bind(&StickySquare::applyPostMoveUpdates, StickySquare, _1, _2, _3);


        // Variables to intialise VMMC object; these shouldn't change
        double maxTrialTranslation = 1.5;
        double maxTrialRotation = 0.0;
        double probTranslate = 1.0;
        double referenceRadius = 0.5;
        bool isRepulsive = false;

        // Initialise VMMC object. 
        vmmc::VMMC vmmc(nParticles, dimension, coordinates, orientations,
            maxTrialTranslation, maxTrialRotation, probTranslate, referenceRadius, 
            maxInteractions, &boxSize[0], isIsotropic, isRepulsive, callbacks, isLattice);



        /* ----------  Create output file  ---------- */

        // Create output file & log initial condition
        InputOutput io;
        io.appendXyzTrajectory(dimension, particles, box, true, n0, description,trajfile);  // "true" is for the first line
        
        // Initalise statistics and write to file
        stats = {0,StickySquare.getEnergy()*nParticles};
        nfrag = StickySquare.computeFragmentHistogram(n0,fragmenthist);
        stats.insert(stats.end(), fragmenthist.begin(), fragmenthist.end());
        // write to file (erase old contents)
        io.appendStats(stats,true,description,statfile);   



        /* ----------  Run the simulation!  ---------- */
        clock_t start_time = clock();  // time the loop
        for (int i=0;i<nsteps;i++)
        {
            // Increment simulation by nmc Monte Carlo Sweeps.
            vmmc += nmc*nParticles; 

            // Append particle coordinates to an xyz trajectory.
            io.appendXyzTrajectory(dimension, particles, false,trajfile);

            // Compute statistics
            stats = {(double)(i+1.0), StickySquare.getEnergy()*nParticles};
            nfrag = StickySquare.computeFragmentHistogram(n0,fragmenthist);
            stats.insert(stats.end(), fragmenthist.begin(), fragmenthist.end());
            io.appendStats(stats,false,"",statfile);
            std::cout << "  i = " << i << ", ncomplete = " << stats.at(f0-1+2) << endl;
        }
        // save time
        double time = (clock() - start_time ) / (double) CLOCKS_PER_SEC;

        /* ----------  Report stuff  ---------- */
        double efinal = StickySquare.getEnergy()*nParticles;

        std::cout << "  Time = " << time << " seconds" << endl;
        std::cout << "  Number of fully completed fragments = " << stats.at(f0-1+2) << endl;
        std::cout << "    Fragments: ";
        for (int x : fragmenthist)  std::cout << x << " "; std::cout << endl;
    
    std::cout << "Complete!";
    // We're done!
    return (EXIT_SUCCESS);
}
