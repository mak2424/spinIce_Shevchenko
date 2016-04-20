#include <iostream>
#include "squarespinicearray.h"
#include "clustermachine.h"
#include "StateMachineFree.h"

using namespace std;

int main()
{
    SquareSpinIceArray sys;
    sys.dropSpinIce(20,20);
    //sys.groundState();
    cout << "\nNum of spins = " << sys.size() << endl;
    
    
    ClusterMachine clusters(&sys,0.8);
    
    cout << "\nclusters.maximalSize = " << clusters.maximalSize() << endl;
    
    //параметр порядка
    double eta = clusters.maximalSize() / (double)sys.size();
    cout << "\neta = " << eta << endl;
    
    
    
    
    sys.save("1.mfsys");
    
    cout << "Finish him!" << endl;
    return 0;
}

