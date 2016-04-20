#include <float.h>

#include "squarespinicearray.h"
#include "clustermachine.h"

using namespace std;

unsigned start_time;
//#define MEASURE(x,y) cout<<y<<": "; start_time =  clock(); x; cout<< setprecision(10)<<(clock()-start_time)/1000.<<endl;

int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;
    
    config::Instance()->set2D(); //размерность системы
    SquareSpinIceArray sys; //объявляем систему типа Square Spin Ice
    sys.dropSpinIce(70,70);
    sys.setInteractionRange(0.8); // установить расстояние до соседей
    
    
    sys.state = sys.groundState(); //приводим к Ground State
    ClusterMachine clusters(&sys,0.8);
    
    double initialSystemEnergy = sys.E();
    
    int nIntervals = 90000;//число интервалов энергий
    double stepOfInterval = (2 * abs(initialSystemEnergy)) / nIntervals;
    
    cout << "Number of spins = " << sys.size() << endl;
    cout << "Initial System Energy (GS) = " << initialSystemEnergy << endl;
    
    vector <unsigned long long> arrIntervals(nIntervals, 0);
    vector <double> arrIntervalsEnergy(nIntervals, 0.);
    double eTotal = 0;
    unsigned long long iterations = 0;
    bool chaotic = true;
    double M_Ei = 0;
    double M_E2 = 0;
    double M2_E = 0;
    double previous_M_Ei = 0;
    double sumMx = 0;
    double sumMx2 = 0;
    double sumMaxCluster = 0;
    
    double sysEnergy = initialSystemEnergy;
    
    unsigned long long nIters = 10000; //число Монте-Карло шагов
    
    double T=20; //температура
    
    ///-------------------------------------------------------------------------------------------------
    start_time =  clock();
    while(chaotic){
        //** Time = 0.004 для 10000 шагов Number of spins = 9940
        double oldE = sysEnergy;
        int num = sys.state.randomize();
        double newE = sys.E();
        
        if (newE > oldE){
            double probability = exp((oldE-newE)/T);
            if (probability<Random::Instance()->nextDouble()){
                sys[num]->rotate(true);
                sysEnergy = oldE;
            }
            else sysEnergy = newE;
        }
        else sysEnergy = newE;
        //*/
        
        eTotal+=sysEnergy;
        
        ///sumMaxCluster += clusters.maximalSize(); ///времязатратный!!!
        
        ///double Mx = sys.M().x; ///времязатратный!!! 0.499 для 10000 шагов Number of spins = 9940
        ///sumMx += Mx;
        ///sumMx2 += Mx*Mx;
        
        iterations++;
        
        //**
        //распределяем энергию по интервалам Time = 0.005 для 10000 шагов Number of spins = 9940
        for (int i = 0; i < nIntervals; i++) {
            if (i == 0 && sysEnergy < initialSystemEnergy) {
                arrIntervals[i]++;
                arrIntervalsEnergy[i] += sysEnergy;
                break;
            }
            else if ((sysEnergy >= (initialSystemEnergy + (i * stepOfInterval))) && ((sysEnergy < (initialSystemEnergy + ((i + 1) * stepOfInterval))))) {
                arrIntervals[i]++;
                arrIntervalsEnergy[i] += sysEnergy;
                break;
            } else if ((i == (nIntervals - 1)) &&(sysEnergy > (initialSystemEnergy + (i * stepOfInterval)))) {
                arrIntervals[i]++;
                arrIntervalsEnergy[i] += sysEnergy;
                break;
            }
        }
        //*/
        
        //** Time = 0.098 из 0.118 для 10000 шагов Number of spins = 9940
        previous_M_Ei = M_Ei;
        M_Ei = 0;
        M_E2 = 0;
        
        for (int i = 0; i < nIntervals; i++)
            if (arrIntervals[i] > 0) {
                double E_i = arrIntervalsEnergy[i] / arrIntervals[i];
                double P_E = (double)arrIntervals[i] / iterations;
                if (P_E>1) cout << "WARNING!!! p = " << P_E << endl;
                M_Ei += E_i * P_E;
                M_E2 += E_i*E_i * P_E;
            }
        //*/
        
        double condition_A = abs(M_Ei / previous_M_Ei);
        double condition_B = abs(previous_M_Ei / M_Ei);
        
        ///if (iterations >= nIters)
        if (iterations >= nIters && (((condition_A >= 0.999999) && (condition_A <= 1.000001)) || ((condition_B >= 0.999999) && (condition_B <= 1.000001))))
            chaotic = false;
    }
    
    cout << "\niterations = " << iterations << endl;
    
    cout << "\nTime(while) = " << (clock()-start_time)/1000. << " sec." << endl;
    ///-------------------------------------------------------------------------------------------------------
    
    M2_E = M_Ei*M_Ei;
    
    double C = ((M_E2-M2_E)/(T*T)) / (double)sys.size();
    
    cout << "Runtime: " << clock() / 1000. << " sec." << endl;
    cout << "\nFinish him!" << endl << endl;
    return 0;
}
