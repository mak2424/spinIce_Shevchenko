#include <float.h>

#include "squarespinicearray.h"
#include "clustermachine.h"

//#include <iomanip>

//
//#include <ctime>
//#include <sstream>
//#include <cmath>
//#include <queue>
//#include <fstream>
//#include <string>

//#include "config.h"
//#include "StateMachine.h"
//#include "PartArray.h"
//#include "honeycombspinicearray.h"
//#include "squarelattice.h"
//#include "StateMachineFree.h"
//#include "random.h"

//#define DBL_EPSILON 2.2204460492503131e-16

#include <mpi.h>


using namespace std;

//Запись информации о системе в файл
void sampleInfo(string type, int latticeSize, unsigned long long numSpins, double interactionRange, unsigned long long nb, double initialSystemEnergy, unsigned long long nIters, int itersC);
void gnuplotScript(int latticeSize, unsigned long long numSpins);
void taskKillScript();
void combineFilesScript(int latticeSize);

unsigned start_time;
//#define MEASURE(x,y) cout<<y<<": "; start_time =  clock(); x; cout<< setprecision(10)<<(clock()-start_time)/1000.<<endl;

int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;
    
    //----mpi
    int rank, size;
    MPI_Status status;
    (void)status;
    MPI_Init (&argc, &argv);	/* starts MPI */
    MPI_Comm_size (MPI_COMM_WORLD, &size);	/* get number of processes */
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
    //printf ("current process id %d \n", rank);
    MPI_Barrier(MPI_COMM_WORLD); //запуск процессов одновременно
    //----
    
    config::Instance()->set2D(); //размерность системы
    SquareSpinIceArray sys; //объявляем систему типа Square Spin Ice
    
    int latticeSize = 70; //размер решетки latticeSize x latticeSize // 22 - 1012, 224 - 100800
    //sys.dropSpinIce(1,1,2); //параметр решетки 2 для SSI 1x1
    sys.dropSpinIce(latticeSize,latticeSize); //набрасываем спины на решетку размером (latticeSize,latticeSize)
    //sys.dropSpinIce(latticeSize,latticeSize,3);
    
    sys.setInteractionRange(0.8); // установить расстояние до соседей
    //sys.interactionRange(); // длина взаимодействия
    
    sys.state = sys.groundState(); //приводим к Ground State
    //sys.state.randomize(100); //случайные направления у спинов ("количество разворотов")
    
    ClusterMachine clusters(&sys,0.8);
    
    double initialSystemEnergy = sys.E();
    
    int nIntervals = 90000;//число интервалов энергий
    double stepOfInterval = (2 * abs(initialSystemEnergy)) / nIntervals;
    unsigned long long nIters = 10000; //число Монте-Карло шагов
    double T=0;
    double sysEnergy = initialSystemEnergy;
    
    if (rank == 0){
        cout << "Number of spins = " << sys.size() << endl;
        cout << "Initial System Energy (GS) = " << initialSystemEnergy << endl;
        //*
        unsigned long long nb = 0;
        for (signed long long i=0; i<sys.size(); i++) {
            if(sys.neighbourSize(i) > nb) nb = sys.neighbourSize(i);
            //cout << sys.neighbourSize(i) << endl; //число соседей у i-го спина
        }
        cout << "Neighbours = " << nb <<endl; //число соседей
        //*/
        
        gnuplotScript(latticeSize, sys.size());
        taskKillScript();
        combineFilesScript(latticeSize);
        sampleInfo(sys.type().toStdString(), latticeSize, sys.size(), sys.interactionRange(), nb, initialSystemEnergy, nIters, itersC);
    }
    
    if (rank !=0){
        stringstream ss;
        ss << rank;
        string str;
        ss >> str;
        
        T = 20;
        ///T = (double)rank/10.0;
        
        unsigned long long counter = 0;
        
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
            
            //sysEnergy = sys.E();
            eTotal+=sysEnergy;
            
            ///sumMaxCluster += clusters.maximalSize();
            
            ///double Mx = sys.M().x; ///времязатратный!!! 0.499 для 10000 шагов Number of spins = 9940
            ///sumMx += Mx;
            ///sumMx2 += Mx*Mx;
            
            iterations++;
            counter++;
            
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
            
            //cout << "\nM_Ei = " << M_Ei << "\nprevious_M_Ei = " << previous_M_Ei << endl;
            double condition_A = abs(M_Ei / previous_M_Ei);
            double condition_B = abs(previous_M_Ei / M_Ei);
            
            //if (iterations >= 10000)
            if (iterations >= nIters && (((condition_A >= 0.999999) && (condition_A <= 1.000001)) || ((condition_B >= 0.999999) && (condition_B <= 1.000001))))
            //if (iterations >= 10000 && (fabs(M_Ei - previous_M_Ei) <= DBL_EPSILON * fmax(fabs(M_Ei), fabs(previous_M_Ei))))
                chaotic = false;
            
            /*
                if (rank == 1)
                    if (counter%10000000 == 0){
                        cout << endl << counter << " of " << nIters << " steps\n";
                        cout << "Runtime: " << clock() / 1000.0 << " sec." << endl;
                    }
            //*/
        }
        
        cout << "Time(while) = " << (clock()-start_time)/1000. << endl;
        ///-------------------------------------------------------------------------------------------------------
        
        
        //cout << "iterations = " << iterations;
        
        /**
        sumMaxCluster /= iterations;
        sumMx /= iterations;
        sumMx2 /= iterations;
        
        //магнитная восприимчивость по оси X
        string address2 = "__Mx.dat1";
        address2.insert ( strlen("_"), str);
        ofstream fMx(address2.c_str());
        double hi_Mx = ((sumMx2 - sumMx*sumMx)/T) / (double)sys.size();
        fMx << hi_Mx <<"\t" << T << endl;
        fMx.close();
        
        ///параметр порядка
        string address7 = "__Eta.dat6";
        address7.insert ( strlen("_"), str);
        ofstream fEta(address7.c_str());
        //double eta = clusters.maximalSize() / (double)sys.size();
        double eta = sumMaxCluster / (double)sys.size();
        //cout << "\neta = " << eta << endl;
        fEta << eta <<"\t" << T << endl;
        fEta.close();
        //*/
        
        M2_E = M_Ei*M_Ei;
        
        //double C = (sigma/(T*T)) /(double)sys.size();
        double C = ((M_E2-M2_E)/(T*T)) / (double)sys.size();
        
        string address = "__capacity.dat";
        address.insert ( strlen("_"), str);
        ofstream fCapacity(address.c_str());
        fCapacity << C <<"\t" << T << endl;
        fCapacity.close();
    }
    sys.save("1.mfsys");
    
    //*
	MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 1) {
        cout << "Runtime: " << clock() / 1000. << " sec." << endl; // время работы программы
        ofstream info("info.txt", ios::app);
        info << "Runtime: " << clock() / 1000. << " sec." << endl;
        info.close();
    }
    //*/
    
    cout << "\n- " << rank << " Finish him!" << endl << endl;
    
    MPI_Finalize();
    return 0;
}


///functions

//Запись информации о системе в файл
void sampleInfo(string type, int latticeSize, unsigned long long numSpins, double interactionRange, unsigned long long nb, double initialSystemEnergy, unsigned long long nIters, int itersC) {
    //(maxC, Tc, maxHi_Mx, Tm)
    ofstream info("info.txt");
    
    info << "Type: " << type << endl;
    info << "Dimension: " << config::Instance()->dimensions() << "D" << endl;
    info << "Size of lattice: " << latticeSize << "x" << latticeSize << endl;
    info << "Number of spins: " << numSpins << endl;
    if (interactionRange == 0) info << "Interaction Range: " << "All-to-All" << endl;
    else info << "Interaction Range: " << interactionRange << endl;
    info << "Neighbours = " << nb << endl;
    //info << "Step of temperature: " << dT << endl;
    info << "Initial System Energy (GS) = " << initialSystemEnergy << endl;
    info << "Initial System Energy (GS) / numSpins = " << initialSystemEnergy / numSpins << endl;
    info << "Monte Carlo Steps: " << nIters*itersC << endl;
    //info << "Max heat capacity = " << maxC << " at T = " << Tc << endl;
    //info << "Max hi_Mx" << maxHi_Mx << " at T = " << Tm << endl;
    //info << "Time: " << clock()/1000.0 << " sec." << endl; // время работы программы
    info.close();
}

//параметры для gnuplot
void gnuplotScript(int latticeSize, unsigned long long numSpins) {
    
    //ofstream gnuplot("_Plot.gp");
    
    stringstream ss;
    ss << latticeSize;
    string str;
    ss >> str;
    
    string name = "_Plot_All.gp";
    name.insert ( strlen(""), str);
    name.insert ( strlen(""), "x");
    name.insert ( strlen(""), str);
    ofstream gnuplot(name.c_str());
    
    
    //_Capacity.txt
    //.png
    gnuplot << "set terminal png font \"Verdana,16\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_Capacity.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"C / (N k)\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_Capacity.txt' using 2:1 notitle, 'heat_" << numSpins << ".txt' using 2:1 with lines notitle" << endl;
    //с барами
    //gnuplot << "plot [*:*] [*:*] '_Capacity.dat' using 2:1:3 with yerrorbars notitle, 'heat_" << numSpins << ".txt' using 2:1 with lines notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 16" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_Capacity.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"C / (N k)\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_Capacity.txt' using 2:1 notitle, 'heat_" << numSpins << ".txt' using 2:1 with lines notitle" << endl;
    //gnuplot << "plot [*:*] [*:*] '_Capacity.txt' using 2:1 notitle
    
    
    /*
    //_XMx.txt
    //.png
    gnuplot << "set terminal png font \"Verdana,16\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_XMx.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"X_Mx / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_XMx.txt' using 2:1 notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 16" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_XMx.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"X_Mx / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_XMx.txt' using 2:1 notitle" << endl;
    
    
    //_EnergyParticle.txt
    //.png
    gnuplot << "set terminal png font \"Verdana,16\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_EnergyParticle.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"E / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_EnergyParticle.txt' using 2:1 notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 16" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_EnergyParticle.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"E / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_EnergyParticle.txt' using 2:1 notitle" << endl;
    
    
    //_EnergySystem.txt
    //.png
    gnuplot << "set terminal png font \"Verdana,16\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_EnergySystem.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"E\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_EnergySystem.txt' using 2:1 notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 16" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_EnergySystem.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"E\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_EnergySystem.txt' using 2:1 notitle" << endl;
    
    
    //_Eta.txt
    //.png
    gnuplot << "set terminal png font \"Verdana,16\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_Eta.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"Eta\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [0:1.05] '_Eta.txt' using 2:1 notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 16" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_Eta.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"Eta\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [0:1.05] '_Eta.txt' using 2:1 notitle" << endl;
    
    
    //_XMy.txt
    //.png
    gnuplot << "set terminal png font \"Verdana,16\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_XMy.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"X_My / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_XMy.txt' using 2:1 notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 16" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_XMy.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"X_My / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_XMy.txt' using 2:1 notitle" << endl;
    
    
    //_XMmodule.txt
    //.png
    gnuplot << "set terminal png font \"Verdana,16\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_XMmodule.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"X_|M| / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_XMmodule.txt' using 2:1 notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 16" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_XMmodule.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"X_|M| / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_XMmodule.txt' using 2:1 notitle" << endl;
    
    
    //_StaggeredEta.txt
    //.png
    gnuplot << "set terminal png font \"Verdana,16\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_StaggeredEta.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"Eta Staggered\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_StaggeredEta.txt' using 2:1 notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 16" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_StaggeredEta.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"Eta Staggered\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "#set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_StaggeredEta.txt' using 2:1 notitle" << endl;
    */
    
    gnuplot << "call \"kill_wgnuplot.exe.bat\"" << endl;
}

//bat-файл для снятия задачи wgnuplot.exe
void taskKillScript() {
    ofstream taskKill("kill_wgnuplot.exe.bat");
    taskKill << "taskkill /f /im wgnuplot.exe" << endl;
}

//bat-файл для объединения файлов
void combineFilesScript(int latticeSize) {
    //ofstream combineFiles("combineFiles.bat");
    ofstream combineFiles("run.bat");
    combineFiles << "@echo off" << endl;
    combineFiles << "copy /b *.dat _Capacity.txt" << endl;
    combineFiles << "copy /b *.dat1 _XMx.txt" << endl;
    combineFiles << "copy /b *.dat2 _XMy.txt" << endl;
    combineFiles << "copy /b *.dat3 _XMmodule.txt" << endl;
    combineFiles << "copy /b *.dat4 _EnergyParticle.txt" << endl;
    combineFiles << "copy /b *.dat5 _EnergySystem.txt" << endl;
    combineFiles << "copy /b *.dat6 _Eta.txt" << endl;
    combineFiles << "copy /b *.dat7 _StaggeredEta.txt" << endl;
    combineFiles << "start "<< latticeSize << "x" << latticeSize << "_Plot_All.gp\"" << endl;
}
