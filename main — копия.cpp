#include <ctime>
#include <sstream>
#include <cmath>
#include <queue>
#include <fstream>
#include <string>
#include <float.h>

#include "config.h"
#include "StateMachine.h"
#include "PartArray.h"
#include "honeycombspinicearray.h"
#include "squarelattice.h"
#include "squarespinicearray.h"
#include "clustermachine.h"
#include "StateMachineFree.h"
#include "random.h"

//#define DBL_EPSILON 2.2204460492503131e-16

#include <mpi.h>


using namespace std;

//Запись информации о системе в файл
void sampleInfo(string type, int latticeSize, unsigned long long numSpins, double interactionRange, unsigned long long nb, double initialSystemEnergy, unsigned long long nIters, int itersC);
void gnuplotScript(int latticeSize, unsigned long long numSpins);
void taskKillScript();
void combineFilesScript(int latticeSize);

int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;
    
    //mpi
    int rank, size;
    MPI_Status status;
    (void)status;
    
    MPI_Init (&argc, &argv);	/* starts MPI */
    
    MPI_Comm_size (MPI_COMM_WORLD, &size);	/* get number of processes */
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
    
    //printf ("current process id %d \n", rank);
    
    MPI_Barrier(MPI_COMM_WORLD); //запуск процессов одновременно
    //
    
    config::Instance()->set2D(); //размерность системы
    
    SquareSpinIceArray sys; //объявляем систему типа Square Spin Ice
    
    //sys.setInteractionRange(0.8); // установить расстояние до соседей
    //sys.interactionRange(); // длина взаимодействия
    
    int latticeSize = 7; //размер решетки latticeSize x latticeSize
    //sys.dropSpinIce(1,1,2); //параметр решетки 2
    sys.dropSpinIce(latticeSize,latticeSize); //набрасываем спины на решетку размером (latticeSize,latticeSize)
    
    sys.state = sys.groundState(); //приводим к Ground State
    //sys.state.randomize(100); //случайные направления у спинов ("количество разворотов")
    
    ClusterMachine clusters(&sys,0.8);
    
    double initialSystemEnergy = sys.E();
    
    int nIntervals = 1000;//число интервалов энергий
    double stepOfInterval = (2 * abs(initialSystemEnergy)) / nIntervals;
    unsigned long long nIters = 2000000; //размер окна (число Монте-Карло шагов / itersC)
    int itersC = 5;
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
    
    //ofstream f("_capacity.dat");
    
    if (rank !=0){
    //for (T=1; T<20; T+=0.5){
    //for (T=0.0001; T<20; T+=1){
    //for (T=20; T>0; T--){
        
        //printf ("current process id %d \n", rank);
        
        stringstream ss;
        ss << rank;
        string str;
        ss >> str;
        
        T = (double)rank/10.0;
        
        //cout << endl << T << endl;
        queue<double> cArray;
        double cTotal = 0;
        
        vector<double> arrIntervalsCapacity;
        vector<double> C_arrIntervals;
        
        unsigned long long counter = 0;
        
        for (int l=0; l < itersC; l++) {
            //vector <unsigned long long> arrIntervals;
            //vector <double> arrIntervalsEnergy;
            vector <unsigned long long> arrIntervals(nIntervals, 0);
            vector <double> arrIntervalsEnergy(nIntervals, 0.);
            double eTotal = 0;
            unsigned long long iterations = 0;
            bool chaotic = true;
            //queue<double> eArray;
            
            double M_Ei = 0;
            double M_E2 = 0;
            double M2_E = 0;
            
            double previous_M_Ei = 0;
            
            double sumMaxCluster = 0;
            
            ///double sumModuleMx = 0;
            ///double sumModuleMy = 0;
            
            double sumMx = 0;
            ///double sumMy = 0;
            double sumMx2 = 0;
            ///double sumMy2 = 0;
            
            ///double sumMlength = 0;
            ///double sumMlength2 = 0;
            
            int o = 0;
            while(chaotic){
                o++;
                double oldE = sysEnergy;
                int num = sys.state.randomize();
                double newE = sys.E();
                if (newE > oldE){
                    double probability = exp((oldE-newE)/T);
                    if (probability<Random::Instance()->nextDouble()){
                        sys[num]->rotate();
                    }
                }
                
                sysEnergy = sys.E();
                eTotal+=sysEnergy;
                
                //sys.M().x; //sum(mx)
                //sys.M().y; //sum(my)
                //sys.M().length(); //|M|=sqrt(mx^2+my^2)
                
                //eArray.push(sys.E());
                
                //MxArray.push(sys.M().x);
                //MyArray.push(sys.M().y);
                //MlengthArray.push(sys.M().length());
                
                sumMaxCluster += clusters.maximalSize();
                
                ///sumModuleMx += abs(sys.M().x);
                ///sumModuleMy += abs(sys.M().y);
                
                //sumModuleMx = 0;
                //sumModuleMy = 0;
                /*for (unsigned long long l=0; l < sys.size(); l++){
                    sumModuleMx += sys.parts[l]->m.x;
                    sumModuleMy += sys.parts[l]->m.y;
                }*/
                
                //cout << "sumModuleMx = " << sumModuleMx << " sys.M().x = " << sys.M().x << endl;
                //cout << "sumModuleMy = " << sumModuleMy << " sys.M().y = " << sys.M().y << endl;
                
                sumMx += sys.M().x;
                ///sumMy += sys.M().y;
                sumMx2 += sys.M().x*sys.M().x;
                ///sumMy2 += sys.M().y*sys.M().y;
                
                ///sumMlength += sys.M().length();
                ///sumMlength2 += sys.M().length()*sys.M().length();
                
                iterations++;
                counter++;
                
                //распределяем энергию по интервалам
                for (int i = 0; i < nIntervals; i++) {
                    if (i == 0 && sysEnergy < initialSystemEnergy) {
                        arrIntervals[i]++;
                        arrIntervalsEnergy[i] += sysEnergy;
                    }
                    else if ((sysEnergy >= (initialSystemEnergy + (i * stepOfInterval))) && ((sysEnergy < (initialSystemEnergy + ((i + 1) * stepOfInterval))))) {
                        arrIntervals[i]++;
                        arrIntervalsEnergy[i] += sysEnergy;
                    } else if ((i == (nIntervals - 1)) &&(sysEnergy > (initialSystemEnergy + (i * stepOfInterval)))) {
                        arrIntervals[i]++;
                        arrIntervalsEnergy[i] += sysEnergy;
                    }
                }
                
                if (o==1){
                    previous_M_Ei = M_Ei;
                    M_Ei = 0;
                    M_E2 = 0;
                    
                    for (int i = 0; i < nIntervals; i++) //{ //arrIntervalsEnergy.size(); i++)
                        if (arrIntervals[i] > 0) {
                            //sigma += (eArray.front()-eAver)*(eArray.front()-eAver);
                            //double E_i = (initialSystemEnergy + ((i + 0.5) * stepOfInterval)); //вычисляем середину интервала
                            double E_i = arrIntervalsEnergy[i] / arrIntervals[i];
                            double P_E = (double)arrIntervals[i] / iterations;
                            if (P_E>1) cout << "WARNING!!! p = " << P_E << endl;
                            M_Ei += E_i * P_E;
                            M_E2 += E_i*E_i * P_E;
                            //eArray.pop();
                        }
                    
                    //cout << "\nM_Ei = " << M_Ei << "\nprevious_M_Ei = " << previous_M_Ei << endl;
                    double condition_A = abs(M_Ei / previous_M_Ei);
                    double condition_B = abs(previous_M_Ei / M_Ei);
                    
                    //if (iterations >= 10000)
                    if (iterations >= nIters && (((condition_A >= 0.999999) && (condition_A <= 1.000001)) || ((condition_B >= 0.999999) && (condition_B <= 1.000001))))
                        //if (iterations >= 10000 && (fabs(M_Ei - previous_M_Ei) <= DBL_EPSILON * fmax(fabs(M_Ei), fabs(previous_M_Ei))))
                        chaotic = false;
                    o = 0;
                }
                if (rank == 1)
                    if (counter%100000 == 0){
                        cout << endl << counter << " of " << nIters*itersC << " steps\n";
                        double programTime = clock();
                        programTime /= 1000.0;
                        cout << "Time: " << programTime << " sec." << endl;
                    }
            }
            
            //cout << "iterations = " << iterations;
            
            sumMaxCluster /= iterations;
            ///sumModuleMx /= iterations;
            ///sumModuleMy /= iterations;
            
            sumMx /= iterations;
            ///sumMy /= iterations;
            sumMx2 /= iterations;
            ///sumMy2 /= iterations;
            
            ///sumMlength /= iterations;
            ///sumMlength2 /= iterations;
            
            ///cout << endl << "Mx = " << sumMx << endl << "My = " << sumMy << endl;
            
            //энергия на одну частицу
            string address5 = "__EnergyParticle.dat4";
            address5.insert ( strlen("_"), str);
            ofstream fEnergyPart(address5.c_str());
            fEnergyPart << M_Ei / (double)sys.size() << "\t" << T << endl;
            fEnergyPart.close();
            
            //энергия системы
            string address6 = "__EnergySystem.dat5";
            address6.insert ( strlen("_"), str);
            ofstream fEnergySys(address6.c_str());
            fEnergySys << M_Ei << "\t" << T << endl;
            fEnergySys.close();
            
            //магнитная восприимчивость по оси X
            string address2 = "__Mx.dat1";
            address2.insert ( strlen("_"), str);
            ofstream fMx(address2.c_str());
            double hi_Mx = ((sumMx2 - sumMx*sumMx)/T) / (double)sys.size();
            fMx << hi_Mx <<"\t" << T << endl;
            fMx.close();
            
            ///магнитная восприимчивость по оси Y
            /*
            string address3 = "__My.dat2";
            address3.insert ( strlen("_"), str);
            ofstream fMy(address3.c_str());
            double hi_My = ((sumMy2 - sumMy*sumMy)/T) / (double)sys.size();
            fMy << hi_My <<"\t" << T << endl;
            fMy.close();
            */
            
            ///магнитная восприимчивость по диагонали (модуль магнитной восприимчивости)
            /*
            string address4 = "__Mmodule.dat3";
            address4.insert ( strlen("_"), str);
            ofstream fMmodule(address4.c_str());
            double hi_Mlength = ((sumMlength2 - sumMlength*sumMlength)/T) / (double)sys.size();
            fMmodule << hi_Mlength <<"\t" << T << endl;
            fMmodule.close();
            */
            
            ///параметр порядка
            string address7 = "__Eta.dat6";
            address7.insert ( strlen("_"), str);
            ofstream fEta(address7.c_str());
            //double eta = clusters.maximalSize() / (double)sys.size();
            double eta = sumMaxCluster / (double)sys.size();
            //cout << "\neta = " << eta << endl;
            fEta << eta <<"\t" << T << endl;
            fEta.close();
            
            ///шахматная намагниченность (параметр порядка)
            /*
            string address8 = "__StaggeredEta.dat7";
            address8.insert ( strlen("_"), str);
            ofstream fStaggeredEta(address8.c_str());
            double staggeredEta = (sumModuleMx + sumModuleMy) / 2.;
            fStaggeredEta << staggeredEta <<"\t" << T << endl;
            fStaggeredEta.close();
            */
            
            //cout << "\niterations = " << iterations << endl;
            //cout << "\nM_Ei = " << M_Ei << "\nprevious_M_Ei = " << previous_M_Ei << endl;
            //system("pause");
            
            //double eAver = eTotal/iterations; //remake average
            //double sigma=0;
            
            /*
             * while(!eArray.empty()){
             * M_Ei += eArray.front() * P_E;
             * M_E2 += eArray.front()*eArray.front() * P_E;
             * eArray.pop();
             * }
             */
            
            //sigma /= (double) iterations;
            
            M2_E = M_Ei*M_Ei;
            
            //double C = (sigma/(T*T)) /(double)sys.size();
            double C = ((M_E2-M2_E)/(T*T)) / (double)sys.size();
            
            //распределяем теплоемкость по интервалам
            if (arrIntervalsCapacity.size() == 0) {
                arrIntervalsCapacity.push_back(C);
                C_arrIntervals.push_back(1);
            }
            else {
                for (unsigned long long i = 0; i < arrIntervalsCapacity.size(); i++){
                    if (fabs(C - arrIntervalsCapacity[i]) <= DBL_EPSILON * fmax(fabs(C), fabs(arrIntervalsCapacity[i]))){
                        C_arrIntervals[i]++;
                        break;
                    }
                    else {
                        arrIntervalsCapacity.push_back(C);
                        C_arrIntervals.push_back(1);
                        break;
                    }
                }
            }
            
            cArray.push(C);
            cTotal += C;
            //f << C <<"\t" << T << endl;
            //cout << "\nC = " << C << endl;
            //sys.save("1.mfsys");
            
            string address = "__capacity.dat";
            address.insert ( strlen("_"), str);
            ofstream fCapacity(address.c_str());
            fCapacity << C <<"\t" << T << endl;
            fCapacity.close();
        }
        
        double M_Ci = 0;
        //ищем наиболее вероятную теплоемкость
        for (unsigned long long i = 0; i < arrIntervalsCapacity.size(); i++){
            double P_C = (double)C_arrIntervals[i] / itersC;
            //cout << "P_C = " << P_C << endl;
            //cout << "arrIntervalsCapacity[i] = " << arrIntervalsCapacity[i] << endl;
            if (P_C>1) cout << "WARNING!!! pС = " << P_C << endl;
            M_Ci += arrIntervalsCapacity[i] * P_C;
            //M_C2 += C_i*C_i * P_C;
        }
        
        //cout << "\nM_Ci = " << M_Ci << endl;
        
        double sigmaC = 0;
        double cAver = cTotal/itersC;
        while(!cArray.empty()){
            sigmaC += (cArray.front()-cAver)*(cArray.front()- cAver);
            cArray.pop();
        }
        
        sigmaC /= (double) itersC;
        sigmaC /= sqrt(sigmaC);
        //f << cAver <<"\t" << T << "\t" << sigmaC << endl;
        string address = "__capacity.dat";
        address.insert ( strlen("_"), str);
        ofstream fCapacity(address.c_str());
        fCapacity << M_Ci <<"\t" << T << endl;
        fCapacity.close();
        //sys.save("1.mfsys");
    }
    sys.save("1.mfsys");
    
    if (rank == 1) {
        double programTime = 0;
        programTime = clock(); //ms
        programTime /= 1000.0; //sec
        cout << "Runtime: " << programTime << " sec." << endl; // время работы программы
        ofstream info("info.txt", ios::app);
        info << "Time: " << programTime << " sec." << endl;
        info.close();
    }
    
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
    gnuplot << "set terminal png font \"Verdana,14\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_Capacity.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"C / (N k)\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_Capacity.txt' using 2:1 notitle, 'heat_" << numSpins << ".txt' using 2:1 with lines notitle" << endl;
    //с барами
    //gnuplot << "plot [*:*] [*:*] '_Capacity.dat' using 2:1:3 with yerrorbars notitle, 'heat_" << numSpins << ".txt' using 2:1 with lines notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 14" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_Capacity.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"C / (N k)\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_Capacity.txt' using 2:1 notitle, 'heat_" << numSpins << ".txt' using 2:1 with lines notitle" << endl;
    //gnuplot << "plot [*:*] [*:*] '_Capacity.txt' using 2:1 notitle
    
    
    //_XMx.txt
    //.png
    gnuplot << "set terminal png font \"Verdana,14\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_XMx.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"X_Mx / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_XMx.txt' using 2:1 notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 14" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_XMx.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"X_Mx / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_XMx.txt' using 2:1 notitle" << endl;
    
    
    //_EnergyParticle.txt
    //.png
    gnuplot << "set terminal png font \"Verdana,14\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_EnergyParticle.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"E / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_EnergyParticle.txt' using 2:1 notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 14" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_EnergyParticle.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"E / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_EnergyParticle.txt' using 2:1 notitle" << endl;
    
    
    //_EnergySystem.txt
    //.png
    gnuplot << "set terminal png font \"Verdana,14\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_EnergySystem.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"E\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_EnergySystem.txt' using 2:1 notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 14" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_EnergySystem.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"E\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_EnergySystem.txt' using 2:1 notitle" << endl;
    
    
    //_Eta.txt
    //.png
    gnuplot << "set terminal png font \"Verdana,14\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_Eta.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"Eta\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [0:1.05] '_Eta.txt' using 2:1 notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 14" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_Eta.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"Eta\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [0:1.05] '_Eta.txt' using 2:1 notitle" << endl;
    
    /*
    //_XMy.txt
    //.png
    gnuplot << "set terminal png font \"Verdana,14\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_XMy.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"X_My / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_XMy.txt' using 2:1 notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 14" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_XMy.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"X_My / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_XMy.txt' using 2:1 notitle" << endl;
    
    
    //_XMmodule.txt
    //.png
    gnuplot << "set terminal png font \"Verdana,14\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_XMmodule.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"X_|M| / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_XMmodule.txt' using 2:1 notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 14" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_XMmodule.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"X_|M| / N\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_XMmodule.txt' using 2:1 notitle" << endl;
    
    
    //_StaggeredEta.txt
    //.png
    gnuplot << "set terminal png font \"Verdana,14\" size 1000, 800" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_StaggeredEta.png\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"Eta Staggered\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_StaggeredEta.txt' using 2:1 notitle" << endl;
    
    //.ps
    gnuplot << "set terminal postscript enhanced lw 2 14" << endl;
    gnuplot << "set output \"" << latticeSize << "x" << latticeSize << "_StaggeredEta.ps\"" << endl;
    gnuplot << "set style data points" << endl;
    gnuplot << "set style fill solid 1.00 border -1" << endl;
    gnuplot << "set ylabel \"Eta Staggered\"" << endl;
    gnuplot << "set xlabel \"Kb T/D\"" << endl;
    gnuplot << "set xtic offset character 1" << endl;
    gnuplot << "plot [*:*] [*:*] '_StaggeredEta.txt' using 2:1 notitle" << endl;
    */
    
    gnuplot << "call \"kill_wgnuplot.exe.bat\"" << endl;
}

//bat-файл для снятия задачи wgnuplot.exe
void taskKillScript() {
    ofstream taskKill("kill_wgnuplot.exe.bat");
    taskKill << "taskkill / F / IM wgnuplot.exe" << endl;
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
