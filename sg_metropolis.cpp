/* 
 * Author: nlare
 * Prog: spin_glass
 * Date: 03.03.15
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include "mersenne.cpp"
#include <omp.h>
#include <boost/filesystem.hpp>

// #define DEBUG

using namespace std;

/*  
    Задается матрица взаимодействия, которая должна быть шестимерной в соответствии с размерностью матрицы спинов
    Параметр порядка Паризи - (s_ij)^2 на нормировочный коэффициент для данной решетки
*/

int ***spin1;
int ***spin2;
// J - обменный интеграл. Здесь задается матрица взаимодействия
int ***J1;
int ***J2;

int L;
// int mult;

int init_spin_matrix();
int init_J();
int neighbour_spins1(int,int,int);
int neighbour_spins2(int,int,int);

double energy();
double pp_parisi();

int write_to_file(string,double);

int main(int argc, char **argv)  {

    double dE, W, T1, T2, r;
    double exp_mcs, max_mcs, max_stat;
    double EE, MM, EE2, PP2, PP, C_v, X;

    int spin_inverted;
    int ci, cj, ck, cl, cm, cn; // random var's
    long seed;

    double time_temp_b, time_temp_e, time_all_b, time_all_e;
    double time_run;

    stringstream ss;

    ofstream f_out_pp, f_out_ee, f_out_x, f_out_cv;
    ofstream script_f;
    ofstream plot_f;

    srand(time(NULL));

    seed = 1 + rand() % 1000;

    L = 4;
    T1 = 0.5;
    T2 = 1.5;
    max_mcs = 10000;
    max_stat = 20;

    boost::filesystem::create_directories("result");
    boost::filesystem::create_directories("graph");
    boost::filesystem::create_directories("temp");

    // spin = new int **[L];
    // for(int i=0;i<N;i++)    {
    //     spin[i] = new int *[L];
    //     for(int j=0;j<N;j++)    {
    //         spin[i][j] = new int [L];
    //     }
    // }

    if(init_spin_matrix()!=0) cerr << "Spins initialization failed!" << endl;
    if(init_J()!=0) cerr << "J-matrix initialization failed!" << endl;

    ss << "result/PP_PARISI_L=" << L;
    f_out_pp.open(ss.str().c_str());
    if(!f_out_pp) cout << "Cannot open \"" << ss.str().c_str() <<"\" Check free space or permission." << endl;
    
    ss.str("");
    ss << "result/ENERGY_L=" << L;
    f_out_ee.open(ss.str().c_str());
    if(!f_out_ee) cout << "Cannot open \"" << ss.str().c_str() <<"\" Check free space or permission." << endl;

    ss.str("");
    ss << "result/C_v_L=" << L;
    f_out_cv.open(ss.str().c_str());
    if(!f_out_cv) cout << "Cannot open \"" << ss.str().c_str() <<"\" Check free space or permission." << endl;

    ss.str("");
    ss << "result/X_L=" << L;
    f_out_x.open(ss.str().c_str());
    if(!f_out_x) cout << "Cannot open \"" << ss.str().c_str() <<"\" Check free space or permission." << endl;

    #ifdef DEBUG
    cout << "In Debug Area" << endl;

    // ci = Mersenne.IRandomX(0,L-1);
    // cj = Mersenne.IRandomX(0,L-1);
    // ck = Mersenne.IRandomX(0,L-1);
    // cl = Mersenne.IRandomX(0,L-1);
    // cm = Mersenne.IRandomX(0,L-1);
    // cn = Mersenne.IRandomX(0,L-1);

    // spin[ci][cj][ck][cl][cm][cn] = 1;
    // cout << "spin[" << ci << "][" << cj << "]["<< ck << "]["<< cl << "][" << cm << "][" << cn << "]=" << spin1[ci][cj][ck][cl][cm][cn] << endl;
    #endif

    // omp_set_num_threads(2);

    time_all_b = omp_get_wtime();

    double T = T1;

    while(T < T2)  {

    if(T > 1.09 && T < 1.20) T+=0.01; else T+=0.1;

    cout << "Begin: T = " << T << "\n";

    time_temp_b = omp_get_wtime();

    for(int stat=0;stat<max_stat;stat++)    {

    for(int mcs=0;mcs<max_mcs;mcs++)    {

    #pragma omp parallel private(ci,cj,ck,r,spin_inverted,dE,W,seed) num_threads(2)
    {
    #pragma omp sections
    {
    #pragma omp section
    {
    seed = time(0) + 1 + random() % 1000;
    CRandomMersenne Mersenne(seed);
    Mersenne.RandomInit(seed);

    for(int i=0;i<L;i++)    {
        for(int j=0;j<L;j++)    {
            for(int k=0;k<L;k++)    {

                ci = Mersenne.IRandom(0,L-1);
                cj = Mersenne.IRandom(0,L-1);
                ck = Mersenne.IRandom(0,L-1);

                r = Mersenne.IRandom(1,10000)/10000.0;

                if(r < 0.5) spin_inverted = spin1[ci][cj][ck];
                else spin_inverted = -spin1[ci][cj][ck];

                dE = (spin1[ci][cj][ck] - spin_inverted)*neighbour_spins1(ci,cj,ck);
                if(dE<0.0) { spin1[ci][cj][ck] = spin_inverted; } 
                else     W=1.0/exp(dE/T);

                if((r = Mersenne.IRandom(1,10000)/10000.0)<=W) spin1[ci][cj][ck] = spin_inverted;

            }
        }
    }

    }

    #pragma omp section
    {
    seed = time(0) + 1 + random() % 1000;
    CRandomMersenne Mersenne(seed);
    Mersenne.RandomInit(seed);

    for(int i=0;i<L;i++)    {
        for(int j=0;j<L;j++)    {
            for(int k=0;k<L;k++)    {

                ci = Mersenne.IRandom(0,L-1);
                cj = Mersenne.IRandom(0,L-1);
                ck = Mersenne.IRandom(0,L-1);

                r = Mersenne.IRandom(1,10000)/10000.0;

                if(r < 0.5) spin_inverted = spin2[ci][cj][ck];
                else spin_inverted = -spin2[ci][cj][ck];
                //cout << T << endl;
                dE = (spin2[ci][cj][ck] - spin_inverted)*neighbour_spins2(ci,cj,ck);
                // cout << dE << endl;
                if(dE<0.0) spin2[ci][cj][ck] = spin_inverted; 
                else     W=1.0/exp(dE/T);

                if((r = Mersenne.IRandom(1,10000)/10000.0<=W)) spin2[ci][cj][ck] = spin_inverted;

            }
        }
    }
    
    }

    }

    }

    if(mcs > 3600)   {
        PP += pp_parisi();
        PP2 += PP*PP;
        EE += energy();
        EE2 += EE*EE;
    }
    // PP += pp_parisi();
    // cout << fixed << setprecision(6) << PP << "\t" << T << endl;

    }

    }

    time_temp_e = omp_get_wtime();

    // cout << "Time_temp: " << time_temp_e - time_temp_b << endl;

    f_out_pp << fixed << setprecision(6) << PP/(max_mcs*max_stat) << "\t" << T << endl;
    f_out_ee << fixed << setprecision(6) << EE/(max_mcs*max_stat) << "\t" << T << endl;

    C_v = ((EE2/(max_mcs*max_stat)) - pow((EE/(max_mcs*max_stat)),2.0))/(T*T);
    X = ((PP2/(max_mcs*max_stat)) - pow((PP/(max_mcs*max_stat)),2.0))/T;

    f_out_cv << fixed << setprecision(6) << C_v << "\t" << T << endl;
    f_out_x << fixed << setprecision(6) << X << "\t" << T << endl;

    std::cout << std::fixed << std::setprecision(5)\
              << "End: T = " << T <<  ", EE2 = " << EE2/max_mcs/max_stat << ", PP2 = " << PP2/max_mcs/max_stat\
              << ", EE = " << EE/max_mcs/max_stat << ", PP = " << PP/max_mcs/max_stat\
              << ", time = " << time_temp_e - time_temp_b << std::endl; 

    PP = 0.0;
    EE = 0.0;
    PP2 = 0.0;
    EE2 = 0.0;
    C_v = 0.0;
    X = 0.0;

    }

    time_all_e = omp_get_wtime();

    time_run = time_all_e - time_all_b;

    cout << "Time_run: " << time_run << "'s" << endl;

    ss.str("");
    ss << "temp/" << "/PP_PARISI_L=" << L << ".plot";
    plot_f.open(ss.str().c_str());

    plot_f << "#!/usr/bin/gnuplot -persist\n" << \
                 "set terminal jpeg font arial 12 size 800,600\n" << \
                 "set output \"graph/PP_PARISI_L=" << L << ".jpg\"\n" << \
                 "set grid x y\n" << \
                 "set xlabel \"T\"\n" << \
                 "set ylabel \"PP\"\n" << \
                 "plot \"result/PP_PARISI_L=" << L <<"\" using 2:1 title \"PP_PARISI_L-" << L << "\" with lines lt rgb \"red\"";

    plot_f.close();

    ss.str("");
    ss << "temp/" << "ENERGY_L=" << L << ".plot"; 
    plot_f.open(ss.str().c_str());

    plot_f << "#!/usr/bin/gnuplot -persist\n" << \
                 "set terminal jpeg font arial 12 size 800,600\n" << \
                 "set output \"graph/ENERGY_L=" << L << ".jpg\"\n" << \
                 "set grid x y\n" << \
                 "set xlabel \"T\"\n" << \
                 "set ylabel \"E\"\n" << \
                 "plot \"result/ENERGY_L=" << L <<"\" using 2:1 title \"ENERGY_L-" << L << "\" with lines lt rgb \"red\"";

    plot_f.close();

    ss.str("");
    ss << "temp/" << "C_v_L=" << L << ".plot"; 
    plot_f.open(ss.str().c_str());

    plot_f << "#!/usr/bin/gnuplot -persist\n" << \
                 "set terminal jpeg font arial 12 size 800,600\n" << \
                 "set output \"graph/C_v_L=" << L << ".jpg\"\n" << \
                 "set grid x y\n" << \
                 "set xlabel \"T\"\n" << \
                 "set ylabel \"C_v\"\n" << \
                 "plot \"result/C_v_L=" << L <<"\" using 2:1 title \"C_v_L-" << L << "\" with lines lt rgb \"red\"";

    plot_f.close();

    ss.str("");
    ss << "temp/" << "X_L=" << L << ".plot"; 
    plot_f.open(ss.str().c_str());

    plot_f << "#!/usr/bin/gnuplot -persist\n" << \
                 "set terminal jpeg font arial 12 size 800,600\n" << \
                 "set output \"graph/X_L=" << L << ".jpg\"\n" << \
                 "set grid x y\n" << \
                 "set xlabel \"T\"\n" << \
                 "set ylabel \"X\"\n" << \
                 "plot \"result/X_L=" << L <<"\" using 2:1 title \"X_L-" << L << "\" with lines lt rgb \"red\"";

    plot_f.close();

    ss.str("");
    ss << "PP_PARISI_L=" << L << "-plot.sh";

    script_f.open(ss.str().c_str());

    ss.str("");
    ss << "#!/bin/bash\n";
    ss << "gnuplot temp/PP_PARISI_L=" << L << ".plot\n";

    script_f << ss.str(); 

    script_f.close();

    ss.str("");
    ss << "ENERGY_L=" << L << "-plot.sh";

    script_f.open(ss.str().c_str());
    
    ss.str("");

    ss << "#!/bin/bash\n";
    ss << "gnuplot temp/ENERGY_L=" << L << ".plot\n";

    script_f << ss.str();

    script_f.close();

    ss.str("");
    ss << "C_v_L=" << L << "-plot.sh";

    script_f.open(ss.str().c_str());
    
    ss.str("");

    ss << "#!/bin/bash\n";
    ss << "gnuplot temp/C_v_L=" << L << ".plot\n";

    script_f << ss.str();

    script_f.close();

    ss.str("");
    ss << "X_L=" << L << "-plot.sh";

    script_f.open(ss.str().c_str());
    
    ss.str("");

    ss << "#!/bin/bash\n";
    ss << "gnuplot temp/X_L=" << L << ".plot\n";

    script_f << ss.str();

    script_f.close();

    f_out_pp.close();
    f_out_ee.close();
    f_out_cv.close();
    f_out_x.close();

    delete spin1;
    delete spin2;

    delete J1;
    delete J2;

    return 0;
}

int init_spin_matrix()    {

    // int seed = 1 + rand() % 10000;
    // CRandomMersenne Mersenne(seed);

    spin1 = new int **[L];
    for(int i=0;i<L;i++)    {
        spin1[i] = new int *[L];
        for(int j=0;j<L;j++)    {
            spin1[i][j] = new int [L];
        }
    }

    spin2 = new int **[L];
    for(int i=0;i<L;i++)    {
        spin2[i] = new int *[L];
        for(int j=0;j<L;j++)    {
            spin2[i][j] = new int [L];
        }
    }

    for(int i=0;i<L;i++)    {
        for(int j=0;j<L;j++)    {
            for(int k=0;k<L;k++)    {

                spin1[i][j][k] = 1;

            }
        }
    }

    for(int i=0;i<L;i++)    {
        for(int j=0;j<L;j++)    {
            for(int k=0;k<L;k++)    {

                spin2[i][j][k] = 1;

            }
        }
    }    

    return 0;
}

int init_J()    {

    int seed = 1 + rand() % 10000;
    CRandomMersenne Mersenne(seed);

    J1 = new int **[L];
    for(int i=0;i<L;i++)    {
        J1[i] = new int *[L];
        for(int j=0;j<L;j++)    {
            J1[i][j] = new int [L];
        }
    }

    J2 = new int **[L];
    for(int i=0;i<L;i++)    {
        J2[i] = new int *[L];
        for(int j=0;j<L;j++)    {
            J2[i][j] = new int [L];
        }
    }

    for(int i=0;i<L;i++)    {
        for(int j=0;j<L;j++)    {
            for(int k=0;k<L;k++)    {
                if((double)(Mersenne.IRandom(1,10000)/10000.0)<0.5) J1[i][j][k] = 1;
                else J1[i][j][k] = -1;
            }
        }
    }

    for(int i=0;i<L;i++)    {
        for(int j=0;j<L;j++)    {
            for(int k=0;k<L;k++)    {
                if((double)(Mersenne.IRandom(1,10000)/10000.0)<0.5) J2[i][j][k] = 1;
                else J2[i][j][k] = -1;
            }
        }
    }

    return 0;
}

inline int neighbour_spins1(int i, int j, int k)    {

    int result;

    if(i==0)    result= J1[L-1][j][k]*spin1[L-1][j][k]; 
    else        result= J1[i-1][j][k]*spin1[i-1][j][k]; 
    if(i==L-1)  result+=J1[0][j][k]  *spin1[0][j][k];  
    else        result+=J1[i+1][j][k]*spin1[i+1][j][k];
    if(j==0)    result+=J1[i][L-1][k]*spin1[i][L-1][k];
    else        result+=J1[i][j-1][k]*spin1[i][j-1][k];
    if(j==L-1)  result+=J1[i][0][k]  *spin1[i][0][k];  
    else        result+=J1[i][j+1][k]*spin1[i][j+1][k];
    if(k==0)    result+=J1[i][j][L-1]*spin1[i][j][L-1];
    else        result+=J1[i][j][k-1]*spin1[i][j][k-1];
    if(k==L-1)  result+=J1[i][j][0]  *spin1[i][j][0];
    else        result+=J1[i][j][k+1]*spin1[i][j][k+1];

    // if(i==0)    result= spin1[L-1][j][k]; 
    // else        result= spin1[i-1][j][k]; 
    // if(i==L-1)  result+=spin1[0][j][k];  
    // else        result+=spin1[i+1][j][k];
    // if(j==0)    result+=spin1[i][L-1][k];
    // else        result+=spin1[i][j-1][k];
    // if(j==L-1)  result+=spin1[i][0][k];  
    // else        result+=spin1[i][j+1][k];
    // if(k==0)    result+=spin1[i][j][L-1];
    // else        result+=spin1[i][j][k-1];
    // if(k==L-1)  result+=spin1[i][j][0];
    // else        result+=spin1[i][j][k+1];

    // cout << "result = " << result << endl; 

    return result;

}

inline int neighbour_spins2(int i, int j, int k)   {

    int result;

    if(i==0)    result= J1[L-1][j][k]*spin2[L-1][j][k]; 
    else        result= J1[i-1][j][k]*spin2[i-1][j][k]; 
    if(i==L-1)  result+=J1[0][j][k]  *spin2[0][j][k];  
    else        result+=J1[i+1][j][k]*spin2[i+1][j][k];
    if(j==0)    result+=J1[i][L-1][k]*spin2[i][L-1][k];
    else        result+=J1[i][j-1][k]*spin2[i][j-1][k];
    if(j==L-1)  result+=J1[i][0][k]  *spin2[i][0][k];  
    else        result+=J1[i][j+1][k]*spin2[i][j+1][k];
    if(k==0)    result+=J1[i][j][L-1]*spin2[i][j][L-1];
    else        result+=J1[i][j][k-1]*spin2[i][j][k-1];
    if(k==L-1)  result+=J1[i][j][0]  *spin2[i][j][0];
    else        result+=J1[i][j][k+1]*spin2[i][j][k+1];

    // if(i==0)    result= spin2[L-1][j][k]; 
    // else        result= spin2[i-1][j][k]; 
    // if(i==L-1)  result+=spin2[0][j][k];  
    // else        result+=spin2[i+1][j][k];
    // if(j==0)    result+=spin2[i][L-1][k];
    // else        result+=spin2[i][j-1][k];
    // if(j==L-1)  result+=spin2[i][0][k];  
    // else        result+=spin2[i][j+1][k];
    // if(k==0)    result+=spin2[i][j][L-1];
    // else        result+=spin2[i][j][k-1];
    // if(k==L-1)  result+=spin2[i][j][0];
    // else        result+=spin2[i][j][k+1];

    return result;

}

double pp_parisi()  {

    double sum = 0.0;

    for(int i=0;i<L;i++)    {
        for(int j=0;j<L;j++)    {
            for(int k=0;k<L;k++)    {
                sum += (spin1[i][j][k] * spin2[i][j][k]);
                // cout << "spin1: " << spin1[i][j][k] << "\tspin2: " << spin2[i][j][k] << "\tmult = " << sum << endl;
            }
        }
    }

    // cout << "sum = " << sum/(L*L*L) << endl;
    // cout << "sum = " << sum << endl;
    return sum/(L*L*L);
}

double energy() {

    double energy = 0.0;

    for(int i=0;i<L;i++)    {
        for(int j=0;j<L;j++)    {
            for(int k=0;k<L;k++)    {
                energy -= (spin1[i][j][k]*neighbour_spins1(i,j,k)*spin2[i][j][k]*neighbour_spins2(i,j,k)); 
            }
        }
    }

    return energy/(6*L*L*L);
}

// double energy2() {
    
//     double energy_per_spin = 0.0;

//     for(int i=0;i<L;i++)    {
//         for(int j=0;j<L;j++)    {
//             for(int k=0;k<L;k++)    {
//                 energy_per_spin -= (spin2[i][j][k]*neighbour_spins2(i,j,k)); 
//             }
//         }
//     }

//     return energy_per_spin/(6*L*L*L);
// }

// double magnetic1() {
    
//     double magnetic = 0.0;

//     for(int i=0;i<L;i++)    {
//         for(int j=0;j<L;j++)    {
//             for(int k=0;k<L;k++)    {
//                 magnetic += spin1[i][j][k]; 
//             }
//         }
//     }

//     return magnetic/(L*L*L);
// }

// double magnetic2() {

//     double magnetic = 0.0;

//     for(int i=0;i<L;i++)    {
//         for(int j=0;j<L;j++)    {
//             for(int k=0;k<L;k++)    {
//                 magnetic += spin2[i][j][k]; 
//             }
//         }
//     }

//     return magnetic/(L*L*L);
// }
