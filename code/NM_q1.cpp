#include<iostream>
#include<assert.h>
#include<math.h>
using namespace std;
#define MAXN 10

// This function computes the bessel function in forward direction
double Jf(int x, int n){
    
    assert(x==1 || x==5 || x==50);
    assert(n>=0 && n<=MAXN);
    double J0, J1, J2;
    if(x==1){
        J0 = 7.65198e-01;
        J1 = 4.40051e-01;
    }
    else if(x==5){
        J0 = -1.77597e-01;
        J1 = -3.27579e-01;
    }
    else{
        J0 = 5.58123e-02;
        J1 = -9.75118e-02;
    }
    
    if(n==0) return J0;
    else if(n==1) return J1;
    
    for(int i=1; i<=n-1; i++){
        J2 =(2.0*(i)/(double)x)*J1 - J0;
        J0 = J1;
        J1 = J2;
    }
    return J2;
}

// This function computes the bessel function in backward direction
double Jb(int x, int n){
    assert(x==1 || x==5 || x==50);
    assert(n>=0 && n<=10);
    
    double J10, J9, J8;
    if(x==1){
        J10 = 2.63062e-10;
        J9 = 5.24925e-09;
    }
    else if(x==5){
        J10 = 1.46780e-03;
        J9 = 5.52028e-03;
    }
    else{
        J10 = -1.138478e-01;
        J9 = -2.71924e-02;
    }
    if(n==10) return J10;
    else if(n==9) return J9;
    
    for(int i=MAXN-1; i>n; i--){
        J8 = (2.0*(i)/(double)x)*J9 - J10;
        J10 = J9;
        J9 = J8;
    }
    return J8;
}

double rel_error(double p_act, double p_eval){
    return abs(p_act - p_eval)/abs(p_eval);
}

double abs_error(double p_act, double p_eval){
    return abs(p_act - p_eval);
}

int main(){
    //test the forward computation
    double J10_1_act = 2.6306151237e-10;
    double J10_5_act = 1.4678026473e-03;
    double J10_50_act = -1.1384784915e-01;
    double J10_1_eval = Jf(1,10);
    double J10_5_eval = Jf(5,10);
    double J10_50_eval = Jf(50,10);
    printf("Calculated values for x = 1, 5, 50 resp :\n");
    printf("%1.10e %1.10e %1.10e\n",J10_1_eval, J10_5_eval, J10_50_eval);
    printf("J10 : x=1 abs_error=%1.10e rel_error=%1.10e\n",abs_error(J10_1_act,J10_1_eval), rel_error(J10_1_act, J10_1_eval));
    printf("J10 : x=5 abs_error=%1.10e rel_error=%1.10e\n",abs_error(J10_5_act,J10_5_eval), rel_error(J10_5_act, J10_5_eval));
    printf("J10 : x=50 abs_error=%1.10e rel_error=%1.10e\n",abs_error(J10_1_act,J10_50_eval), rel_error(J10_50_act, J10_50_eval));
    
    //test the backward computation
    double J0_1_act = 7.6519768656e-01;
    double J0_5_act = -1.7759677131e-01;
    double J0_50_act = 5.5812327669e-02;
    double J0_1_eval = Jb(1,0);
    double J0_5_eval = Jb(5,0);
    double J0_50_eval = Jb(50,0);
    printf("Calculated values for x = 1, 5, 50 resp :\n");
    printf("%1.10e %1.10e %1.10e\n",J0_1_eval, J0_5_eval, J0_50_eval);
    printf("J0 : x=1 abs_error=%1.10e rel_error=%1.10e\n",abs_error(J0_1_act,J0_1_eval), rel_error(J0_1_act, J0_1_eval));
    printf("J0 : x=5 abs_error=%1.10e rel_error=%1.10e\n",abs_error(J0_5_act,J0_5_eval), rel_error(J0_5_act, J0_5_eval));
    printf("J0 : x=50 abs_error=%1.10e rel_error=%1.10e\n",abs_error(J0_1_act,J0_50_eval), rel_error(J0_50_act, J0_50_eval));
    
}
    