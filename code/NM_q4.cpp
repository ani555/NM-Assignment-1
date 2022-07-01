#include<iostream>
#include<math.h>
#include<utility>
using namespace std;
#define inf 1e8
#define MAXITERS 100000

double f1(double x, double y){
    return sin(x) + 3*cos(y) -2 ; 
}

double f2(double x, double y){
    return cos(x) - sin(y) + 0.2;
}

double df1x(double x, double y){
    return cos(x);
}

double df1y(double x, double y){
    return -3*sin(y);
}

double df2x(double x, double y){
    return -sin(x);
}

double df2y(double x, double y){
    return -cos(y);
}

double rel_error(pair<double, double> p_act, pair<double, double> p_eval){
    double x0 = p_act.first, y0 = p_act.second;
    double x1 = p_eval.first, y1 = p_eval.second;
    return sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))/sqrt(x0*x0 + y0*y0);
}

pair<pair<double,double>, int> newtons_method(double (*f1)(double, double), double (*f2)(double, double), double (*df1x)(double, double), double (*df1y)(double, double), double (*df2x)(double, double), double (*df2y)(double, double), pair<double, double> p0, double tol = 1e-8){
    double x0 = p0.first;
    double y0 = p0.second;
    double x1, y1;
    double invJ11, invJ12, invJ21, invJ22, det;
    int niters=0;
    double rel_err = inf;
    while(rel_err > tol){
        
        invJ11 = (*df2y)(x0, y0);
        invJ22 = (*df1x)(x0, y0);
        invJ12 = -(*df1y)(x0, y0);
        invJ21 = -(*df2x)(x0, y0);
        det = invJ11*invJ22 - invJ12*invJ21;
        
        x1 = x0 - (invJ11 * (*f1)(x0, y0) + invJ12*(*f2)(x0, y0))/det;
        y1 = y0 - (invJ21 * (*f1)(x0, y0) + invJ22*(*f2)(x0, y0))/det;
        rel_err = rel_error(make_pair(x0,y0), make_pair(x1, y1));
        x0 = x1;
        y0 = y1;
        niters++;
        if(niters == MAXITERS){
            printf("Solution did not converge in %d iterations",MAXITERS);
            break;
        }
    }
    return make_pair(make_pair(x1, y1), niters);
}

int main(){
    //test newtons method for multivariable non linear eqn
    pair<pair<double, double>, int> res = newtons_method(f1, f2, df1x, df1y, df2x, df2y, make_pair(1.0,1.0));
    printf("Root: x = %.8lf y = %.8lf \n", res.first.first, res.first.second);
    printf("#iters = %d\n", res.second);
}
