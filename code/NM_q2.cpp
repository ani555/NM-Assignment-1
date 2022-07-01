#include<iostream>
#include<math.h>
#include<utility>
using namespace std;
#define inf 1e8
#define MAXITERS 100000

// Functions
double f1(double x){
    return x*sin(x) + 3*cos(x) - x;
}

double f2(double x){
    return sin(x) - 0.1*x;
}

double f3(double x){
    return (x+2) * pow(x-1,4);
}

// Derivatives
double df1(double x){
    return x*cos(x) - 2*sin(x) -1;
}

double df2(double x){
    return cos(x) - 0.1;
}

double df3(double x){
    return pow(x-1,4) + 4*(x+2)*pow(x-1,3);
}

// Double derivatives
double d2f1(double x){
    return -x*sin(x) - cos(x);
}

double d2f2(double x){
    return -sin(x);
}

double d2f3(double x){
    return 8*pow(x-1,3) + 12*(x+2)*pow(x-1,2); 
}

double rel_error(double p_act, double p_eval){
    return abs(p_act - p_eval)/abs(p_eval);
}

/*
This function calculates the root of a function f using newtons method. It takes input the pointer to function f, its derivative df, intial guess x0 and tolerance value tol which is optional defaulting to 1e-06, it then returns the root and niters as a pair
*/
pair<double,int> newton_raphson(double (*f)(double), double (*df)(double),double x0, double tol=1e-06){
    double rel_err = inf;
    double x1;
    int niters=0;
    while(rel_err > tol){
        x1 = x0 - (*f)(x0)/(*df)(x0);
        rel_err = rel_error(x1, x0);
        x0 = x1;
        niters++;
        if(niters == MAXITERS){
            printf("Solution did not converge in %d iterations",MAXITERS);
            break;
        }
        //printf("iter# %d Error: %1.6e\n",niters,rel_err);
    }
    //printf("\n");
    return make_pair(x1, niters);
}

/* 
This function calculates the root of a function f using secant method. It takes input the pointer to function f, intial guesses x0 and x1 and tolerance value tol which is optional defaulting to 1e-06, it then returns the root and niters as a pair 
*/
pair<double, int> secant_method(double (*f)(double), double x0, double x1, double tol=1e-6){
    double rel_err = inf;
    double x2;
    int niters = 0;
    while(rel_err > tol){
        x2 = x1 - (x1 - x0)*(*f)(x1)/((*f)(x1)-(*f)(x0));
        rel_err = rel_error(x2, x1);
        x0 = x1;
        x1 = x2;
        niters++;
        if(niters == MAXITERS){
            printf("Solution did not converge in %d iterations",MAXITERS);
            break;
        }
    }
    return make_pair(x2, niters);
}

/*
This function calculates the root of a function f using newtons method. It takes input the pointer to function f, its derivative df and double derivative d2f, intial guess x0 and tolerance value tol which is optional defaulting to 1e-06, it then returns the root and niters as a pair
*/
pair<double, int> modified_newton(double (*f)(double), double (*df)(double), double (*d2f)(double), double x0, double tol=1e-6){
    double rel_err = inf;
    double x1;
    int niters=0;
    while(rel_err > tol){
        x1 = x0 - (*f)(x0)*(*df)(x0)/((*df)(x0)*(*df)(x0) - (*f)(x0)*(*d2f)(x0));
        rel_err = rel_error(x1, x0);
        x0 = x1;
        niters++;
        if(niters == MAXITERS){
            printf("Solution did not converge in %d iterations",MAXITERS);
            break;
        }
    }
    return make_pair(x1, niters);
}
int main(){
    
    printf("-----Newton Raphson-----\n");
    //test newton raphson to find roots of f1 between (-6,6)
    pair<double, int> root1 = newton_raphson(f1, df1, 1.5);
    pair<double, int> root2 = newton_raphson(f1, df1, -2);
    pair<double, int> root3 = newton_raphson(f1, df1, -5);
    printf("Roots of x*sin(x) + 3*cos(x) - x between (-6,6)\n");    
    printf("Root 1 = %.6lf #iters = %d\n",root1.first, root1.second);
    printf("Root 2 = %.6lf #iters = %d\n",root2.first, root2.second);
    printf("Root 3 = %.6lf #iters = %d\n",root3.first, root3.second);
    
    //test newton raphson to find all positive roots of f2
    root1 = newton_raphson(f2, df2, 3);
    root2 = newton_raphson(f2, df2, 7);
    root3 = newton_raphson(f2, df2, 8);
    printf("Positive roots of sin(x) - 0.1*x\n");
    printf("Root 1 = %.6lf #iters = %d\n",root1.first, root1.second);
    printf("Root 2 = %.6lf #iters = %d\n",root2.first, root2.second);
    printf("Root 3 = %.6lf #iters = %d\n",root3.first, root3.second);
    
    printf("\n-----Secant Method-----\n");
    //test secant method to find roots of f1 between (-6,6)
    root1 = secant_method(f1, 1.5, 1.6);
    root2 = secant_method(f1, -2.5, -2);
    root3 = secant_method(f1,  -4.5, -5);
    printf("Roots of x*sin(x) + 3*cos(x) - x between (-6,6)\n");    
    printf("Root 1 = %.6lf #iters = %d\n",root1.first, root1.second);
    printf("Root 2 = %.6lf #iters = %d\n",root2.first, root2.second);
    printf("Root 3 = %.6lf #iters = %d\n",root3.first, root3.second);
    
    //test the secant method to find all positive roots of f2
    root1 = secant_method(f2, 3, 2.5);
    root2 = secant_method(f2, 7.2, 7);
    root3 = secant_method(f2, 8.5, 8);
    printf("Positive roots of sin(x) - 0.1*x\n");
    printf("Root 1 = %.6lf #iters = %d\n",root1.first, root1.second);
    printf("Root 2 = %.6lf #iters = %d\n",root2.first, root2.second);
    printf("Root 3 = %.6lf #iters = %d\n",root3.first, root3.second);

    printf("\n-----Modified Newton-----\n");
    //test moified newton method to find roots of f1 between (-6,6)
    root1 = modified_newton(f1, df1, d2f1, 1.5);
    root2 = modified_newton(f1, df1, d2f1, -2);
    root3 = modified_newton(f1, df1, d2f1, -5);
    printf("Roots of x*sin(x) + 3*cos(x) - x between (-6,6)\n");    
    printf("Root 1 = %.6lf #iters = %d\n",root1.first, root1.second);
    printf("Root 2 = %.6lf #iters = %d\n",root2.first, root2.second);
    printf("Root 3 = %.6lf #iters = %d\n",root3.first, root3.second);
    
    //test the modified newton method to find all positive roots of f2
    root1 = modified_newton(f2, df2, d2f2, 3);
    root2 = modified_newton(f2, df2, d2f2, 7);
    root3 = modified_newton(f2, df2, d2f2, 8);
    printf("Positive roots of sin(x) - 0.1*x\n");
    printf("Root 1 = %.6lf #iters = %d\n",root1.first, root1.second);
    printf("Root 2 = %.6lf #iters = %d\n",root2.first, root2.second);
    printf("Root 3 = %.6lf #iters = %d\n",root3.first, root3.second);

    // newtons method vs modified newtons method in case of multiple roots
    printf("\nNewton vs modified newton for f = (x+2)*(x-1)^4\n");
    pair<double, int> root = newton_raphson(f3, df3, 1.2);
    printf("Newtons method :  Root = %.6lf #iters = %d\n",root.first, root.second);
    root = modified_newton(f3, df3, d2f3, 1.2);
    printf("Modified newton : Root = %.6lf #iters = %d\n",root.first, root.second);    
    
    
}
