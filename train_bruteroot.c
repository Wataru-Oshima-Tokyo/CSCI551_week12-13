// This is the simplest brute force stepping method possible to find a root
// on an interval with a fixed step size equal to error tolerance.
//
// It simply detects a sign change and declares the root to be half way between
// the step before the sign change and the step after and quits.
//
// Inspired by questions from Section 02 CSCI 551, 2022
//
// Why not just search by making small steps for a zero crossing.  This definitely works, but
// does require small step sizes, potentially many steps to search a large interval, and error
// will simply be stepsize / 2.0
//
// Makes an intersting comparison to Newton, Bisection, and Regula Falsi, all which base stepping
// on an estimate of where the root is on an interval.
//
// Generally, Newton should converge fastest, followed by Regula Falsi, then Bisection, and the Bruteroot will
// get the root stepping in one direction looking for a simple sign change.
//
#include<stdio.h>
#include<math.h>
#include <limits.h>
#include <omp.h>
#include <time.h>
#define NUM_THREADS 2
#define OMP
#define _pi 3.141592
// model a function with desmos for example to "see" roots
// 
// https://www.desmos.com/calculator
double f(double x){
    return 0.2361225 * sin((2*_pi*x)/1800);
    // y = 1E-15x6 - 2E-10x5 + 1E-06x4 - 0.0017x3 + 1.0981x2 - 179.85x + 6478.7
    // return (1*pow(10,-15)*(x*x*x*x*x*x));
    // y = 3E-15x6 - 1E-11x5 + 3E-08x4 - 2E-05x3 + 0.0091x2 - 0.8738x + 23.82
    
}

double Left_Riemann_sum(double velocity,long double delta){
    long double width, height, area;
    width = delta;
    height = velocity;
    area = width*height;
    if (area<0) area*=-1;
    // printf("%Lf\n", area);
    return area;
}

int main(void)
{
    int itr, maxitr, rootfound=0;
    double step, x0, x1, x2, sign, start;
    struct timespec _start, _stop;
    double fstart, fstop;
    printf("\nEnter x0 start, x1 end, step size (allowed error is 1/2 step size)\n");
    scanf("%lf %lf %lf", &x0, &x2, &step);
    start=x0;

    printf("Will step at %lf from %lf to %lf to find root\n", step, x0, x2);
    maxitr = fabs(x2-x0)/step;
    clock_gettime(CLOCK_MONOTONIC, &_start); fstart=(double)_start.tv_sec + ((double)_start.tv_nsec/1000000000.0);
#ifdef OMP
    #pragma omp parallel for num_threads(NUM_THREADS) default(none) reduction(+:rootfound) private(x1, x0, itr, sign) shared(maxitr, step)
#endif
    for (itr=1; itr<=maxitr; itr++)
    {

        x1=x0+step;

        // Any negative value x positive or vice versa will result in a sign less than zero - a crossing
        sign = f(x0) * f(x1);

        if(sign < 0.0)
        {
            // if(x1 >x2) break;
            printf("Sign change at Iteration no. %3d, x = %20.15f, root estimated at %20.15lf\n", itr, x1, (x1+x0)/2.0);
            rootfound++;

            // We can exit the loop on first root (zero crossing found) or keep looking
            // break;
            
        }

        x0=x1;
    }

    if(!rootfound){
        printf("After %d iterations: No solution (zero crossing) found on interval %lf to %lf\n", maxitr, start, x0);    
    }
    clock_gettime(CLOCK_MONOTONIC, &_stop); fstop=(double)_stop.tv_sec + ((double)_stop.tv_nsec/1000000000.0);
    printf("The total number of roots is %d\n",rootfound);
    printf("Finished in %lf seconds\n",(fstop-fstart));
    return 1;
}
