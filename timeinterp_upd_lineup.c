#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>
// For values between 1 second indexed data, use linear interpolation to determine profile value at any "t".
//
// Wikipedia - https://en.wikipedia.org/wiki/Linear_interpolation
//
// This assumes that the profile is "piecewise linear", which is accurate for linear or constant functions, and 
// reasonably accurate for most non-linear functions over small intervals of 1 second.
//
// This test program simply interpolates values between 1 second intervals to produce a value every 1/10th of 
// a second.

// #include "ex3accel.h"
// #include "ex6_linear.h"
// #define LINEAR
// #define OMP


#include "_ex6_linear.h"
#include "_ex6_non_linear.h"


#define DELTA 0.1
#define NUM_THREADS 4
#define OMP
// #define VELOCITY
struct timespec start, stop;
double fstart, fstop;
double lin_velocities[1801];
double nonlin_velocities[1801];
// table look-up for acceleration profile given and velocity profile determined
//
// Note: for 2 functions (2 trains) we would want to make 2 different versions of this
//       function or better yet, pass in the table to use.
//
double table_accel(int timeidx, int lin);
double table_velocity(int timeidx, int lin);
double faccel(double time, int lin);
double fvelocity(double time, int lin);


double Left_Riemann_sum(double velocity, long double start, long double delta){
    long double width, height, area;
    width = delta;
    height = velocity;
    area = width*height;
    return area;
}


int main(int argc, char *argv[])
{
    int idx;
    double dt = DELTA;
    double time, steps_per_sec;
    double lvnot=0, nlvnot=0;
    double lin_left_Riemann_sum = 0, nlin_left_Riemann_sum = 0;;

    if(argc == 2)
    {
        sscanf(argv[1], "%le", &dt);
        // printf("Will use %d steps per sec\n", steps_per_sec);
    }
    #ifdef OMP
        printf("With OpenMP\n");
    #else
        printf("Without OpenMP\n");
    #endif
    steps_per_sec = 1.0 /dt ;
    int maxitr = steps_per_sec*1800;
    
    //change accel to velocity so that I can calculate the area under the curve by left_Riemann sum
    for(int i=0; i<1801; i++){
        lin_velocities[i] = faccel(i,1)*1+lvnot;
        nonlin_velocities[i] = faccel(i,0)*1+nlvnot;
        lvnot = lin_velocities[i];
        nlvnot = nonlin_velocities[i];
    }

    
    clock_gettime(CLOCK_MONOTONIC, &start); fstart=(double)start.tv_sec + ((double)start.tv_nsec/1000000000.0);
    #ifdef OMP
        #pragma omp parallel for num_threads(NUM_THREADS) default(none) reduction(+:lin_left_Riemann_sum,nlin_left_Riemann_sum) private(time) shared(dt, idx,maxitr)
    #endif
    for(idx=0; idx <= maxitr; idx++)
    {
        // time you would use in your integrator and faccel(time) is the fuction to integrate
        time = 0.0 + (dt*(double)idx);
        lin_left_Riemann_sum += Left_Riemann_sum(fvelocity(time,1), time, dt);
        nlin_left_Riemann_sum += Left_Riemann_sum(fvelocity(time,0), time, dt);
        #ifdef VELOCITY
            if((fabs(fvelocity(time,1) - fvelocity(time,0))<1) && time > 400 && time <1400){
                #ifdef OMP
                    printf("The two trains line up when the velocites are %lf(linear) m/sec and %lf(non linear) m/sec at %lf seconds\n",fvelocity(time,1),fvelocity(time,0), time);
                #else
                    printf("The two trains line up for the first time on the invterval %d to %d when the velocites are %lf(linear) m/sec and %lf(non linear) m/sec at %lf seconds\n", 400, 1400, fvelocity(time,1),fvelocity(time,0), time);
                    break;
                #endif
            }
            
        #else
            if((fabs(lin_left_Riemann_sum - nlin_left_Riemann_sum)<1) && time > 600 && time <800){
                #ifdef OMP
                    printf("The two trains line up when the location is %lf meters at %lf seconds\n", lin_left_Riemann_sum, time);
                #else
                    printf("The two trains line up for the first time on the invterval %d to %d when the location is %lf meters at %lf seconds\n", 600, 800, lin_left_Riemann_sum, time);
                    break;
                #endif
             }
        #endif

    }
    clock_gettime(CLOCK_MONOTONIC, &stop); fstop=(double)stop.tv_sec + ((double)stop.tv_nsec/1000000000.0);
    printf("The toatal time is %f\n", (fstop-fstart));
    printf("The linear left riemann sum is %f\n", lin_left_Riemann_sum);
    printf("The non linear left riemann sum is %f\n", nlin_left_Riemann_sum);
    // printf("The trapozoidal sum is %f\n", trapezoidal_sum);
    return 0;
}


// Simple look-up in accleration profile array
//
// Added array bounds check for known size of train arrays
//
double table_accel(int timeidx, int lin)
{
    long unsigned int tsize = sizeof(Lin_DefaultProfile) / sizeof(double);

    // Check array bounds for look-up table
    if(timeidx > tsize)
    {
        printf("timeidx=%d exceeds table size = %lu and range %d to %lu\n", timeidx, tsize, 0, tsize-1);
        return 0;
        // exit(-1);
    }
    if(lin==1)
        return Lin_DefaultProfile[timeidx];
    else   
        return NLin_DefaultProfile[timeidx];
}


// Simple linear interpolation example for table_accel(t) for any floating point t value
// for a table of accelerations that are 1 second apart in time, evenly spaced in time.
//
// accel[timeidx] <= accel[time] < accel[timeidx_next]
//
//
double faccel(double time, int lin)
{
    // The timeidx is an index into the known acceleration profile at a time <= time of interest passed in
    //
    // Note that conversion to integer truncates double to next lowest integer value or floor(time)
    //
    int timeidx = (int)time;

    // The timeidx_next is an index into the known acceleration profile at a time > time of interest passed in
    //
    // Note that the conversion to integer truncates double and the +1 is added for ceiling(time)
    //
    int timeidx_next = ((int)time)+1;

    // delta_t = time of interest - time at known value < time
    //
    // For more general case
    // double delta_t = (time - (double)((int)time)) / ((double)(timeidx_next - timeidx);
    //
    // If time in table is always 1 second apart, then we can simplify since (timeidx_next - timeidx) = 1.0 by definition here
    double delta_t = time - (double)((int)time);

    return ( 
               // The accel[time] is a linear value between accel[timeidx] and accel[timeidx_next]
               // 
               // The accel[time] is a value that can be determined by the slope of the interval and accel[timedix] 
               //
               // I.e. accel[time] = accel[timeidx] + ( (accel[timeidx_next] - accel[timeidx]) / ((double)(timeidx_next - timeidx)) ) * delta_t
               //
               //      ((double)(timeidx_next - timeidx)) = 1.0
               // 
               //      accel[time] = accel[timeidx] + (accel[timeidx_next] - accel[timeidx]) * delta_t
               //
               table_accel(timeidx, lin) + ( (table_accel(timeidx_next, lin) - table_accel(timeidx, lin)) * delta_t)
           );
}


double table_velocity(int timeidx, int lin)
{
    long unsigned int tsize = sizeof(lin_velocities) / sizeof(double);

    // Check array bounds for look-up table
    if(timeidx > tsize)
    {
        printf("timeidx=%d exceeds table size = %lu and range %d to %lu\n", timeidx, tsize, 0, tsize-1);
        return 0;
        // exit(-1);
    }
    if(lin==1)
        return lin_velocities[timeidx];
    else
        return  nonlin_velocities[timeidx];
}


double fvelocity(double time, int lin)
{
    // The timeidx is an index into the known acceleration profile at a time <= time of interest passed in
    //
    // Note that conversion to integer truncates double to next lowest integer value or floor(time)
    //
    int timeidx = (int)time;

    // The timeidx_next is an index into the known acceleration profile at a time > time of interest passed in
    //
    // Note that the conversion to integer truncates double and the +1 is added for ceiling(time)
    //
    int timeidx_next = ((int)time)+1;

    // delta_t = time of interest - time at known value < time
    //
    // For more general case
    // double delta_t = (time - (double)((int)time)) / ((double)(timeidx_next - timeidx);
    //
    // If time in table is always 1 second apart, then we can simplify since (timeidx_next - timeidx) = 1.0 by definition here
    double delta_t = time - (double)((int)time);

    return ( 
               // The accel[time] is a linear value between accel[timeidx] and accel[timeidx_next]
               // 
               // The accel[time] is a value that can be determined by the slope of the interval and accel[timedix] 
               //
               // I.e. accel[time] = accel[timeidx] + ( (accel[timeidx_next] - accel[timeidx]) / ((double)(timeidx_next - timeidx)) ) * delta_t
               //
               //      ((double)(timeidx_next - timeidx)) = 1.0
               // 
               //      accel[time] = accel[timeidx] + (accel[timeidx_next] - accel[timeidx]) * delta_t
               //
               
               table_velocity(timeidx, lin) + ( (table_velocity(timeidx_next, lin) - table_velocity(timeidx, lin)) * delta_t)
           );
}