#include <math.h>
#include <stdlib.h>
#include <stdio.h>
// #include <sys/sem.h>
// #include <sys/types.h>
// #include <sys/wait.h>
// #include <semaphore.h>

// #include <unistd.h> 
// #include <sys/mman.h>


//#include "mkl.h"

double * timeDifference(double *, double*, int, int, int);
double * timeDifference_fast(double *, double*, int, int, int);
double * timeDifference_fast_double(double *, double*, int, double, int);
void timeDifference_inPlace(double *, double*, double*, int, int, int);
void cleanup(void *);
// double * timeDifference_multi(double *, double *, int, int, int, int, int);
// double * timeDifference_multi(double *, double *, int, int, int, int, int);

int fftSize(int, int);
int fftSize_double(int, double);
// double * linearFit(double *, int);
void wrap_qsort(float *, int);
int cmpfunc(const void *, const void *);


// A function to release malloc'ed memory.
void cleanup(void * object)
{
    free(object);
    return;
}

// We need a function to pass to qsort. This is it.
// I got this from stack exchange
// https://stackoverflow.com/questions/27284185/how-does-the-compare-function-in-qsort-work
int cmpfunc (const void * a, const void * b)
{
   // qsort() passes in `void*` types because it can't know the actual types being sorted
   // convert those pointers to pointers to int and deref them to get the actual int values

   float val1 = *(float*)a;
   float val2 = *(float*)b;

   // qsort() expects the comparison function to return:
   // 
   //    a negative result if val1 < val2
   //    0 if val1 == val2
   //    a positive result if val1 > val2

   return ( val1 - val2 ); 
}

// Wrap qsort so it can be easily called from python.
void wrap_qsort (float * input, int size){

    // Not much to do here. Just call qsort
    qsort(input, size, sizeof(float), cmpfunc);
}


// Perform a linear fit to a sorted power array
// But, don't actually use this function. It is much faster to calculate
// the slope and Y-intercept from the correlation and standard deviations
// of the datapoints.

// // Also, this is WRONG. B needs to fill with the log of ii or something.
// double * linearFit(double * sortPower, int powerSize){

//     // A fast way of performing least squares linear regression: 
//     // Ax = B
//     // Where I want to end up with a solution of the form y = mx + b
//     // A is of the form [1, power1,
//     //                   1, power2, 
//     //                   1, power3]
//     // x is of the form [b,
//     //                   m]
//     // B is of the form [0,
//     //                   1,
//     //                   2]

//     printf("Calling linearFit");

//     // Initialize MKL integers for use in the mkl function
//     // rows: The number of rows in A, the same as the number of power values
//     // columns: the number of columns in A, which will be 2 [1, POWER_VALUE]
//     // nrhs: The number of columns in B, which will be 1
//     // lda: The leading dimension of A, equal to max(1, columns) for row major, or 2
//     // ldb: THe leading dimension of B, equal 1.
//     MKL_INT rows = powerSize, columns = 2, nrhs = 1, lda = 2, ldb = 1, info;

//     // A place to store the input and output matrices
//     double *A;
//     double *B;

//     // Add a leading column of 1's
//     // Even though in our example A has 2 columns, the lapacke function
//     // takes in a single long vector and wraps it appropriately
//     A = (double *)calloc(2 * powerSize, sizeof(double*));

//     // This is how to create a powerSize array
//     B = (double *)calloc(powerSize, sizeof(double));

//     // This may take a lot of memory. I should check that it is all okay. 
//     // if memory cannot be allocated
//     if(A == NULL)                     
//     {
//         printf("Creation of B didn't work.");
//         exit(0);
//     }
//     if(B == NULL)                     
//     {
//         printf("Creation of B didn't work.");
//         exit(0);
//     }

//     printf("Filling B");
//     // Fill the arrays.
//     for (int ii = 0; ii < powerSize; ii++){
//         // Just a list of ascending values.
//         B[ii] = ii;
//     }

//     printf("Filling A");
//     // Create the vector for a that will be appropriately wrapped
//     for (int ii = 0; ii < 2*powerSize; ii+=2){
//         A[ii] = 1;
//         A[ii + 1] = sortPower[ii];

//     }


//     printf("Calling the fitter.");
//     // Solve the linear system
//     info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', rows, columns, nrhs, A, lda, B, ldb);


//     /* Check for the full rank */
//     if( info > 0 ) {
//             printf( "The diagonal element %i of the triangular factor ", info );
//             printf( "of A is zero, so that A does not have full rank;\n" );
//             printf( "the least squares solution could not be computed.\n" );
//             exit( 1 );
//     }

//     return B;
// }

// Calculate the size of an FFT
int fftSize(int windowSize, int maxFreq){
    return 2 * floor(windowSize * maxFreq);
}

//Calculate the size of an FFT, but allow for a double maxFreq
int fftSize_double(int windowSize, double maxFreq){
    return 2 * floor(windowSize * maxFreq);
}

// Calculates a histogram of time differences
double * timeDifference(double *photonTimes, double *photonWeights,
                        int windowSize, int maxFreq, int nPhotons){

    double *outHistogram;
    int nbins;
    float timeResol;
    double photonDiff;
    int freqBin;

    //Calculate the time resolution
    timeResol = 0.5 / maxFreq;

    //calculate the size of the fft
    nbins = fftSize(windowSize, maxFreq);

    // Create an array where we will store the output
    outHistogram = (double *)calloc(nbins, sizeof(double));

    // if memory cannot be allocated
    if(outHistogram == NULL)                     
    {
        printf("Creation of outHistogram didn't work.");
        exit(0);
    }

    // Loop through all the photons, except for the last one
    for (int ii = 0; ii < nPhotons - 1; ii++)
    {

        for (int jj = ii + 1; jj < nPhotons; jj++)
        {
            // Calculate the time difference  
            photonDiff = photonTimes[jj] - photonTimes[ii];

            // Exit this second loop if the time difference is too large
            if (photonDiff >= windowSize)
                break;

            // Otherwise, calculate the frequency bin
            freqBin = floor(photonDiff / timeResol);

            // Store the data in the output
            outHistogram[freqBin] += photonWeights[ii] * photonWeights[jj];
        }
    }

    return outHistogram;

}

// Calculates a histogram of time differences... quickly.
void timeDifference_inPlace(double *photonTimes, double *photonWeights,
                            double *outHistogram,
                            int windowSize, int maxFreq, int nPhotons){

    int nbins;
    float timeResol;
    double photonDiff;
    int freqBin;
    int skip=1;

    //Calculate the time resolution
    timeResol = 0.5 / maxFreq;

    // Loop through all the photons, except for the last one
    for (int ii = 0; ii < nPhotons - 1; ii++)
    {
        // printf("\n");
        // printf("ii is: %d\n", ii);
        // printf("Skip at start of ii: %d\n", skip);

        // If photon 0 differenced up to photon 100 before hitting the
        // window size, then I know for sure that photon 1 will also
        // difference up to photon 100 before hitting the window size.
        // Thus, I can skip checking photonDiff >= windowSize for the 
        // first 99 photons of the photon 1 time differencing.
        // "skip" impliments this skipping.

        // Quickly loop through the photons that I already know are okay.
        for (int jj = ii + 1; jj < (int)fmin(ii + skip, nPhotons); jj++)
        {
            // printf("In first loop. Skip is %d\n", skip);
            // printf("jj is: %d\n", jj);
            // Calculate the time difference  
            photonDiff = photonTimes[jj] - photonTimes[ii];

            // Otherwise, calculate the frequency bin
            freqBin = floor(photonDiff / timeResol);

            // Store the data in the output
            outHistogram[freqBin] += photonWeights[ii] * photonWeights[jj];
        }

        // Loop through additional photons.
        for (int jj = ii + skip; jj < nPhotons; jj++)
        {
            // printf("In second loop. Skip is %d\n", skip);
            // printf("jj is %d\n", jj);

            // Calculate the time difference  
            photonDiff = photonTimes[jj] - photonTimes[ii];
            // printf("photonDiff is: %f\n", photonDiff);

            // Exit this second loop if the time difference is too large
            if (photonDiff >= windowSize)
            {
                // printf("Exiting. Photon diff is %d, window size is %d.\n", photonDiff, windowSize);
                skip = jj - ii - 1;
                //printf("skip at exit is %d\n", skip);
                break;
            }

            // Otherwise, calculate the frequency bin
            freqBin = floor(photonDiff / timeResol);

            // Store the data in the output
            outHistogram[freqBin] += photonWeights[ii] * photonWeights[jj];

        }
    }

    return;

}


// Calculates a histogram of time differences... quickly.
double * timeDifference_fast(double *photonTimes, double *photonWeights,
                             int windowSize, int maxFreq, int nPhotons){

    double *outHistogram;
    int nbins;
    float timeResol;
    double photonDiff;
    int freqBin;
    int skip=1;

    //Calculate the time resolution
    timeResol = 0.5 / maxFreq;

    //calculate the size of the fft
    nbins = fftSize(windowSize, maxFreq);

    // Create an array where we will store the output
    // Actually, let C create this, and I will just modify it in place.
    // I think this will help with memory leaks.
    outHistogram = (double *)calloc(nbins, sizeof(double));

    // if memory cannot be allocated
    if(outHistogram == NULL)                     
    {
        printf("Creation of outHistogram didn't work.");
        exit(0);
    }

    // Loop through all the photons, except for the last one
    for (int ii = 0; ii < nPhotons - 1; ii++)
    {
        // printf("\n");
        // printf("ii is: %d\n", ii);
        // printf("Skip at start of ii: %d\n", skip);

        // If photon 0 differenced up to photon 100 before hitting the
        // window size, then I know for sure that photon 1 will also
        // difference up to photon 100 before hitting the window size.
        // Thus, I can skip checking photonDiff >= windowSize for the 
        // first 99 photons of the photon 1 time differencing.
        // "skip" impliments this skipping.

        // Quickly loop through the photons that I already know are okay.
        for (int jj = ii + 1; jj < (int)fmin(ii + skip, nPhotons); jj++)
        {
            // printf("In first loop. Skip is %d\n", skip);
            // printf("jj is: %d\n", jj);
            // Calculate the time difference  
            photonDiff = photonTimes[jj] - photonTimes[ii];

            // Otherwise, calculate the frequency bin
            freqBin = floor(photonDiff / timeResol);

            // Store the data in the output
            outHistogram[freqBin] += photonWeights[ii] * photonWeights[jj];
        }

        // Loop through additional photons.
        for (int jj = ii + skip; jj < nPhotons; jj++)
        {
            // printf("In second loop. Skip is %d\n", skip);
            // printf("jj is %d\n", jj);

            // Calculate the time difference  
            photonDiff = photonTimes[jj] - photonTimes[ii];
            // printf("photonDiff is: %f\n", photonDiff);

            // Exit this second loop if the time difference is too large
            if (photonDiff >= windowSize)
            {
                // printf("Exiting. Photon diff is %d, window size is %d.\n", photonDiff, windowSize);
                skip = jj - ii - 1;
                //printf("skip at exit is %d\n", skip);
                break;
            }

            // Otherwise, calculate the frequency bin
            freqBin = floor(photonDiff / timeResol);

            // Store the data in the output
            outHistogram[freqBin] += photonWeights[ii] * photonWeights[jj];

        }
    }

    return outHistogram;

}

// Calculates a histogram of time differences... quickly.
// Allows for a double max frequency
double * timeDifference_fast_double(double *photonTimes, double *photonWeights,
                             int windowSize, double maxFreq, int nPhotons){

    double *outHistogram;
    int nbins;
    float timeResol;
    double photonDiff;
    int freqBin;
    int skip=1;

    //Calculate the time resolution
    timeResol = 0.5 / maxFreq;

    //calculate the size of the fft
    nbins = fftSize_double(windowSize, maxFreq);

    // Create an array where we will store the output
    // Actually, let C create this, and I will just modify it in place.
    // I think this will help with memory leaks.
    outHistogram = (double *)calloc(nbins, sizeof(double));

    // if memory cannot be allocated
    if(outHistogram == NULL)                     
    {
        printf("Creation of outHistogram didn't work.");
        exit(0);
    }

    // Loop through all the photons, except for the last one
    for (int ii = 0; ii < nPhotons - 1; ii++)
    {
        // printf("\n");
        // printf("ii is: %d\n", ii);
        // printf("Skip at start of ii: %d\n", skip);

        // If photon 0 differenced up to photon 100 before hitting the
        // window size, then I know for sure that photon 1 will also
        // difference up to photon 100 before hitting the window size.
        // Thus, I can skip checking photonDiff >= windowSize for the 
        // first 99 photons of the photon 1 time differencing.
        // "skip" impliments this skipping.

        // Quickly loop through the photons that I already know are okay.
        for (int jj = ii + 1; jj < (int)fmin(ii + skip, nPhotons); jj++)
        {
            // printf("In first loop. Skip is %d\n", skip);
            // printf("jj is: %d\n", jj);
            // Calculate the time difference  
            photonDiff = photonTimes[jj] - photonTimes[ii];

            // Otherwise, calculate the frequency bin
            freqBin = floor(photonDiff / timeResol);

            // Store the data in the output
            outHistogram[freqBin] += photonWeights[ii] * photonWeights[jj];
        }

        // Loop through additional photons.
        for (int jj = ii + skip; jj < nPhotons; jj++)
        {
            // printf("In second loop. Skip is %d\n", skip);
            // printf("jj is %d\n", jj);

            // Calculate the time difference  
            photonDiff = photonTimes[jj] - photonTimes[ii];
            // printf("photonDiff is: %f\n", photonDiff);

            // Exit this second loop if the time difference is too large
            if (photonDiff >= windowSize)
            {
                // printf("Exiting. Photon diff is %d, window size is %d.\n", photonDiff, windowSize);
                skip = jj - ii - 1;
                //printf("skip at exit is %d\n", skip);
                break;
            }

            // Otherwise, calculate the frequency bin
            freqBin = floor(photonDiff / timeResol);

            // Store the data in the output
            outHistogram[freqBin] += photonWeights[ii] * photonWeights[jj];

        }
    }

    return outHistogram;

}

// // Calculates a histogram of time differences... using multiple cores.
// double * timeDifference_multi2(double *photonTimes, double *photonWeights,
//                               int windowSize, int maxFreq, int nPhotons,
//                               int nCores, int nSem)
// {

//     double *outHistogram;
//     int nbins = fftSize(windowSize, maxFreq);
//     float timeResol;
//     double photonDiff;
//     int freqBin;
//     int skip=1;
//     int breakArray[nCores];
//     int unitLength;
//     int pid = 0;
//     int start;
//     int stop;
//     int semCheck;
//     int bufferTime[buffer];
//     int bufferWeight[buffer];
//     sem_t* sem = mmap(NULL, sizeof(sem_t)*nSem,
//                       PROT_READ|PROT_WRITE, 
//                       MAP_SHARED|MAP_ANONYMOUS,-1, 0);

//     // Initalize the semaphore
//     for (int ii = 0; ii < nSem; ii++){
//             semCheck = sem_init(&sem[ii], 1, 1);
//     }

//     // Check to make sure the semaphore was created successfully.
//     if (semCheck != 0){
//         printf("Creating of semaphore failed.");
//         exit(0);
//     } 

//     // Calculate the unit length
//     unitLength = (int)floor(nPhotons / nCores);


//     //Calculate the time resolution
//     timeResol = 0.5 / maxFreq;

//     //calculate the size of the fft
//     //nbins = fftSize(windowSize, maxFreq);

//     // Create an array where we will store the output
//     // Actually, it is probably better if I let python handle the creation of
//     // this, and I just take the pointer as an input
//     outHistogram = mmap(0, nbins * sizeof(double), PROT_READ | PROT_WRITE,
//                         MAP_ANONYMOUS | MAP_SHARED, -1, 0);
//     if (outHistogram == MAP_FAILED) {perror("mmap");exit(1);}

//     // if memory cannot be allocated
//     if(outHistogram == NULL)                     
//     {
//         printf("Creation of outHistogram didn't work.");
//         exit(0);
//     }

//     // What if different cores try to write to the same outHist bin at the
//     // same time? Bad things will happen. There are way too many bins to lock
//     // each one with a semaphore. I can either use a single semaphore to lock
//     // a range of output bins, or I can store what I want to write to
//     // outHistogram in a buffer that gets read out at once. So, a core will
//     // lock outHistogram, write its buffer, then unlock out histogram. It also
//     // seems possible to do a mixture of the two: sort the buffer, unlock a
//     // portion of outHistogram, write, unlock, and move onto the next segment.

//     // Calculate the start and stop points for splitting the photon list
//     // between processors.
//     // The points are stored in breakArray
//     // Each process N will start at breakArray[N] and run to 
//     // (breakArray[N+1] - 1)
//     // unitLength is floor(nPhotons / nCores), and shows how many photons each
//     // process will iterate through. "floor" is needed because you can't run
//     // through a fractional number of photons, e.g. photons 0-3.333. 
//     // Note that breakArray is of length nCores + 1, and the below function
//     // only loops throgh nCores. The last index is automatically filled with
//     // the length of the photon array. In such a manner, the last process will
//     // pickup whatever photons are remaining.


//     // Loop through each core
//     for (int ii = 0; ii < nCores; ii++){
//         // Calculate the break points
//         breakArray[ii] = ii * unitLength;
//     }

//     // Set the last break point to ensure all photons are looked at
//     breakArray[nCores] = nPhotons;

//     // Begin setting up each individual process
//     for (int proc = 0; proc < nCores; proc++){

//         // Make sure the process can be created successfully
//         if ((pid=fork()) < 0){
//             printf("There was a fork error.");
//             perror("The fork did not complete successfully.");
//             exit(1);
//         }
//         //printf("The PID is %d\n", pid);
//         // If the process was created successfully...
//         if (pid == 0){
//             // The start index is the first value of the break array
//             start = breakArray[proc];

//             // The stop index is the next value, minus 1
//             stop  = breakArray[proc + 1] - 1;

//             // Loop through all the photons
//             for (int ii = start; ii <= stop; ii++)
//             {

//                 // If photon 0 differenced up to photon 100 before hitting the
//                 // window size, then I know for sure that photon 1 will also
//                 // difference up to photon 100 before hitting the window size.
//                 // Thus, I can skip checking photonDiff >= windowSize for the 
//                 // first 99 photons of the photon 1 time differencing.
//                 // "skip" impliments this skipping.

//                 // Quickly loop through the photons that I already know are okay.
//                 for (int jj = ii + 1; jj < (int)fmin(ii + skip, nPhotons); jj++)
//                 {
//                     // Calculate the time difference  
//                     photonDiff = photonTimes[jj] - photonTimes[ii];

//                     // Otherwise, calculate the frequency bin
//                     freqBin = floor(photonDiff / timeResol);

//                     // Store the data

//                     // Lock the semaphore
//                     semCheck = sem_wait(&sem[freqBin%nSem]);

//                     // Add the element directly
//                     outHistogram[freqBin] += photonWeights[ii] *
//                                              photonWeights[jj];

//                     sem_post(&sem[freqBin%nSem]);
//                 }

//                 // Loop through additional photons.
//                 for (int jj = ii + skip; jj <= nPhotons; jj++)
//                 {
//                     // printf("jj is %d\n", jj);

//                     // Calculate the time difference  
//                     photonDiff = photonTimes[jj] - photonTimes[ii];
//                     // printf("photonDiff is: %f\n", photonDiff);

//                     // Exit this second loop if the time difference is too large
//                     if (photonDiff >= windowSize)
//                     {
//                         // printf("Exiting. Photon diff is %d, window size is %d.\n", photonDiff, windowSize);
//                         skip = jj - ii - 1;
//                         //printf("skip at exit is %d\n", skip);

//                         break;
//                     }

//                     // Otherwise, calculate the frequency bin
//                     freqBin = floor(photonDiff / timeResol);

//                     // Store the data

//                     // Try and lock the semaphore
//                     semCheck = sem_wait(&sem[freqBin%nSem]);

//                     // Add the element directly
//                     outHistogram[freqBin] += photonWeights[ii] *
//                                              photonWeights[jj];

//                     sem_post(&sem[freqBin%nSem]);
//                 }

//            }
//         _exit(0);
//         }
//     while(wait(NULL) != -1);
//     }

//     // remove the semaphore
//     for (int ii = 0; ii < nSem; ii++){
//         semCheck = sem_destroy(&sem[ii]);
//     }
//     return outHistogram;
// }

// // Calculates a histogram of time differences... using multiple cores.
// double * timeDifference_multi(double *photonTimes, double *photonWeights,
//                               int windowSize, int maxFreq, int nPhotons,
//                               int nCores, int buffer)
// {

//     double *outHistogram;
//     int nbins = fftSize(windowSize, maxFreq);
//     float timeResol;
//     double photonDiff;
//     int freqBin;
//     int skip=1;
//     int breakArray[nCores];
//     int unitLength;
//     int pid = 0;
//     int start;
//     int stop;
//     int bufferTime[buffer];
//     double bufferWeight[buffer];
//     int inBuffer = 0;
//     int semCheck;
//     sem_t sem;

//     // Initalize the semaphore
//     semCheck = sem_init(&sem, 1, 1);

//     // Check to make sure the semaphore was created successfully.
//     if (semCheck != 0){
//         printf("Creating of semaphore failed.");
//         exit(0);
//     } 

//     // Calculate the unit length
//     unitLength = (int)floor(nPhotons / nCores);


//     //Calculate the time resolution
//     timeResol = 0.5 / maxFreq;

//     //calculate the size of the fft
//     //nbins = fftSize(windowSize, maxFreq);

//     // Create an array where we will store the output
//     // Actually, it is probably better if I let python handle the creation of
//     // this, and I just take the pointer as an input
//     outHistogram = mmap(0, nbins * sizeof(double), PROT_READ | PROT_WRITE,
//                         MAP_ANONYMOUS | MAP_SHARED, -1, 0);
//     if (outHistogram == MAP_FAILED) {perror("mmap");exit(1);}

//     // if memory cannot be allocated
//     if(outHistogram == NULL)                     
//     {
//         printf("Creation of outHistogram didn't work.");
//         exit(0);
//     }

//     // What if different cores try to write to the same outHist bin at the
//     // same time? Bad things will happen. There are way too many bins to lock
//     // each one with a semaphore. I can either use a single semaphore to lock
//     // a range of output bins, or I can store what I want to write to
//     // outHistogram in a buffer that gets read out at once. So, a core will
//     // lock outHistogram, write its buffer, then unlock out histogram. It also
//     // seems possible to do a mixture of the two: sort the buffer, unlock a
//     // portion of outHistogram, write, unlock, and move onto the next segment.

//     // Calculate the start and stop points for splitting the photon list
//     // between processors.
//     // The points are stored in breakArray
//     // Each process N will start at breakArray[N] and run to 
//     // (breakArray[N+1] - 1)
//     // unitLength is floor(nPhotons / nCores), and shows how many photons each
//     // process will iterate through. "floor" is needed because you can't run
//     // through a fractional number of photons, e.g. photons 0-3.333. 
//     // Note that breakArray is of length nCores + 1, and the below function
//     // only loops throgh nCores. The last index is automatically filled with
//     // the length of the photon array. In such a manner, the last process will
//     // pickup whatever photons are remaining.


//     // Loop through each core
//     for (int ii = 0; ii < nCores; ii++){
//         // Calculate the break points
//         breakArray[ii] = ii * unitLength;
//     }

//     // Set the last break point to ensure all photons are looked at
//     breakArray[nCores] = nPhotons;

//     // Begin setting up each individual process
//     for (int proc = 0; proc < nCores; proc++){

//         // Make sure the process can be created successfully
//         if ((pid=fork()) < 0){
//             printf("There was a fork error.");
//             perror("The fork did not complete successfully.");
//             exit(1);
//         }
//         //printf("The PID is %d\n", pid);
//         // If the process was created successfully...
//         if (pid == 0){
//             // The start index is the first value of the break array
//             start = breakArray[proc];

//             // The stop index is the next value, minus 1
//             stop  = breakArray[proc + 1] - 1;

//             // Loop through all the photons
//             for (int ii = start; ii <= stop; ii++)
//             {

//                 // If photon 0 differenced up to photon 100 before hitting the
//                 // window size, then I know for sure that photon 1 will also
//                 // difference up to photon 100 before hitting the window size.
//                 // Thus, I can skip checking photonDiff >= windowSize for the 
//                 // first 99 photons of the photon 1 time differencing.
//                 // "skip" impliments this skipping.

//                 // Quickly loop through the photons that I already know are okay.
//                 for (int jj = ii + 1; jj < (int)fmin(ii + skip, nPhotons); jj++)
//                 {
//                     // Calculate the time difference  
//                     photonDiff = photonTimes[jj] - photonTimes[ii];

//                     // Otherwise, calculate the frequency bin
//                     freqBin = floor(photonDiff / timeResol);

//                     // Store the data

//                     // Try and lock the semaphore
//                     semCheck = sem_trywait(&sem);

//                     // If we did not successfully lock the semaphore, and there
//                     // is still space in our buffer, just add the photon to the
//                     // buffer
//                     if (semCheck == -1 && inBuffer <= buffer - 1){

//                         // Store the frequency bin and associated weight
//                         // in the buffer
//                         bufferTime[inBuffer] = freqBin;
//                         bufferWeight[inBuffer] = photonWeights[ii] * photonWeights[jj];

//                         // Keep track of how many elements we have in the
//                         // buffer. 
//                         inBuffer += 1;

//                     } 

//                     // If the output is not locked, add the element directly
//                     // to the output. Also dump the buffer.
//                     // Additionally, if the output is locked, but the buffer
//                     // is full, we are going to have to wait in line and 
//                     // do this anyways.  
//                     else if (semCheck == 0)
//                     {

//                         // Add the element directly
//                         outHistogram[freqBin] += photonWeights[ii] *
//                                                  photonWeights[jj];

//                         // Dump the buffer
//                         for (int kk = 0; kk < inBuffer; kk++){
//                             outHistogram[bufferTime[kk]] = bufferWeight[kk];

//                             // Reset the in buffer counter to zero
//                             inBuffer = 0;
//                         }

//                         // Set the semaphore to unlock
//                         sem_post(&sem);

//                     }

//                     // If we get here, it means that the outHistogram is 
//                     // locked and our buffer is full. We have no choice but
//                     // to wait and write when possible. 
//                     else if (semCheck == -1 && inBuffer > buffer - 1)
//                     {
//                         sem_wait(&sem);

//                         // Add the element directly
//                         outHistogram[freqBin] += photonWeights[ii] *
//                                                  photonWeights[jj];

//                         // Dump the buffer
//                         for (int kk = 0; kk < inBuffer; kk++){
//                             outHistogram[bufferTime[kk]] = bufferWeight[kk];

//                             // Reset the in buffer counter to zero
//                             inBuffer = 0;
//                         }

//                         // Set the semaphore to unlock
//                         sem_post(&sem);
//                     }
//                 }

//                 // Loop through additional photons.
//                 for (int jj = ii + skip; jj <= nPhotons; jj++)
//                 {
//                     // printf("jj is %d\n", jj);

//                     // Calculate the time difference  
//                     photonDiff = photonTimes[jj] - photonTimes[ii];
//                     // printf("photonDiff is: %f\n", photonDiff);

//                     // Exit this second loop if the time difference is too large
//                     if (photonDiff >= windowSize)
//                     {
//                         // printf("Exiting. Photon diff is %d, window size is %d.\n", photonDiff, windowSize);
//                         skip = jj - ii - 1;
//                         //printf("skip at exit is %d\n", skip);

//                         // Dump the buffer
//                         // Actually... I can keep adding to the buffer just
//                         // fine, right?
//                         // semCheck = sem_wait(&sem);
//                         // for (int kk = 0; kk < inBuffer; kk++){
//                         //     outHistogram[bufferTime[kk]] = bufferWeight[kk];

//                         //     // Reset the in buffer counter to zero
//                         //     inBuffer = 0;
//                         // }
//                         // semCheck = sem_post(&sem);

//                         break;
//                     }

//                     // Otherwise, calculate the frequency bin
//                     freqBin = floor(photonDiff / timeResol);

//                     // Store the data

//                     // Try and lock the semaphore
//                     semCheck = sem_trywait(&sem);

//                     // If outHist is currently locked...
//                     // And we have not yet maxed out our buffer...
//                     if (semCheck == -1 && inBuffer <= buffer - 1){

//                         // Store the frequency bin and associated weight
//                         // in the buffer
//                         bufferTime[inBuffer] = freqBin;
//                         bufferWeight[inBuffer] = photonWeights[ii] * photonWeights[jj];

//                         // Keep track of how many elements we have in the
//                         // buffer. 
//                         inBuffer += 1;

//                     } 

//                     // If the output is not locked, add the element directly
//                     // to the output. Also dump the buffer.
//                     // Additionally, if the output is locked, but the buffer
//                     // is full, we are going to have to wait in line and 
//                     // do this anyways.  
//                     else if (semCheck == 0)
//                     {                        

//                         // Add the element directly
//                         outHistogram[freqBin] += photonWeights[ii] *
//                                                  photonWeights[jj];

//                         // Dump the buffer
//                         for (int kk = 0; kk < inBuffer; kk++){
//                             outHistogram[bufferTime[kk]] = bufferWeight[kk];

//                             // Reset the in buffer counter to zero
//                             inBuffer = 0;
//                         }

//                         // Set the semaphore to unlock
//                         sem_post(&sem);
//                     }
//                     else if (semCheck == -1 && inBuffer > buffer - 1)
//                     {
//                     sem_wait(&sem);

//                     // Add the element directly
//                     outHistogram[freqBin] += photonWeights[ii] *
//                                              photonWeights[jj];

//                     // Dump the buffer
//                     for (int kk = 0; kk < inBuffer; kk++){
//                         outHistogram[bufferTime[kk]] = bufferWeight[kk];

//                         // Reset the in buffer counter to zero
//                         inBuffer = 0;
//                     }

//                     // Set the semaphore to unlock
//                     sem_post(&sem);
//                     } 
//                 }
//             // Dump the buffer
//             semCheck = sem_wait(&sem);
//             for (int kk = 0; kk < inBuffer; kk++){
//                 outHistogram[bufferTime[kk]] = bufferWeight[kk];

//                 // Reset the in buffer counter to zero
//                 inBuffer = 0;
//             }
//             semCheck = sem_post(&sem);

//            }
//         _exit(0);
//         }
//     while(wait(NULL) != -1);
//     }

//     // remove the semaphore
//     semCheck = sem_destroy(&sem);
//     return outHistogram;
// }



void main(void){}

// void main(void){
//     double testArray[] = {1, 1.01, 1.02, 1.5, 2, 3, 4, 5, 6};
//     double testWeights[] = {1, 1, 1, 1, 2, 3, 4, 5, 6};
//     int windowSize = 2;
//     double maxFreq = 1;
//     int nphotons = 9;
//     int size;
//     double *histogram;
//     double *histogram2;

//     printf("Working on histogram 1\n");
//     histogram = timeDifference(testArray, testWeights, 
//                                windowSize, maxFreq, nphotons);

//     printf("Working on histogram 2\n");
//     histogram2 = timeDifference_multi(testArray, testWeights, 
//                                     windowSize, maxFreq, nphotons, 4, 125000);

//     size = fftSize(windowSize, maxFreq);

//     for (int ii=0; ii < size; ii++)
//     {
//         printf("%f vs %f\n", histogram[ii], histogram2[ii]);
//     }

//     return;
// }


// void main(void){
//     double testX[] = {1, 2, 3, 4, 5, 6};
//     int length = 5;
//     double * output;

//     output = linearFit(testX, length);

//     return;
// }