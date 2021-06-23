#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#include "mkl.h"

double * timeDifference(double *, double*, int, int, int);
double * timeDifference_fast(double *, double*, int, int, int);
int fftSize(int, int);
// double * linearFit(double *, int);
void wrap_qsort(float *, int);
int cmpfunc(const void *, const void *);


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
double * timeDifference_fast(double *photonTimes, double *photonWeights,
                        int windowSize, int maxFreq, int nPhotons){

    double *outHistogram;
    int nbins;
    float timeResol;
    double photonDiff;
    int freqBin;
    int skip=1;
    printf("test print.");

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

//     histogram = timeDifference(testArray, testWeights, 
//                                windowSize, maxFreq, nphotons);

//     histogram2 = timeDifference_fast(testArray, testWeights, 
//                                     windowSize, maxFreq, nphotons);

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