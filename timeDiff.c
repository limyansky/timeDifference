#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double * timeDifference(double *, double*, int, int, int);
double * timeDifference_fast(double *, double*, int, int, int);
int fftSize(int, int);

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