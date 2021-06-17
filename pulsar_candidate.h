/******************************************************************************
 * file: pulsar_candidate.h
 * version: v0_black_widows
 * date: 11-11-2013
 * author: mario (mario@piffio.org)
 * description:
 *   describe the class PulsarCandidate that collects the information about
 *   a pulsar (spin parameters) and its significance in the search
 * changelog:
 *   v0_bw: got rid of the refinement parameters
 ******************************************************************************/
#ifndef BLINDSEARCH_pulsar_candidate_H_
#define BLINDSEARCH_pulsar_candidate_H_

class PulsarCandidate {
  public:
    PulsarCandidate(double f0, double f1);
    PulsarCandidate(void);

    void SetFFTPower(float power, float significance);

    double FFTF0(void);
    double FFTF1(void);
    float FFTPower(void);
    float FFTSignif(void);
  private:
    // variables coming from the FFT spectrum analysis
    double FFT_f0_, FFT_f1_;
    float FFT_power_, FFT_significance_;
};

#endif // BLINDSEARCH_pulsar_candidate_H_
