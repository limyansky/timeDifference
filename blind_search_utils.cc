/******************************************************************************
 * file: blind_search_utils.cc
 * version: v4B_black_widows
 * date: 4-12-2013
 * author: mario (mario@piffio.org)
 * description:
 *   here are some of the core functions of the blind search program:
 *   CorrectTimes() that compresses the TOAs time series to compensate for the
 *   pulsar spin-down, TimeDiff() that applies the time-differencing transform
 *   to the series of TOAs and generates that of time differences.
 * changelog:
 *   v0_bw: moved functions around and dropped f2 support
 *   v2_bw: use spindown_step instead of step
 *   v3_bw: added orbital demodulation
 *   v4B_bw: alternative B for picking the best candidates
 ******************************************************************************/
#include <stdio.h> // Basic I/O
#include <math.h> // sqrt
#include <time.h>
#include <sys/time.h> // it provides gettimeofday() on linux
#include <vector> // Vector utilities

#include "pulsar_candidate.h"
#include "benchmark_tool.h"
#include "spindown_step.h"
#include "power_vector.h"

#include "blind_search_utils.h"

/******************************************************************************
 * CorrectTimes()
 *   correct the photon times of arrival for F1 before running the FFT.
 *   In an "almost periodic" time series, the probability to get an event
 *   in a unit time has the property:
 * P(p(t + 2pi*T)) == P(p(t))  with  p(t) = F0*t + (F1*t^2)/2
 *   This is  expected from a millisecond pulsar, in its frame of reference.
 *   In order to take advantage of an FFT, the time series should be periodic:
 * P(F0*(t' + 2pi*T)) == P(F0*t')  or  p(t') = F0*t'
 *   This goal can be achieved by rescaling the time t into:
 * t'(t) = t + [(F1/F0)*t^2]/2
 *   This linearization does not require to know in advance the values of both
 *   F0 and F1, but only their ratio F1/F0. (young pulsars are more complicated)
 *   F0 can then be found through an FFT on the linearized time of arrivals.
 ******************************************************************************/
std::vector<double> CorrectTimes(std::vector<double> &photon_toas,
    double f1_f0) {
  std::vector<double> corrected_photon_times;

  unsigned int number_of_photons = static_cast<unsigned int>(
      photon_toas.size());
  double f1_corr = 0.0;
  double toa = 0.0;

  for (unsigned int phot = 0; phot < number_of_photons; phot++) {
    toa = photon_toas[phot];
    f1_corr = toa * toa * f1_f0 / 2.0;
    /* For millisecond pulsars higher order derivatives are superfluous */
    corrected_photon_times.push_back(toa + f1_corr);
  }
  return corrected_photon_times;
}

/******************************************************************************
 * TimeDiff()
 *   given the toas of the photons, extract the differences, up to a maximum
 *   value, time_window_for_diffs.
 ******************************************************************************/
std::vector<double> TimeDiff(
    std::vector<double> &photon_toas, std::vector<float> &photon_weights,
    double time_window_for_diffs, std::vector<float> &diff_weights) {
  std::vector<double> time_differences;

  unsigned int number_of_photons = static_cast<unsigned int>(
      photon_toas.size());
  /* difference between the time of arrivals of the two photons */
  double time_diff = 0.0;
  float diff_weight = 0.0;
  float total_weight = 0.0;
  float scale = 0.0;

  diff_weights.clear();
  unsigned int phot_2 = 0;
  for (unsigned int phot_1 = 0; phot_1 < number_of_photons - 1; phot_1++) {
    for (phot_2 = phot_1 + 1; phot_2 < number_of_photons; phot_2++) {
      time_diff = photon_toas[phot_2] - photon_toas[phot_1];
      if (time_diff > time_window_for_diffs) break;
      if (time_diff < 0.0) break;
      /* The weight of a difference is the Bayesian probability that it
       * belongs to the source. For this to happen both events must come
       * from the source so the probability of a time difference is the
       * composition, and thus the product, of the weights of the photons */
      diff_weight = photon_weights[phot_2] * photon_weights[phot_1];

      time_differences.push_back(time_diff);
      diff_weights.push_back(diff_weight);
      total_weight += diff_weight;
    }
  }
  if (total_weight > 0.0) {
    /* It is not obvious that this is the best choice for the normalization */
    scale = static_cast<float>(time_differences.size()) / total_weight;
    for (unsigned int diff = 0; diff < time_differences.size(); diff++)
      diff_weights[diff] *= scale;
  }
  return time_differences;
}

/******************************************************************************
 * PickTheBest_B()
 *   parse the power vector and select the best candidates, filling reults.
 *   The definition of best candidates rely on a figure of merit, in this case
 *   the raw power as it appears in the power vector entries.
 *   The selection of the best candidates passes through the extraction of a
 *   tail, that is a sorted array of the best candidates. The tail is then
 *   directly fit minimizing the chi2
 ******************************************************************************/
int PickTheBest_B(PowerVector &PV, const double f1_f0,
    const int min_num, const int max_num, const float sig_threshold,
    PulsarCandidate *results) {
  int ncandidates = 0;
  /* XXX the size of the tail is key to the performance and robustness 
   * of the algorithm that extracts the candidates.
   * the cpu complexity scales as (N_tail^2 + x*N_PV)*log(N_PV) 
   * where x is some positive constant. if N_tail=sqrt(x*N_PV)
   * the system is well balanced between the two parts. beyond
   * that the complexity becomes a burden. on the other side, the
   * longer the tail, the more robust the fit. the value of x
   * represents the ratio between the complexity of checking whether
   * a value should be considered for insertion, and the insertion
   * itself. it could be worked out analytically, but I rather get
   * this value from direct testing */
  int tail_size = 10000; // should be proportional to sqrt(len(PV))
  unsigned long index_tail[tail_size];
  for (int j=0; j<tail_size; j++)
    index_tail[j] = 0;

  // STEP 1: get a suitable tail out of the power vector (DONE)
  PV.GetExponentialTail(tail_size, index_tail);
  // STEP 2: bin the tail, fit, and get constant+slope (DONE)
  PV.FitExponentialTail_without_ROOT(tail_size, index_tail);
  // STEP 3: get the best candidates out of the tail itself (DONE)
  ncandidates = PV.GetNCandidates(sig_threshold, min_num, max_num,
      tail_size, index_tail);

  PV.FillResults(ncandidates, tail_size, index_tail, f1_f0, results);
  return ncandidates;
}

/******************************************************************************
 * PickTheBest()
 *   parse the power vector and select the best candidates, filling reults.
 *   The definition of best candidates rely on a figure of merit and the
 *   conversion of this figure of merit into significance can be performed
 *   in many different ways. Here is the implemented strategy:
 *   2) or B) extracts from the power vector a "tail" with only the best
 *      candidates and analizes that tail to extract the best ones and to
 *      assign them a significance
 ******************************************************************************/
int PickTheBest(PowerVector &PV, SpindownStep &f1f0_step,
    const int min_num, const int max_num, const float sig_threshold,
    PulsarCandidate *results) {
  int ncand = 0;

  double f1_f0 = f1f0_step.F1F0();
  ncand = PickTheBest_B(PV, f1_f0, min_num, max_num, sig_threshold, results);
  return ncand;
}

/******************************************************************************
 * RemoveCircularOrbitModulation()
 * compensate for the Roemer delay associated to a circular binary orbit.
 * in the Newtonian regime, the orbit is described by 3 parameters:
 * A1 (lt sec) semiaxis of the pulsar orbit projected along the line of site,
 * PB (days) orbital period, and T0 (MJD) epoch of ascending nodes.
 * the event times are coded as toas, in seconds from epoch (MET)
 ******************************************************************************/
std::vector<double> RemoveCircularOrbitModulation(
    std::vector<double> &photon_toas, double epoch, 
    double PB_d, double A1_lts, double T0_MJD) {
  std::vector<double> demodulated_toas;
  double orbital_phase = 0.0;
  double toa = 0.0;

  /* this is the definition of MET, Mission Elapsed Time */
  double T0_MET = (T0_MJD - 51910.0) * 86400.0 - 64.184;
  double Omega_Hz = M_PI / (43200.0 * PB_d);

  for (int i = 0; i < photon_toas.size(); i++) {
    toa = photon_toas[i];
    orbital_phase = Omega_Hz * (toa + epoch - T0_MET);
    demodulated_toas.push_back(toa - A1_lts * sin(orbital_phase));
  }
  return demodulated_toas;
}
