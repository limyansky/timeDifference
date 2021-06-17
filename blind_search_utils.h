/******************************************************************************
 * file: blind_search_utils.h
 * version: v3_black_widows
 * date: 12-11-2013
 * author: mario (mario@piffio.org)
 * description:
 *   describe some of the core functions of the blind search program
 * changelog:
 *   v0_bw: moved functions around and dropped f2 support
 *   v2_bw: use spindown_step instead of step
 *   v3_bw: added orbital demodulation
 ******************************************************************************/
#ifndef BLINDSEARCH_blind_search_utils_H_
#define BLINDSEARCH_blind_search_utils_H_

std::vector<double> TimeDiff(
    std::vector<double> &photon_toas, std::vector<float> &photon_weights,
    double time_window_for_diffs, std::vector<float> &diff_weights);

std::vector<double> CorrectTimes(std::vector<double> &photon_times,
    double f1_f0);

int PickTheBest(PowerVector &PV, SpindownStep &f1f0_step,
    const int min_num, const int max_num, const float sig_threshold,
    PulsarCandidate *results);
void append_sorted(int array_size, unsigned long index, float value,
    unsigned long *index_list, float *value_list);

std::vector<double> RemoveCircularOrbitModulation(
    std::vector<double> &photon_toas, double epoch, 
    double PB, double A1, double T0);

#endif // BLINDSEARCH_blind_search_utils_H_
