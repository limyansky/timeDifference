/******************************************************************************
 * file: power_vector.h
 * version: v4B_black_widows
 * date: 4-12-2013
 * author: mario (mario@piffio.org)
 * description:
 *   implements the class PowerVector that stores the power and interprets it
 * changelog:
 *   v4B_bw: alternative B in the conversion of the power into significance
 ******************************************************************************/
#ifndef BLINDSEARCH_power_vector_H
#define BLINDSEARCH_power_vector_H

class PowerVector {
  public:
    PowerVector(unsigned long size, int harm_num);
    PowerVector(void);
    void Clear(void);
    ~PowerVector(void);

    void SetFreqRange(const double min_f0, const double max_f0);
    void Set(unsigned long index, float value);
    float Get(unsigned long index);
    double GetFreq(unsigned long index);
    float GetProb(unsigned long index);
    float GetProb(float power);
    float GetPower(float significance);
    float GetSignif(unsigned long index);
    float GetSignif(float power);
    unsigned long MinIndex(void);
    unsigned long MaxIndex(void);
    int HarmNum(void);
    
    // NEW FIT USING HISTOGRAM + PDF
    void GetExponentialTail(const int tail_size, unsigned long *tail_index);
    void FitExponentialTail_without_ROOT(const int tail_size,
        unsigned long *tail_index);

    int GetNCandidates(const float sig_threshold,
        const int min_num, const int max_num,
        const int tail_size, unsigned long *index_tail);
    void FillResults(const int ncandidates,
        const int tail_size, unsigned long *index_tail,
        const float f1_f0, PulsarCandidate *results);
  private:
    float InsertSorted(const int array_size, unsigned long index,
        unsigned long *index_list);

    double size_, min_f0_, max_f0_, harm_num_;
    unsigned long min_index_, max_index_;
    float constant_, slope_, start_, end_;
    float *data_;
};
#endif // BLINDSEARCH_power_vector_H
