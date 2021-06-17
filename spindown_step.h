/******************************************************************************
 * file: spindown_step.h
 * version: v2_black_widows
 * date: 11-11-2013
 * author: mario (mario@piffio.org)
 * description:
 *   describe the class Step that describes the search settings in the
 *   spindown and harmonics combinations parameters grid
 * changelog:
 *   v0_bw: got rid of f2, unrelevant parameter for black widows
 *   v2_bw: renamed, dropped identifier, turned into a linked list
 ******************************************************************************/
#ifndef BLINDSEARCH_spindown_step_H
#define BLINDSEARCH_spindown_step_H

class SpindownStep {
  public:
    SpindownStep(const int step_num, const double start,
        const double step_size, const int N_steps, const int chat_level);
    SpindownStep(void);

    SpindownStep *Next(void);
    void Open(void);
    void SetHarmonic(const int harm_num);
    void Close(void);

    int Num(void);
    int HarmNum(void);
    double F1F0(void);
  private:
    int harm_, index_, N_steps_;
    SpindownStep *next_;
    double f1_f0_;
    int debug_;
    BenchmarkTool *timer_;
};

#endif // BLINDSEARCH_spindown_step_H
