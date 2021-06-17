/******************************************************************************
 * benchmark_tool.h
 * defines a BenchmarkTool class to measure how long an operation takes
 * example:
 *
 * #include <my_function.h>
 * #include benchmark_tool.h
 *
 * int main(int argc, char* argv[]) {
 *   int chat_level = 2;
 *   int high_precision = 1;
 *   int severity = 3;
 *   BenchmarkTool b(high_precision,chat_level);
 *   b.Tic();
 *   my_function();
 *   b.Tac("running my_function", severity);
 * }
 * Author: Mario 20/11/10
 ******************************************************************************/
#ifndef BLINDSEARCH_benchmark_tool_H
#define BLINDSEARCH_benchmark_tool_H

class BenchmarkTool {
  public:
    BenchmarkTool(const int high_precision,const int chat_level);
/* chat_levels [default DEBUG]:
 * 0 [SILENT]  = Don't print anything
 * 1 [QUIET]   = Print only ERRORs
 * 2 [NORMAL]  = Print ERRORs and WARNINGs 
 * 3 [VERBOSE] = Print ERRORs, WARNINGs, and NOTICEs
 * 4 [DEBUG]   = Print everything
 */
    explicit BenchmarkTool(const int high_precision);
    BenchmarkTool(void);

    void Tic(void);
    double Tac(const char* operation_name, int severity);
/* severity levels [default ERROR]:
 * 0 [IGNORE]  = Never print
 * 1 [DEBUG]   = Print only if DEBUG is turned on
 * 2 [NOTICE]  = Print if VERBOSE or DEBUG are set
 * 3 [WARNING/IMPORTANT] = Print also for the NORMAL case
 * 4 [ERROR/FUNDAMENTAL] = Print even if the QUIET option is set
 */
    double Tac(const char* operation_name);
  private:
    const char* format;
    bool subsec;
    int severity_threshold;

    double start;
    clock_t cstart;
    time_t tstart;
    time_t tstop;
};

#endif // BLINDSEARCH_benchmark_tool_H
