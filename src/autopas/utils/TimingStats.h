#ifndef TIMINGSTATS_H
#define TIMINGSTATS_H

#include <cmath>
#include <limits>
#include <spdlog/spdlog.h>
#include <vector>
class TimingStats{
    public:
        TimingStats(const std::string& name) : operationName(name), sumTimings(0), minimumTiming(std::numeric_limits<double>::max()), maximumTiming(0), count(0) {}
        ~TimingStats() {
            printStats(operationName);
        }

        // no copying
        TimingStats(const TimingStats&) = delete;
        TimingStats& operator=(const TimingStats&) = delete;

        void addTiming(double timing) {
            sumTimings += timing;
            if (timing < minimumTiming) {
                minimumTiming = timing;
            }
            if (timing > maximumTiming) {
                maximumTiming = timing;
            }
            times.push_back(timing);
            ++count;
        }
        
        void printStats(const std::string& operationName) const {
            if (count > 0) {
                double averageTiming = sumTimings / count;
                double stdev = 0;
                for (const auto& time : times) {
                    stdev += (time - averageTiming) * (time - averageTiming);
                }
                stdev = std::sqrt(stdev / count);
                spdlog::info("{} - Average: {} s, Min: {} s, Max: {} s, StdDev: {} s over {} instances", operationName, averageTiming, minimumTiming, maximumTiming, stdev, count);
            } else {
                spdlog::info("{} - No timings recorded.", operationName);
            }
        }

    private:
        std::string operationName;
        std::vector<double> times;
        double sumTimings;
        double minimumTiming;
        double maximumTiming;
        std::size_t count;

};

struct KernelTimings{
    std::string _operationName{};
    TimingStats _preparation{"preparation before Kernel"};
    TimingStats _kernel{"actual Kernel/operation"};
    TimingStats _cleanup{"anything after Kernel"};
    TimingStats _total{"total"};

    KernelTimings(const std::string& name) : _operationName(name){}
    ~KernelTimings(){
        spdlog::info("\n\n\n Timing Stats of ({}) \n",_operationName);
    }
};

#endif  // TIMINGSTATS_H