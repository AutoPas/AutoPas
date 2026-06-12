#ifndef TIMINGSTATS_H
#define TIMINGSTATS_H

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

#endif  // TIMINGSTATS_H