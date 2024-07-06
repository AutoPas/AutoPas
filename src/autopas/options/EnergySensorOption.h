/**
 * @file EnergySensorOption.h
 * @author Maximilian Praus
 * @date 06.07.2024
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {

/**
 * Class representing the different energy sensor options
 */

class EnergySensorOption : public Option<EnergySensorOption> {
    public:
    /**
     * Possible choices for energy sensor
     */
    enum Value {
        /**
         * No energy sensor. Used when energy measurement not available
         */
        none,
        /**
         * Use NVML
         */
        nvml,
        /**
         * Use NVIDIA
         */
        nvidia,
        /**
         * Use ROCM
         */
        rocm,
        /**
         * Use Cray
         */
        cray,
        /**
         * Use LIKWID
         */
        likwid,
        /**
         * Use RAPL
         */
        rapl,
        /**
         * Use TEGRA
         */
        tegra,
        /**
         * Use Xilinx
         */
        xilinx,
        /**
         * Use Powersensor2
         */
        powersensor2,
        /**
         * Use PowerSensor3
         */
        powersensor3,
    };

    /**
     * Cosntructor
     */
    EnergySensorOption() = default;

    /**
     * Cosntructor with selected option
     * @param option
     */
    constexpr EnergySensorOption(Value option) : _value(option) {}

    /**
     * Cast to value
     * @return
     */
    constexpr operator Value() const {return _value;}

    /**
     * Provide a way to iterate over the options of EnergySensorOption
     * @return map option -> string representation
     */
    static std::map<EnergySensorOption, std::string> getOptionNames() {
        return {
            {EnergySensorOption::none, "none"},
            {EnergySensorOption::nvml, "nvml"},
            {EnergySensorOption::nvidia, "nvidia"},
            {EnergySensorOption::rocm, "rocm"},
            {EnergySensorOption::cray, "cray"},
            {EnergySensorOption::likwid, "likwid"},
            {EnergySensorOption::rapl, "rapl"},
            {EnergySensorOption::tegra, "tgera"},
            {EnergySensorOption::xilinx, "xilinx"},
            {EnergySensorOption::powersensor2, "powersensor2"},
            {EnergySensorOption::powersensor3, "powersensor3"},
        };
    };


    private:
    Value _value{Value(-1)};
};

}