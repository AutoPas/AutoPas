#ifndef KMP_STATS_TIMING_H
#define KMP_STATS_TIMING_H

/** @file kmp_stats_timing.h
 * Access to real time clock and timers.
 */

//===----------------------------------------------------------------------===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is dual licensed under the MIT and the University of Illinois Open
// Source Licenses. See LICENSE.txt for details.
//
//===----------------------------------------------------------------------===//

#include "kmp_os.h"
#include <limits>
#include <stdint.h>
#include <string>
#if KMP_HAVE_X86INTRIN_H
#include <x86intrin.h>
#endif

class tsc_tick_count {
private:
  int64_t my_count;

public:
  class tsc_interval_t {
    int64_t value;
    explicit tsc_interval_t(int64_t _value) : value(_value) {}

  public:
    tsc_interval_t() : value(0) {} // Construct 0 time duration
#if KMP_HAVE_TICK_TIME
    double seconds() const; // Return the length of a time interval in seconds
#endif
    double ticks() const { return double(value); }
    int64_t getValue() const { return value; }
    tsc_interval_t &operator=(int64_t nvalue) {
      value = nvalue;
      return *this;
    }

    friend class tsc_tick_count;

    friend tsc_interval_t operator-(const tsc_tick_count &t1,
                                    const tsc_tick_count &t0);
    friend tsc_interval_t operator-(const tsc_tick_count::tsc_interval_t &i1,
                                    const tsc_tick_count::tsc_interval_t &i0);
    friend tsc_interval_t &operator+=(tsc_tick_count::tsc_interval_t &i1,
                                      const tsc_tick_count::tsc_interval_t &i0);
  };

#if KMP_HAVE___BUILTIN_READCYCLECOUNTER
  tsc_tick_count()
      : my_count(static_cast<int64_t>(__builtin_readcyclecounter())) {}
#elif KMP_HAVE___RDTSC
  unsigned int dummy;
  tsc_tick_count() : my_count(static_cast<int64_t>(__rdtscp(&dummy))) {}  //use rdtscp instead of rdtsc  AUTO by Ali
#else
#error Must have high resolution timer defined
#endif
  tsc_tick_count(int64_t value) : my_count(value) {}
  int64_t getValue() const { return my_count; }
  tsc_tick_count later(tsc_tick_count const other) const {
    return my_count > other.my_count ? (*this) : other;
  }
  tsc_tick_count earlier(tsc_tick_count const other) const {
    return my_count < other.my_count ? (*this) : other;
  }
#if KMP_HAVE_TICK_TIME
  static double tick_time(); // returns seconds per cycle (period) of clock
#endif
  static tsc_tick_count now() {
    return tsc_tick_count();
  } // returns the rdtsc register value
  friend tsc_tick_count::tsc_interval_t operator-(const tsc_tick_count &t1,
                                                  const tsc_tick_count &t0);
};

inline tsc_tick_count::tsc_interval_t operator-(const tsc_tick_count &t1,
                                                const tsc_tick_count &t0) {
  return tsc_tick_count::tsc_interval_t(t1.my_count - t0.my_count);
}

// added here the tick time definition to enable using it without enabling LIBOMP_STATS by Ali
#if KMP_HAVE_TICK_TIME
#if KMP_MIC
double tsc_tick_count::tick_time() {
  // pretty bad assumption of 1GHz clock for MIC
  return 1 / ((double)1000 * 1.e6);
}
#elif KMP_ARCH_X86 || KMP_ARCH_X86_64
#include <string.h>
// Extract the value from the CPUID information
double tsc_tick_count::tick_time() {
  static double result = 0.0;

  if (result == 0.0) {
    kmp_cpuid_t cpuinfo;
    char brand[256];

    __kmp_x86_cpuid(0x80000000, 0, &cpuinfo);
    memset(brand, 0, sizeof(brand));
    int ids = cpuinfo.eax;

    for (unsigned int i = 2; i < (ids ^ 0x80000000) + 2; i++)
      __kmp_x86_cpuid(i | 0x80000000, 0,
                      (kmp_cpuid_t *)(brand + (i - 2) * sizeof(kmp_cpuid_t)));

    char *start = &brand[0];
    for (; *start == ' '; start++)
      ;

    char *end = brand + KMP_STRLEN(brand) - 3;
    uint64_t multiplier;

    if (*end == 'M')
      multiplier = 1000LL * 1000LL;
    else if (*end == 'G')
      multiplier = 1000LL * 1000LL * 1000LL;
    else if (*end == 'T')
      multiplier = 1000LL * 1000LL * 1000LL * 1000LL;
    else {
      std::cout << "Error determining multiplier '" << *end << "'\n";
      exit(-1);
    }
    *end = 0;
    while (*end != ' ')
      end--;
    end++;

    double freq = strtod(end, &start);
    if (freq == 0.0) {
      std::cout << "Error calculating frequency " << end << "\n";
      exit(-1);
    }

    result = ((double)1.0) / (freq * multiplier);
  }
  return result;
}
#endif
#endif

inline tsc_tick_count::tsc_interval_t
operator-(const tsc_tick_count::tsc_interval_t &i1,
          const tsc_tick_count::tsc_interval_t &i0) {
  return tsc_tick_count::tsc_interval_t(i1.value - i0.value);
}

inline tsc_tick_count::tsc_interval_t &
operator+=(tsc_tick_count::tsc_interval_t &i1,
           const tsc_tick_count::tsc_interval_t &i0) {
  i1.value += i0.value;
  return i1;
}

#if KMP_HAVE_TICK_TIME
inline double tsc_tick_count::tsc_interval_t::seconds() const {
  return value * tick_time();
}
#endif

extern std::string formatSI(double interval, int width, char unit);

inline std::string formatSeconds(double interval, int width) {
  return formatSI(interval, width, 'S');
}

inline std::string formatTicks(double interval, int width) {
  return formatSI(interval, width, 'T');
}

#endif // KMP_STATS_TIMING_H
