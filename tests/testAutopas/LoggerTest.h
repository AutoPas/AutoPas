/**
 * @file LoggerTest.h
 * @author seckler
 * @date 16.04.18
 */

#pragma once

#include <ostream>

// Helper:

class ScopedRedirect {
 public:
  ScopedRedirect(std::ostream& inOriginal, std::ostream& inRedirect)
      : mOriginal(inOriginal), mRedirect(inRedirect) {
    mOriginal.rdbuf(mRedirect.rdbuf(mOriginal.rdbuf()));
  }

  ~ScopedRedirect() { mOriginal.rdbuf(mRedirect.rdbuf(mOriginal.rdbuf())); }

 private:
  ScopedRedirect(const ScopedRedirect&);
  ScopedRedirect& operator=(const ScopedRedirect&);

  std::ostream& mOriginal;
  std::ostream& mRedirect;
};