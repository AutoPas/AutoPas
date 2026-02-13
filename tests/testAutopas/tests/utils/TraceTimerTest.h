//
// Created by mmueh on 2/12/2026.
//
#pragma once

#include "AutoPasTestBase.h"
#include "autopas/utils/TraceTimer.h"

class TraceTimerTest : public AutoPasTestBase {
 protected:
  autopas::utils::TraceTimer _timer{};
};