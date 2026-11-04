#pragma once

enum kmp_agent_type {
  Q_LEARNING_OLD = 0,
  SARSA_OLD = 1,
  Q_LEARNING = 2,
  DOUBLEQ = 3,
  SARSA = 4,
  EXPECTED_SARSA = 5,
  QV_LEARNING = 6,
  R_LEARNING = 7,
  DEEPQ_LEARNING = 8,
  CHUNK_LEARNING = 9
};