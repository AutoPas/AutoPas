/**
 * @file CycleCounterQuery.cpp
 * @author MehdiHachicha
 * @date 06.05.2024
 * Script to check if system provides a cycle counter register.
 */

#include <iostream>

int main() {
  std::cout << __has_builtin(__builtin_readcyclecounter) << std::endl;
  return 0;
}
