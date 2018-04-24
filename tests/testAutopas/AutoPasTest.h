/**
 * @file AutoPasTest.h
 * @author seckler
 * @date 24.04.18
 */

#include <gtest/gtest.h>
#include <utils/Logger.h>

#pragma once

class AutoPasTest : public testing::Test {
 public:
  AutoPasTest(){
    autopas::logger::create();
  }

  virtual ~AutoPasTest(){
    autopas::logger::unregister();
  }
};



