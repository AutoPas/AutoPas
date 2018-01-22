//
// Created by seckler on 22.01.18.
//

#include "SPHTest.h"
#include "sph/autopassph.h"
#include "autopas.h"

TEST(SPHTest, testW) {
    double value = autopas::sph::W({1., 1., 1.}, 1.);
    double should_be_value = 0.00944773;
    EXPECT_NEAR(value, should_be_value, 1e-8);

    value = autopas::sph::W({1.,.5,.25},.5);
    should_be_value = 0.00151727;
    EXPECT_NEAR(value, should_be_value, 1e-8);
}

TEST(SPHTest, testGradW) {
    std::array<double,3> value = autopas::sph::gradW({1., 1., 1.}, 1.);
    std::array<double,3> should_be_value = {-0.0213086,   -0.0213086,    -0.0213086};
    for(int i=0;i<3;i++) {
        EXPECT_NEAR(value[i], should_be_value[i], 1e-7);
    }

    value = autopas::sph::gradW({1.,.5,.25},.5);
    should_be_value = {-0.038073,   -0.0190365,    -0.00951825};
    for(int i=0;i<3;i++) {
        EXPECT_NEAR(value[i], should_be_value[i], 1e-7);
    }
}
