//
// Created by Damitha on 3/21/2018.
//
#include "gtest/gtest.h"

TEST(basic_check, test_eq) {
    for(int i;i<4;i++){
        EXPECT_EQ(1, 1);
    }
    EXPECT_EQ(1, 1);
}

TEST(basic_check, test_neq) {
    EXPECT_EQ(1, 0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
