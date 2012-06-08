#include <divsim/stdinc.h>
#include <divsim/KDTreeModel.h>

#include <mathic/KDTree.h>
#include <mathic/DivList.h>

#include <gtest/gtest.h>

TEST(DivFinder, NoOp) {
  KDTreeModel<1,1,1,1,1> model(1, 1, 0, 0, 1.0, 1000);
};
