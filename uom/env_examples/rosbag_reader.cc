// Copyright @2023 UOM project. All rights reserved.

#include "glog/logging.h"
#include "fmt/printf.h"
#include "gflags/gflags.h"
#include "rosbag/bag.h"
#include "rosbag/view.h"

DEFINE_string(input, "", "path to the rosbag file.");

int main(int argc, char* argv[]) {
  google::ParseCommandLineFlags(&argc, &argv, true);
  rosbag::Bag bag(FLAGS_input);
  rosbag::View view(bag);
  for (const auto& connect : view.getConnections()) {
     LOG(INFO) << fmt::format("topic: {}, datatype: {}\n", connect->topic, connect->datatype);
  }
}
