// Streaming demo: runs RaPlace-Light over a directory of polar PNGs,
// reporting loop-closure candidates and per-stage timing.
//
//   ./raplace_demo <polar_dir> [max_frames]

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <numeric>
#include <string>
#include <vector>

#include "raplace.hpp"

namespace fs = std::filesystem;
using clk = std::chrono::steady_clock;

static double ms(clk::time_point a, clk::time_point b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <polar_dir> [max_frames]\n", argv[0]);
    return 1;
  }
  const std::string dir = argv[1];
  const size_t max_frames = argc > 2 ? std::stoul(argv[2]) : SIZE_MAX;

  std::vector<std::string> files;
  for (const auto& e : fs::directory_iterator(dir)) {
    if (e.path().extension() == ".png") files.push_back(e.path().string());
  }
  std::sort(files.begin(), files.end());
  if (files.size() > max_frames) files.resize(max_frames);
  std::printf("%zu frames\n", files.size());

  raplace::Params p;
  raplace::Descriptor desc_gen(p);
  raplace::Matcher matcher(p);

  std::vector<double> t_desc, t_query;
  for (size_t i = 0; i < files.size(); ++i) {
    cv::Mat polar = cv::imread(files[i], cv::IMREAD_GRAYSCALE);
    auto t0 = clk::now();
    cv::Mat d = desc_gen.fromPolar(polar);
    auto t1 = clk::now();
    int best_idx = -1;
    double min_dist = 0;
    bool ok = matcher.query(d, &best_idx, &min_dist);
    auto t2 = clk::now();
    matcher.add(d);
    t_desc.push_back(ms(t0, t1));
    t_query.push_back(ms(t1, t2));
    if (ok && (i + 1) % 500 == 0) {
      std::printf("frame %zu: best=%d dist=%.5f desc=%.2fms query=%.2fms\n",
                  i, best_idx, min_dist, t_desc.back(), t_query.back());
    }
  }
  auto mean = [](const std::vector<double>& v) {
    return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
  };
  std::printf("mean desc %.2f ms, mean query %.2f ms (db=%d)\n",
              mean(t_desc), mean(t_query), matcher.size());
  return 0;
}
