// RaPlace-Light: real-time radar place recognition (RaPlace, IROS 2023,
// lightweight C++ port).
//
// Pipeline per scan:  polar (range x azimuth, uint8)
//   -> optional range crop + noise-floor subtraction
//   -> polar-to-cartesian remap at low resolution (cached LUT)
//   -> Radon transform (rotate + column sum) at n_theta angles
//   -> area-resize to n_theta_out columns (low-pass over theta)
//   -> |FFT| along rho, keep half spectrum, z-normalize  => descriptor
//
// Matching: circular cross-correlation along theta in the frequency
// domain (rotation invariant), distance = (autocorr - best corr) / 1000.

#pragma once

#include <opencv2/opencv.hpp>

#include <vector>

namespace raplace {

struct Params {
  int work_size = 69;       // cartesian image size (pixels)
  int n_theta = 180;        // radon angles computed
  int n_theta_out = 18;     // descriptor theta columns
  // crop_frac 0.5 + denoise 20 improves F1 for tight (20 m) revisit
  // criteria; raw (1.0 / 0) is better at 100 m scale. Default: raw.
  double crop_frac = 1.0;   // fraction of range bins kept
  int denoise = 0;          // subtracted from polar intensities
  int exclude_recent = 50;  // frames excluded around the query
  // Minimum scan signal fraction (signalFraction) for a frame to take
  // part in matching; gated-out frames stay in the DB (index-aligned)
  // but are never candidates, and a gated-out query yields no decision.
  // 0.0 disables. Recommended ~0.005 for marine sequences.
  double signal_gate = 0.0;
};

// Fraction of strong-return pixels — information score of a polar scan.
double signalFraction(const cv::Mat& polar, int intensity_th = 60);

class Descriptor {
 public:
  explicit Descriptor(const Params& p = Params());

  // polar: CV_8UC1, rows = range bins, cols = azimuth.
  // returns CV_32FC1 descriptor (rho_half rows x n_theta_out cols).
  cv::Mat fromPolar(const cv::Mat& polar);

 private:
  cv::Mat fromCart(const cv::Mat& cart);
  void buildMaps(int n_range, int n_az, int cart_size);

  Params p_;
  cv::Mat map_x_, map_y_;          // polar->cart LUT
  int map_range_ = -1, map_az_ = -1;
  std::vector<cv::Mat> rot_mats_;  // per-angle 2x3 rotation matrices
  int diag_ = 0;
};

class Matcher {
 public:
  explicit Matcher(const Params& p = Params());

  // signal: the scan's signalFraction; only used when signal_gate > 0.
  void add(const cv::Mat& desc, double signal = 1.0);

  // Search valid frames older than exclude_recent. Returns best index and
  // the RaPlace distance via out params; false if the database is too
  // small or the query/candidates are gated out.
  bool query(const cv::Mat& desc, int* best_idx, double* min_dist,
             double signal = 1.0) const;

  int size() const { return static_cast<int>(db_.size()); }

 private:
  cv::Mat thetaFFT(const cv::Mat& desc) const;  // CV_32FC2, per-row FFT

  Params p_;
  std::vector<cv::Mat> db_;  // theta-axis FFT of each descriptor
  std::vector<bool> valid_;  // passed the signal gate at add() time
};

}  // namespace raplace
