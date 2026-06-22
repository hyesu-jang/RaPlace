#include "raplace.hpp"

#include <cmath>
#include <mutex>

namespace raplace {

Descriptor::Descriptor(const Params& p) : p_(p) {
  diag_ = static_cast<int>(std::ceil(p_.work_size * std::sqrt(2.0))) + 2;
  cv::Point2f center((diag_ - 1) / 2.0f, (diag_ - 1) / 2.0f);
  rot_mats_.reserve(p_.n_theta);
  for (int i = 0; i < p_.n_theta; ++i) {
    double th = i * 180.0 / p_.n_theta;
    rot_mats_.push_back(cv::getRotationMatrix2D(center, th, 1.0));
  }
}

void Descriptor::buildMaps(int n_range, int n_az, int cart_size) {
  const int s = cart_size;
  map_x_.create(s, s, CV_32FC1);
  map_y_.create(s, s, CV_32FC1);
  const float c = (s - 1) / 2.0f;
  const float r_scale = n_range / (s / 2.0f);
  for (int y = 0; y < s; ++y) {
    float* mx = map_x_.ptr<float>(y);
    float* my = map_y_.ptr<float>(y);
    for (int x = 0; x < s; ++x) {
      float dx = x - c, dy = y - c;
      float az = std::atan2(dy, dx) / static_cast<float>(2.0 * CV_PI);
      az -= std::floor(az);  // [0,1)
      mx[x] = az * n_az;
      my[x] = std::hypot(dx, dy) * r_scale;
    }
  }
  map_range_ = n_range;
  map_az_ = n_az;
}

cv::Mat Descriptor::fromPolar(const cv::Mat& polar) {
  cv::Mat pol = polar;
  if (p_.crop_frac < 1.0) {
    int rows = static_cast<int>(polar.rows * p_.crop_frac);
    pol = polar.rowRange(0, rows);
  }
  if (p_.denoise > 0) {
    cv::Mat tmp;
    cv::subtract(pol, cv::Scalar(p_.denoise), tmp);  // saturates at 0
    pol = tmp;
  }
  // Render the cart image supersampled, then area-average down: the cart
  // grid point-samples polar, and thousands of range bins mapped onto
  // ~work_size/2 pixels would alias away thin radar returns (validated
  // against the dataset-provided cart images).
  const int kSupersample = 3;
  const int render = kSupersample * p_.work_size;
  const int pre_rows = 2 * render;
  if (pol.rows > 2 * pre_rows) {
    cv::Mat shrunk;
    cv::resize(pol, shrunk, cv::Size(pol.cols, pre_rows), 0, 0,
               cv::INTER_AREA);
    pol = shrunk;
  }
  if (pol.rows != map_range_ || pol.cols != map_az_) {
    buildMaps(pol.rows, pol.cols, render);
  }
  cv::Mat cart;
  cv::remap(pol, cart, map_x_, map_y_, cv::INTER_LINEAR,
            cv::BORDER_CONSTANT, 0);
  cv::resize(cart, cart, cv::Size(p_.work_size, p_.work_size), 0, 0,
             cv::INTER_AREA);
  return fromCart(cart);
}

cv::Mat Descriptor::fromCart(const cv::Mat& cart) {
  cv::Mat img;
  cart.convertTo(img, CV_32F);

  // radon: zero-pad so rotations keep the full image inside the frame
  cv::Mat padded = cv::Mat::zeros(diag_, diag_, CV_32F);
  int off_y = (diag_ - img.rows) / 2, off_x = (diag_ - img.cols) / 2;
  img.copyTo(padded(cv::Rect(off_x, off_y, img.cols, img.rows)));

  cv::Mat sino(diag_, p_.n_theta, CV_32F);
  cv::parallel_for_(cv::Range(0, p_.n_theta), [&](const cv::Range& range) {
    cv::Mat rot, col_sum;
    for (int i = range.start; i < range.end; ++i) {
      cv::warpAffine(padded, rot, rot_mats_[i], padded.size(),
                     cv::INTER_LINEAR, cv::BORDER_CONSTANT, 0);
      cv::reduce(rot, col_sum, 0, cv::REDUCE_SUM, CV_32F);  // 1 x diag
      col_sum.reshape(1, diag_).copyTo(sino.col(i));
    }
  });

  double mx;
  cv::minMaxLoc(sino, nullptr, &mx);
  if (mx > 0) sino /= mx;
  cv::resize(sino, sino, cv::Size(p_.n_theta_out, sino.rows),
             0, 0, cv::INTER_AREA);

  // |FFT| along rho: per-row FFT of the transposed sinogram
  cv::Mat sino_t, spec_c;
  cv::transpose(sino, sino_t);  // n_theta_out x rho
  cv::dft(sino_t, spec_c, cv::DFT_ROWS | cv::DFT_COMPLEX_OUTPUT);
  std::vector<cv::Mat> planes;
  cv::split(spec_c, planes);
  cv::Mat mag;
  cv::magnitude(planes[0], planes[1], mag);
  cv::Mat desc;
  cv::transpose(mag, desc);                       // rho x n_theta_out
  desc = desc.rowRange(0, sino.rows / 2).clone();  // half spectrum

  cv::Scalar mean, stddev;
  cv::meanStdDev(desc, mean, stddev);
  if (stddev[0] > 0) desc = (desc - mean[0]) / stddev[0];
  return desc;
}

double signalFraction(const cv::Mat& polar, int intensity_th) {
  return static_cast<double>(cv::countNonZero(polar > intensity_th)) /
         (static_cast<double>(polar.rows) * polar.cols);
}

Matcher::Matcher(const Params& p) : p_(p) {}

cv::Mat Matcher::thetaFFT(const cv::Mat& desc) const {
  cv::Mat f;
  cv::dft(desc, f, cv::DFT_ROWS | cv::DFT_COMPLEX_OUTPUT);
  return f;  // rho_half x n_theta_out, CV_32FC2
}

void Matcher::add(const cv::Mat& desc, double signal) {
  db_.push_back(thetaFFT(desc));
  valid_.push_back(signal >= p_.signal_gate);
}

bool Matcher::query(const cv::Mat& desc, int* best_idx, double* min_dist,
                    double signal) const {
  if (p_.signal_gate > 0 && signal < p_.signal_gate) return false;
  const int n_cand = static_cast<int>(db_.size()) - p_.exclude_recent + 1;
  if (n_cand <= 0) return false;

  cv::Mat fq = thetaFFT(desc);
  const int rows = fq.rows, T = fq.cols;

  double best = -1e300;
  int best_i = -1;
  std::mutex mu;
  cv::parallel_for_(cv::Range(0, n_cand), [&](const cv::Range& range) {
    double loc_best = -1e300;
    int loc_i = -1;
    cv::Mat cross(1, T, CV_32FC2), corr;
    for (int i = range.start; i < range.end; ++i) {
      if (!valid_[i]) continue;
      const cv::Mat& fn = db_[i];
      cross.setTo(0);
      auto* acc = cross.ptr<cv::Vec2f>(0);
      for (int r = 0; r < rows; ++r) {
        const auto* a = fq.ptr<cv::Vec2f>(r);
        const auto* b = fn.ptr<cv::Vec2f>(r);
        for (int t = 0; t < T; ++t) {  // a * conj(b)
          acc[t][0] += a[t][0] * b[t][0] + a[t][1] * b[t][1];
          acc[t][1] += a[t][1] * b[t][0] - a[t][0] * b[t][1];
        }
      }
      cv::dft(cross, corr,
              cv::DFT_INVERSE | cv::DFT_SCALE | cv::DFT_ROWS |
                  cv::DFT_REAL_OUTPUT);
      double mx;
      cv::minMaxLoc(corr, nullptr, &mx);
      if (mx > loc_best) {
        loc_best = mx;
        loc_i = i;
      }
    }
    std::lock_guard<std::mutex> lock(mu);
    if (loc_best > best) {
      best = loc_best;
      best_i = loc_i;
    }
  });
  if (best_i < 0) return false;  // every candidate was gated out
  const double autoc = cv::norm(desc, cv::NORM_L2SQR);
  *best_idx = best_i;
  *min_dist = (autoc - best) / 1000.0;
  return true;
}

}  // namespace raplace
