#pragma once

#include <exception>
#include <iomanip>

template<class T = double>
class ConfusionMatrix {
public:
  T tp = 0.0, tn = 0.0, fp = 0.0, fn = 0.0;
  /* get recall */
  T recall() const {
    if (tp + tn == 0) {
      throw std::runtime_error("tp + tn == 0");
    }
    return tp / (tp + fn);
  }
  /* get precision */
  T precision() const {
    if (tp + fp == 0) {
      throw std::runtime_error("tp + fp == 0");
    }
    return tp / (tp + fp);
  }
  /* get f1 score */
  T f1() const {
    auto p = precision();
    auto r = recall();
    if (p + r == 0) {
      throw std::runtime_error("precision + recall == 0");
    }
    return 2 * p * r / (p + r);
  }
  /* get accuracy */
  T accuracy() const {
    if (tp + tn + fp + fn == 0) {
      throw std::runtime_error("tp + tn + fp + fn == 0");
    }
    return (tp + tn) / (tp + tn + fp + fn);
  }
  /* get sensitivity */
  T sensitivity() const {
    return recall();
  }
  /* get specificity */
  T specificity() const {
    if (tn + fp == 0) {
      throw std::runtime_error("tn + fp == 0");
    }
    return tn / (tn + fp);
  }
  /* get false positive rate */
  T fpr() const {
    return 1 - specificity();
  }
  /* get false negative rate */
  T fnr() const {
    return 1 - sensitivity();
  }
  /* get false discovery rate */
  T fdr() const {
    return 1 - precision();
  }

  /* operator + */
  ConfusionMatrix operator+ (const ConfusionMatrix& mat) const {
    ConfusionMatrix ret;
    ret.tp = tp + mat.tp;
    ret.tn = tn + mat.tn;
    ret.fp = fp + mat.fp;
    ret.fn = fn + mat.fn;
    return ret;
  }

  /* operator += */
  ConfusionMatrix& operator+= (const ConfusionMatrix& mat) {
    tp += mat.tp;
    tn += mat.tn;
    fp += mat.fp;
    fn += mat.fn;
    return *this;
  }

  friend auto operator<< (std::ostream& os, const ConfusionMatrix& mat) {
    os << std::setw(6);
    os << "tp: " << mat.tp << '\t';
    os << "tn: " << mat.tn << '\n';
    os << "fp: " << mat.fp << '\t';
    os << "fn: " << mat.fn << '\n';
    return os;
  }

};

/* create operator << for output confusion matrix */

