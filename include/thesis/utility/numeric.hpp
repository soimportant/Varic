#pragma once
#include <concepts>
#include <numeric>

/**
 * @brief C++26 feature
 * @link https://en.cppreference.com/w/cpp/numeric/add_sat
 * @link https://en.cppreference.com/w/cpp/numeric/sub_sat
 * @tparam T
 */

template<class T>
concept addable = requires(T a, T b) {
  { a + b } -> std::convertible_to<T>;
};

template<class T>
concept subtractable = requires(T a, T b) {
  { a - b } -> std::convertible_to<T>;
};

template<addable T>
constexpr T add_sat(T a, T b) {
  if (a > std::numeric_limits<T>::max() - b) {
    return std::numeric_limits<T>::max();
  }
  return a + b;
}

template<subtractable T>
constexpr T sub_sat(T a, T b) {
  if (a < std::numeric_limits<T>::min() + b) {
    return std::numeric_limits<T>::min();
  }
  return a - b;
}

