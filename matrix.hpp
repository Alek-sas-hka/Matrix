#pragma once

#include <algorithm>
#include <vector>

template <size_t N, size_t M, typename T = int64_t>
class Matrix {
 public:
  size_t GetStr() const { return str_; }

  size_t GetCol() const { return col_; }

  Matrix() { matrix_.resize(N, std::vector<T>(M)); }

  Matrix(const std::vector<std::vector<T>>& src) {
    matrix_.resize(N, std::vector<T>(M));
    matrix_ = src;
  }

  Matrix(const T& elem) { matrix_.resize(N, std::vector<T>(M, elem)); }

  Matrix(const Matrix<N, M, T>& src) : Matrix(src.matrix_) {}

  Matrix<N, M, T>& operator=(const Matrix<N, M, T>& src) {
    Matrix<N, M, T> tmp(src);
    this->matrix_.swap(tmp.matrix_);
    return *this;
  }

  friend Matrix<N, M, T> operator+(const Matrix<N, M, T>& left,
                                   const Matrix<N, M, T>& right) {
    Matrix<N, M, T> answer = left;
    answer += right;
    return answer;
  }

  friend Matrix<N, M, T> operator-(const Matrix<N, M, T>& left,
                                   const Matrix<N, M, T>& right) {
    Matrix<N, M, T> answer = left;
    answer -= right;
    return answer;
  }

  Matrix<N, M, T>& operator+=(const Matrix<N, M, T>& src) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        matrix_[i][j] += src.matrix_[i][j];
      }
    }
    return *this;
  }

  Matrix<N, M, T>& operator-=(const Matrix<N, M, T>& src) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        matrix_[i][j] -= src.matrix_[i][j];
      }
    }
    return *this;
  }

  friend Matrix<N, M, T> operator*(const Matrix<N, M, T>& multiplier,
                                   const T& elem) {
    Matrix<N, M, T> answer;
    for (size_t i = 0; i < N; i++) {
      for (size_t j = 0; j < M; j++) {
        answer.matrix_[i][j] = multiplier.matrix_[i][j] * elem;
      }
    }
    return answer;
  }

  template <size_t U, size_t K>
  friend bool CheckRightSize(const Matrix<N, M, T>& left,
                             const Matrix<U, K, T>& right) {
    return left.GetCol() == right.GetStr();
  }

  template <size_t U, size_t K>
  friend Matrix<N, K, T> operator*(const Matrix<N, M, T>& left,
                                   const Matrix<U, K, T>& right) {
    if (!CheckRightSize(left, right)) {
      exit(1);
    }
    Matrix<N, K, T> answer;
    for (size_t string = 0; string < left.GetStr(); ++string) {
      for (size_t column = 0; column < right.GetCol(); ++column) {
        answer(string, column) = 0;
        for (size_t i = 0; i < left.GetCol(); i++) {
          answer(string, column) += left(string, i) * right(i, column);
        }
      }
    }
    return answer;
  }

  Matrix<M, N, T> Transposed() {
    Matrix<M, N, T> answer;
    for (size_t string = 0; string < N; ++string) {
      for (size_t column = 0; column < M; ++column) {
        answer(column, string) = matrix_[string][column];
      }
    }
    return answer;
  }

  T operator()(const size_t& str_idx, const size_t& col_idx) const {
    return matrix_[str_idx][col_idx];
  }

  T& operator()(const size_t& str_idx, const size_t& col_idx) {
    return matrix_[str_idx][col_idx];
  }

  bool operator==(const Matrix<N, M, T>& right) {
    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++i) {
        if (*this(i, j) != right(i, j)) {
          return false;
        }
      }
    }
    return true;
  }

  T Trace() {
    T trace = 0;
    for (size_t i = 0; i < N; i++) {
      trace += matrix_[i][i];
    }
    return trace;
  }

 private:
  std::vector<std::vector<T>> matrix_;
  size_t str_ = N;
  size_t col_ = M;
};
