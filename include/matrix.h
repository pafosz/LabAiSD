#pragma once

#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <utility>

template <typename T>
class Matrix {
private:
	T** data_;
	int rows_;
	int cols_;

public:
	Matrix(int rows, int cols, T init_value);
	Matrix(int rows, int cols, T lower_bound, T upper_bound);
	Matrix(const Matrix<T>& other);

	void Swap(Matrix<T>& other) noexcept;
	Matrix<T>& operator=(const Matrix<T> other);

	~Matrix();

	T& operator()(int row, int col);
	const T& operator()(int row, int col) const;

	Matrix<T> operator+=(const Matrix<T>& other);
	Matrix<T> operator+(const Matrix<T>& other) const;
	Matrix<T> operator-(const Matrix<T>& other) const;
	Matrix<T> operator*(const Matrix<T>& other);
	Matrix<T> operator*(const T scalar) const;
	Matrix<T> operator/(const T scalar);

	T Trace() const;
	void Print() const;
};


template <typename T>
Matrix<T>::Matrix(int rows, int cols, T init_value) : rows_(rows), cols_(cols) {
	data_ = new T*[rows];
	for (size_t i = 0; i < rows; i++) {
		data_[i] = new T[cols];
		for (size_t j = 0; j < cols; j++) {
			data_[i][j] = init_value;
		}
	}
}

template <typename T>
Matrix<T>::Matrix(int rows, int cols, T lower_bound, T upper_bound) : rows_(rows), cols_(cols) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<T> dist(lower_bound, upper_bound);

	data_ = new T * [rows];
	for (int i = 0; i < rows; ++i) {
		data_[i] = new T[cols];
		for (int j = 0; j < cols; ++j) {
			data_[i][j] = dist(gen);
		}
	}
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& other) : rows_(other.rows_), cols_(other.cols_) {
	data_ = new T * [rows_];
	for (size_t i = 0; i < rows_; ++i) {
		data_ = new T [cols_];
		for (size_t j = 0; j < cols_; ++j) {
			data_[i][j] = other.data_[i][j];
		}
	}
}

template <typename T>
Matrix<T>::~Matrix() {
	for (size_t i = 0; i < rows_; ++i) {
		delete[] data_[i];
	}
	delete[] data_;
	data_ = nullptr;
}

template <typename T>
void Matrix<T>::Swap(Matrix<T>& other) noexcept {
	std::swap(data_, other.data_);
	std::swap(rows_, other.rows_);
	std::swap(cols_, other.cols_);
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T> other) {
	Swap(other);
	return *this;
}

template <typename T>
T& Matrix<T>::operator()(int row, int col) {
	return data_[row][col];
}

template <typename T>
const T& Matrix<T>::operator()(int row, int col) const {
	return data_[row][col];
}

template <typename T>
Matrix<T> Matrix<T>::operator+=(const Matrix<T>& other) {
	if (rows_ != other.rows_ || cols_ != other.cols_) {
		throw std::invalid_argument("Matrix dimensions do not match");
	}
	for (size_t i = 0; i < rows_; ++i) {
		for (size_t j = 0; j < cols_; ++j) {
			data_[i][j] += other.data_[i][j];
		}
	}
	return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) const {
	if (rows_ != other.rows_ || cols_ != other.cols_) {
		throw std::invalid_argument("Matrix dimensions do not match");
	}
	Matrix<T> res = *this;
	res += other;	
	return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& other) const {
	if (rows_ != other.rows_ || cols_ != other.cols_) {
		throw std::invalid_argument("Matrix dimensions do not match");
	}
	Matrix res(rows_, cols_, 0);
	for (size_t i = 0; i < rows_; ++i) {
		for (size_t j = 0; j < cols_; ++j) {
			res(i, j) = data_[i][j] - other.data_[i][j];
		}
	}
	return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& other) {
	if (cols_ != other.rows_) {
		throw std::invalid_argument("Matrix dimensions do not match");
	}
	Matrix res(rows_, other.cols_, 0);
	for (size_t i = 0; i < rows_; ++i) {
		for (size_t j = 0; j < other.cols_; ++j) {
			for (size_t k = 0; k < cols_; ++k) {
				res(i, j) += data_[i][k] * other.data_[k][j];
			}
		}
	}
	return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const T scalar) const {
	Matrix<T> res(rows_, cols_, 0);
	for (size_t i = 0; i < rows_; ++i) {
		for (size_t j = 0; j < cols_; ++j) {
			res(i, j) = data_[i][j] * scalar;
		}
	}
	return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator/(const T scalar) {
	Matrix<T> res(rows_, cols_, 0);
	for (size_t i = 0; i < rows_; ++i) {
		for (size_t j = 0; j < cols_; ++j) {
			res(i, j) = data_[i][j] / scalar;
		}
	}
	return res;
}

template <typename T>
T Matrix<T>::Trace() const {
	if (rows_ != cols_)
		throw std::invalid_argument("Matrix is not square");
	T res = 0;
	for (size_t i = 0; i < rows_; ++i) {
		for (size_t j = 0; j < cols_; ++j) {
			if (i == j) res += data_[i][j];
		}
	}
	return res;
}

template <typename T>
void Matrix<T>::Print() const {
	for (size_t i = 0; i < rows_; ++i) {
		for (size_t j = 0; j < cols_; ++j) {
			std::cout << std::setw(5) << data_[i][j] << " ";
		}
		std::cout << "\n";
	}
}