#pragma once

#include <iostream>
#include <cmath>
#include<complex>
#include<random>
#include <utility>
#include <stdexcept>

template<typename T>
class Matrix{
	
	T** data_;
	int rows_, cols_;

public:
	Matrix(int rows, int cols, T init_value);

	Matrix(int rows, int cols, T lower_bound, T upper_bound);
	
	Matrix(const Matrix<T>& other);

	void swap(Matrix<T>& other) noexcept;

	Matrix<T>& operator=(const Matrix<T>& other);

	~Matrix();

	T& operator()(int row, int col);

	Matrix<T> operator+(const Matrix<T>& other);

	Matrix<T> operator-(const Matrix<T>& other);

	Matrix<T> operator*(const Matrix<T>& other);

	Matrix<T> operator*(const T scalar) const;

	friend Matrix<T> operator*(const T scalar, const Matrix<T>& matrix);

	Matrix<T> operator/(const T scalar);

	T Trace() const;

};

template<typename T>
Matrix<T>::Matrix(int rows, int cols, T init_value) : rows_(rows), cols_(cols), data_(){

}