#include <complex>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <type_traits>
#include <utility>


template <typename T>
class Matrix {
public:
	Matrix();
	Matrix(int rows, int cols, T init_value);
	Matrix(int rows, int cols, T lower_bound, T upper_bound);
	Matrix(const Matrix<T>& other);

	void Swap(Matrix<T>& other) noexcept;
	Matrix<T>& operator=(const Matrix<T> other);

	~Matrix();

	T& operator()(int row, int col);
	const T& operator()(int row, int col) const;

	Matrix<T>& operator+=(const Matrix<T>& other);
	Matrix<T> operator+(const Matrix<T>& other) const;
	Matrix<T> operator-(const Matrix<T>& other) const;
	Matrix<T> operator*(const Matrix<T>& other) const;
	Matrix<T> operator*(const T scalar) const;
	Matrix<double> operator/(const T scalar) const;

	int get_rows() const;
	int get_cols() const;

	T Trace() const;
	void Print() const;

private:
	T** data_;
	int rows_;
	int cols_;
};

template <typename T>
Matrix<T>::Matrix() : data_(nullptr), rows_(0), cols_(0) {}

template <typename T>
Matrix<T>::Matrix(int rows, int cols, T init_value) : rows_(rows), cols_(cols) {
	data_ = new T * [rows];
	for (int i = 0; i < rows; i++) {
		data_[i] = new T[cols];
		for (int j = 0; j < cols; j++) {
			data_[i][j] = init_value;
		}
	}
}

template <typename T>
Matrix<T>::Matrix(int rows, int cols, T lower_bound, T upper_bound) : rows_(rows), cols_(cols) {
	std::random_device rd;
	std::mt19937 gen(rd());

	if constexpr (std::is_integral_v<T>) {
		std::uniform_int_distribution<T> dist(lower_bound, upper_bound);

		data_ = new T * [rows];
		for (int i = 0; i < rows; ++i) {
			data_[i] = new T[cols];
			for (int j = 0; j < cols; ++j) {
				data_[i][j] = dist(gen);
			}
		}
	}
	else if constexpr (std::is_floating_point_v<T>) {
		std::uniform_real_distribution<T> dist(lower_bound, upper_bound);

		data_ = new T * [rows];
		for (int i = 0; i < rows; ++i) {
			data_[i] = new T[cols];
			for (int j = 0; j < cols; ++j) {
				data_[i][j] = dist(gen);
			}
		}
	}
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T>& other) : rows_(other.rows_), cols_(other.cols_) {
	data_ = new T * [rows_];
	for (int i = 0; i < rows_; ++i) {
		data_[i] = new T[cols_];
		for (int j = 0; j < cols_; ++j) {
			data_[i][j] = other.data_[i][j];
		}
	}
}

template <typename T>
Matrix<T>::~Matrix() {
	for (int i = 0; i < rows_; ++i) {
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
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& other) {
	if (rows_ != other.rows_ || cols_ != other.cols_) {
		throw std::invalid_argument("Matrix dimensions do not match");
	}
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < cols_; ++j) {
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
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < cols_; ++j) {
			res(i, j) = data_[i][j] - other.data_[i][j];
		}
	}
	return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& other) const {
	if (cols_ != other.rows_) {
		throw std::invalid_argument("Matrix dimensions do not match");
	}
	Matrix res(rows_, other.cols_, 0);
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < other.cols_; ++j) {
			for (int k = 0; k < cols_; ++k) {
				res(i, j) += data_[i][k] * other.data_[k][j];
			}
		}
	}
	return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const T scalar) const {
	Matrix<T> res(rows_, cols_, 0);
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < cols_; ++j) {
			res(i, j) = data_[i][j] * scalar;
		}
	}
	return res;
}

template <typename T>
Matrix<double> Matrix<T>::operator/(const T scalar) const {
	if (!scalar) {
		throw std::invalid_argument("Divide by zero!");
	}
	Matrix<double> res(rows_, cols_, 0);
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < cols_; ++j) {
			res(i, j) = static_cast<double>(data_[i][j]) / scalar;
		}
	}
	return res;
}

template <typename T>
int Matrix<T>::get_rows() const {
	return rows_;
}

template <typename T>
int Matrix<T>::get_cols() const {
	return cols_;
}

template <typename T>
bool operator==(const Matrix<T>& matrix_first, const Matrix<T>& matrix_second) {
	if (matrix_first.get_rows() != matrix_second.get_rows() || matrix_first.get_cols() != matrix_second.get_cols()) {
		return false;
	}
	const double kEpsilon = 1.0E-5;
	for (int i = 0; i < matrix_first.get_rows(); ++i) {
		for (int j = 0; j < matrix_first.get_cols(); ++j) {
			if (matrix_first(i, j) - matrix_second(i, j) > kEpsilon) 
				return false;
		}
	}
	return true;
}

template <>
bool operator==(const Matrix<std::complex<double>>& matrix_first, const Matrix<std::complex<double>>& matrix_second) {
	if (matrix_first.get_rows() != matrix_second.get_rows() || matrix_first.get_cols() != matrix_second.get_cols()) {
		return false;
	}
	const double kEpsilon = 1.0E-5;
	for (int i = 0; i < matrix_first.get_rows(); ++i) {
		for (int j = 0; j < matrix_first.get_cols(); ++j) {
			if ((abs(matrix_first(i, j)) - abs(matrix_second(i, j))) > kEpsilon)
				return false;
		}
	}
	return true;
}

template <typename T>
bool operator!=(const Matrix<T>& matrix_first, const Matrix<T>& matrix_second) {
	return (!(matrix_first == matrix_second));
}

template <>
bool operator!=(const Matrix<std::complex<double>>& matrix_first, const Matrix<std::complex<double>>& matrix_second) {
	return (!(matrix_first == matrix_second));
}

template <typename T>
T Matrix<T>::Trace() const {
	if (rows_ != cols_)
		throw std::invalid_argument("Matrix is not square");
	T res = 0;
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < cols_; ++j) {
			if (i == j) res += data_[i][j];
		}
	}
	return res;
}

template <>
std::complex<double> Matrix<std::complex<double>>::Trace() const {
	if (rows_ != cols_)
		throw std::invalid_argument("Matrix is not square");
	std::complex<double> res = 0;
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < cols_; ++j)
			if (i == j) res += data_[i][j];
	}		
	return res;
}
		
template <typename T>
void Matrix<T>::Print() const {
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < cols_; ++j) {
			std::cout.precision(4);
			std::cout << std::setw(5) << data_[i][j] << " ";
		}
		std::cout << "\n";
	}
}
