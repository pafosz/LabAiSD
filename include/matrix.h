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
	T Det() const;
	void Print() const;
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
bool operator==(const Matrix<T>& A, const Matrix<T>& B) {
	if (A.get_rows() != B.get_rows() || A.get_cols() != B.get_cols()) {
		return false;
	}
	const double kEpsilon = 1.0E-5;
	for (int i = 0; i < A.get_rows(); ++i) {
		for (int j = 0; j < A.get_cols(); ++j) {
			if (A(i, j) - B(i, j) > kEpsilon) return false;
		}
	}
	return true;
}

template <typename T>
bool operator!=(const Matrix<T>& A, const Matrix<T>& B) {
	return (!(A == B));
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

template <typename T>
T Matrix<T>::Det() const {
	if(rows_ != cols_)
		throw  std::invalid_argument("The determinant exists only for square matrices");
	int N = rows_;
	if(N == 2) return (data_[0][0] * data_[1][1]) - (data_[0][1] * data_[1][0]);
	T det = 0;
	
	// Рекурсивно вычисляем определитель
	for (int col = 0; col < N; col++) {
		double** subMatrix = new double* [N - 1];
		for (int i = 0; i < N - 1; i++) {
			subMatrix[i] = new double[N - 1];
		}

		for (int i = 1; i < N; i++) {
			int subMatrixCol = 0;
			for (int j = 0; j < N; j++) {
				if (j != col) {
					subMatrix[i - 1][subMatrixCol] = data_[i][j];
					subMatrixCol++;
				}
			}
		}

		double sign = (col % 2 == 0) ? 1 : -1;
		det += sign * data_[0][col] * subMatrix.calculateDeterminant(); //subMatrix, N - 1

		// Don't forget to free the memory
		for (int i = 0; i < N - 1; i++) {
			delete[] subMatrix[i];
		}
		delete[] subMatrix;
	}

	return det;




}

template <typename T>
void Matrix<T>::Print() const {
	for (int i = 0; i < rows_; ++i) {
		for (int j = 0; j < cols_; ++j) {
			std::cout << std::setw(5) << data_[i][j] << " ";
		}
		std::cout << "\n";
	}
}

//template <typename T>
//double* solveLinearSystem(const Matrix<T>& A, const double* b) {
//	int len = sizeof(b)/sizeof(int);
//	static double x[len];
//	double y[len];
//
//	for (int i = 0; i < len; ++i) {
//		double sum = 0.0;
//		for (int j = 0; j < len; ++j) {
//			sum += A(i, j) * x[j];
//		}
//		y[i] = b[i] - sum;
//	}
//
//	for (int i = len - 1; i >= 0; --i) {
//		double sum = 0.0;
//		for (int j = i + 1; j < len; ++j) {
//			sum += A(i, j) * x[j];
//		}
//		x[i] = (1.0 / A(i, i)) * (y[i] - sum);
//	}
//
//	return x;
//}

//template <typename T>
//std::vector<double> solveLinearSystem(const Matrix<T>& A, const std::vector<double>& b) {
//	int n = A.get_rows();
//	std::vector<double> x(n);
//	std::vector<double> y(n);
//
//	for (int i = 0; i < n; ++i) {
//		double sum = 0.0;
//		for (int j = 0; j < n; ++j) {
//			sum += A(i, j) * x[j];
//		}
//		y[i] = b[i] - sum;
//	}
//
//	for (int i = n - 1; i >= 0; --i) {
//		double sum = 0.0;
//		for (int j = i + 1; j < n; ++j) {
//			sum += A(i, j) * x[j];
//		}
//		x[i] = (1.0 / A(i, i)) * (y[i] - sum);
//	}
//
//	return x;
//}
template <typename T>
double* solveLinearSystem(const Matrix<T>& A, const double* vx, const int size_vx) {
	if (A.get_cols() != size_vx)
		throw std::invalid_argument("The dimensions of the matrix and vector do not match");

	double* vb = new double[size_vx];
	std::fill(vb, vb + size_vx, 0);  // initialize the array with zeros

	for (int i = 0; i < A.get_rows(); ++i) {
		for (int j = 0; j < A.get_cols(); ++j) {
			vb[i] += A(i, j) * vx[j];
		}
	}
	return vb;
}