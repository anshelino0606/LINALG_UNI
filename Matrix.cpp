//
// Created by Anhelina Modenko on 24.09.2023.
//

#include <stdexcept>
#include <iostream>
#include <fstream>
#include "Matrix.h"

Matrix::Matrix()
    : rows(0), cols(0), data(nullptr) {}

Matrix::Matrix(unsigned int rows, unsigned int cols)
    : rows(rows), cols(cols) {
    data = new double*[rows];
    for (int i = 0; i < rows; ++i) {
        data[i] = new double[cols];
    }
}

Matrix Matrix::operator+(const Matrix& other) const {
    if (rows != other.rows || cols != other.cols) {
        throw std::invalid_argument("ERROR::DIFFERENT SIZED MATRICES");
    }

    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.data[i][j] = data[i][j] + other.data[i][j];
        }
    }

    return result;
}

Matrix Matrix::operator-(const Matrix& other) const {
    if (rows != other.rows || cols != other.cols) {
        throw std::invalid_argument("ERROR::DIFFERENT SIZED MATRICES");
    }

    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.data[i][j] = data[i][j] - other.data[i][j];
        }
    }

    return result;
}

Matrix Matrix::operator*(const Matrix& other) const {
    if (cols != other.rows) {
        throw std::invalid_argument("ERROR::DIFFERENT SIZED MATRICES");
    }

    Matrix result(rows, other.cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < other.cols; ++j) {
            for (int k = 0; k < cols; ++k) {
                result.data[i][j] += data[i][k] * other.data[k][j];
            }
        }
    }

    return result;
}

Matrix Matrix::identity(unsigned int size) {
    Matrix result(size, size);
    for (int i = 0; i < size; ++i) {
        result.data[i][i] = 1;
    }

    return result;
}

Matrix Matrix::zero(unsigned int rows, unsigned int cols) {
    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        result.data[i][i] = 0;
    }

    return result;
}

Matrix Matrix::fromArray(double* arr, unsigned int rows, unsigned int cols) {
    Matrix result(rows, cols);
    int k = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j, ++k) {
            result.data[i][j] = arr[k];
        }
    }

    return result;
}

void Matrix::print() const {
    std::cout << "Matrix: " << std::endl;
    for (int i = 0; i < rows; ++i) {
        std::cout << "[ ";
        for (int j = 0; j < cols; ++j) {
            std::cout << data[i][j] << " ";
        }
        std::cout << "]" << std::endl;
    }
}

Matrix Matrix::fromFile(const char *filename) {
    // get matrix from a file

    return Matrix();
}

Matrix Matrix::multiply(const Matrix& other) const {
    return *this * other;
}

void Matrix::add(const Matrix& other) const {
    *this + other;
}

Matrix Matrix::subtract(const Matrix &other) const {
    return *this - other;
}

Matrix Matrix::fromCol(double *arr, unsigned int rows) {
    return fromArray(arr, rows, 1);
}

Matrix Matrix::fromRow(double *arr, unsigned int cols) {
    return fromArray(arr, 1, cols);
}

bool Matrix::equals(const Matrix &other) const {
    if (rows != other.rows || cols != other.cols) {
        return false;
    }

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < other.cols; ++j) {
            if (data[i][j] != other.data[i][j]) {
                return false;
            }
        }
    }

    return true;
}

Matrix Matrix::multiplyByScalar(double scalar) const {
    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            result.data[i][j] = data[i][j] * scalar;
        }
    }

    return result;
}

Matrix Matrix::transpose() const {
    Matrix result(cols, rows);
    for (int i = 0; i < cols; ++i) {
        for (int j = 0; j < rows; ++j) {
            result.data[i][j] = data[j][i];
        }
    }

    return result;
}

Matrix Matrix::gaussianEliminationREF() const {
    unsigned int currentRow = 0;
    for (int col = 0; col < cols; col++) {
        unsigned int row = currentRow;

        if (row >= this->rows) break;

        // non-zero entry
        for (;row < this->rows; row++) {
            if (data[row][col] != 0.0f) break;
        }

        if (row == this->rows) continue;

        swapRows(row, currentRow);

        float factor = 1 / this->data[currentRow][col];
        for (int i = col; i < this->cols; i++) {
            this->data[currentRow][i] *= factor;
        }

        for (int i = currentRow + 1; i < this->rows; i++) {
            float factor = -this->data[i][col];
            for (int j = col; j < this->cols; j++) {
                this->data[i][j] += factor * this->data[currentRow][j];
            }
        }

        currentRow++;
    }
    return *this;
}

void Matrix::swapRows(unsigned int row1, unsigned int row2) const {
    double* temp = data[row1];
    data[row1] = data[row2];
    data[row2] = temp;
}

void Matrix::swapCols(unsigned int col1, unsigned int col2) const {
    for (int i = 0; i < rows; ++i) {
        double temp = data[i][col1];
        data[i][col1] = data[i][col2];
        data[i][col2] = temp;
    }
}

Matrix Matrix::gaussianEliminationRREF() const {

    unsigned int currentRow = 0;
    for (int col = 0; col < cols; col++) {
        unsigned int row = currentRow;

        if (row >= this->rows) break;

        // non-zero entry
        for (;row < this->rows; row++) {
            if (data[row][col] != 0.0f) break;
        }

        if (row == this->rows) continue;

        swapRows(row, currentRow);

        float factor = 1 / this->data[currentRow][col];
        for (int i = col; i < this->cols; i++) {
            this->data[currentRow][i] *= factor;
        }

        for (int i = currentRow + 1; i < this->rows; i++) {
            float factor = -this->data[i][col];
            for (int j = col; j < this->cols; j++) {
                this->data[i][j] += factor * this->data[currentRow][j];
            }
        }

        currentRow++;
    }
    return *this;
}

Matrix Matrix::augment(const Matrix &augment) const {
    if (rows != augment.rows) {
        throw std::invalid_argument("ERROR::DIFFERENT SIZED MATRICES");
    }

    Matrix result(rows, cols + augment.cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols + augment.cols; ++j) {
            if (j < cols) {
                result.data[i][j] = data[i][j];
            } else {
                result.data[i][j] = augment.data[i][j - cols];
            }
        }
    }

    return result;
}

Matrix Matrix::subMatrix(unsigned int row, unsigned int col, unsigned int rows, unsigned int cols) const {
    Matrix result(rows, cols);
    for (int i = row; i < row + rows; ++i) {
        for (int j = col; j < col + cols; ++j) {
            result.data[i - row][j - col] = data[i][j];
        }
    }

    return result;
}

double Matrix::determinant() const {
    if (rows != cols) {
        throw std::invalid_argument("ERROR::DIFFERENT SIZED MATRICES");
    }

    if (rows == 1) {
        return data[0][0];
    }

    if (rows == 2) {
        return data[0][0] * data[1][1] - data[0][1] * data[1][0];
    }

    double result = 0;
    for (int i = 0; i < cols; ++i) {
        result += data[0][i] * cofactor().data[0][i];
    }

    return result;
}

Matrix Matrix::inverse() const {
    if (rows != cols) {
        throw std::invalid_argument("ERROR::DIFFERENT SIZED MATRICES");
    }

    if (!isInvertible()) {
        throw std::invalid_argument("ERROR::NOT INVERTIBLE");
    }

    return adjugate().multiplyByScalar(1 / determinant());
}

Matrix Matrix::cofactor() const {
    if (rows != cols) {
        throw std::invalid_argument("ERROR::DIFFERENT SIZED MATRICES");
    }

    Matrix result(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; j++) {
            result.data[i][j] = subMatrix(i, j, rows - 1, cols - 1).determinant();
        }
    }

    return result;
}

Matrix Matrix::adjugate() const {
    if (rows != cols) {
        throw std::invalid_argument("ERROR::DIFFERENT SIZED MATRICES");
    }

    return cofactor().transpose();
}

bool Matrix::isInvertible() const {
    return determinant() != 0;
}

Matrix Matrix::inverseGaussianElimination() {
    if (rows != cols) {
        throw std::invalid_argument("ERROR::DIFFERENT SIZED MATRICES");
    }

    if (!isInvertible()) {
        throw std::invalid_argument("ERROR::NOT INVERTIBLE");
    }

    Matrix result = augment(identity(rows)).gaussianEliminationRREF();
    return result.subMatrix(0, rows, rows, rows);
}

Matrix::~Matrix() {
    for (int i = 0; i < rows; ++i) {
        delete[] data[i];
    }
    delete[] data;
}

Matrix Matrix::strassenMultiply(const Matrix &other) const {
if (cols != other.rows) {
        throw std::invalid_argument("ERROR::DIFFERENT SIZED MATRICES");
    }

    if (rows == 1 && cols == 1 && other.rows == 1 && other.cols == 1) {
        Matrix result(1, 1);
        result.data[0][0] = data[0][0] * other.data[0][0];
        return result;
    }

    unsigned int newSize = rows / 2;
    Matrix a11 = subMatrix(0, 0, newSize, newSize);
    Matrix a12 = subMatrix(0, newSize, newSize, newSize);
    Matrix a21 = subMatrix(newSize, 0, newSize, newSize);
    Matrix a22 = subMatrix(newSize, newSize, newSize, newSize);

    Matrix b11 = other.subMatrix(0, 0, newSize, newSize);
    Matrix b12 = other.subMatrix(0, newSize, newSize, newSize);
    Matrix b21 = other.subMatrix(newSize, 0, newSize, newSize);
    Matrix b22 = other.subMatrix(newSize, newSize, newSize, newSize);

    Matrix p1 = a11.strassenMultiply(b12 - b22);
    Matrix p2 = (a11 + a12).strassenMultiply(b22);
    Matrix p3 = (a21 + a22).strassenMultiply(b11);
    Matrix p4 = a22.strassenMultiply(b21 - b11);
    Matrix p5 = (a11 + a22).strassenMultiply(b11 + b22);
    Matrix p6 = (a12 - a22).strassenMultiply(b21 + b22);
    Matrix p7 = (a11 - a21).strassenMultiply(b11 + b12);

    Matrix c11 = p5 + p4 - p2 + p6;
    Matrix c12 = p1 + p2;
    Matrix c21 = p3 + p4;
    Matrix c22 = p5 + p1 - p3 - p7;

    Matrix result(rows, cols);
    // insert c11
    for (int i = 0; i < newSize; ++i) {
        for (int j = 0; j < newSize; ++j) {
            result.data[i][j] = c11.data[i][j];
        }
    }

    // insert c12
    for (int i = 0; i < newSize; ++i) {
        for (int j = newSize; j < cols; ++j) {
            result.data[i][j] = c12.data[i][j - newSize];
        }
    }

    // insert c21
    for (int i = newSize; i < rows; ++i) {
        for (int j = 0; j < newSize; ++j) {
            result.data[i][j] = c21.data[i - newSize][j];
        }
    }

    // insert c22
    for (int i = newSize; i < rows; ++i) {
        for (int j = newSize; j < cols; ++j) {
            result.data[i][j] = c22.data[i - newSize][j - newSize];
        }
    }

    return result;
}

Matrix Matrix::add(const Matrix &other, Matrix &result) const {

    // check if matrices are of the same size
    if (rows != other.rows || cols != other.cols) {
        throw std::invalid_argument("ERROR::DIFFERENT SIZED MATRICES");
    }

    // add two matrices (to the current)
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < other.cols; ++j) {
            result.data[i][j] = data[i][j] + other.data[i][j];
        }
    }

    return result;
}

void Matrix::setValues(double **values) {

    // set values of the matrix
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            data[i][j] = values[i][j];
        }
    }

}

Matrix &Matrix::operator=(const Matrix &other) {
    if (this == &other) {
        return *this;
    }

    // delete old data
    this->clear();

    // =
    rows = other.rows;
    cols = other.cols;

    // allocate new data
    this->allocateMemory(rows, cols);

    // copy data
    setValues(other.data);

    return *this;
}

void Matrix::allocateMemory(unsigned int row, unsigned int col) {

    data = new double*[row];
    for (int i = 0; i < row; ++i) {
        data[i] = new double[col];
    }

    rows = row;
    cols = col;
}

void Matrix::clear() {

    // clear matrix
    for (int i = 0; i < rows; ++i) {
        delete[] data[i];
    }
    delete[] data;

    rows = 0;
    cols = 0;
    data = nullptr;

}

bool Matrix::operator==(const Matrix& other) const {
    if (rows != other.rows || cols != other.cols) {
        return false;
    }

    // check if matrices are equal
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < other.cols; ++j) {
            if (data[i][j] != other.data[i][j]) {
                return false;
            }
        }
    }

    return true;
}

bool operator!=(const Matrix& other1, const Matrix& other2) {
    return !(other1 == other2);
}

std::ostream& operator << (std::ostream& os, const Matrix& matrix) {

    matrix.print();
    return os;

}


