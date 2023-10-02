//
// Created by Anhelina Modenko on 24.09.2023.
//

#ifndef LINALG_UNI_MATRIX_H
#define LINALG_UNI_MATRIX_H


class Matrix {

public:
    /*
     * CONSTRUCTORS
     */

    Matrix();
    Matrix(unsigned int rows, unsigned int cols);

    /*
     * OPERATORS
     */

    Matrix operator+(const Matrix& other) const;
    Matrix operator-(const Matrix& other) const;
    Matrix operator*(const Matrix& other) const;


    /*
     * FUNCTIONALITY
     */

    /*
     * allocate memory for new matrix
     */
    void allocateMemory(unsigned int row, unsigned int col);

    /*
     * Identity matrix
     */
    Matrix identity(unsigned int size);

    /*
     * Zero matrix
     */
    Matrix zero(unsigned int rows, unsigned int cols);

    /*
     * New matrix from array
     */
    Matrix fromArray(double* arr, unsigned int rows, unsigned int cols);

    /*
     * Print matrix
     */
    void print() const;

    /*
     * Get matrix from file
     */
    Matrix fromFile(const char* filename);

    /*
     * Set matrix values
     */
    void setValues(double** values);

    /*
     * Multiplication of matrices
     * function
     */
    Matrix multiply(const Matrix& other) const;

    /*
     * Strassen multiplication
     * function
     */
    Matrix strassenMultiply(const Matrix& other) const;

    /*
     * Addition of matrices
     */
    void add(const Matrix& other) const;

    /*
     * Addition of matrices (result is new matrix)
     */
    Matrix add(const Matrix& other, Matrix& result) const;

    /*
     * Subtraction of matrices
     */
    Matrix subtract(const Matrix& other) const;

    /*
     * Matrix from col
     */
    Matrix fromCol(double* arr, unsigned int rows);

    /*
     * Matrix from row
     */
    Matrix fromRow(double* arr, unsigned int cols);

    /*
     * Are matrices equal
     */
    bool equals(const Matrix& other) const;

    /*
     * Multiple matrix by scalar
     */
    Matrix multiplyByScalar(double scalar) const;

    /*
     * Transpose matrix
     */
    Matrix transpose() const;

    /*
     * Swap rows
     */
    void swapRows(unsigned int row1, unsigned int row2) const;

    /*
     * Swap columns
     */
    void swapCols(unsigned int col1, unsigned int col2) const;

    /*
     * performs Gaussian elimination to transform the matrix
     * to row echelon form REF
     */
    Matrix gaussianEliminationREF() const;

    /*
     * performs Gaussian elimination to transform the matrix
     * to reduced row echelon form RREF
     */
    Matrix gaussianEliminationRREF() const;

    /*
     * Augment matrix on the end of a matrix
     * @param matrix the original matrix
     * @param augment the matrix to augment
     */
    Matrix augment(const Matrix& augment) const;

    /*
     * Get submatrix
     */
    Matrix subMatrix(unsigned int row, unsigned int col, unsigned int rows, unsigned int cols) const;

    /*
     * Get determinant
     */
    double determinant() const;

    /*
     * Get inverse matrix
     */
    Matrix inverse() const;

    /*
     * Get cofactor
     */
    Matrix cofactor() const;

    /*
     * Get adjugate
     */
    Matrix adjugate() const;

    /*
     * Bool is invertible
     */
    bool isInvertible() const;

    /*
     * calculate the inverse of a matrix using Gaussian Elimination
     */
    Matrix inverseGaussianElimination();

    /*
     * clearing up memory
     */
    void clear();

    /*
     * OPERATORS
     */
    Matrix& operator=(const Matrix& other);

    /*
     * DESTRUCTOR
     */

    ~Matrix();

private:

    unsigned int rows;
    unsigned int cols;

    double** data;
};


#endif //LINALG_UNI_MATRIX_H
