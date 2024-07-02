#include <cmath>
#include "Matrix.h"

Matrix Matrix::transpose(Matrix& matrix) {
    Matrix matrix2(matrix.cols(),matrix.rows());
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            matrix2.at(i,j) =  matrix.at(j,i);
        }
    }  
    return matrix2;
}

Matrix Matrix::algebraic_complement_transp(Matrix& matrix) {
    Matrix matrix2(matrix.cols(),matrix.rows());
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            matrix2.at(i,j) =  pow((-1),(i+j))*det_minor(matrix,i,j);
        }
    }    
    Matrix matrix3(matrix2.cols(),matrix2.rows());
    for (int i = 0; i < matrix2.rows(); ++i) {
        for (int j = 0; j < matrix2.cols(); ++j) {
            matrix3.at(i,j) =  matrix2.at(j,i);
        }
    }
    return matrix3;
}

double Matrix::det_minor(Matrix& matrix, int i, int j) {
    double epsl = 0.0000001;
    Matrix matrix2(matrix.rows()-1, matrix.cols()-1);
    int rk;
    for (int k = 0; k < i; ++k) {
        for (int l = 0; l < j; ++l) {
            matrix2.at(k,l) = matrix.at(k,l);
        }
        for (int l = j+1; l < matrix.cols(); ++l) {
            matrix2.at(k,l-1) = matrix.at(k,l);
        }
    }
    for (int k = i+1; k < matrix.rows(); ++k) {
        for (int l = 0; l < j; ++l) {
            matrix2.at(k-1,l) = matrix.at(k,l);
        }
        for (int l = j+1; l < matrix.cols(); ++l) {
            matrix2.at(k-1,l-1) = matrix.at(k,l);
        }
    }
    rk = matrix2.gauss_method(epsl, matrix2);
    double det = 1;
    for (int i = 0; i < matrix2.cols(); ++i) {
        det *= matrix2.at(i,i);
    }
    return det;
}

int Matrix::gauss_method(double eps, Matrix& matrix) {
    double a;
    int i,j,k,l;
    i = 0;
    j = 0;
    int m = matrix.rows();
    int n = matrix.cols();

    while (i < m && j < n) {
        a=0;
        for (k = i; k < m; ++k){
            if (fabs(matrix.at(k,j)) > a) {
                l = k;
                a = fabs(matrix.at(k,j));
            }
        }
        if (a <= eps) {
            for (k = i; k < m; ++k) {
                matrix.at(k,j) = 0;
            }
            ++j;
            continue;
        }

        if (l != i) {
            for (k = j; k < n; ++k) {
                a = matrix.at(i,k);
                matrix.at(i,k) = matrix.at(l,k);
                matrix.at(l,k) = (-a);
            }

        }

        a = matrix.at(i,j);
        assert(fabs(a) > eps);
        for (k = i + 1; k < m; ++k) {
            double c = (-matrix.at(k,j)) / a;
            matrix.at(k,j) = 0;
            for (l = j + 1; l < n; ++l) {
                matrix.at(k,l) += c * matrix.at(i,l);
            }
        }
        ++i; 
        ++j;
    }
    return i;
}

Matrix Matrix::lin_syst_solve(Matrix& matrix) {
    Matrix X_vector(1, matrix.cols() - 1);
    double s = 0;
    for (int i = 0; i < matrix.cols() - 1; ++i) {
        X_vector.at(0,matrix.cols() - 2 - i) = (matrix.at(matrix.rows()-i-1,matrix.cols()-1) - s)/matrix.at(matrix.rows()-i -1,matrix.cols()-2-i);
        s = 0;

        //s = matrix.at(matrix.rows()-i-2,matrix.cols()-2-i)*X_vector.at(0,matrix.cols() - 2 - i);

        for (int j = 0; j <= i; ++j) {
            s += matrix.at(matrix.rows()-i-2,matrix.cols()-2-j)*X_vector.at(0,matrix.cols() - 2 - j);
        }
    }
    return X_vector;
}

