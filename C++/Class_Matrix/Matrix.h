#include <cmath>
#include <iostream>

class Matrix {
private:
    int m;
    int n;
    double* elems;
public:
    //Конструктор по умолчанию:
    Matrix(){
        m=1;
        n=1;
        elems=new double [m*n];
        elems[0]=0;
    }
    double rows() const { return m; }     
    double cols() const { return n; }
    double* data() const { return elems;}

    //Конструктор
    Matrix(int rows, int cols):
        m(rows),
        n(cols),
        elems(new double[m*n])
    {
        for(int i=0; i<m*n;++i)
        {
            elems[i]=0;
        }
    }
    // Copy constructor
    Matrix(const Matrix& a){
        m=a.m;
        n=a.n;
        elems=new double[m*n];
        for(int i=0; i<m*n; ++i){
            elems[i]=a.elems[i];
        }
    }
    ~Matrix()
    {
        delete[] elems;
    }
    double at(int i, int j) const { // x = a.at(i, j);
        return elems[i*n + j];
    }
    double& at(int i, int j) {      // a.at(i, j) = 2.73;
        return elems[i*n + j];
    }
    const double* operator[](int i) const { // x = a[i][j];
        return elems + i*n;
    }
    double* operator[](int i) { // a[i][j] = 1.5;
        return elems + i*n;
    }
    Matrix operator+(const Matrix& matrix){
        Matrix result(m,n);
        if (m == matrix.m && n == matrix.n){
            for (int i = 0; i < m * n; ++i)
            {
                result.elems[i] = elems[i] + matrix.elems[i];
            }
            return result;
        }
         else {
            std::cout<<"Operation impossible ";
            return result;
        }
    }
    Matrix operator-(const Matrix& matrix){
        Matrix result(m,n);
        if (m == matrix.m && n == matrix.n){
            for (int i = 0; i < m * n; ++i)
            {
                result.elems[i] = elems[i] - matrix.elems[i];
            }
            return result;
        }
         else {
            std::cout<<"Operation impossible ";
            return result;
        }
    }
    Matrix operator*(const Matrix& matrix){
        Matrix result(m,matrix.n);
        if (n == matrix.m)
        {
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < matrix.n; j++)
                {
                    result.elems[i * result.n +j] = 0;
                    for (int k = 0; k < n; k++)
                    {
                        result.elems[i * result.n + j] += elems[i * n + k] * matrix.elems[k * matrix.n + j];
                    }
                }
            }
            return result;      
        }
        else {
            std::cout<<"Operation impossible ";
            return result;
        }
    }
    int gauss_method(double eps, Matrix& matrix);
    Matrix transpose(Matrix& matrix);
    Matrix algebraic_complement_transp(Matrix& matrix);
    double det_minor(Matrix& matrix, int i, int j);
    Matrix lin_syst_solve(Matrix& matrix);
};
/*
    Complex& operator+=(const Complex& z) {
        a += z.a;
        b += z.b;
        return *this;
    }
    
    Complex& operator-=(const Complex& z) {
        a -= z.a;
        b -= z.b;
        return *this;
    }
    
    Complex& operator*=(const Complex& z) {
        double r = a*z.a - b*z.b;
        double i = a*z.b + b*z.a;
        a = r;
        b = i;
        return *this;
    }
*/

inline std::istream& operator>>(std::istream& s, Matrix& matrix) {
    for (int i =0; i < matrix.rows(); ++i)
    {
        for (int j =0; j < matrix.cols(); ++j)
        {
            s >> matrix.at(i, j);
        }
    }
    return s;
}

inline std::ostream& operator<<(std::ostream& s, const Matrix& matrix) {
    for (int i =0; i<matrix.rows(); ++i)
    {
        for(int j=0; j<matrix.cols(); ++j)
        {
            if (j > 0)
                s << " ";
            s << matrix.at(i, j);
        }
        s << std::endl;
    }
    return s;
}

// Bad style!
// using namespace std;



/*
class Complex {
private:
    double a;
    double b;
    // z = r*e^(i*phi)
public:
    Complex(double x = 0., double y = 0.):
        a(x),
        b(y)
    {}
    
    // Bad style:
    // Complex(double x = 0., double y = 0.) {
    //     a = x; b = y;
    // }
    
    // Destructor
    // ~Complex() {}
    
    // Copy constructor
    Complex(const Complex& z):
        a(z.a),
        b(z.b)
    {}
    
    double real() const { return a; }   // double r = z.real();
    double& real() { return a; }        // z.real() = 0.5;
    double imag() const { return b; }
    double& imag() { return b; }
    
    double abs() const {
        return sqrt(a*a + b*b);
    }
    
    double abs2() const {
        return a*a + b*b;
    }
    
    double arg() const {
        return atan2(b, a);
    }
    ...
    // z = u + v;  z = u.operator+(v);
    Complex operator+(const Complex& z) const {
        return Complex(a + z.a, b + z.b);
    }
    ...

    Complex& operator+=(const Complex& z) {
        a += z.a;
        b += z.b;
        return *this;
    }
    
    Complex& operator-=(const Complex& z) {
        a -= z.a;
        b -= z.b;
        return *this;
    }
    
    Complex& operator*=(const Complex& z) {
        double r = a*z.a - b*z.b;
        double i = a*z.b + b*z.a;
        a = r;
        b = i;
        return *this;
    }
    
    Complex conjugate() const {
        return Complex(a, -b);
    }
    
    Complex inverse() const {
        double mod2 = abs2();
        double r = a/mod2;
        double i = (-b)/mod2;
        return Complex(r, i);
    }
    
    Complex& operator/=(const Complex& z) {
        *this *= z.inverse();
        return *this;
    }
    
    void roots(int n, Complex* root) const;
};

inline Complex operator+(
    const Complex& u, const Complex& v
) {
    Complex res = u;
    res += v;
    return res;
}

inline Complex operator-(
    const Complex& u, const Complex& v
) {
    Complex res = u;
    res -= v;
    return res;
}

inline Complex operator*(
    const Complex& u, const Complex& v
) {
    Complex res = u;
    res *= v;
    return res;
}

inline Complex operator/(
    const Complex& u, const Complex& v
) {
    Complex res = u;
    res /= v;
    return res;
}

inline std::ostream& operator<<(std::ostream& s, const Complex& z) {
    s << z.real() << " + " << z.imag() << "*i";
    return s;
}

inline std::istream& operator>>(std::istream& s, Complex& z) {
    s >> z.real() >> z.imag();
    return s;
}

#endif
*/
