#include <iostream>
#include <cmath>
#include "Matrix.h"

using namespace std;

int main() {
    int m, n, rk;
    double eps;
    std::cout<<"Enter m,n for matrix A: ";
    std::cin>>m>>n;
    Matrix A(m,n);
    std::cout<<"Enter elements of matrix A: "<<std::endl;
    std::cin>>A;
    std::cout<<"Matrix:"<<std::endl;
    std::cout<<A;

    std::cout<<"Enter m,n for matrix B: ";
    std::cin>>m>>n;
    Matrix B(m,n);
    std::cout<<"Enter elements of matrix B: "<<std::endl;
    std::cin>>B;

    Matrix C = A + B;
    std::cout<<"A + B:"<<std::endl;
    std::cout<<C;
    Matrix D = A - B;
    std::cout<<"A - B:"<<std::endl;
    std::cout<<D;
    Matrix F = A * B;
    std::cout<<"A * B:"<<std::endl;
    std::cout<<F;

    std::cout << "Введите точность eps" << std::endl;
    std::cin>>eps;

    Matrix G = B;
    rk = G.gauss_method(eps, G);
    std::cout<<"Ступенчатый вид матрицы B: " << std::endl;
    std::cout<<G<<std::endl;
    std::cout<<"ранг матрицы B = ";
    std::cout<<rk<<std::endl;

    if (m == n) {
        double det = 1;
        for (int i = 0; i < n; ++i) {
            det *= G.at(i,i);
        }
        std::cout<<"Определитель матрицы B = ";
        std::cout<<det<<std::endl;

        if (det !=0){
            F = B.algebraic_complement_transp(B);
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < n; ++j) {
                    F.at(i,j) = F.at(i,j)/det;
                }
            }
            std::cout<<"Обратная матрица к B:"<<std::endl;
            std::cout<<F<<std::endl;
        }
    }

    Matrix H = G.lin_syst_solve(G);
    std::cout<<"Решение системы B = ";
    std::cout<<H<<std::endl;

    return 0;
}
