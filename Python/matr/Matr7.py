from fractions import *
class Matrix:
    '''Matrix n x m'''
    def __init__(self, n = 0, m = 0):
        self.n = n
        self.m = m
        self.matrix = [[Fraction(0)]*m for i in range(n)]

    def set_matrix(self, list):
        for i in range(self.n):
            for j in range(self.m):
                self.matrix[i][j] = Fraction(list[i*self.m+j])
        return self
    
    def __str__(self):
        s = "\n" + "\n".join([str(i) for i in self.matrix]) + "\n"
        return s
    
    def __add__(self, B):
        assert type(B) == Matrix
        C=Matrix(self.n,self.m)
        if type(B) == Matrix:
            for i in range(self.n):
                for j in range(self.m):
                    C.matrix[i][j]=Fraction(self.matrix[i][j]+B.matrix[i][j])
            return C
        
    def __iadd__(self, B):
        assert type(B) == Matrix
        if type(B) == Matrix:
            for i in range(self.n):
                for j in range(self.m):
                    self.matrix[i][j]=Fraction(self.matrix[i][j]+B.matrix[i][j])
        return self
    
    def __sub__(self, B):
        assert type(B) == Matrix
        C=Matrix(self.n,self.m)
        if type(B) == Matrix:
            for i in range(self.n):
                for j in range(self.m):
                    C.matrix[i][j]=Fraction(self.matrix[i][j]-B.matrix[i][j])
            return C
        
    def __isub__(self, B):
        assert type(B) == Matrix
        if type(B) == Matrix:
            for i in range(self.n):
                for j in range(self.m):
                    self.matrix[i][j]=Fraction(self.matrix[i][j]-B.matrix[i][j])
        return self
    
    def __mul__(self, k):
        if type(k)==Matrix:
            return(self.multiplyMatrices(k))
        '''Multiply a Matrix by a number'''
        B=Matrix(self.n,self.m)
        for i in range(self.n):
            for j in range(self.m):
                B.matrix[i][j] = Fraction(self.matrix[i][j]*k)
        return B
    
    def __matmul__(self, b):
        return(self.multiplyMatrices(b))
    
    def __imul__(self, k):
        '''Multiply a Matrix by a number'''
        B=Matrix(self.n,self.m)
        for i in range(self.n):
            for j in range(self.m):
                self.matrix[i][j] = Fraction(self.matrix[i][j]*k)
        return self
    
    def __getitem__(self, i):
        assert i[0] <= self.n and i[1] <= self.m
        if (i[0] <= self.n and i[1] <= self.m):
            return self.matrix[i[0]][i[1]]
        else:
            raise IndexError("index out of range")
        
    def __setitem__(self, i, a):
        assert i[0] <= self.n and i[1] <= self.m
        if (i[0] <= self.n and i[1] <= self.m):
            self.matrix[i[0]][i[1]]=Fraction(a)
        else:
            raise IndexError("index out of range") 
         
    def copyRationalMatrix(self):
        m = self.m
        assert(m > 0)
        n = self.n
        assert(n > 0)
        b = Matrix(n,m)
        for i in range(n):
            for j in range(m):
                b.matrix[i][j]=self.matrix[i][j]
        return b
    
    def showRationalMatrix(self):
        res = ""
        m = self.m
        n = self.n
        maxLen = 0
        for i in range(m):
            for j in range(n):
                txt = str(self.matrix[i][j])
                l = len(txt)
                if l > maxLen:
                    maxLen = l
    
        for i in range(m):
            if i > 0:
                res += '\n'     # Line delimiter
            res += '['
            for j in range(n):
                if j > 0:
                    res += ' '  # Element delimiter
                txt = str(self.matrix[i][j])
                l = len(txt)
                for k in range(maxLen - l):
                    res += ' '
                res += txt
            res += ']'
        return res
    
    def swapMatrixRows(self, i, k):
        assert(0 <= i and i < self.n and 0 <= k and k < self.n)
        if i == k:
            return
        m = self.m
        assert(m == len(self.matrix[k]))
        for j in range(m):
            (self.matrix[i][j], self.matrix[k][j]) = (self.matrix[k][j], -self.matrix[i][j])
        
    def addMatrixRows(self, k, i, r):
        assert(0 <= i and i < self.n and 0 <= k and k < self.n)
        m = self.m
        assert(len(self.matrix[k]) == m)
        for j in range(m):
            self.matrix[k][j] += self.matrix[i][j]*r
        return self
  
    def echelonFormOfRationalMatrix(self):
        b = self.copyRationalMatrix()
        m = self.m      # number of rows
        n = self.n   # number of columns
        i = 0           # number of rows processed    
        j = 0           # number of columns processed
        while i < n and j < m:
            # 1. Looking for nonzero element in j-th column
            #    starting from i-th row
            iMax = i
            for k in range(i, n):
                if b.matrix[k][j] != 0:
                    iMax = k
                    break
            if b.matrix[iMax][j] == 0:
                # Zero column
                j += 1
                continue
            
            assert(b.matrix[iMax][j] != 0)
            if iMax != i:
                b.swapMatrixRows(i, iMax)
            assert(b.matrix[i][j] != 0)
                
            # 2. Annihilate elements of j-th column,
            #    starting from i+1-th row
            c = b.matrix[i][j]
            for k in range(i+1, n):
                b.addMatrixRows(k, i, -b.matrix[k][j]/c)
            i += 1
            j += 1    
        return (b,i)          

    def determinantOfRationalMatrix(self):
        if self.n == 0 or self.n != self.m:
            raise ValueError("Determinant of non-square matrix")
        (b, rank) = self.echelonFormOfRationalMatrix()
        d = b.matrix[0][0]
        for i in range(1, b.m):
            d *= b.matrix[i][i]
        return d

    def multiplyMatrices(self, b):
        if self.n == 0 or self.m == 0:
            raise ValueError("Empty matrix a") 
        if self.n == 0 or self.m == 0:
            raise ValueError("Empty matrix b")
        rows_a = self.n
        cols_a = self.m
        rows_b = b.n
        cols_b = b.m
        if cols_a != rows_b:
            raise ValueError("Matrices cannot be multiplied")
    
        # Create a zero matrix of size rows_a * cols_b
        C=Matrix(rows_a,cols_b)
        for i in range(rows_a):
            for j in range(cols_b):
                for k in range(cols_a):
                    C.matrix[i][j] += self.matrix[i][k]*b.matrix[k][j]
        return C
    
    def det_minor(self, i, j):
        B=Matrix(self.n-1, self.m-1)
        for k in range(0,i):
            for l in range(0,j):
                B.matrix[k][l]=self.matrix[k][l]
            for l in range(j+1, self.m):
                B.matrix[k][l-1] = self.matrix[k][l]
        for k in range(i+1,self.n):
            for l in range(0,j):
                B.matrix[k-1][l] = self.matrix[k][l]
            for l in range(j+1, self.m):
                B.matrix[k-1][l-1] = self.matrix[k][l]
        b = B.echelonFormOfRationalMatrix()[0]
        d = b.matrix[0][0]
        for i in range(1, b.m):
            d = d*b.matrix[i][i]
        return d


    def algebraic_complement_transp(self):
        B=Matrix(self.n,self.m)
        for i in range(0,self.n):
            for j in range(0,self.m):
                B.matrix[i][j]=((-1)**(i+j))*(self.det_minor(i, j))
        C=Matrix(B.m,B.n)
        for i in range(B.n):
            for j in range(B.m):
                C.matrix[i][j]=B.matrix[j][i]
        return C
    
    def transpose(self):
        B=Matrix(self.m,self.n)
        for i in range(B.n):
            for j in range(B.m):
                B.matrix[i][j]=self.matrix[j][i]
        return B
    
    def lin_syst_solve(self):
        if self.n + 1 != self.m:
            raise ValueError("число строк и столбцов определены не корректно")
        B=Matrix(1,self.m-1)
        s=Fraction(0)
        for i in range(0,self.m-1):
            B.matrix[0][self.m-2-i]=Fraction(self.matrix[self.n-i-1][self.m-1] - s)/Fraction(self.matrix[self.n-i-1][self.m-2-i])
            s=Fraction(0)
            for j in range(0,i+1):
                s+=self.matrix[self.n-i-2][self.m-2-j]*B.matrix[0][self.m-2-j]
        return B


    def inv_matrix(self):
        if (self.determinantOfRationalMatrix()==0):
            raise ValueError("Вырожденная матрица")
        if (self.m!=self.n):
            raise ValueError("Не квадратная матрица")
        F = self.algebraic_complement_transp()
        for i in range(self.n):
            for j in range(self.m):
                F.matrix[i][j] = F.matrix[i][j]/self.determinantOfRationalMatrix()
        return F
       

        


    '''
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
'''



    


    
    