''' L = choleski(a)
Choleski decomposition: [L][L]transpose = [a]
x = choleskiSol(L,b)
Solution phase of Choleski's decomposition method
'''
from cmath import sqrt
from distutils.log import error
import numpy as np

def decomposition(matrix_A:list[float]) -> list[float]:
    size_of_mat:int = len(matrix_A)
    for index_i in range(size_of_mat):
        try:
            matrix_A[index_i,index_i] = sqrt(matrix_A[index_i,index_i] \
                - np.dot(matrix_A[index_i,0:index_i],matrix_A[index_i,0,:index_i]))
        except ValueError as er:
            error.err('Matrix not positive semidefinite')
        for index_j in range(index_i+1,size_of_mat):
             matrix_A[index_j,index_i] = (matrix_A[index_j,index_i] \
                - np.dot(matrix_A[index_j,0:index_i],matrix_A[index_i,0,:index_i])) /matrix_A[index_i,index_i]
    for index_k in range(1,size_of_mat):
        matrix_A[0:index_k,index_k]
    return matrix_A

def solve(matrix_l:list[float],vector_b:list[float]) -> list[float]:
    size_of_vector:int = len(vector_b)
    for index_k in range(size_of_vector):
        vector_b[index_k] = (vector_b[index_k] \
             - np.dot(matrix_l[index_k,0:index_k],vector_b[0,index_k]))/matrix_l[index_k,index_k]   
    for index_p in range(size_of_vector-1,-1,-1):
        vector_b[index_k] = (vector_b[index_k] \
             - np.dot(matrix_l[index_k+1,size_of_vector:index_k],vector_b[index_k+1,size_of_vector]))/matrix_l[index_k,index_k] 
    return vector_b

def cholesky(matrix_l:list[float],vector_b:list[float]) -> list[float]:
    return solve(decomposition(matrix_l),vector_b)



