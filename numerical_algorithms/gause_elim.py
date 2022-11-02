## module liear_systems
''' x = gaussElimin(a,b).
Solves [a]{b} = {x} by Gauss elimination.
'''
import numpy as np
def gauss_elimination(matrix_A: list[float],vector_b:list[float]) -> list[float]:
    size_of_vector: int = len(vector_b)
    # Elimmination Phase
    for index_k in range(0,size_of_vector-1):
        for index_i in range(index_k+1,size_of_vector):
            if matrix_A[index_i,index_k] != 0.0:
                lam : float = matrix_A[index_i,index_k] / matrix_A[index_k,index_k]
                matrix_A[index_i,index_k+1:size_of_vector] = matrix_A[index_i,index_k+1:size_of_vector] -lam*matrix_A[index_i,index_k+1:size_of_vector]
                vector_b[index_i] = vector_b[index_i] - lam*vector_b[index_i]
    # Backward Substitution
    for index_kk in range(size_of_vector-1,-1,-1):
        vector_b[index_kk] = (vector_b[index_kk] - np.dot(matrix_A[index_kk,index_kk+1:size_of_vector], vector_b[index_kk+1:size_of_vector]))/(matrix_A[index_kk,index_kk])
    return vector_b

def vandemonde(vector_bb: list[float]) -> list[float]:
    size_of_vector_bb : int = len(vector_bb)
    matrix_aa = np.zeros((size_of_vector_bb,size_of_vector_bb))
    for index_j in range(size_of_vector_bb):
        matrix_aa[:,index_j] = vector_bb**size_of_vector_bb-index_j-1
    return matrix_aa

v = np.array([2.0,1.8,1.6,1.4,1.2,1.0])
b = np.array([0.0, 1.0, 0.0])

#a = vandemonde(v)

a = np.array([[4.00,-2.00,1.00],
      [-2.00,4.00,-2.00],
      [1.00,-2.00,4.00]])
aOrig = a.copy()
bOrg = b.copy()

x = gauss_elimination(a,b)

print("x is {}".format(x))
print()
print("the residual is {}".format(np.dot(aOrig,x) - bOrg))
print()
print("determinant is {}".format(np.prod(np.diagonal(a))))



    
