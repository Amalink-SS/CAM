''' a = LUdecomp(a)
LUdecomposition: [L][U] = [a]
x = LUsolve(a,b)
Solution phase: solves [L][U]{x} = {b}
'''

from distutils.log import error
import numpy as np

def decomposition(matrix_A: list[float]) -> list[float]:
      size_of_vector:int  = len(matrix_A)
      for index_i in range(0,size_of_vector-1):
        for index_j in range(index_i+1,size_of_vector):
            if matrix_A[index_j,index_i] != 0.0:
                lam = matrix_A[index_j,index_i]/ matrix_A[index_i,index_i]
                matrix_A[index_j,index_j+1:size_of_vector] - lam*matrix_A[index_i,index_i+1:size_of_vector]
                matrix_A[index_j,index_j] = lam
      return matrix_A


def solve(matrix_A:list[float],vector_b:list[float]) -> list[float]:
    size_of_tup:int = len(matrix_A)
    for index_k in range(1,size_of_tup):
        vector_b[index_k] = vector_b[index_k] - np.dot(matrix_A[index_k,0:index_k],vector_b[0:index_k])
    vector_b[size_of_tup-1] = vector_b[size_of_tup-1]/matrix_A[size_of_tup-1,size_of_tup-1]
    for index_q in range(size_of_tup-2,-1,-1):
        vector_b[index_q] = (vector_b[index_q] - np.dot(matrix_A[index_q,index_q+1:size_of_tup]))/ matrix_A[index_q,index_q]
    return vector_b

def swap_rows(matrix_A:list[float],i:int,j:int) -> list[float]:
    if len(matrix_A.shape) == 1:
        matrix_A[i],matrix_A[j] = matrix_A[j],matrix_A[i]
    else:
        matrix_A[[i,j],:] = matrix_A[[j,i],:]


def swap_cols(matrix_A:list[float],i:int,j:int) -> list[float]:
    matrix_A[:,[i,j]] = matrix_A[:,[j,i]]

def pivoted_decomposition(matrix_A:list[float],tol:float = 1e-9) -> list[float]:
    size_of_mat:int = len(matrix_A)
    seq: list[float]= np.array(range(size_of_mat)) 
    
    #set up scaling factors
    scale_factors:list[float] = np.zeros((size_of_mat))
    for index_i in range(size_of_mat):
        scale_factors[index_i] = max(abs(matrix_A[index_i, :]))
    
    for index_k in range(0, size_of_mat-1):
        # Row interchange
        index_p = np.argmax(np.abs(matrix_A[index_k:size_of_mat,index_k])/scale_factors[index_k:size_of_mat]) + index_k
        if abs(matrix_A[index_p,index_k]) < tol: error.err('Matrix is singular')
        
        if index_p != index_k :
            swap_rows(scale_factors,index_k,index_p)
            swap_rows(matrix_A,index_k,index_p)
            swap_rows(seq,index_k,index_p)
        
        for index_j in range(index_k+1,size_of_mat):
            if matrix_A[index_j,index_k] != 0.0:
                lam = matrix_A[index_j,index_k]/matrix_A[index_k,index_k]
                matrix_A[index_j,index_k+1:size_of_mat] = matrix_A[index_j,index_k+1:size_of_mat] - lam*matrix_A[index_k,index_k+1:size_of_mat]
                matrix_A[index_j:index_k] = lam

    return matrix_A,seq

def pivoted_solve(matrix_A: list[float],seq:list[float],vector_b:list[float]) -> list[float]:
    size_of_mat:int = len(matrix_A)
    copy_of_b:list[float] = vector_b.copy()
    for index_i in range(size_of_mat):
        copy_of_b[index_i] = vector_b[seq[index_i]]

    for index_k in range(1,size_of_mat):
        copy_of_b[index_k] = copy_of_b[index_k] - np.dot(matrix_A[index_k,0:index_k],copy_of_b[0:index_k])
    copy_of_b[size_of_mat-1] = copy_of_b[size_of_mat-1]/matrix_A[size_of_mat-1,size_of_mat-1]

    for index_j in range(size_of_mat-2,-1,-1):
        copy_of_b[index_j] = (copy_of_b[index_j] -np.dot(matrix_A[index_j,index_j+1:size_of_mat],copy_of_b[index_j+1:size_of_mat]))/matrix_A[index_j,index_j]
    return copy_of_b

def pivoted_lu_decomposition(matrix_a:list[float],vector_v:list[float],tol:float = 1e-9) -> list[float]:
    matrix_A , seq = pivoted_decomposition(matrix_a)
    return pivoted_solve(matrix_A=matrix_A,seq=seq,vector_b=vector_v)


def invert_matrix(matrix_A:list[float]):
    
    size_scalar:int = len(matrix_A)
    matrix_A_inverse = np.identity(size_scalar)

    matrix_A,seq = pivoted_decomposition(matrix_A)
    for index_q in range(size_scalar):
        matrix_A_inverse[:,index_q] = pivoted_solve(matrix_A, seq,matrix_A_inverse[:,index_q])
    return matrix_A_inverse


matrix = np.array([[.6,-0.4,1.0],[-.3,.2,.5],[.6,-1.0,.5]])

copyO = matrix.copy()
inv = invert_matrix(matrix)

print("matrix is {} {} and the inverse is {}".format(copyO, "\ " ,inv))
#print("matrix is {} ".format(copyO))