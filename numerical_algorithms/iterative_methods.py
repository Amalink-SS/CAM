from mpmath import sqrt
import numpy as np
import math


def gauss_seidel(iterator, matrix_x: list[float], relaxation_coef: float, number_iter: int, extra_iter: int, tol: float = 1e-9):

    for index_i in range(1, 1000):
        Copy_of_matrix_x = matrix_x.copy()
        matrix_x = iterator(matrix_x, relaxation_coef)
        dx = math.sqrt(np.dot(matrix_x-Copy_of_matrix_x,
                       matrix_x-Copy_of_matrix_x))

        if dx < tol:
            return matrix_x, index_i, relaxation_coef
        if index_i == number_iter:
            dy = dx
        if index_i == number_iter + extra_iter:
            dy = dx
            relaxation_coef = 2.0 / \
                (1.0 + sqrt(1 - pow((dy/dx), (1/extra_iter))))
    else:
        print("Method failed to converge")


def iterator(matrix_x: list[float], relaxation_coef: float):
    size_scalar: int = len(matrix_x)
    matrix_x[0] = relaxation_coef*(matrix_x[1] - matrix_x[size_scalar-1]) / \
        2.0 + (1.0 - relaxation_coef)*matrix_x[0]
    for index_i in range(1, size_scalar-1):
        matrix_x[index_i] = relaxation_coef * \
            (matrix_x[index_i-1] + matrix_x[index_i+1]) / \
            2.0 - (1.0 - relaxation_coef)*matrix_x[index_i]
    matrix_x[size_scalar-1] = relaxation_coef*(1.0 - matrix_x[0] + matrix_x[size_scalar - 2])/2.0 \
        + (1.0 - relaxation_coef)*matrix_x[size_scalar-1]
    return matrix_x


matrix_x = np.array([[6, -2, 1], [-2, 7, 2], [1, 2, -5]])

vector_y = np.array([11, 5, -1])

matrix_x, number_iter, relaxation_coef = gauss_seidel(
    iterator, matrix_x, 1.0, 10, 1)

print(
    "solution is {}".format(matrix_x),
    "Number of Iterations is {}".format(number_iter),
    "Relaxation coefficient is {}".format(relaxation_coef)

)
