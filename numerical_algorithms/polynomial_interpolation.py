import numpy as np


def newton_polynomial_coefficients(x_data: list[float], y_data: list[float]) -> list[float]:
    numberOf_dataPoints: int = len(x_data)
    coeff_vector = y_data.copy()
    for k in range(1, numberOf_dataPoints):
        coeff_vector[k:numberOf_dataPoints] = (
            (np.array(coeff_vector[k:numberOf_dataPoints]) - coeff_vector[k-1])/(np.array(x_data[k:numberOf_dataPoints]) - x_data[k-1]))
    return coeff_vector


def evaluateNewton_polynomial(coefficients: list[float], x_data: list[float], x: float) -> float:
    polynomial_degree = len(x_data) - 1
    polynomial = coefficients[polynomial_degree]
    for k in range(1, polynomial_degree + 1):
        polynomial = coefficients[polynomial_degree-k] + \
            (x - x_data[polynomial_degree-k])*polynomial
    return polynomial


def show_newton_poly(x_data: list[float], y_data: list[float], x: float):
    print(evaluateNewton_polynomial(
        newton_polynomial_coefficients(x_data, y_data), x_data, 1.0))


def neville_interpolant(x_data: list[float], y_data: list[float], point: float) -> list[float]:
    numberOf_dataPoints = len(x_data)
    copy_of_y_data = y_data.copy()
    for i in range(1, numberOf_dataPoints):
        copy_of_y_data[0:numberOf_dataPoints-i] = ((
            point - np.array(x_data[i:numberOf_dataPoints]))*copy_of_y_data[0:numberOf_dataPoints-i] +
            (np.array(x_data[0:numberOf_dataPoints-i]) - point)*copy_of_y_data[1:numberOf_dataPoints-i+1]) / \
            (np.array(x_data[0:numberOf_dataPoints-i]) -
             np.array(x_data[i:numberOf_dataPoints]))
    return copy_of_y_data


def show_neville(x_data: list[float], y_data: list[float], point: float):
    print(neville_interpolant(x_data, y_data, point))


def rational_interpolation(x_data: list[float], y_data: list[float], x: float, tolerance: float) -> list[float]:
    vector_length = len(x_data)
    r_vector = y_data.copy()
    r_vector_old = np.zeros(vector_length)
    for k in range(vector_length):
        for i in range(vector_length-k-1):
            if abs(x - x_data[i+k+1]) < tolerance:
                return y_data[i+k+1]
            else:
                const_1 = r_vector[i+1] - r_vector[i]
                const_2 = r_vector[i+1] - r_vector_old[i+1]
                const_3 = (x - np.array(x_data[i])) / \
                    (x - np.array(x_data[i+k+1]))
                r_vector[i] = r_vector[i+1] + const_1 / \
                    (const_3*(1.0 - (const_1/const_3)) - 1.0)
                r_vector_old[i+1] = r_vector[i+1]
    return r_vector


def show_rational(x_data: list[float], y_data: list[float], x: float, tolerance: float):
    print(rational_interpolation(x_data, y_data, x, tolerance))


a = [0.1, 0.2, 0.5, 0.6, 0.8, 1.2, 1.5]
b = [-1.5342, -1.0811, -0.4445, -0.3085, -0.0868, 0.2281, 0.3824]

show_rational(a, b, 0.1, 1e-9)
