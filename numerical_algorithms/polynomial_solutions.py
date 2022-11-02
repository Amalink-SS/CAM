from numpy import sign
import math
from sympy import symbols


def bisection(function_to_eval, lower_int_bound: float, upper_int_bound: float, tolerance=1e-9):
    first_output = function_to_eval(lower_int_bound)
    second_output = function_to_eval(upper_int_bound)
    if first_output == 0.0:
        return lower_int_bound
    if second_output == 0.0:
        return second_output
    if sign(first_output) == sign(second_output):
        error.err('Root is not bracketed')
    number_of_iterations: int = int(math.ceil(
        math.log(abs(upper_int_bound-lower_int_bound)/tolerance)/math.log(2.0)))
    for i in range(number_of_iterations):
        mid_point = (upper_int_bound + lower_int_bound)/(2.0)
        f_at_mid_point = function_to_eval(mid_point)
        if (abs(f_at_mid_point) > abs(first_output)) and (abs(f_at_mid_point) > abs(second_output)):
            return None
        elif f_at_mid_point == 0.0:
            return mid_point
        elif sign(second_output) != sign(f_at_mid_point):
            lower_int_bound = mid_point
            first_output = f_at_mid_point
        else:
            upper_int_bound = mid_point
            second_output = f_at_mid_point
    return (lower_int_bound + upper_int_bound) / (2.0)


def function_to_eval(x):
    return x**3 - 10.0*x**2 + 5.0


print(bisection(function_to_eval,
      0.0, 1.0))


# more to be implemented
