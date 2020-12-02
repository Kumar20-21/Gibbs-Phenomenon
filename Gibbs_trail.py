from matplotlib import pyplot as plt
import numpy as np
import math as ma


"""
from sumpy.interactive import printing
printing.init_printing(use_latex = true)
"""


def function_triangle_wave(a, b, x):
    L = b-a
    x = x - (x//L+1)*L-a
    if x <= 0:
        return x+1
    else:
        return 1-x


def function_square_wave(a, b, x):
    L = b - a
    x = x - (x//L+1)*L-a
    if x < 0:
        return 1
    else:
        return 0


def sawtooth_function(a, b, x):
    L = b - a
    x = x - (x//L+1)*L-a
    return x


def x_square(a, b, x):
    L = b - a
    x = x - (x//L+1)*L
    return x*x


def midpoint_rule_integration(a, b, n, comp_int_point):
    L = b-a
    h = L/comp_int_point
    sum1 = 0
    for jj in range(comp_int_point):
        x_star = ((a+jj*h)+(a+(jj+1)*h))/2.0
        sum1 += sawtooth_function(a, b, x_star)*np.exp(-2*np.pi*1j*n*x_star/L)
    sum1 = h*sum1/L
    return sum1


def complex_to_real(comp_arr, real_arr, m):
    if abs(comp_arr[m].real) < 10**(-12):  # Setting degree of precision to be 10^(-11)
        real_arr += [0.0]
    else:  # Truncating the number after 11 digits
        real_arr += [ma.trunc(comp_arr[m].real * (10 ** 11)) / (10 ** 11)]

    for alp in range(m):
        if abs(comp_arr[m+1+alp].real) < 10**(-12):
            real_arr += [0.0]
        else:
            real_arr += [2*ma.trunc(comp_arr[m+1+alp].real*(10**11))/(10**11)]
        if abs(comp_arr[m+1+alp].imag) < 10**(-12):
            real_arr += [0.0]
        else:
            real_arr += [-2*ma.trunc(comp_arr[m+1+alp].imag*(10**11))/(10**11)]
    real_arr = np.array(real_arr)
    return real_arr.transpose()


def plotting_graph(a, b, real_arr, n, x):
    L = b-a
    sum1 = real_arr[0]
    for jj in range(n):
        sum1 += real_arr[2*jj+1]*np.cos(2*np.pi*(jj+1)*x/L)
        sum1 += real_arr[2*(jj+1)]*np.sin(2*np.pi*(jj+1)*x/L)
    return sum1


def error_at_point(function_value_at_x, sum_value):
    er = abs(function_value_at_x-sum_value)*100
    return er


l_i_p = -1  # For storing left end point of the interval
r_i_p = 1  # For storing right end point of the interval

l_i_pp = -2  # The left end point for graph plotting
r_i_pp = 2  # The right end point for graph plotting

no_plotting = 2000  # Number of points for graph plotting
no_integration = 1000  # Number of points for composite numerical integration
no_terms = 100  # Number of terms in the Fourier series expansion

x_array = []  # Array for storing x values for graph plotting
y_sawtooth = []  # Array for storing y values for sawtooth function plotting
y_x_square = []  # Array for storing y values for x**2 function plotting
y_triangle = []  # Array for storing y values for triangular wave function plotting
y_square_wave = []  # Array for storing y values for square wave function 

comp_Four_arr = []  # Array for storing complex Fourier coefficients
real_Four_arr = []  # Array for storing real Fourier coefficients

sum_value_arr = []  # Array for storing function value obtained by Fourier sum at a point
err_arr = []  # Array for storing error between actual function value and Fourier sum at a point

for ii in range(-no_terms, no_terms+1):
    comp_Four_arr += [midpoint_rule_integration(l_i_p, r_i_p, ii, no_integration)]

real_Four_arr = complex_to_real(comp_Four_arr, real_Four_arr, no_terms)  # Array for storing real Fourier coefficients

for ii in range(no_plotting+1):
    h_val = (r_i_pp - l_i_pp) / no_plotting
    x_prime = l_i_pp + h_val * ii
    x_array += [x_prime]
    y_sawtooth += [sawtooth_function(l_i_p, r_i_p, x_prime)]
    sum_value_arr += [plotting_graph(l_i_p, r_i_p, real_Four_arr, no_terms, x_prime)]
    err_arr += [error_at_point(y_sawtooth[ii], sum_value_arr[ii])]


plt.plot(x_array, y_sawtooth, 'red', linestyle='solid', linewidth=1, label='Sawtooth function')
plt.plot(x_array, sum_value_arr, 'black', linestyle='dashed', linewidth=1, label='Fourier Reconstruction')
plt.title('Gibbs Phenomenon for sawtooth')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.legend()
plt.show()
