from math import exp

import Methods
import matplotlib.pyplot as plt
import numpy as np
from scipy import *

# input your slope function as a string
# i.e. for the DE y'(x)=f(x,y). In this linear DE example
# f(x,y)-1-3x+y but you can enter whatever function you want.
"""
#x_0 = 1
#y_0 = 0
x_0 = 0
y_0 = 0

#b = 2  # max x-value
b = 1
n = 10  # number of steps

#actual = np.log(2)
actual = np.pi
print(actual)
done = False
count = 0

def math_function(x, y):
    #return 1 / x
    return 4 / (x**2 + 1)


while not done:
    estimate = Methods.runge_kutta(x_0, y_0, b, n, math_function)[1][n]
    print(n, estimate)

    error = np.abs(estimate - actual) / np.sqrt(estimate ** 2 + actual ** 2)
    if error < 0.00000001:
        done = True

    n = n * 2



def diff(x, y):
    return 9.8 - 0.2 * y

x = np.linspace(0, 50, 40)
y = np.linspace(0, 60, 40)

# use x,y
for j in x:
    for k in y:
        slope = diff(j, k)
        domain = np.linspace(j-0.1, j+0.1, 2)
        plt.plot(domain, slope*(domain - j) + k, solid_capstyle='projecting', solid_joinstyle='bevel')

plt.title("Slope field y'")
plt.grid(True)
plt.show()

"""
# Draw the slope field here
fig = plt.figure(num=1, figsize=(8, 8))
ax = fig.add_subplot(111)


# Vector field function
def vf(t, v):
    return v


def integrand(t):
    return np.e**t

'''
# Solution curves
t0 = 0
tEnd = 10
dt = 0.01

ic = [[0, 0]]
color = ['b']

for k in range(len(ic)):
    Y = []
    T = []

    step = 0.1
    for i in range(0, 51):
        T.append(i * step)
        Y.append(integrand(i * step))

    ax.plot(T, Y, color=color[k], lw=1.25)

# Vector field
T, V = np.meshgrid(np.linspace(0, 1, 20), np.linspace(0, 4, 30))
W = 12
X = V
# Normalize arrows
N = np.sqrt(W**2+X**2)
U2, V2 = W/N, X/N
ax.quiver(T, V, U2, V2)
'''


def dy_dx(x, y):
    return y


a = 0
b = 1
y_0 = 1

euler_errors = []
trapeze_errors = []
simpson_errors = []
runge_errors = []

steps = range(3, 16)

''' Euler's method '''
for step in steps:
    n = 2**step

    ''' Get the y-terms '''
    euler = Methods.eulers(a, b, y_0, n, dy_dx)[1]
    runge = Methods.runge_kutta(a, b, y_0, n, dy_dx)[1]
    trapazoid = Methods.trapazoid(a, b, y_0, n)[1]
    simpson = Methods.simpson(a, b, y_0, exp(a + 2**-step), n)[1]

    euler_errors.append(np.abs(np.e - euler[n]))
    runge_errors.append(np.abs(np.e - runge[n]))
    trapeze_errors.append(np.abs(np.e - trapazoid[n]))
    simpson_errors.append(np.abs(np.e - simpson[n]))


print(euler_errors)
print(runge_errors)
print(trapeze_errors)
print(simpson_errors)
# plt.yscale('log')
ax.scatter(steps, np.array(euler_errors), color='red', label='Euler')
ax.scatter(steps, np.array(runge_errors), color='blue', label='Runge-Kutta')
ax.scatter(steps, np.array(trapeze_errors), color='orange', label='Trapezoidal')
ax.scatter(steps, np.array(simpson_errors), color='green', label='Simpson’s-Rule')
ax.legend()

plt.title('Error for Varying Methods')
plt.xlabel('2^n steps')
plt.ylabel('Error: y(1)-y_n')
plt.savefig('euler-runge.png')

'''
    ax.plot(euler[0], euler[1], color='red', linestyle='dashdot', label='Euler')
    ax.legend()

    plt.title('Euler\'s Method: y\'=y')
    plt.xlim([0, 1])
    plt.ylim([1, 3])
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    plt.savefig('euler.png')
'''
