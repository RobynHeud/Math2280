import numpy as np
from scipy import optimize


def runge_kutta(a, b, y_0, n, f):
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = np.zeros(n + 1)

    y[0] = y_0

    for i in range(1, n + 1):
        k1 = f(x[i - 1], y[i - 1])
        k2 = f(x[i - 1] + h / 2, y[i - 1] + k1 * h / 2)
        k3 = f(x[i - 1] + h / 2, y[i - 1] + k2 * h / 2)
        k4 = f(x[i], y[i - 1] + h * k3)
        y[i] = y[i - 1] + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6

    return x, y


def eulers(a, b, y_0, n, f):
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = np.zeros(n + 1)

    y[0] = y_0

    for i in range(1, n + 1):
        y[i] = y[i - 1] + h * f(x[i - 1], y[i - 1])

    return x, y

def simpson(a, b, y_0, y_1, n):

    h = (b - a) / n

    x = np.linspace(a, b, n + 1)
    y = np.zeros(n + 1)

    y[0] = y_0
    y[1] = y_1

    for i in range(2, n + 1):
        y[i] = -(h * y[i - 2] + 4 * h * y[i - 1] + 3 * y[i - 2]) / (h - 3)

    return x, y


def trapazoid(a, b, y_0, n):
    h = (b - a) / n

    x = np.linspace(a, b, n + 1)
    y = np.zeros(n + 1)

    y[0] = y_0

    for i in range(1, n + 1):
        y[i] = (2 * y[i - 1] + h * y[i - 1]) / (2 - h)

    return x, y

def improved_eulers(a, b, y_0, n, f):
    h = (b - a) / n
    x = np.linspace(a, b, n + 1)
    y = np.zeros(n + 1, dy)

    y[0] = y_0

    for i in range(1, n + 1):
        k1 = f(x[i - 1], y[i - 1])
        u = y[i - 1] + h * k1
        k2 = f(x[i], u)
        y[i] = y[i - 1] + h * (k1 + k2) / 2

    return x, y
