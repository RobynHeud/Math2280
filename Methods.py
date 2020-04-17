import numpy as np
from scipy import optimize


def runge_kutta(a, b, y_0, n, f):
    h = (b - a) / n
    x = np.linspace(a, b, n+1)
    y = np.zeros(n+1)

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
    x = np.linspace(a, b, n+1)
    y = np.zeros(n+1)

    y[0] = y_0

    for i in range(1, n + 1):
        y[i] = y[i - 1] + h * f(x[i - 1], y[i - 1])

    return x, y


def improved_eulers(a, b, y_0, n, f):
    h = (b - a) / n
    x = np.linspace(a, b, n+1)
    y = np.zeros(n+1)

    y[0] = y_0

    for i in range(1, n + 1):
        k1 = f(x[i - 1], y[i - 1])
        u = y[i - 1] + h * k1
        k2 = f(x[i], u)
        y[i] = y[i - 1] + h * (k1 + k2) / 2

    return x, y

