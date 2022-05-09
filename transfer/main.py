import numpy as np
import matplotlib.pyplot as plt

# Локальные функции
import methods
from methods import leftCorner
from methods import rightCorner
from methods import laxWendroff
from methods import errorMethod

from methods import plotMethod

# Коэффициент в уравнение
a = 2
# Границы
x_1 = 0
x_2 = 1
T = 1

b = 0.5
def phi(x):
  if np.abs(x - 0.4) <= b**2:
    return np.exp(-b**2 / (b**2 - (x - 0.4)**2))
  return 0

def mu(t):
  return 0

# Аналитическое решение
def solution(x, t):
  if np.abs(x - a * t - 0.4) <= b**2:
    return np.exp(-b**2 / (b**2 - (x - a * t - 0.4)**2))
  return 0


N = 200
J = (a + 1) * N

h = (x_2 - x_1) / (N - 1)
tau = T / (J - 1)
current = a * tau / h
print("Current equals " + str(current))
plotMethod(leftCorner, solution, a, phi, mu, x_1, x_2, T, N, J, "левый уголок")
errorMethod(leftCorner, solution, a, phi, mu, x_1, x_2, T, "левый уголок")
plotMethod(rightCorner, solution, a, phi, mu, x_1, x_2, T, N, J, "правый уголок")
errorMethod(rightCorner, solution, a, phi, mu, x_1, x_2, T, "правый уголок")
plotMethod(laxWendroff, solution, a, phi, mu, x_1, x_2, T, N, J, "схема Лакса-Вендроффа")
errorMethod(laxWendroff, solution, a, phi, mu, x_1, x_2, T, "схема Лакса-Вендроффа")

J = N
h = (x_2 - x_1) / (N - 1)
tau = T / (J - 1)
current = a * h / tau

print("Current equals " + str(current))
plotMethod(leftCorner, solution, a, phi, mu, x_1, x_2, T, N, J, "левый уголок")
plotMethod(rightCorner, solution, a, phi, mu, x_1, x_2, T, N, J, "правый уголок")
plotMethod(laxWendroff, solution, a, phi, mu, x_1, x_2, T, N, J, "схема Лакса-Вендроффа")

def phi(x):
  if x <= 0.3:
    return 0.5
  else:
    return 0

def mu(t):
  return 0.5

# Аналитическое решение
def solution(x, t):
  if x - a * t <= 0.3:
    return 0.5
  else:
    return 0

plotMethod(leftCorner, solution, a, phi, mu, x_1, x_2, T, N, J, "левый уголок")
plotMethod(rightCorner, solution, a, phi, mu, x_1, x_2, T, N, J, "правый уголок")
plotMethod(laxWendroff, solution, a, phi, mu, x_1, x_2, T, N, J, "схема Лакса-Вендроффа")

def phi(x):
  if x <= 1 and x >= 0.4:
    return 0.5
  else:
    return 0

def mu(t):
  return 0

# Аналитическое решение
def solution(x, t):
  if x - a * t <= 1 and x - a * t >= 0.4:
    return 0.5
  else:
    return 0

plotMethod(leftCorner, solution, a, phi, mu, x_1, x_2, T, N, J, "левый уголок")
plotMethod(rightCorner, solution, a, phi, mu, x_1, x_2, T, N, J, "правый уголок")
plotMethod(laxWendroff, solution, a, phi, mu, x_1, x_2, T, N, J, "схема Лакса-Вендроффа")

def phi(x):
  if x <= 0.4 and x >= 0.1:
    return 0.5
  else:
    return 0

def mu(t):
  return 0

# Аналитическое решение
def solution(x, t):
  if x - a * t <= 0.4 and x - a * t >= 0.1:
    return 0.5
  else:
    return 0

plotMethod(leftCorner, solution, a, phi, mu, x_1, x_2, T, N, J, "левый уголок")
plotMethod(rightCorner, solution, a, phi, mu, x_1, x_2, T, N, J, "правый уголок")
plotMethod(laxWendroff, solution, a, phi, mu, x_1, x_2, T, N, J, "схема Лакса-Вендроффа")