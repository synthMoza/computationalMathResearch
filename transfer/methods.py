import numpy as np
import matplotlib.pyplot as plt
import random

from mpl_toolkits import mplot3d
from matplotlib.widgets import Slider, Button

# Левый уголок
# N, J - кол-во отсчетов x и t соответственно
#          (n,j + 1)  *
#                     |
#      (n, j - 1) *---* (n, j)
def leftCorner(a, phi, mu, x_1, x_2, T, N, J):
  h = (x_2 - x_1) / (N - 1)
  tau = T / (J - 1)
  # Инициализируем ответ
  y = np.zeros(shape = (J, N))
  for i in range(0, N):
    y[0][i] = phi(x_1 + h * i)
  for i in range(0, J):
    y[i][0] = mu(tau * i)

  for i in range(1, J):
    for n in range(1, N):
      y[i][n] = y[i - 1][n] - (a * tau / h) * (y[i - 1][n] - y[i - 1][n - 1])

  # Посчитали значение в каждом узле, сгенерируем три массива
  size = N * J
  x_ans = np.zeros(size)
  y_ans = np.zeros(size)
  t_ans = np.zeros(size)

  for i in range(0, size):
    t_index = int(i / J)
    x_index = int(i % N)

    y_ans[i] = y[t_index][x_index]
    x_ans[i] = x_1 + h * x_index
    t_ans[i] = tau * t_index
  
  return x_ans, t_ans, y_ans

# Правый уголок
# N, J - кол-во отсчетов x и t соответственно
#        *---* (n,j + 1)
#            |
# (n - 1, j) * (n, j)
def rightCorner(a, phi, mu, x_1, x_2, T, N, J):
  h = (x_2 - x_1) / (N - 1)
  tau = T / (J - 1)
  # Инициализируем ответ
  y = np.zeros(shape = (J, N))
  for i in range(0, N):
    y[0][i] = phi(x_1 + h * i)
  for i in range(0, J):
    y[i][0] = mu(tau * i)

  for i in range(1, J):
    for n in range(1, N):
      y[i][n] = (h * y[i - 1][n] + a * tau * y[i][n - 1]) / (h + a * tau)

  # Посчитали значение в каждом узле, сгенерируем три массива
  size = N * J
  x_ans = np.zeros(size)
  y_ans = np.zeros(size)
  t_ans = np.zeros(size)

  for i in range(0, size):
    t_index = int(i / J)
    x_index = int(i % N)

    y_ans[i] = y[t_index][x_index]
    x_ans[i] = x_1 + h * x_index
    t_ans[i] = tau * t_index
  
  return x_ans, t_ans, y_ans

# Лакса - Вендроффа
# N, J - кол-во отсчетов x и t соответственно
#                 *(n,j + 1)
#                 |
# (n - 1, j)  *---*---* (n + 1, j)
#               (n, j)
def laxWendroff(a, phi, mu, x_1, x_2, T, N, J):
  h = (x_2 - x_1) / (N - 1)
  tau = T / (J - 1)
  # Инициализируем ответ
  y = np.zeros(shape = (J + N - 2, J + N - 2))
  for i in range(0, J + N - 2):
    y[0][i] = phi(x_1 + h * i)
  for i in range(0, J + N - 2):
    y[i][0] = mu(tau * i)

  for i in range(1, J):
    for n in range(1, J + N - 2 - i):
      y[i][n] = (a * a * tau * tau / (2 * h * h) - a * tau / (2 * h)) * y[i - 1][n + 1] + (1 - a * a * tau * tau / (h * h)) * y[i - 1][n] + (a * tau / (2 * h) + a * a * tau * tau / (2 * h * h)) * y[i - 1][n - 1]

  # Посчитали значение в каждом узле, сгенерируем три массива
  size = N * J
  x_ans = np.zeros(size)
  y_ans = np.zeros(size)
  t_ans = np.zeros(size)

  for i in range(0, size):
    t_index = int(i / J)
    x_index = int(i % N)

    y_ans[i] = y[t_index][x_index]
    x_ans[i] = x_1 + h * x_index
    t_ans[i] = tau * t_index
  
  return x_ans, t_ans, y_ans

def getOrder(name):
    if name[0] == "с":
      return 2.1
    else:
      return 0.9

# Функция для построения графика численного и аналитического решений с ползунков
def plotMethod(method, solution, a, phi, mu, x_1, x_2, T, N, J, name):
  x, t, y = method(a, phi, mu, x_1, x_2, T, N, J)

  plt.figure(figsize=(11.7,8.3))
  plt.title("Сравнение аналитического и численного решений, " + name)
  plt.grid(which='both')
  plt.grid(which='minor', alpha=0.2)
  plt.grid(which='major', alpha=0.5)
  plt.minorticks_on()
  plt.autoscale()
  plt.xlabel("$x$", fontsize=10)
  plt.ylabel("$y(x)$", fontsize=10)
  # Изначально строим при t = 0
  lineNum, = plt.plot(x[0 : N - 1], y[0 : N - 1], 'r', label='численное решение')
  trueY = [solution(x[i], t[i]) for i in range(0, N)]
  lineTrue, = plt.plot(x[0 : N], trueY, 'g--', label='аналитическое решение')
  plt.legend()

  # Сдвинуть график так, чтобы ползунок поместился
  plt.subplots_adjust(bottom = 0.25)

  # Ползунок для времени
  axtime = plt.axes([0.25, 0.1, 0.65, 0.03])
  timeSlider = Slider(
      ax = axtime,
      label = 'Время (в отсчетах)',
      valmin = 0,
      valmax = J - 1,
      valstep = 1,
      valinit = 0,
  )

  def update(val):
    lineNum.set_xdata(x[val * N : (val + 1) * N])
    lineNum.set_ydata(y[val * N : (val + 1) * N])
    lineTrue.set_xdata(x[val * N : (val + 1) * N])
    trueY = [solution(x[i], t[i]) for i in range(val * N, (val + 1) * N)]
    lineTrue.set_ydata(trueY)

  timeSlider.on_changed(update)
  plt.show()

def errorMethod(method, solution, a, phi, mu, x_1, x_2, T, name):
  # h - var, tau - const
  N_start = 100
  N_end = 150
  N_step = 2
  N_amount = int((N_end - N_start) / N_step)
  error = np.zeros(N_amount)
  h = np.zeros(N_amount)
  J = a * N_end

  for i in range(0, N_amount):
    N = N_start + N_step * i
    h[i] = (x_2 - x_1) / (N - 1)
    x, t, y = method(a, phi, mu, x_1, x_2, T, N, J)
    error[i] = np.max([np.abs(y[j] - solution(x[j], t[j])) for j in range(0, N * J)])

  error = [error[0] / h[0] * np.power(h[i], getOrder(name)) for i in range(0, N_amount)]
  plt.figure(figsize=(11.7,8.3))
  plt.title("Зависимость ошибки от шага при $\\tau = const$, " + name)
  plt.grid(which='both')
  plt.grid(which='minor', alpha=0.2)
  plt.grid(which='major', alpha=0.5)
  plt.minorticks_on()
  plt.autoscale()
  plt.xlabel("$ln(h)$", fontsize=10)
  plt.ylabel("$ln(error)$", fontsize=10)
  # Коэффициент
  coef = np.polyfit(np.log(h), np.log(error), 1)
  plt.plot(np.log(h), np.log(error), 'r', label='k = ' + str(round(coef[0], 1)))
  plt.legend()
  plt.show()