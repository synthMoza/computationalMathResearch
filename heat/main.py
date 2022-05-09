import numpy as np
import matplotlib.pyplot as plt

from method import solveWeighted
from method import solveWeightedHighOrder
from mpl_toolkits import mplot3d
from matplotlib.widgets import Slider, Button

# Equation parameters
a = 1

def f(x, t):
	return x

# Initial conditions
L = 1
T = 1

# u(x, 0) = phi(x)
# u(0, t) = psi(t)
# u(L, t) = mu(t)

def phi(x):
	return np.sin(3 * np.pi * x / 2)
def psi(t):
	return 0
def mu(t):
	return t - np.exp(-(3 * np.pi / 2)**2 * t)

# Solution
def solution(x, t):
	return x * t + np.exp(-(3 * np.pi / 2)**2 * t) * np.sin(3 * np.pi * x / 2)

def plot(sigma, N, J):
	h = L / (N - 1)
	tau = T / (J - 1)

	x, t, y = solveWeighted(sigma, N, J, T, L, phi, psi, mu, f)

	plt.figure(figsize=(11.7,8.3))
	plt.grid(which='both')
	plt.grid(which='minor', alpha=0.2)
	plt.grid(which='major', alpha=0.5)
	plt.title("$N = " + str(N) + ", J = " + str(J) + ", \\sigma = " + str(sigma) + "$")
	plt.minorticks_on()
	plt.autoscale()
	plt.xlabel("$x$", fontsize=10)
	plt.ylabel("$y(x)$", fontsize=10)

	valStart = 0
	start = valStart * N
	end = (valStart + 1) * N

	# Изначально строим при t = 0
	lineNum, = plt.plot(x[start : end], y[start : end], 'r', label='численное решение')

	trueY = [solution(x[i], tau * valStart) for i in range(start, end)]
	lineTrue, = plt.plot(x[start : end], trueY, 'g--', label='аналитическое решение')
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
		start = val * N
		end = start + N

		lineNum.set_xdata(x[start : end])
		lineNum.set_ydata(y[start : end])
		
		lineTrue.set_xdata(x[start : end])
		trueY = [solution(x[i], val * tau) for i in range(start, end)]
		lineTrue.set_ydata(trueY)

	timeSlider.on_changed(update)
	plt.show()

def plotError(sigma, order_h, order_tau):
	N_start = 4
	N_end = 60
	N_step = 1

	length = int((N_end - N_start) / N_step)
	error = np.zeros(length)
	h = np.zeros(length)
	tau = np.zeros(length)

	countError = 0
	for N in range(N_start, N_end, N_step):
		J = 2 * N * N
		x, t, y = solveWeighted(sigma, N, J, T, L, phi, psi, mu, f)

		realY = solution(x, t)
		error[countError] = np.max(np.abs(y - realY))
		h[countError] = L / (N - 1)
		tau[countError] = T / (J - 1)

		countError += 1

	# Plot error
	plt.figure(figsize=(11.7,8.3))
	plt.grid(which='both')
	plt.grid(which='minor', alpha=0.2)
	plt.grid(which='major', alpha=0.5)
	plt.title("Error in log scale, $\\sigma = " + str(sigma) + "$")
	plt.minorticks_on()
	plt.autoscale()
	tmp = h**order_h + tau**order_tau
	k = np.polyfit(np.log(tmp), np.log(error), 1)
	plt.plot(np.log(tmp), np.log(error), 'r-', label='k = ' + str(round(k[0], 2)))
	plt.xlabel("$ln(h^" +str(order_h) + " + \\tau^" + str(order_tau) + ")$", fontsize=10)
	plt.ylabel("$ln(error)$", fontsize=10)
	plt.legend()
	plt.show()


def plotHighOrder(N, J):
	h = L / (N - 1)
	tau = T / (J - 1)

	x, t, y = solveWeightedHighOrder(N, J, T, L, phi, psi, mu, f)

	plt.figure(figsize=(11.7,8.3))
	plt.grid(which='both')
	plt.grid(which='minor', alpha=0.2)
	plt.grid(which='major', alpha=0.5)
	plt.title("$N = " + str(N) + ", J = " + str(J) + ", \\sigma = \\sigma^*$")
	plt.minorticks_on()
	plt.autoscale()
	plt.xlabel("$x$", fontsize=10)
	plt.ylabel("$y(x)$", fontsize=10)

	valStart = 0
	start = valStart * N
	end = (valStart + 1) * N

	# Изначально строим при t = 0
	lineNum, = plt.plot(x[start : end], y[start : end], 'r', label='численное решение')

	trueY = [solution(x[i], tau * valStart) for i in range(start, end)]
	lineTrue, = plt.plot(x[start : end], trueY, 'g--', label='аналитическое решение')
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
		start = val * N
		end = start + N

		lineNum.set_xdata(x[start : end])
		lineNum.set_ydata(y[start : end])
		
		lineTrue.set_xdata(x[start : end])
		trueY = [solution(x[i], val * tau) for i in range(start, end)]
		lineTrue.set_ydata(trueY)

	timeSlider.on_changed(update)
	plt.show()

def plotHighOrderError():
	N_start = 6
	N_end = 40
	N_step = 1

	order_h = 4
	order_tau = 2

	length = int((N_end - N_start) / N_step)
	error = np.zeros(length)
	h = np.zeros(length)
	tau = np.zeros(length)

	countError = 0
	for N in range(N_start, N_end, N_step):
		J = 2 * N * N + 5
		x, t, y = solveWeightedHighOrder(N, J, T, L, phi, psi, mu, f)

		realY = solution(x, t)
		error[countError] = np.max(np.abs(y - realY))
		h[countError] = L / (N - 1)
		tau[countError] = T / (J - 1)

		countError += 1

	# Plot error
	plt.figure(figsize=(11.7,8.3))
	plt.grid(which='both')
	plt.grid(which='minor', alpha=0.2)
	plt.grid(which='major', alpha=0.5)
	plt.title("Error in log scale, $\\sigma = \\sigma^*$")
	plt.minorticks_on()
	plt.autoscale()
	k = np.polyfit(np.log(h**order_h + tau**order_tau), np.log(error), 1)
	plt.plot(np.log(h**order_h + tau**order_tau), np.log(error), 'r-', label='k = ' + str(round(k[0], 2)))
	plt.xlabel("$ln(h^" +str(order_h) + " + \\tau^" + str(order_tau) + ")$", fontsize=10)
	plt.ylabel("$ln(error)$", fontsize=10)
	plt.legend()
	plt.show()

# ==========================================================================

# Mesh parameters
N = 5
J = 2 * (N * N) + 5

plot(0, N, N)
plot(0, N, J)
plot(1, N, J)
plot(0.5, N, J)
plotHighOrder(N, J)

# sigma = 0, O(h^2 + tau)
plotError(sigma=0, order_h=2, order_tau=1)
# sigma = 1, O(h^2 + tau)
plotError(sigma=1, order_h=2, order_tau=1)
# sigma = 0.5, O(h^2 + tau^2)
plotError(sigma=0.5, order_h=2, order_tau=2)
# sigma = sigma*, O(h^4 + tau^2)
plotHighOrderError()