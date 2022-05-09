import numpy as np
import matplotlib.pyplot as plt

def solveWeighted(sigma, N, J, T, L, phi, psi, mu, f):
	h = L / (N - 1)
	tau = T / (J - 1)

	a = sigma * tau / (h * h)
	b = (1 - sigma) * tau / (h * h)

	# Mesh
	y = np.zeros(shape = (J, N))
	# Boundary conditions
	for n in range(0, N):
		y[0][n] = phi(h * n)
	for j in range(0, J):
		y[j][0] = psi(tau * j)
		y[j][-1] = mu(tau * j)

	# Generate matrix
	A = np.zeros(shape = (N - 2, N - 2))
	A[0][0] = 2 * a + 1
	A[0][1] = -a

	for i in range(1, N - 3):
		A[i][i - 1] = -a
		A[i][i] = 2 * a + 1
		A[i][i + 1] = -a

	A[-1][-2] = -a
	A[-1][-1] = 2 * a + 1

	# Solve line by line
	for j in range(1, J):
		# Generate column
		B = np.zeros(N - 2)
		B[0] = a * y[j][0]
		B[-1] = a * y[j][-1]

		for n in range(1, N - 1):
			value = y[j - 1][n] * (1 - 2 * b) + b * y[j - 1][n + 1] + b * y[j - 1][n - 1] + tau * f(n * h, tau * (j - 1))
			B[n - 1] += value

		x = np.linalg.solve(A, B)
		for n in range(1, N - 1):
			y[j][n] = x[n - 1]

	# Посчитали значение в каждом узле, сгенерируем три массива
	size = N * J
	x_ans = np.zeros(size)
	y_ans = np.zeros(size)
	t_ans = np.zeros(size)

	counter = 0
	for j in range(0, J):
		for n in range(0, N):
			x_ans[counter] = n * h
			t_ans[counter] = j * tau
			y_ans[counter] = y[j][n]
			counter += 1
  	
	return x_ans, t_ans, y_ans

def solveWeightedHighOrder(N, J, T, L, phi, psi, mu, f):
	h = L / (N - 1)
	tau = T / (J - 1)

	sigma = 0.5 - h * h / (12 * tau)

	a = sigma * tau / (h * h)
	b = (1 - sigma) * tau / (h * h)

	# Mesh
	y = np.zeros(shape = (J, N))
	# Boundary conditions
	for n in range(0, N):
		y[0][n] = phi(h * n)
	for j in range(0, J):
		y[j][0] = psi(tau * j)
		y[j][-1] = mu(tau * j)

	# Generate matrix
	A = np.zeros(shape = (N - 2, N - 2))
	A[0][0] = 2 * a + 1
	A[0][1] = -a

	for i in range(1, N - 3):
		A[i][i - 1] = -a
		A[i][i] = 2 * a + 1
		A[i][i + 1] = -a

	A[-1][-2] = -a
	A[-1][-1] = 2 * a + 1

	# Solve line by line
	for j in range(1, J):
		# Generate column
		B = np.zeros(N - 2)
		B[0] = a * y[j][0]
		B[-1] = a * y[j][-1]

		for n in range(1, N - 1):
			func = f(n * h, tau * (j - 1)) + (f((n + 1) * h, tau * (j - 1)) - 2 * f(n * h, tau * (j - 1)) + f((n - 1) * h, tau * (j - 1))) / 24
			value = y[j - 1][n] * (1 - 2 * b) + b * y[j - 1][n + 1] + b * y[j - 1][n - 1] + tau * func
			B[n - 1] += value

		x = np.linalg.solve(A, B)
		for n in range(1, N - 1):
			y[j][n] = x[n - 1]

	# Посчитали значение в каждом узле, сгенерируем три массива
	size = N * J
	x_ans = np.zeros(size)
	y_ans = np.zeros(size)
	t_ans = np.zeros(size)

	counter = 0
	for j in range(0, J):
		for n in range(0, N):
			x_ans[counter] = n * h
			t_ans[counter] = j * tau
			y_ans[counter] = y[j][n]
			counter += 1
  	
	return x_ans, t_ans, y_ans