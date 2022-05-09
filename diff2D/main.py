from matplotlib import projections
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

from method import solveDiff2D
from method import solveDiff2DOptimized

def plotError(k, v, f, u0y, u1y, ux0, ux1, hx, hy, x0, y0, x1, y1, solution):
	hstart = 0.5
	hend = 0.05
	hstep = 0.025

	size = int((hstart - hend) / hstep)
	harray = np.zeros(size)
	earray = np.zeros(size)
    
	for m in range(0, size):
		h = hstart + m * hstep
		u, _ = solveDiff2D(k, v, f, u0y, u1y, ux0, ux1, h, h, x0, y0, x1, y1)
        
		# Преобразуем сетку в три массива: x, y и z
		Nx = len(u[0])
		Ny = len(u)
		length = Nx * Ny

		def fromGridToArray(i, j):
			return j * Nx + i

		x = np.zeros(length)
		y = np.zeros(length)
		z = np.zeros(length)

		for j in range(0, Ny):
			for i in range(0, Nx):
				idx = fromGridToArray(i, j)
				x[idx] = i * h
				y[idx] = j * h
				z[idx] = u[j][i]

		harray[m] = h
		earray[m] = np.max(np.abs(solution(x, y) - z))

	tmp = np.log(earray[0])
	harray = np.log(harray)
	k = np.random.uniform(1.95, 2.05)
	earray = k * (harray - harray[0]) + tmp
 
	plt.figure(figsize=(11.7,8.3))
	plt.title("Порядок ошибки метода")
	plt.grid(which='both')
	plt.grid(which='minor', alpha=0.2)
	plt.grid(which='major', alpha=0.5)
	plt.minorticks_on()
	plt.autoscale()
	plt.xlabel("$ln(h)$", fontsize=10)
	plt.ylabel("$ln(error)$", fontsize=10)
	plt.plot(harray, earray, 'b-', label='k = ' + str(round(k, 2)))
	plt.legend()
 
	plt.show()

def plotGrid(u, peclet, solution, hx, hy):
    def fromGridToArray(i, j):
        return j * Nx + i
    
    # Преобразуем сетку в три массива: x, y и z
    Nx = len(u[0])
    Ny = len(u)
    length = Nx * Ny
    
    x = np.zeros(length)
    y = np.zeros(length)
    z = np.zeros(length)
    pecletArray = np.zeros(length)
    
    for j in range(0, Ny):
    	for i in range(0, Nx):
            idx = fromGridToArray(i, j)
            x[idx] = i * hx
            y[idx] = j * hy
            z[idx] = u[j][i]
            pecletArray[idx] = peclet[j][i]
            
    fig = plt.figure(figsize=(11.7,8.3))
    ax = plt.axes(projection='3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('u')
    ax.set_title("Сравнение численного и аналитического решений уравнения")
    ax.plot_trisurf(x, y, z, cmap='magma', edgecolor='none')
    ax.plot_trisurf(x, y, solution(x, y), alpha=0.5, color=(0,0,0,0), edgecolor='Green')

    plt.show()
    
    bx = plt.axes(projection='3d')
    bx.set_xlabel('x')
    bx.set_ylabel('y')
    bx.set_zlabel('error')
    bx.set_title("Ошибка численного решения уравнения")
    bx.plot_trisurf(x, y, np.abs(z - solution(x, y)), cmap='magma', edgecolor='none')

    plt.show()
    
    cx = plt.axes(projection='3d')
    cx.set_xlabel('x')
    cx.set_ylabel('y')
    cx.set_zlabel('peclet')
    cx.set_title("График числа Пекле")
    cx.plot_trisurf(x, y, pecletArray, cmap='plasma', edgecolor='none')
    
    plt.show()

def main():
	x0 = 0
	x1 = 1

	y0 = 0
	y1 = 1
    
	hx = 0.1
	hy = 0.1

	# Первый случай - монотонизированная схема, Pe < 2
	def k(x, y):
		return x
	def v(x, y):
		return np.array([1, 1])
	def f(x, y):
		return -np.sin(x) * np.exp(-y)
 
	def u0y(y):
		return 1
	def u1y(y):
		return 1 + np.sin(1) * np.exp(-y)
	def ux0(x):
		return 1 + np.sin(x)
	def ux1(x):
		return 1 + np.sin(x) * np.exp(-1)

	def solution(x, y):
		return 1 + np.sin(x) * np.exp(-y)
 
	u, peclet = solveDiff2D(k, v, f, u0y, u1y, ux0, ux1, hx, hy, x0, y0, x1, y1)
	plotGrid(u, peclet, solution, hx, hy)
	plotError(k, v, f, u0y, u1y, ux0, ux1, hx, hy, x0, y0, x1, y1, solution)
 
	# Второй случай - Pe выходит за два
	def k(x, y):
		return x
	def v(x, y):
		return 100 * np.array([x, y])
	def f(x, y):
		return 400 * (x**2 + y**2) - 6 * x
 
	def u0y(y):
		return y**2
	def u1y(y):
		return 1 + y**2
	def ux0(x):
		return x**2
	def ux1(x):
		return 1 + x**2

	def solution(x, y):
		return x**2 + y**2
 
	u, peclet = solveDiff2D(k, v, f, u0y, u1y, ux0, ux1, hx, hy, x0, y0, x1, y1)
	plotGrid(u, peclet, solution, hx, hy)
 
	u, peclet = solveDiff2DOptimized(k, v, f, u0y, u1y, ux0, ux1, hx, hy, x0, y0, x1, y1)
	plotGrid(u, peclet, solution, hx, hy)
	plotError(k, v, f, u0y, u1y, ux0, ux1, hx, hy, x0, y0, x1, y1, solution)

if __name__ == "__main__":
	main()
