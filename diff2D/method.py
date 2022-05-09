from operator import eq
import numpy as np
import matplotlib.pyplot as plt

def solveDiff2D(k, v, f, u0y, u1y, ux0, ux1, hx, hy, x0, y0, x1, y1):
    Nx = int((x1 - x0) / hx) + 1
    Ny = int((y1 - y0) / hy) + 1
    length = Nx * Ny
    
    A = np.zeros(shape=(length, length))
    B = np.zeros(length)
    
    u = np.zeros(shape=(Ny, Nx))
    # Дополнительно посчитаем число пекле в каждом узле
    peclet = np.zeros(shape=(Ny, Nx))
    
    def fromGridToArray(i, j):
        return j * Nx + i
    def insertIntoMatrix(i, j, row, coef):
        idx = fromGridToArray(i, j)
        A[row][idx] = coef
    def calculateCoefUp(i, j):
        coef = -k(i * hx, (j + 0.5) * hy) / (hy * hy) + v(i * hx, (j + 0.5) * hy)[1] / (2 * hy)
        #print("coef up:", coef)
        return coef
    def calculateCoefRight(i, j):
        coef = -k((i + 0.5) * hx, j * hy) / (hx * hx) + v((i + 0.5) * hx, j * hy)[0] / (2 * hx)
        #print("coef right:", coef)
        return coef
    def calculateCoefDown(i, j):
        coef = -k(i * hx, (j - 0.5) * hy) / (hy * hy) - v(i * hx, (j - 0.5) * hy)[1] / (2 * hy)
        #print("coef down", coef)
        return coef
    def calculateCoefLeft(i, j):
        coef = -k((i - 0.5) * hx, j * hy) / (hx * hx) - v((i - 0.5) * hx, j * hy)[0] / (2 * hx)
        #print("coef left:", coef)
        return coef
    def calculateCoefMiddle(i, j): 
        coef = (k((i + 0.5) * hx, j * hy) + k((i - 0.5) * hx, j * hy)) / (hx * hx) + (k(i * hx, (j + 0.5) * hy) + k(i * hx, (j - 0.5) * hy)) / (hy * hy) + (v((i + 0.5) * hx, j * hy)[0] - v((i - 0.5) * hx, j * hy)[0]) / (2 * hx) + (v(i * hx, (j + 0.5) * hy)[1] - v(i * hx, (j - 0.5) * hy)[1]) / (2 * hy)
        # print("coef middle:", coef)
        return coef
    def addIntoColumn(row, value):
        B[row] += value
    
    # Какое уравнение (строчку в матрице и столбце) заполняем
    eqCounter = 0
    
    # Заполним граничные условия
    for i in range(0, Nx):
        u[0][i] = ux0(i * hx)
        insertIntoMatrix(i, 0, eqCounter, 1)
        addIntoColumn(eqCounter, u[0][i])
        eqCounter += 1
        
        u[-1][i] = ux1(i * hx)
        insertIntoMatrix(i, Ny - 1, eqCounter, 1)
        addIntoColumn(eqCounter, u[-1][i])
        eqCounter += 1
        
    for j in range(1, Ny - 1):
        u[j][0] = u0y(j * hy)
        insertIntoMatrix(0, j, eqCounter, 1)
        addIntoColumn(eqCounter, u[j][0])
        eqCounter += 1
               
        u[j][-1] = u1y(j * hy)
        insertIntoMatrix(Nx - 1, j, eqCounter, 1)
        addIntoColumn(eqCounter, u[j][-1])
        eqCounter += 1
    
    # Построим матрицу A и столбец B
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            # print("i=", i, "j=", j, "eq=", eqCounter)
            # Middle
            insertIntoMatrix(i, j, eqCounter, calculateCoefMiddle(i, j))
            # Up
            if j + 1 != Ny - 1:
                insertIntoMatrix(i, j + 1, eqCounter, calculateCoefUp(i, j))
            else:
                addIntoColumn(eqCounter, -calculateCoefUp(i, j) * u[j + 1][i])
            # Down
            if j - 1 != 0:
                insertIntoMatrix(i, j - 1, eqCounter, calculateCoefDown(i, j))
            else:
                addIntoColumn(eqCounter, -calculateCoefDown(i, j) * u[j - 1][i])
            # Left
            if i - 1 != 0:
                insertIntoMatrix(i - 1, j, eqCounter, calculateCoefLeft(i, j))
            else:
                addIntoColumn(eqCounter, -calculateCoefLeft(i, j) * u[j][i - 1])
            # Right
            if i + 1 != Nx - 1:
                insertIntoMatrix(i + 1, j, eqCounter, calculateCoefRight(i, j))
            else:
                addIntoColumn(eqCounter, -calculateCoefRight(i, j) * u[j][i + 1])
            
            addIntoColumn(eqCounter, f(i * hx, j * hy))
            eqCounter += 1
    
    # Решаем получившуюся линейную систему
    x = np.linalg.solve(A, B)

    # Переносим решение в сетку
    for j in range(0, Ny):
        for i in range(0, Nx):
            u[j][i] = x[fromGridToArray(i, j)]
            if k(i * hx, j * hy) != 0:
                peclet[j][i] = np.max([hx, hy]) * np.max(v(i * hx, j * hy)) / k(i * hx, j * hy)
            else:
                peclet[j][i] = -1
            
    return u, peclet

def solveDiff2DOptimized(k, v, f, u0y, u1y, ux0, ux1, hx, hy, x0, y0, x1, y1):
    Nx = int((x1 - x0) / hx) + 1
    Ny = int((y1 - y0) / hy) + 1
    length = Nx * Ny
    
    A = np.zeros(shape=(length, length))
    B = np.zeros(length)
    
    u = np.zeros(shape=(Ny, Nx))
    # Дополнительно посчитаем число пекле в каждом узле
    peclet = np.zeros(shape=(Ny, Nx))

    def kOptimized(x, y):
        # Сначала посчитаем Пекле в данном узле
        curPeclet = np.max([hx, hy]) * np.max(v(x, y)) / k(x, y)
        if curPeclet < 2:
            return k(x, y)
        else:
            return k(x, y) * curPeclet / 2

    def fromGridToArray(i, j):
        return j * Nx + i
    def insertIntoMatrix(i, j, row, coef):
        idx = fromGridToArray(i, j)
        A[row][idx] = coef
    def calculateCoefUp(i, j):
        coef = -kOptimized(i * hx, (j + 0.5) * hy) / (hy * hy) + v(i * hx, (j + 0.5) * hy)[1] / (2 * hy)
        #print("coef up:", coef)
        return coef
    def calculateCoefRight(i, j):
        coef = -kOptimized((i + 0.5) * hx, j * hy) / (hx * hx) + v((i + 0.5) * hx, j * hy)[0] / (2 * hx)
        #print("coef right:", coef)
        return coef
    def calculateCoefDown(i, j):
        coef = -kOptimized(i * hx, (j - 0.5) * hy) / (hy * hy) - v(i * hx, (j - 0.5) * hy)[1] / (2 * hy)
        #print("coef down", coef)
        return coef
    def calculateCoefLeft(i, j):
        coef = -kOptimized((i - 0.5) * hx, j * hy) / (hx * hx) - v((i - 0.5) * hx, j * hy)[0] / (2 * hx)
        #print("coef left:", coef)
        return coef
    def calculateCoefMiddle(i, j): 
        coef = (kOptimized((i + 0.5) * hx, j * hy) + kOptimized((i - 0.5) * hx, j * hy)) / (hx * hx) + (kOptimized(i * hx, (j + 0.5) * hy) + kOptimized(i * hx, (j - 0.5) * hy)) / (hy * hy) + (v((i + 0.5) * hx, j * hy)[0] - v((i - 0.5) * hx, j * hy)[0]) / (2 * hx) + (v(i * hx, (j + 0.5) * hy)[1] - v(i * hx, (j - 0.5) * hy)[1]) / (2 * hy)
        # print("coef middle:", coef)
        return coef
    def addIntoColumn(row, value):
        B[row] += value
    
    # Какое уравнение (строчку в матрице и столбце) заполняем
    eqCounter = 0
    
    # Заполним граничные условия
    for i in range(0, Nx):
        u[0][i] = ux0(i * hx)
        insertIntoMatrix(i, 0, eqCounter, 1)
        addIntoColumn(eqCounter, u[0][i])
        eqCounter += 1
        
        u[-1][i] = ux1(i * hx)
        insertIntoMatrix(i, Ny - 1, eqCounter, 1)
        addIntoColumn(eqCounter, u[-1][i])
        eqCounter += 1
        
    for j in range(1, Ny - 1):
        u[j][0] = u0y(j * hy)
        insertIntoMatrix(0, j, eqCounter, 1)
        addIntoColumn(eqCounter, u[j][0])
        eqCounter += 1
               
        u[j][-1] = u1y(j * hy)
        insertIntoMatrix(Nx - 1, j, eqCounter, 1)
        addIntoColumn(eqCounter, u[j][-1])
        eqCounter += 1
    
    # Построим матрицу A и столбец B
    pecletBase = np.max([hx, hy])
    for j in range(1, Ny - 1):
        for i in range(1, Nx - 1):
            # Посчитаем число Пекле у каждого из соседей
            pecletUp = pecletBase * np.max(v(i * hx, (j + 0.5) * hy)) / kOptimized(i * hx, (j + 0.5) * hy)
            pecletDown = pecletBase * np.max(v(i * hx, (j - 0.5) * hy)) / kOptimized(i * hx, (j - 0.5) * hy)
            pecletRight = pecletBase * np.max(v((i + 0.5) * hx, j * hy)) / kOptimized((i + 0.5) * hx, j * hy)
            pecletLeft = pecletBase * np.max(v((i - 0.5) * hx, j * hy)) / kOptimized((i - 0.5) * hx, j * hy)
            # В качестве числа Пекле данного узла возьмем максимальный из соседей
            peclet[j][i] = np.max([pecletUp, pecletDown, pecletRight, pecletLeft])
            
            # print("i=", i, "j=", j, "eq=", eqCounter)
            # Middle
            insertIntoMatrix(i, j, eqCounter, calculateCoefMiddle(i, j))
            # Up
            if j + 1 != Ny - 1:
                insertIntoMatrix(i, j + 1, eqCounter, calculateCoefUp(i, j))
            else:
                addIntoColumn(eqCounter, -calculateCoefUp(i, j) * u[j + 1][i])
            # Down
            if j - 1 != 0:
                insertIntoMatrix(i, j - 1, eqCounter, calculateCoefDown(i, j))
            else:
                addIntoColumn(eqCounter, -calculateCoefDown(i, j) * u[j - 1][i])
            # Left
            if i - 1 != 0:
                insertIntoMatrix(i - 1, j, eqCounter, calculateCoefLeft(i, j))
            else:
                addIntoColumn(eqCounter, -calculateCoefLeft(i, j) * u[j][i - 1])
            # Right
            if i + 1 != Nx - 1:
                insertIntoMatrix(i + 1, j, eqCounter, calculateCoefRight(i, j))
            else:
                addIntoColumn(eqCounter, -calculateCoefRight(i, j) * u[j][i + 1])
            
            addIntoColumn(eqCounter, f(i * hx, j * hy))
            eqCounter += 1
    
    # Решаем получившуюся линейную систему
    x = np.linalg.solve(A, B)

    # Переносим решение в сетку
    for j in range(0, Ny):
        for i in range(0, Nx):
            u[j][i] = x[fromGridToArray(i, j)]
            
    return u, peclet