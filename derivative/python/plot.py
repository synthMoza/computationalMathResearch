import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import numpy as np
import csv

# Original function and its derivatives for plotting
def f(x):
    return np.sin(x) * x
def f1d(x):
    return np.sin(x) + x * np.cos(x)
def f2d(x):
    return 2 * np.cos(x) - np.sin(x) * x

# Initialize arrays for filling with data from csv
steps = np.array([])
r1 = np.array([])
r2 = np.array([])
r3 = np.array([])

deviation_file_path = "../deviation.csv"
first_derivative_file_path = "../first_derivative.csv"
second_derivative_file_path = "../second_derivative.csv"

# Read deviation
with open(deviation_file_path) as File:
    reader = csv.reader(File, delimiter='\t')
    for row in reader:
        steps = np.append(steps, float(row[0]))
        r1 = np.append(r1, float(row[1]))
        r2 = np.append(r2, float(row[2]))
        r3 = np.append(r3, float(row[3]))

# Read first derivative
x_f1d = np.array([])
y_f1d = np.array([])

with open(first_derivative_file_path) as File:
    reader = csv.reader(File, delimiter='\t')
    for row in reader:
        x_f1d = np.append(x_f1d, float(row[0]))
        y_f1d = np.append(y_f1d, float(row[1]))

# Read second derivative
x_f2d = np.array([])
y_f2d = np.array([])

with open(second_derivative_file_path) as File:
    reader = csv.reader(File, delimiter='\t')
    for row in reader:
        x_f2d = np.append(x_f2d, float(row[0]))
        y_f2d = np.append(y_f2d, float(row[1]))

# Polyfit deviation in log scale
r1_coef = np.polyfit(np.log(steps), np.log(r1), 1)
r2_coef = np.polyfit(np.log(steps), np.log(r2), 1)
r3_coef = np.polyfit(np.log(steps), np.log(r3), 1)

# Plot points and polyfitted lines
plt.figure(figsize=(8,6))

plt.subplot()
plt.grid()
plt.minorticks_on()
plt.autoscale()
plt.xlabel("$h$, step")
plt.ylabel("$r$, deviation")
plt.title("Deviation-step dependency on derivatives in log scale")
plt.loglog(steps, r1, 'ro', label = "first derivative, first order deviation, k = " + str(round(r1_coef[0], 3)))
plt.loglog(steps, r2, 'go', label = "first derivative, second order deviation, k = " + str(round(r2_coef[0], 3)))
plt.loglog(steps, r3, 'bo', label = "second derivative, second order deviation, k = " + str(round(r3_coef[0], 3)))
plt.legend()
plt.tight_layout()
plt.savefig('deviation.png')

plt.clf()
# Plot original function
plt.figure(figsize=(24,6))
plt.grid()
plt.minorticks_on()
plt.autoscale()
plt.xlabel("$x$, argument")
plt.ylabel("$y$, value")
plt.title("Original function")
y = f(x_f1d)
plt.plot(x_f1d, y, 'b-', label = "$y=sin(x)*x$")
plt.tight_layout()
plt.legend()
plt.savefig('original_function.png')

plt.clf()
# Plot first derivative
plt.figure(figsize=(24,6))
plt.grid()
plt.minorticks_on()
plt.autoscale()
plt.xlabel("$x$, argument")
plt.ylabel("$y$, value")
plt.title("Numerical first derivative compared to the real one")
plt.plot(x_f1d, y_f1d, 'go', label="numerical derivative")

y = f1d(x_f1d)
plt.plot(x_f1d, y, 'r-', label = "real derivative")
plt.legend()
plt.tight_layout()
plt.savefig('first_derivative.png')

plt.clf()
# Plot second derivative
plt.figure(figsize=(24,6))
plt.grid()
plt.minorticks_on()
plt.autoscale()
plt.xlabel("$x$, argument")
plt.ylabel("$y$, value")
plt.title("Numerical second derivative compared to the real one")
plt.plot(x_f2d, y_f2d, 'go', label="numerical derivative")

y = f2d(x_f2d)
plt.plot(x_f2d, y, 'r-', label = "real derivative")
plt.legend()
plt.tight_layout()
plt.savefig('second_derivative.png')