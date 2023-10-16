import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy.optimize import bisect
import os


def getCls():
    path = "results/"
    files = os.listdir(path)
    cls_values=[]
    energie_values=[]
    for file in files:
        with open(path+file, 'r') as f:
            lines = f.readlines()
            
        line_12 = lines[11]
        line_4 = lines[3]
        cls_values.append(float(line_12[13:]))
        energie_values.append(float(line_4[19:23]))
    return([energie_values, cls_values])

def plot_and_fit(x_values, k_values, y_values):

    # Fit linéaire
    coeffs = np.polyfit(x_values, y_values, 1)

    # Pour obtenir les coeffs de la droite
    slope, intercept = coeffs

    # Trouver la valeur de x lorsque CLs = 0.05
    x_value_for_f_0_05 = (0.05 - intercept) / slope

    # Tracer
    plt.figure(figsize=(8, 6))
    plt.scatter(x_values, y_values, color='blue', marker='o', s=50, label='CLs values')

    # Tracer la droite
    x_fit = np.linspace(min(x_values), max(x_values), 100)
    y_fit = np.polyval(coeffs, x_fit)
    plt.plot(x_fit, y_fit, color='red', label='Linear fit')

    plt.xlabel(r'$E_{tr}$')
    plt.ylabel('CLs')
    plt.title('Limit on Etr')
    plt.legend()
    plt.grid(True)
    plt.axhline(y=0.05, color='gray', linestyle='--')
    plt.axvline(x=x_value_for_f_0_05, color='green', linestyle='--')
    plt.legend()
    plt.savefig('plot_fit_E.png')
    plt.show()

    coeffs = np.polyfit(k_values, y_values, 1)
    
    # Pour obtenir les coefficients de la droite
    slope, intercept = coeffs

    # Trouver la valeur de x lorsque CLs = 0.05
    x_value_for_f_0_05 = (0.05 - intercept) / slope

    # Tracer
    plt.figure(figsize=(8, 6))
    plt.scatter(k_values, y_values, color='blue', marker='o', s=50, label='CLs values')

    # Tracer la droite
    k_fit = np.linspace(min(k_values), max(k_values), 100)
    y_fit = np.polyval(coeffs, k_fit)
    plt.plot(k_fit, y_fit, color='red', label='Linear fit')

    plt.xlabel(r'$\kappa$')
    plt.ylabel('CLs')
    plt.title('Limit on kappa')
    plt.legend()
    #plt.grid(True)
    plt.axhline(y=0.05, color='gray', linestyle='--')
    plt.axvline(x=x_value_for_f_0_05, color='green', linestyle='--')
    plt.text(-1.135e-13, 0.06, r'$\kappa > {:.3e}$'.format(x_value_for_f_0_05), fontsize=12, color='purple', fontweight='bold')
    plt.legend()
    plt.savefig('plot_fit_K.png')
    plt.show()
    # Show the first cls point beneath 0.05
    for x, y, k in zip(x_values, y_values, k_values):
        if y < 0.05:
            print("For CLs =", y, "(first point below 0.05), we have Eth =", x, "and kappa =", k)
            break

def plot_and_interpolate(x, y):
    interpolation = interpolate.interp1d(x, y, kind='linear')
    x_new = np.linspace(x.min(), x.max(), 100)
    k = np.array([1 / (1 - ((xi ** 2) / (2 * 511e-6 ** 2))) for xi in x])
    k_new=np.linspace(k.min(), k.max(), 100)
    y_new = interpolation(x_new)

    def find_x(x):
        return interpolation(x) - 0.05
    
    #For E

    x_solution = bisect(find_x, x.min(), x.max())

    plt.plot(x, y, 'o', label='CLs values')
    plt.plot(x_new, y_new, '-', label='Interpolation')
    plt.xlabel(r'$E_{tr}$')
    plt.ylabel('CLs')
    plt.title('Limit on Etr')
    plt.legend()
    plt.axhline(y=0.05, color='gray', linestyle='--')
    plt.axvline(x=x_solution, color='green', linestyle='--')

    plt.text(2100, 0.06, r'$E > {:.3e}$'.format(x_solution), fontsize=12, color='purple', fontweight='bold')
    plt.legend()
    plt.savefig('plot_interpol_E.png')
    plt.show()

    # For Kappa
    interpolation = interpolate.interp1d(k, y, kind='linear')
    y_new = interpolation(k_new)

    k_solution = 1 / (1 - ((x_solution ** 2) / (2 * 511e-6 ** 2)))

    plt.plot(k, y, 'o', label='CLs values')
    plt.plot(k_new, y_new, '-', label='Interpolation')
    plt.xlabel(r'$\kappa$')
    plt.ylabel('CLs')
    plt.title('Limit on kappa')
    plt.legend()
    plt.axhline(y=0.05, color='gray', linestyle='--')
    plt.axvline(x=k_solution, color='green', linestyle='--')

    plt.text(-1.305e-13, 0.06, r'$\kappa > {:.3e}$'.format(k_solution), fontsize=12, color='purple', fontweight='bold')

    #plt.text(2200, 0.06, r'$E_{tr} > {:.3e}$'.format(x_solution), fontsize=12, color='purple', fontweight='bold')
    plt.legend()
    plt.savefig('plot_interpol_K.png')
    plt.show()




# x_values = np.array([2300, 2250, 2225, 2210, 2200, 2175, 2150]) # First values
x_values = np.array([2400, 2300, 2250, 2237, 2230, 2225, 2200, 2100]) # Second values (madgraph)
x_values = np.array([2400, 2300, 2250, 2237, 2230, 2225, 2200, 2100, 2050, 2025, 2000]) # Third values (Data uncertainties)
x_values = np.array(getCls()[0])

# y_values = np.array([0.258267, 0.186887, 0.168647, 0.126563, 0.0177835, 0.00588068, 0.00470678]) # First presentation
# y_values = np.array([0.176181, 0.154514, 0.151199, 0.11079, 0.0112652, 0.00278396, 0.00165911]) # First (sherpa 900GeV ?) velue with poisson continuous
y_values = np.array([0.374682, 0.129847, 0.0634117, 0.050172, 0.045412, 0.0425831, 0.0251836, 0.000102955]) # Second values (madgraph 120/1000/2000)
y_values = np.array([0.390181, 0.151625, 0.0928074, 0.0662857, 0.0552424, 0.048478, 0.0281215, 0.0020202, 0.00102249, 0.00104254, 0]) # Third values (Data uncertainties)
y_values = np.array(getCls()[1])
plot_and_interpolate(x_values, y_values)
