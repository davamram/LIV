import matplotlib.pyplot as plt
import numpy as np

def plot_and_fit():
    # Les coordonnées des points mesurés
    x_values = np.array([2300, 2250, 2225, 2210, 2200, 2175, 2150])
    k_values = [1 / (1 - ((x ** 2) / (2 * 511e-6 ** 2))) for x in x_values]
    y_values = np.array([0.258267, 0.186887, 0.168647, 0.126563, 0.0177835, 0.00588068, 0.00470678])

    # Ajustement linéaire (fit d'une droite)
    coeffs = np.polyfit(x_values, y_values, 1)

    # Extraire les coefficients de la droite
    slope, intercept = coeffs

    # Trouver la valeur de x lorsque f(x) = 0.05 (y = 0.05)
    x_value_for_f_0_05 = (0.05 - intercept) / slope

    # Tracer les points
    plt.figure(figsize=(8, 6))
    plt.scatter(x_values, y_values, color='blue', marker='o', s=50, label='CLs values')

    # Tracer la droite ajustée
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
    
    # Extraire les coefficients de la droite
    slope, intercept = coeffs

    # Trouver la valeur de x lorsque f(x) = 0.05 (y = 0.05)
    x_value_for_f_0_05 = (0.05 - intercept) / slope

    # Tracer les points
    plt.figure(figsize=(8, 6))
    plt.scatter(k_values, y_values, color='blue', marker='o', s=50, label='CLs values')

    # Tracer la droite ajustée
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

def plot_and_interpolate():
    # Les coordonnées des points mesurés
    x_values = np.array([2500, 2300, 2250, 2225, 2210, 2150])
    y_values = np.array([0.764851, 0.258267, 0.186887, 0.168647, 0.126563, 0.00470678])

    # Créer une liste de valeurs x pour l'interpolation
    x_fit = np.linspace(min(x_values), max(x_values), 100)

    # Interpolation linéaire entre chaque paire de points
    y_fit = np.interp(x_fit, x_values, y_values)

    # Tracer les points mesurés et la courbe interpolée
    plt.figure(figsize=(8, 6))
    plt.scatter(x_values, y_values, color='blue', marker='o', s=50, label='Points mesurés')
    plt.plot(x_fit, y_fit, color='red', label='Courbe interpolée')
    
    plt.xlabel('X Values')
    plt.ylabel('Y Values')
    plt.title('Valeurs Mesurées et Interpolation Linéaire')
    plt.legend()
    plt.grid(True)
    plt.show()

#plot_and_interpolate()
plot_and_fit()
