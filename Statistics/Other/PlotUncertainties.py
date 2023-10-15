import matplotlib.pyplot as plt
from math import sqrt

# Données
Atlas_down = [-7.827176e-02, -3.551370e-02, -1.697738e-02, -7.084197e-03, -2.575804e-03, -1.157590e-03, -5.549091e-04, -2.714763e-04, -1.185141e-04, -5.092944e-05, -2.129382e-05, -8.701203e-06, -2.891808e-06, -5.895342e-07, -7.166077e-08, -8.281072e-09]
Atlas_up = [5.233613e-02, 2.433757e-02, 1.214272e-02, 5.235209e-03, 1.963962e-03, 8.571483e-04, 4.382003e-04, 2.079510e-04, 9.104776e-05, 3.976086e-05, 1.729810e-05, 7.207799e-06, 2.486231e-06, 5.904051e-07, 7.280856e-08, 9.336763e-09]
Bins_values = [1.375000e+02, 1.625000e+02, 1.875000e+02, 2.250000e+02, 2.750000e+02, 3.250000e+02, 3.750000e+02, 4.350000e+02, 5.100000e+02, 6.000000e+02, 7.000000e+02, 8.250000e+02, 1.000000e+03, 1.300000e+03, 1.750000e+03, 2.250000e+03]
Bins_values_stairs = [1.250000e+02, 1.500000e+02, 1.750000e+02, 2.000000e+02, 2.500000e+02, 3.000000e+02, 3.500000e+02, 4.000000e+02, 4.700000e+02, 5.500000e+02, 6.500000e+02, 7.500000e+02, 9.000000e+02, 1.100000e+03, 1.500000e+03, 2.000000e+03, 2.500000e+03]
Atlas_stat = [0.00295715, 0.00177455, 0.000834419, 0.000268998, 0.000112464, 6.68152e-05, 3.31831e-05, 1.32362e-05, 6.20673e-06, 2.69945e-06, 1.42323e-06, 5.45151e-07, 1.94116e-07, 3.17007e-08, 4.10542e-09, 5.15698e-10]
Local_stat = [7.074313e-02, 4.497576e-02, 2.768121e-02, 1.581030e-02, 7.730480e-03, 7.547453e-03, 3.877605e-03, 2.543693e-03, 1.195816e-03, 8.785437e-04, 1.057928e-03, 3.362813e-04, 5.268441e-05, 1.021554e-07, 2.813266e-08, 2.251434e-09]
Down = [-sqrt(Adown**2 - Astat**2 + Lstat**2) for Adown, Astat, Lstat in zip(Atlas_down, Atlas_stat, Local_stat)]
Up = [sqrt(Aup**2 - Astat**2 + Lstat**2) for Aup, Astat, Lstat in zip(Atlas_up, Atlas_stat, Local_stat)]
Values = [2.224240e+00, 1.017430e+00, 5.141328e-01, 2.227860e-01, 8.262708e-02, 3.580690e-02, 1.729010e-02, 8.139827e-03, 3.521509e-03, 1.464270e-03, 6.108828e-04, 2.379368e-04, 7.349850e-05, 1.409705e-05, 1.395380e-06, 9.125426e-08]

Ratio_up = [up/val for up, val in zip(Up, Values)]
Ratio_down = [down/val for down, val in zip(Down, Values)]
Ratio_stat_u = [stat/val for stat, val in zip(Local_stat, Values)]
Ratio_stat_d = [-stat/val for stat, val in zip(Local_stat, Values)]
Ratio_Atlas_u = [up/val for up, val in zip(Atlas_up, Values)]
Ratio_Atlas_d = [up/val for up, val in zip(Atlas_down, Values)]

# Tracer un graphique à barres sans remplissage de couleur
plt.stairs(Ratio_stat_u[-3:], Bins_values_stairs[-4:], color='red', edgecolor='red', label = "Stat")
plt.stairs(Ratio_stat_d[-3:], Bins_values_stairs[-4:], color='red', edgecolor='red')
plt.stairs(Ratio_down[-3:], Bins_values_stairs[-4:], color='none', edgecolor='blue', label = "Stat+Sys")
plt.stairs(Ratio_up[-3:], Bins_values_stairs[-4:], color='none', edgecolor='blue')
plt.stairs(Ratio_Atlas_u[-3:], Bins_values_stairs[-4:], color='none', edgecolor='green', label = "Atlas stat+sys")
plt.stairs(Ratio_Atlas_d[-3:], Bins_values_stairs[-4:], color='none', edgecolor='green')
plt.xscale('log')
plt.xlabel('Energy')
plt.ylabel('Uncertainties')
plt.title('Stat and Syst')
plt.grid(True)
plt.legend()
plt.show()
