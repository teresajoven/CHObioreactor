import pandas as pd
import os
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt

ex_data = pd.read_excel('./heatmap_excel.xlsx', sheet_name = "simulations")

Dmin, Dmax, DNsamps = 0.020, 0.0314, 2  # (h-1) Minimal, maximal and number of dilution rates 
X0, Xf, XNsamps = 1.0, 3.0, 5
round_y_axis = []
round_x_axis = []

x_axis_labels = np.linspace(X0, Xf, XNsamps)
y_axis_labels = np.linspace(Dmin, Dmax, DNsamps)

col1 = x_axis_labels
col2 = y_axis_labels

for i in col1:
    label = str(round(i, 2))
    round_x_axis.append(label)

for i in col2:
    label = str(round(i, 4))
    round_y_axis.append(label)
 
x_axis_labels= round_x_axis
y_axis_labels = round_y_axis

N=["GLC","GLN","PHE","ARG","ASN","ASP","LYS","MET","HIS","ILE","LEU","PRO","VAL","SER","THR","TRP","TYR"]

for i in N: 
 col3 = ex_data[i+ '(mM)'].tolist()
 array_flux = np.array(col3)
 good_ar = array_flux.reshape(len(y_axis_labels), len(x_axis_labels))
 fig, ax = plt.subplots(figsize=(10,10))
 sns.heatmap(good_ar, annot=True, xticklabels=x_axis_labels, yticklabels=y_axis_labels, linewidths=.5, ax=ax, cmap = "RdYlBu_r", fmt='.3g')
 ax.set_title('Concentration of '+i+' in tank (g L-1)', size = 18)
 ax.invert_yaxis()

 plt.xlabel('Biomass (gDW L-1)', fontsize = 15)
 plt.ylabel('Dilution rate (h-1)', fontsize = 15)
 plt.yticks(rotation=360)
 plt.show() 


