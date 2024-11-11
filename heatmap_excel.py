#!/usr/bin/python
# Transform cobra model to Flexible Net (FN)
from __future__ import division, print_function
import numpy as np
import cobra.test
from CHOBiorFN_funcion import comProductivity, loadCHOmodel
import xlwt
import time


# Measure time
start = time.time()

# Creating xls file
row = 0
xlsfile = "heatmap_excel.xlsx"
wb = xlwt.Workbook()
bstyle = xlwt.easyxf('font: bold on')
ws = wb.add_sheet("simulations")

N=["GLC"] #concatenate before GLC (glucose) the amino acids present in the medium

ws.write(row, 0, 'D (h-1)', bstyle)
ws.write(row, 1, 'Biomass (gDW L-1)', bstyle)   #units change with the appropiate conversion
ws.write(row, 2, 'Flux (pg/cell day-1)', bstyle) #units change with the appropiate conversion
j=3
for i in N: 
 ws.write(row, j, i+ '(mM)', bstyle)
 j=j+1

wb.save(xlsfile)


# Parameters
Dmin, Dmax, DNsamps = 0.020, 0.0314, 2  # example of minimal, maximal and number of dilution rates (h-1)
X0, Xf, XNsamps = 1.0, 3.0, 5 # parameters to the cell densities to which the cell density will be constrained
#Xint = (Xf - X0)/(XNsamps - 1)
                                         # See function comProductivity() for details

# Load model
fnet = loadCHOmodel(filename = "iCHOv1", name = "CHOFN", solver = "cplex") #loading an example of genome-scale metabolic model

# Optimize the net for each dilution rate and each biomass
for D in np.linspace(Dmin, Dmax, DNsamps):
    for X in np.linspace(X0, Xf, XNsamps):
       row = row + 1
       flux, T = comProductivity(fnet, D, Glc=35.156197964437965, X0 = X, Xf = X, XNsamps = 1, Xint = 1) #concatenate before GLC (glucose) the amino acids present in the medium
       ws.write(row, 0, D)
       ws.write(row, 1, X)
       ws.write(row, 2, flux*0.001*24*150000*300/X) #converting from mM h-1 to pg/cell day-1
       j=3
       for i in T:   
        ws.write(row, j, i*150)
        j=j+1
       
wb.save(xlsfile)  

print(time.time()-start)
