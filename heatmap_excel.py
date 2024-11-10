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

N=["GLC","GLN","PHE","ARG","ASN","ASP","LYS","MET","HIS","ILE","LEU","PRO","VAL","SER","THR","TRP","TYR"]

ws.write(row, 0, 'D (h-1)', bstyle)
ws.write(row, 1, 'Biomass (gDW L-1)', bstyle)
ws.write(row, 2, 'Flux (pg/cell day-1)', bstyle)
j=3
for i in N: 
 ws.write(row, j, i+ '(mM)', bstyle)
 j=j+1

wb.save(xlsfile)


# Parameters
Dmin, Dmax, DNsamps = 0.020, 0.0314, 2  # (h-1) Minimal, maximal and number of dilution rates 
X0, Xf, XNsamps = 1.0, 3.0, 5 # Parameters to the cell densities to which the cell density will be constrained

#Xint = (Xf - X0)/(XNsamps - 1)
                                         # See function comProductivity() for details

# Load model
fnet = loadCHOmodel(filename = "iCHOv1", name = "CHOFN", solver = "cplex")

# Optimize the net for each dilution rate and each biomass
for D in np.linspace(Dmin, Dmax, DNsamps):
    for X in np.linspace(X0, Xf, XNsamps):
       row = row + 1
       flux, T = comProductivity(fnet, D, Glc=35.15619796,Gln=7.872825297,Phe=1.275823895,Arg=1.946111508,Asn=5.122945878,Asp=1.354860924,Lys=2.380196919,Met=1.010750472,His=1.216191003,Ile=2.944437965,Leu=3.964891601,Pro=4.621036929,Val=2.946237191,Ser=4.956440488,Thr=2.774941185,Trp=0.926515767,Tyr=0.978250505, X0 = X, Xf = X, XNsamps = 1, Xint = 1)
       ws.write(row, 0, D)
       ws.write(row, 1, X)
       ws.write(row, 2, flux*0.001*24*150000*300/X) #convierte mM h-1 en pg/cell day-1
       j=3
       for i in T:   
        ws.write(row, j, i*150)
        j=j+1
       
wb.save(xlsfile)  

print(time.time()-start)
