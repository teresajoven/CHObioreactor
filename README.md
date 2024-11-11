# CHObioreactor
<p><div class="text-justify">
This repository contains the Python script, <em>CHObioreactor.py</em>, that performs the computational model corresponding to the bioreactor of the figure presented in the document <em>bioreactor.pdf</em>.  
  
Synthesis reactions from a generic file _reactions.csv_(for example, for antibodies production) can be added to a model extracted from [BiGG Models], (http://bigg.ucsd.edu/models/), which is an open source knowledgebase of genome-scale metabolic models (GEMs). This modified model is converted into a Flexible Net (FN) by utilizing the open-source Python package _fnyzer_ (Flexible Nets analYZER). An FN allows for modeling the multiscale system with the integration of a metabolic network and the dynamics of the bioreactor. Then, adding flux bounds, imposing steady state conditions and introducing an objective function (as maximizing the antibodies production) the problem to be addressed is a linear programming problem (LPP). 

The files <em>heatmap_excel</em> and <em>heatmap_plot</em> enable the visualization of the production on a heatmap.

### Requirements
To run the script, it will be necessary to install:

*Python package *cobrapy* that provides a simple interface to metabolic constraint-based model (CBM) reconstruction and analysis.

*Python libraries: NumPy, to create multidimensional vectors and matrices, and Pandas, which will allow for data manipulation and storage. 

*Optimization software for fnyzer to solve the LPP associated with the FN. For instance, solvers such as GLPK, Gurobi or CPLEX are appropiate.


## Step-by-Step Guide: Modeling development process

**1. Adding the necessary synthesis reactions to the GEM.**

**2. Modeling the bioreactor dynamics.**

**2.1. Including uptake rates of medium metabolites.**

**2.2. Modeling the nutrient uptake from the tank into the cell.**

**2.3. Transporting the nutrient out of the tank to the effluent.**

**2.4. Biomass production.**

**2.5. Analyzing the productivity (as antibodies production).**
</div></p>
