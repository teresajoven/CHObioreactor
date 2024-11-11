# -*- coding: utf-8 -*-

# Functions to compute the maximum productivity of a Flexible Net integrating a metabolic network with synthesis reactions and bioreactor variables.

from __future__ import division, print_function
import numpy as np
import cobra
from fnyzer import FNFactory, cobra2fn
import pandas as pd
import numpy as np 



def comProductivity(fnet,D,Glc,X0,Xf,XNsamps,Xint): #concatenate before Glc (glucose) the amino acids present in the medium
 
    
    for X in np.linspace(X0, Xf, XNsamps):
        genCHOBiorFN(fnet,D,Glc,Xmin=X0,Xmax=Xf)
        netobj = FNFactory(fnet) # Build net object
        T=[]
        try:
            netobj.optimize()
            print (f"{netobj.objval:.15f}")
            optX = netobj.places['X'].avm
            D=dict(glc=Glc)                       #concatenate before Glc (glucose) the amino acids present in the medium
            flux = netobj.trans['taout'].avl
            for i in D:     
              nut_tank =netobj.places[i.upper()].avm
              nut_out = netobj.trans['t'+i+'out'].avl
              nut_in = netobj.trans['t'+i+'in'].avl
              nut_incell = netobj.trans['t'+i+'t'].avl
              print(str(i.upper())+' in tank: '+str(nut_tank)+" mM")
              print('Flux of '+str(i)+' entering the cell: '+str(nut_incell)+" mM h-1")
              print('Flux of '+str(i)+' leaving through the effluent: '+str(nut_out)+" mM h-1")
              T.append(nut_tank)
            print("Flux of iggM1 with X in [",f"{X:.2f}",",",f"{Xf:.2f}","] gDW L-1:", f"{flux:.15f}", "mM h-1")       
        except:
            print("  Problem NOT feasible")
            flux=1000
    return flux, T
    
def loadCHOmodel(filename = "iCHOv1", name = "CHOFN", solver = "cplex"): #loading an example of genome-scale metabolic model
   
    model = cobra.io.read_sbml_model(filename+'.xml')
    exch = list(model.exchanges)

    for i in exch: 
     if i.lower_bound != 0.0:
          i.lower_bound = -1000
          
    # adding syntesis reactions
    antibody = wpd.read_csv("reactions.csv") #example of antibodies synthesis reactions
    ab_names = ["antiCD20", "iggM1", "iggM2", "iggM3", "mAb", "mAb2"]
    for ab in ab_names:
        reaction = cobra.Reaction(ab)
        reaction.name = "%s synthesis" % ab
        reaction.subsystem = "PROTEIN PRODUCTION"
        reaction.lower_bound = 0
        reaction.upper_bound = 0
        print(reaction.name)
        r_dict = {}
        for i in range(len(antibody)):
            met = model.metabolites.get_by_id(antibody.loc[i]["name"][2:]) 
            r_dict[met] = round(antibody.loc[i]["%s_coef" % ab], 3)         
        reaction.add_metabolites(r_dict)
        model.add_reactions([reaction])
    fnet = cobra2fn(model) # Build Flexible Net
    fnet['name'] =  name
    fnet['solver'] = solver
    
    return fnet


def genCHOBiorFN(fnet,D,Glc,Xmin,Xmax): #concatenate before Glc (glucose) the amino acids present in the medium
    
    ### Bioreactor(tank) concentrations
    fnet['places']['A'] = 0.0
    fnet['places']['X'] = 1.9512
    fnet['mbounds'] = [str(Xmin)+"<=m['X']", "m['X']<="+str(Xmax)]
   
    ### Bioreactor reactions (transitions and handlers)
    fnet['shandlers'] = {}
    
    N=dict(glc=Glc) #concatenate before Glc (glucose) the amino acids present in the medium
    
    for i in N:    
      fnet['places'][i.upper()] = 0.0
      # Nutrient feed
      fnet['trans']['t'+i+'in'] = {'l0': D*N[i], 'a0': 0}
      fnet['vhandlers']['v'+i+'in'] = [{'a':('v'+i+'in',i.upper()), 'v':('t'+i+'in','v'+i+'in')},
                                  'a == v']

      # Nutrient uptake (from the tank into the cell)
      fnet['trans']['t'+i+'t'] = {'l0': 0, 'a0': 0}
      fnet['vhandlers']['v'+i+'t'] = [{'a':(i.upper(),'v'+i+'t'), 'v':('t'+i+'t','v'+i+'t')},
                                 'a == v']
      if i=='glc':
       nut_ex = 'EX_glc__D_e' # Reaction ID of glucose exchange in metabolic network
      else:
       nut_ex = 'EX_'+i+'__L_e'   #ID of exchange reactions
        
      nut_in = 't_'+nut_ex+'_b'  # b is for backward reaction and f is for forward reaction
      fnet['trans'][nut_in]['l0'] = 0 # Nutrient uptake determined by its intensity handler
      nut_out = 't_'+nut_ex+'_f'
      fnet['trans'][nut_out]['l0'] = 0 # No nutrient out of the cell allowed
      fnet['shandlers']['h'+i] = [{i+'g':('h'+i,nut_in), i+'t':('h'+i,'t'+i+'t')}, i+'g*'+str(Xmin)+'<= '+i+'t', i+'t <= '+i+'g*'+str(Xmax)] 
      
      # Nutrient out of the tank (to effluent)
      fnet['trans']['t'+i+'out'] = {'l0': 0, 'a0': 0}
      fnet['vhandlers']['v'+i+'out'] = [{'a':(i.upper(),'v'+i+'out'), 'v':('t'+i+'out','v'+i+'out')},
                                  'a == v']
      fnet['shandlers']['s'+i+'out'] = [{'g':(i.upper(),'s'+i+'out'), 'r':('s'+i+'out','t'+i+'out')},
                                   'r == g*'+str(D)]

    # Cell growth
    fnet['trans']['txt'] = {'l0': 0, 'a0': 0}  
    fnet['vhandlers']['vxt'] = [{'a':('vxt','X'), 'v':('txt','vxt')},
                                'a == v']
    biomass_reaction = 'BIOMASS_cho' # Reaction ID of biomass in metabolic network
    tBiomass = 't_'+biomass_reaction+ '_f' # transition for biomass
    fnet['trans'][tBiomass]['l0'] = 0 # Biomass production (i.e. growth rate) determined by its intensity handler. It should be equal to D in steady state.
    fnet['shandlers']['hr'] = [{'r':('hr',tBiomass), 'rt':('hr','txt')},
                                'r*'+str(Xmin)+'<= rt', 'rt <= r*'+str(Xmax)]

    # Cells out of the tank (to effluent)
    fnet['trans']['txout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vxout'] = [{'a':('X','vxout'), 'v':('txout','vxout')},
                                  'a == v']
    fnet['shandlers']['sxout'] = [{'x':('X','sxout'), 'r':('sxout','txout')},
                                   'r == x*'+str(D)]
    
    # Antibodies production (from cell to tank)
    fnet['trans']['tat'] = {'l0': 0, 'a0': 0} 
    fnet['vhandlers']['vat'] = [{'a':('vat','A'), 'v':('tat','vat')},
                               'a == v']
    tab = 't_iggM1_f'
    fnet['trans'][tab]['l0'] = 0 
    fnet['shandlers']['ha'] = [{'z':('ha',tab), 'at':('ha','tat')},
                                'z*'+str(Xmin)+'<= at', 'at <= z*'+str(Xmax)]

    # Antibodies out of the tank (to effluent)
    fnet['trans']['taout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vaout'] = [{'a':('A','vaout'), 'v':('taout','vaout')},
                                  'a == v']
    fnet['shandlers']['saout'] = [{'c':('A','saout'), 'r':('saout','taout')},
                                   'r == c*'+str(D)]    
      

    # Other net options
    L=[]
    for i in N: 
     L.append(i.upper())
    L.append('X')
    L.append('A')
    fnet['obj'] = {'f': "avl['tat']", 'sense': 'max'}
    fnet['exavtrans'] = 'all'
    fnet['extrans'] = 'all'
    fnet['actavplaces'] = L
    fnet['actplaces'] = L
    fnet['options'] = {
            'antype': 'cst', #mpc
            'savenet': False,            
            'printres': False,
            'printmodel': False,
            'writevars': {
                'avm': L,
                'avl':'all'},
            'plotres': False,
            'writexls': True,
            }

fnet = loadCHOmodel(filename = "../BIGG/iCHOv1", name = "CHOFN", solver = "cplex") #take a genome-scale metabolic model of BIGG models
comProductivity(fnet,D=0.0314,Glc=35.15619796,X0=1.0,Xf=1.0,XNsamps=1,Xint=1): #concatenate before Glc (glucose) the amino acids present in the medium

