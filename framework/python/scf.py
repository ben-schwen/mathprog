# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 10:41:29 2019

@author: ben
"""

import networkx as nx
import sys
sys.path =['','/home1/e01225371/kmst', '/home1/e01225371/anaconda3/lib/python35.zip', '/home1/e01225371/anaconda3/lib/python3.5', '/home1/e01225371/anaconda3/lib/python3.5/plat-linux', '/home1/e01225371/anaconda3/lib/python3.5/lib-dynload', '/home1/e01225371/.local/lib/python3.5/site-packages', '/home1/e01225371/anaconda3/lib/python3.5/site-packages', '/home1/e01225371/anaconda3/lib/python3.5/site-packages/Sphinx-1.4.6-py3.5.egg', '/home1/e01225371/anaconda3/lib/python3.5/site-packages/setuptools-27.2.0-py3.5.egg']

import time
import math
from gurobipy import *

def parse_graph(edges_str):
  n_nodes = int(edges_str.pop(0))
  edges_str.pop(0)
  
  G = nx.empty_graph(n_nodes)
  
  for e_str in edges_str:
    if(e_str == ''):
      continue
    e = e_str.split(' ')
    G.add_edge(int(e[1]), int(e[2]), weight=int(e[3]))
  
  return(G)

def scf(G,k):
  weights = [data['weight'] for (_,_,data) in G.edges(data=True)]
  edges = list(G.edges())
  di_edges = [(i,j) for (i,j) in G.edges()] + [(j,i) for (i,j) in G.edges()]

  try:

    model = Model("k_mst")
    model.setParam( 'OutputFlag', False )
    n = G.number_of_nodes()
    
    # CREATE VARIABLES
    y = model.addVars(n, vtype=GRB.BINARY, name='y')
    d = model.addVars(G.number_of_edges()*2, vtype=GRB.BINARY, name='d')
    f = model.addVars(G.number_of_edges()*2, vtype=GRB.INTEGER, name='f')
      
    # ADD CONSTRAINTS
    # Constraint 2
    model.addConstr(d.sum() == (k+1), "c_0")
    # Constraint 3  
    model.addConstr(y.sum() == (k+1), "c_1")
  
    # Constraint 4
    for i in range(1,G.number_of_nodes()):
      model.addConstr(d.sum([di_edges.index(x) for x in [(a,_) for (a,_) in di_edges if a==i]]) == y[i])
      #Constraint 7
      model.addConstr(
        f.sum([di_edges.index(x) for x in [(a,_) for (a,_) in di_edges if a==i]]) -
        f.sum([di_edges.index(x) for x in [(_,a) for (_,a) in di_edges if a==i]]) == y[i])
    
    # Constraint 5  
    #ingoing to 0
    model.addConstr(d.sum([di_edges.index(x) for x in [(_,a) for (_,a) in di_edges if a==0]])==1)
    #outgoing from 0
    model.addConstr(d.sum([di_edges.index(x) for x in [(a,_) for (a,_) in di_edges if a==0]])==1)  

    #Constraint 6 
    for v in range(G.number_of_edges()*2):
      (i,j) = di_edges[v]
      model.addConstr(d[v] <= y[i])
      model.addConstr(d[v] <= y[j])
      # Constraint 8
      model.addConstr(f[v] <= (n-1)*d[v])
      model.addConstr(f[v] >= 0)    
    
    # ADD OBJECTIVE
    expr = LinExpr()
    for i in range(G.number_of_edges()*2):
      if (weights[i % G.number_of_edges()] != 0):
        expr += weights[i % G.number_of_edges()]*d[i]
    model.setObjective(expr, GRB.MINIMIZE)
  
    model.write("scf.lp")
  
    start_time = time.time()
    model.optimize()
    end_time = time.time()
    
#==============================================================================
#     nodeSol = model.getAttr('X', y)
#     diSol = model.getAttr('X', d)
#   
#     y_ls = []
#     d_ls = []
#   
#   
#     for i in range(n):
#       if(nodeSol[i] > 0.5):
#         y_ls.append(i)
#       
#     for i in range(G.number_of_edges()*2):
#       if(diSol[i] > 0.5):
#         d_ls.append(di_edges[i])
#==============================================================================
      
    #print('Obj: %g' % model.objVal)
    #print(y_ls)  
    #print(d_ls)
  
    #G1 = nx.empty_graph(G.number_of_nodes())
    #G1.add_edges_from(d_ls)
    #plt.figure(2)
    #nx.draw_circular(G1, with_labels=True)
    #plt.show()
    return(model, end_time-start_time)
#==============================================================================
#     return("Objective = " + str(model.objVal) + 
#            "\tNodes explored = " + str(model.getAttr("NodeCount")) +
#            "\tTime = " + str(end_time-start_time))
#==============================================================================
    
  except GurobiError:
    print('Error reported')
    return("Error_X")


print(sys.argv)
if len(sys.argv) == 2:
  input_arg = sys.argv[1]
  fname = 'g0' + str(input_arg) + '.dat'
  filepath = '/home1/e01225371/kmst/data/' + fname
  
edges_str = []
#filepath = "/home/ben/Schreibtisch/mathprog/framework/data/g04.dat"
with open(filepath) as f:
        for l in f:
            edges_str += (l.strip().split('\n'))
            
G = parse_graph(edges_str)
k1 = math.ceil((G.number_of_nodes()-1)/5) 
k2 = math.ceil((G.number_of_nodes()-1)/2)

[model1, t1] = scf(G,k1)
[model2, t2] = scf(G,k2)

o0 = str("Model Inst \t k \t Obj \t Nodes \t Time")
o1 = ("SCF \t" + str(input_arg) +"\t"+ str(k1) +"\t"+ str(model1.objVal) 
      +"\t"+ str(model1.getAttr("NodeCount")) +"\t"+  str(t1))
o2 = ("SCF \t" + str(input_arg) +"\t"+ str(k2) +"\t"+ str(model2.objVal) 
      +"\t"+ str(model2.getAttr("NodeCount")) +"\t"+  str(t2))
#print(o0)
print(o1)
print(o2)

    