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
import matplotlib.pyplot as plt
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
  
def dcut(model, where):
  if where == GRB.Callback.MIPSOL:
    n = G.number_of_nodes()
    edges = list(G.edges())
    di_edges = [(i,j) for (i,j) in G.edges()] + [(j,i) for (i,j) in G.edges()]
    
    nodesSol = model.cbGetSolution(model._y)
    archesSol = model.cbGetSolution(model._d)
    edgesSol = model.cbGetSolution(model._e)
  
    n_ls = []
    a_ls = []
    e_ls = []
  
    for i in range(n):
      if(nodesSol[i] > 0.5):
        n_ls.append(i)

    for i in range(G.number_of_edges()*2):
      if(archesSol[i] > 0.5):
        a_ls.append(di_edges[i])
        
    for i in range(G.number_of_edges()):
      if(edgesSol[i] > 0.5):
        e_ls.append(di_edges[i])
        
    G1 = nx.empty_graph()
    G1.add_nodes_from(n_ls)
    G1.add_edges_from(a_ls)
    G1.remove_node(0)
    
    n_ls.sort()
    a_ls.sort()
    
    components = list(nx.connected_components(G1))
    if(len(components) > 1):
      for S in components:
        good_edges = [(u,v) for (u,v) in edges if ((u in S and v not in S) or (v in S and u not in S)) and (u*v != 0)]
        good_index = [edges.index(x) for x in good_edges]

        if(len(good_index) > 0):
          for i in S:
            model.cbLazy(model._e.sum(good_index) >= model._y[i])
            break
        global counter
        counter += 1
      

def dcc(G,k):
  weights = [data['weight'] for (_,_,data) in G.edges(data=True)]
  edges = list(G.edges())
  di_edges = [(i,j) for (i,j) in G.edges()] + [(j,i) for (i,j) in G.edges()]

  try:

    model = Model("k_mst")
    model.setParam( 'OutputFlag', False )
    n = G.number_of_nodes()
    model.Params.lazyConstraints = 1
    global counter    
    counter = 0
    # CREATE VARIABLES
    y = model.addVars(n, vtype=GRB.BINARY, name='y')
    d = model.addVars(G.number_of_edges()*2, vtype=GRB.BINARY, name='d')
    e = model.addVars(G.number_of_edges(), vtype=GRB.BINARY, name='e')
    
    model.update()
    # ADD CONSTRAINTS

    for i in range(G.number_of_edges()):
      (u,v) = di_edges[i]      
      #if(u*v != 0):
      model.addConstr(d[i]+d[i+G.number_of_edges()] == e[i])
      model.addConstr(d[i]+d[i+G.number_of_edges()] <= 1)

    # Constraint 2
    model.addConstr(d.sum() == (k+1), "c_0")
    # Constraint 3  
    model.addConstr(y.sum() == (k+1), "c_1")
  
    # Constraint 4
    for i in range(G.number_of_nodes()):
      model.addConstr(d.sum([di_edges.index(x) for x in [(_,a) for (_,a) in di_edges if a==i]]) <= y[i])
  
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
      
      
    # ADD OBJECTIVE
    expr = LinExpr()
    for i in range(G.number_of_edges()*2):
      if (weights[i % G.number_of_edges()] != 0):
        expr += weights[i % G.number_of_edges()]*d[i]
    model.setObjective(expr, GRB.MINIMIZE)
  
    model.write("dcc.lp")
  
    model._y = y
    model._d = d
    model._e = e 
    start_time = time.time()
    model.optimize(dcut)
    end_time = time.time()
 
    return(model, end_time-start_time, counter)

  except GurobiError:
    print('Error reported')
    return("Error_X")

if len(sys.argv) == 3:
  input_arg = sys.argv[1]
  k = int(sys.argv[2])
  fname = 'g0' + str(input_arg) + '.dat'
  filepath = '/home1/e01225371/kmst/data/' + fname  
  k1=k 
#  
#input_arg = 3
#
#fname = 'g0' + str(input_arg) + '.dat'
#filepath = '/home/ben/Schreibtisch/mathprog/framework/data/' + fname
  
edges_str = []
with open(filepath) as f:
        for l in f:
            edges_str += (l.strip().split('\n'))
            
G = parse_graph(edges_str)
#k1 = math.ceil((G.number_of_nodes()-1)/5) 
#k2 = math.ceil((G.number_of_nodes()-1)/2)

[model1, t1, c1] = dcc(G,k1)
#[model2, t2, c2] = dcc(G,k2)
#
#o0 = str("Model Inst \t k \t Obj \t Nodes \t Time")
o1 = ("DCC \t" + str(input_arg) +"\t"+ str(k1) +"\t"+ str(int(model1.objVal)) 
      +"\t"+ str(model1.getAttr("NodeCount")) +"\t"+  str(t1) + "\t" + str(c1))
#o2 = ("DCC \t" + str(input_arg) +"\t"+ str(k2) +"\t"+ str(int(model2.objVal)) 
#      +"\t"+ str(model2.getAttr("NodeCount")) +"\t"+  str(t2) + "\t" + str(c2))
#print(o0)
print(o1)
#print(o2)