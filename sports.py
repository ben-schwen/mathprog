from gurobipy import *

n = 20
k = 3
teams = list(range(0,n))

m = Model("sports")

xij = [(i,j,k) for i in range(n) for j in range(n) for k in range(3) if i != j]

plays = m.addVars(xij, vtype=GRB.BINARY)
points = m.addVars(teams, vtype=GRB.INTEGER)

m.update()
m.modelSense = GRB.MAXIMIZE

for a in range(0,n):
  for b in range(0,n):
    if(a != b):
      s = "Constraints for play" + str(a) + "|" + str(b)
      m.addConstr(sum([plays[x] for x in [(i,j,k) for (i,j,k) in xij if (i==a and j==b)]]), GRB.EQUAL, 1, s)

for t in range(0,(n-1)):
  m.addConstr(points[t+1], GRB.LESS_EQUAL, points[t], "order")
  
for t in range(0,n):
  s = "Constraints for points" + str(t)
  m.addConstr((sum([plays[x] for x in [(i,j,k) for (i,j,k) in xij if (i==t or j==t) and k==0]]) +
    3*sum([plays[x] for x in [(i,j,k) for (i,j,k) in xij if (i==t) and k==1]]) + 
    3*sum([plays[x] for x in [(i,j,k) for (i,j,k) in xij if (j==t) and k==2]])), GRB.EQUAL, points[t], s)
m.setObjective(points[n-k-1], GRB.MAXIMIZE)
m.write("sports.lp")

m.optimize()

