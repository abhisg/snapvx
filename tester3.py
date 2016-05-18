import time
from snapvx import *
import numpy as np

np.random.seed(1)
num_nodes = 2
num_edges = 1
n = 50
snapGraph = GenRndGnm(PUNGraph, num_nodes, num_edges)
gvx = TGraphVX(snapGraph,use_proximal_updates=True)

mat = np.random.randn(n+1,100*(n+1))
mat = np.cov(mat)
#mat = (mat+mat.T)/2
print mat
rho = 1.0
lamb = 10.0
"""Theta = np.zeros((n+1,n+1))
Z = np.zeros((n+1,n+1))
U = np.zeros((n+1,n+1))
for i in xrange(10):
	mt = rho*(Z+Z.T-U-U.T)/2 - mat
	w,v = np.linalg.eigh(mt)
	w = w + np.sqrt(w**2 + 4*rho)
	Theta = np.dot(v,np.dot(np.diag(w),v.T))*1.0/(2*rho)
	for j in xrange(n+1):
		Z[0,j] = Theta[0,j] + U[0,j]
		Z[j,0] = Theta[0,j] + U[0,j]
	for j in xrange(1,n+1):
		for k in xrange(1,n+1):
			if j == k:
				Z[j,k] = Theta[j,k] + U[j,k]
			else:
				theta = max(0,1-lamb/(rho*(Theta[j,k]+U[j,k])))
				Z[j,k] = theta * (Theta[j,k] + U[j,k])
	U += Theta - Z
	print "Theta",Theta
	print "Z",Z
	print "U",U"""
variables = []
constants = {}
var = Variable(1,name='nu')
variables.append(var)
constants[var] = {'a':[mat[0,0]]}
for i in xrange(n):
	var = Variable(1,name='x'+str(i))
	variables.append(var)
	constants[var] = {'a':[mat[i+1,0]]}

for i in xrange(n):
	for j in xrange(n):
		var = Variable(1,name='x'+str(3+i*3+j))
		variables.append(var)
		constants[var] = {'a':[mat[i+1,j+1]]}
gvx.SetNodeProximalArgs(0,variables,proximalArgs = constants,proximalOperator="NETLAPLACE")
x = Variable(n,name='x')
gvx.SetNodeProximalArgs(1,[x],proximalArgs={x:{'a':'foobar'}},proximalOperator="FOOBAR")
gvx.SetEdgeProximalArgs(0,1,proximalOperator="EDGELASSO")
	
#For each node, add an objective (using random data)
"""for i in range(num_nodes):
	x = Variable(n,name='x') #Each node has its own variable named 'x'
  	b = np.random.randn(n)
	mu = np.random.rand()*0.01
	y = np.random.randint(0,2)
  	#gvx.SetNodeObjective(i, square(norm(x-a)))
	#gvx.SetNodeProximalArgs(i,[x],proximalArgs={x:{'a':b}},proximalOperator=SQUARE)
	gvx.SetNodeProximalArgs(i,[x],proximalArgs={x:{'b':b,'mu':mu,'y':y}},proximalOperator="MOD_SQUARE")
#def netLasso(src, dst, data):
#	return (norm(src['x'] - dst['x'],2), [])
#gvx.AddEdgeObjectives(netLasso)
for edge in gvx.Edges():
	etup = edge.GetSrcNId(), edge.GetDstNId()
	gvx.SetEdgeProximalArgs(etup[0],etup[1],proximalOperator="NETLASSO")
	

#Create new graph
gvx = TGraphVX()

#Use CVXPY syntax to define a problem
x1 = Variable(2, name='x1')
obj = square(x1)
#Add Node 1 with the given objective, with the constraint that x1 <= 10
gvx.AddNode(1, Objective=obj, Constraints=[],proximalargs=[np.array([1.5,2.5])])
#gvx.AddNode(1, obj, [])

#Similarly, add Node 2 with objective |x2 + 3|
x2 = Variable(2, name='x2')
#obj2 = abs(x2+3)
obj2 = square(x2)
gvx.AddNode(2, obj2, [],proximalargs=[np.array([-3,-4])])

#Similarly, add Node 3 with objective |x2 + 3|
x3 = Variable(2, name='x3')
#obj2 = abs(x2+3)
obj2 = square(x3)
gvx.AddNode(3, obj2, [],proximalargs=[np.array([-3,-3])])

#Add an edge between the two nodes, 
#Define an objective, constraints using CVXPY syntax
gvx.AddEdge(1, 2, Objective=square(norm(x1 - x2)), Constraints=[])
gvx.AddEdge(1, 3, Objective=square(norm(x1 - x3)), Constraints=[])
gvx.AddEdge(2, 3, Objective=square(norm(x2 - x3)), Constraints=[])
"""
#Solve the problem, and print the solution
if __name__ == '__main__':
#	start = time.time()
	gvx.Solve(Rho=1.0,Lambda=lamb,EpsAbs=0.000001,EpsRel=0.000001)
#	print time.time() - start
	print gvx.PrintSolution()

