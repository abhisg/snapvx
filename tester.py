import time
from snapvx import *
import numpy as np

np.random.seed(1)
num_nodes = 25000
num_edges = 2000000
n = 50
snapGraph = GenRndGnm(PUNGraph, num_nodes, num_edges)
gvx = TGraphVX(snapGraph,use_proximal_updates=True)

#For each node, add an objective (using random data)
for i in range(num_nodes):
	x = Variable(n,name='x') #Each node has its own variable named 'x'
  	a = numpy.random.randn(n)
  	#gvx.SetNodeObjective(i, square(norm(x-a)))
	gvx.SetNodeProximalArgs(i,[x],proximalArgs={x:a})
#def netLasso(src, dst, data):
#	return (norm(src['x'] - dst['x'],2), [])
#gvx.AddEdgeObjectives(netLasso)
for edge in gvx.Edges():
	etup = edge.GetSrcNId(), edge.GetDstNId()
	gvx.SetEdgeProximalArgs(etup[0],etup[1])
	

"""#Create new graph
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
	start = time.time()
	gvx.Solve(MaxIters=1)
	print time.time() - start
	print gvx.PrintSolution()

