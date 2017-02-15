from dolfin import * 
import numpy as np 

# number of elements
N = 4
p = 1

mesh = IntervalMesh( N, 0.0, 1.0 )
plot(mesh, interactive=True)


V = functionSpace( mesh, 'CG', p )

f = Expression('1.0')

def Boundary(x,on_boundary):
   return on_boundary 

bc = DirichletBC( V, 0.0, boundary ) 

#functions 
u = TrialFunction(V)
v = TestFunction(V) 

a = inner( grad(u) , grad(v) ) *dx  
l = f*v*dx 

# Solve problem 
u = Function(V) 
problem = LinearVariationalProblem( a, l, u, bc)
solver = LinearVariationalSolver( problem ) 
solver.solve() 

plot(u, title='Solution') 

#h = CellSize(mesh) 
#n = FacetNormal(mesh) 
      
