from dolfin import *
import mshr
import numpy as np

### mozna by bylo dobre zaokrouhlit souradnice ?

gridName ='square01-c14.grid'
subGridName = 'c14.subgrid'
# GEOMETRY
p1= Point(0.0, 0.0)
p2= Point(1.0, 1.0)


#p5 = Point(0.75, 0.75) 
#p6 = Point( 0.76, 0.76) 
p3 = Point( 0.45 , 0.15) 
p4 = Point( 0.55, 0.25)

#triangle_vertices = [p3,
#                     p4,
#                     p5]
#p4 = Point(0.755, 0.755)
##p3 = Point(0.49, 0.11)
##p4 = Point(0.51, 0.13)

#p5 = Point(0.245, 0.245)
#p6 = Point(0.255, 0.255) 



# #elems ~
N = 8

def distance(x,y): 
   vx = x[0] - y[0] 
   vy = x[1] - y[1] 
   rng = sqrt( vx*vx + vy*vy ) 
   return rng

#BOUNDARY CONDITIONS 
class Boundary(SubDomain): 
   def inside(self, x, on_boundary): 
      return on_boundary 

boundary = Boundary() 




x  = []
#triangle 
#x.append( [p3[0], p3[1]] )
#x.append( [p4[0], p4[1]] )
#x.append( [p5[0], p5[1]] )
#x.append( [p3[0], p3[1]] )

#rectangle
x.append( [p3[0], p3[1]] )
x.append( [p3[0], p4[1]] )
x.append( [p4[0], p4[1]] )
x.append( [p4[0], p3[1]] )
x.append( [p3[0], p3[1]] ) 



# SQUAREtriangle
#domain = mshr.Rectangle(p1,p2)
#domain.set_subdomain(1, mshr.Polygon(triangle_vertices) )
#domain.set_subdomain(2, mshr.Rectangle(p1,p2) - mshr.Polygon(triangle_vertices) ) 

# Square with two subdomains
#domain = mshr.Rectangle(p1,p2)

#domain.set_subdomain(1, mshr.Rectangle(p3, p4))
#domain.set_subdomain(2, mshr.Rectangle(p5, p6))
#domain.set_subdomain(3, mshr.Rectangle(p1,p2) - mshr.Rectangle(p3, p4) - mshr.Rectangle(p5, p6) ) 

# Square subdomain
domain = mshr.Rectangle(p1,p2)

domain.set_subdomain(1, mshr.Rectangle(p3, p4))
domain.set_subdomain(2, mshr.Rectangle(p1,p2) - mshr.Rectangle(p3, p4) ) 

#print 'Before mesh generation!'
mesh2d = mshr.generate_mesh(domain, N)

#mesh2d = UnitSquareMesh( N, N, 'crossed')

# coordinates
coords = mesh2d.coordinates()
n = FacetNormal(mesh2d)

plot(mesh2d, interactive = True) 


# Initialize mesh function for boundary domains
boundMark = FacetFunction("bool", mesh2d)
boundMark.set_all(False)
boundary.mark(boundMark, True)



# COORDINATES OF THE VERTICES 
coordinates = [] 
for v in vertices(mesh2d): 
   coordinates.append( coords[v.index()] )
#   print('coord: ',v.index(), coords[v.index()] )


# COORDINATES OF THE VERTICES 
vert = []
for cell in cells(mesh2d):
#   print "cell", cell.index(), "has edges :", cell.entities(1) 
#   print "cell", cell.index(), "has vertices :", cell.entities(0) 
   

   vertex = []
   firstTime = True 
   normala = np.array([ cell.normal(0).x(), cell.normal(0).y() ])

   for e in facets(cell): 
      if (firstTime):
         firstTime = False
#         print 'vertices?',  e.entities(0)[0] , e.entities(0)[1]  
         
#         normala = np.array([e.normal().x() , e.normal().y()])
         vec1 = coords[ e.entities(0)[1]  ] - coords[ e.entities(0)[0] ] 
#         print 'vec1:', vec1
         vec2 = np.array([-vec1[1],vec1[0]])  
         vec2 = vec2 / np.linalg.norm(vec2)
         if ( np.linalg.norm(vec2 + normala) > 1E-9 ): 
#            print 'if - stejne normaly - otocit'
            
            vertex = [ cell.entities(0)[1], cell.entities(0)[0], cell.entities(0)[2]  ]
#            print 'vertex', vertex
#            print 'vec2', vec2 
#            print 'normala', normala
#            print 'norma:' , np.linalg.norm(vec2 + normala) 

         else: 
#            print 'else'
            vertex = [ cell.entities(0)[0], cell.entities(0)[1], cell.entities(0)[2]  ]
#            print 'vertex', vertex
   vert.append(vertex)




edges = []
# BOUNDARY EDGES indices of the vertices
for e in facets(mesh2d): 
   edge = []
   if ( boundMark[ e.index() ] ):
      for v in vertices(e):
         edge.append( v.index() ) 
         #print('v', v.index() )
#      print 'edge:', edge , coords[ edge[0] ], coords[edge[1]]
#      print 'normala:' ,  e.normal().x(), e.normal().y()
      
      normala = np.array([e.normal().x() , e.normal().y()])
      vec1 = coords[ edge[1] ] - coords[ edge[0] ] 
#      print 'vec1' , vec1
#       normala - otoceni smeroveho vektoru strany o 90 stupnu proti smeru hod rucicek
      vec2 = np.array([-vec1[1],vec1[0]])  
      vec2 = vec2 / np.linalg.norm(vec2)
#       vec2 ma tedy opacny smer nez vnejsi normala
#      print 'normala:', normala
#      print 'vec:',vec2
      # musime otocit poradi vrcholu na strane
      if ( np.linalg.norm(vec2 + normala) > 0 ):
#         print 'edge:' , edge
#         print 'normala:', normala
#         print 'vec:',vec2
         edge.reverse()
         if ( np.linalg.norm(vec2 + normala) < 0.5 ):
            print 'value of the sum:', np.linalg.norm(vec2 + normala)
      #   print 'edge again:' , edge
       
      edges.append(edge)   


# BASIC parameters    
npoin = mesh2d.num_vertices()    
nelem = mesh2d.num_cells()
nbelm = len(edges)
#print 'nbelm:', nbelm
nbc = 1
ibc = nbc
triangle = 3 

# WRITE INTO THE FILE

# File open 
with open(gridName, 'w') as f:
   f.write(str(npoin)+ ' ') 
   f.write(str(nelem) + ' ')
   f.write(str(nbelm) + ' ')
   f.write(str(nbc) + '\n')
   #periodicity  
   f.write('0.  0. 0 0  0.  0. 0 0 \n' )
   #vertices 
   for c in coordinates: 
      for cc in c: 
         f.write( str(cc) + ' ' )
      f.write('\n')
   #elems
   for v in vert:
      f.write( str(triangle) + ' ' )
      for w in v:
         f.write( str( w + 1 ) )
         f.write(' ')
      f.write('\n')
   for e in edges: 
      for ie in e: 
         f.write( str(ie + 1) + ' ')
      f.write( str(ibc) )
      f.write('\n')
   f.close()
   

# SUBGRID

with open(subGridName, 'w') as f:
   f.write(str(4)+ '\n') 
   f.write('convex \n')
   for i in range(len(x) -1): 
      f.write(str(i+1) + ' ')
      for j in range(len(x[i])): 
         f.write(str(x[i][j]) + ' ')
      for j in range(len(x[i+1])): 
         f.write(str(x[i+1][j]) + ' ')    
      f.write('\n')   
   f.close()
   
   
print 'DONE with ', nelem, 'elems!'    
