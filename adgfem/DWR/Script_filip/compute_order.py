import numpy as np 
import sys

# draws the line denoting given order in log scale
# can be used in comparison with real results

if len(sys.argv) != 4:
   sys.stderr.write('compute_order needs exactly 3  arguments! ')
   sys.stderr.write('order.dat k p')
   sys.stderr.write('k - k_th column of order.dat will be used (together with 4th - DOF)')
   sys.stderr.write('p = degree')
   sys.stderr.write('\n')
   sys.exit(1) 
   


#dof1 = float(sys.argv[1])
#dof2 = float(sys.argv[2])
k = int(sys.argv[2])
p = int(sys.argv[3])
fileName = 'order-' + str(p) + '.dat'
inputName = sys.argv[1]

# e2 = e1 * (N2/N1)^p
#e2 = e1 * np.power( dof1/dof2 ,p)
#print('e2:', e2) 
x = []

with open(inputName, 'r') as f:

   for line in f:
        float_list = [float(i) for i in line.split()]
        x.append( float_list )
#        print line.split(" ")
   f.close()
   
#print(x)   

last_row = len( x) - 1
#print( last_row ) 

dof1 = x[0][3] 
dof2 = x[last_row][3]
e1 = x[0][k-1]

#e2 = e1 * (N2/N1)^p
e2 = e1 * np.power( dof1/dof2 ,p)

if (len( x[last_row] ) < k) | (len( x[last_row] ) < 4): 
   sys.stderr.write('rows of the file are shorter than given k!\n')
   sys.exit(1)
else:
   with open(fileName, 'w') as f:
      f.write( str(dof1) + ' ' + str(e1) +'\n') 
      f.write( str(dof2) + ' ' + str(e2) +'\n') 
      f.close()

print('Data were written to file', fileName)
