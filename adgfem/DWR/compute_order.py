import numpy as np 
import sys


if len(sys.argv) != 5:
   sys.stderr.write('compute_order needs exactly 4 arguments! ')
   sys.stderr.write('dof1 dof2 error1 p')
   sys.stderr.write('\n')
   sys.exit(1) 
   


dof1 = float(sys.argv[1])
dof2 = float(sys.argv[2])
e1 = float(sys.argv[3])
p = int(sys.argv[4])
fileName = 'order-' + str(p) + '.dat'

# e2 = e1 * (N2/N1)^p
e2 = e1 * np.power( dof1/dof2 ,p)
print('e2:', e2) 

with open(fileName, 'w') as f:
   f.write(str(dof1) + ' ' + str(e1) +'\n') 
   f.write(str(dof2) + ' ' + str(e2) +'\n') 
   f.close() 
   
print('Data were written to file', fileName)
