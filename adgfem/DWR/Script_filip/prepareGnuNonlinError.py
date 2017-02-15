import numpy as np 
import sys 
from decimal import Decimal


def closest(num, li):
    higher = []

    for number in li:
        if number > num:
            higher.append(number)

    if higher:
        lowest = sorted(higher)[0]

        return lowest
    else:
        lowest = li[-1]
        return lowest
        
        

def setOutputFile( inputName, outName, tableName ):

   def addTexTable(tableName, adaptLevels, nelem, iterations, Etot, etaAll, etaS, etaN, etaA, etaAD, Itot):
      ''' 
         f - the TeX file to write into 
      ''' 

      table_beg = r'''
{\scriptsize             %  mensi
%{\tiny                   %  nejmensi
\begin{tabular}{cc|rrrr|r}
\hline
$N_{h}$ & iter & J(e) & 
 $\eta$ & $\eta_S $ & 
 $\eta_A $ & $i_{eff}$ \\ 
\hline  
         
         ''' 
      table_end = r''' 
\hline
\end{tabular}
} %%% end of font size  
\end{document}
         '''
      with open(tableName, 'write') as f: 
         f.write( table_beg ) 
         for i in adaptLevels: 
             f.write( '{:5}{:2}{:5}{:2}{:.2E}{:2}{:.2E}{:2}{:.2E}{:2}{:.2E}{:2}{:.2f}{:6}'.format( \
               int(nelem[i]), ' &',int(iterations[i]), ' &',Etot[i],' &', etaAll[i],' &', etaS[i],' &', etaN[i] + etaA[i] + etaAD[i],' &', \
               Itot[i], ' \\\\ \n' ) )
#             f.write('\\\\ a \n')
         f.write( table_end )                 
         f.close()
      return  
      
      

   x = []

   with open(inputName, 'r') as f: 
      for line in f:
           float_list = [float(i) for i in line.split()]
           x.append( float_list )
   #        print line.split(" ")
      f.close() 
      
   #print(x)

   last_row = len( x) -1

   Ju = x[last_row][5]  

   # compute Juh - always last row for each adapt level
   adaptLevels = []
   for i in range(1,len(x)): 
      if (x[i][0] > x[i-1][0]):
         adaptLevels.append(i-1)  
   adaptLevels.append(last_row)
         
   # set J(u_h)       
   Juh = []
   for i in range(len(adaptLevels)): 
      Juh.append( x[adaptLevels[i]][6] ) 

         
   # Eh = J(u) - J(uh)
   Eh = []
   for i in range(len(x)): 
      p = closest( i, adaptLevels )
   #   print p 
      Eh.append( Ju - x[p][6] )  
   #print( Eh) 


   # En = J(uh) - J(uhn) 
   En = []
   for i in range(len(x)): 
      p = closest( i, adaptLevels )
      En.append( x[p][6] - x[i][6])  
   #print( En) 

   # Etot = J(u) - J(uhn) 
   Etot = []
   for i in range(len(x)): 
      Etot.append( Ju - x[i][6]) 
      
   # nelem 
   nelem = [ x[i][1] for i in range(len(x)) ]  
   
   # iter 
   iterations = [ (x[i][3]+x[i][4]) for i in range(len(x)) ]    
   
   # eta_S 
   etaS = [ x[i][7] for i in range(len(x)) ] 
   #print 'etaS' , etaS 

   # eta_S 
   etaN = [ x[i][8] for i in range(len(x)) ] 
    

   # eta_S 
   etaA = [ x[i][9] for i in range(len(x)) ] 


   # eta_S 
   etaAD = [ x[i][10] for i in range(len(x)) ] 
   #print 'etaAD' , etaAD 
   
   #eta_all 
   etaAll = [ (x[i][7]+x[i][8]+x[i][9]+x[i][10]) for i in range(len(x)) ] 
     

   # Is = etaS / Eh  
   Is = [] 
   for i in range(len(etaS)): 
      Is.append( min(etaS[i] / Eh[i],100) )



   # In = etaN + etaA +etaAD / En 
   In = [] 
   for i in range(len(etaS)): 
      In.append( min((etaN[i] + etaA[i] + etaAD[i]) / Eh[i], 100) )


   # Itot = etaS + etaN + etaA +etaAD / Eh + En 
   Itot = [] 
   for i in range(len(etaS)): 
      Itot.append( min( abs( (etaS[i] + etaN[i] + etaA[i] + etaAD[i]) / (Eh[i] + En[i] )),100.) ) 
   
      

#with open(outName, 'w') as f: 
#      for i in range(len(x)): 
#          f.write( '{:5}{:15.8f}{:2}{:15.8f}{:2}{:15.8f}{:2}{:15.8f}{:2}{:15.8f}{:2}{:15.8f}{:2}{:15.8f}'.format( \
#            i+1, Eh[i],'  ', En[i],'  ', etaS[i],'  ', etaN[i] + etaA[i] + etaAD[i],'  ', \
#            Is[i],'  ', In[i],'  ', Itot[i] ) )
#          f.write('\n')
#          
#      f.close()
#   return 

   
   # nelem iter Etot eta, etaS, etaA, Itot
#   with open(tableName, 'w') as f: 
#      for i in adaptLevels: 
#          f.write( '{:5}{:5}{:15.8f}{:2}{:15.8f}{:2}{:15.8f}{:2}{:15.8f}{:2}{:15.8f}'.format( \
#            int(nelem[i]), int(iterations[i]), Etot[i],'  ', etaAll[i],'  ', etaS[i],'  ', etaN[i] + etaA[i] + etaAD[i],'  ', \
#            Itot[i] ) )
#          f.write('\n')
#          
#      f.close()
#   return      
      
      
   addTexTable(tableName, adaptLevels, nelem, iterations, Etot, etaAll, etaS, etaN, etaA, etaAD, Itot) 
   
   return 
      
   
   

inputName = '../SNA_P2_adapt_rezL2/aDWR_nlErrors' 
# not used
outName = 'effectivity' 
tableName = 'tableTexREZ_P2_NEW'

setOutputFile( inputName, outName, tableName)  


