#!/usr/bin/env python

import bitstring
import re
import math
import numpy
import sys
import time
import openbabel
import pybel
import os
import random
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from optparse import OptionParser
import pylab


__DEBUG__=False
#Changes SPE iterations from 10000 to 100 to test new changes
#to code



def distance(pair1,pair2):

  square_sum=0.0
  for i in range(len(pair1)):
    square_sum+=(pair1[i]-pair2[i])**2

  return math.sqrt(square_sum)

###############################################################################
def tanimoto(a,b):
  #Calculates the Tanimoto dissimilarity for two bitstrings

  diffs=(a^b).count(True)  
  #Counts the number of bits that are different in the two bitstrings

  ones=(a&b).count(True)
  #Counts the number of bits that are both 1 in the bitstring

  
  return float(diffs)/(diffs+ones+2.2250738585072014e-308)
  #Must add smallest machine number to denominator to avoid dividing by zero 
  #in case two bitstrings are exactly the same

###############################################################################

###############################################################################
def get_bitstring_bin(unformatted_bitstring):

  return bitstring.BitStream(bin=re.sub(r'[^\w]','',unformatted_bitstring))
  #Get rid of any crap that might be in bitstring, JChem's mdgenerate inserts
  #pipe character "|" every 8 bits    
  #BitStream object

###############################################################################

###############################################################################
def get_bitstring_int(unformatted_intstring):

  first_line=True
  for byte in unformatted_intstring.split('\t'):
    if first_line:
      cf=bitstring.BitArray(int=int(byte),length=32)
      first_line=False
    else:
      cf+=bitstring.BitArray(int=int(byte),length=32)
   
  return cf
  

###############################################################################


###############################################################################
def calculate_tanimoto_matrix(molfile):


  for i in range(len(molfile)):

    for j in range(len(molfile)-1,i-1,-1):
      print'{0:.2f}'.format(tanimoto(get_bitstring_int(molfile[i].data['CF']), \
                                     get_bitstring_int(molfile[j].data['CF'])))
      #tanimoto(get_bitstring_int(molfile[i].data['CF']),get_bitstring_int(molfile[j].data['CF']))


    

###############################################################################

###############################################################################

def read_sdf_file(filename):

  basename,extname=os.path.splitext(str(filename))
  #uses os.path to split the filename into its base and extenstion

  return pybel.readfile(re.sub(r'[^\w]','',extname),filename)
  #return list of mol objects which can be manipulated in the main program.

###############################################################################

###############################################################################
def SPE_Embed(mollist,dim,loop_ratio):

  dim=int(dim)
  global __DEBUG__

  cutoff_length=0.4

  mollist_size=len(mollist)
  inner_loop_iterations=loop_ratio*mollist_size
  outer_loop_iterations=10000

  if __DEBUG__: outer_loop_iterations=10
  
  initial_learning_rate=1.0
  
  learning_rate=initial_learning_rate
  
  if dim == 2:
    coords=numpy.random.random((mollist_size,2))
    #Initialize initial 2D values for each molecule
  elif dim == 3:
    coords=numpy.random.random((mollist_size,3))

  tanimoto_matrix=numpy.ones((mollist_size,mollist_size))



 
  for olindex in range(outer_loop_iterations):
  # Outer loop
    ilindex=0
    for ilindex in range(int(inner_loop_iterations)):
  
      mol1_index=random.randint(0,mollist_size-1)
      mol2_index=random.randint(0,mollist_size-1)
      #Select two numbers within domain of mollist


      while mol1_index==mol2_index:
        mol2_index=random.randint(0,mollist_size-1)
        # Make sure we pick two different molecules

      if mol2_index > mol1_index: 
        temp=mol2_index
        mol2_index=mol1_index
        mol1_index=temp
        #We want mol1_index to larger so we only have to look at lower triangle
        #of tanimoto_matrix

      mol1=mollist[mol2_index]
      mol2=mollist[mol1_index]
      # pick the two randomol_listm molecules

      if tanimoto_matrix[mol1_index][mol2_index] < 1:
      #Check if the tanimoto distance has been calculated

        tanimoto_distance=tanimoto_matrix[mol1_index][mol2_index]
        #Retrieve value as tanimoto_distance          
  
      else:
        tanimoto_distance=tanimoto(get_bitstring_int(mol1.data['CF']), \
                                   get_bitstring_int(mol2.data['CF']))
        #Elsewise, calculate it 

        tanimoto_matrix[mol1_index][mol2_index]=tanimoto_distance
        #And store it in matrix


      embedded_distance=distance(coords[mol1_index],coords[mol2_index])
      if (tanimoto_distance > cutoff_length and embedded_distance >= tanimoto_distance):
        pass
      else: 

     
        # We don't want to update any coords whose distance above the cutoff length
        # and whose distance is greater than the original distance

      

        scale_factor=0.5*learning_rate*(tanimoto_distance - embedded_distance)/ \
                                 (embedded_distance + 2.2250738585072014e-308)

        coords[mol1_index]+=scale_factor*(coords[mol1_index]-coords[mol2_index])
        coords[mol2_index]+=scale_factor*(coords[mol2_index]-coords[mol1_index])                                 
        #Update coordinates 
        
              
    
    learning_rate-= initial_learning_rate/outer_loop_iterations

    
    if olindex%(outer_loop_iterations/10)==0: print float(olindex)/(outer_loop_iterations)*100,"% "

  for molindex in range(len(mollist)):
    mollist[molindex].data['SPE.x']=coords[molindex][0]
    mollist[molindex].data['SPE.y']=coords[molindex][1]
    
    if dim == 3:
      mollist[molindex].data['SPE.z']=coords[molindex][2]

  return mollist

###############################################################################

##############################################################################
def plot_spe(mol_list,dim):

  dim=int(dim)

  cluster_coords=[]

  #color_array=['b','g','r','c','m','y','k']

  number_clusters=0
  for molindex in range(len(mol_list)):
    if int(mol_list[molindex].data['clid']) > number_clusters:
      number_clusters=int(mol_list[molindex].data['clid'])
  # Need to find out how many clusters there are      
  

  color_array=[]
  cm=pylab.get_cmap('Paired')  
  for i in range(number_clusters):
    color_array.append(cm(1.*i/number_clusters))


  for i in range(number_clusters):
    cluster_coords.append([])
  #Give each cluster a list member, to plot for later  

  cluster_coords.append([])
  #For singletons


  for molindex in range(len(mol_list)):

    mol=mol_list[molindex]

    if int(mol.data['clid']) > 0:
    # If compound is a cluster member, put it into its respective container
      if dim==2:
        cluster_coords[int(mol.data['clid'])-1].append([mol.data['SPE.x'],mol.data['SPE.y']])
      elif dim==3:
        cluster_coords[int(mol.data['clid'])-1].append([mol.data['SPE.x'],mol.data['SPE.y'],mol.data['SPE.z']])  
   
    else: 
    #Else, put into the list of singletons at the end of cluster_coords container
      if dim==2:
        cluster_coords[number_clusters].append([mol.data['SPE.x'],mol.data['SPE.y']])
      elif dim==3:
        cluster_coords[number_clusters].append([mol.data['SPE.x'],mol.data['SPE.y'],mol.data['SPE.z']])

  fig = plt.figure(figsize=(8, 6))

  if dim==2:
    ax = fig.add_subplot(111, axisbg='#FFFFCC')
  elif dim==3:
    ax = Axes3D(fig)
  
  for i in range(len(cluster_coords)):
    xarray=numpy.ones((len(cluster_coords[i])))
    yarray=numpy.ones((len(cluster_coords[i])))
  
    if dim==3:  
      zarray=numpy.ones((len(cluster_coords[i])))

    #cluster_identifier=color_array[i%len(color_array)]+'o'
    cluster_identifier='o'    

    if i == len(cluster_coords)-1 or len(cluster_coords[i]) < 5:
      cluster_identifier='0.2'


    for j in range(len(cluster_coords[i])):
      xarray[j]=cluster_coords[i][j][0]
      yarray[j]=cluster_coords[i][j][1]
      if dim == 3:
        zarray[j]=cluster_coords[i][j][2]     

      if i == len(cluster_coords)-1 or len(cluster_coords[i]) < 5:
        if dim == 2:  
          ax.plot(xarray,yarray,cluster_identifier,marker=',',linewidth=0)      
        
        elif dim == 3:
          ax.plot(xarray,yarray,zarray,cluster_identifier,marker=',',linewidth=0)

      else:
        if dim == 2:
          ax.plot(xarray, yarray, 'o', color=color_array[i])
        elif dim == 3:
          ax.plot(xarray,yarray,zarray,'o',color=color_array[i])
  




  plt.show()

  



def __main__():
  
  plot=True
  embed=True
  tanimoto_matrix=False  

  parser=OptionParser()

  parser.add_option("-i","--input",action="store",dest="input_file")
  parser.add_option("-o","--output",action="store",dest="output_file")
  parser.add_option("-p","--only-plot",action="store_true",dest="only_plot",default=False)
  parser.add_option("-e","--only-embed",action="store_true",dest="only_embed",default=False)
  parser.add_option("-d","--dimension",type="int",dest="dimension",default=int(2))
  parser.add_option("-l","--loops",action="store",type="float",dest="loop_ratio",default=float(0.1))
  parser.add_option("--log",action="store_true",dest="write_log",default=False)
  parser.add_option("-t","--tanimoto-matrix",action="store_true",dest="only_tanimoto",default=False)

  (options,args)=parser.parse_args()


  if options.only_plot: embed=False
  if options.only_embed: plot=False
  if options.only_tanimoto: embed=False; plot=False; tanimoto_matrix=True;

  mol_list=list(read_sdf_file(options.input_file))  
  dim=options.dimension
  
  if tanimoto_matrix:
    
    calculate_tanimoto_matrix(mol_list)



  if embed:

    start_time=time.time()
    mol_list=SPE_Embed(mol_list,dim,options.loop_ratio)
    embed_time=time.time()-start_time
    

    

    outfile=pybel.Outputfile(format='sdf',filename=options.output_file,overwrite=True)
    for mol in mol_list:  
      outfile.write(mol)
    outfile.close()
    
    if options.write_log: 
      outname=options.output_file + '.log'
      f=open(outname,'w')
      f.write('Wallclock time: '+str(embed_time)+' seconds for '+str(len(mol_list))+' fingerprints')
      f.close()



  if plot:
    plot_spe(mol_list,dim)
    embed_time=0

  return len(mol_list),embed_time


data_size,timings=__main__()
print "Wallclock time:",timings,"seconds for",data_size,"fingerprints."
  












