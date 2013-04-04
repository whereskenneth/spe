#!/usr/bin/env python

import re
import math
import numpy
import sys
import openbabel
import pybel
import os
import random
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from optparse import OptionParser
import pylab


def find_nearest(array,value):
  idx=(numpy.abs(array-value)).argmin()
  if value-array[idx] == 0 : return idx
  else: return -1



###############################################################################

def read_sdf_file(filename):

  basename,extname=os.path.splitext(str(filename))
  #uses os.path to split the filename into its base and extenstion

  return pybel.readfile(re.sub(r'[^\w]','',extname),filename)
  #return list of mol objects which can be manipulated in the main program.

###############################################################################


##############################################################################
def plot_spe(mol_list,dim,output):

  save_png=False

  if output: save_png=True


  min_clusters_to_plot=5
  
  
  first_cluster=True

  cluster_array=[]
  cluster_array.append([])
  cluster_array[0]=[numpy.array([]),numpy.array([])]
  
  singleton_array=[numpy.array([]),numpy.array([])]
  
  if dim==3: 
    singleton_array.append(numpy.array([]))
    cluster_array.append(nump.array([]))

  cluster_list=numpy.array([])

        
  for mol in mol_list: 
    
    x=float(mol.data['SPE.x'])
    y=float(mol.data['SPE.y'])
    if dim==3: z=float(mol.data['SPE.z'])

    try: 
      cluster_id=int(mol.data['clid'])
    except: 
      cluster_id=-1
    
    
    
    if cluster_id<0 :
      singleton_array[0]=numpy.append(singleton_array[0],x)
      singleton_array[1]=numpy.append(singleton_array[1],y)
      if dim==3: singleton_array[2]=numpy.append(singleton_array[2],z)

    else :

      if first_cluster:
        cluster_list=numpy.append(cluster_list,cluster_id)  
        cluster_array[0][0]=numpy.append(cluster_array[0][0],x)
        cluster_array[0][1]=numpy.append(cluster_array[0][1],y)
        if dim==3: cluster_array[0][2]=numpy.append(cluster_array[0][2],z)
        cluster_array.append([])
        first_cluster=False
          
      else:
        clindex=find_nearest(cluster_list,cluster_id)
        if clindex >= 0:
          cluster_array[clindex][0]=numpy.append(cluster_array[clindex][0],x)
          cluster_array[clindex][1]=numpy.append(cluster_array[clindex][1],y)          
          if dim == 3:cluster_array[clindex][2]=numpy.append(cluster_array[clindex][2],z)      
        else:          

          cluster_array[cluster_list.size]=[numpy.array([]),numpy.array([])]
          cluster_array[cluster_list.size][0]=numpy.append(cluster_array[cluster_list.size][0],x)
          cluster_array[cluster_list.size][1]=numpy.append(cluster_array[cluster_list.size][1],y)
          if dim == 3: 
            cluster_array[cluster_list.size].append(numpy.array([]))
            cluster_array[cluster_list.size][2]=numpy.append(cluster_array[cluster_list.size][2],z)
          cluster_list=numpy.append(cluster_list,cluster_id) 
          cluster_array.append([])






  # Now there is a cluster_array, and a singleton array, each with a list of x,y,z coordinates
  # in [0],[1],[2] repectively. cluster_array will be an array of clusters with numpy lists for coordinates

  cluster_array.pop()

#  color_array=numpy.array((1,1))
  color_array=[]
  cm=pylab.get_cmap('hsv')  
  cluster_counter=0
  for i in range(len(cluster_array)):
    if len(cluster_array[i][0]) > min_clusters_to_plot-1:
      cluster_counter+=1

  for i in range(cluster_counter):
    color_array.append(cm(1.*i/cluster_counter))
  

  fig = plt.figure(figsize=(16, 16))

  if dim==2:
    ax = fig.add_subplot(111, axisbg='#FFFFCC')
   # ax=fig.add_subplot(111)
  elif dim==3:
    ax = Axes3D(fig)


  xarray=singleton_array[0]  
  yarray=singleton_array[1]
  cluster_identifier='0.5'

  if dim==3: zarray=singleton_array[2]
 
  if dim==2:    
    ax.plot(xarray,yarray,cluster_identifier,marker=',',linewidth=0)      
        
  elif dim == 3:
    ax.plot(xarray,yarray,zarray,cluster_identifier,marker=',',linewidth=0)

  color_counter=0
  for i in range(len(cluster_array)):
    xarray=cluster_array[i][0]
    yarray=cluster_array[i][1]
    if dim==3: zarray=cluster_array[i][2]
    
    if len(cluster_array[i][0]) > min_clusters_to_plot -1 :
      cluster_identifier=color_array[color_counter]
      color_counter+=1
    else:
      cluster_identifier='0.5' 
      
    if dim==2:
      ax.plot(xarray,yarray,color=cluster_identifier,marker=',',linewidth=0)
    if dim==3:
      ax.plot(xarray,yarray,zarray,color=cluster_identifier,marker=',',linewidth=0)


  #pylab.ylim([-100,100])
  #pylab.xlim([-100,100])  

  if save_png:
    plt.savefig(output,dpi=600,orientation='landscape',format='png')

  else:
    plt.show()


def __main__():
  
  dim=0

  parser=OptionParser()

  parser.add_option("-i","--input",action="store",dest="input_file")
  parser.add_option("-o","--output",action="store",dest="output_file")
  parser.add_option("-d","--dimension",type="int",dest="dimension",default=int(2))
  parser.add_option("-p","--png",type="string",dest="png_outfile",default="")

  (options,args)=parser.parse_args()


  mol_list=read_sdf_file(options.input_file)


  if options.dimension > 3: dim = 3
  elif options.dimension < 2: dim=2



  print "Writing coordinates to ", options.png_outfile, "."
  
  plot_spe(mol_list,2,options.png_outfile)

  return 0


__main__()










