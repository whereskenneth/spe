#!/usr/bin/env python

import re
import openbabel
import pybel
import os
from optparse import OptionParser

###############################################################################

def read_sdf_file(filename):

  basename,extname=os.path.splitext(str(filename))
  #uses os.path to split the filename into its base and extenstion

  return pybel.readfile(re.sub(r'[^\w]','',extname),filename)
  #return list of mol objects which can be manipulated in the main program.

###############################################################################

  

def __main__():
  
  plot=True
  embed=True
  tanimoto_matrix=False  

  parser=OptionParser()

  parser.add_option("-i","--input",action="store",dest="input_file")
  parser.add_option("-o","--output",action="store",dest="output_file")
  parser.add_option("-t","--tag",action="store",type="string",dest="SDtag",default="")
  parser.add_option("-v","--value",action="store",dest="value",type="string",default="")
  

  (options,args)=parser.parse_args()

  mol_list=read_sdf_file(options.input_file)  

  

  outfile=pybel.Outputfile(format='sdf',filename=options.output_file,overwrite=True)
  for mol in mol_list:  
    mol.data[options.SDtag]=options.value    
    outfile.write(mol)
 

  outfile.close()
    


  return 0


__main__()












