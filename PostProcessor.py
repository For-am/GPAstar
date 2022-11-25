'''
python PostProcessor.py 

'''
#!/usr/bin/env python
import sys
import os
import time
import shelve
import hashlib
import heapq
import argparse
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

parser = argparse.ArgumentParser(description = 'test')

parser.add_argument('-g','--g',type = str,required = 'True') 


args = parser.parse_args()   
  
def main():
    
    sTree = shelve.open("tree", writeback = True)
    #goal_str = sTree['goal']['key']
    goal_str = args.g
    print(goal_str)
    
    keys = goal_str.split('_')
    key = ''
    t = []
    for k in keys:
        key = key + k + '_'
        t.append(key)
    
    for key in t:
        key = key.rstrip('_')
        print(key)
        #this gives us ech key in sequence
        
    
    
    #gmx trjconv -center -pbc mol -f md_0_10.gro -o unfolded.gro -s md_0_10.tpr
    sTree.close()
    
main()



