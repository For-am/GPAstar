'''
python precomb.py -u md.mdp -g unfc.gro -t topol.top -p  targetu.gro -s 4 -a 0.1
'''
#!/usr/bin/env python
import sys
import os
import time
import shelve
import hashlib
import heapq
import argparse
import functools
import shutil
import numpy as np
from functools import total_ordering
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

parser = argparse.ArgumentParser(description = 'Files required to run simulations')

parser.add_argument('-u','--mdp',type = str,required = 'True') 
parser.add_argument('-g','--gro',type = str,required = 'True') 
parser.add_argument('-t','--top',type = str,required = 'True') 
parser.add_argument('-p','--pdb',type = str,required = 'True')     
parser.add_argument('-s','--seed',type = int,required = 'True')   
parser.add_argument('-a','--alpha',type = float,required = 'True')  

args = parser.parse_args()   
  

def g_calcDist(initial,target):
    
    #add mda rmsd here
    m = mda.Universe(initial)
    ref = mda.Universe(target)

    mobile = m.atoms.select_atoms('backbone')
    reference = ref.atoms.select_atoms('backbone')
    align.alignto(mobile,reference)
    
    d = rmsd(mobile.atoms.positions,reference.atoms.positions)

    return d

def m_calcDist(initial,target):

    mobile = initial.select_atoms('backbone')
    reference = target.select_atoms('backbone')
    align.alignto(mobile,reference) 
    d = rmsd(mobile.atoms.positions,reference.atoms.positions)
    return d
    '''
    force = "Backbone\nBackbone"
    os.system("echo 'Backbone\nBackbone' | gmx rms -s " +  initial + " -f " + target + " -o rmsd.xvg")
    time.sleep(5)
    d = out("cat rmsd.xvg|grep -v '#'|grep -v '@'|tr -s ' '|cut -f3 -d ' ' | tail -1")
    os.remove("rmsd.xvg")
    return d
    '''

def getHash(mystr):
    coded =  mystr.encode()
    return hashlib.sha256(coded).hexdigest()

def addDict(atoms,g,h,f,depth):
    keys = ['atoms','g','h','f','depth']
    values = [atoms, g, h, f,depth]
    tdict = dict(zip(keys,values))
    return tdict


def main():

    sTree = shelve.open("tree", writeback = True)
    
    src = args.gro
    dest = r"template.gro"
    path = shutil.copyfile(src,dest)
    sTree['template'] = 'template.gro'
    a = mda.Universe(args.gro)
    myatoms = a.atoms    
    
    
    #input from user, part of shelve but not tree
    #storing positions and box vectors separately
    sTree['initialgrocoordinates'] = (myatoms.positions,myatoms.dimensions)
    sTree['initialgroname'] = args.gro
    sTree['topology'] = args.top
    
    t = mda.Universe(args.pdb)
    targ = t.atoms
   
    sTree['targetcoordinates'] =  (targ.positions,targ.dimensions)
    sTree['targetname'] = args.pdb
    sTree['seed'] = args.seed
    sTree['alpha'] = args.alpha
    a = sTree['alpha']
    
    g = g_calcDist(sTree['initialgroname'], sTree['initialgroname'])
    time.sleep(3)
    h = g_calcDist(sTree['initialgroname'], sTree['targetname'])
    f = (g*a) +h
    depth = 1
    grohash = getHash('R')
    
    initname = sTree['initialgroname']
    finalname = 'g_' + grohash + '.gro'
    os.rename(initname,finalname)
    sTree['initialgroname'] = finalname
    
    grotarget = getHash(sTree['targetname'])
    
    initname = sTree['targetname']
    finalname = 'g_' + grotarget + '.gro'
    os.rename(initname,finalname)
    sTree['targetname'] = finalname
    
    #u = mda.Universe(sTree['initialgro'])
    #atoms = u.atoms            
    name = sTree['initialgroname'];
    atoms = (sTree['initialgrocoordinates'])  
    d = addDict(atoms,g,h,f,depth)
    sTree['R'] =  d
                     
        
    #The following code generates the mdp files according to provided seeds
    file = open(args.mdp)
    lines = [""]
    text = str(file.readline())
    while text:
        lines.append(text)
        text = str(file.readline())
   
    foundv = False
    founds = False
    
    for i in range(sTree['seed']):
        fname = 'fold' + str(i) + '.mdp'
        f = open(fname, 'w')
        for l in lines:
            if l.startswith('gen_vel'):
                l = l.replace("gen_vel                 = no ","gen_vel                 = yes  \n")
                foundv = True
            if l.startswith('gen_seed'):
                l = l.replace(l,"gen_seed                = " + str(i)+ "\n")
                founds = True
            f.write(l)
        if foundv == False:
            f.write("gen_vel                 = yes  \n")
        if founds == False:
            f.write("gen_seed                = " + str(i)+ "\n")

           
    sTree.close()
    
main()



