'''
python Preprocessor.py -u fold.mdp -g unfc.gro -t topol.top -p  targetu.gro -s 4 -a 0.1
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
  

def calcDist(initial,target):
    
    #add mda rmsd here
    m = mda.Universe(initial)
    ref = mda.Universe(target)

    mobile = m.atoms.select_atoms('backbone')
    reference = ref.atoms.select_atoms('backbone')
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

def addDict(grofile,g,h,f,depth):
    keys = ['gro','g','h','f','depth']
    values = [grofile, g, h, f,depth]
    tdict = dict(zip(keys,values))
    
    return tdict

def main():
    
    sTree = shelve.open("tree", writeback = True)
    
    #changing code to work with mda,keeping original for reference:
    

    #input from user, part of shelve but not tree
    #THIS OR UNIVERSE? 
    sTree['initialgro'] = args.gro
    
    sTree['topology'] = args.top
    
    #pdb didn't work with groreader, however it does work as a universe
    sTree['target'] = args.pdb
    sTree['seed'] = args.seed
    sTree['alpha'] = args.alpha
    a = sTree['alpha']
    
    #rmsd works with universes, so should we store the universe obj in shelve?
    g = calcDist(sTree['initialgro'], sTree['initialgro'])
    time.sleep(3)
    h = calcDist(sTree['initialgro'], sTree['target'])
    f = (g*a) +h
    depth = 1
    grohash = getHash(sTree['initialgro'])
    
    initname = sTree['initialgro']
    finalname = grohash + '.gro'
    os.rename(initname,finalname)
    sTree['initialgro'] = finalname
    
    grotarget = getHash(sTree['target'])
    
    initname = sTree['target']
    finalname = grotarget + '.gro'
    os.rename(initname,finalname)
    sTree['target'] = finalname
    
    #u = mda.Universe(sTree['initialgro'])
    #atoms = u.atoms            
    atoms = sTree['initialgro'];
    d = addDict(atoms,g,h,f,depth)
    sTree['R'] =  d
    print(sTree['R']['gro'])
                     
        
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



