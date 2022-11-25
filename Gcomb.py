#!/usr/bin/env python

# python Gcomb.py -t topol.top -r 1.5 -d 3

import sys
import os
import time
import shelve
import hashlib
import heapq
import argparse
import os.path
import functools
from functools import total_ordering
import shutil
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

class Node:
    def __init__(self, key, f):
        self.key = key  
        self.f = f 
    
    def __lt__(self, obj):
        return ((self.f) > (obj.f))
  
    def __gt__(self, obj):
        return ((self.f) < (obj.f))
  
    def __le__(self, obj):
        return ((self.f) >= (obj.f))
  
    def __ge__(self, obj):
        return ((self.f) <= (obj.f))
  
    def __eq__(self, obj):
        return (self.f == obj.f)

    
parser = argparse.ArgumentParser(description = 'Files required to run simulations')
parser.add_argument('-t','--top',type = str,required = 'True')
parser.add_argument('-r','--tol',type = float,required = 'True')
parser.add_argument('-d','--dep',type = int,required = 'True')

args = parser.parse_args()   

maxDepth = args.dep
tolerance = args.tol

def m_GenTpr(init,numseed,keyhash): 
    #gmx grompp -f unfold.mdp -c md_0_3.gro  -p topol.top -o md_0_4.tpr
    ######md version for gro
    os.system('gmx grompp -f fold'  + str(numseed) + '.mdp -c m_' + init +' -r m_' + init +' -p '+ args.top + ' -o m_' + keyhash +'.tpr   ')
  

def g_GenTpr(init,numseed,keyhash): 
    ########gmx version of gro
    os.system('gmx grompp -f fold'  + str(numseed) + '.mdp -c ' + init +' -r ' + init +' -p '+ args.top + ' -o g_' + keyhash +'.tpr   ')
    #gropath = os.path.realpath(init + '.gro')
    #os.system('gmx grompp -f fold'  + str(numseed) + '.mdp -c ' + gropath +' -r ' + gropath +' -p '+ args.top + ' -o ' + keyhash +'.tpr')

    
def getHash(mystr):
    coded =  mystr.encode()
    return hashlib.sha256(coded).hexdigest()

def addDict(name,atoms,g,h,f,depth):
    keys = ['name','atoms','g','h','f','depth']
    values = [name,atoms, g, h, f,depth]
    tdict = dict(zip(keys,values))
    return tdict

def out(command):
    result = os.popen(command).read()
    return float(result)

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

def Sim(parent,child,sTree,tpr): 
    
    #os.system("srun -T 8 singularity exec /nfshome/jphillips/images/gromacs.sif gmx mdrun -v -deffnm "+ tpr  +" -nice 19")
    os.system("gmx mdrun -v -deffnm g_"+ tpr  +" -nice 19 -nt 8")
    #remove tpr.log,edr,cpt,trr,
    os.remove('g_'+tpr + ".log")
    os.remove('g_'+tpr + ".edr")
    os.remove('g_'+tpr + ".cpt")
    os.remove('g_'+tpr + ".xtc")
    os.remove('g_'+tpr + ".tpr")
    tprfile = 'g_'+ tpr + '.gro'
    
    #u = mda.Universe(tprfile)
    #atoms = u.atoms
    #just keep tpr.gro
    #g = sTree[parent]['g'] + calcDist(sTree[parent]['name'], tprfile)
    
    #v = mda.Universe(sTree['target'])
    #targ = v.atoms
    
    #h = calcDist(tprfile,sTree['target'])
    #f = (g * sTree['alpha']) + h
    #depth =sTree[parent]['depth'] + 1
    grohash = getHash(tpr)     
   
    initname = 'g_' + tpr + '.gro'
    finalname = 'g_' + grohash + '.gro'
    os.rename(initname,finalname)
    #d = addDict(finalname,g,h,f,depth)

    #if tpr isn't there, skip run
    tprfile = 'm_'+ tpr + '.gro'
    os.system("gmx mdrun -v -deffnm m_"+ tpr  +" -nice 19")
    #remove tpr.log,edr,cpt,trr,
    os.remove('m_'+tpr + ".log")
    os.remove('m_'+tpr + ".edr")
    os.remove('m_'+tpr + ".cpt")
    os.remove('m_'+tpr + ".xtc")
    os.remove('m_'+tpr + ".tpr")
    
    u = mda.Universe(tprfile)
    atoms = u.atoms
    
    univ = mda.Universe(sTree['template'])
    p = univ.atoms
    
    temp = sTree[parent]['atoms']
    p.positions = temp[0]
    p.dimensions = temp[1]
    
    g = sTree[parent]['g'] + m_calcDist(p, atoms)
    
    univt = mda.Universe(sTree['template'])
    t = univt.atoms.select_atoms('protein')
    #t.positions = np.zeros(t.positions.shape)
    #t.dimensions = np.zeros(t.dimensions.shape)
    tempt = sTree['targetcoordinates']
    #print(t.positions.shape, tempt.shape)
    t.positions = tempt[0]
    t.dimensions = tempt[1]
    h = m_calcDist(atoms,t)
    f = (g * sTree['alpha']) + h
    depth =sTree[parent]['depth'] + 1
    #grohash = getHash(tpr)     
   
    initname = tpr
    finalname = grohash + '.gro'
    print("key", tprfile
    os.rename(initname,finalname)
    tree = (atoms.positions,atoms.dimensions)
    
    d = addDict(finalname,tree,g,h,f,depth)
    sTree[child] = d

def getChildrenKeys(parent_key,seed):
    return [parent_key + '_' + str(i)  for i in range(seed)]

def getChildrenHash(keys):               
    return [gethash(str(i)) for i in keys]


def getChildrenhval(children):
    return [m_calcDist(sTree[i]['atoms'] ,sTree['targetcoordinates']) for i in children]
  
def getChildrengval(parent,childrenkeys):
    return [(sTree[parent]['g'] + m_calcDist(sTree[parent]['atoms'], sTree[i]['atoms'])) for i in childrenkeys]


def addSortFrontier(pq,newNode):
    pq.append(newNode)
    temp = sorted(pq)
    return temp

def goalCheck(sTree,key):
    if sTree[key]['h'] <= tolerance:         #is this ok?
        return True
    else:
        return False
    
def setGoal(sTree,goalKey,tolerance):
    keys = ['key','tolerance']
    values = [goalKey,tolerance]
    tdict = dict(zip(keys,values))
    sTree['goal'] = tdict

def checkExist(children,sTree,i = 0):
    print(i)
    if i == sTree['seed']:
        return
    else:
        return (children[i] in sTree.keys() and checkExist(children,sTree, i+1))
    
def main():
    sTree = shelve.open("tree", writeback = True)
   
    curr_key = 'R'
    
    curr_depth = sTree[curr_key]['depth']
    frontier = []
    currNode = Node(curr_key,sTree[curr_key]['f'])
    frontier.append(currNode)
    foundGoal = False
    seed = sTree['seed']
        
    while len(frontier) != 0:   
        
        universe = mda.Universe(sTree['template'])
        myatoms = universe.atoms
        
        currNode = frontier.pop()
        curr_key = currNode.key
        curr_depth = sTree[curr_key]['depth']
        #goal check dist tolerance val h
        if goalCheck(sTree,curr_key):
            foundGoal = True
            goalKey = curr_key
            break
        if curr_depth > maxDepth:
            continue
        if currNode.key not in sTree.keys():
            print("Error")    
        children = getChildrenKeys(curr_key,seed)
        if(checkExist(children,sTree)):
            os.remove(sTree[curr_key]['gro'] )
        #print("children",children)
        for c,t in zip(children,range(seed)):
            if c in sTree.keys():
                N = Node(c,sTree[c]['f'])
                frontier = addSortFrontier(frontier,N)
            else:
                h = getHash(c)
                #writing required gro 
                gfilename = sTree[curr_key]['name'] 
                mfile = getHash(curr_key) + '.gro'
                #myatoms = sTree[curr_key]['gro']
                #myatoms.write(filename)
                #Generating tpr files for production
                mfilename ='m_'+ mfile
                temp =  sTree[curr_key]['atoms']
                myatoms.positions = temp[0]
                myatoms.dimensions = temp[1]
                myatoms.write(mfilename)
                g_GenTpr(gfilename,t,h)  
                m_GenTpr(mfile,t,h)  
                #removing gro
                #os.remove(mfilename)
                
                #Production step     
                Sim(curr_key,c,sTree,h)
                N = Node(c,sTree[c]['f'])
                frontier = addSortFrontier(frontier,N)
                print(c)

    for i in sTree.keys():
        print(i)
        
    if foundGoal == True:
        setGoal(sTree,goalKey,tolerance)
        #stop here then postproc
    #PRINTING TO CHECK HOW FRONTIER IS WORKING:
    curr_key = 'R'
    frontier = []
    currNode = Node(curr_key, sTree[curr_key]['f'])
    frontier.append(currNode)
    while len(frontier) != 0 :     
        currNode = frontier.pop()
        print(currNode.key, "--->\n", currNode.f, sTree[currNode.key]['g'],sTree[currNode.key]['h'])
        curr_key = currNode.key
        #print("H: ",sTree[curr_key]['gro'] )
        curr_depth = sTree[curr_key]['depth']
        if curr_depth > maxDepth:
            continue
        if currNode.key not in sTree.keys():
            print("Error") 
        children = getChildrenKeys(curr_key,seed)
        for c in children:
            if c in sTree.keys():
                N = Node(c,sTree[c]['f'])
                frontier = addSortFrontier(frontier,N)

    sTree.close()       

main()

