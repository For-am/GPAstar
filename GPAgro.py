#!/usr/bin/env python

# python GPASearch.py -t topol.top -r 1.5 -d 3
# python GPASearch.py -t topol.top -r 3 -d 4
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

def GenTpr(init,numseed,keyhash): 
    print("11111",init, keyhash)
    #gmx grompp -f unfold.mdp -c md_0_3.gro  -p topol.top -o md_0_4.tpr
    os.system('gmx grompp -f fold'  + str(numseed) + '.mdp -c ' + init +' -r ' + init +' -p '+ args.top + ' -o ' + keyhash +'.tpr   ')
    #gropath = os.path.realpath(init + '.gro')
    #os.system('gmx grompp -f fold'  + str(numseed) + '.mdp -c ' + gropath +' -r ' + gropath +' -p '+ args.top + ' -o ' + keyhash +'.tpr')

def getHash(mystr):
    coded =  mystr.encode()
    return hashlib.sha256(coded).hexdigest()

def addDict(grofile,g,h,f,depth):
    keys = ['gro','g','h','f','depth']
    values = [grofile, g, h, f,depth]
    tdict = dict(zip(keys,values))
    return tdict

def out(command):
    result = os.popen(command).read()
    return float(result)

def calcDist(initial,target):
    
    #add mda rmsd here
    m = mda.Universe(initial)
    ref = mda.Universe(target)

    mobile = m.atoms.select_atoms('backbone')
    reference = ref.atoms.select_atoms('backbone')
    align.alignto(mobile,reference)
    
    d = rmsd(mobile.atoms.positions,reference.atoms.positions)

    return d

def Sim(parent,child,sTree,tpr): 
    
    #os.system("srun -T 8 singularity exec /nfshome/jphillips/images/gromacs.sif gmx mdrun -v -deffnm "+ tpr  +" -nice 19")
    os.system("gmx mdrun -v -deffnm "+ tpr  +" -nice 19 -nt 8")
    #remove tpr.log,edr,cpt,trr,
    os.remove(tpr + ".log")
    os.remove(tpr + ".edr")
    os.remove(tpr + ".cpt")
    os.remove(tpr + ".xtc")
    os.remove(tpr + ".tpr")
    tprfile = tpr + '.gro'
    
    #u = mda.Universe(tprfile)
    #atoms = u.atoms
    #just keep tpr.gro
    g = sTree[parent]['g'] + calcDist(sTree[parent]['gro'], tprfile)
    
    #v = mda.Universe(sTree['target'])
    #targ = v.atoms
    
    h = calcDist(tprfile,sTree['target'])
    f = (g * sTree['alpha']) + h
    depth =sTree[parent]['depth'] + 1
    grohash = getHash(tpr)     
   
    initname = tpr + '.gro'
    finalname = grohash + '.gro'
    os.rename(initname,finalname)
    d = addDict(finalname,g,h,f,depth)
    
    sTree[child] =  d
    '''
    cwd = os.getcwd()
    parent = os.path.realpath(cwd)
    
    directory = "groFiles"
    
    destination = os.path.join(parent, directory)
    source = os.path.realpath(finalname)
    shutil.move(source, destination)'''
    
def getChildrenKeys(parent_key,seed):
    return [parent_key + '_' + str(i)  for i in range(seed)]

def getChildrenHash(keys):               
    return [gethash(str(i)) for i in keys]


def getChildrenhval(children):
    return [calcDist(sTree[i]['gro'] ,sTree['target']) for i in children]
  
def getChildrengval(parent,childrenkeys):
    return [(sTree[parent]['g'] + calcDist(sTree[parent]['gro'], sTree[i]['gro'])) for i in childrenkeys]


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
    print(sTree[curr_key]['gro'])
    curr_depth = sTree[curr_key]['depth']
    frontier = []
    currNode = Node(curr_key,sTree[curr_key]['f'])
    frontier.append(currNode)
    foundGoal = False
    seed = sTree['seed']
        
    while len(frontier) != 0:   
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
                filename = sTree[curr_key]['gro'] 
                #myatoms = sTree[curr_key]['gro']
                #myatoms.write(filename)
                #Generating tpr files for production
                #print("00000", filename, curr_key)
                GenTpr(filename,t,h)  
                #removing gro
                #os.remove(filename)
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

