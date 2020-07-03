import math
import numpy as np
from numpy import dtype
from itertools import *


    


    

def drawGraph(r,adj):
    #turns adjacency matrix into latex graph drawing
    #returns a string, which can be copied and pasted into latex
    #drawing is of all graph vertices arranged on a circle, equally spaced
    #adj is numpy matrix of graph to be drawn
    #r is radius of circle

    n = len(adj)
    s= "\\begin{tikzpicture}\n"
    locs = [0 for i in range(0,n)]
    for i in range(0,n):
        x = round(r*math.cos(2*math.pi * i / n),4)
        y = round(r*math.sin(2*math.pi * i / n),4)
        locs[i] = [str(x),str(y)]
        s1 = "\draw[fill=black] ("+str(x) + "," + str(y) + ") circle (3pt); \n"
        s = s + s1
        
    for i in range(0,n):
        for j in range(i,n):
            if(adj[i][j] == 1):
                s1 = "\draw[thin] (" + locs[i][0] + "," + locs[i][1] + ") -- (" + locs[j][0] + "," + locs[j][1] + ") ;\n"
                s = s + s1
    s = s + "\\end{tikzpicture}"
    return s

def testSrg(a):
    #a is numpy matrix, adjacency matrix of a graph
    #test if A satisfies A^2 = xA + yI + zJ
    #I is identity matrix, J is matrix of all ones
    #this condition is the strongly regular graph condition
    

    b = a @ a 
    
    x = 0
    y = 0
    z = 0
    i = 0
    j = 0
    while((i < len(a)) &(j < len(a))):
        
        if (i != j):
            if (a[i][j] == 0):
                z = b[i][j]
                break
        if (j == len(a)-1):
            i,j = i+1,0
        else:
            i,j = i,j+1
    i = 0
    j = 0
    while((i < len(a)) &(j < len(a))):
        if (i != j):
            if (a[i][j] == 1):
                x = b[i][j] - z
                break
        if (j == len(a)-1):
            i,j = i+1,0
        else:
            i,j = i,j+1
    y = b[0][0] - z
    ar = [y+z,x+z,z]
    
    return( [np.array_equal(b ,x*a + y * np.identity(len(a)) + z * np.ones((len(a),len(a)))), ar])

def testThreeSRG(a):
    #a is numpy matrix, adjacency matrix of a graph

    #test if A satisfies A^3 = xA + yI + zJ
    #I is identity matrix, J is matrix of all ones
    #this condition is the 3-strongly regular graph condition that I have been exploring
    b = a @ a @ a 
    
    x = 0
    y = 0
    z = 0
    i = 0
    j = 0
    while((i < len(a)) &(j < len(a))):
        
        if (i != j):
            if (a[i][j] == 0):
                
                z = b[i][j]
                break
        if (j == len(a)-1):
            i,j = i+1,0
        else:
            i,j = i,j+1
    i = 0
    j = 0
    while((i < len(a)) &(j < len(a))):
        if (i != j):
            if (a[i][j] == 1):
                
                x = b[i][j] - z
                break
        if (j == len(a)-1):
            i,j = i+1,0
        else:
            i,j = i,j+1
    y = b[0][0] - z
    
    ar = [y+z,x+z,z]
    return( [np.array_equal(b ,x*a + y * np.identity(len(a)) + z * np.ones((len(a),len(a)))), ar])

def testnSRG(n,a):
    #a is numpy matrix, adjacency matrix of a graph
    #n is integer

    #test if A satisfies A^n = xA + yI + zJ
    #I is identity matrix, J is matrix of all ones
    #this condition is a more general version of SRG and 3SRG conditions
    b = np.linalg.matrix_power(a,n)
    
    x = 0
    y = 0
    z = 0
    i = 0
    j = 0
    while((i < len(a)) &(j < len(a))):
        
        if (i != j):
            if (a[i][j] == 0):
                
                z = b[i][j]
                break
        if (j == len(a)-1):
            i,j = i+1,0
        else:
            i,j = i,j+1
    i = 0
    j = 0
    while((i < len(a)) &(j < len(a))):
        if (i != j):
            if (a[i][j] == 1):
                
                x = b[i][j] - z
                break
        if (j == len(a)-1):
            i,j = i+1,0
        else:
            i,j = i,j+1
    y = b[0][0] - z
    
    ar = [y+z,x+z,z]
    return( [np.array_equal(b ,x*a + y * np.identity(len(a)) + z * np.ones((len(a),len(a)))), ar])


def g6toAdjMat(i):
    #i is string of graph represented in g6 format
    #takes in graph in g6 format and returns numpy adjacency matrix of graph
    l = i
    s = ""
    nn = i[0]
    n = ord(nn)-63
    
    for j in l[1:]:
        l1 = (str(bin(ord(j)-63))[2:])
        while(len(l1) < 6):
            l1 = "0"+l1
        s = s + l1
    bitstr = (s[0:n*(n-1)//2])
    mat=  np.zeros((n,n), np.int32)
    for r in range(1,n):
        for c in range(0,r):
            mat[r][c] = int(bitstr[0])
            mat[c][r] = int(bitstr[0])
            bitstr = bitstr[1:]
    return mat

def adjMatTog6(mat):
    #mat is numpy matrix, adjacency matrix of graph
    #takes in adjacency matrix and returns a string of g6 format
    n = len(mat)
    N = chr(n+63)
    s = []
    for i in range(1,n):
        for j in range(0,i):
            s.append(str(mat[j][i]))
    while(len(s)%6 != 0):
        s.append("0")
        
    st = N
    while(len(s) > 0):
        t = s[:6]
        l = "".join(t)
        i = int(l,2)+63
        st = st + chr(i)
        s = s[6:]
    return st
        
        
    
    

def testConnected(a):
    #test if graph is connected
    #by setting diagonal equal to 1, we add a path from each vertex to itself
    #meaning powers of adjacency matrix now count all paths of length less than 
    #or equal to instead of equal to
    #this method is naive and inefficient, but was quick to code
    for i in range(0,len(a)):
        a[i][i] = 1
    b = np.linalg.matrix_power(a,len(a))
    for i in range(0,len(a)):
        a[i][i] = 0
    if (0 in b):
        return False
    return True

def spanTree(mat):
    #test if mat is not connected
    if (not testConnected(mat)):
        return -1
    #takes as input adjacency matrix of connected graph
    #returns adjacency matrix of spanning tree of graph
    
    n = len(mat)
    
    for r in range(1,n):
        for c in range(0,r):
            #find edge
            if(mat[r][c] == 1):
                #test if graph is still connected after removing edge rc
                mat[r][c] = 0
                mat[c][r] = 0
                nghbrs = [r]
                check2 = True
                success = False
                res = []
                while(check2):
                    
                    
                    temp = []
                    res = res + nghbrs
                    for i in nghbrs:
                        for j in range(0,n):
                            
                            if((mat[i][j] == 1) & (not(j in res))):
                                
                                temp.append(j)
                            
                    
                    nghbrs = temp
                    if(c in nghbrs):
                        check2 = False
                        success = True
                    elif(len(nghbrs) == 0):
                        check2 = False
                if(success):
                    # if graph is still connected, recurse on graph with edge removed
                    return spanTree(mat)
                else:
                    #if graph is not still connected, put edge back and search for different edge
                    mat[c][r] = 1
                    mat[r][c] = 1
    # base case of recursion, return graph which is minimally connected
    return mat

def floydWarshall(mat):
    #standard floyd-warshall, uses dp to find length of shortest path between any two points
    dists = [[100 for i in range(len(mat))] for j in range(len(mat))]
    for i in range(0,len(mat)):
        dists[i][i] = 0
        for j in range(0,len(mat)):
            if(mat[i][j]):
                dists[i][j] = 1
    for k in range(len(mat)):
        for i in range(len(mat)):
            for j in range(len(mat)):
                if(dists[i][j] > dists[i][k] + dists[k][j]):
                    dists[i][j] = dists[i][k] + dists[k][j]
    return dists

def testDistanceRegular(mat, m=3):
    #checks if a graph is distance regular
    #concept of distance regular explained on https://en.wikipedia.org/wiki/Distance-regular_graph
    
    dist = floydWarshall(mat)
    arr = [[-1] for k in range(m+1)]
    for i in range(len(mat)):
        for j in range(len(mat)):
            a = dist[i][j]
            tempdis = [[0 for m in range(m+1)] for n in range(m+1)]
            for k in range(len(mat)):
                
                b = dist[i][k]
                c = dist[j][k]
                
                tempdis[b][c] = tempdis[b][c] + 1
            if (arr[a] == [-1]):
                arr[a] = tempdis
            elif(arr[a] != tempdis):
                return False
    return True
                
def graphComplement(mat):
    #returns graph complement of a graph
    #graph complement is explained on:https://en.wikipedia.org/wiki/Complement_graph
    mat1 = np.zeros((len(mat),len(mat)), np.int32)
    for i in range(0,len(mat)):
        for j in range(0,len(mat)):
            if(i!=j):
                if(mat[i][j] == 0):
                    mat1[i][j] = 1
               
    return mat1
                
def testRegular(mat):
    #tests if a graph is regular
    #a regular graph is one where every vertex has the same degree
    s = sum(mat[0])
    for i in range(1,len(mat)):
        if(sum(mat[i]) != s):
            return False
    return True

def testVVregular(mat):
    if(len(mat)%2 == 1):
        return [False,[]]
    #test if a graph can be partitioned into v1 and v2, sets of vertices
    #such that elements of v1 have m neighbors in v1 and n neighbors in v2
    #and vice verse for elements in v2
    #this is equivalent to a simplified form of the concept of equitable partitions
    n = len(mat)
    for i in allSubsets(n):
        if(len(i) == n//2):
            check = True
            ct1 = 0
            ct2 = 0
            temp = -1000
            for j in range(0,n):
                ct1 = 0
                ct2 = 0
                
                for k in range(0,n):
                    if(mat[j][k]):
                        if((j in i) == (k in i)):
                            ct1 = ct1 + 1
                        else:
                            ct2 = ct2 + 1
                #print(ct1,ct2)
                #print(temp)
                if(temp != -1000):
                    if(ct1-ct2 != temp):
                        check = False
                        break
                else:
                    
                    temp = ct1-ct2
            if(check):
                return [True, i]
    return [False, []]


def countVVregular(mat):
    if(len(mat)%2 == 1):
        return 0
    #counts how many ways the graph can be partitioned into v1 and v2, as in testVVregular
    ct = 0
    n = len(mat)
    for i in allSubsets(n):
        if(len(i) == n//2):
            check = True
            ct1 = 0
            ct2 = 0
            temp = -1000
            for j in range(0,n):
                ct1 = 0
                ct2 = 0
                
                for k in range(0,n):
                    if(mat[j][k]):
                        if((j in i) == (k in i)):
                            ct1 = ct1 + 1
                        else:
                            ct2 = ct2 + 1
                #print(ct1,ct2)
                #print(temp)
                if(temp != -1000):
                    if(ct1-ct2 != temp):
                        check = False
                        break
                else:
                    
                    temp = ct1-ct2
            if(check):
                ct = ct + 1
    return ct

                    
            
def relabelAdjMat(mat, tup):
    #takes in an adjacency matrix, and returns a relabelled adjacency matrix where the vertices included in tup, a tuple, are the first rows of the matrix.
    tem = tup[:]
    for i in range(0,len(mat)):
        if(not(i in tem)):
            tem = tem+[i]
    ret = np.zeros((len(mat),len(mat)), np.int32)
    for i in range(0,len(mat)):
        for j in range(0,len(mat)):
            ret[i][j] = mat[tem[i]][tem[j]]
    return ret
                       
                            
def graphExpansionOperation(mat1,mat2,mat3):
    #returns a block matrix constructed from mat1,2,3 equal to:
    #[mat1   mat2]
    #[mat2^T mat3]
    n = len(mat1)
    ret = np.zeros((2*n,2*n), np.int32)
    for i in range(0,n):
        for j in range(0,n):
            if(mat1[i][j]):
                ret[i][j] = 1
                
            if(mat3[i][j]):
                
                ret[i+n][j+n] = 1
    for i in range(0,n):
        for j in range(0,n):
            if(mat2[i][j]):
                
                ret[i][j+n] = 1
                ret[j+n][i] = 1
    return ret

def complete(n):
    # n is integer
    # returns adjacency matrix of complete graph on n vertices
    ret = np.ones((n,n), np.int32)
    for i in range(n):
        ret[i][i] = 0
    return ret

def findIsomorphism(m1,m2):
    #m1, m2 are numpy adjacency matrices of graphs
    #attempts to find isomorphism between the two graphs m1 and m2
    #returns the a tuple representing the permutation that is an isomorphism
    n = len(m1)
    s = [i for i in range(0,n)]
    for i in permutations(s):
        p1 = np.zeros((n,n), np.int32)
        p2 = np.zeros((n,n), np.int32)
        for j in range(n):
            p1[j][i[j]] = 1
            p2[i[j]][j] = 1
        if(np.array_equal(m1, p1@m2@p2)):
            return(i)
    return []

def HammingCodeMat(n):
    #generates adjacency matrix of n-dimensional hamming code
    #this code is explained on https://en.wikipedia.org/wiki/Hamming_code
    #notice that this returns the hamming code minus the zero vector
    gens = []
    for i in range(0,n):
        v = [0 for i in range(0,2**n-1)]
        jum = 2**(n-1-i)
        ind = 0 
        par = 1
        while(ind < len(v)):
            for j in range(0,jum):
                if(ind + j < len(v)):
                    v[ind+j] = par
            if(par == 1):
                par = 0
            else:
                par = 1
            ind = ind + jum
        gens.append(v)
        
    
    m = []
    for i in allBitStrings(n)[1:]:
        t = []
        for j in range(len(i)):
            if(i[j]):
                t.append(gens[j])
        m.append(addVec(t))
        
    return np.array(m, np.int32)
    
def LatexAppropriateAdjmat(mat):
    #mat is numpy matrix
    #outputs a string of LaTeX code to represent the matrix
    s= "\\begin{pmatrix} "
    for i in mat:
        for j in i:
            s = s + str(j) + "&"
        s = s[:-1]
        s = s + "\\\\"
    s = s[:-1]
    s = s + "\\end{pmatrix}"
    return s
                    


    
    

                          
                        
def allBitStrings(n):
    #returns a list of all bitstrings of length n
    l = [[]]
    for i in range(0,n):
        l = [i + [0] for i in l] + [i + [1] for i in l]
    return l

def addVec(vs):
    #sums a list of vectors, considered as lists of numbers, of all the same length
    k = [0 for i in range(len(vs[0]))]
    for i in range(len(vs[0])):
        sum = 0
        for j in range(len(vs)):
            sum = sum + vs[j][i]
        k[i] = sum%2
    return k

def allSubsets(n):
    l = [[]]
    for i in range(0,n):
        #print(i)
        l1 = l[:]
        for j in range(len(l1)):
            if(len(l1[j]) < n/2):
                l1[j] = l1[j] + [i]
        l = l + l1
        
    return l
