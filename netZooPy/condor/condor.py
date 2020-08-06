import numpy as np
import pandas as pd
from igraph import *
import time
"""
Description:
    Python implementation of the BRIM algorithm for bipartite community structure detection as described in 
    "Modularity and community detection in bipartite networks" by Michael J. Barber." This package is somewhat translated from the R package CONDOR https://github.com/jplatig/condor.

Usage:
    To use the package you can import it using: ```import condor```. And then access the different functions implemented with ```condor.function()```.
    Usage instructions:
    Supose you have a network (weighted or not) as an edgelist loaded into a pandas dataframe.
    ```
    import condor
    co = condor.condor_object(net)
    ```
    Returns a condor object (dictionary with a graph, node names, community membership dataframes...)
    ```
    co = condor.initial_community(co)
    ```
    Computes the initial community structure and updates the condor object.
    ```
    co = condor.brim(co)
    ```
    Runs the iterative modularity optimization algorithm (BRIM) and updates the condor object with the final membership.
    To see the results type:
    ```
    co["tar_memb"]
    co["reg_memb"]
    ```
    To compute the qscores for the vertices type:
    ```
    co = condor.qscores(co) # Computes the qscores
    co["qscores"]["reg_qscores"] # Dataframe containing the qscores for the regulators.
    co["qscores"]["tar_qscores"] # Dataframe containing the qscores for the targets.
"""
def condor_object(net):
    """
    Description:
        Initialization of the condor object. The function gets a network in edgelist format encoded in a pandas dataframe.
        Returns a dictionary with an igraph network, names of the targets and regulators, list of edges, modularity, and vertex memberships.   

    Inputs:
        net: Input adjacency matrix.

    Outputs:
        CO: CONDOR initial object.
    """
    t = time.time()
        #Error flags.
    assert len(set(net.iloc[:,0]).intersection(net.iloc[:,1]))==0, "The network must be bipartite."
    assert not net.isnull().any().any(), "NaN values detected."
    assert not ("" in list(net.iloc[:,0]) or "" in list(net.iloc[:,1])), "Empty strings detected."
    
    #Builds graph object.
    if (net.shape[1] == 3): 
        print("Weights detected")
        edges = list(zip(net.iloc[:,0],net.iloc[:,1],net.iloc[:,2]))
        Gr = Graph.TupleList(edges,weights=True)
    else:
        print("Unweighted network. Weights initialized as 1.")
        edges = list(zip(net.iloc[:,0],net.iloc[:,1],[1 for i in net.iloc[:,1]]))
        Gr = Graph.TupleList(edges,weights=True)
    
    #Assigns color names (bipartite sets).
    reg_names = sorted(set(net.iloc[:,1]))
    tar_names = sorted(set(net.iloc[:,0]))
    Gr.vs["type"] = 0 #Tar
    for j in [i.index for i in Gr.vs if i["name"] in reg_names]:
        Gr.vs[j]["type"] = 1 #Reg
    
    index_dict = {k.index:k["name"] for k in Gr.vs  }
    
    print("Condor object built in",time.time()-t)
    CO = {"G":Gr,"tar_names":tar_names,"reg_names":reg_names,"index_dict":index_dict,"edges":edges,"modularity":None,"reg_memb":None,"Qcoms":None}
    return CO

def bipartite_modularity(B,m,R,T,CO):  
    """
    Description:
        Computation of the bipartite modularity as described in ""Modularity and community detection in bipartite networks" by Michael J. Barber." 
        
    Inputs:
        B : Adjacency matrix of the network - adjacency matrix of the null model
        m : Modularity matrix.
        R : Indices of the target nodes.
        T : Indices of the source nodes.
        CO: CONDOR object with initial community.

    Outputs:
        Q : Modularity.
        CO: Updated CONDOR object.
    """
    RtBT = R.transpose().dot(B.dot(T))
    Qcoms = (1/m)*(np.diagonal(RtBT))
    Q = sum(Qcoms)
    Qcoms = Qcoms[Qcoms>0]
    CO["Qcoms"] = Qcoms
    return Q,CO

def initial_community(CO,method="LCS",project=False):
    """
    Description:
        Computation of the initial community structure based on unipartite methods.
        The implementation using bipartite projection is not yet available, but project=False performs better modularity-wise (at least with the networks I worked on).

    Inputs:
        CO     : CONDOR object.
        method : Method to determine intial community assignment.
                 "LCS": Multilevel method.
                 "LEC": Leading Eigenvector method.
                 "FG" : Fast Greedy method.

    Outputs:
        CO: Updated CONDOR object.
    """
    if project: print("Not yet implemented.")
    
    t = time.time()
    
    #Computes initial community assignement based on different methods. Default method is Louvain clustering which is really fast. The others take several times more.
    if not project:
        if method=="LCS":
            vc = Graph.community_multilevel(CO["G"],weights="weight")
        if method=="LEC":
            vc = Graph.community_leading_eigenvector(CO["G"],weights="weight")
        if method=="FG":
            vc = Graph.community_fastgreedy(CO["G"],weights="weight").as_clustering()
    print("Initial community structure computed in ",time.time()-t,". Modularity = ",vc.modularity)
     
    CO["modularity"] = vc.modularity  
    
    #Stores initial community in the condor object.
    reg_index = [i.index for i in CO["G"].vs.select(type_in=[1])]
    reg_memb = [vc.membership[i] for i in reg_index]
    T0 = pd.DataFrame(zip(reg_index,reg_memb))
    T0.columns = ["index","community"]    
    CO["reg_memb"] = T0
    
    return CO 

def brim(CO,deltaQmin="def",c=25):
    """
    Description:
        Implementation of the BRIM algorithm to iteratively maximize bipartite modularity.
        Note that c is the maximum number of communities. Dynamic choice of c is not yet implemented.

    Inputs:
        CO       : Initial CONDOR Object with initial community assignment.
        deltaQmin: Modularity stopping criterion.
        c        : The maximum number of communities.

    Outputs:
        CO: Updated CONDOR object with community assignment.
    """    
    #Gets modularity matrix, initial community matrix and index dictionary.
    B,m,T0,R0,gn,rg = matrices(CO,c)
    
    #Default deltaQmin.
    if(deltaQmin == "def"):
        deltaQmin = min(1/len(CO["edges"]),1e-5)
    
    #BRIM iterative process
    Qnow = 0
    deltaQ = 1
    p,q = B.shape
    while(deltaQ>deltaQmin):
        #Right sweep
        Tp = B.dot(T0)
        R = np.zeros((p,c))
        am = np.array(np.argmax(Tp, axis=1))
        for i in range(0,len(am)):
            R[i,am[i][0]]=1
        #Left sweep
        Rp = B.transpose().dot(R)
        T = np.zeros((q,c))
        am = np.array(np.argmax(Rp, axis=1))
        for i in range(0,len(am)):
            T[i,am[i][0]]=1
        T0 = T
        
        Qthen = Qnow
        Qnow,CO = bipartite_modularity(B,m,R,T,CO)
        deltaQ = Qnow-Qthen
        print(Qnow)
        
    #Update modularity attribute   
    CO["modularity"] = Qnow
    
    #Update membership dataframes.   
    CO["tar_memb"] = pd.DataFrame(list(zip(list(gn),[R[i,:].argmax() for i in range(0,len(gn))])))
    CO["reg_memb"] = pd.DataFrame(list(zip(list(rg),[T[i,:].argmax() for i in range(0,len(rg))]))) 
    CO["tar_memb"].columns = ["tar","com"]
    CO["reg_memb"].columns = ["reg","com"]
        
    return CO

def matrices(CO,c):
    """
    Description:
        Computation of modularity matrix and initial community matrix.

    Inputs:
        CO: Initial CONDOR object.
        c : The maximum number of communities.

    Outputs:
        B : Adjacency matrix of the network - adjacency matrix of the null model
        m : Modularity matrix.
        T0: Indices of the source nodes, for igraph.
        R0: Indices of the target nodes, for igraph.
        gn: Indices of the target node labels.
        rg: Indices of the source node labels.
    """
    t = time.time()
    
    #Dimensions of the matrix
    p = len(CO["tar_names"])
    q = len(CO["reg_names"])
    
    #Index dictionaries for the matrix. Note that this set of indices is different of that in the CO object (that one is for the igraph network.)
    rg = {CO["reg_names"][i]: i for i in range(0,q)}
    gn = {CO["tar_names"][i]: i for i in range(0,p)}
    
    #Computes bipartite adjacency matrix.
    A = np.matrix(np.zeros((p,q)))
    for edge in CO["edges"]:
        A[gn[edge[0]],rg[edge[1]]] = edge[2]
    
    #Computes bipartite modularity matrix.
    ki = A.sum(1)
    dj = A.sum(0)   
    m = float(sum(ki))
    B = A - ((ki@dj)/m)
    
    #Creates initial community T0 matrix.
    d = CO["index_dict"]
    if ("index" in CO["reg_memb"].columns): 
        ed = zip([rg[j] for j in [d[i] for i in CO["reg_memb"].iloc[:,0]]],CO["reg_memb"].iloc[:,1])
    else: 
        ed = zip([rg[j] for j in CO["reg_memb"].iloc[:,0]],CO["reg_memb"].iloc[:,1])
    T0 = np.zeros((q,c))
    for edge in ed:
        T0[edge]=1
    
    if ("tar_memb" not in CO):
        print("Matrices computed in",time.time()-t)
        return B,m,T0,0,gn,rg
    ed = zip([gn[j] for j in CO["tar_memb"].iloc[:,0]],CO["tar_memb"].iloc[:,1])
    R0 = np.zeros((p,c))
    for edge in ed:
        R0[edge]=1
    print("Matrices computed in",time.time()-t)
    return B,m,T0,R0,gn,rg

def qscores(CO):
    """
    Description:
        Computes the qscores (contribution of a vertex to its community modularity)
        for each vertex in the network.

    Inputs:
        CO: Initial CONDOR object.

    Outputs:
        CO: CONDOR object with qscores field.
    """
    c = 1 + max(CO["reg_memb"]["com"])
    B,m,T,R,gn,rg = matrices(CO,6)
    CO["Qscores"]={"reg_qscores":None,"tar_qscores":None}
    
    #Qscores for the regulators:
    Rq = R.transpose().dot(B)/(2*m)
    Qj = list()
    for j,r in CO["reg_memb"].iterrows():
        Qjh = Rq[r["com"],j]/CO["Qcoms"][r["com"]]
        Qj.append(Qjh)
    CO["Qscores"]["reg_qscores"]=CO["reg_memb"].copy()
    CO["Qscores"]["reg_qscores"]["qscore"] = Qj
    
    #Qscores for the targets:
    Tq = B.dot(T)/(2*m)
    Qi = list()
    for i,r in CO["tar_memb"].iterrows():
        Qih = Tq[i,r["com"]]/CO["Qcoms"][r["com"]]
        Qi.append(Qih)
    CO["Qscores"]["tar_qscores"]=CO["tar_memb"].copy()
    CO["Qscores"]["tar_qscores"]["qscore"] = Qi
    
    return CO

def condor(filename,c=25,deltaQmin="def"):
    """
    Description:
        Default settings run of condor with output of the membership dataframes.
        Reads a network in csv format and index_col=0.

    Inputs:
        filename : path of file to save community assignment.
        c        : The maximum number of communities.
        deltaQmin: Modularity stopping criterion.

    Reference:
        Platig, John, et al. "Bipartite community structure of eQTLs." PLoS computational biology 12.9 (2016): e1005033.
    """
    t = time.time()
    net = pd.read_csv(filename,sep=",",index_col=0)
    CO = condor_object(net)
    CO = initial_community(CO)
    CO = brim(CO,deltaQmin,c)
    
    CO["tar_memb"].to_csv("tar_memb.txt")
    CO["reg_memb"].to_csv("reg_memb.txt")
    
    return "Runtime "+str(time.time()-t)+"s"