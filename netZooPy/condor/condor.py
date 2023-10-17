import numpy as np
import pandas as pd
from igraph import *
from .timer import Timer

class condor_object:
    """
    Initialization of the condor object. The function gets a network
    in edgelist format as a path to a file or encoded in a pandas dataframe.
    Builds a condor_object with an edgelist,an igraph network, names of the
    targets and regulators.

    Note: The edgelist is assumed to contain a bipartite network.
    The program will relabel the nodes so that the edgelist represents a bipartite network anyway.
    It is on the user to know that the network they are using is suitable for the method.

    Parameters
    ------------

        network_file: str
            Path to file encoding an edgelist.
        sep: str
            Separator used in the file.
        index_col: int
            Column that stores the index of the edgelist. E.g. None, 0...
        header: int
            Row that stores the header of the edgelist. E.g. None, 0...
        dataframe: dataFrame
            Pandas DataFrame containing the edgelist. Use as alternative
            to filename
        silent: Silent mode, will not print anything on console.

    Returns
    ---------

        condor_object : object
            The object has the following attributes:
            - net: Contains the network edgelist as a pandas DataFrame.
            - graph: iGraph object containing the bipartite network.
            - reg_names: list of the nodes in the first column.
            - tar_names: list of the nodes in the second column.
            - index_dict: dictionary keeping track of the indices in the "graph" variable and the actual names of the nodes.
    
        warning:
            Condor uses iGraph for the node assignment. For the initialization of the
            assignment, there is some stochasticity involved, that can be somehow controlled 
            by setting the python random seed.
            
            For instance, running condor twice
            on the same dataset, in the same python process, might result in slightly different assignments. 
            To avoid this behavior, you can set the seed ( random.seed(0) ) before calling 
            condor.

                >>> import random
                >>> random.seed(1)
                >>> c1 = condor(...)
                >>> random.seed(1)
                >>> c2 = condor(...)

            In this case c1 and c2 are exactly the same.

            On the contrary, if condor is called twice during two different python calls, you will 
            have the exact same results, as the random seed will have resetted. 
            Instead, the stochasticity of the initial assignment can be kept by setting 
            a random seed at the beginning

                >>> import random
                >>> random.seed(random.randint(1,10000000)) 
                >>> condor(...)
    
    """

    def __init__(
        self, network_file=None, sep=",", index_col=0, header=0, dataframe=None,silent=False
    ):

        # Checks that either a dataframe or a path to an edgelist are provided, and not to both at the same time.
        assert (
            dataframe is not None or network_file is not None
        ), "Could not create condor object. Edgelist must be passed either as a pandas DataFrame or path to file."
        assert (
            dataframe is None or network_file is None
        ), "Could not create condor object. Either pass a DataFrame or a file path. Not both."

        if network_file is not None:
            # Reads the edgelist and creates a DataFrame
            self.net = pd.read_csv(
                network_file, sep=sep, index_col=index_col, header=header
            )

        if dataframe is not None:
            # Checks that the dataframe has the proper format.
            assert (
                dataframe.shape[1] < 4
            ), "DataFrame must have at most three columns. Two containing vertices and an additional one containing edgeweights."
            if dataframe.shape[1] == 3:
                assert pd.api.types.is_numeric_dtype(
                    dataframe.iloc[:, 2]
                ), "Weight column must contain numeric values."
            self.net = dataframe

        # Forces the edgelist to induce a bipartite network by renaming the columns.
        cnames = self.net.columns
        self.net = self.net.astype({cnames[0]: str, cnames[1]: str})
        self.net.iloc[:, 0] = "reg_" + self.net.iloc[:, 0]
        self.net.iloc[:, 1] = "tar_" + self.net.iloc[:, 1]
        # self.net.columns = ["V1","V2","weight"]
        self.silent = silent

        with Timer("Object creation:",self.silent):
            # Checks that the DataFrame gives rise to a well-defined bipartite network.
            assert not self.net.isnull().any().any(), "NaN values detected."
            assert not (
                "" in list(self.net.iloc[:, 0]) or "" in list(self.net.iloc[:, 1])
            ), "Empty strings detected."

            # Checks if weights are provided. If not, a column with 1s is added to the DataFrame.
            if self.net.shape[1] != 3:
                print("Unweighted network. Weights initialized as 1.")
                self.net["weight"] = 1

            self.net.columns = ["V1", "V2", "weight"]
            # Creates iGraph object from the DataFrame.
            print('hello')
            print(self.net)
            self.graph = Graph.DataFrame(self.net, directed=False, use_vids=False)
            # from igraph 0.10 add parameter: use_vids=False

            self.reg_names = sorted(set(self.net.iloc[:, 0]))
            self.tar_names = sorted(set(self.net.iloc[:, 1]))

            # By construction of the graph object, and the fact that reg_ is sorted in front of tar_
            # we have that the graph nodes are sorted by first the reg nodes and then the tar nodes.
            #types = [0] * len(self.reg_names)
            #types.extend([1] * len(self.tar_names))
            types = [0 if self.graph.vs[i]['name'].startswith('reg') else 1 for i in range(len(self.graph.vs))]
            
            self.graph.vs["type"] = types

            # Dictionary to keep track of node indices and node names should they be rearranged.
            # TO-DO: Check if this is really necessary.
            self.index_dict = {node.index: node["name"] for node in self.graph.vs}

            self.modularity = None
            self.reg_memb = None
            self.tar_memb = None
            self.Qcoms = None

    def initial_community(self, method="LDN", project=False,resolution=1):
        """
            Computation of the initial community structure based on unipartite methods.

        Parameters
        -----------

            method : str
                Method to determine intial community assignment.
                    - "LCS": Multilevel method.
                    - "LDN": Leiden method
            project: bool
                Whether to apply the initial community structure on the bipartite network
                disregarding the bipartite structure or apply it to the unipartite network resulting from the projection onto
                one of these nodes.
        Outputs:
            self     : updates condor object with
                     tar_memb: DataFrame of initial target node membership.
                     reg_memb: DataFrame of initial reg node membership.
        """

        weights_id = self.net.columns[2]
        if project:
            with Timer("Initial community structure with projection:",self.silent):

                # Obtains bipartite projection onto target subset.
                projected_graph = self.graph.bipartite_projection(which=1)

                if method == "LCS":
                    vc = Graph.community_multilevel(projected_graph, weights=weights_id)
                if method == "LDN":
                    vc = Graph.community_leiden(
                        projected_graph, objective_function='modularity',resolution_parameter=resolution, weights=weights_id
                    )

                self.modularity = vc.modularity
                if not self.silent: print("Initial modularity: ", self.modularity)

                tar_index = [i.index for i in self.graph.vs.select(type_in=[1])]
                # By the ordering on the indices, the target nodes indices begin at len(co.reg_names)
                # in order not to mess with the indices in vc we substract that starting value.
                tar_memb = [vc.membership[i - len(self.reg_names)] for i in tar_index]
                T0 = pd.DataFrame(zip(tar_index, tar_memb))
                T0.columns = ["index", "community"]
                self.tar_memb = T0

                # Here we have only computed the community structure for the target nodes. We initialize a properly sized dataframe for the reg nodes.
                reg_index = [i.index for i in self.graph.vs.select(type_in=[0])]
                # By the ordering on the indices, the target nodes indices begin at len(co.reg_names)
                # in order not to mess with the indices in vc we substract that starting value.
                reg_memb = [0 for i in reg_index]
                R0 = pd.DataFrame(zip(reg_index, reg_memb))
                R0.columns = ["index", "community"]
                self.reg_memb = R0

                # Sort of the same as above but without projecting. The bipartite network is treated as a unipartite network.
        else:
            with Timer("Initial community structure without projection:",self.silent):
                if method == "LCS":
                    vc = Graph.community_multilevel(self.graph, weights=weights_id)
                if method == "LDN":
                    vc = Graph.community_leiden(
                        self.graph,objective_function='modularity',resolution_parameter=resolution, weights=weights_id
                    )

                self.modularity = vc.modularity
                if not self.silent: print("Initial modularity: ", self.modularity)

                tar_index = [i.index for i in self.graph.vs.select(type_in=[1])]
                tar_memb = [vc.membership[i] for i in tar_index]
                tar_index = [self.index_dict[i] for i in tar_index]
                T0 = pd.DataFrame(zip(tar_index, tar_memb))
                T0.columns = ["index", "community"]
                self.tar_memb = T0

                reg_index = [i.index for i in self.graph.vs.select(type_in=[0])]
                reg_memb = [vc.membership[i] for i in reg_index]
                reg_index = [self.index_dict[i] for i in reg_index]
                R0 = pd.DataFrame(zip(reg_index, reg_memb))
                R0.columns = ["index", "community"]
                self.reg_memb = R0


    def bipartite_modularity(self, B, m, R, T):
        """ Computation of the bipartite modularity as described in 
        Michael J. Barber. Modularity and community detection in bipartite networks.

        Parameters
        ------------
            B        : array
                modularity matrix.
            m        : array
                sum of the weights (or number of edges in the unweighted case).
            R        : array
                community assignement matrix for reg nodes.
            T        : array
                community assignement matrix for tar nodes.
        
        Returns
        --------
            Q: _
                Modularity score.
        Notes
        --------
            self.Qcoms: _
                Modularity contribution by each community.
            self.modularity: _
                Modularity score.

        References
        -----------
        .. [1]__ Michael J. Barber. Modularity and community detection in bipartite networks.

        """

        RtBT = T.transpose().dot(B.dot(R))
        Qcoms = (1 / m) * (np.diagonal(RtBT))
        Q = sum(Qcoms)
        self.Qcoms = Qcoms[Qcoms > 0]
        self.modularity = Q
        return Q

    def matrices(self, c,resolution):
        """ Computation of modularity matrix and initial community matrix.

        Parameters
        ------------
            c        : int
                max number of communities.
        Returns
        ----------
            B        : array
                Modularity matrix.
            m        : int
                Sum of the weights (or number of edges in the unweighted case.
            T0       : array
                Initial community structure matrix for target nodes.
            R0       : array
                Initial community structure matrix for reg nodes.
            gn       : dict
                Index dictionary for tar node names.
            rg       : dict
                Index dictionary for reg node names.

        """

        with Timer("Matrix computation:",self.silent):

            # Dimensions of the matrix
            p = len(self.tar_names)
            q = len(self.reg_names)

            # Index dictionaries for the matrix. Note that this set of indices is different of that in the condor object (that one is for the igraph network.)
            rg = {self.reg_names[i]: i for i in range(0, q)}
            gn = {self.tar_names[i]: i for i in range(0, p)}

            # Computes weighted biadjacency matrix.
            A = np.matrix(np.zeros((p, q)))
            for edge in self.net.iterrows():
                A[gn[edge[1][1]], rg[edge[1][0]]] = edge[1][2]

            # Computes node degrees for the nodesets.
            ki = A.sum(1)
            dj = A.sum(0)
            # Computes sum of edges and bimodularity matrix.
            m = float(sum(ki))
            B = A - resolution*((ki @ dj) / m)

            # d = self.index_dict

            # Computation of initial modularity matrix for tar and reg nodes from the membership dataframe.
            T_ed = zip(
                [gn[j] for j in [i for i in self.tar_memb.iloc[:, 0]]],
                self.tar_memb.iloc[:, 1],
            )
            T0 = np.zeros((p, c))
            for edge in T_ed:
                T0[edge] = 1

            R_ed = zip(
                [rg[j] for j in [i for i in self.reg_memb.iloc[:, 0]]],
                self.reg_memb.iloc[:, 1],
            )
            R0 = np.zeros((q, c))
            for edge in R_ed:
                R0[edge] = 1

        return B, m, T0, R0, gn, rg

    def brim(self, deltaQmin="def", c="def", resolution=1):

        """ Implementation of the BRIM algorithm to iteratively maximize 
        bipartite modularity. Note that c is the maximum number of communities. 
        Dynamic choice of c is not yet implemented.

        Parameters
        ------------
            deltaQmin: str
                Difference modularity threshold for stopping the iterative process
            c        : int
                max number of communities.
            resolution: float
                
        
        Notes
        ---------
            Updates condor object with the following:
            - self.modularity: Modularity score for the final assignement.
            - self.tar_memb: Final community assignement for tar nodes.
            - self.reg_memb: Final community assignement for reg nodes.

        Note:
            c        : Has to be bigger than the number of communities given by the initial community assignement. Otherwise the program will crash. The default option gives room for 20% more communities which rarely fails.

        """

        if c == "def":
            c = int(len(self.tar_memb["community"].unique()) * 1.2)

        B, m, T0, R0, gn, rg = self.matrices(c,resolution)


        # Default deltaQmin.
        if deltaQmin == "def":
            deltaQmin = min(1 / m, 1e-5)

        with Timer("BRIM: ",self.silent):
            Qnow = 0
            deltaQ = 1
            p, q = B.shape
            while deltaQ > deltaQmin:
                # Right sweep
                Tp = T0.transpose().dot(B)
                R = np.zeros((q, c))
                am = np.array(np.argmax(Tp.transpose(), axis=1))
                for i in range(0, len(am)):
                    R[i, am[i][0]] = 1

                # Left sweep
                Rp = B.dot(R)
                T = np.zeros((p, c))
                am = np.array(np.argmax(Rp, axis=1))

                for i in range(0, len(am)):
                    T[i, am[i][0]] = 1
                T0 = T

                Qthen = Qnow
                Qnow = self.bipartite_modularity(B, m, R, T)
                deltaQ = Qnow - Qthen
                if not self.silent: print(Qnow)

            self.modularity = Qnow

        self.tar_memb = pd.DataFrame(
            list(zip(list(gn), [T0[i, :].argmax() for i in range(0, len(gn))]))
        )
        self.reg_memb = pd.DataFrame(
            list(zip(list(rg), [R0[i, :].argmax() for i in range(0, len(rg))]))
        )
        self.tar_memb.columns = ["tar", "community"]
        self.reg_memb.columns = ["reg", "community"]

    def qscores(self):
        """
            Computes the qscores (contribution of a vertex to its community modularity)
            for each vertex in the network.
        """

        c = 1 + max(self.reg_memb["com"])
        B, m, T, R, gn, rg = self.matrices(c)
        self.Qscores = {"reg_qscores": None, "tar_qscores": None}

        # Qscores for the targets:
        Rq = B.dot(R) / (2 * m)
        Qj = list()
        for j, r in self.tar_memb.iterrows():
            Qjh = Rq[j, r["com"]] / self.Qcoms[r["com"]]
            Qj.append(Qjh)
        self.Qscores["tar_qscores"] = self.tar_memb.copy()
        self.Qscores["tar_qscores"]["qscore"] = Qj

        # Qscores for the regulators:
        Tq = T.transpose().dot(B) / (2 * m)
        Qi = list()
        for i, r in self.reg_memb.iterrows():
            Qih = Tq[r["com"], i] / self.Qcoms[r["com"]]
            Qi.append(Qih)
        self.Qscores["reg_qscores"] = self.reg_memb.copy()
        self.Qscores["reg_qscores"]["qscore"] = Qi


def run_condor(
    network_file,
    sep=",",
    index_col=0,
    header=0,
    initial_method="LDN",
    initial_project=False,
    com_num="def",
    deltaQmin="def",
    resolution=1,
    return_output=False,
    tar_output="tar_memb.txt",
    reg_output="reg_memb.txt",
    silent = False
):
    """
        Computation of the whole condor process. It creates a condor object and runs all the steps of BRIM on it. The function outputs

        Note: The edgelist is assumed to contain a bipartite network. The program will relabel the nodes so that the edgelist represents a bipartite network anyway.
        It is on the user to know that the network they are using is suitable for the method.
    
    Parameters
    -----------
        network_file: str
            Path to file encoding an edgelist.
        sep: str
            Separator used in the file.
        index_col: int
            Column that stores the index of the edgelist. E.g. None, 0...
        header: int
            Row that stores the header of the edgelist. E.g. None, 0...
        initial_method: str
            Method to determine intial community assignment. (By default Leiden method).
        initial_project: bool
            Whether to project the network onto one of the bipartite sets for the initial community detection.
        com_num: int
            Max number of communities. It is recomended to leave this to default, otherwise if the initial community assignement is bigger the program will crash.
        deltaQmin: float
            Difference modularity threshold for stopping the iterative process.
        resolution: int
            Not yet implemented.
        return_output:  bool
            Whether the function returns the created condor object.
        tar_output: str
            Filename for saving the tar node final membership.
        reg_output: str
            Filename for saving the reg node final membership.
        silent: Run in silent mode

    Returns
    --------
        Files "tar_memb.txt" and "reg_memb.txt" encoding the final tar and reg node membership.
    """

    co = condor_object(network_file, sep, index_col, header,dataframe=None,silent=silent)

    co.initial_community(method=initial_method, project=initial_project,resolution=resolution)


    co.brim(deltaQmin, c=com_num, resolution=resolution)
    co.tar_memb.to_csv(tar_output)
    co.reg_memb.to_csv(reg_output)
    

    if return_output == True:
        return co
