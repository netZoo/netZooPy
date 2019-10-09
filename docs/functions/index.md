# Functions

## Panda

### run_panda

    Description:

        Run PANDA algorithm from the command line.

    Usage:

        -h, --help: help
        -e, --expression: expression values
        -m, --motif: pair file of motif edges, or Pearson correlation matrix when not provided
        -p, --ppi: pair file of PPI edges
        -o, --out: output file
        -r, --rm_missing
        -q, --lioness: output for Lioness single sample networks

    Example:

        python run_panda.py -e ./ToyData/ToyExpressionData.txt -m ./ToyData/ToyMotifData.txt -p ./ToyData/ToyPPIData.txt -o test_panda.txt -q output_panda.txt

### run_lioness

    Description:

        Run LIONESS algorithm from the command line.

    Usage:

        -h, --help: help
        -e, --expression: expression matrix (.npy)
        -m, --motif: motif matrix, normalized (.npy)
        -p, --ppi: ppi matrix, normalized (.npy)
        -n, --npy: PANDA network (.npy)
        -o, --out: output folder
        -f, --format: output format (txt, npy, or mat)
        start: to start from nth sample (optional)
        end: to end at nth sample (optional, must with start)

    Example:

        python run_lioness.py -e expression.npy -m motif.npy -p ppi.npy -n panda.npy -o /tmp -f npy 1 100

### panda

#### class Panda

    Panda(object)

    Description:

    	Using PANDA to infer gene regulatory network.

    Usage:

    	1. Reading in input data (expression data, motif prior, TF PPI data)
    	2. Computing coexpression network
    	3. Normalizing networks
    	4. Running PANDA algorithm
    	5. Writing out PANDA network

    Authors: 

	cychen, davidvi, alessandromarin

##### init
	init__(self, expression_file, motif_file, ppi_file, save_memory = False, save_tmp=True, remove_missing=False, keep_expression_matrix = False):
        # =====================================================================
        # Data loading
        # =====================================================================
       
##### remove_missing
	remove_missing(self)
        Remove genes and tfs not present in all files.

##### normalize_network
	normalize_network(self, x)

##### panda_loop
	panda_loop(self, correlation_matrix, motif_matrix, ppi_matrix)
        Panda algorithm.

###### t_function
        t_function(x, y=None)
	T function.

###### update_diagonal
	update_diagonal(diagonal_matrix, num, alpha, step)
        Update diagonal.

##### pearson_results_data_frame
	pearson_results_data_frame(self)
        Results to data frame.

##### save_panda_results
	save_panda_results(self, path='panda.npy')

##### top_network_plot
	top_network_plot(self, top = 100, file = 'panda_top_100.png')
        Select top genes.

##### shape_plot_network
	shape_plot_network(self, subset_panda_results, file = 'panda.png')
        Create plot.

##### create_plot
	create_plot(self, unique_genes, links, file = 'panda.png')
        Run plot.

###### split_label
	split_label(label)

##### return_panda_indegree
	return_panda_indegree(self)
        Return Panda indegree.

##### return_panda_outdegree
	return_panda_outdegree(self)
        Return Panda outdegree.

## Lioness

### lioness

#### class Lioness

    Lioness(Panda)
 
    Description:

       Using LIONESS to infer single-sample gene regulatory networks.

    Usage:

       1. Reading in PANDA network and preprocessed middle data
       2. Computing coexpression network
       3. Normalizing coexpression network
       4. Running PANDA algorithm
       5. Writing out LIONESS networks

    Authors:

       cychen, davidvi

##### init

    init__(self, obj, start=1, end=None, save_dir='lioness_output', save_fmt='npy')

##### lioness_loop

    lioness_loop(self)

##### save_lioness_results

    save_lioness_results(self, file='lioness.txt')

### AnalyzeLioness

#### class AnalyzeLioness

    AnalyzeLioness(Lioness)

##### init

    init__(self, lioness_data)
    Load variables from lioness.

##### top_network_plot

    top_network_plot(self, column = 0, top = 100, file = 'lioness_top_100.png')
    Select top genes.

### lioness_for_puma

#### class LionessPuma

    LionessPuma(Puma)

    Description:
         Using LIONESS to infer single-sample gene regulatory networks.

    Usage:
        1. Reading in PUMA network and preprocessed middle data
        2. Computing coexpression network
        3. Normalizing coexpression network
        4. Running PUMA algorithm
        5. Writing out LIONESS networks

    Authors:
        cychen, davidvi

##### init

    init__(self, obj, start=1, end=None, save_dir='lioness_output', save_fmt='npy')

##### lioness_loop

     lioness_loop(self)

##### save_lioness_results

     save_lioness_results(self, file='lioness.txt')

## Puma

### run_puma

    Description:
        Run PUMA algorithm from the command line.

    Usage:
        run_puma
        -h, --help: help
        -e, --expression: expression values
        -m, --motif: pair file of motif edges, or Pearson correlation matrix when not provided
        -p, --ppi: pair file of PPI edges
        -i, --mir (required): miR file
        -o, --out: output file
        -r, --rm_missing
        -q, --lioness: output for Lioness single sample networks

     Example:
         python run_puma.py -e ./ToyData/ToyExpressionData.txt -m ./ToyData/ToyMotifData.txt -p ./ToyData/ToyPPIData.txt -i ToyData/ToyMiRList.txt -o test_puma.txt -q output_lioness.txt

### puma

#### class Puma

    Puma(object)

    Description:
        Using PUMA to infer gene regulatory network.

    Usage:
        1. Reading in input data (expression data, motif prior, TF PPI data, miR)
        2. Computing coexpression network
        3. Normalizing networks
        4. Running PUMA algorithm
        5. Writing out PUMA network

    Authors:
        cychen, davidvi, alessandromarin

##### init

    init__(self, expression_file, motif_file, ppi_file, mir_file, save_memory = False, save_tmp=True, remove_missing=False, keep_expression_matrix = False)
    # =====================================================================
    # Data loading
    # =====================================================================

##### remove_missing

    remove_missing(self)
    Remove genes and tfs not present in all files.

##### normalize_network

    normalize_network(self, x)

##### puma_loop

    puma_loop(self, correlation_matrix, motif_matrix, ppi_matrix)
    Puma algorithm

###### t_function

    t_function(x, y=None)
    T function

###### update_diagonal

    update_diagonal(diagonal_matrix, num, alpha, step)
    Update diagoanl

##### pearson_results_data_frame

    pearson_results_data_frame(self)
    Results to data frame

##### pearson_results_data_frame

    pearson_results_data_frame(self)
    Results to data frame.

##### save_puma_results

    save_puma_results(self, path='puma.npy')

##### top_network_plot

    top_network_plot(self, top = 100, file = 'puma_top_100.png')
    Select top genes.

##### shape_plot_network

    shape_plot_network(self, subset_puma_results, file = 'puma.png')
    Create plot

##### create_plot

    create_plot(self, unique_genes, links, file = 'puma.png')
    Run plot

###### split_label

    split_label(label)

##### return_puma_indegree

    return_puma_indegree(self)
    Return Puma indegree

##### return_puma_outdegree

     return_puma_outdegree(self)
     Return Puma outdegree
