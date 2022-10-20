import numpy as np
import numpy.random as ra
import pandas as pd
N = 30
P = 10

layer1_data = ra.multivariate_normal(mean = np.zeros(shape=P), cov = np.identity(P), size = N)
layer2_data = ra.multivariate_normal(mean = np.zeros(shape=P), cov = np.identity(P), size = N)

# create id column
ids = np.arange(N).reshape((N, 1))

# creat colnames
colnames = ["gene" + str(i) for i in range(P)]

layer1_with_id = pd.DataFrame(data=layer1_data,
             index=range(N),
             columns=colnames)

layer2_with_id = pd.DataFrame(data=layer2_data,
             index=range(N),
             columns=colnames)

layer1_with_id.index.name = "id"
layer2_with_id.index.name = "id"

layer1_with_id.to_csv('toy_data_layer1.tsv',sep='\t')
layer2_with_id.to_csv('toy_data_layer2.tsv',sep='\t')
