"""Regenerate the DRAGON parity fixtures.

Run this from the netZooPy repo root after intentional changes to DRAGON's
covariance / precision / partial-correlation pipeline. The same files must
then be copied into netZooR's tests/testthat/dragon_parity/ to keep the
two repos in sync.

    python tests/dragon_parity/regenerate.py
"""
import os
import numpy as np
from netZooPy import dragon
from netZooPy.dragon.dragon import Scale

OUT = os.path.dirname(os.path.abspath(__file__))

n, p1, p2 = 50, 20, 15
SEED = 20250510

X1_raw, X2_raw, _, _ = dragon.simulate_dragon_data(
    eta11=0.05, eta12=0.05, eta22=0.05,
    p1=p1, p2=p2, epsilon=[0.1, 0.1], n=n, seed=SEED,
)
X1 = Scale(X1_raw)
X2 = Scale(X2_raw)

lambdas, _ = dragon.estimate_penalty_parameters_dragon(X1, X2)
cov     = dragon.get_shrunken_covariance_dragon(X1, X2, lambdas)
prec, _ = dragon.get_precision_matrix_dragon(X1, X2, lambdas)
ggm     = dragon.get_partial_correlation_dragon(X1, X2, lambdas)

np.savetxt(os.path.join(OUT, "X1.csv"),   X1,   delimiter=",", fmt="%.17g")
np.savetxt(os.path.join(OUT, "X2.csv"),   X2,   delimiter=",", fmt="%.17g")
np.savetxt(os.path.join(OUT, "cov.csv"),  cov,  delimiter=",", fmt="%.17g")
np.savetxt(os.path.join(OUT, "prec.csv"), prec, delimiter=",", fmt="%.17g")
np.savetxt(os.path.join(OUT, "ggm.csv"),  ggm,  delimiter=",", fmt="%.17g")
with open(os.path.join(OUT, "lambdas.txt"), "w") as f:
    f.write("%.17g\n%.17g\n" % (lambdas[0], lambdas[1]))

print("seed:", SEED, "n,p1,p2:", n, p1, p2, "lambdas:", lambdas)
print("wrote:", OUT)
