# Refinement of Interval Discrete Mixed-Monotonic Equations (RIDMME)

Author: Alexandra Catano-Lopez 
- [Orcid](https://orcid.org/0000-0003-4370-3382)
- [Google Scholar](https://scholar.google.com/citations?user=MSjCw2MAAAAJ&hl=es)

## Description

The RIDMME is an algorithm created for solving the forward problem using Interval Discrete Equations and handling initial conditions and parameters with uncertainty. This library includes reachability analysis using the following methods:

1. The Mixed-Monotonic method proposed and developed by Meyer et al [1]. You can find the (TIRA Toolbox) [https://www.mathworks.com/matlabcentral/fileexchange/105115-tira-toolbox-for-interval-reachability-analysis] for this method.
2. Improvement of the Mixed-Monotonic method through adaptive subdivision of interval-valued initial conditions and parameters by implementing automatic differentiation (RIDMME algorithm).
3. Montecarlo simulation.

For the interval-valued arithmetic implemented in the automatic differentiation process, we use the (PROFIL/BIAS library) [https://www.tuhh.de/ti3/keil/profil/index_e.html], considered the fastest library for interval-valued arithmetic.


## Installation

To install the PROFIL/BIAS library, use the following command:

> ./install_bias.sh

We select the configuration number 11 (_x86-64-Linux-compat-gcc_), making only avaliable the for **Linux** distribution.

Then, create the build folder
> mkdir build

and compile it
> cd build

> cmake .. .

> cmake --build . --config Release | Debug

## Citation 

Cite this library as:
> Catano-Lopez, A. 2023. Refinement of Interval Discrete Mixed-Monotonic Equations (RIDMME). Available at: https://github.com/alexacl95/RIDMME 

```
@software{Acatano2023,
  author = {Catano-Lopez, Alexandra},
  title = {{Refinement of Interval Discrete Mixed-Monotonic Equations (RIDMME)}},
  url = {https://github.com/alexacl95/RIDMME},
  version = {1.0},
  year = {2023}
}
```

## References

[1] Meyer, P.J., Devonport, A., and Arcak, M. 2021. Interval Reachability Analysis. In SpringerBriefs in Electrical and Computer Engineering. Springer International Publishing. https://doi.org/10.1007/978-3-030-65110-7
