---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.18.1
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Algorithms

We need to compute

$$
I \propto \frac{1}{\mathcal{Z}(T)}\sum_{i}e^{- E_{i}/(k_\mathrm{B}T)} \\
\times \sum_f |M_{fi}|^2 \delta(E_f + \hbar\omega_{\boldsymbol{k}^\prime} - E_i - \hbar\omega_{\boldsymbol{k}} )
$$

$$
M_{fi} = \sum_n \frac{\bra{f} {\cal D}^\dagger_{\boldsymbol{k}^\prime\hat{\epsilon}^\prime}\ket{n}\bra{n} {\cal D}^{\phantom\dagger}_{\boldsymbol{k}\hat{\epsilon}}\ket{i}}{E_n - E_i - \hbar\omega_{\boldsymbol{k}}+\mathrm{i}\Gamma_n}.
$$

## Direct solution
EDRIXS hosts a pure Python implementation of the cross-section that is suitable for performance insensitive cases and testing. This was what we used in the [atomic model](../atomic/atomic_model.md) section.
See the diagonalization step in [edrixs.solvers.ed_1v1c_py](https://github.com/EDRIXS/edrixs/blob/master/edrixs/solvers.py#L28) and spectrum construction in [edrixs.solvers.rixs_1v1c_py](https://github.com/EDRIXS/edrixs/blob/master/edrixs/solvers.py#L419). The following is a _crude_ illustration of what happens

* Build `emat` and `umat` in $Y_l^m$ single-particle basis

```{code-cell} ipython3
:tags: [remove_output]

import edrixs
import scipy
import numpy as np
```

```{code-cell} ipython3
ten_dq = 1.0
emat_i = np.zeros((14, 14), dtype=complex)
emat_i[:10, :10] = edrixs.cf_cubic_d(ten_dq) # plus SOC and hopping

F0, F2, F4 = 4.0, 1.0, 0.3
umat_i = np.zeros((14, 14, 14, 14), dtype=complex)
umat_i = edrixs.get_umat_slater('d', F0, F2, F4) # plus core hole terms 

# For L3-edge RIXS of d-shell material
print(f"emat is 14x14")
print(f"umat is 14x14x14x14")
```

* Construct [Fock basis](https://github.com/EDRIXS/edrixs/blob/master/edrixs/fock_basis.py#L54) for $\ket{i}$ and $\ket{n}$.

```{code-cell} ipython3
v_norb = 10
v_noccu = 8
c_norb = 4
basis_i = edrixs.get_fock_bin_by_N(v_norb, v_noccu, c_norb, c_norb)
basis_n = edrixs.get_fock_bin_by_N(v_norb, v_noccu+1, c_norb, c_norb - 1)
print(np.array(basis_i)[:3, :10])
```

* Build Hamiltonian for initial and intermediate state using [edrixs.two_fermion](https://github.com/EDRIXS/edrixs/blob/master/edrixs/manybody_operator.py#L109) and [edrixs.four_fermion](https://github.com/EDRIXS/edrixs/blob/master/edrixs/manybody_operator.py#L52).

```{code-cell} ipython3
ncfg_i, ncfg_n = len(basis_i), len(basis_n)
print(f"{ncfg_i} ground state vectors")
print(f"{ncfg_n} intermediate state vectors")
hmat_i = np.zeros((ncfg_i, ncfg_i), dtype=complex)
hmat_n = np.zeros((ncfg_n, ncfg_n), dtype=complex)
hmat_i[:, :] += edrixs.two_fermion(emat_i, basis_i, basis_i)
hmat_i[:, :] += edrixs.four_fermion(umat_i, basis_i)
# hmat_n left out for brevity. 
```

* Diagonalize Hamiltonian for initial and intermediate state using scipy, which calls LAPACK under the hood.

```{code-cell} ipython3
eval_i, evec_i = scipy.linalg.eigh(hmat_i)
eval_n, evec_n = scipy.linalg.eigh(hmat_n)
```

* Load relevant dipole matrix operators

```{code-cell} ipython3
trans_op_Ylm = edrixs.get_trans_oper('dp32')
```

* Build the many-body transition operators for each polarization and express them in the eigenbasis

```{code-cell} ipython3
npol = 3
trans_op = np.zeros((npol, ncfg_n, ncfg_i), dtype=complex)
for i in range(npol):
    trans_op[i] = edrixs.cb_op2(edrixs.two_fermion(trans_op_Ylm[i], basis_n, basis_i),
                                evec_n, evec_i)
print(trans_op.shape)
```

* RIXS spectra can be built from `trans_op`, `eval_i`, and `eval_n`.

## Current EDRIXS solver
Currently, EDRIXS calls into its [Fortran layer](https://github.com/EDRIXS/edrixs/tree/master/src). This is sparsely documented and suffers from serious limitations in portability, maintainability, and flexibility. See [edrixs.ed_siam_fort](https://github.com/EDRIXS/edrixs/blob/master/edrixs/solvers.py#L1819) and [edrixs.rixs_siam_fort](https://github.com/EDRIXS/edrixs/blob/master/edrixs/solvers.py#L2478), which calls [ed_driver.f90](https://github.com/EDRIXS/edrixs/blob/master/src/ed_driver.f90) and [rixs_driver.f90](https://github.com/EDRIXS/edrixs/blob/master/src/rixs_driver.f90). We used this in the [Anderson impurity model section](../AIM/AIM.md).

It is useful to re-express the cross-section as

$$
M_{fi} = \bra{f} {\cal D}^\dagger_{\boldsymbol{k}^\prime,\hat{\epsilon}^\prime} \frac{1}{{\cal H}_n - E_i - \hbar\omega_{\boldsymbol{k}}+i\Gamma_n} {\cal D}^{\phantom\dagger}_{\boldsymbol{k},\hat{\epsilon}}\ket{i}
$$

The general process is:

* Fock basis is generated in integer representation

```{code-cell} ipython3
for binary_string in ["011", "101", "110"]:
    print(f"{binary_string} is {int(binary_string, 2)}")
```

* Hamiltonian is constructed in CSR format. For example, the matrix

$$
\begin{bmatrix}
10 & 20 &  0 &  0 &  0 &  0 \\
 0 & 30 &  0 & 40 &  0 &  0 \\
 0 &  0 & 50 & 60 & 70 &  0 \\
 0 &  0 &  0 &  0 &  0 & 80
\end{bmatrix}
$$

is encoded as

```{code-cell} ipython3
V         = [10, 20, 30, 40, 50, 60, 70, 80]
COL_INDEX = [0,  1,  1,  3,  2,  3,  4,  5]
ROW_INDEX = [0,  2,  4,  7,  8]
```

* A handful of the eigenvectors of the ground state Hamiltonian ${\cal H}$ are obtained via the Lanczos method within the parallel ARPACK library. 

$$
{\cal H} \ket{i} = E_i \ket{i}
$$

* Apply the absorption transition operator to the ground state

$$
\ket{b} = {\cal D}_i \ket{i}.
$$

* Solve the following linear equation, involving the intermediate state Hamiltontian ${\cal H}^\prime$ via sparse [MINRES](https://en.wikipedia.org/wiki/Minimal_residual_method) methods 

$$
(\frac{1}{{\cal H}^\prime - E_i - \hbar\omega_{\boldsymbol{k}}+i\Gamma_n}) \ket{x} = \ket{b}
$$

* Apply emission operator

$$
\ket{F} = {\cal D}_f^\dagger \ket{x}
$$

* The spectrum can then be represented as

$$
I \propto \sum_{i}e^{- E_{i}/(k_\mathrm{B}T)} \Im \bra{F} \frac{1}{{\cal H} - E_i - \hbar\omega_{\boldsymbol{k}}+i \Gamma_f} \ket{F}
$$

* The continued fraction technique is then used to construct the spectrum


## Plans
We are still finalizing the planned approach. 

### Algorithm changes

* Improved basis lookup methods for constructing Hamiltonian based on Ref. [^1].

* Implement methods to re-use the shared structure of Krylov subspaces for different incident x-ray energies based on Ref. [^2]. Kipton recommends this preprint [^3].


### Infrastructure changes
The current plan is to call into more modern platform-agnostic CPU/GPU packages 
[PETSc](https://petsc.org/release/)+[SLEPc](https://slepc.upv.es/). These, which are primarily written in C, have Python interfaces.

## References

[^1]: C.J. Jia, Y. Wang, C.B. Mendl, B. Moritz, T.P. Devereaux, Paradeisos: A perfect hashing algorithm for many-body eigenvalue problems
       [Computer Physics Communications
 224, 81-89 (2017)](https://doi.org/10.1016/j.cpc.2017.11.011)

[^2]: Prakash Sharma, Luogen Xu, Fei Xue, and Yao Wang, Paradeisos: Accelerating resonant spectroscopy simulations using multishifted biconjugate gradient
       [Phys. Rev. B 112, 115113 (2025)](https://doi.org/10.1103/gfdn-pyr2)

[^3]: Tyler Chen, The Lanczos algorithm for matrix functions: a handbook for scientists [arXiv:2410.11090](https://arxiv.org/abs/2410.11090)
