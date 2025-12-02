---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.17.2
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# RIXS calculations for an Anderson Impurity Model
An Anderson Impurity Model (AIM) refers to a Hamiltonian with a set of correlated orbitals, often referred to as the impurity or metal states, that hybridize with a set of uncorrelated orbitals, commonly known as the ligands or bath states. 

This example examines RIXS simulations for NiPS₃ based on Ref. [^1].

![AIM](./levels.png)

In this example, we will focus on some properties of the model. If desired, the function we will run can be examined [here](https://github.com/EDRIXS/MAIQMag_workshop/blob/main/tutorials/AIM/helper_function.py).


```{code-cell} ipython3
:tags: [remove_output]

import numpy as np
import matplotlib.pyplot as plt
from helper_function import make_rixs

%matplotlib inline
```

```{code-cell} ipython3
:tags: [remove_input, remove_output]

# Secretly make figures look nice!
plt.rcParams.update({'figure.dpi': 150, 'savefig.dpi': 150,
                     'font.size': 8})
```

## Original model

```{code-cell} ipython3
---
tags: [remove_output]
jupyter:
  outputs_hidden: true
---
ominc, eloss, rixs, impurity_occupation = make_rixs()
```

```{code-cell} ipython3
:tags: [hide-output]

fig, ax = plt.subplots(figsize=(3, 4))
art = ax.pcolormesh(ominc, eloss, rixs.sum(-1).T, cmap='terrain', shading='gouraud',
                   vmin=0, vmax=1.5)
plt.colorbar(art, ax=ax, orientation='horizontal')
ax.set_xlabel('Incident energy (eV)')
ax.set_ylabel('Energy loss (eV)')
ax.set_title(f"Model from Nat. Comm. 15, 3496 (2024)")
plt.tight_layout()
plt.show()
```
Let's discuss how we obtained this model. (If desired, the full gory details are in Ref. [^2].)

## Magnetic field behavior
What is the spin of the ground state? And what will happen when a magnetic field is applied?

```{code-cell} ipython3
:tags: [remove_output]

ominc, eloss, rixs, impurity_occupation = make_rixs(ext_B=np.array([0.1, 0, 0]))
```

```{code-cell} ipython3
:tags: [hide-output]

fig, ax = plt.subplots(figsize=(3, 4))
art = ax.pcolormesh(ominc, eloss, rixs.sum(-1).T, cmap='terrain', shading='gouraud',
                   vmin=0, vmax=0.5)
plt.colorbar(art, ax=ax, orientation='horizontal')
ax.set_xlabel('Incident energy (eV)')
ax.set_ylabel('Energy loss (eV)')
ax.set_title(f"{impurity_occupation:.1f} electrons on Ni")
ax.set_ylim(-0.1, 0.5)
plt.tight_layout()
plt.show()
```

## Charge-transfer energy behavior
What happens to the charge on Ni and the spin-flip transition energy when the charge transfer energy is made very large?

```{code-cell} ipython3
---
tags: [remove_output]
jupyter:
  outputs_hidden: true
---
ominc, eloss, rixs, impurity_occupation = make_rixs(ext_B=np.array([0.1, 0, 0]),
                                                    Delta=10, c_level=-696)
```

```{code-cell} ipython3
:tags: [hide-output]

fig, ax = plt.subplots(figsize=(3, 4))
art = ax.pcolormesh(ominc, eloss, rixs.sum(-1).T, cmap='terrain', shading='gouraud',
                   vmin=0, vmax=0.5)
plt.colorbar(art, ax=ax, orientation='horizontal')
ax.set_xlabel('Incident energy (eV)')
ax.set_ylabel('Energy loss (eV)')
ax.set_title(f"{impurity_occupation:.1f} electrons on Ni")
ax.set_ylim(-0.1, 0.5)
plt.tight_layout()
plt.show()
```

## Quality of match
This plot, taken from our recent paper [^1], illustrates the level of precision in a good theory-experiment match.
![Match](./match.png)


## Challenge
Parameter extraction is hard because (a) the experiment theory matche is imperfect (b) degeneracies and the solution including an anticorrelation between charge transfer and Hund's exchange. Several different solutions exist in the literature.[^1],[^3],[^4],[^5]

[^1]: Wei He et al.,
       [Nature Communications 15, 3496 (2024)](https://doi.org/10.1038/s41467-024-47852-x).

[^2]: Wei He et al.,
       Nature Communications 15, 3496 (2024) [Supplementary information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-47852-x/MediaObjects/41467_2024_47852_MOESM1_ESM.pdf)

[^3]: S. Kang et al.,
      [Nature 583, 785–789 (2020)](https://doi.org/10.1038/s41586-020-2520-5)

[^4]: Y. Cao et al., [arXiv:2407.01881](https://arxiv.org/abs/2407.01881)

[^5]: L. Schue et al., [arXiv:2506.13381](https://arxiv.org/abs/2506.13381)
