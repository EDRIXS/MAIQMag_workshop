import edrixs
import numpy as np
from mpi4py import MPI
comm = MPI.COMM_WORLD

F0_dd = 7.88
F2_dd = 10.68
F4_dd = 6.68

F0_dp = 7.45
F2_dp = 6.56
G1_dp = 4.92
G3_dp = 2.80

Delta = 1.15 # Controls d-ligand splitting (many different definitions in literature)
ten_dq = 0.42
ten_dq_bath = 2.0

Veg = 1.65
Vt2g = 0.95 

ext_B = np.array([0, 0, 0])

c_level = -678

def make_rixs(F0_dd=F0_dd,
              F2_dd=F2_dd,
              F4_dd=F4_dd,
              F0_dp=F0_dp,
              F2_dp=F2_dp,
              G1_dp=G1_dp,
              G3_dp=G3_dp,
              Delta=Delta,
              ten_dq=ten_dq,
              ten_dq_bath=ten_dq_bath,
              Veg=Veg,
              Vt2g=Vt2g,
              ext_B=ext_B,
              c_level=c_level):
    """
    Wrapper function for RIXS calculation.

    See:
    edrixs.solvers.ed_siam_fort
    and
    edrixs.solvers.rixs_siam_fort

    for parameter defintions and
    https://edrixs.github.io/edrixs/auto_examples/example_3_AIM_XAS.html#sphx-glr-auto-examples-example-3-aim-xas-py
    for explanation.
    """
    
    ## Electrons and active shells
    nd = 8
    norb_d = 10
    norb_bath = 10
    nbath = 1
    v_noccu  = nd + nbath*norb_d
    shell_name = ('d', 'p') 
    
    ## Coulomb interactions
    slater = ([F0_dd, F2_dd, F4_dd],  # initial
              [F0_dd, F2_dd, F4_dd, F0_dp, F2_dp, G1_dp, G3_dp])  # with core hole
    
    ## Energies for different shells
    U_dd = F0_dd - edrixs.get_F0('d', F2_dd, F4_dd)
    U_dp = F0_dp - edrixs.get_F0('dp', G1_dp, G3_dp)
    E_d, E_L = edrixs.CT_imp_bath(U_dd, Delta, nd)
    E_dc, E_Lc, E_p = edrixs.CT_imp_bath_core_hole(U_dd, U_dp, Delta, nd)
    
    
    ## d-electron crystal field
    CF = np.zeros((norb_d, norb_d), dtype=complex)
    diagonal_indices = np.arange(norb_d)
    
    orbital_energies = np.array([e for orbital_energy in
                                 [+0.6 * ten_dq, # dz2
                                  -0.4 * ten_dq, # dzx
                                  -0.4 * ten_dq, # dzy
                                  +0.6 * ten_dq, # dx2-y2
                                  -0.4 * ten_dq] # dxy)
                                 for e in [orbital_energy]*2])
    CF[diagonal_indices, diagonal_indices] = orbital_energies
    
    ## Spin orbit coupling
    zeta_d_i = 0.083
    zeta_d_n = 0.102
    c_soc = 11.4
    trans_c2n = edrixs.tmat_c2r('d',True)
    soc_i = edrixs.cb_op(edrixs.atom_hsoc('d', zeta_d_i), trans_c2n)
    soc_n = edrixs.cb_op(edrixs.atom_hsoc('d', zeta_d_n), trans_c2n)
    
    ## Assemble ``emat`` for Ni
    imp_mat = CF + soc_i + E_d*np.eye(norb_d)
    imp_mat_n = CF + soc_n + E_dc*np.eye(norb_d)
    
    ## Ligand crystal field and hopping
    bath_CF = np.full((nbath, norb_d), 0, dtype=complex)
    bath_CF[0, :2] += ten_dq_bath*.6  # 3z2-r2
    bath_CF[0, 2:6] -= ten_dq_bath*.4  # zx/yz
    bath_CF[0, 6:8] += ten_dq_bath*.6  # x2-y2
    bath_CF[0, 8:] -= ten_dq_bath*.4  # xy
    
    bath_level = bath_CF + np.full((nbath, norb_d), E_L)
    bath_level_n = bath_CF + np.full((nbath, norb_d), E_Lc)
    
    hyb = np.zeros((nbath, norb_d), dtype=complex)
    hyb[0, :2] = Veg  # 3z2-r2
    hyb[0, 2:6] = Vt2g  # zx/yz
    hyb[0, 6:8] = Veg  # x2-y2
    hyb[0, 8:] = Vt2g  # xy
    
    eval_i, denmat, noccu_gs = edrixs.ed_siam_fort(
        comm, shell_name, nbath, siam_type=0, imp_mat=imp_mat, imp_mat_n=imp_mat_n,
        bath_level=bath_level, bath_level_n=bath_level_n, hyb=hyb, c_level=c_level,
        c_soc=c_soc, slater=slater,  trans_c2n=trans_c2n, v_noccu=v_noccu,
        do_ed=1, ed_solver=2, neval=50, nvector=3, ncv=100, idump=True,
        ext_B=ext_B, on_which='spin')  

    impurity_occupation = np.sum(denmat[0].diagonal()[0:norb_d]).real

    ominc = np.arange(850.5, 856, .25)
    eloss = np.arange(-0.25, 2, 0.01)

    gamma_f = 0.025
    gamma_c = 0.6

    temperature = 40
    thin = np.deg2rad(22.6)
    thout = np.deg2rad(150 - 22.6)
    pol_type = [('linear', 0, 'linear', 0), ('linear', 0, 'linear', np.pi/2)]

    rixs, poles = edrixs.rixs_siam_fort(
        comm,
        shell_name,
        nbath,
        ominc,
        eloss,
        gamma_c=gamma_c,
        gamma_f=gamma_f,
        v_noccu=v_noccu,
        thin=thin,
        thout=thout,
        pol_type=pol_type,
        num_gs=3,
        temperature=temperature
    )
    return ominc, eloss, rixs, impurity_occupation
