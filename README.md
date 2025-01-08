# URSA

Release version: 1.0

URSA is an electronic structure code that uses a real-space Finite Element Method (FEM) discretization to perform GW approximation calculations using all-electron ground-state DFT or Hartree-Fock as the starting points. The code is taking advantage of the state-of-the-art FEAST [1,2] as the non-linear eigenvalue solver for the GW quasiparticle equations [3]. URSA is an extension of an 3D FEM legacy code - NESSIE [4]. The real space mesh is generated using the software tetgen [5]. Remark 1: Current atom database is limited to two rows of the periodic table (but it is easy to add more); Only isolated systems (non-periodic) are considered in this release (no bandsructrure calculations). Remark 2: Only graphical solution and Spectral function method of G0W0@DFT and G0W0@HF are included in this current URSA release while G0W0@DFT and G0W0@HF with non-linear eigenvalue solver will be added in the next release.

[1] [E. Polizzi, "Density-Matrix-Based Algorithms for Solving Eigenvalue Problems",in Phys. Rev. B. Vol. 79, 115112, 2009](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.79.115112)

[2] [FEAST Eigenvalue Solver](http://www.feast-solver.org/)

[3] [D. Li, E. Polizzi, "FEAST nonlinear eigenvalue algorithm for GW quasiparticle equations", arXiv:2409.06119, 2024](https://arxiv.org/abs/2409.06119)

[4] [NESSIE](http://www.nessie-code.org/)

[5] [Tetgen](https://wias-berlin.de/software/index.jsp?id=TetGen)



## Installation

__1__- Install the Intel-MKL library.

__2__- Install the [FEAST library](http://www.feast-solver.org/) (BSD License)--see corresponding installation instructions.

__3__- Install the [Tetgen library](https://wias-berlin.de/software/index.jsp?id=TetGen) (GNU-GPL License)-- install the executable __tetgen__ in your PATH or copy it in the /bin URSA directory.

__4__- Installation of [LibXC library](https://libxc.gitlab.io/) is optional. Please make sure it is compiled with the same compiler. Without LibXC, only LDA XC functional of DFT is supported. 

__5__- Deﬁne the Shell variable URSAROOT, e.g. export URSAROOT=<URSA directory> or set URSAROOT=<URSA directory> respectively for the BASH or CSH shells.
One of this command can be placed in the appropriate shell startup ﬁle in $HOME (i.e .bashrc or .cshrc).

__6__- In your shell startup file you also need to include this line (using BASH shell): export PATH="$URSAROOT/bin:$PATH"

__7__- Navigate to the /src directory and execute

```bash
make all
```
By default, the compilation is using ifx and intel MPI. If you wish to use a different compiler or MPI directives, just enter "__make__" to see all the possible options.

__8__- done


## Usage and Example

__1__- Create a standard .xyz atomic file which contains the atoms and their coordinates in Angstrum. The /simulations directory contains multiple sub-directories with their
example xyz files:


__He__ for Helium,  __H2__ for diHydrogen, __CO__ for Carbon Monoxide, __H2O__  for Water, __N2__ for Nitrogen molecule, etc

In the following, let us take the example of __H2__.

__2__- To run G0W0 simulations, we need to create the input file: __H2.gw__
This input file can be automatically generated by executing (for example):
```bash
ursa_configure H2 dft cd p2
```
Here:

-__dft__ stands for DFT simulations (!Options: dft, hf),

-__cd__ stands for the Contour-Deformation calculation method for the GW self-energy (Options: cd, casida, none (no GW)),

-__p2__ stands for P2-FEM real-space mesh discretization (p3 is also available for higher accuracy). We note that running G0W0 with p2 provides good enough results for the HOMO energies.

Before executing the G0W0 simulations, you can edit the file __H2.gw__ and change various options. Many of the these options are self-explanatory (such as p2 vs p3) but many others are expert options. A comprehensive URSA documentation is not yet available. 

__3__- After executing the __ursa_configure__ file, the simulation can be run using the following line
```bash
ursa_compute H2
```


Many information will be displayed on screen during the G0W0 calculations. Here is the screenshot of some results. The last list of (complex) numbers are the correlation energies ($<\psi_{KS}^{HOMO}|\Sigma_c(\omega)|\psi_{KS}^{HOMO}>$) in atomic units calculated by Contour-Deformation method with 20 (N_grid=20) energy grids -- between \[Emin=-17.5 eV, Emax=2.5 eV\].

```bash
 Nn -        4824 Ne -        3459
 # of occupied orbitals           1
  -------------------------------------------
 | compute initial guess of electron density |
 | & diagonalize Kohn-Sham equation          |
  -------------------------------------------
 ***
  ------------------------------------------
 |           SCF-iteration starts           |
  ------------------------------------------
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SCF-Cycle ----           0
 --- Energy states up to occupied + LUMO (eV) ---
           1  -13.5149078487864     
           2 -0.473829703823659     
 SCF-loop residual -------   0.710946383567262     
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SCF-Cycle ----           1
 --- Energy states up to occupied + LUMO (eV) ---
           1  -11.7032590040045     
           2  0.294520615075446     
 SCF-loop residual -------   0.240422824425118     
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SCF-Cycle ----           2
 --- Energy states up to occupied + LUMO (eV) ---
           1  -10.0675793945623     
           2  0.696014651915786     
 SCF-loop residual -------   5.741681215220480E-002
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SCF-Cycle ----           3
 --- Energy states up to occupied + LUMO (eV) ---
           1  -10.4119842346510     
           2  0.568924542027834     
 SCF-loop residual -------   8.864128718577743E-003
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SCF-Cycle ----           4
 --- Energy states up to occupied + LUMO (eV) ---
           1  -10.4223724892283     
           2  0.582736135007934     
 SCF-loop residual -------   2.721093419270812E-003
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SCF-Cycle ----           5
 --- Energy states up to occupied + LUMO (eV) ---
           1  -10.4216693053465     
           2  0.595110999481273     
 SCF-loop residual -------   5.086691923474552E-004
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SCF-Cycle ----           6
 --- Energy states up to occupied + LUMO (eV) ---
           1  -10.4216946789234     
           2  0.595162615603790     
 SCF-loop residual -------   3.760992403227131E-004
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SCF-Cycle ----           7
 --- Energy states up to occupied + LUMO (eV) ---
           1  -10.4204027406174     
           2  0.595533313433651     
 SCF-loop residual -------   9.826690430502450E-005
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SCF-Cycle ----           8
 --- Energy states up to occupied + LUMO (eV) ---
           1  -10.4200786349540     
           2  0.595551701010133     
 SCF-loop residual -------   4.267419464840010E-006
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SCF-Cycle ----           9
 --- Energy states up to occupied + LUMO (eV) ---
           1  -10.4200508289930     
           2  0.595555406533187     
 SCF-loop residual -------   2.982778552984982E-006
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SCF-Cycle ----          10
 --- Energy states up to occupied + LUMO (eV) ---
           1  -10.4200549726954     
           2  0.595553391471003     
 SCF-loop residual -------   1.497352175115428E-006
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SCF-Cycle ----          11
 --- Energy states up to occupied + LUMO (eV) ---
           1  -10.4200595676844     
           2  0.595553259971305     
 SCF-loop residual -------   1.823926367985015E-007
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SCF-Cycle ----          12
 --- Energy states up to occupied + LUMO (eV) ---
           1  -10.4200600773861     
           2  0.595553278859117     
 SCF-loop residual -------   3.078541208248705E-008
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SCF-Cycle ----          13
 --- Energy states up to occupied + LUMO (eV) ---
           1  -10.4200602244323     
           2  0.595553276152077     
 SCF-loop residual -------   2.945577497462028E-009
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 Ground State Convergence Reached!!!
 ==============================================
 -----  Orbital,  DFT_Ekinetic (eV) -----
           1   15.4531746077192     
           2   5.71143440447254     
 -----  Orbital,  DFT_Eext (eV) -----
           1  -49.1971316953183     
           2  -17.9081505880650     
 -----  Orbital,  DFT_Eh (eV) -----
           1   35.2792082059234     
           2   16.2386747820364     
 -----  Orbital,  DFT_Ex (eV) -----
           1  -9.62703427631154     
           2  -3.05776786301122     
 -----  Orbital,  DFT_Ec (eV) -----
           1  -1.42748640418873     
           2 -0.263486593096187     
 -----  Orbital,  DFT_Ex+Ec pbe_g (eV) -----
           1 -0.900790662256304     
           2 -0.125150866184476     
 ==============================================
 --- solve more unoccupied states ---
 --- Energy states up to all solved unoccupied states (eV) ---
           1  -10.4200602244323     
           2  0.595553276149457     
           3  0.985167452736554     
           4   1.32807546060004     
           5   1.38774468798066     
           6   3.29998811274398     
           7   3.36450136229456     
           8   3.50035296198647     
           9   3.71031164019140     
          10   3.76774824910850     
          11   3.97818611709550     
          12   4.66779729548959     
          13   4.72581990609421     
          14   4.96632407222522     
          15   5.07625849939762     
          16   5.51121011729669     
          17   5.61426746940361     
          18   5.86851008576786     
          19   5.90329068421243     
          20   5.98535176200943     
          21   6.16105873985066     
          22   6.46387683357575     
          23   6.74013755613410     
  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!  GW starts  !!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 *** construct bare Coulomb potential -- 1/|r-r'| ***
 *** construct GW Exchange Self-energy -- Sigma_x ***
 ==============================================
 ----  Orbital,  GW <psi|Sigma_x|psi> (eV) ----
           1  -17.7834868694120     
           2 -0.988475335417110     
 ==============================================
 *** compute the energy grid points ***
 *** compute imaginary axis energy integral gauss nodes ***
 *** compute dynamically screened Coulomb interaction - W ***
 (2.656711212987296E-002,0.000000000000000E+000)
 (1.574689817913091E-002,0.000000000000000E+000)
 (7.197646933419861E-003,0.000000000000000E+000)
 (1.186865375138612E-004,0.000000000000000E+000)
 (-5.946502572104874E-003,0.000000000000000E+000)
 (-1.128195609957185E-002,0.000000000000000E+000)
 (-1.609014391825136E-002,0.000000000000000E+000)
 (-2.045063381777931E-002,0.000000000000000E+000)
 (-2.452695506960541E-002,0.000000000000000E+000)
 (-2.835224871320701E-002,0.000000000000000E+000)
 (-3.199462104384110E-002,0.000000000000000E+000)
 (-3.549983892163585E-002,0.000000000000000E+000)
 (-3.890686250391783E-002,0.000000000000000E+000)
 (-4.225053643746908E-002,0.000000000000000E+000)
 (-4.556216099315275E-002,0.000000000000000E+000)
 (-4.887139935935044E-002,0.000000000000000E+000)
 (-5.220801981582089E-002,0.000000000000000E+000)
 (-5.578439578417040E-002,0.000000000000000E+000)
 (-5.916480901781986E-002,0.000000000000000E+000)
 (-6.271383117651150E-002,0.000000000000000E+000)
```

Eigenvectors are all saved in output files. The code is also generating a mesh point .vtk file that can be used for plotting the ground state density and potential.



## License
[BSD](https://opensource.org/licenses/BSD-3-Clause)
The URSA code is distributed under the BSD software license.

Copyright (c) 2022-2024, The Regents of the University of Massachusetts, Amherst.
Dongming Li
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the University nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
