var documenterSearchIndex = {"docs":
[{"location":"charged_particle_3d/#Charged-Particles-in-3D","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"","category":"section"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"The motion of charged particles in an electromagnetic field (EB) is governed by the Lorentz force,","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"ddotx (t) = fracem big E(x(t)) + dotx (t) times B(x(t)) big ","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"where m and e denote the particle's mass and charge, respectively.","category":"page"},{"location":"charged_particle_3d/#Canonical-Formulation","page":"Charged Particles in 3D","title":"Canonical Formulation","text":"","category":"section"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"The canonical form of the equations can be obtained from the Hamiltonian","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"H (xp) = frac12m (p-A(x))^2 + e phi(x)","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"as","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"beginaligned\ndotx (t) = fracpartial Hpartial p (x(t)p(t)) = frac1m (p(t) - A(x(t))) \ndotp (t) = - fracpartial Hpartial x (x(t)p(t)) = frac1m nabla A(x(t)) cdot (p(t)-A(x(t))) - nabla phi(x(t)) \nendaligned","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"where the fields (EB) are related to the potentials (phi A) by","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"beginaligned\nE (x) = - nabla phi (x)  \nB (x) = nabla times A (x) \nendaligned","category":"page"},{"location":"charged_particle_3d/#Noncanonical-Formulation","page":"Charged Particles in 3D","title":"Noncanonical Formulation","text":"","category":"section"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"The noncanonical form of the equations can be obtained from the phasespace Lagrangian","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"beginaligned\nL (xdotxvdotv) = (e A(x) + mv) cdot dotx - H(xv)  \nH(xv) = fracm2 v^2 + e phi(x)\nendaligned","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"as","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"beginaligned\ndotx (t) = v (t)  \ndotv (t) = fracem big nabla A (x(t)) cdot dotx(t) - dotA (x(t)) - nabla phi(x(t)) big \nendaligned","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"Computing the time derivative of A and using the relation between the potentials (phi A) and the fields (EB), this can be rewritten as","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"beginaligned\ndotx (t) = v (t)  \ndotv (t) = fracem big E(x(t)) + v (t) times B(x(t)) big \nendaligned","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"This constitutes a noncanonical Hamiltonian system of the form","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"dotz (t) = Omega^-T (z(t)) nabla H(z(t)) ","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"with z = (xv) and the symplectic matrix Omega given by","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"Omega = frac1m beginpmatrix\n  mathbb0  mathbb1 \n- mathbb1  e hatB \nendpmatrix ","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"and","category":"page"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"hatB = beginpmatrix\n0  -B_3  B_2 \nB_3  0  - B_1 \n- B_2  B_1  0 \nendpmatrix ","category":"page"},{"location":"charged_particle_3d/#Modules","page":"Charged Particles in 3D","title":"Modules","text":"","category":"section"},{"location":"charged_particle_3d/","page":"Charged Particles in 3D","title":"Charged Particles in 3D","text":"Modules = [ChargedParticleDynamics.ChargedParticle3d.SingularField,\n           ChargedParticleDynamics.ChargedParticle3d.SymmetricField,\n           ChargedParticleDynamics.ChargedParticle3d.ThetaPinchCanonical,\n           ChargedParticleDynamics.ChargedParticle3d.ThetaPinchNoncanonical,\n           ChargedParticleDynamics.ChargedParticle3d.TokamakSmallCartesian,\n           ChargedParticleDynamics.ChargedParticle3d.TokamakSmallCylindrical,\n           ChargedParticleDynamics.ChargedParticle3d.TokamakSmallToroidal,\n           ChargedParticleDynamics.ChargedParticle3d.TokamakSmallNoncanonical,\n           ChargedParticleDynamics.ChargedParticle3d.TokamakIterCylindrical,\n           ChargedParticleDynamics.ChargedParticle3d.SolovevIter,\n           ChargedParticleDynamics.ChargedParticle3d.SolovevIterXpoint]","category":"page"},{"location":"charged_particle_3d/#ChargedParticleDynamics.ChargedParticle3d.SingularField","page":"Charged Particles in 3D","title":"ChargedParticleDynamics.ChargedParticle3d.SingularField","text":"Charged Particle in a singular magnetic field of the form B(xyz) = (x^2 + y^2)^-32 e_z.\n\n\n\n\n\n","category":"module"},{"location":"charged_particle_3d/#ChargedParticleDynamics.ChargedParticle3d.SymmetricField","page":"Charged Particles in 3D","title":"ChargedParticleDynamics.ChargedParticle3d.SymmetricField","text":"Charged Particle in an axisymmetric magnetic field of the form B(xyz) = (1 + x^2 + y^2) e_z.\n\n\n\n\n\n","category":"module"},{"location":"charged_particle_3d/#ChargedParticleDynamics.ChargedParticle3d.ThetaPinchCanonical","page":"Charged Particles in 3D","title":"ChargedParticleDynamics.ChargedParticle3d.ThetaPinchCanonical","text":"Charged Particle in an uniform magnetic field of the form B(xyz) = B_0 e_z.\n\n\n\n\n\n","category":"module"},{"location":"normalization/#Normalization","page":"Normalization","title":"Normalization","text":"","category":"section"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"In order to normalize the equations of the various models implemented in this package, we start at the level of the Lagrangian. In the following, this is explained exemplary for the charged particle Lagrangian, however, it generalizes straightforwardly to other systems like the Pauli particle or the guiding center system.","category":"page"},{"location":"normalization/#Charged-Particle-Lagrangian","page":"Normalization","title":"Charged Particle Lagrangian","text":"","category":"section"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"Consider the phasespace Lagrangian","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"L (q dotq v) = ( mv + A (x) ) cdot dotx  - fracm2 vert v vert^2  - e phi (x) ","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"and, in full generality, introduce the following normalizations:","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"t = hatt t  quad\nx = hatx x  quad\nv = hatv v  quad\nA = hatA A  quad\nphi = hatphi phi  quad\nL = hatL L  ","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"This leads us to","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"L\n= fracLhatL\n= frac1hatL  big( m hatv v + e hatA A big) cdot frachatxhatt dotx - fracm hatvhatL fracvert v vert^22 - frace hatphihatL phi ","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"Let us choose the following normalizations (and note that others are possible and may be more appropriate, depending on the problem at hand):","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"beginaligned\nhatv = sqrtfrachatWm  \nhatA = hatl hatB  \nhatphi = hatl hatE  \nhatL = m hatv^2 = e hatphi = hatW \nendaligned","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"with the particle energy W. The normalized Lagrangian becomes","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"L = bigg( frachatxhatt hatv  v + underbracefrace hatBm_hatomega_c frachatx hatl hatt hatv^2  A bigg) cdot dotx - fracvert v vert^22 - frachatlhatx  phi ","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"where omega_c = e B  m is the gyration frequency. This suggests to set","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"beginaligned\nhatt = omega_c^-1  \nhatx = hatt hatv  \nhatl = hatx \nendaligned","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"such that","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"L = ( v + A ) cdot dotx - fracvert v vert^22 - phi ","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"We thus have obtained the normalized Lagrangian.","category":"page"},{"location":"normalization/#Alternative-Normalization","page":"Normalization","title":"Alternative Normalization","text":"","category":"section"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"Often, especially when simulating an ensemble of particles, it is more appropriate to normalize the velocity to the thermal velocity and to choose different normalizations for v and dotx, specifically","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"beginaligned\nhatv = v_mathrmth  \nhatt = omega_c^-1  \nhatA = hatl hatB  \nhatphi = hatl hatE  \nhatL = m hatv^2 = e hatphi = hatW \nendaligned","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"where W = k_B T now denotes the thermal energy and v_mathrmth = sqrt 2 W  m. Consequently,","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"beginaligned\nhatx = hatv hatt = fracv_mathrmthomega_c = rho_mathrmth \nendaligned","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"and the normalized Lagrangian becomes","category":"page"},{"location":"normalization/","page":"Normalization","title":"Normalization","text":"L = bigg( v + frachatlrho_mathrmth  A bigg) cdot dotx - fracvert v vert^22 - frachatlrho_mathrmth  phi ","category":"page"},{"location":"initialization/#Initialization","page":"Initialization","title":"Initialization","text":"","category":"section"},{"location":"initialization/","page":"Initialization","title":"Initialization","text":"In the following, we will discuss how to initialise the various models (charged particle, Pauli particle, guiding center) from a common set of initial data (particle position, total energy, pitch angle).","category":"page"},{"location":"pauli_particle_3d/#Pauli-Particles-in-3D","page":"Pauli Particles in 3D","title":"Pauli Particles in 3D","text":"","category":"section"},{"location":"pauli_particle_3d/#Modules","page":"Pauli Particles in 3D","title":"Modules","text":"","category":"section"},{"location":"pauli_particle_3d/","page":"Pauli Particles in 3D","title":"Pauli Particles in 3D","text":"Modules = [ChargedParticleDynamics.PauliParticle3d.SymmetricField,\n           ChargedParticleDynamics.PauliParticle3d.ThetaPinchField,\n           ChargedParticleDynamics.PauliParticle3d.TokamakSmallCartesian,\n           ChargedParticleDynamics.PauliParticle3d.TokamakSmallCylindrical,\n           ChargedParticleDynamics.PauliParticle3d.TokamakSmallToroidal,\n           ChargedParticleDynamics.PauliParticle3d.TokamakIterCylindrical,\n           ChargedParticleDynamics.PauliParticle3d.SolovevIter,\n           ChargedParticleDynamics.PauliParticle3d.SolovevIterXpoint]","category":"page"},{"location":"pauli_particle_3d/#ChargedParticleDynamics.PauliParticle3d.SymmetricField","page":"Pauli Particles in 3D","title":"ChargedParticleDynamics.PauliParticle3d.SymmetricField","text":"Charged Particle in an axisymmetric magnetic field of the form B(xyz) = (1 + x^2 + y^2) e_z.\n\n\n\n\n\n","category":"module"},{"location":"pauli_particle_3d/#ChargedParticleDynamics.PauliParticle3d.ThetaPinchField","page":"Pauli Particles in 3D","title":"ChargedParticleDynamics.PauliParticle3d.ThetaPinchField","text":"Charged Particle in an uniform magnetic field of the form B(xyz) = B_0 e_z.\n\n\n\n\n\n","category":"module"},{"location":"#ChargedParticleDynamics.jl","page":"Overview","title":"ChargedParticleDynamics.jl","text":"","category":"section"},{"location":"#Models","page":"Overview","title":"Models","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"Pages = [\"charged_particle_3d.md\",\n         \"guiding_center_4d.md\"\n]","category":"page"},{"location":"#License","page":"Overview","title":"License","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"Copyright (c) Michael Kraus <michael.kraus@ipp.mpg.de>Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.","category":"page"},{"location":"guiding_center_4d/#Guiding-Center-Dynamics-in-4D","page":"Guiding Center Dynamics in 4D","title":"Guiding Center Dynamics in 4D","text":"","category":"section"},{"location":"guiding_center_4d/","page":"Guiding Center Dynamics in 4D","title":"Guiding Center Dynamics in 4D","text":"Guiding centre dynamics is a reduced version of charged particle dynamics, where the motion of the particle in a strong magnetic field B is reduced to the motion of the guiding centre, that is the centre of the gyro motion of the particle about a magnetic field line. The dynamics of the guiding centre can be described in terms of only four coordinates (as compared to six for the full motion of the charged particle), the position of the guiding centre r = (xyz) and the parallel velocity u, where parallel refers to the direction of the magnetic field.","category":"page"},{"location":"guiding_center_4d/#Lagrangian-Formulation","page":"Guiding Center Dynamics in 4D","title":"Lagrangian Formulation","text":"","category":"section"},{"location":"guiding_center_4d/","page":"Guiding Center Dynamics in 4D","title":"Guiding Center Dynamics in 4D","text":"The simplest form of the guiding centre equations can be obtained from the Lagrangian","category":"page"},{"location":"guiding_center_4d/","page":"Guiding Center Dynamics in 4D","title":"Guiding Center Dynamics in 4D","text":"beginaligned\nL = (A (r) + u b (r)) cdot dotr - H(ru)  \nH = tfrac12 u^2 + mu vert B (r) vert \nendaligned","category":"page"},{"location":"guiding_center_4d/","page":"Guiding Center Dynamics in 4D","title":"Guiding Center Dynamics in 4D","text":"where b = B  vert B vert is the unit vector of the magnetic field B = nabla times A with A the magnetic vector potential and mu is the magnetic moment. The Euler-Lagrange equations are computed as","category":"page"},{"location":"guiding_center_4d/","page":"Guiding Center Dynamics in 4D","title":"Guiding Center Dynamics in 4D","text":"beginaligned\nnabla vartheta^T ( r(t) u(t) ) cdot dotr (t) - dotvartheta ( r(t) u(t) ) = nabla H ( r(t) u(t) )  \nb ( r(t) ) cdot dotr (t) = u (t) \nendaligned","category":"page"},{"location":"guiding_center_4d/","page":"Guiding Center Dynamics in 4D","title":"Guiding Center Dynamics in 4D","text":"with vartheta(ru) = A (r) + u  b (r) and the gradient nabla denoting the derivative with respect to r.","category":"page"},{"location":"guiding_center_4d/#Hamiltonian-Formulation","page":"Guiding Center Dynamics in 4D","title":"Hamiltonian Formulation","text":"","category":"section"},{"location":"guiding_center_4d/","page":"Guiding Center Dynamics in 4D","title":"Guiding Center Dynamics in 4D","text":"Computing the time derivative of vartheta, the Euler-Lagrange equations can be rewritten in an explicit form as","category":"page"},{"location":"guiding_center_4d/","page":"Guiding Center Dynamics in 4D","title":"Guiding Center Dynamics in 4D","text":"beginaligned\ndotr (t) = dfracu (t)  beta (r (t))b (r (t)) cdot beta (r (t)) + dfracB (r (t))B (r (t)) cdot beta (r (t)) times nabla H (r (t)u (t))  \ndotu (t) = - dfracbeta (r (t))b (r (t)) cdot beta (r (t)) cdot nabla H (r (t) u (t)) \nendaligned","category":"page"},{"location":"guiding_center_4d/","page":"Guiding Center Dynamics in 4D","title":"Guiding Center Dynamics in 4D","text":"where beta = nabla times vartheta. This constitutes a noncanonical Hamiltonian system of the form","category":"page"},{"location":"guiding_center_4d/","page":"Guiding Center Dynamics in 4D","title":"Guiding Center Dynamics in 4D","text":"dotq = Omega^-T (q) nabla H(q) ","category":"page"},{"location":"guiding_center_4d/","page":"Guiding Center Dynamics in 4D","title":"Guiding Center Dynamics in 4D","text":"with q = (xyzu) and the symplectic matrix Omega given by","category":"page"},{"location":"guiding_center_4d/","page":"Guiding Center Dynamics in 4D","title":"Guiding Center Dynamics in 4D","text":"Omega_ij = dfracpartial vartheta_jpartial q^i - dfracpartial vartheta_ipartial q^j ","category":"page"},{"location":"guiding_center_4d/#Modules","page":"Guiding Center Dynamics in 4D","title":"Modules","text":"","category":"section"},{"location":"guiding_center_4d/","page":"Guiding Center Dynamics in 4D","title":"Guiding Center Dynamics in 4D","text":"Modules = [ChargedParticleDynamics.GuidingCenter4d.SymmetricField,\n           ChargedParticleDynamics.GuidingCenter4d.ThetaPinchField,\n           ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCartesian,\n           ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCylindrical,\n           ChargedParticleDynamics.GuidingCenter4d.TokamakSmallToroidal,\n           ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCartesian,\n           ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCylindrical,\n           ChargedParticleDynamics.GuidingCenter4d.TokamakIterCylindrical,\n           ChargedParticleDynamics.GuidingCenter4d.SolovevIter,\n           ChargedParticleDynamics.GuidingCenter4d.SolovevIterXpoint]","category":"page"},{"location":"guiding_center_4d/#ChargedParticleDynamics.GuidingCenter4d.SymmetricField","page":"Guiding Center Dynamics in 4D","title":"ChargedParticleDynamics.GuidingCenter4d.SymmetricField","text":"First and second Poincaré invariant for a guiding center particle in an axisymmetric magnetic field of the form B(xyz) = B_0 (1 + x^2 + y^2) e_z.\n\nThe loop for the first invariant is initialized by\n\nq (tau) = beginpmatrix\nr_x cos (2pi tau) \nr_y sin (2pi tau) \nz_0 + z_1 sin (2pi tau) \nu_0 + u_1 cos (2pi tau) \nendpmatrix\n\nwith parameters\n\nB_0 = 1 quad\nr_x = 05 quad\nr_y = 03 quad\nz_0 = 00 quad\nz_1 = 01 quad\nu_0 = 05 quad\nu_1 = 005 quad\nmu = 001 \n\nThe surface for the second invariant is initialized by\n\nq (tau) = beginpmatrix\nr_0 (sigma - 05) \nr_0 (tau   - 05) \nz_0 + z_1 cos (2pi sigma) cos (2pi tau) \nu_0 + u_1 sin (2pi sigma) sin (2pi tau) \nendpmatrix\n\nwith parameters\n\nB_0 = 1 quad\nr_0 = 05 quad\nz_0 = 00 quad\nz_1 = 01 quad\nu_0 = 05 quad\nu_1 = 001 quad\nmu = 001 \n\n\n\n\n\n","category":"module"},{"location":"guiding_center_4d/#ChargedParticleDynamics.GuidingCenter4d.ThetaPinchField","page":"Guiding Center Dynamics in 4D","title":"ChargedParticleDynamics.GuidingCenter4d.ThetaPinchField","text":"First Poincaré invariant for a guiding center particle in an θ-pinch magnetic field of the form B(xyz) = B_0  e_z.\n\nThe loop for the first Poincaré invariant is initialized by\n\nq (tau) = beginpmatrix\nr_x cos (2pi tau) \ny_0 \nr_z sin (2pi tau) \nendpmatrix\n\nwith parameters\n\nB_0 = 1 quad\nr_x = 05 quad\nr_z = 03 quad\ny_0 = 00 quad\nu_0 = 05 quad\nmu = 25 times 10^-6 \n\n\n\n\n\n","category":"module"},{"location":"guiding_center_4d/#ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCartesian","page":"Guiding Center Dynamics in 4D","title":"ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCartesian","text":"Analytic axisymmetric small tokamak equilibrium in cartesian coordinates.\n\n\n\n\n\n","category":"module"},{"location":"guiding_center_4d/#ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCylindrical","page":"Guiding Center Dynamics in 4D","title":"ChargedParticleDynamics.GuidingCenter4d.TokamakSmallCylindrical","text":"Analytic axisymmetric small tokamak equilibrium in cylindrical coordinates.\n\n\n\n\n\n","category":"module"},{"location":"guiding_center_4d/#ChargedParticleDynamics.GuidingCenter4d.TokamakSmallToroidal","page":"Guiding Center Dynamics in 4D","title":"ChargedParticleDynamics.GuidingCenter4d.TokamakSmallToroidal","text":"Analytic axisymmetric small tokamak equilibrium in circular coordinates.\n\n\n\n\n\n","category":"module"},{"location":"guiding_center_4d/#ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCartesian","page":"Guiding Center Dynamics in 4D","title":"ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCartesian","text":"Analytic axisymmetric medium-size tokamak equilibrium in cartesian coordinates.\n\n\n\n\n\n","category":"module"},{"location":"guiding_center_4d/#ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCylindrical","page":"Guiding Center Dynamics in 4D","title":"ChargedParticleDynamics.GuidingCenter4d.TokamakMediumCylindrical","text":"Analytic axisymmetric medium-size tokamak equilibrium in cartesian coordinates.\n\n\n\n\n\n","category":"module"},{"location":"guiding_center_4d/#ChargedParticleDynamics.GuidingCenter4d.TokamakIterCylindrical","page":"Guiding Center Dynamics in 4D","title":"ChargedParticleDynamics.GuidingCenter4d.TokamakIterCylindrical","text":"Analytic ITER-like Solov'ev equilibrium with X-point.\n\n\n\n\n\n","category":"module"},{"location":"guiding_center_4d/#ChargedParticleDynamics.GuidingCenter4d.SolovevIter","page":"Guiding Center Dynamics in 4D","title":"ChargedParticleDynamics.GuidingCenter4d.SolovevIter","text":"Analytic ITER-like Solov'ev equilibrium.\n\n\n\n\n\n","category":"module"},{"location":"guiding_center_4d/#ChargedParticleDynamics.GuidingCenter4d.SolovevIterXpoint","page":"Guiding Center Dynamics in 4D","title":"ChargedParticleDynamics.GuidingCenter4d.SolovevIterXpoint","text":"Analytic ITER-like Solov'ev equilibrium with X-point.\n\n\n\n\n\n","category":"module"}]
}
