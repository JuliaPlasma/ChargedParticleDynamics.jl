
# ChargedParticleDynamics.jl

![CI](https://github.com/JuliaPlasma/ChargedParticleDynamics.jl/workflows/CI/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/JuliaPlasma/ChargedParticleDynamics.jl/badge.svg)](https://coveralls.io/github/JuliaPlasma/ChargedParticleDynamics.jl)
[![codecov](https://codecov.io/gh/JuliaPlasma/ChargedParticleDynamics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaPlasma/ChargedParticleDynamics.jl)


## General

```@contents
Pages = ["normalization.md",
         "initialization.md",
]
```


## Models

```@contents
Pages = ["charged_particle_3d.md",
         "guiding_center_3d.md",
         "guiding_center_4d.md",
         "pauli_particle_3d.md",
         "audit.md",
]
```


## References

- Xinjie Li, Ruili Zhang, Jian Liu, *Approximately symplectic Runge-Kutta methods for the guiding
  center dynamics*, Journal of Computational Physics **563**, 115069 (2026),
  [doi:10.1016/j.jcp.2026.115069](https://doi.org/10.1016/j.jcp.2026.115069).
- Ruili Zhang, Jian Liu, Tong Liu, Wenxiang Li, Xiaogang Wang, *Canonical Hamiltonian guiding
  center theory and classical intrinsic magnetic moment*, Frontiers of Physics **21**(2), 026200
  (2026).
- Jianyuan Xiao, Hong Qin, *Slow manifolds of classical Pauli particle enable structure-preserving
  geometric algorithms for guiding center dynamics*,
  [arXiv:2006.03818](https://arxiv.org/abs/2006.03818) (2020).
- Joshua W. Burby, *Guiding center dynamics as motion on a formal slow manifold in loop space*,
  Journal of Mathematical Physics **61**, 012703 (2020),
  [doi:10.1063/1.5119801](https://doi.org/10.1063/1.5119801).


## Acknowledgements

Stefan Possaner contributed vitally to sorting out the appalling details of normalization and initialization.


## License

> Copyright (c) Michael Kraus <michael.kraus@ipp.mpg.de>
>
> Permission is hereby granted, free of charge, to any person obtaining a copy
> of this software and associated documentation files (the "Software"), to deal
> in the Software without restriction, including without limitation the rights
> to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
> copies of the Software, and to permit persons to whom the Software is
> furnished to do so, subject to the following conditions:
>
> The above copyright notice and this permission notice shall be included in all
> copies or substantial portions of the Software.
>
> THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
> IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
> FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
> AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
> LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
> OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
> SOFTWARE.
