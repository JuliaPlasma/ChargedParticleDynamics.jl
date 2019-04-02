#!/bin/bash

if [ ${TRAVIS_OS_NAME} = "linux" ]; then
    julia --project --color=yes -e 'using Pkg; Pkg.instantiate(); import ChargedParticleDynamics; include(joinpath(dirname(pathof(ChargedParticleDynamics)), "..", "docs", "make.jl"))';
fi
