#!/bin/sh

#Run three scripts in sequence, but clear memory between them
Rscript /Integratomics/R_integration_code_dockerVersion.r "$@"
Rscript -e "gc()"

Rscript /Integratomics/Integration_bootstrapping_dockerVersion.r "$@"
Rscript -e "gc()"

Rscript /Integratomics/Integration_bootstrap_analysis_dockerVersion.r "$@"