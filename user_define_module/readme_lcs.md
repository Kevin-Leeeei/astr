# pre/post-process files  written by LCS
Updated by C.S. Luo @ Beihang University, 1/11/2024

All pp modules allow using parallel computing. We use [FFTW](www.fftw.org/) to realize FFT, with self-module [fftwlink](../src/fftwlink.F90) used for special parallel allocations.

## Using method
`mpirun  -np $NP $ASTR_DIR/astr pp $MODULE $SUBMODULE $PARAMETER $INPUTFILE`

example: `mpirun  -np 16 $ASTR_DIR/astr pp hitgen 2Dherring ./datin/input.2d`

## Modules
The user-defined-pp is separated for 4 principal modules.
### Modules and submodules
1. [`spectra`](./udf_pp_spectra.F90): spectrum-related calculation
    - `instant3D` or `instant2D` : calculation of velocity spectrum $E^s(k)$ and $E^c(k)$ for solenoidal(incompressible) and dilational(compressible) kinetic energy of one snapshot. `$PARAMETER1`: output file step. `$PARAMETER2`: spectral calculation method, basing on [[1]](#ref1), method 2 is recommended. (Method 1/2: $\langle E(\boldsymbol{k}) \rangle_{|\boldsymbol{k}| \in S_k} 2 \pi k$ or $\langle E(\boldsymbol{k}) \rangle_{|\boldsymbol{k}| \in S_k} 4 \pi k^2$,  Method 3: $\sum_{|\boldsymbol{k}| \in S_k} E(\boldsymbol{k})$. Method 1: $S_{log,k}$, Method 2/3: $S_{k}$)
    - `skewness2D`: skewness factor calculation using kinetic energy spectrum. `$PARAMETER`: output file step. 
    - `triad3D` or `triad2D` : calculation of kinetic energy transfer spectrum (of triads) $T^{\bullet}(k)$. `$PARAMETER1`: output file step. `$PARAMETER2`: spectral calculation method, same as `instant3D`.
    - `initparam3D` or `initparam2D`: calculation of initial field parameter(initial integral length, large-eddy-turnover time *etc.*) using spectral method. No `$PARAMETER`. 
2. [`SGS`](./udf_pp_SGS.F90): energy transfer calculation in physical space
    The calculation of this part requires a [`datin/SGSinput`](#datinsgsinput-explanation)
    - `Pi3Dint` or `Pi2Dint`: calculation of energy transfer using integral equivalent formula, cf. Eq.(2.9). `$PARAMETER`: output file step.  
    - `stress3D` or `stress2D`: calculation of SGS stress using integral equivalent formula, cf. Eq.(2.9). `$PARAMETER`: output file step. 
    - `Pi3Dtot` or `Pi2Dtot` : calculation of energy transfer using definition, cf. Eq.(2.1). `$PARAMETER`: output file step. 
    - `Pi3Dlocal` or `Pi2Dlocal` : calculation of energy transfer using local formula, cf. Eq.(2.13). `$PARAMETER`: output file step. 
    -  `LES3D`: calculation of LES judgement formula, cf. section 4.4. `$PARAMETER`:output file step. 
    - `T3D`: calculation of energy transfer term. `$PARAMETER`: output file step.  
    - `PiOmega2D`: calculation of SGS enstrophy transfer. `$PARAMETER`: output file step. 
3. [`velgradient`](./udf_pp_velgrad.F90): velocity gradient calculation
    - `instant`: calculation of velocity gradient of a snapshot. `$PARAMETER`: output file step.  
    - `ScaleLen` : calculation of scale length. `$PARAMETER`: output file step.
    - `vortex2D`: `$PARAMETER`: output file step.
4. [`hitgen`](./udf_pp_hitgen.F90): initial field generation for homogenous isotropic turbulence(HIT). No `$PARAMETER`. 
    - `3DIC4` or `2DIC4`
    - `2Dpic`
    - `2Dherring`

### `$PARAMETER` explanation
- output file step: integer. 
    - `flowfield****.h5`: `$PARAMETER=****`. 
    - *example: `flowfield0022.h5`: `$PARAMETER=22`*
    - `flowfield.h5`, `$PARAMETER=0`

### `datin/SGSinput` explanation
The file should be like following: parameter details cf. [./udf_pp_SGS.F90](./udf_pp_SGS.F90).

```Bash
# Input
#
# num_l,num_alpha,num_alphamin
10,20,15

# ratio_max, ratio_min
32.d0,0.8d0

# loutput
t
```

## Reference
[1] <span id="ref1">[Stepanov, Rodion, et al. "Systematic bias in the calculation of spectral density from a three-dimensional spatial grid." Physical Review E 90.5 (2014): 053309.](https://doi.org/10.1103/PhysRevE.90.053309)</span>