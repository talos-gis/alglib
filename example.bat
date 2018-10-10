@echo off
IF NOT "%~1"=="" GOTO lbl_1_end
type _internal\example-bat-help
EXIT /B 1
:lbl_1_end
IF "%~4"=="" GOTO lbl_2_end
echo Too many parameters specified
echo Did you enclose parameters to be passed to the compiler in double quotes?
EXIT /B 1
:lbl_2_end
IF NOT "%~1"=="list" GOTO lbl_3_end
IF "%~2"=="dforest" GOTO lbl_4_dforest
IF "%~2"=="tsort" GOTO lbl_4_tsort
IF "%~2"=="descriptivestatistics" GOTO lbl_4_descriptivestatistics
IF "%~2"=="bdss" GOTO lbl_4_bdss
IF "%~2"=="kmeans" GOTO lbl_4_kmeans
IF "%~2"=="blas" GOTO lbl_4_blas
IF "%~2"=="lda" GOTO lbl_4_lda
IF "%~2"=="hblas" GOTO lbl_4_hblas
IF "%~2"=="reflections" GOTO lbl_4_reflections
IF "%~2"=="creflections" GOTO lbl_4_creflections
IF "%~2"=="sblas" GOTO lbl_4_sblas
IF "%~2"=="ablasf" GOTO lbl_4_ablasf
IF "%~2"=="ablas" GOTO lbl_4_ablas
IF "%~2"=="ortfac" GOTO lbl_4_ortfac
IF "%~2"=="rotations" GOTO lbl_4_rotations
IF "%~2"=="hsschur" GOTO lbl_4_hsschur
IF "%~2"=="evd" GOTO lbl_4_evd
IF "%~2"=="hqrnd" GOTO lbl_4_hqrnd
IF "%~2"=="matgen" GOTO lbl_4_matgen
IF "%~2"=="trfac" GOTO lbl_4_trfac
IF "%~2"=="trlinsolve" GOTO lbl_4_trlinsolve
IF "%~2"=="safesolve" GOTO lbl_4_safesolve
IF "%~2"=="rcond" GOTO lbl_4_rcond
IF "%~2"=="matinv" GOTO lbl_4_matinv
IF "%~2"=="linreg" GOTO lbl_4_linreg
IF "%~2"=="gammafunc" GOTO lbl_4_gammafunc
IF "%~2"=="normaldistr" GOTO lbl_4_normaldistr
IF "%~2"=="igammaf" GOTO lbl_4_igammaf
IF "%~2"=="bdsvd" GOTO lbl_4_bdsvd
IF "%~2"=="svd" GOTO lbl_4_svd
IF "%~2"=="logit" GOTO lbl_4_logit
IF "%~2"=="mlpbase" GOTO lbl_4_mlpbase
IF "%~2"=="xblas" GOTO lbl_4_xblas
IF "%~2"=="densesolver" GOTO lbl_4_densesolver
IF "%~2"=="mlpe" GOTO lbl_4_mlpe
IF "%~2"=="linmin" GOTO lbl_4_linmin
IF "%~2"=="minlbfgs" GOTO lbl_4_minlbfgs
IF "%~2"=="mlptrain" GOTO lbl_4_mlptrain
IF "%~2"=="pca" GOTO lbl_4_pca
IF "%~2"=="odesolver" GOTO lbl_4_odesolver
IF "%~2"=="conv" GOTO lbl_4_conv
IF "%~2"=="ftbase" GOTO lbl_4_ftbase
IF "%~2"=="fft" GOTO lbl_4_fft
IF "%~2"=="corr" GOTO lbl_4_corr
IF "%~2"=="fht" GOTO lbl_4_fht
IF "%~2"=="autogk" GOTO lbl_4_autogk
IF "%~2"=="gq" GOTO lbl_4_gq
IF "%~2"=="gkq" GOTO lbl_4_gkq
IF "%~2"=="lsfit" GOTO lbl_4_lsfit
IF "%~2"=="minlm" GOTO lbl_4_minlm
IF "%~2"=="polint" GOTO lbl_4_polint
IF "%~2"=="ratinterpolation" GOTO lbl_4_ratinterpolation
IF "%~2"=="ratint" GOTO lbl_4_ratint
IF "%~2"=="taskgen" GOTO lbl_4_taskgen
IF "%~2"=="spline2d" GOTO lbl_4_spline2d
IF "%~2"=="spline3" GOTO lbl_4_spline3
IF "%~2"=="spline1d" GOTO lbl_4_spline1d
IF "%~2"=="idwint" GOTO lbl_4_idwint
IF "%~2"=="nearestneighbor" GOTO lbl_4_nearestneighbor
IF "%~2"=="matdet" GOTO lbl_4_matdet
IF "%~2"=="sdet" GOTO lbl_4_sdet
IF "%~2"=="ldlt" GOTO lbl_4_ldlt
IF "%~2"=="spdgevd" GOTO lbl_4_spdgevd
IF "%~2"=="sinverse" GOTO lbl_4_sinverse
IF "%~2"=="inverseupdate" GOTO lbl_4_inverseupdate
IF "%~2"=="srcond" GOTO lbl_4_srcond
IF "%~2"=="ssolve" GOTO lbl_4_ssolve
IF "%~2"=="estnorm" GOTO lbl_4_estnorm
IF "%~2"=="schur" GOTO lbl_4_schur
IF "%~2"=="mincg" GOTO lbl_4_mincg
IF "%~2"=="minasa" GOTO lbl_4_minasa
IF "%~2"=="airyf" GOTO lbl_4_airyf
IF "%~2"=="bessel" GOTO lbl_4_bessel
IF "%~2"=="betaf" GOTO lbl_4_betaf
IF "%~2"=="chebyshev" GOTO lbl_4_chebyshev
IF "%~2"=="dawson" GOTO lbl_4_dawson
IF "%~2"=="elliptic" GOTO lbl_4_elliptic
IF "%~2"=="expintegrals" GOTO lbl_4_expintegrals
IF "%~2"=="fresnel" GOTO lbl_4_fresnel
IF "%~2"=="hermite" GOTO lbl_4_hermite
IF "%~2"=="ibetaf" GOTO lbl_4_ibetaf
IF "%~2"=="jacobianelliptic" GOTO lbl_4_jacobianelliptic
IF "%~2"=="laguerre" GOTO lbl_4_laguerre
IF "%~2"=="legendre" GOTO lbl_4_legendre
IF "%~2"=="psif" GOTO lbl_4_psif
IF "%~2"=="trigintegrals" GOTO lbl_4_trigintegrals
IF "%~2"=="binomialdistr" GOTO lbl_4_binomialdistr
IF "%~2"=="nearunityunit" GOTO lbl_4_nearunityunit
IF "%~2"=="chisquaredistr" GOTO lbl_4_chisquaredistr
IF "%~2"=="correlation" GOTO lbl_4_correlation
IF "%~2"=="fdistr" GOTO lbl_4_fdistr
IF "%~2"=="correlationtests" GOTO lbl_4_correlationtests
IF "%~2"=="studenttdistr" GOTO lbl_4_studenttdistr
IF "%~2"=="jarquebera" GOTO lbl_4_jarquebera
IF "%~2"=="mannwhitneyu" GOTO lbl_4_mannwhitneyu
IF "%~2"=="poissondistr" GOTO lbl_4_poissondistr
IF "%~2"=="stest" GOTO lbl_4_stest
IF "%~2"=="studentttests" GOTO lbl_4_studentttests
IF "%~2"=="variancetests" GOTO lbl_4_variancetests
IF "%~2"=="wsr" GOTO lbl_4_wsr
GOTO lbl_4___error
:lbl_4_dforest
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_tsort
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_descriptivestatistics
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_bdss
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_kmeans
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_blas
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_lda
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_hblas
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_reflections
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_creflections
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_sblas
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_ablasf
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_ablas
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_ortfac
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_rotations
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_hsschur
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_evd
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_hqrnd
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_matgen
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_trfac
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_trlinsolve
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_safesolve
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_rcond
echo rcond_1             condition number of the Vandermonde matrix
GOTO lbl_4___end
:lbl_4_matinv
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_linreg
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_gammafunc
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_normaldistr
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_igammaf
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_bdsvd
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_svd
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_logit
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_mlpbase
echo mlp_process         MLPProcess() example
echo mlp_process_cls     MLPProcess() example
echo mlp_randomize       MLPRandomize() example
echo mlp_serialize       MLPSerialize()/MLPUnserialize() example
GOTO lbl_4___end
:lbl_4_xblas
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_densesolver
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_mlpe
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_linmin
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_minlbfgs
echo minlbfgs_1          L-BFGS optimization #1
echo minlbfgs_2          L-BFGS optimization #2
GOTO lbl_4___end
:lbl_4_mlptrain
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_pca
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_odesolver
echo ode_example1        ODE solver: example #1
echo ode_example2        ODE solver: example #2
GOTO lbl_4___end
:lbl_4_conv
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_ftbase
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_fft
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_corr
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_fht
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_autogk
echo autogk_singular     1-dimensional integration of singular functions
echo autogk_smooth       1-dimensional integration
GOTO lbl_4___end
:lbl_4_gq
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_gkq
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_lsfit
echo lsfit_linear        Linear least squares fitting
echo lsfit_nonlinear     Nonlinear least squares fitting, Hessian-free
echo lsfit_nonlinear2    Nonlinear least squares fitting, using Hessian
GOTO lbl_4___end
:lbl_4_minlm
echo minlm_fgh           Levenberg-Marquardt optimization, FGH scheme
echo minlm_fgj           Levenberg-Marquardt optimization, FGJ scheme
echo minlm_fj            Levenberg-Marquardt optimization, FJ scheme
echo minlm_fj2           Levenberg-Marquardt optimization, FJ scheme
GOTO lbl_4___end
:lbl_4_polint
echo polint_cheb1        1-dimensional Chebyshev polynomial interpolation
echo polint_cheb2        1-dimensional Chebyshev polynomial interpolation
echo polint_eqdist       1-dimensional equidistant polynomial interpolation
echo polint_fit          1-dimensional polynomial fitting
echo polint_gen          1-dimensional polynomial interpolation
GOTO lbl_4___end
:lbl_4_ratinterpolation
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_ratint
echo ratint_fit          Rational fitting
GOTO lbl_4___end
:lbl_4_taskgen
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_spline2d
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_spline3
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_spline1d
echo spline1d_calc       1-dimensional spline operations
echo spline1d_cubic      cubic spline interpolation
echo spline1d_fit        1-dimensional spline fitting
echo spline1d_fitc       1-dimensional constrained spline fitting
echo spline1d_hermite    1-dimensional Hermite spline interpolation
echo spline1d_linear     1-dimensional linear spline interpolation
GOTO lbl_4___end
:lbl_4_idwint
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_nearestneighbor
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_matdet
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_sdet
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_ldlt
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_spdgevd
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_sinverse
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_inverseupdate
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_srcond
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_ssolve
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_estnorm
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_schur
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_mincg
echo mincg_2             CG optimization #2 (StpMax demo)
echo mincg_1             CG optimization #1
GOTO lbl_4___end
:lbl_4_minasa
echo minasa_1            Bound constrained optimization #1
echo minasa_2            Bound constrained optimization #2
GOTO lbl_4___end
:lbl_4_airyf
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_bessel
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_betaf
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_chebyshev
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_dawson
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_elliptic
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_expintegrals
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_fresnel
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_hermite
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_ibetaf
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_jacobianelliptic
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_laguerre
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_legendre
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_psif
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_trigintegrals
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_binomialdistr
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_nearunityunit
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_chisquaredistr
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_correlation
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_fdistr
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_correlationtests
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_studenttdistr
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_jarquebera
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_mannwhitneyu
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_poissondistr
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_stest
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_studentttests
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_variancetests
echo No examples for this unit
GOTO lbl_4___end
:lbl_4_wsr
echo No examples for this unit
GOTO lbl_4___end
:lbl_4___error
echo unknown unit
EXIT /B 1
:lbl_4___end
EXIT /B 0
:lbl_3_end
IF NOT "%~1"=="view" GOTO lbl_5_end
IF "%~2"=="autogk_singular" GOTO lbl_6_autogk_singular
IF "%~2"=="autogk_smooth" GOTO lbl_6_autogk_smooth
IF "%~2"=="minlbfgs_1" GOTO lbl_6_minlbfgs_1
IF "%~2"=="lsfit_linear" GOTO lbl_6_lsfit_linear
IF "%~2"=="lsfit_nonlinear" GOTO lbl_6_lsfit_nonlinear
IF "%~2"=="lsfit_nonlinear2" GOTO lbl_6_lsfit_nonlinear2
IF "%~2"=="minlm_fgh" GOTO lbl_6_minlm_fgh
IF "%~2"=="minlm_fgj" GOTO lbl_6_minlm_fgj
IF "%~2"=="minlm_fj" GOTO lbl_6_minlm_fj
IF "%~2"=="minlm_fj2" GOTO lbl_6_minlm_fj2
IF "%~2"=="mlp_process" GOTO lbl_6_mlp_process
IF "%~2"=="mlp_process_cls" GOTO lbl_6_mlp_process_cls
IF "%~2"=="mlp_randomize" GOTO lbl_6_mlp_randomize
IF "%~2"=="mlp_serialize" GOTO lbl_6_mlp_serialize
IF "%~2"=="ode_example1" GOTO lbl_6_ode_example1
IF "%~2"=="ode_example2" GOTO lbl_6_ode_example2
IF "%~2"=="polint_cheb1" GOTO lbl_6_polint_cheb1
IF "%~2"=="polint_cheb2" GOTO lbl_6_polint_cheb2
IF "%~2"=="polint_eqdist" GOTO lbl_6_polint_eqdist
IF "%~2"=="polint_fit" GOTO lbl_6_polint_fit
IF "%~2"=="polint_gen" GOTO lbl_6_polint_gen
IF "%~2"=="ratint_fit" GOTO lbl_6_ratint_fit
IF "%~2"=="rcond_1" GOTO lbl_6_rcond_1
IF "%~2"=="spline1d_calc" GOTO lbl_6_spline1d_calc
IF "%~2"=="spline1d_cubic" GOTO lbl_6_spline1d_cubic
IF "%~2"=="spline1d_fit" GOTO lbl_6_spline1d_fit
IF "%~2"=="spline1d_fitc" GOTO lbl_6_spline1d_fitc
IF "%~2"=="spline1d_hermite" GOTO lbl_6_spline1d_hermite
IF "%~2"=="spline1d_linear" GOTO lbl_6_spline1d_linear
IF "%~2"=="minlbfgs_2" GOTO lbl_6_minlbfgs_2
IF "%~2"=="mincg_2" GOTO lbl_6_mincg_2
IF "%~2"=="mincg_1" GOTO lbl_6_mincg_1
IF "%~2"=="minasa_1" GOTO lbl_6_minasa_1
IF "%~2"=="minasa_2" GOTO lbl_6_minasa_2
GOTO lbl_6___error
:lbl_6_autogk_singular
type examples\_demo_autogk_singular.*
GOTO lbl_6___end
:lbl_6_autogk_smooth
type examples\_demo_autogk_smooth.*
GOTO lbl_6___end
:lbl_6_minlbfgs_1
type examples\_demo_minlbfgs_1.*
GOTO lbl_6___end
:lbl_6_lsfit_linear
type examples\_demo_lsfit_linear.*
GOTO lbl_6___end
:lbl_6_lsfit_nonlinear
type examples\_demo_lsfit_nonlinear.*
GOTO lbl_6___end
:lbl_6_lsfit_nonlinear2
type examples\_demo_lsfit_nonlinear2.*
GOTO lbl_6___end
:lbl_6_minlm_fgh
type examples\_demo_minlm_fgh.*
GOTO lbl_6___end
:lbl_6_minlm_fgj
type examples\_demo_minlm_fgj.*
GOTO lbl_6___end
:lbl_6_minlm_fj
type examples\_demo_minlm_fj.*
GOTO lbl_6___end
:lbl_6_minlm_fj2
type examples\_demo_minlm_fj2.*
GOTO lbl_6___end
:lbl_6_mlp_process
type examples\_demo_mlp_process.*
GOTO lbl_6___end
:lbl_6_mlp_process_cls
type examples\_demo_mlp_process_cls.*
GOTO lbl_6___end
:lbl_6_mlp_randomize
type examples\_demo_mlp_randomize.*
GOTO lbl_6___end
:lbl_6_mlp_serialize
type examples\_demo_mlp_serialize.*
GOTO lbl_6___end
:lbl_6_ode_example1
type examples\_demo_ode_example1.*
GOTO lbl_6___end
:lbl_6_ode_example2
type examples\_demo_ode_example2.*
GOTO lbl_6___end
:lbl_6_polint_cheb1
type examples\_demo_polint_cheb1.*
GOTO lbl_6___end
:lbl_6_polint_cheb2
type examples\_demo_polint_cheb2.*
GOTO lbl_6___end
:lbl_6_polint_eqdist
type examples\_demo_polint_eqdist.*
GOTO lbl_6___end
:lbl_6_polint_fit
type examples\_demo_polint_fit.*
GOTO lbl_6___end
:lbl_6_polint_gen
type examples\_demo_polint_gen.*
GOTO lbl_6___end
:lbl_6_ratint_fit
type examples\_demo_ratint_fit.*
GOTO lbl_6___end
:lbl_6_rcond_1
type examples\_demo_rcond_1.*
GOTO lbl_6___end
:lbl_6_spline1d_calc
type examples\_demo_spline1d_calc.*
GOTO lbl_6___end
:lbl_6_spline1d_cubic
type examples\_demo_spline1d_cubic.*
GOTO lbl_6___end
:lbl_6_spline1d_fit
type examples\_demo_spline1d_fit.*
GOTO lbl_6___end
:lbl_6_spline1d_fitc
type examples\_demo_spline1d_fitc.*
GOTO lbl_6___end
:lbl_6_spline1d_hermite
type examples\_demo_spline1d_hermite.*
GOTO lbl_6___end
:lbl_6_spline1d_linear
type examples\_demo_spline1d_linear.*
GOTO lbl_6___end
:lbl_6_minlbfgs_2
type examples\_demo_minlbfgs_2.*
GOTO lbl_6___end
:lbl_6_mincg_2
type examples\_demo_mincg_2.*
GOTO lbl_6___end
:lbl_6_mincg_1
type examples\_demo_mincg_1.*
GOTO lbl_6___end
:lbl_6_minasa_1
type examples\_demo_minasa_1.*
GOTO lbl_6___end
:lbl_6_minasa_2
type examples\_demo_minasa_2.*
GOTO lbl_6___end
:lbl_6___error
echo unknown example
EXIT /B 1
:lbl_6___end
EXIT /B 0
:lbl_5_end
CHDIR _tmp
IF EXIST * DEL /F /Q *
CHDIR ..
IF "%~1"=="dcc32" GOTO lbl_7_dcc32
GOTO lbl_7___error
:lbl_7_dcc32
COPY src\*.pas _tmp > NUL 2> NUL
IF NOT ERRORLEVEL 1 GOTO lbl_8
echo Error copying ALGLIB source.
EXIT /B 1
:lbl_8
GOTO lbl_7___end
:lbl_7___error
echo unknown compiler
EXIT /B 1
:lbl_7___end
IF "%~2"=="autogk_singular" GOTO lbl_9_autogk_singular
IF "%~2"=="autogk_smooth" GOTO lbl_9_autogk_smooth
IF "%~2"=="minlbfgs_1" GOTO lbl_9_minlbfgs_1
IF "%~2"=="lsfit_linear" GOTO lbl_9_lsfit_linear
IF "%~2"=="lsfit_nonlinear" GOTO lbl_9_lsfit_nonlinear
IF "%~2"=="lsfit_nonlinear2" GOTO lbl_9_lsfit_nonlinear2
IF "%~2"=="minlm_fgh" GOTO lbl_9_minlm_fgh
IF "%~2"=="minlm_fgj" GOTO lbl_9_minlm_fgj
IF "%~2"=="minlm_fj" GOTO lbl_9_minlm_fj
IF "%~2"=="minlm_fj2" GOTO lbl_9_minlm_fj2
IF "%~2"=="mlp_process" GOTO lbl_9_mlp_process
IF "%~2"=="mlp_process_cls" GOTO lbl_9_mlp_process_cls
IF "%~2"=="mlp_randomize" GOTO lbl_9_mlp_randomize
IF "%~2"=="mlp_serialize" GOTO lbl_9_mlp_serialize
IF "%~2"=="ode_example1" GOTO lbl_9_ode_example1
IF "%~2"=="ode_example2" GOTO lbl_9_ode_example2
IF "%~2"=="polint_cheb1" GOTO lbl_9_polint_cheb1
IF "%~2"=="polint_cheb2" GOTO lbl_9_polint_cheb2
IF "%~2"=="polint_eqdist" GOTO lbl_9_polint_eqdist
IF "%~2"=="polint_fit" GOTO lbl_9_polint_fit
IF "%~2"=="polint_gen" GOTO lbl_9_polint_gen
IF "%~2"=="ratint_fit" GOTO lbl_9_ratint_fit
IF "%~2"=="rcond_1" GOTO lbl_9_rcond_1
IF "%~2"=="spline1d_calc" GOTO lbl_9_spline1d_calc
IF "%~2"=="spline1d_cubic" GOTO lbl_9_spline1d_cubic
IF "%~2"=="spline1d_fit" GOTO lbl_9_spline1d_fit
IF "%~2"=="spline1d_fitc" GOTO lbl_9_spline1d_fitc
IF "%~2"=="spline1d_hermite" GOTO lbl_9_spline1d_hermite
IF "%~2"=="spline1d_linear" GOTO lbl_9_spline1d_linear
IF "%~2"=="minlbfgs_2" GOTO lbl_9_minlbfgs_2
IF "%~2"=="mincg_2" GOTO lbl_9_mincg_2
IF "%~2"=="mincg_1" GOTO lbl_9_mincg_1
IF "%~2"=="minasa_1" GOTO lbl_9_minasa_1
IF "%~2"=="minasa_2" GOTO lbl_9_minasa_2
GOTO lbl_9___error
:lbl_9_autogk_singular
COPY examples\_demo_autogk_singular.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_autogk_smooth
COPY examples\_demo_autogk_smooth.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_minlbfgs_1
COPY examples\_demo_minlbfgs_1.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_lsfit_linear
COPY examples\_demo_lsfit_linear.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_lsfit_nonlinear
COPY examples\_demo_lsfit_nonlinear.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_lsfit_nonlinear2
COPY examples\_demo_lsfit_nonlinear2.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_minlm_fgh
COPY examples\_demo_minlm_fgh.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_minlm_fgj
COPY examples\_demo_minlm_fgj.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_minlm_fj
COPY examples\_demo_minlm_fj.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_minlm_fj2
COPY examples\_demo_minlm_fj2.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_mlp_process
COPY examples\_demo_mlp_process.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_mlp_process_cls
COPY examples\_demo_mlp_process_cls.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_mlp_randomize
COPY examples\_demo_mlp_randomize.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_mlp_serialize
COPY examples\_demo_mlp_serialize.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_ode_example1
COPY examples\_demo_ode_example1.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_ode_example2
COPY examples\_demo_ode_example2.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_polint_cheb1
COPY examples\_demo_polint_cheb1.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_polint_cheb2
COPY examples\_demo_polint_cheb2.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_polint_eqdist
COPY examples\_demo_polint_eqdist.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_polint_fit
COPY examples\_demo_polint_fit.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_polint_gen
COPY examples\_demo_polint_gen.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_ratint_fit
COPY examples\_demo_ratint_fit.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_rcond_1
COPY examples\_demo_rcond_1.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_spline1d_calc
COPY examples\_demo_spline1d_calc.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_spline1d_cubic
COPY examples\_demo_spline1d_cubic.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_spline1d_fit
COPY examples\_demo_spline1d_fit.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_spline1d_fitc
COPY examples\_demo_spline1d_fitc.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_spline1d_hermite
COPY examples\_demo_spline1d_hermite.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_spline1d_linear
COPY examples\_demo_spline1d_linear.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_minlbfgs_2
COPY examples\_demo_minlbfgs_2.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_mincg_2
COPY examples\_demo_mincg_2.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_mincg_1
COPY examples\_demo_mincg_1.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_minasa_1
COPY examples\_demo_minasa_1.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9_minasa_2
COPY examples\_demo_minasa_2.pas _tmp\_demo.pas > NUL 2> NUL
GOTO lbl_9___end
:lbl_9___error
echo unknown example
EXIT /B 1
:lbl_9___end
IF "%~1"=="dcc32" GOTO lbl_10_dcc32
GOTO lbl_10___error
:lbl_10_dcc32
CHDIR _tmp
dcc32 -$C+ -CC -$Q+ -$R+ -$X+ %~3 _demo.pas >> ../log.txt 2>&1 2>&1
IF NOT ERRORLEVEL 1 GOTO lbl_11
echo Error while compiling (see ../log.txt for more info)
CHDIR ..
EXIT /B 1
:lbl_11
CHDIR ..
pushd _tmp
.\_demo
popd
GOTO lbl_10___end
:lbl_10___error
echo unknown compiler
EXIT /B 1
:lbl_10___end
EXIT /B 0
