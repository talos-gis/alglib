@echo off
IF NOT "%~1"=="" GOTO lbl_1_end
type _internal\check-bat-help
EXIT /B 1
:lbl_1_end
IF "%~4"=="" GOTO lbl_2_end
echo Too many parameters specified
echo Did you enclose parameters to be passed to the compiler in double quotes?
EXIT /B 1
:lbl_2_end
IF NOT "%~2"=="all" GOTO lbl_3_end
pushd .
CALL .\check  "%~1"  "dforest_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_4
EXIT /B 1
:lbl_4
pushd .
CALL .\check  "%~1"  "tsort_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_5
EXIT /B 1
:lbl_5
pushd .
CALL .\check  "%~1"  "descriptivestatistics_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_6
EXIT /B 1
:lbl_6
pushd .
CALL .\check  "%~1"  "bdss_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_7
EXIT /B 1
:lbl_7
pushd .
CALL .\check  "%~1"  "kmeans_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_8
EXIT /B 1
:lbl_8
pushd .
CALL .\check  "%~1"  "blas_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_9
EXIT /B 1
:lbl_9
pushd .
CALL .\check  "%~1"  "lda_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_10
EXIT /B 1
:lbl_10
pushd .
CALL .\check  "%~1"  "hblas_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_11
EXIT /B 1
:lbl_11
pushd .
CALL .\check  "%~1"  "reflections_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_12
EXIT /B 1
:lbl_12
pushd .
CALL .\check  "%~1"  "creflections_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_13
EXIT /B 1
:lbl_13
pushd .
CALL .\check  "%~1"  "sblas_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_14
EXIT /B 1
:lbl_14
pushd .
CALL .\check  "%~1"  "ablasf_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_15
EXIT /B 1
:lbl_15
pushd .
CALL .\check  "%~1"  "ablas_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_16
EXIT /B 1
:lbl_16
pushd .
CALL .\check  "%~1"  "ortfac_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_17
EXIT /B 1
:lbl_17
pushd .
CALL .\check  "%~1"  "rotations_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_18
EXIT /B 1
:lbl_18
pushd .
CALL .\check  "%~1"  "hsschur_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_19
EXIT /B 1
:lbl_19
pushd .
CALL .\check  "%~1"  "evd_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_20
EXIT /B 1
:lbl_20
pushd .
CALL .\check  "%~1"  "hqrnd_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_21
EXIT /B 1
:lbl_21
pushd .
CALL .\check  "%~1"  "matgen_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_22
EXIT /B 1
:lbl_22
pushd .
CALL .\check  "%~1"  "trfac_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_23
EXIT /B 1
:lbl_23
pushd .
CALL .\check  "%~1"  "trlinsolve_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_24
EXIT /B 1
:lbl_24
pushd .
CALL .\check  "%~1"  "safesolve_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_25
EXIT /B 1
:lbl_25
pushd .
CALL .\check  "%~1"  "rcond_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_26
EXIT /B 1
:lbl_26
pushd .
CALL .\check  "%~1"  "matinv_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_27
EXIT /B 1
:lbl_27
pushd .
CALL .\check  "%~1"  "linreg_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_28
EXIT /B 1
:lbl_28
pushd .
CALL .\check  "%~1"  "gammafunc_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_29
EXIT /B 1
:lbl_29
pushd .
CALL .\check  "%~1"  "normaldistr_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_30
EXIT /B 1
:lbl_30
pushd .
CALL .\check  "%~1"  "igammaf_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_31
EXIT /B 1
:lbl_31
pushd .
CALL .\check  "%~1"  "bdsvd_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_32
EXIT /B 1
:lbl_32
pushd .
CALL .\check  "%~1"  "svd_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_33
EXIT /B 1
:lbl_33
pushd .
CALL .\check  "%~1"  "logit_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_34
EXIT /B 1
:lbl_34
pushd .
CALL .\check  "%~1"  "mlpbase_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_35
EXIT /B 1
:lbl_35
pushd .
CALL .\check  "%~1"  "xblas_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_36
EXIT /B 1
:lbl_36
pushd .
CALL .\check  "%~1"  "densesolver_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_37
EXIT /B 1
:lbl_37
pushd .
CALL .\check  "%~1"  "mlpe_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_38
EXIT /B 1
:lbl_38
pushd .
CALL .\check  "%~1"  "lbfgs_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_39
EXIT /B 1
:lbl_39
pushd .
CALL .\check  "%~1"  "mlptrain_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_40
EXIT /B 1
:lbl_40
pushd .
CALL .\check  "%~1"  "pca_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_41
EXIT /B 1
:lbl_41
pushd .
CALL .\check  "%~1"  "odesolver_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_42
EXIT /B 1
:lbl_42
pushd .
CALL .\check  "%~1"  "conv_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_43
EXIT /B 1
:lbl_43
pushd .
CALL .\check  "%~1"  "ftbase_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_44
EXIT /B 1
:lbl_44
pushd .
CALL .\check  "%~1"  "fft_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_45
EXIT /B 1
:lbl_45
pushd .
CALL .\check  "%~1"  "corr_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_46
EXIT /B 1
:lbl_46
pushd .
CALL .\check  "%~1"  "fht_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_47
EXIT /B 1
:lbl_47
pushd .
CALL .\check  "%~1"  "autogk_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_48
EXIT /B 1
:lbl_48
pushd .
CALL .\check  "%~1"  "gq_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_49
EXIT /B 1
:lbl_49
pushd .
CALL .\check  "%~1"  "gkq_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_50
EXIT /B 1
:lbl_50
pushd .
CALL .\check  "%~1"  "lsfit_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_51
EXIT /B 1
:lbl_51
pushd .
CALL .\check  "%~1"  "minlm_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_52
EXIT /B 1
:lbl_52
pushd .
CALL .\check  "%~1"  "polint_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_53
EXIT /B 1
:lbl_53
pushd .
CALL .\check  "%~1"  "ratinterpolation_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_54
EXIT /B 1
:lbl_54
pushd .
CALL .\check  "%~1"  "ratint_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_55
EXIT /B 1
:lbl_55
pushd .
CALL .\check  "%~1"  "taskgen_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_56
EXIT /B 1
:lbl_56
pushd .
CALL .\check  "%~1"  "spline2d_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_57
EXIT /B 1
:lbl_57
pushd .
CALL .\check  "%~1"  "spline3_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_58
EXIT /B 1
:lbl_58
pushd .
CALL .\check  "%~1"  "spline1d_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_59
EXIT /B 1
:lbl_59
pushd .
CALL .\check  "%~1"  "idwint_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_60
EXIT /B 1
:lbl_60
pushd .
CALL .\check  "%~1"  "nearestneighbor_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_61
EXIT /B 1
:lbl_61
pushd .
CALL .\check  "%~1"  "matdet_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_62
EXIT /B 1
:lbl_62
pushd .
CALL .\check  "%~1"  "sdet_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_63
EXIT /B 1
:lbl_63
pushd .
CALL .\check  "%~1"  "ldlt_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_64
EXIT /B 1
:lbl_64
pushd .
CALL .\check  "%~1"  "spdgevd_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_65
EXIT /B 1
:lbl_65
pushd .
CALL .\check  "%~1"  "sinverse_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_66
EXIT /B 1
:lbl_66
pushd .
CALL .\check  "%~1"  "inverseupdate_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_67
EXIT /B 1
:lbl_67
pushd .
CALL .\check  "%~1"  "srcond_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_68
EXIT /B 1
:lbl_68
pushd .
CALL .\check  "%~1"  "ssolve_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_69
EXIT /B 1
:lbl_69
pushd .
CALL .\check  "%~1"  "estnorm_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_70
EXIT /B 1
:lbl_70
pushd .
CALL .\check  "%~1"  "schur_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_71
EXIT /B 1
:lbl_71
pushd .
CALL .\check  "%~1"  "airyf_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_72
EXIT /B 1
:lbl_72
pushd .
CALL .\check  "%~1"  "bessel_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_73
EXIT /B 1
:lbl_73
pushd .
CALL .\check  "%~1"  "betaf_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_74
EXIT /B 1
:lbl_74
pushd .
CALL .\check  "%~1"  "chebyshev_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_75
EXIT /B 1
:lbl_75
pushd .
CALL .\check  "%~1"  "dawson_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_76
EXIT /B 1
:lbl_76
pushd .
CALL .\check  "%~1"  "elliptic_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_77
EXIT /B 1
:lbl_77
pushd .
CALL .\check  "%~1"  "expintegrals_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_78
EXIT /B 1
:lbl_78
pushd .
CALL .\check  "%~1"  "fresnel_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_79
EXIT /B 1
:lbl_79
pushd .
CALL .\check  "%~1"  "hermite_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_80
EXIT /B 1
:lbl_80
pushd .
CALL .\check  "%~1"  "ibetaf_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_81
EXIT /B 1
:lbl_81
pushd .
CALL .\check  "%~1"  "jacobianelliptic_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_82
EXIT /B 1
:lbl_82
pushd .
CALL .\check  "%~1"  "laguerre_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_83
EXIT /B 1
:lbl_83
pushd .
CALL .\check  "%~1"  "legendre_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_84
EXIT /B 1
:lbl_84
pushd .
CALL .\check  "%~1"  "psif_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_85
EXIT /B 1
:lbl_85
pushd .
CALL .\check  "%~1"  "trigintegrals_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_86
EXIT /B 1
:lbl_86
pushd .
CALL .\check  "%~1"  "binomialdistr_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_87
EXIT /B 1
:lbl_87
pushd .
CALL .\check  "%~1"  "nearunityunit_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_88
EXIT /B 1
:lbl_88
pushd .
CALL .\check  "%~1"  "chisquaredistr_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_89
EXIT /B 1
:lbl_89
pushd .
CALL .\check  "%~1"  "correlation_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_90
EXIT /B 1
:lbl_90
pushd .
CALL .\check  "%~1"  "fdistr_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_91
EXIT /B 1
:lbl_91
pushd .
CALL .\check  "%~1"  "correlationtests_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_92
EXIT /B 1
:lbl_92
pushd .
CALL .\check  "%~1"  "studenttdistr_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_93
EXIT /B 1
:lbl_93
pushd .
CALL .\check  "%~1"  "jarquebera_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_94
EXIT /B 1
:lbl_94
pushd .
CALL .\check  "%~1"  "mannwhitneyu_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_95
EXIT /B 1
:lbl_95
pushd .
CALL .\check  "%~1"  "poissondistr_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_96
EXIT /B 1
:lbl_96
pushd .
CALL .\check  "%~1"  "stest_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_97
EXIT /B 1
:lbl_97
pushd .
CALL .\check  "%~1"  "studentttests_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_98
EXIT /B 1
:lbl_98
pushd .
CALL .\check  "%~1"  "variancetests_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_99
EXIT /B 1
:lbl_99
pushd .
CALL .\check  "%~1"  "wsr_short"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_100
EXIT /B 1
:lbl_100
EXIT /B 0
:lbl_3_end
IF NOT "%~2"=="all_silent" GOTO lbl_101_end
pushd .
CALL .\check  "%~1"  "dforest_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_102
EXIT /B 1
:lbl_102
pushd .
CALL .\check  "%~1"  "tsort_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_103
EXIT /B 1
:lbl_103
pushd .
CALL .\check  "%~1"  "descriptivestatistics_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_104
EXIT /B 1
:lbl_104
pushd .
CALL .\check  "%~1"  "bdss_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_105
EXIT /B 1
:lbl_105
pushd .
CALL .\check  "%~1"  "kmeans_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_106
EXIT /B 1
:lbl_106
pushd .
CALL .\check  "%~1"  "blas_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_107
EXIT /B 1
:lbl_107
pushd .
CALL .\check  "%~1"  "lda_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_108
EXIT /B 1
:lbl_108
pushd .
CALL .\check  "%~1"  "hblas_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_109
EXIT /B 1
:lbl_109
pushd .
CALL .\check  "%~1"  "reflections_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_110
EXIT /B 1
:lbl_110
pushd .
CALL .\check  "%~1"  "creflections_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_111
EXIT /B 1
:lbl_111
pushd .
CALL .\check  "%~1"  "sblas_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_112
EXIT /B 1
:lbl_112
pushd .
CALL .\check  "%~1"  "ablasf_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_113
EXIT /B 1
:lbl_113
pushd .
CALL .\check  "%~1"  "ablas_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_114
EXIT /B 1
:lbl_114
pushd .
CALL .\check  "%~1"  "ortfac_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_115
EXIT /B 1
:lbl_115
pushd .
CALL .\check  "%~1"  "rotations_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_116
EXIT /B 1
:lbl_116
pushd .
CALL .\check  "%~1"  "hsschur_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_117
EXIT /B 1
:lbl_117
pushd .
CALL .\check  "%~1"  "evd_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_118
EXIT /B 1
:lbl_118
pushd .
CALL .\check  "%~1"  "hqrnd_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_119
EXIT /B 1
:lbl_119
pushd .
CALL .\check  "%~1"  "matgen_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_120
EXIT /B 1
:lbl_120
pushd .
CALL .\check  "%~1"  "trfac_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_121
EXIT /B 1
:lbl_121
pushd .
CALL .\check  "%~1"  "trlinsolve_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_122
EXIT /B 1
:lbl_122
pushd .
CALL .\check  "%~1"  "safesolve_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_123
EXIT /B 1
:lbl_123
pushd .
CALL .\check  "%~1"  "rcond_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_124
EXIT /B 1
:lbl_124
pushd .
CALL .\check  "%~1"  "matinv_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_125
EXIT /B 1
:lbl_125
pushd .
CALL .\check  "%~1"  "linreg_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_126
EXIT /B 1
:lbl_126
pushd .
CALL .\check  "%~1"  "gammafunc_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_127
EXIT /B 1
:lbl_127
pushd .
CALL .\check  "%~1"  "normaldistr_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_128
EXIT /B 1
:lbl_128
pushd .
CALL .\check  "%~1"  "igammaf_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_129
EXIT /B 1
:lbl_129
pushd .
CALL .\check  "%~1"  "bdsvd_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_130
EXIT /B 1
:lbl_130
pushd .
CALL .\check  "%~1"  "svd_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_131
EXIT /B 1
:lbl_131
pushd .
CALL .\check  "%~1"  "logit_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_132
EXIT /B 1
:lbl_132
pushd .
CALL .\check  "%~1"  "mlpbase_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_133
EXIT /B 1
:lbl_133
pushd .
CALL .\check  "%~1"  "xblas_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_134
EXIT /B 1
:lbl_134
pushd .
CALL .\check  "%~1"  "densesolver_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_135
EXIT /B 1
:lbl_135
pushd .
CALL .\check  "%~1"  "mlpe_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_136
EXIT /B 1
:lbl_136
pushd .
CALL .\check  "%~1"  "lbfgs_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_137
EXIT /B 1
:lbl_137
pushd .
CALL .\check  "%~1"  "mlptrain_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_138
EXIT /B 1
:lbl_138
pushd .
CALL .\check  "%~1"  "pca_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_139
EXIT /B 1
:lbl_139
pushd .
CALL .\check  "%~1"  "odesolver_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_140
EXIT /B 1
:lbl_140
pushd .
CALL .\check  "%~1"  "conv_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_141
EXIT /B 1
:lbl_141
pushd .
CALL .\check  "%~1"  "ftbase_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_142
EXIT /B 1
:lbl_142
pushd .
CALL .\check  "%~1"  "fft_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_143
EXIT /B 1
:lbl_143
pushd .
CALL .\check  "%~1"  "corr_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_144
EXIT /B 1
:lbl_144
pushd .
CALL .\check  "%~1"  "fht_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_145
EXIT /B 1
:lbl_145
pushd .
CALL .\check  "%~1"  "autogk_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_146
EXIT /B 1
:lbl_146
pushd .
CALL .\check  "%~1"  "gq_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_147
EXIT /B 1
:lbl_147
pushd .
CALL .\check  "%~1"  "gkq_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_148
EXIT /B 1
:lbl_148
pushd .
CALL .\check  "%~1"  "lsfit_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_149
EXIT /B 1
:lbl_149
pushd .
CALL .\check  "%~1"  "minlm_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_150
EXIT /B 1
:lbl_150
pushd .
CALL .\check  "%~1"  "polint_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_151
EXIT /B 1
:lbl_151
pushd .
CALL .\check  "%~1"  "ratinterpolation_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_152
EXIT /B 1
:lbl_152
pushd .
CALL .\check  "%~1"  "ratint_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_153
EXIT /B 1
:lbl_153
pushd .
CALL .\check  "%~1"  "taskgen_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_154
EXIT /B 1
:lbl_154
pushd .
CALL .\check  "%~1"  "spline2d_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_155
EXIT /B 1
:lbl_155
pushd .
CALL .\check  "%~1"  "spline3_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_156
EXIT /B 1
:lbl_156
pushd .
CALL .\check  "%~1"  "spline1d_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_157
EXIT /B 1
:lbl_157
pushd .
CALL .\check  "%~1"  "idwint_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_158
EXIT /B 1
:lbl_158
pushd .
CALL .\check  "%~1"  "nearestneighbor_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_159
EXIT /B 1
:lbl_159
pushd .
CALL .\check  "%~1"  "matdet_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_160
EXIT /B 1
:lbl_160
pushd .
CALL .\check  "%~1"  "sdet_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_161
EXIT /B 1
:lbl_161
pushd .
CALL .\check  "%~1"  "ldlt_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_162
EXIT /B 1
:lbl_162
pushd .
CALL .\check  "%~1"  "spdgevd_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_163
EXIT /B 1
:lbl_163
pushd .
CALL .\check  "%~1"  "sinverse_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_164
EXIT /B 1
:lbl_164
pushd .
CALL .\check  "%~1"  "inverseupdate_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_165
EXIT /B 1
:lbl_165
pushd .
CALL .\check  "%~1"  "srcond_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_166
EXIT /B 1
:lbl_166
pushd .
CALL .\check  "%~1"  "ssolve_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_167
EXIT /B 1
:lbl_167
pushd .
CALL .\check  "%~1"  "estnorm_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_168
EXIT /B 1
:lbl_168
pushd .
CALL .\check  "%~1"  "schur_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_169
EXIT /B 1
:lbl_169
pushd .
CALL .\check  "%~1"  "airyf_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_170
EXIT /B 1
:lbl_170
pushd .
CALL .\check  "%~1"  "bessel_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_171
EXIT /B 1
:lbl_171
pushd .
CALL .\check  "%~1"  "betaf_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_172
EXIT /B 1
:lbl_172
pushd .
CALL .\check  "%~1"  "chebyshev_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_173
EXIT /B 1
:lbl_173
pushd .
CALL .\check  "%~1"  "dawson_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_174
EXIT /B 1
:lbl_174
pushd .
CALL .\check  "%~1"  "elliptic_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_175
EXIT /B 1
:lbl_175
pushd .
CALL .\check  "%~1"  "expintegrals_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_176
EXIT /B 1
:lbl_176
pushd .
CALL .\check  "%~1"  "fresnel_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_177
EXIT /B 1
:lbl_177
pushd .
CALL .\check  "%~1"  "hermite_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_178
EXIT /B 1
:lbl_178
pushd .
CALL .\check  "%~1"  "ibetaf_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_179
EXIT /B 1
:lbl_179
pushd .
CALL .\check  "%~1"  "jacobianelliptic_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_180
EXIT /B 1
:lbl_180
pushd .
CALL .\check  "%~1"  "laguerre_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_181
EXIT /B 1
:lbl_181
pushd .
CALL .\check  "%~1"  "legendre_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_182
EXIT /B 1
:lbl_182
pushd .
CALL .\check  "%~1"  "psif_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_183
EXIT /B 1
:lbl_183
pushd .
CALL .\check  "%~1"  "trigintegrals_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_184
EXIT /B 1
:lbl_184
pushd .
CALL .\check  "%~1"  "binomialdistr_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_185
EXIT /B 1
:lbl_185
pushd .
CALL .\check  "%~1"  "nearunityunit_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_186
EXIT /B 1
:lbl_186
pushd .
CALL .\check  "%~1"  "chisquaredistr_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_187
EXIT /B 1
:lbl_187
pushd .
CALL .\check  "%~1"  "correlation_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_188
EXIT /B 1
:lbl_188
pushd .
CALL .\check  "%~1"  "fdistr_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_189
EXIT /B 1
:lbl_189
pushd .
CALL .\check  "%~1"  "correlationtests_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_190
EXIT /B 1
:lbl_190
pushd .
CALL .\check  "%~1"  "studenttdistr_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_191
EXIT /B 1
:lbl_191
pushd .
CALL .\check  "%~1"  "jarquebera_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_192
EXIT /B 1
:lbl_192
pushd .
CALL .\check  "%~1"  "mannwhitneyu_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_193
EXIT /B 1
:lbl_193
pushd .
CALL .\check  "%~1"  "poissondistr_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_194
EXIT /B 1
:lbl_194
pushd .
CALL .\check  "%~1"  "stest_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_195
EXIT /B 1
:lbl_195
pushd .
CALL .\check  "%~1"  "studentttests_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_196
EXIT /B 1
:lbl_196
pushd .
CALL .\check  "%~1"  "variancetests_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_197
EXIT /B 1
:lbl_197
pushd .
CALL .\check  "%~1"  "wsr_silent"  "%~3" 
popd
IF NOT ERRORLEVEL 1 GOTO lbl_198
EXIT /B 1
:lbl_198
EXIT /B 0
:lbl_101_end
CHDIR _tmp
IF EXIST * DEL /F /Q *
CHDIR ..
IF "%~1"=="dcc32" GOTO lbl_199_dcc32
GOTO lbl_199___error
:lbl_199_dcc32
COPY src\*.pas _tmp > NUL 2> NUL
IF NOT ERRORLEVEL 1 GOTO lbl_200
echo Error copying ALGLIB source files.
EXIT /B 1
:lbl_200
GOTO lbl_199___end
:lbl_199___error
echo unknown compiler
EXIT /B 1
:lbl_199___end
IF "%~2"=="dforest" GOTO lbl_201_dforest
IF "%~2"=="dforest_silent" GOTO lbl_201_dforest_silent
IF "%~2"=="dforest_short" GOTO lbl_201_dforest_short
IF "%~2"=="tsort" GOTO lbl_201_tsort
IF "%~2"=="tsort_silent" GOTO lbl_201_tsort_silent
IF "%~2"=="tsort_short" GOTO lbl_201_tsort_short
IF "%~2"=="descriptivestatistics" GOTO lbl_201_descriptivestatistics
IF "%~2"=="descriptivestatistics_silent" GOTO lbl_201_descriptivestatistics_silent
IF "%~2"=="descriptivestatistics_short" GOTO lbl_201_descriptivestatistics_short
IF "%~2"=="bdss" GOTO lbl_201_bdss
IF "%~2"=="bdss_silent" GOTO lbl_201_bdss_silent
IF "%~2"=="bdss_short" GOTO lbl_201_bdss_short
IF "%~2"=="kmeans" GOTO lbl_201_kmeans
IF "%~2"=="kmeans_silent" GOTO lbl_201_kmeans_silent
IF "%~2"=="kmeans_short" GOTO lbl_201_kmeans_short
IF "%~2"=="blas" GOTO lbl_201_blas
IF "%~2"=="blas_silent" GOTO lbl_201_blas_silent
IF "%~2"=="blas_short" GOTO lbl_201_blas_short
IF "%~2"=="lda" GOTO lbl_201_lda
IF "%~2"=="lda_silent" GOTO lbl_201_lda_silent
IF "%~2"=="lda_short" GOTO lbl_201_lda_short
IF "%~2"=="hblas" GOTO lbl_201_hblas
IF "%~2"=="hblas_silent" GOTO lbl_201_hblas_silent
IF "%~2"=="hblas_short" GOTO lbl_201_hblas_short
IF "%~2"=="reflections" GOTO lbl_201_reflections
IF "%~2"=="reflections_silent" GOTO lbl_201_reflections_silent
IF "%~2"=="reflections_short" GOTO lbl_201_reflections_short
IF "%~2"=="creflections" GOTO lbl_201_creflections
IF "%~2"=="creflections_silent" GOTO lbl_201_creflections_silent
IF "%~2"=="creflections_short" GOTO lbl_201_creflections_short
IF "%~2"=="sblas" GOTO lbl_201_sblas
IF "%~2"=="sblas_silent" GOTO lbl_201_sblas_silent
IF "%~2"=="sblas_short" GOTO lbl_201_sblas_short
IF "%~2"=="ablasf" GOTO lbl_201_ablasf
IF "%~2"=="ablasf_silent" GOTO lbl_201_ablasf_silent
IF "%~2"=="ablasf_short" GOTO lbl_201_ablasf_short
IF "%~2"=="ablas" GOTO lbl_201_ablas
IF "%~2"=="ablas_silent" GOTO lbl_201_ablas_silent
IF "%~2"=="ablas_short" GOTO lbl_201_ablas_short
IF "%~2"=="ortfac" GOTO lbl_201_ortfac
IF "%~2"=="ortfac_silent" GOTO lbl_201_ortfac_silent
IF "%~2"=="ortfac_short" GOTO lbl_201_ortfac_short
IF "%~2"=="rotations" GOTO lbl_201_rotations
IF "%~2"=="rotations_silent" GOTO lbl_201_rotations_silent
IF "%~2"=="rotations_short" GOTO lbl_201_rotations_short
IF "%~2"=="hsschur" GOTO lbl_201_hsschur
IF "%~2"=="hsschur_silent" GOTO lbl_201_hsschur_silent
IF "%~2"=="hsschur_short" GOTO lbl_201_hsschur_short
IF "%~2"=="evd" GOTO lbl_201_evd
IF "%~2"=="evd_silent" GOTO lbl_201_evd_silent
IF "%~2"=="evd_short" GOTO lbl_201_evd_short
IF "%~2"=="hqrnd" GOTO lbl_201_hqrnd
IF "%~2"=="hqrnd_silent" GOTO lbl_201_hqrnd_silent
IF "%~2"=="hqrnd_short" GOTO lbl_201_hqrnd_short
IF "%~2"=="matgen" GOTO lbl_201_matgen
IF "%~2"=="matgen_silent" GOTO lbl_201_matgen_silent
IF "%~2"=="matgen_short" GOTO lbl_201_matgen_short
IF "%~2"=="trfac" GOTO lbl_201_trfac
IF "%~2"=="trfac_silent" GOTO lbl_201_trfac_silent
IF "%~2"=="trfac_short" GOTO lbl_201_trfac_short
IF "%~2"=="trlinsolve" GOTO lbl_201_trlinsolve
IF "%~2"=="trlinsolve_silent" GOTO lbl_201_trlinsolve_silent
IF "%~2"=="trlinsolve_short" GOTO lbl_201_trlinsolve_short
IF "%~2"=="safesolve" GOTO lbl_201_safesolve
IF "%~2"=="safesolve_silent" GOTO lbl_201_safesolve_silent
IF "%~2"=="safesolve_short" GOTO lbl_201_safesolve_short
IF "%~2"=="rcond" GOTO lbl_201_rcond
IF "%~2"=="rcond_silent" GOTO lbl_201_rcond_silent
IF "%~2"=="rcond_short" GOTO lbl_201_rcond_short
IF "%~2"=="matinv" GOTO lbl_201_matinv
IF "%~2"=="matinv_silent" GOTO lbl_201_matinv_silent
IF "%~2"=="matinv_short" GOTO lbl_201_matinv_short
IF "%~2"=="linreg" GOTO lbl_201_linreg
IF "%~2"=="linreg_silent" GOTO lbl_201_linreg_silent
IF "%~2"=="linreg_short" GOTO lbl_201_linreg_short
IF "%~2"=="gammafunc" GOTO lbl_201_gammafunc
IF "%~2"=="gammafunc_silent" GOTO lbl_201_gammafunc_silent
IF "%~2"=="gammafunc_short" GOTO lbl_201_gammafunc_short
IF "%~2"=="normaldistr" GOTO lbl_201_normaldistr
IF "%~2"=="normaldistr_silent" GOTO lbl_201_normaldistr_silent
IF "%~2"=="normaldistr_short" GOTO lbl_201_normaldistr_short
IF "%~2"=="igammaf" GOTO lbl_201_igammaf
IF "%~2"=="igammaf_silent" GOTO lbl_201_igammaf_silent
IF "%~2"=="igammaf_short" GOTO lbl_201_igammaf_short
IF "%~2"=="bdsvd" GOTO lbl_201_bdsvd
IF "%~2"=="bdsvd_silent" GOTO lbl_201_bdsvd_silent
IF "%~2"=="bdsvd_short" GOTO lbl_201_bdsvd_short
IF "%~2"=="svd" GOTO lbl_201_svd
IF "%~2"=="svd_silent" GOTO lbl_201_svd_silent
IF "%~2"=="svd_short" GOTO lbl_201_svd_short
IF "%~2"=="logit" GOTO lbl_201_logit
IF "%~2"=="logit_silent" GOTO lbl_201_logit_silent
IF "%~2"=="logit_short" GOTO lbl_201_logit_short
IF "%~2"=="mlpbase" GOTO lbl_201_mlpbase
IF "%~2"=="mlpbase_silent" GOTO lbl_201_mlpbase_silent
IF "%~2"=="mlpbase_short" GOTO lbl_201_mlpbase_short
IF "%~2"=="xblas" GOTO lbl_201_xblas
IF "%~2"=="xblas_silent" GOTO lbl_201_xblas_silent
IF "%~2"=="xblas_short" GOTO lbl_201_xblas_short
IF "%~2"=="densesolver" GOTO lbl_201_densesolver
IF "%~2"=="densesolver_silent" GOTO lbl_201_densesolver_silent
IF "%~2"=="densesolver_short" GOTO lbl_201_densesolver_short
IF "%~2"=="mlpe" GOTO lbl_201_mlpe
IF "%~2"=="mlpe_silent" GOTO lbl_201_mlpe_silent
IF "%~2"=="mlpe_short" GOTO lbl_201_mlpe_short
IF "%~2"=="lbfgs" GOTO lbl_201_lbfgs
IF "%~2"=="lbfgs_silent" GOTO lbl_201_lbfgs_silent
IF "%~2"=="lbfgs_short" GOTO lbl_201_lbfgs_short
IF "%~2"=="mlptrain" GOTO lbl_201_mlptrain
IF "%~2"=="mlptrain_silent" GOTO lbl_201_mlptrain_silent
IF "%~2"=="mlptrain_short" GOTO lbl_201_mlptrain_short
IF "%~2"=="pca" GOTO lbl_201_pca
IF "%~2"=="pca_silent" GOTO lbl_201_pca_silent
IF "%~2"=="pca_short" GOTO lbl_201_pca_short
IF "%~2"=="odesolver" GOTO lbl_201_odesolver
IF "%~2"=="odesolver_silent" GOTO lbl_201_odesolver_silent
IF "%~2"=="odesolver_short" GOTO lbl_201_odesolver_short
IF "%~2"=="conv" GOTO lbl_201_conv
IF "%~2"=="conv_silent" GOTO lbl_201_conv_silent
IF "%~2"=="conv_short" GOTO lbl_201_conv_short
IF "%~2"=="ftbase" GOTO lbl_201_ftbase
IF "%~2"=="ftbase_silent" GOTO lbl_201_ftbase_silent
IF "%~2"=="ftbase_short" GOTO lbl_201_ftbase_short
IF "%~2"=="fft" GOTO lbl_201_fft
IF "%~2"=="fft_silent" GOTO lbl_201_fft_silent
IF "%~2"=="fft_short" GOTO lbl_201_fft_short
IF "%~2"=="corr" GOTO lbl_201_corr
IF "%~2"=="corr_silent" GOTO lbl_201_corr_silent
IF "%~2"=="corr_short" GOTO lbl_201_corr_short
IF "%~2"=="fht" GOTO lbl_201_fht
IF "%~2"=="fht_silent" GOTO lbl_201_fht_silent
IF "%~2"=="fht_short" GOTO lbl_201_fht_short
IF "%~2"=="autogk" GOTO lbl_201_autogk
IF "%~2"=="autogk_silent" GOTO lbl_201_autogk_silent
IF "%~2"=="autogk_short" GOTO lbl_201_autogk_short
IF "%~2"=="gq" GOTO lbl_201_gq
IF "%~2"=="gq_silent" GOTO lbl_201_gq_silent
IF "%~2"=="gq_short" GOTO lbl_201_gq_short
IF "%~2"=="gkq" GOTO lbl_201_gkq
IF "%~2"=="gkq_silent" GOTO lbl_201_gkq_silent
IF "%~2"=="gkq_short" GOTO lbl_201_gkq_short
IF "%~2"=="lsfit" GOTO lbl_201_lsfit
IF "%~2"=="lsfit_silent" GOTO lbl_201_lsfit_silent
IF "%~2"=="lsfit_short" GOTO lbl_201_lsfit_short
IF "%~2"=="minlm" GOTO lbl_201_minlm
IF "%~2"=="minlm_silent" GOTO lbl_201_minlm_silent
IF "%~2"=="minlm_short" GOTO lbl_201_minlm_short
IF "%~2"=="polint" GOTO lbl_201_polint
IF "%~2"=="polint_silent" GOTO lbl_201_polint_silent
IF "%~2"=="polint_short" GOTO lbl_201_polint_short
IF "%~2"=="ratinterpolation" GOTO lbl_201_ratinterpolation
IF "%~2"=="ratinterpolation_silent" GOTO lbl_201_ratinterpolation_silent
IF "%~2"=="ratinterpolation_short" GOTO lbl_201_ratinterpolation_short
IF "%~2"=="ratint" GOTO lbl_201_ratint
IF "%~2"=="ratint_silent" GOTO lbl_201_ratint_silent
IF "%~2"=="ratint_short" GOTO lbl_201_ratint_short
IF "%~2"=="taskgen" GOTO lbl_201_taskgen
IF "%~2"=="taskgen_silent" GOTO lbl_201_taskgen_silent
IF "%~2"=="taskgen_short" GOTO lbl_201_taskgen_short
IF "%~2"=="spline2d" GOTO lbl_201_spline2d
IF "%~2"=="spline2d_silent" GOTO lbl_201_spline2d_silent
IF "%~2"=="spline2d_short" GOTO lbl_201_spline2d_short
IF "%~2"=="spline3" GOTO lbl_201_spline3
IF "%~2"=="spline3_silent" GOTO lbl_201_spline3_silent
IF "%~2"=="spline3_short" GOTO lbl_201_spline3_short
IF "%~2"=="spline1d" GOTO lbl_201_spline1d
IF "%~2"=="spline1d_silent" GOTO lbl_201_spline1d_silent
IF "%~2"=="spline1d_short" GOTO lbl_201_spline1d_short
IF "%~2"=="idwint" GOTO lbl_201_idwint
IF "%~2"=="idwint_silent" GOTO lbl_201_idwint_silent
IF "%~2"=="idwint_short" GOTO lbl_201_idwint_short
IF "%~2"=="nearestneighbor" GOTO lbl_201_nearestneighbor
IF "%~2"=="nearestneighbor_silent" GOTO lbl_201_nearestneighbor_silent
IF "%~2"=="nearestneighbor_short" GOTO lbl_201_nearestneighbor_short
IF "%~2"=="matdet" GOTO lbl_201_matdet
IF "%~2"=="matdet_silent" GOTO lbl_201_matdet_silent
IF "%~2"=="matdet_short" GOTO lbl_201_matdet_short
IF "%~2"=="sdet" GOTO lbl_201_sdet
IF "%~2"=="sdet_silent" GOTO lbl_201_sdet_silent
IF "%~2"=="sdet_short" GOTO lbl_201_sdet_short
IF "%~2"=="ldlt" GOTO lbl_201_ldlt
IF "%~2"=="ldlt_silent" GOTO lbl_201_ldlt_silent
IF "%~2"=="ldlt_short" GOTO lbl_201_ldlt_short
IF "%~2"=="spdgevd" GOTO lbl_201_spdgevd
IF "%~2"=="spdgevd_silent" GOTO lbl_201_spdgevd_silent
IF "%~2"=="spdgevd_short" GOTO lbl_201_spdgevd_short
IF "%~2"=="sinverse" GOTO lbl_201_sinverse
IF "%~2"=="sinverse_silent" GOTO lbl_201_sinverse_silent
IF "%~2"=="sinverse_short" GOTO lbl_201_sinverse_short
IF "%~2"=="inverseupdate" GOTO lbl_201_inverseupdate
IF "%~2"=="inverseupdate_silent" GOTO lbl_201_inverseupdate_silent
IF "%~2"=="inverseupdate_short" GOTO lbl_201_inverseupdate_short
IF "%~2"=="srcond" GOTO lbl_201_srcond
IF "%~2"=="srcond_silent" GOTO lbl_201_srcond_silent
IF "%~2"=="srcond_short" GOTO lbl_201_srcond_short
IF "%~2"=="ssolve" GOTO lbl_201_ssolve
IF "%~2"=="ssolve_silent" GOTO lbl_201_ssolve_silent
IF "%~2"=="ssolve_short" GOTO lbl_201_ssolve_short
IF "%~2"=="estnorm" GOTO lbl_201_estnorm
IF "%~2"=="estnorm_silent" GOTO lbl_201_estnorm_silent
IF "%~2"=="estnorm_short" GOTO lbl_201_estnorm_short
IF "%~2"=="schur" GOTO lbl_201_schur
IF "%~2"=="schur_silent" GOTO lbl_201_schur_silent
IF "%~2"=="schur_short" GOTO lbl_201_schur_short
IF "%~2"=="airyf" GOTO lbl_201_airyf
IF "%~2"=="airyf_silent" GOTO lbl_201_airyf_silent
IF "%~2"=="airyf_short" GOTO lbl_201_airyf_short
IF "%~2"=="bessel" GOTO lbl_201_bessel
IF "%~2"=="bessel_silent" GOTO lbl_201_bessel_silent
IF "%~2"=="bessel_short" GOTO lbl_201_bessel_short
IF "%~2"=="betaf" GOTO lbl_201_betaf
IF "%~2"=="betaf_silent" GOTO lbl_201_betaf_silent
IF "%~2"=="betaf_short" GOTO lbl_201_betaf_short
IF "%~2"=="chebyshev" GOTO lbl_201_chebyshev
IF "%~2"=="chebyshev_silent" GOTO lbl_201_chebyshev_silent
IF "%~2"=="chebyshev_short" GOTO lbl_201_chebyshev_short
IF "%~2"=="dawson" GOTO lbl_201_dawson
IF "%~2"=="dawson_silent" GOTO lbl_201_dawson_silent
IF "%~2"=="dawson_short" GOTO lbl_201_dawson_short
IF "%~2"=="elliptic" GOTO lbl_201_elliptic
IF "%~2"=="elliptic_silent" GOTO lbl_201_elliptic_silent
IF "%~2"=="elliptic_short" GOTO lbl_201_elliptic_short
IF "%~2"=="expintegrals" GOTO lbl_201_expintegrals
IF "%~2"=="expintegrals_silent" GOTO lbl_201_expintegrals_silent
IF "%~2"=="expintegrals_short" GOTO lbl_201_expintegrals_short
IF "%~2"=="fresnel" GOTO lbl_201_fresnel
IF "%~2"=="fresnel_silent" GOTO lbl_201_fresnel_silent
IF "%~2"=="fresnel_short" GOTO lbl_201_fresnel_short
IF "%~2"=="hermite" GOTO lbl_201_hermite
IF "%~2"=="hermite_silent" GOTO lbl_201_hermite_silent
IF "%~2"=="hermite_short" GOTO lbl_201_hermite_short
IF "%~2"=="ibetaf" GOTO lbl_201_ibetaf
IF "%~2"=="ibetaf_silent" GOTO lbl_201_ibetaf_silent
IF "%~2"=="ibetaf_short" GOTO lbl_201_ibetaf_short
IF "%~2"=="jacobianelliptic" GOTO lbl_201_jacobianelliptic
IF "%~2"=="jacobianelliptic_silent" GOTO lbl_201_jacobianelliptic_silent
IF "%~2"=="jacobianelliptic_short" GOTO lbl_201_jacobianelliptic_short
IF "%~2"=="laguerre" GOTO lbl_201_laguerre
IF "%~2"=="laguerre_silent" GOTO lbl_201_laguerre_silent
IF "%~2"=="laguerre_short" GOTO lbl_201_laguerre_short
IF "%~2"=="legendre" GOTO lbl_201_legendre
IF "%~2"=="legendre_silent" GOTO lbl_201_legendre_silent
IF "%~2"=="legendre_short" GOTO lbl_201_legendre_short
IF "%~2"=="psif" GOTO lbl_201_psif
IF "%~2"=="psif_silent" GOTO lbl_201_psif_silent
IF "%~2"=="psif_short" GOTO lbl_201_psif_short
IF "%~2"=="trigintegrals" GOTO lbl_201_trigintegrals
IF "%~2"=="trigintegrals_silent" GOTO lbl_201_trigintegrals_silent
IF "%~2"=="trigintegrals_short" GOTO lbl_201_trigintegrals_short
IF "%~2"=="binomialdistr" GOTO lbl_201_binomialdistr
IF "%~2"=="binomialdistr_silent" GOTO lbl_201_binomialdistr_silent
IF "%~2"=="binomialdistr_short" GOTO lbl_201_binomialdistr_short
IF "%~2"=="nearunityunit" GOTO lbl_201_nearunityunit
IF "%~2"=="nearunityunit_silent" GOTO lbl_201_nearunityunit_silent
IF "%~2"=="nearunityunit_short" GOTO lbl_201_nearunityunit_short
IF "%~2"=="chisquaredistr" GOTO lbl_201_chisquaredistr
IF "%~2"=="chisquaredistr_silent" GOTO lbl_201_chisquaredistr_silent
IF "%~2"=="chisquaredistr_short" GOTO lbl_201_chisquaredistr_short
IF "%~2"=="correlation" GOTO lbl_201_correlation
IF "%~2"=="correlation_silent" GOTO lbl_201_correlation_silent
IF "%~2"=="correlation_short" GOTO lbl_201_correlation_short
IF "%~2"=="fdistr" GOTO lbl_201_fdistr
IF "%~2"=="fdistr_silent" GOTO lbl_201_fdistr_silent
IF "%~2"=="fdistr_short" GOTO lbl_201_fdistr_short
IF "%~2"=="correlationtests" GOTO lbl_201_correlationtests
IF "%~2"=="correlationtests_silent" GOTO lbl_201_correlationtests_silent
IF "%~2"=="correlationtests_short" GOTO lbl_201_correlationtests_short
IF "%~2"=="studenttdistr" GOTO lbl_201_studenttdistr
IF "%~2"=="studenttdistr_silent" GOTO lbl_201_studenttdistr_silent
IF "%~2"=="studenttdistr_short" GOTO lbl_201_studenttdistr_short
IF "%~2"=="jarquebera" GOTO lbl_201_jarquebera
IF "%~2"=="jarquebera_silent" GOTO lbl_201_jarquebera_silent
IF "%~2"=="jarquebera_short" GOTO lbl_201_jarquebera_short
IF "%~2"=="mannwhitneyu" GOTO lbl_201_mannwhitneyu
IF "%~2"=="mannwhitneyu_silent" GOTO lbl_201_mannwhitneyu_silent
IF "%~2"=="mannwhitneyu_short" GOTO lbl_201_mannwhitneyu_short
IF "%~2"=="poissondistr" GOTO lbl_201_poissondistr
IF "%~2"=="poissondistr_silent" GOTO lbl_201_poissondistr_silent
IF "%~2"=="poissondistr_short" GOTO lbl_201_poissondistr_short
IF "%~2"=="stest" GOTO lbl_201_stest
IF "%~2"=="stest_silent" GOTO lbl_201_stest_silent
IF "%~2"=="stest_short" GOTO lbl_201_stest_short
IF "%~2"=="studentttests" GOTO lbl_201_studentttests
IF "%~2"=="studentttests_silent" GOTO lbl_201_studentttests_silent
IF "%~2"=="studentttests_short" GOTO lbl_201_studentttests_short
IF "%~2"=="variancetests" GOTO lbl_201_variancetests
IF "%~2"=="variancetests_silent" GOTO lbl_201_variancetests_silent
IF "%~2"=="variancetests_short" GOTO lbl_201_variancetests_short
IF "%~2"=="wsr" GOTO lbl_201_wsr
IF "%~2"=="wsr_silent" GOTO lbl_201_wsr_silent
IF "%~2"=="wsr_short" GOTO lbl_201_wsr_short
GOTO lbl_201___error
:lbl_201_dforest
COPY _internal\_run_testforestunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testforestunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_dforest_short
COPY _internal\_run_short_testforestunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testforestunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_dforest_silent
COPY _internal\_run_silent_testforestunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testforestunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_tsort
COPY _internal\_run_testtsortunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testtsortunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_tsort_short
COPY _internal\_run_short_testtsortunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testtsortunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_tsort_silent
COPY _internal\_run_silent_testtsortunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testtsortunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_descriptivestatistics
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_descriptivestatistics_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_descriptivestatistics_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_bdss
COPY _internal\_run_testbdssunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testbdssunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_bdss_short
COPY _internal\_run_short_testbdssunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testbdssunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_bdss_silent
COPY _internal\_run_silent_testbdssunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testbdssunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_kmeans
COPY _internal\_run_testkmeansunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testkmeansunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_kmeans_short
COPY _internal\_run_short_testkmeansunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testkmeansunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_kmeans_silent
COPY _internal\_run_silent_testkmeansunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testkmeansunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_blas
COPY _internal\_run_testblasunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testblasunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_blas_short
COPY _internal\_run_short_testblasunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testblasunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_blas_silent
COPY _internal\_run_silent_testblasunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testblasunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_lda
COPY _internal\_run_testldaunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testldaunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_lda_short
COPY _internal\_run_short_testldaunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testldaunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_lda_silent
COPY _internal\_run_silent_testldaunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testldaunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_hblas
COPY _internal\_run_testhblasunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testhblasunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_hblas_short
COPY _internal\_run_short_testhblasunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testhblasunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_hblas_silent
COPY _internal\_run_silent_testhblasunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testhblasunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_reflections
COPY _internal\_run_testreflectionsunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testreflectionsunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_reflections_short
COPY _internal\_run_short_testreflectionsunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testreflectionsunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_reflections_silent
COPY _internal\_run_silent_testreflectionsunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testreflectionsunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_creflections
COPY _internal\_run_testcreflunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testcreflunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_creflections_short
COPY _internal\_run_short_testcreflunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testcreflunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_creflections_silent
COPY _internal\_run_silent_testcreflunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testcreflunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_sblas
COPY _internal\_run_testsblasunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testsblasunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_sblas_short
COPY _internal\_run_short_testsblasunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testsblasunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_sblas_silent
COPY _internal\_run_silent_testsblasunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testsblasunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_ablasf
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_ablasf_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_ablasf_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_ablas
COPY _internal\_run_testablasunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testablasunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_ablas_short
COPY _internal\_run_short_testablasunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testablasunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_ablas_silent
COPY _internal\_run_silent_testablasunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testablasunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_ortfac
COPY _internal\_run_testortfacunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testortfacunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_ortfac_short
COPY _internal\_run_short_testortfacunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testortfacunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_ortfac_silent
COPY _internal\_run_silent_testortfacunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testortfacunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_rotations
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_rotations_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_rotations_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_hsschur
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_hsschur_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_hsschur_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_evd
COPY _internal\_run_testevdunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testevdunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_evd_short
COPY _internal\_run_short_testevdunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testevdunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_evd_silent
COPY _internal\_run_silent_testevdunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testevdunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_hqrnd
COPY _internal\_run_testhqrndunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testhqrndunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_hqrnd_short
COPY _internal\_run_short_testhqrndunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testhqrndunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_hqrnd_silent
COPY _internal\_run_silent_testhqrndunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testhqrndunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_matgen
COPY _internal\_run_testmatgenunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testmatgenunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_matgen_short
COPY _internal\_run_short_testmatgenunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testmatgenunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_matgen_silent
COPY _internal\_run_silent_testmatgenunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testmatgenunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_trfac
COPY _internal\_run_testtrfacunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testtrfacunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_trfac_short
COPY _internal\_run_short_testtrfacunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testtrfacunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_trfac_silent
COPY _internal\_run_silent_testtrfacunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testtrfacunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_trlinsolve
COPY _internal\_run_testsstunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testsstunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_trlinsolve_short
COPY _internal\_run_short_testsstunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testsstunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_trlinsolve_silent
COPY _internal\_run_silent_testsstunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testsstunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_safesolve
COPY _internal\_run_testsafesolveunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testsafesolveunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_safesolve_short
COPY _internal\_run_short_testsafesolveunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testsafesolveunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_safesolve_silent
COPY _internal\_run_silent_testsafesolveunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testsafesolveunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_rcond
COPY _internal\_run_testrcondunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testrcondunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_rcond_short
COPY _internal\_run_short_testrcondunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testrcondunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_rcond_silent
COPY _internal\_run_silent_testrcondunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testrcondunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_matinv
COPY _internal\_run_testmatinvunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testmatinvunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_matinv_short
COPY _internal\_run_short_testmatinvunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testmatinvunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_matinv_silent
COPY _internal\_run_silent_testmatinvunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testmatinvunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_linreg
COPY _internal\_run_testregressunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testregressunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_linreg_short
COPY _internal\_run_short_testregressunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testregressunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_linreg_silent
COPY _internal\_run_silent_testregressunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testregressunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_gammafunc
COPY _internal\_run_testgammaunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testgammaunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_gammafunc_short
COPY _internal\_run_short_testgammaunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testgammaunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_gammafunc_silent
COPY _internal\_run_silent_testgammaunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testgammaunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_normaldistr
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_normaldistr_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_normaldistr_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_igammaf
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_igammaf_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_igammaf_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_bdsvd
COPY _internal\_run_testbdsvdunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testbdsvdunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_bdsvd_short
COPY _internal\_run_short_testbdsvdunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testbdsvdunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_bdsvd_silent
COPY _internal\_run_silent_testbdsvdunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testbdsvdunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_svd
COPY _internal\_run_testsvdunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testsvdunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_svd_short
COPY _internal\_run_short_testsvdunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testsvdunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_svd_silent
COPY _internal\_run_silent_testsvdunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testsvdunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_logit
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_logit_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_logit_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_mlpbase
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_mlpbase_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_mlpbase_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_xblas
COPY _internal\_run_testxblasunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testxblasunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_xblas_short
COPY _internal\_run_short_testxblasunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testxblasunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_xblas_silent
COPY _internal\_run_silent_testxblasunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testxblasunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_densesolver
COPY _internal\_run_testdensesolverunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testdensesolverunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_densesolver_short
COPY _internal\_run_short_testdensesolverunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testdensesolverunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_densesolver_silent
COPY _internal\_run_silent_testdensesolverunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testdensesolverunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_mlpe
COPY _internal\_run_testmlpeunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testmlpeunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_mlpe_short
COPY _internal\_run_short_testmlpeunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testmlpeunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_mlpe_silent
COPY _internal\_run_silent_testmlpeunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testmlpeunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_lbfgs
COPY _internal\_run_testlbfgs.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testlbfgs.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_lbfgs_short
COPY _internal\_run_short_testlbfgs.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testlbfgs.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_lbfgs_silent
COPY _internal\_run_silent_testlbfgs.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testlbfgs.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_mlptrain
COPY _internal\_run_testmlpunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testmlpunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_mlptrain_short
COPY _internal\_run_short_testmlpunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testmlpunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_mlptrain_silent
COPY _internal\_run_silent_testmlpunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testmlpunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_pca
COPY _internal\_run_testpcaunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testpcaunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_pca_short
COPY _internal\_run_short_testpcaunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testpcaunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_pca_silent
COPY _internal\_run_silent_testpcaunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testpcaunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_odesolver
COPY _internal\_run_testodesolverunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testodesolverunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_odesolver_short
COPY _internal\_run_short_testodesolverunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testodesolverunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_odesolver_silent
COPY _internal\_run_silent_testodesolverunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testodesolverunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_conv
COPY _internal\_run_testconvunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testconvunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_conv_short
COPY _internal\_run_short_testconvunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testconvunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_conv_silent
COPY _internal\_run_silent_testconvunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testconvunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_ftbase
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_ftbase_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_ftbase_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_fft
COPY _internal\_run_testfftunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testfftunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_fft_short
COPY _internal\_run_short_testfftunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testfftunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_fft_silent
COPY _internal\_run_silent_testfftunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testfftunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_corr
COPY _internal\_run_testcorrunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testcorrunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_corr_short
COPY _internal\_run_short_testcorrunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testcorrunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_corr_silent
COPY _internal\_run_silent_testcorrunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testcorrunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_fht
COPY _internal\_run_testfhtunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testfhtunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_fht_short
COPY _internal\_run_short_testfhtunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testfhtunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_fht_silent
COPY _internal\_run_silent_testfhtunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testfhtunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_autogk
COPY _internal\_run_testautogk.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testautogk.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_autogk_short
COPY _internal\_run_short_testautogk.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testautogk.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_autogk_silent
COPY _internal\_run_silent_testautogk.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testautogk.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_gq
COPY _internal\_run_testgq.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testgq.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_gq_short
COPY _internal\_run_short_testgq.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testgq.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_gq_silent
COPY _internal\_run_silent_testgq.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testgq.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_gkq
COPY _internal\_run_testgkq.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testgkq.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_gkq_short
COPY _internal\_run_short_testgkq.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testgkq.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_gkq_silent
COPY _internal\_run_silent_testgkq.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testgkq.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_lsfit
COPY _internal\_run_llstestunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\llstestunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_lsfit_short
COPY _internal\_run_short_llstestunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\llstestunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_lsfit_silent
COPY _internal\_run_silent_llstestunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\llstestunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_minlm
COPY _internal\_run_testlm.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testlm.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_minlm_short
COPY _internal\_run_short_testlm.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testlm.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_minlm_silent
COPY _internal\_run_silent_testlm.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testlm.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_polint
COPY _internal\_run_testpolintunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testpolintunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_polint_short
COPY _internal\_run_short_testpolintunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testpolintunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_polint_silent
COPY _internal\_run_silent_testpolintunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testpolintunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_ratinterpolation
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_ratinterpolation_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_ratinterpolation_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_ratint
COPY _internal\_run_testratinterpolation.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testratinterpolation.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_ratint_short
COPY _internal\_run_short_testratinterpolation.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testratinterpolation.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_ratint_silent
COPY _internal\_run_silent_testratinterpolation.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testratinterpolation.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_taskgen
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_taskgen_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_taskgen_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_spline2d
COPY _internal\_run_testinterpolation2dunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testinterpolation2dunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_spline2d_short
COPY _internal\_run_short_testinterpolation2dunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testinterpolation2dunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_spline2d_silent
COPY _internal\_run_silent_testinterpolation2dunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testinterpolation2dunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_spline3
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_spline3_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_spline3_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_spline1d
COPY _internal\_run_testsplineinterpolationunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testsplineinterpolationunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_spline1d_short
COPY _internal\_run_short_testsplineinterpolationunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testsplineinterpolationunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_spline1d_silent
COPY _internal\_run_silent_testsplineinterpolationunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testsplineinterpolationunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_idwint
COPY _internal\_run_testidwunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testidwunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_idwint_short
COPY _internal\_run_short_testidwunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testidwunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_idwint_silent
COPY _internal\_run_silent_testidwunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testidwunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_nearestneighbor
COPY _internal\_run_testnearestneighborunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testnearestneighborunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_nearestneighbor_short
COPY _internal\_run_short_testnearestneighborunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testnearestneighborunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_nearestneighbor_silent
COPY _internal\_run_silent_testnearestneighborunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testnearestneighborunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_matdet
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_matdet_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_matdet_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_sdet
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_sdet_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_sdet_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_ldlt
COPY _internal\_run_testldltunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testldltunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_ldlt_short
COPY _internal\_run_short_testldltunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testldltunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_ldlt_silent
COPY _internal\_run_silent_testldltunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testldltunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_spdgevd
COPY _internal\_run_testspdgevdunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testspdgevdunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_spdgevd_short
COPY _internal\_run_short_testspdgevdunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testspdgevdunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_spdgevd_silent
COPY _internal\_run_silent_testspdgevdunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testspdgevdunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_sinverse
COPY _internal\_run_testinvldltunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testinvldltunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_sinverse_short
COPY _internal\_run_short_testinvldltunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testinvldltunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_sinverse_silent
COPY _internal\_run_silent_testinvldltunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testinvldltunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_inverseupdate
COPY _internal\_run_testinverseupdateunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testinverseupdateunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_inverseupdate_short
COPY _internal\_run_short_testinverseupdateunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testinverseupdateunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_inverseupdate_silent
COPY _internal\_run_silent_testinverseupdateunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testinverseupdateunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_srcond
COPY _internal\_run_testrcondldltunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testrcondldltunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_srcond_short
COPY _internal\_run_short_testrcondldltunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testrcondldltunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_srcond_silent
COPY _internal\_run_silent_testrcondldltunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testrcondldltunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_ssolve
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_ssolve_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_ssolve_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_estnorm
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_estnorm_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_estnorm_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_schur
COPY _internal\_run_testschurunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testschurunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_schur_short
COPY _internal\_run_short_testschurunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testschurunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_schur_silent
COPY _internal\_run_silent_testschurunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testschurunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_airyf
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_airyf_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_airyf_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_bessel
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_bessel_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_bessel_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_betaf
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_betaf_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_betaf_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_chebyshev
COPY _internal\_run_testchebyshevunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testchebyshevunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_chebyshev_short
COPY _internal\_run_short_testchebyshevunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testchebyshevunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_chebyshev_silent
COPY _internal\_run_silent_testchebyshevunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testchebyshevunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_dawson
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_dawson_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_dawson_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_elliptic
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_elliptic_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_elliptic_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_expintegrals
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_expintegrals_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_expintegrals_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_fresnel
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_fresnel_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_fresnel_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_hermite
COPY _internal\_run_testhermiteunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testhermiteunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_hermite_short
COPY _internal\_run_short_testhermiteunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testhermiteunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_hermite_silent
COPY _internal\_run_silent_testhermiteunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testhermiteunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_ibetaf
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_ibetaf_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_ibetaf_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_jacobianelliptic
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_jacobianelliptic_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_jacobianelliptic_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_laguerre
COPY _internal\_run_testlaguerreunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testlaguerreunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_laguerre_short
COPY _internal\_run_short_testlaguerreunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testlaguerreunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_laguerre_silent
COPY _internal\_run_silent_testlaguerreunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testlaguerreunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_legendre
COPY _internal\_run_testlegendreunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testlegendreunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_legendre_short
COPY _internal\_run_short_testlegendreunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testlegendreunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_legendre_silent
COPY _internal\_run_silent_testlegendreunit.pas _tmp\_test.pas > NUL 2> NUL
COPY tests\testlegendreunit.* _tmp > NUL 2> NUL
GOTO lbl_201___end
:lbl_201_psif
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_psif_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_psif_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_trigintegrals
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_trigintegrals_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_trigintegrals_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_binomialdistr
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_binomialdistr_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_binomialdistr_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_nearunityunit
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_nearunityunit_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_nearunityunit_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_chisquaredistr
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_chisquaredistr_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_chisquaredistr_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_correlation
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_correlation_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_correlation_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_fdistr
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_fdistr_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_fdistr_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_correlationtests
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_correlationtests_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_correlationtests_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_studenttdistr
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_studenttdistr_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_studenttdistr_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_jarquebera
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_jarquebera_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_jarquebera_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_mannwhitneyu
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_mannwhitneyu_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_mannwhitneyu_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_poissondistr
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_poissondistr_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_poissondistr_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_stest
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_stest_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_stest_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_studentttests
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_studentttests_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_studentttests_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_variancetests
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_variancetests_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_variancetests_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201_wsr
echo No separate test file for this unit
echo It is tested somewhere else
EXIT /B 0
GOTO lbl_201___end
:lbl_201_wsr_short
EXIT /B 0
GOTO lbl_201___end
:lbl_201_wsr_silent
EXIT /B 0
GOTO lbl_201___end
:lbl_201___error
echo unknown unit
EXIT /B 1
:lbl_201___end
IF "%~1"=="dcc32" GOTO lbl_202_dcc32
GOTO lbl_202___error
:lbl_202_dcc32
CHDIR _tmp
dcc32 -$C+ -CC -$Q+ -$R+ -$X+ %~3 _test.pas >> ../log.txt 2>&1 2>&1
IF NOT ERRORLEVEL 1 GOTO lbl_203
echo Error while compiling (see ../log.txt for more info)
CHDIR ..
EXIT /B 1
:lbl_203
CHDIR ..
pushd _tmp
.\_test
popd
GOTO lbl_202___end
:lbl_202___error
echo unknown compiler
EXIT /B 1
:lbl_202___end
EXIT /B 0
