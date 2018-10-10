(*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
      pseudocode.

See subroutines comments for additional copyrights.

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************)
unit nsevd;
interface
uses Math, Sysutils, Ap, blas, reflections, rotations, hsschur, hessenberg;

function RMatrixEVD(A : TReal2DArray;
     N : AlglibInteger;
     VNeeded : AlglibInteger;
     var WR : TReal1DArray;
     var WI : TReal1DArray;
     var VL : TReal2DArray;
     var VR : TReal2DArray):Boolean;
function NonSymmetricEVD(A : TReal2DArray;
     N : AlglibInteger;
     VNeeded : AlglibInteger;
     var WR : TReal1DArray;
     var WI : TReal1DArray;
     var VL : TReal2DArray;
     var VR : TReal2DArray):Boolean;

implementation

procedure InternalTREVC(const T : TReal2DArray;
     N : AlglibInteger;
     SIDE : AlglibInteger;
     HOWMNY : AlglibInteger;
     VSELECT : TBoolean1DArray;
     var VL : TReal2DArray;
     var VR : TReal2DArray;
     var M : AlglibInteger;
     var INFO : AlglibInteger);forward;
procedure InternalHSEVDLALN2(const LTRANS : Boolean;
     const NA : AlglibInteger;
     const NW : AlglibInteger;
     const SMIN : Double;
     const CA : Double;
     const A : TReal2DArray;
     const D1 : Double;
     const D2 : Double;
     const B : TReal2DArray;
     const WR : Double;
     const WI : Double;
     var RSWAP4 : TBoolean1DArray;
     var ZSWAP4 : TBoolean1DArray;
     var IPIVOT44 : TInteger2DArray;
     var CIV4 : TReal1DArray;
     var CRV4 : TReal1DArray;
     var X : TReal2DArray;
     var SCL : Double;
     var XNORM : Double;
     var INFO : AlglibInteger);forward;
procedure InternalHSEVDLADIV(const A : Double;
     const B : Double;
     const C : Double;
     const D : Double;
     var P : Double;
     var Q : Double);forward;


(*************************************************************************
Finding eigenvalues and eigenvectors of a general matrix

The algorithm finds eigenvalues and eigenvectors of a general matrix by
using the QR algorithm with multiple shifts. The algorithm can find
eigenvalues and both left and right eigenvectors.

The right eigenvector is a vector x such that A*x = w*x, and the left
eigenvector is a vector y such that y'*A = w*y' (here y' implies a complex
conjugate transposition of vector y).

Input parameters:
    A       -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N       -   size of matrix A.
    VNeeded -   flag controlling whether eigenvectors are needed or not.
                If VNeeded is equal to:
                 * 0, eigenvectors are not returned;
                 * 1, right eigenvectors are returned;
                 * 2, left eigenvectors are returned;
                 * 3, both left and right eigenvectors are returned.

Output parameters:
    WR      -   real parts of eigenvalues.
                Array whose index ranges within [0..N-1].
    WR      -   imaginary parts of eigenvalues.
                Array whose index ranges within [0..N-1].
    VL, VR  -   arrays of left and right eigenvectors (if they are needed).
                If WI[i]=0, the respective eigenvalue is a real number,
                and it corresponds to the column number I of matrices VL/VR.
                If WI[i]>0, we have a pair of complex conjugate numbers with
                positive and negative imaginary parts:
                    the first eigenvalue WR[i] + sqrt(-1)*WI[i];
                    the second eigenvalue WR[i+1] + sqrt(-1)*WI[i+1];
                    WI[i]>0
                    WI[i+1] = -WI[i] < 0
                In that case, the eigenvector  corresponding to the first
                eigenvalue is located in i and i+1 columns of matrices
                VL/VR (the column number i contains the real part, and the
                column number i+1 contains the imaginary part), and the vector
                corresponding to the second eigenvalue is a complex conjugate to
                the first vector.
                Arrays whose indexes range within [0..N-1, 0..N-1].

Result:
    True, if the algorithm has converged.
    False, if the algorithm has not converged.

Note 1:
    Some users may ask the following question: what if WI[N-1]>0?
    WI[N] must contain an eigenvalue which is complex conjugate to the
    N-th eigenvalue, but the array has only size N?
    The answer is as follows: such a situation cannot occur because the
    algorithm finds a pairs of eigenvalues, therefore, if WI[i]>0, I is
    strictly less than N-1.

Note 2:
    The algorithm performance depends on the value of the internal parameter
    NS of the InternalSchurDecomposition subroutine which defines the number
    of shifts in the QR algorithm (similarly to the block width in block-matrix
    algorithms of linear algebra). If you require maximum performance
    on your machine, it is recommended to adjust this parameter manually.


See also the InternalTREVC subroutine.

The algorithm is based on the LAPACK 3.0 library.
*************************************************************************)
function RMatrixEVD(A : TReal2DArray;
     N : AlglibInteger;
     VNeeded : AlglibInteger;
     var WR : TReal1DArray;
     var WI : TReal1DArray;
     var VL : TReal2DArray;
     var VR : TReal2DArray):Boolean;
var
    A1 : TReal2DArray;
    VL1 : TReal2DArray;
    VR1 : TReal2DArray;
    WR1 : TReal1DArray;
    WI1 : TReal1DArray;
    I : AlglibInteger;
    MX : Double;
begin
    A := DynamicArrayCopy(A);
    Assert((VNeeded>=0) and (VNeeded<=3), 'RMatrixEVD: incorrect VNeeded!');
    SetLength(A1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        APVMove(@A1[I][0], 1, N, @A[I-1][0], 0, N-1);
        Inc(I);
    end;
    Result := NonSymmetricEVD(A1, N, VNeeded, WR1, WI1, VL1, VR1);
    if Result then
    begin
        SetLength(WR, N-1+1);
        SetLength(WI, N-1+1);
        APVMove(@WR[0], 0, N-1, @WR1[0], 1, N);
        APVMove(@WI[0], 0, N-1, @WI1[0], 1, N);
        if (VNeeded=2) or (VNeeded=3) then
        begin
            SetLength(VL, N-1+1, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                APVMove(@VL[I][0], 0, N-1, @VL1[I+1][0], 1, N);
                Inc(I);
            end;
        end;
        if (VNeeded=1) or (VNeeded=3) then
        begin
            SetLength(VR, N-1+1, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                APVMove(@VR[I][0], 0, N-1, @VR1[I+1][0], 1, N);
                Inc(I);
            end;
        end;
    end;
end;


function NonSymmetricEVD(A : TReal2DArray;
     N : AlglibInteger;
     VNeeded : AlglibInteger;
     var WR : TReal1DArray;
     var WI : TReal1DArray;
     var VL : TReal2DArray;
     var VR : TReal2DArray):Boolean;
var
    S : TReal2DArray;
    Tau : TReal1DArray;
    SEL : TBoolean1DArray;
    I : AlglibInteger;
    INFO : AlglibInteger;
    M : AlglibInteger;
begin
    A := DynamicArrayCopy(A);
    Assert((VNeeded>=0) and (VNeeded<=3), 'NonSymmetricEVD: incorrect VNeeded!');
    if VNeeded=0 then
    begin
        
        //
        // Eigen values only
        //
        ToUpperHessenberg(A, N, Tau);
        InternalSchurDecomposition(A, N, 0, 0, WR, WI, S, INFO);
        Result := INFO=0;
        Exit;
    end;
    
    //
    // Eigen values and vectors
    //
    ToUpperHessenberg(A, N, Tau);
    UnpackQFromUpperHessenberg(A, N, Tau, S);
    InternalSchurDecomposition(A, N, 1, 1, WR, WI, S, INFO);
    Result := INFO=0;
    if  not Result then
    begin
        Exit;
    end;
    if (VNeeded=1) or (VNeeded=3) then
    begin
        SetLength(VR, N+1, N+1);
        I:=1;
        while I<=N do
        begin
            APVMove(@VR[I][0], 1, N, @S[I][0], 1, N);
            Inc(I);
        end;
    end;
    if (VNeeded=2) or (VNeeded=3) then
    begin
        SetLength(VL, N+1, N+1);
        I:=1;
        while I<=N do
        begin
            APVMove(@VL[I][0], 1, N, @S[I][0], 1, N);
            Inc(I);
        end;
    end;
    InternalTREVC(A, N, VNeeded, 1, SEL, VL, VR, M, INFO);
    Result := INFO=0;
end;


procedure InternalTREVC(const T : TReal2DArray;
     N : AlglibInteger;
     SIDE : AlglibInteger;
     HOWMNY : AlglibInteger;
     VSELECT : TBoolean1DArray;
     var VL : TReal2DArray;
     var VR : TReal2DArray;
     var M : AlglibInteger;
     var INFO : AlglibInteger);
var
    ALLV : Boolean;
    BOTHV : Boolean;
    LEFTV : Boolean;
    OVER : Boolean;
    PAIR : Boolean;
    RIGHTV : Boolean;
    SOMEV : Boolean;
    I : AlglibInteger;
    IERR : AlglibInteger;
    II : AlglibInteger;
    IP : AlglibInteger;
    IIS : AlglibInteger;
    J : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
    JNXT : AlglibInteger;
    K : AlglibInteger;
    KI : AlglibInteger;
    N2 : AlglibInteger;
    BETA : Double;
    BIGNUM : Double;
    EMAX : Double;
    OVFL : Double;
    REC : Double;
    REMAX : Double;
    SCL : Double;
    SMIN : Double;
    SMLNUM : Double;
    ULP : Double;
    UNFL : Double;
    VCRIT : Double;
    VMAX : Double;
    WI : Double;
    WR : Double;
    XNORM : Double;
    X : TReal2DArray;
    WORK : TReal1DArray;
    TEMP : TReal1DArray;
    TEMP11 : TReal2DArray;
    TEMP22 : TReal2DArray;
    TEMP11B : TReal2DArray;
    TEMP21B : TReal2DArray;
    TEMP12B : TReal2DArray;
    TEMP22B : TReal2DArray;
    SkipFlag : Boolean;
    K1 : AlglibInteger;
    K2 : AlglibInteger;
    K3 : AlglibInteger;
    K4 : AlglibInteger;
    VT : Double;
    RSWAP4 : TBoolean1DArray;
    ZSWAP4 : TBoolean1DArray;
    IPIVOT44 : TInteger2DArray;
    CIV4 : TReal1DArray;
    CRV4 : TReal1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    VSELECT := DynamicArrayCopy(VSELECT);
    SetLength(X, 2+1, 2+1);
    SetLength(TEMP11, 1+1, 1+1);
    SetLength(TEMP11B, 1+1, 1+1);
    SetLength(TEMP21B, 2+1, 1+1);
    SetLength(TEMP12B, 1+1, 2+1);
    SetLength(TEMP22B, 2+1, 2+1);
    SetLength(TEMP22, 2+1, 2+1);
    SetLength(WORK, 3*N+1);
    SetLength(TEMP, N+1);
    SetLength(RSWAP4, 4+1);
    SetLength(ZSWAP4, 4+1);
    SetLength(IPIVOT44, 4+1, 4+1);
    SetLength(CIV4, 4+1);
    SetLength(CRV4, 4+1);
    if HOWMNY<>1 then
    begin
        if (SIDE=1) or (SIDE=3) then
        begin
            SetLength(VR, N+1, N+1);
        end;
        if (SIDE=2) or (SIDE=3) then
        begin
            SetLength(VL, N+1, N+1);
        end;
    end;
    
    //
    // Decode and test the input parameters
    //
    BOTHV := SIDE=3;
    RIGHTV := (SIDE=1) or BOTHV;
    LEFTV := (SIDE=2) or BOTHV;
    ALLV := HOWMNY=2;
    OVER := HOWMNY=1;
    SOMEV := HOWMNY=3;
    INFO := 0;
    if N<0 then
    begin
        INFO := -2;
        Exit;
    end;
    if  not RIGHTV and  not LEFTV then
    begin
        INFO := -3;
        Exit;
    end;
    if  not ALLV and  not OVER and  not SOMEV then
    begin
        INFO := -4;
        Exit;
    end;
    
    //
    // Set M to the number of columns required to store the selected
    // eigenvectors, standardize the array SELECT if necessary, and
    // test MM.
    //
    if SOMEV then
    begin
        M := 0;
        PAIR := False;
        J:=1;
        while J<=N do
        begin
            if PAIR then
            begin
                PAIR := False;
                VSELECT[J] := False;
            end
            else
            begin
                if J<N then
                begin
                    if AP_FP_Eq(T[J+1,J],0) then
                    begin
                        if VSELECT[J] then
                        begin
                            M := M+1;
                        end;
                    end
                    else
                    begin
                        PAIR := True;
                        if VSELECT[J] or VSELECT[J+1] then
                        begin
                            VSELECT[J] := True;
                            M := M+2;
                        end;
                    end;
                end
                else
                begin
                    if VSELECT[N] then
                    begin
                        M := M+1;
                    end;
                end;
            end;
            Inc(J);
        end;
    end
    else
    begin
        M := N;
    end;
    
    //
    // Quick return if possible.
    //
    if N=0 then
    begin
        Exit;
    end;
    
    //
    // Set the constants to control overflow.
    //
    UNFL := MinRealNumber;
    OVFL := 1/UNFL;
    ULP := MachineEpsilon;
    SMLNUM := UNFL*(N/ULP);
    BIGNUM := (1-ULP)/SMLNUM;
    
    //
    // Compute 1-norm of each column of strictly upper triangular
    // part of T to control overflow in triangular solver.
    //
    WORK[1] := 0;
    J:=2;
    while J<=N do
    begin
        WORK[J] := 0;
        I:=1;
        while I<=J-1 do
        begin
            WORK[J] := WORK[J]+AbsReal(T[I,J]);
            Inc(I);
        end;
        Inc(J);
    end;
    
    //
    // Index IP is used to specify the real or complex eigenvalue:
    // IP = 0, real eigenvalue,
    //      1, first of conjugate complex pair: (wr,wi)
    //     -1, second of conjugate complex pair: (wr,wi)
    //
    N2 := 2*N;
    if RIGHTV then
    begin
        
        //
        // Compute right eigenvectors.
        //
        IP := 0;
        IIS := M;
        KI:=N;
        while KI>=1 do
        begin
            SkipFlag := False;
            if IP=1 then
            begin
                SkipFlag := True;
            end
            else
            begin
                if KI<>1 then
                begin
                    if AP_FP_Neq(T[KI,KI-1],0) then
                    begin
                        IP := -1;
                    end;
                end;
                if SOMEV then
                begin
                    if IP=0 then
                    begin
                        if  not VSELECT[KI] then
                        begin
                            SkipFlag := True;
                        end;
                    end
                    else
                    begin
                        if  not VSELECT[KI-1] then
                        begin
                            SkipFlag := True;
                        end;
                    end;
                end;
            end;
            if  not SkipFlag then
            begin
                
                //
                // Compute the KI-th eigenvalue (WR,WI).
                //
                WR := T[KI,KI];
                WI := 0;
                if IP<>0 then
                begin
                    WI := SQRT(AbsReal(T[KI,KI-1]))*SQRT(AbsReal(T[KI-1,KI]));
                end;
                SMIN := Max(ULP*(AbsReal(WR)+AbsReal(WI)), SMLNUM);
                if IP=0 then
                begin
                    
                    //
                    // Real right eigenvector
                    //
                    WORK[KI+N] := 1;
                    
                    //
                    // Form right-hand side
                    //
                    K:=1;
                    while K<=KI-1 do
                    begin
                        WORK[K+N] := -T[K,KI];
                        Inc(K);
                    end;
                    
                    //
                    // Solve the upper quasi-triangular system:
                    //   (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK.
                    //
                    JNXT := KI-1;
                    J:=KI-1;
                    while J>=1 do
                    begin
                        if J>JNXT then
                        begin
                            Dec(J);
                            Continue;
                        end;
                        J1 := J;
                        J2 := J;
                        JNXT := J-1;
                        if J>1 then
                        begin
                            if AP_FP_Neq(T[J,J-1],0) then
                            begin
                                J1 := J-1;
                                JNXT := J-2;
                            end;
                        end;
                        if J1=J2 then
                        begin
                            
                            //
                            // 1-by-1 diagonal block
                            //
                            TEMP11[1,1] := T[J,J];
                            TEMP11B[1,1] := WORK[J+N];
                            InternalHSEVDLALN2(False, 1, 1, SMIN, 1, TEMP11, 1.0, 1.0, TEMP11B, WR, 0.0, RSWAP4, ZSWAP4, IPIVOT44, CIV4, CRV4, X, SCL, XNORM, IERR);
                            
                            //
                            // Scale X(1,1) to avoid overflow when updating
                            // the right-hand side.
                            //
                            if AP_FP_Greater(XNORM,1) then
                            begin
                                if AP_FP_Greater(WORK[J],BIGNUM/XNORM) then
                                begin
                                    X[1,1] := X[1,1]/XNORM;
                                    SCL := SCL/XNORM;
                                end;
                            end;
                            
                            //
                            // Scale if necessary
                            //
                            if AP_FP_Neq(SCL,1) then
                            begin
                                K1 := N+1;
                                K2 := N+KI;
                                APVMul(@WORK[0], K1, K2, SCL);
                            end;
                            WORK[J+N] := X[1,1];
                            
                            //
                            // Update right-hand side
                            //
                            K1 := 1+N;
                            K2 := J-1+N;
                            K3 := J-1;
                            VT := -X[1,1];
                            i1_ := (1) - (K1);
                            for i_ := K1 to K2 do
                            begin
                                WORK[i_] := WORK[i_] + VT*T[i_+i1_,J];
                            end;
                        end
                        else
                        begin
                            
                            //
                            // 2-by-2 diagonal block
                            //
                            TEMP22[1,1] := T[J-1,J-1];
                            TEMP22[1,2] := T[J-1,J];
                            TEMP22[2,1] := T[J,J-1];
                            TEMP22[2,2] := T[J,J];
                            TEMP21B[1,1] := WORK[J-1+N];
                            TEMP21B[2,1] := WORK[J+N];
                            InternalHSEVDLALN2(False, 2, 1, SMIN, 1.0, TEMP22, 1.0, 1.0, TEMP21B, WR, 0, RSWAP4, ZSWAP4, IPIVOT44, CIV4, CRV4, X, SCL, XNORM, IERR);
                            
                            //
                            // Scale X(1,1) and X(2,1) to avoid overflow when
                            // updating the right-hand side.
                            //
                            if AP_FP_Greater(XNORM,1) then
                            begin
                                BETA := Max(WORK[J-1], WORK[J]);
                                if AP_FP_Greater(BETA,BIGNUM/XNORM) then
                                begin
                                    X[1,1] := X[1,1]/XNORM;
                                    X[2,1] := X[2,1]/XNORM;
                                    SCL := SCL/XNORM;
                                end;
                            end;
                            
                            //
                            // Scale if necessary
                            //
                            if AP_FP_Neq(SCL,1) then
                            begin
                                K1 := 1+N;
                                K2 := KI+N;
                                APVMul(@WORK[0], K1, K2, SCL);
                            end;
                            WORK[J-1+N] := X[1,1];
                            WORK[J+N] := X[2,1];
                            
                            //
                            // Update right-hand side
                            //
                            K1 := 1+N;
                            K2 := J-2+N;
                            K3 := J-2;
                            K4 := J-1;
                            VT := -X[1,1];
                            i1_ := (1) - (K1);
                            for i_ := K1 to K2 do
                            begin
                                WORK[i_] := WORK[i_] + VT*T[i_+i1_,K4];
                            end;
                            VT := -X[2,1];
                            i1_ := (1) - (K1);
                            for i_ := K1 to K2 do
                            begin
                                WORK[i_] := WORK[i_] + VT*T[i_+i1_,J];
                            end;
                        end;
                        Dec(J);
                    end;
                    
                    //
                    // Copy the vector x or Q*x to VR and normalize.
                    //
                    if  not OVER then
                    begin
                        K1 := 1+N;
                        K2 := KI+N;
                        i1_ := (K1) - (1);
                        for i_ := 1 to KI do
                        begin
                            VR[i_,IIS] := WORK[i_+i1_];
                        end;
                        II := ColumnIdxAbsMax(VR, 1, KI, IIS);
                        REMAX := 1/AbsReal(VR[II,IIS]);
                        for i_ := 1 to KI do
                        begin
                            VR[i_,IIS] := REMAX*VR[i_,IIS];
                        end;
                        K:=KI+1;
                        while K<=N do
                        begin
                            VR[K,IIS] := 0;
                            Inc(K);
                        end;
                    end
                    else
                    begin
                        if KI>1 then
                        begin
                            for i_ := 1 to N do
                            begin
                                TEMP[i_] := VR[i_,KI];
                            end;
                            MatrixVectorMultiply(VR, 1, N, 1, KI-1, False, WORK, 1+N, KI-1+N, 1.0, TEMP, 1, N, WORK[KI+N]);
                            for i_ := 1 to N do
                            begin
                                VR[i_,KI] := TEMP[i_];
                            end;
                        end;
                        II := ColumnIdxAbsMax(VR, 1, N, KI);
                        REMAX := 1/AbsReal(VR[II,KI]);
                        for i_ := 1 to N do
                        begin
                            VR[i_,KI] := REMAX*VR[i_,KI];
                        end;
                    end;
                end
                else
                begin
                    
                    //
                    // Complex right eigenvector.
                    //
                    // Initial solve
                    //     [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
                    //     [ (T(KI,KI-1)   T(KI,KI)   )               ]
                    //
                    if AP_FP_Greater_Eq(AbsReal(T[KI-1,KI]),AbsReal(T[KI,KI-1])) then
                    begin
                        WORK[KI-1+N] := 1;
                        WORK[KI+N2] := WI/T[KI-1,KI];
                    end
                    else
                    begin
                        WORK[KI-1+N] := -WI/T[KI,KI-1];
                        WORK[KI+N2] := 1;
                    end;
                    WORK[KI+N] := 0;
                    WORK[KI-1+N2] := 0;
                    
                    //
                    // Form right-hand side
                    //
                    K:=1;
                    while K<=KI-2 do
                    begin
                        WORK[K+N] := -WORK[KI-1+N]*T[K,KI-1];
                        WORK[K+N2] := -WORK[KI+N2]*T[K,KI];
                        Inc(K);
                    end;
                    
                    //
                    // Solve upper quasi-triangular system:
                    // (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2)
                    //
                    JNXT := KI-2;
                    J:=KI-2;
                    while J>=1 do
                    begin
                        if J>JNXT then
                        begin
                            Dec(J);
                            Continue;
                        end;
                        J1 := J;
                        J2 := J;
                        JNXT := J-1;
                        if J>1 then
                        begin
                            if AP_FP_Neq(T[J,J-1],0) then
                            begin
                                J1 := J-1;
                                JNXT := J-2;
                            end;
                        end;
                        if J1=J2 then
                        begin
                            
                            //
                            // 1-by-1 diagonal block
                            //
                            TEMP11[1,1] := T[J,J];
                            TEMP12B[1,1] := WORK[J+N];
                            TEMP12B[1,2] := WORK[J+N+N];
                            InternalHSEVDLALN2(False, 1, 2, SMIN, 1.0, TEMP11, 1.0, 1.0, TEMP12B, WR, WI, RSWAP4, ZSWAP4, IPIVOT44, CIV4, CRV4, X, SCL, XNORM, IERR);
                            
                            //
                            // Scale X(1,1) and X(1,2) to avoid overflow when
                            // updating the right-hand side.
                            //
                            if AP_FP_Greater(XNORM,1) then
                            begin
                                if AP_FP_Greater(WORK[J],BIGNUM/XNORM) then
                                begin
                                    X[1,1] := X[1,1]/XNORM;
                                    X[1,2] := X[1,2]/XNORM;
                                    SCL := SCL/XNORM;
                                end;
                            end;
                            
                            //
                            // Scale if necessary
                            //
                            if AP_FP_Neq(SCL,1) then
                            begin
                                K1 := 1+N;
                                K2 := KI+N;
                                APVMul(@WORK[0], K1, K2, SCL);
                                K1 := 1+N2;
                                K2 := KI+N2;
                                APVMul(@WORK[0], K1, K2, SCL);
                            end;
                            WORK[J+N] := X[1,1];
                            WORK[J+N2] := X[1,2];
                            
                            //
                            // Update the right-hand side
                            //
                            K1 := 1+N;
                            K2 := J-1+N;
                            K3 := 1;
                            K4 := J-1;
                            VT := -X[1,1];
                            i1_ := (K3) - (K1);
                            for i_ := K1 to K2 do
                            begin
                                WORK[i_] := WORK[i_] + VT*T[i_+i1_,J];
                            end;
                            K1 := 1+N2;
                            K2 := J-1+N2;
                            K3 := 1;
                            K4 := J-1;
                            VT := -X[1,2];
                            i1_ := (K3) - (K1);
                            for i_ := K1 to K2 do
                            begin
                                WORK[i_] := WORK[i_] + VT*T[i_+i1_,J];
                            end;
                        end
                        else
                        begin
                            
                            //
                            // 2-by-2 diagonal block
                            //
                            TEMP22[1,1] := T[J-1,J-1];
                            TEMP22[1,2] := T[J-1,J];
                            TEMP22[2,1] := T[J,J-1];
                            TEMP22[2,2] := T[J,J];
                            TEMP22B[1,1] := WORK[J-1+N];
                            TEMP22B[1,2] := WORK[J-1+N+N];
                            TEMP22B[2,1] := WORK[J+N];
                            TEMP22B[2,2] := WORK[J+N+N];
                            InternalHSEVDLALN2(False, 2, 2, SMIN, 1.0, TEMP22, 1.0, 1.0, TEMP22B, WR, WI, RSWAP4, ZSWAP4, IPIVOT44, CIV4, CRV4, X, SCL, XNORM, IERR);
                            
                            //
                            // Scale X to avoid overflow when updating
                            // the right-hand side.
                            //
                            if AP_FP_Greater(XNORM,1) then
                            begin
                                BETA := Max(WORK[J-1], WORK[J]);
                                if AP_FP_Greater(BETA,BIGNUM/XNORM) then
                                begin
                                    REC := 1/XNORM;
                                    X[1,1] := X[1,1]*REC;
                                    X[1,2] := X[1,2]*REC;
                                    X[2,1] := X[2,1]*REC;
                                    X[2,2] := X[2,2]*REC;
                                    SCL := SCL*REC;
                                end;
                            end;
                            
                            //
                            // Scale if necessary
                            //
                            if AP_FP_Neq(SCL,1) then
                            begin
                                APVMul(@WORK[0], 1+N, KI+N, SCL);
                                APVMul(@WORK[0], 1+N2, KI+N2, SCL);
                            end;
                            WORK[J-1+N] := X[1,1];
                            WORK[J+N] := X[2,1];
                            WORK[J-1+N2] := X[1,2];
                            WORK[J+N2] := X[2,2];
                            
                            //
                            // Update the right-hand side
                            //
                            VT := -X[1,1];
                            i1_ := (1) - (N+1);
                            for i_ := N+1 to N+J-2 do
                            begin
                                WORK[i_] := WORK[i_] + VT*T[i_+i1_,J-1];
                            end;
                            VT := -X[2,1];
                            i1_ := (1) - (N+1);
                            for i_ := N+1 to N+J-2 do
                            begin
                                WORK[i_] := WORK[i_] + VT*T[i_+i1_,J];
                            end;
                            VT := -X[1,2];
                            i1_ := (1) - (N2+1);
                            for i_ := N2+1 to N2+J-2 do
                            begin
                                WORK[i_] := WORK[i_] + VT*T[i_+i1_,J-1];
                            end;
                            VT := -X[2,2];
                            i1_ := (1) - (N2+1);
                            for i_ := N2+1 to N2+J-2 do
                            begin
                                WORK[i_] := WORK[i_] + VT*T[i_+i1_,J];
                            end;
                        end;
                        Dec(J);
                    end;
                    
                    //
                    // Copy the vector x or Q*x to VR and normalize.
                    //
                    if  not OVER then
                    begin
                        i1_ := (N+1) - (1);
                        for i_ := 1 to KI do
                        begin
                            VR[i_,IIS-1] := WORK[i_+i1_];
                        end;
                        i1_ := (N2+1) - (1);
                        for i_ := 1 to KI do
                        begin
                            VR[i_,IIS] := WORK[i_+i1_];
                        end;
                        EMAX := 0;
                        K:=1;
                        while K<=KI do
                        begin
                            EMAX := Max(EMAX, AbsReal(VR[K,IIS-1])+AbsReal(VR[K,IIS]));
                            Inc(K);
                        end;
                        REMAX := 1/EMAX;
                        for i_ := 1 to KI do
                        begin
                            VR[i_,IIS-1] := REMAX*VR[i_,IIS-1];
                        end;
                        for i_ := 1 to KI do
                        begin
                            VR[i_,IIS] := REMAX*VR[i_,IIS];
                        end;
                        K:=KI+1;
                        while K<=N do
                        begin
                            VR[K,IIS-1] := 0;
                            VR[K,IIS] := 0;
                            Inc(K);
                        end;
                    end
                    else
                    begin
                        if KI>2 then
                        begin
                            for i_ := 1 to N do
                            begin
                                TEMP[i_] := VR[i_,KI-1];
                            end;
                            MatrixVectorMultiply(VR, 1, N, 1, KI-2, False, WORK, 1+N, KI-2+N, 1.0, TEMP, 1, N, WORK[KI-1+N]);
                            for i_ := 1 to N do
                            begin
                                VR[i_,KI-1] := TEMP[i_];
                            end;
                            for i_ := 1 to N do
                            begin
                                TEMP[i_] := VR[i_,KI];
                            end;
                            MatrixVectorMultiply(VR, 1, N, 1, KI-2, False, WORK, 1+N2, KI-2+N2, 1.0, TEMP, 1, N, WORK[KI+N2]);
                            for i_ := 1 to N do
                            begin
                                VR[i_,KI] := TEMP[i_];
                            end;
                        end
                        else
                        begin
                            VT := WORK[KI-1+N];
                            for i_ := 1 to N do
                            begin
                                VR[i_,KI-1] := VT*VR[i_,KI-1];
                            end;
                            VT := WORK[KI+N2];
                            for i_ := 1 to N do
                            begin
                                VR[i_,KI] := VT*VR[i_,KI];
                            end;
                        end;
                        EMAX := 0;
                        K:=1;
                        while K<=N do
                        begin
                            EMAX := Max(EMAX, AbsReal(VR[K,KI-1])+AbsReal(VR[K,KI]));
                            Inc(K);
                        end;
                        REMAX := 1/EMAX;
                        for i_ := 1 to N do
                        begin
                            VR[i_,KI-1] := REMAX*VR[i_,KI-1];
                        end;
                        for i_ := 1 to N do
                        begin
                            VR[i_,KI] := REMAX*VR[i_,KI];
                        end;
                    end;
                end;
                IIS := IIS-1;
                if IP<>0 then
                begin
                    IIS := IIS-1;
                end;
            end;
            if IP=1 then
            begin
                IP := 0;
            end;
            if IP=-1 then
            begin
                IP := 1;
            end;
            Dec(KI);
        end;
    end;
    if LEFTV then
    begin
        
        //
        // Compute left eigenvectors.
        //
        IP := 0;
        IIS := 1;
        KI:=1;
        while KI<=N do
        begin
            SkipFlag := False;
            if IP=-1 then
            begin
                SkipFlag := True;
            end
            else
            begin
                if KI<>N then
                begin
                    if AP_FP_Neq(T[KI+1,KI],0) then
                    begin
                        IP := 1;
                    end;
                end;
                if SOMEV then
                begin
                    if  not VSELECT[KI] then
                    begin
                        SkipFlag := True;
                    end;
                end;
            end;
            if  not SkipFlag then
            begin
                
                //
                // Compute the KI-th eigenvalue (WR,WI).
                //
                WR := T[KI,KI];
                WI := 0;
                if IP<>0 then
                begin
                    WI := SQRT(AbsReal(T[KI,KI+1]))*SQRT(AbsReal(T[KI+1,KI]));
                end;
                SMIN := Max(ULP*(AbsReal(WR)+AbsReal(WI)), SMLNUM);
                if IP=0 then
                begin
                    
                    //
                    // Real left eigenvector.
                    //
                    WORK[KI+N] := 1;
                    
                    //
                    // Form right-hand side
                    //
                    K:=KI+1;
                    while K<=N do
                    begin
                        WORK[K+N] := -T[KI,K];
                        Inc(K);
                    end;
                    
                    //
                    // Solve the quasi-triangular system:
                    // (T(KI+1:N,KI+1:N) - WR)'*X = SCALE*WORK
                    //
                    VMAX := 1;
                    VCRIT := BIGNUM;
                    JNXT := KI+1;
                    J:=KI+1;
                    while J<=N do
                    begin
                        if J<JNXT then
                        begin
                            Inc(J);
                            Continue;
                        end;
                        J1 := J;
                        J2 := J;
                        JNXT := J+1;
                        if J<N then
                        begin
                            if AP_FP_Neq(T[J+1,J],0) then
                            begin
                                J2 := J+1;
                                JNXT := J+2;
                            end;
                        end;
                        if J1=J2 then
                        begin
                            
                            //
                            // 1-by-1 diagonal block
                            //
                            // Scale if necessary to avoid overflow when forming
                            // the right-hand side.
                            //
                            if AP_FP_Greater(WORK[J],VCRIT) then
                            begin
                                REC := 1/VMAX;
                                APVMul(@WORK[0], KI+N, N+N, REC);
                                VMAX := 1;
                                VCRIT := BIGNUM;
                            end;
                            i1_ := (KI+1+N)-(KI+1);
                            VT := 0.0;
                            for i_ := KI+1 to J-1 do
                            begin
                                VT := VT + T[i_,J]*WORK[i_+i1_];
                            end;
                            WORK[J+N] := WORK[J+N]-VT;
                            
                            //
                            // Solve (T(J,J)-WR)'*X = WORK
                            //
                            TEMP11[1,1] := T[J,J];
                            TEMP11B[1,1] := WORK[J+N];
                            InternalHSEVDLALN2(False, 1, 1, SMIN, 1.0, TEMP11, 1.0, 1.0, TEMP11B, WR, 0, RSWAP4, ZSWAP4, IPIVOT44, CIV4, CRV4, X, SCL, XNORM, IERR);
                            
                            //
                            // Scale if necessary
                            //
                            if AP_FP_Neq(SCL,1) then
                            begin
                                APVMul(@WORK[0], KI+N, N+N, SCL);
                            end;
                            WORK[J+N] := X[1,1];
                            VMAX := Max(AbsReal(WORK[J+N]), VMAX);
                            VCRIT := BIGNUM/VMAX;
                        end
                        else
                        begin
                            
                            //
                            // 2-by-2 diagonal block
                            //
                            // Scale if necessary to avoid overflow when forming
                            // the right-hand side.
                            //
                            BETA := Max(WORK[J], WORK[J+1]);
                            if AP_FP_Greater(BETA,VCRIT) then
                            begin
                                REC := 1/VMAX;
                                APVMul(@WORK[0], KI+N, N+N, REC);
                                VMAX := 1;
                                VCRIT := BIGNUM;
                            end;
                            i1_ := (KI+1+N)-(KI+1);
                            VT := 0.0;
                            for i_ := KI+1 to J-1 do
                            begin
                                VT := VT + T[i_,J]*WORK[i_+i1_];
                            end;
                            WORK[J+N] := WORK[J+N]-VT;
                            i1_ := (KI+1+N)-(KI+1);
                            VT := 0.0;
                            for i_ := KI+1 to J-1 do
                            begin
                                VT := VT + T[i_,J+1]*WORK[i_+i1_];
                            end;
                            WORK[J+1+N] := WORK[J+1+N]-VT;
                            
                            //
                            // Solve
                            //    [T(J,J)-WR   T(J,J+1)     ]'* X = SCALE*( WORK1 )
                            //    [T(J+1,J)    T(J+1,J+1)-WR]             ( WORK2 )
                            //
                            TEMP22[1,1] := T[J,J];
                            TEMP22[1,2] := T[J,J+1];
                            TEMP22[2,1] := T[J+1,J];
                            TEMP22[2,2] := T[J+1,J+1];
                            TEMP21B[1,1] := WORK[J+N];
                            TEMP21B[2,1] := WORK[J+1+N];
                            InternalHSEVDLALN2(True, 2, 1, SMIN, 1.0, TEMP22, 1.0, 1.0, TEMP21B, WR, 0, RSWAP4, ZSWAP4, IPIVOT44, CIV4, CRV4, X, SCL, XNORM, IERR);
                            
                            //
                            // Scale if necessary
                            //
                            if AP_FP_Neq(SCL,1) then
                            begin
                                APVMul(@WORK[0], KI+N, N+N, SCL);
                            end;
                            WORK[J+N] := X[1,1];
                            WORK[J+1+N] := X[2,1];
                            VMAX := Max(AbsReal(WORK[J+N]), Max(AbsReal(WORK[J+1+N]), VMAX));
                            VCRIT := BIGNUM/VMAX;
                        end;
                        Inc(J);
                    end;
                    
                    //
                    // Copy the vector x or Q*x to VL and normalize.
                    //
                    if  not OVER then
                    begin
                        i1_ := (KI+N) - (KI);
                        for i_ := KI to N do
                        begin
                            VL[i_,IIS] := WORK[i_+i1_];
                        end;
                        II := ColumnIdxAbsMax(VL, KI, N, IIS);
                        REMAX := 1/AbsReal(VL[II,IIS]);
                        for i_ := KI to N do
                        begin
                            VL[i_,IIS] := REMAX*VL[i_,IIS];
                        end;
                        K:=1;
                        while K<=KI-1 do
                        begin
                            VL[K,IIS] := 0;
                            Inc(K);
                        end;
                    end
                    else
                    begin
                        if KI<N then
                        begin
                            for i_ := 1 to N do
                            begin
                                TEMP[i_] := VL[i_,KI];
                            end;
                            MatrixVectorMultiply(VL, 1, N, KI+1, N, False, WORK, KI+1+N, N+N, 1.0, TEMP, 1, N, WORK[KI+N]);
                            for i_ := 1 to N do
                            begin
                                VL[i_,KI] := TEMP[i_];
                            end;
                        end;
                        II := ColumnIdxAbsMax(VL, 1, N, KI);
                        REMAX := 1/AbsReal(VL[II,KI]);
                        for i_ := 1 to N do
                        begin
                            VL[i_,KI] := REMAX*VL[i_,KI];
                        end;
                    end;
                end
                else
                begin
                    
                    //
                    // Complex left eigenvector.
                    //
                    // Initial solve:
                    //   ((T(KI,KI)    T(KI,KI+1) )' - (WR - I* WI))*X = 0.
                    //   ((T(KI+1,KI) T(KI+1,KI+1))                )
                    //
                    if AP_FP_Greater_Eq(AbsReal(T[KI,KI+1]),AbsReal(T[KI+1,KI])) then
                    begin
                        WORK[KI+N] := WI/T[KI,KI+1];
                        WORK[KI+1+N2] := 1;
                    end
                    else
                    begin
                        WORK[KI+N] := 1;
                        WORK[KI+1+N2] := -WI/T[KI+1,KI];
                    end;
                    WORK[KI+1+N] := 0;
                    WORK[KI+N2] := 0;
                    
                    //
                    // Form right-hand side
                    //
                    K:=KI+2;
                    while K<=N do
                    begin
                        WORK[K+N] := -WORK[KI+N]*T[KI,K];
                        WORK[K+N2] := -WORK[KI+1+N2]*T[KI+1,K];
                        Inc(K);
                    end;
                    
                    //
                    // Solve complex quasi-triangular system:
                    // ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2
                    //
                    VMAX := 1;
                    VCRIT := BIGNUM;
                    JNXT := KI+2;
                    J:=KI+2;
                    while J<=N do
                    begin
                        if J<JNXT then
                        begin
                            Inc(J);
                            Continue;
                        end;
                        J1 := J;
                        J2 := J;
                        JNXT := J+1;
                        if J<N then
                        begin
                            if AP_FP_Neq(T[J+1,J],0) then
                            begin
                                J2 := J+1;
                                JNXT := J+2;
                            end;
                        end;
                        if J1=J2 then
                        begin
                            
                            //
                            // 1-by-1 diagonal block
                            //
                            // Scale if necessary to avoid overflow when
                            // forming the right-hand side elements.
                            //
                            if AP_FP_Greater(WORK[J],VCRIT) then
                            begin
                                REC := 1/VMAX;
                                APVMul(@WORK[0], KI+N, N+N, REC);
                                APVMul(@WORK[0], KI+N2, N+N2, REC);
                                VMAX := 1;
                                VCRIT := BIGNUM;
                            end;
                            i1_ := (KI+2+N)-(KI+2);
                            VT := 0.0;
                            for i_ := KI+2 to J-1 do
                            begin
                                VT := VT + T[i_,J]*WORK[i_+i1_];
                            end;
                            WORK[J+N] := WORK[J+N]-VT;
                            i1_ := (KI+2+N2)-(KI+2);
                            VT := 0.0;
                            for i_ := KI+2 to J-1 do
                            begin
                                VT := VT + T[i_,J]*WORK[i_+i1_];
                            end;
                            WORK[J+N2] := WORK[J+N2]-VT;
                            
                            //
                            // Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2
                            //
                            TEMP11[1,1] := T[J,J];
                            TEMP12B[1,1] := WORK[J+N];
                            TEMP12B[1,2] := WORK[J+N+N];
                            InternalHSEVDLALN2(False, 1, 2, SMIN, 1.0, TEMP11, 1.0, 1.0, TEMP12B, WR, -WI, RSWAP4, ZSWAP4, IPIVOT44, CIV4, CRV4, X, SCL, XNORM, IERR);
                            
                            //
                            // Scale if necessary
                            //
                            if AP_FP_Neq(SCL,1) then
                            begin
                                APVMul(@WORK[0], KI+N, N+N, SCL);
                                APVMul(@WORK[0], KI+N2, N+N2, SCL);
                            end;
                            WORK[J+N] := X[1,1];
                            WORK[J+N2] := X[1,2];
                            VMAX := Max(AbsReal(WORK[J+N]), Max(AbsReal(WORK[J+N2]), VMAX));
                            VCRIT := BIGNUM/VMAX;
                        end
                        else
                        begin
                            
                            //
                            // 2-by-2 diagonal block
                            //
                            // Scale if necessary to avoid overflow when forming
                            // the right-hand side elements.
                            //
                            BETA := Max(WORK[J], WORK[J+1]);
                            if AP_FP_Greater(BETA,VCRIT) then
                            begin
                                REC := 1/VMAX;
                                APVMul(@WORK[0], KI+N, N+N, REC);
                                APVMul(@WORK[0], KI+N2, N+N2, REC);
                                VMAX := 1;
                                VCRIT := BIGNUM;
                            end;
                            i1_ := (KI+2+N)-(KI+2);
                            VT := 0.0;
                            for i_ := KI+2 to J-1 do
                            begin
                                VT := VT + T[i_,J]*WORK[i_+i1_];
                            end;
                            WORK[J+N] := WORK[J+N]-VT;
                            i1_ := (KI+2+N2)-(KI+2);
                            VT := 0.0;
                            for i_ := KI+2 to J-1 do
                            begin
                                VT := VT + T[i_,J]*WORK[i_+i1_];
                            end;
                            WORK[J+N2] := WORK[J+N2]-VT;
                            i1_ := (KI+2+N)-(KI+2);
                            VT := 0.0;
                            for i_ := KI+2 to J-1 do
                            begin
                                VT := VT + T[i_,J+1]*WORK[i_+i1_];
                            end;
                            WORK[J+1+N] := WORK[J+1+N]-VT;
                            i1_ := (KI+2+N2)-(KI+2);
                            VT := 0.0;
                            for i_ := KI+2 to J-1 do
                            begin
                                VT := VT + T[i_,J+1]*WORK[i_+i1_];
                            end;
                            WORK[J+1+N2] := WORK[J+1+N2]-VT;
                            
                            //
                            // Solve 2-by-2 complex linear equation
                            //   ([T(j,j)   T(j,j+1)  ]'-(wr-i*wi)*I)*X = SCALE*B
                            //   ([T(j+1,j) T(j+1,j+1)]             )
                            //
                            TEMP22[1,1] := T[J,J];
                            TEMP22[1,2] := T[J,J+1];
                            TEMP22[2,1] := T[J+1,J];
                            TEMP22[2,2] := T[J+1,J+1];
                            TEMP22B[1,1] := WORK[J+N];
                            TEMP22B[1,2] := WORK[J+N+N];
                            TEMP22B[2,1] := WORK[J+1+N];
                            TEMP22B[2,2] := WORK[J+1+N+N];
                            InternalHSEVDLALN2(True, 2, 2, SMIN, 1.0, TEMP22, 1.0, 1.0, TEMP22B, WR, -WI, RSWAP4, ZSWAP4, IPIVOT44, CIV4, CRV4, X, SCL, XNORM, IERR);
                            
                            //
                            // Scale if necessary
                            //
                            if AP_FP_Neq(SCL,1) then
                            begin
                                APVMul(@WORK[0], KI+N, N+N, SCL);
                                APVMul(@WORK[0], KI+N2, N+N2, SCL);
                            end;
                            WORK[J+N] := X[1,1];
                            WORK[J+N2] := X[1,2];
                            WORK[J+1+N] := X[2,1];
                            WORK[J+1+N2] := X[2,2];
                            VMAX := Max(AbsReal(X[1,1]), VMAX);
                            VMAX := Max(AbsReal(X[1,2]), VMAX);
                            VMAX := Max(AbsReal(X[2,1]), VMAX);
                            VMAX := Max(AbsReal(X[2,2]), VMAX);
                            VCRIT := BIGNUM/VMAX;
                        end;
                        Inc(J);
                    end;
                    
                    //
                    // Copy the vector x or Q*x to VL and normalize.
                    //
                    if  not OVER then
                    begin
                        i1_ := (KI+N) - (KI);
                        for i_ := KI to N do
                        begin
                            VL[i_,IIS] := WORK[i_+i1_];
                        end;
                        i1_ := (KI+N2) - (KI);
                        for i_ := KI to N do
                        begin
                            VL[i_,IIS+1] := WORK[i_+i1_];
                        end;
                        EMAX := 0;
                        K:=KI;
                        while K<=N do
                        begin
                            EMAX := Max(EMAX, AbsReal(VL[K,IIS])+AbsReal(VL[K,IIS+1]));
                            Inc(K);
                        end;
                        REMAX := 1/EMAX;
                        for i_ := KI to N do
                        begin
                            VL[i_,IIS] := REMAX*VL[i_,IIS];
                        end;
                        for i_ := KI to N do
                        begin
                            VL[i_,IIS+1] := REMAX*VL[i_,IIS+1];
                        end;
                        K:=1;
                        while K<=KI-1 do
                        begin
                            VL[K,IIS] := 0;
                            VL[K,IIS+1] := 0;
                            Inc(K);
                        end;
                    end
                    else
                    begin
                        if KI<N-1 then
                        begin
                            for i_ := 1 to N do
                            begin
                                TEMP[i_] := VL[i_,KI];
                            end;
                            MatrixVectorMultiply(VL, 1, N, KI+2, N, False, WORK, KI+2+N, N+N, 1.0, TEMP, 1, N, WORK[KI+N]);
                            for i_ := 1 to N do
                            begin
                                VL[i_,KI] := TEMP[i_];
                            end;
                            for i_ := 1 to N do
                            begin
                                TEMP[i_] := VL[i_,KI+1];
                            end;
                            MatrixVectorMultiply(VL, 1, N, KI+2, N, False, WORK, KI+2+N2, N+N2, 1.0, TEMP, 1, N, WORK[KI+1+N2]);
                            for i_ := 1 to N do
                            begin
                                VL[i_,KI+1] := TEMP[i_];
                            end;
                        end
                        else
                        begin
                            VT := WORK[KI+N];
                            for i_ := 1 to N do
                            begin
                                VL[i_,KI] := VT*VL[i_,KI];
                            end;
                            VT := WORK[KI+1+N2];
                            for i_ := 1 to N do
                            begin
                                VL[i_,KI+1] := VT*VL[i_,KI+1];
                            end;
                        end;
                        EMAX := 0;
                        K:=1;
                        while K<=N do
                        begin
                            EMAX := Max(EMAX, AbsReal(VL[K,KI])+AbsReal(VL[K,KI+1]));
                            Inc(K);
                        end;
                        REMAX := 1/EMAX;
                        for i_ := 1 to N do
                        begin
                            VL[i_,KI] := REMAX*VL[i_,KI];
                        end;
                        for i_ := 1 to N do
                        begin
                            VL[i_,KI+1] := REMAX*VL[i_,KI+1];
                        end;
                    end;
                end;
                IIS := IIS+1;
                if IP<>0 then
                begin
                    IIS := IIS+1;
                end;
            end;
            if IP=-1 then
            begin
                IP := 0;
            end;
            if IP=1 then
            begin
                IP := -1;
            end;
            Inc(KI);
        end;
    end;
end;


procedure InternalHSEVDLALN2(const LTRANS : Boolean;
     const NA : AlglibInteger;
     const NW : AlglibInteger;
     const SMIN : Double;
     const CA : Double;
     const A : TReal2DArray;
     const D1 : Double;
     const D2 : Double;
     const B : TReal2DArray;
     const WR : Double;
     const WI : Double;
     var RSWAP4 : TBoolean1DArray;
     var ZSWAP4 : TBoolean1DArray;
     var IPIVOT44 : TInteger2DArray;
     var CIV4 : TReal1DArray;
     var CRV4 : TReal1DArray;
     var X : TReal2DArray;
     var SCL : Double;
     var XNORM : Double;
     var INFO : AlglibInteger);
var
    ICMAX : AlglibInteger;
    J : AlglibInteger;
    BBND : Double;
    BI1 : Double;
    BI2 : Double;
    BIGNUM : Double;
    BNORM : Double;
    BR1 : Double;
    BR2 : Double;
    CI21 : Double;
    CI22 : Double;
    CMAX : Double;
    CNORM : Double;
    CR21 : Double;
    CR22 : Double;
    CSI : Double;
    CSR : Double;
    LI21 : Double;
    LR21 : Double;
    SMINI : Double;
    SMLNUM : Double;
    TEMP : Double;
    U22ABS : Double;
    UI11 : Double;
    UI11R : Double;
    UI12 : Double;
    UI12S : Double;
    UI22 : Double;
    UR11 : Double;
    UR11R : Double;
    UR12 : Double;
    UR12S : Double;
    UR22 : Double;
    XI1 : Double;
    XI2 : Double;
    XR1 : Double;
    XR2 : Double;
    TMP1 : Double;
    TMP2 : Double;
begin
    ZSWAP4[1] := False;
    ZSWAP4[2] := False;
    ZSWAP4[3] := True;
    ZSWAP4[4] := True;
    RSWAP4[1] := False;
    RSWAP4[2] := True;
    RSWAP4[3] := False;
    RSWAP4[4] := True;
    IPIVOT44[1,1] := 1;
    IPIVOT44[2,1] := 2;
    IPIVOT44[3,1] := 3;
    IPIVOT44[4,1] := 4;
    IPIVOT44[1,2] := 2;
    IPIVOT44[2,2] := 1;
    IPIVOT44[3,2] := 4;
    IPIVOT44[4,2] := 3;
    IPIVOT44[1,3] := 3;
    IPIVOT44[2,3] := 4;
    IPIVOT44[3,3] := 1;
    IPIVOT44[4,3] := 2;
    IPIVOT44[1,4] := 4;
    IPIVOT44[2,4] := 3;
    IPIVOT44[3,4] := 2;
    IPIVOT44[4,4] := 1;
    SMLNUM := 2*MinRealNumber;
    BIGNUM := 1/SMLNUM;
    SMINI := Max(SMIN, SMLNUM);
    
    //
    // Don't check for input errors
    //
    INFO := 0;
    
    //
    // Standard Initializations
    //
    SCL := 1;
    if NA=1 then
    begin
        
        //
        // 1 x 1  (i.e., scalar) system   C X = B
        //
        if NW=1 then
        begin
            
            //
            // Real 1x1 system.
            //
            // C = ca A - w D
            //
            CSR := CA*A[1,1]-WR*D1;
            CNORM := AbsReal(CSR);
            
            //
            // If | C | < SMINI, use C = SMINI
            //
            if AP_FP_Less(CNORM,SMINI) then
            begin
                CSR := SMINI;
                CNORM := SMINI;
                INFO := 1;
            end;
            
            //
            // Check scaling for  X = B / C
            //
            BNORM := AbsReal(B[1,1]);
            if AP_FP_Less(CNORM,1) and AP_FP_Greater(BNORM,1) then
            begin
                if AP_FP_Greater(BNORM,BIGNUM*CNORM) then
                begin
                    SCL := 1/BNORM;
                end;
            end;
            
            //
            // Compute X
            //
            X[1,1] := B[1,1]*SCL/CSR;
            XNORM := AbsReal(X[1,1]);
        end
        else
        begin
            
            //
            // Complex 1x1 system (w is complex)
            //
            // C = ca A - w D
            //
            CSR := CA*A[1,1]-WR*D1;
            CSI := -WI*D1;
            CNORM := AbsReal(CSR)+AbsReal(CSI);
            
            //
            // If | C | < SMINI, use C = SMINI
            //
            if AP_FP_Less(CNORM,SMINI) then
            begin
                CSR := SMINI;
                CSI := 0;
                CNORM := SMINI;
                INFO := 1;
            end;
            
            //
            // Check scaling for  X = B / C
            //
            BNORM := AbsReal(B[1,1])+AbsReal(B[1,2]);
            if AP_FP_Less(CNORM,1) and AP_FP_Greater(BNORM,1) then
            begin
                if AP_FP_Greater(BNORM,BIGNUM*CNORM) then
                begin
                    SCL := 1/BNORM;
                end;
            end;
            
            //
            // Compute X
            //
            InternalHSEVDLADIV(SCL*B[1,1], SCL*B[1,2], CSR, CSI, TMP1, TMP2);
            X[1,1] := TMP1;
            X[1,2] := TMP2;
            XNORM := AbsReal(X[1,1])+AbsReal(X[1,2]);
        end;
    end
    else
    begin
        
        //
        // 2x2 System
        //
        // Compute the real part of  C = ca A - w D  (or  ca A' - w D )
        //
        CRV4[1+0] := CA*A[1,1]-WR*D1;
        CRV4[2+2] := CA*A[2,2]-WR*D2;
        if LTRANS then
        begin
            CRV4[1+2] := CA*A[2,1];
            CRV4[2+0] := CA*A[1,2];
        end
        else
        begin
            CRV4[2+0] := CA*A[2,1];
            CRV4[1+2] := CA*A[1,2];
        end;
        if NW=1 then
        begin
            
            //
            // Real 2x2 system  (w is real)
            //
            // Find the largest element in C
            //
            CMAX := 0;
            ICMAX := 0;
            J:=1;
            while J<=4 do
            begin
                if AP_FP_Greater(AbsReal(CRV4[J]),CMAX) then
                begin
                    CMAX := AbsReal(CRV4[J]);
                    ICMAX := J;
                end;
                Inc(J);
            end;
            
            //
            // If norm(C) < SMINI, use SMINI*identity.
            //
            if AP_FP_Less(CMAX,SMINI) then
            begin
                BNORM := Max(AbsReal(B[1,1]), AbsReal(B[2,1]));
                if AP_FP_Less(SMINI,1) and AP_FP_Greater(BNORM,1) then
                begin
                    if AP_FP_Greater(BNORM,BIGNUM*SMINI) then
                    begin
                        SCL := 1/BNORM;
                    end;
                end;
                TEMP := SCL/SMINI;
                X[1,1] := TEMP*B[1,1];
                X[2,1] := TEMP*B[2,1];
                XNORM := TEMP*BNORM;
                INFO := 1;
                Exit;
            end;
            
            //
            // Gaussian elimination with complete pivoting.
            //
            UR11 := CRV4[ICMAX];
            CR21 := CRV4[IPIVOT44[2,ICMAX]];
            UR12 := CRV4[IPIVOT44[3,ICMAX]];
            CR22 := CRV4[IPIVOT44[4,ICMAX]];
            UR11R := 1/UR11;
            LR21 := UR11R*CR21;
            UR22 := CR22-UR12*LR21;
            
            //
            // If smaller pivot < SMINI, use SMINI
            //
            if AP_FP_Less(AbsReal(UR22),SMINI) then
            begin
                UR22 := SMINI;
                INFO := 1;
            end;
            if RSWAP4[ICMAX] then
            begin
                BR1 := B[2,1];
                BR2 := B[1,1];
            end
            else
            begin
                BR1 := B[1,1];
                BR2 := B[2,1];
            end;
            BR2 := BR2-LR21*BR1;
            BBND := Max(AbsReal(BR1*(UR22*UR11R)), AbsReal(BR2));
            if AP_FP_Greater(BBND,1) and AP_FP_Less(AbsReal(UR22),1) then
            begin
                if AP_FP_Greater_Eq(BBND,BIGNUM*AbsReal(UR22)) then
                begin
                    SCL := 1/BBND;
                end;
            end;
            XR2 := BR2*SCL/UR22;
            XR1 := SCL*BR1*UR11R-XR2*(UR11R*UR12);
            if ZSWAP4[ICMAX] then
            begin
                X[1,1] := XR2;
                X[2,1] := XR1;
            end
            else
            begin
                X[1,1] := XR1;
                X[2,1] := XR2;
            end;
            XNORM := Max(AbsReal(XR1), AbsReal(XR2));
            
            //
            // Further scaling if  norm(A) norm(X) > overflow
            //
            if AP_FP_Greater(XNORM,1) and AP_FP_Greater(CMAX,1) then
            begin
                if AP_FP_Greater(XNORM,BIGNUM/CMAX) then
                begin
                    TEMP := CMAX/BIGNUM;
                    X[1,1] := TEMP*X[1,1];
                    X[2,1] := TEMP*X[2,1];
                    XNORM := TEMP*XNORM;
                    SCL := TEMP*SCL;
                end;
            end;
        end
        else
        begin
            
            //
            // Complex 2x2 system  (w is complex)
            //
            // Find the largest element in C
            //
            CIV4[1+0] := -WI*D1;
            CIV4[2+0] := 0;
            CIV4[1+2] := 0;
            CIV4[2+2] := -WI*D2;
            CMAX := 0;
            ICMAX := 0;
            J:=1;
            while J<=4 do
            begin
                if AP_FP_Greater(AbsReal(CRV4[J])+AbsReal(CIV4[J]),CMAX) then
                begin
                    CMAX := AbsReal(CRV4[J])+AbsReal(CIV4[J]);
                    ICMAX := J;
                end;
                Inc(J);
            end;
            
            //
            // If norm(C) < SMINI, use SMINI*identity.
            //
            if AP_FP_Less(CMAX,SMINI) then
            begin
                BNORM := Max(AbsReal(B[1,1])+AbsReal(B[1,2]), AbsReal(B[2,1])+AbsReal(B[2,2]));
                if AP_FP_Less(SMINI,1) and AP_FP_Greater(BNORM,1) then
                begin
                    if AP_FP_Greater(BNORM,BIGNUM*SMINI) then
                    begin
                        SCL := 1/BNORM;
                    end;
                end;
                TEMP := SCL/SMINI;
                X[1,1] := TEMP*B[1,1];
                X[2,1] := TEMP*B[2,1];
                X[1,2] := TEMP*B[1,2];
                X[2,2] := TEMP*B[2,2];
                XNORM := TEMP*BNORM;
                INFO := 1;
                Exit;
            end;
            
            //
            // Gaussian elimination with complete pivoting.
            //
            UR11 := CRV4[ICMAX];
            UI11 := CIV4[ICMAX];
            CR21 := CRV4[IPIVOT44[2,ICMAX]];
            CI21 := CIV4[IPIVOT44[2,ICMAX]];
            UR12 := CRV4[IPIVOT44[3,ICMAX]];
            UI12 := CIV4[IPIVOT44[3,ICMAX]];
            CR22 := CRV4[IPIVOT44[4,ICMAX]];
            CI22 := CIV4[IPIVOT44[4,ICMAX]];
            if (ICMAX=1) or (ICMAX=4) then
            begin
                
                //
                // Code when off-diagonals of pivoted C are real
                //
                if AP_FP_Greater(AbsReal(UR11),AbsReal(UI11)) then
                begin
                    TEMP := UI11/UR11;
                    UR11R := 1/(UR11*(1+AP_Sqr(TEMP)));
                    UI11R := -TEMP*UR11R;
                end
                else
                begin
                    TEMP := UR11/UI11;
                    UI11R := -1/(UI11*(1+AP_Sqr(TEMP)));
                    UR11R := -TEMP*UI11R;
                end;
                LR21 := CR21*UR11R;
                LI21 := CR21*UI11R;
                UR12S := UR12*UR11R;
                UI12S := UR12*UI11R;
                UR22 := CR22-UR12*LR21;
                UI22 := CI22-UR12*LI21;
            end
            else
            begin
                
                //
                // Code when diagonals of pivoted C are real
                //
                UR11R := 1/UR11;
                UI11R := 0;
                LR21 := CR21*UR11R;
                LI21 := CI21*UR11R;
                UR12S := UR12*UR11R;
                UI12S := UI12*UR11R;
                UR22 := CR22-UR12*LR21+UI12*LI21;
                UI22 := -UR12*LI21-UI12*LR21;
            end;
            U22ABS := AbsReal(UR22)+AbsReal(UI22);
            
            //
            // If smaller pivot < SMINI, use SMINI
            //
            if AP_FP_Less(U22ABS,SMINI) then
            begin
                UR22 := SMINI;
                UI22 := 0;
                INFO := 1;
            end;
            if RSWAP4[ICMAX] then
            begin
                BR2 := B[1,1];
                BR1 := B[2,1];
                BI2 := B[1,2];
                BI1 := B[2,2];
            end
            else
            begin
                BR1 := B[1,1];
                BR2 := B[2,1];
                BI1 := B[1,2];
                BI2 := B[2,2];
            end;
            BR2 := BR2-LR21*BR1+LI21*BI1;
            BI2 := BI2-LI21*BR1-LR21*BI1;
            BBND := Max((AbsReal(BR1)+AbsReal(BI1))*(U22ABS*(AbsReal(UR11R)+AbsReal(UI11R))), AbsReal(BR2)+AbsReal(BI2));
            if AP_FP_Greater(BBND,1) and AP_FP_Less(U22ABS,1) then
            begin
                if AP_FP_Greater_Eq(BBND,BIGNUM*U22ABS) then
                begin
                    SCL := 1/BBND;
                    BR1 := SCL*BR1;
                    BI1 := SCL*BI1;
                    BR2 := SCL*BR2;
                    BI2 := SCL*BI2;
                end;
            end;
            InternalHSEVDLADIV(BR2, BI2, UR22, UI22, XR2, XI2);
            XR1 := UR11R*BR1-UI11R*BI1-UR12S*XR2+UI12S*XI2;
            XI1 := UI11R*BR1+UR11R*BI1-UI12S*XR2-UR12S*XI2;
            if ZSWAP4[ICMAX] then
            begin
                X[1,1] := XR2;
                X[2,1] := XR1;
                X[1,2] := XI2;
                X[2,2] := XI1;
            end
            else
            begin
                X[1,1] := XR1;
                X[2,1] := XR2;
                X[1,2] := XI1;
                X[2,2] := XI2;
            end;
            XNORM := Max(AbsReal(XR1)+AbsReal(XI1), AbsReal(XR2)+AbsReal(XI2));
            
            //
            // Further scaling if  norm(A) norm(X) > overflow
            //
            if AP_FP_Greater(XNORM,1) and AP_FP_Greater(CMAX,1) then
            begin
                if AP_FP_Greater(XNORM,BIGNUM/CMAX) then
                begin
                    TEMP := CMAX/BIGNUM;
                    X[1,1] := TEMP*X[1,1];
                    X[2,1] := TEMP*X[2,1];
                    X[1,2] := TEMP*X[1,2];
                    X[2,2] := TEMP*X[2,2];
                    XNORM := TEMP*XNORM;
                    SCL := TEMP*SCL;
                end;
            end;
        end;
    end;
end;


procedure InternalHSEVDLADIV(const A : Double;
     const B : Double;
     const C : Double;
     const D : Double;
     var P : Double;
     var Q : Double);
var
    E : Double;
    F : Double;
begin
    if AP_FP_Less(AbsReal(D),AbsReal(C)) then
    begin
        E := D/C;
        F := C+D*E;
        P := (A+B*E)/F;
        Q := (B-A*E)/F;
    end
    else
    begin
        E := C/D;
        F := D+C*E;
        P := (B+A*E)/F;
        Q := (-A+B*E)/F;
    end;
end;


end.