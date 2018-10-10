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
unit crcond;
interface
uses Math, Sysutils, Ap, clu, ctrlinsolve;

function CMatrixRCond1(const A : TComplex2DArray; N : AlglibInteger):Double;
function CMatrixLURCond1(const LUDcmp : TComplex2DArray;
     N : AlglibInteger):Double;
function CMatrixRCondInf(const A : TComplex2DArray; N : AlglibInteger):Double;
function CMatrixLURCondInf(const LUDcmp : TComplex2DArray;
     N : AlglibInteger):Double;
function ComplexRCond1(A : TComplex2DArray; N : AlglibInteger):Double;
function ComplexRCond1LU(const LU : TComplex2DArray; N : AlglibInteger):Double;
function ComplexRCondInf(A : TComplex2DArray; N : AlglibInteger):Double;
function ComplexRCondInfLU(const LU : TComplex2DArray;
     N : AlglibInteger):Double;
procedure InternalEstimateComplexRCondLU(const LU : TComplex2DArray;
     const N : AlglibInteger;
     OneNorm : Boolean;
     IsANormProvided : Boolean;
     ANORM : Double;
     var RCOND : Double);

implementation

procedure InternalComplexRCondEstimateNorm(const N : AlglibInteger;
     var V : TComplex1DArray;
     var X : TComplex1DArray;
     var EST : Double;
     var KASE : AlglibInteger;
     var ISAVE : TInteger1DArray;
     var RSAVE : TReal1DArray);forward;
function InternalComplexRCondSCSUM1(const X : TComplex1DArray;
     N : AlglibInteger):Double;forward;
function InternalComplexRCondICMAX1(const X : TComplex1DArray;
     N : AlglibInteger):AlglibInteger;forward;
procedure InternalComplexRCondSaveAll(var ISAVE : TInteger1DArray;
     var RSAVE : TReal1DArray;
     var I : AlglibInteger;
     var ITER : AlglibInteger;
     var J : AlglibInteger;
     var JLAST : AlglibInteger;
     var JUMP : AlglibInteger;
     var ABSXI : Double;
     var ALTSGN : Double;
     var ESTOLD : Double;
     var TEMP : Double);forward;
procedure InternalComplexRCondLoadAll(var ISAVE : TInteger1DArray;
     var RSAVE : TReal1DArray;
     var I : AlglibInteger;
     var ITER : AlglibInteger;
     var J : AlglibInteger;
     var JLAST : AlglibInteger;
     var JUMP : AlglibInteger;
     var ABSXI : Double;
     var ALTSGN : Double;
     var ESTOLD : Double;
     var TEMP : Double);forward;


(*************************************************************************
Estimate of a matrix condition number (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.

Result: 1/LowerBound(cond(A))
*************************************************************************)
function CMatrixRCond1(const A : TComplex2DArray; N : AlglibInteger):Double;
var
    I : AlglibInteger;
    A1 : TComplex2DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert(N>=1, 'CMatrixRCond1: N<1!');
    SetLength(A1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        i1_ := (0) - (1);
        for i_ := 1 to N do
        begin
            A1[I,i_] := A[I-1,i_+i1_];
        end;
        Inc(I);
    end;
    Result := ComplexRCond1(A1, N);
end;


(*************************************************************************
Estimate of the condition number of a matrix given by its LU decomposition (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    LUDcmp      -   LU decomposition of a matrix in compact form. Output of
                    the CMatrixLU subroutine.
    N           -   size of matrix A.

Result: 1/LowerBound(cond(A))
*************************************************************************)
function CMatrixLURCond1(const LUDcmp : TComplex2DArray;
     N : AlglibInteger):Double;
var
    I : AlglibInteger;
    A1 : TComplex2DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert(N>=1, 'CMatrixLURCond1: N<1!');
    SetLength(A1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        i1_ := (0) - (1);
        for i_ := 1 to N do
        begin
            A1[I,i_] := LUDcmp[I-1,i_+i1_];
        end;
        Inc(I);
    end;
    Result := ComplexRCond1LU(A1, N);
end;


(*************************************************************************
Estimate of a matrix condition number (infinity-norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.

Result: 1/LowerBound(cond(A))
*************************************************************************)
function CMatrixRCondInf(const A : TComplex2DArray; N : AlglibInteger):Double;
var
    I : AlglibInteger;
    A1 : TComplex2DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert(N>=1, 'CMatrixRCondInf: N<1!');
    SetLength(A1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        i1_ := (0) - (1);
        for i_ := 1 to N do
        begin
            A1[I,i_] := A[I-1,i_+i1_];
        end;
        Inc(I);
    end;
    Result := ComplexRCondInf(A1, N);
end;


(*************************************************************************
Estimate of the condition number of a matrix given by its LU decomposition
(infinity norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    LUDcmp  -   LU decomposition of a matrix in compact form. Output of
                the CMatrixLU subroutine.
    N       -   size of matrix A.

Result: 1/LowerBound(cond(A))
*************************************************************************)
function CMatrixLURCondInf(const LUDcmp : TComplex2DArray;
     N : AlglibInteger):Double;
var
    I : AlglibInteger;
    A1 : TComplex2DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert(N>=1, 'CMatrixLURCondInf: N<1!');
    SetLength(A1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        i1_ := (0) - (1);
        for i_ := 1 to N do
        begin
            A1[I,i_] := LUDcmp[I-1,i_+i1_];
        end;
        Inc(I);
    end;
    Result := ComplexRCondInfLU(A1, N);
end;


function ComplexRCond1(A : TComplex2DArray; N : AlglibInteger):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Nrm : Double;
    Pivots : TInteger1DArray;
begin
    A := DynamicArrayCopy(A);
    Nrm := 0;
    J:=1;
    while J<=N do
    begin
        V := 0;
        I:=1;
        while I<=N do
        begin
            V := V+AbsComplex(A[I,J]);
            Inc(I);
        end;
        Nrm := Max(Nrm, V);
        Inc(J);
    end;
    ComplexLUDecomposition(A, N, N, Pivots);
    InternalEstimateComplexRCondLU(A, N, True, True, Nrm, V);
    Result := V;
end;


function ComplexRCond1LU(const LU : TComplex2DArray; N : AlglibInteger):Double;
var
    V : Double;
begin
    InternalEstimateComplexRCondLU(LU, N, True, False, 0, V);
    Result := V;
end;


function ComplexRCondInf(A : TComplex2DArray; N : AlglibInteger):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Nrm : Double;
    Pivots : TInteger1DArray;
begin
    A := DynamicArrayCopy(A);
    Nrm := 0;
    I:=1;
    while I<=N do
    begin
        V := 0;
        J:=1;
        while J<=N do
        begin
            V := V+AbsComplex(A[I,J]);
            Inc(J);
        end;
        Nrm := Max(Nrm, V);
        Inc(I);
    end;
    ComplexLUDecomposition(A, N, N, Pivots);
    InternalEstimateComplexRCondLU(A, N, False, True, Nrm, V);
    Result := V;
end;


function ComplexRCondInfLU(const LU : TComplex2DArray;
     N : AlglibInteger):Double;
var
    V : Double;
begin
    InternalEstimateComplexRCondLU(LU, N, False, False, 0, V);
    Result := V;
end;


procedure InternalEstimateComplexRCondLU(const LU : TComplex2DArray;
     const N : AlglibInteger;
     OneNorm : Boolean;
     IsANormProvided : Boolean;
     ANORM : Double;
     var RCOND : Double);
var
    CWORK1 : TComplex1DArray;
    CWORK2 : TComplex1DArray;
    CWORK3 : TComplex1DArray;
    CWORK4 : TComplex1DArray;
    ISAVE : TInteger1DArray;
    RSAVE : TReal1DArray;
    KASE : AlglibInteger;
    KASE1 : AlglibInteger;
    AINVNM : Double;
    SMLNUM : Double;
    CW : Boolean;
    V : Complex;
    I : AlglibInteger;
    i_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Exit;
    end;
    SetLength(CWORK1, N+1);
    SetLength(CWORK2, N+1);
    SetLength(CWORK3, N+1);
    SetLength(CWORK4, N+1);
    SetLength(ISAVE, 4+1);
    SetLength(RSAVE, 3+1);
    RCOND := 0;
    if N=0 then
    begin
        RCOND := 1;
        Exit;
    end;
    SMLNUM := MinRealNumber;
    
    //
    // Estimate the norm of inv(A).
    //
    if  not IsANormProvided then
    begin
        ANORM := 0;
        if OneNorm then
        begin
            KASE1 := 1;
        end
        else
        begin
            KASE1 := 2;
        end;
        KASE := 0;
        repeat
            InternalComplexRCondEstimateNorm(N, CWORK4, CWORK1, ANORM, KASE, ISAVE, RSAVE);
            if KASE<>0 then
            begin
                if KASE=KASE1 then
                begin
                    
                    //
                    // Multiply by U
                    //
                    I:=1;
                    while I<=N do
                    begin
                        V := C_Complex(0.0);
                        for i_ := I to N do
                        begin
                            V := C_Add(V,C_Mul(LU[I,i_],CWORK1[i_]));
                        end;
                        CWORK1[I] := V;
                        Inc(I);
                    end;
                    
                    //
                    // Multiply by L
                    //
                    I:=N;
                    while I>=1 do
                    begin
                        V := C_Complex(0);
                        if I>1 then
                        begin
                            V := C_Complex(0.0);
                            for i_ := 1 to I-1 do
                            begin
                                V := C_Add(V,C_Mul(LU[I,i_],CWORK1[i_]));
                            end;
                        end;
                        CWORK1[I] := C_Add(V,CWORK1[I]);
                        Dec(I);
                    end;
                end
                else
                begin
                    
                    //
                    // Multiply by L'
                    //
                    I:=1;
                    while I<=N do
                    begin
                        CWORK2[I] := C_Complex(0);
                        Inc(I);
                    end;
                    I:=1;
                    while I<=N do
                    begin
                        V := CWORK1[I];
                        if I>1 then
                        begin
                            for i_ := 1 to I-1 do
                            begin
                                CWORK2[i_] := C_Add(CWORK2[i_], C_Mul(V, Conj(LU[I,i_])));
                            end;
                        end;
                        CWORK2[I] := C_Add(CWORK2[I],V);
                        Inc(I);
                    end;
                    
                    //
                    // Multiply by U'
                    //
                    I:=1;
                    while I<=N do
                    begin
                        CWORK1[I] := C_Complex(0);
                        Inc(I);
                    end;
                    I:=1;
                    while I<=N do
                    begin
                        V := CWORK2[I];
                        for i_ := I to N do
                        begin
                            CWORK1[i_] := C_Add(CWORK1[i_], C_Mul(V, Conj(LU[I,i_])));
                        end;
                        Inc(I);
                    end;
                end;
            end;
        until KASE=0;
    end;
    
    //
    // Quick return if possible
    //
    if AP_FP_Eq(ANORM,0) then
    begin
        Exit;
    end;
    
    //
    // Estimate the norm of inv(A).
    //
    AINVNM := 0;
    if OneNorm then
    begin
        KASE1 := 1;
    end
    else
    begin
        KASE1 := 2;
    end;
    KASE := 0;
    repeat
        InternalComplexRCondEstimateNorm(N, CWORK4, CWORK1, AINVNM, KASE, ISAVE, RSAVE);
        if KASE<>0 then
        begin
            if KASE=KASE1 then
            begin
                
                //
                // Multiply by inv(L).
                //
                CW := ComplexSafeSolveTriangular(LU, N, CWORK1, False, 0, True, CWORK2, CWORK3);
                if  not CW then
                begin
                    RCOND := 0;
                    Exit;
                end;
                
                //
                // Multiply by inv(U).
                //
                CW := ComplexSafeSolveTriangular(LU, N, CWORK1, True, 0, False, CWORK2, CWORK3);
                if  not CW then
                begin
                    RCOND := 0;
                    Exit;
                end;
            end
            else
            begin
                
                //
                // Multiply by inv(U').
                //
                CW := ComplexSafeSolveTriangular(LU, N, CWORK1, True, 2, False, CWORK2, CWORK3);
                if  not CW then
                begin
                    RCOND := 0;
                    Exit;
                end;
                
                //
                // Multiply by inv(L').
                //
                CW := ComplexSafeSolveTriangular(LU, N, CWORK1, False, 2, True, CWORK2, CWORK3);
                if  not CW then
                begin
                    RCOND := 0;
                    Exit;
                end;
            end;
        end;
    until KASE=0;
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if AP_FP_Neq(AINVNM,0) then
    begin
        RCOND := 1/AINVNM;
        RCOND := RCOND/ANORM;
    end;
end;


procedure InternalComplexRCondEstimateNorm(const N : AlglibInteger;
     var V : TComplex1DArray;
     var X : TComplex1DArray;
     var EST : Double;
     var KASE : AlglibInteger;
     var ISAVE : TInteger1DArray;
     var RSAVE : TReal1DArray);
var
    ITMAX : AlglibInteger;
    I : AlglibInteger;
    ITER : AlglibInteger;
    J : AlglibInteger;
    JLAST : AlglibInteger;
    JUMP : AlglibInteger;
    ABSXI : Double;
    ALTSGN : Double;
    ESTOLD : Double;
    SAFMIN : Double;
    TEMP : Double;
    i_ : AlglibInteger;
begin
    
    //
    //Executable Statements ..
    //
    ITMAX := 5;
    SAFMIN := MinRealNumber;
    if KASE=0 then
    begin
        I:=1;
        while I<=N do
        begin
            X[I] := C_Complex(AP_Double(1)/N);
            Inc(I);
        end;
        KASE := 1;
        JUMP := 1;
        InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
        Exit;
    end;
    InternalComplexRCondLoadAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
    
    //
    // ENTRY   (JUMP = 1)
    // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
    //
    if JUMP=1 then
    begin
        if N=1 then
        begin
            V[1] := X[1];
            EST := AbsComplex(V[1]);
            KASE := 0;
            InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
            Exit;
        end;
        EST := InternalComplexRCondSCSUM1(X, N);
        I:=1;
        while I<=N do
        begin
            ABSXI := AbsComplex(X[I]);
            if AP_FP_Greater(ABSXI,SAFMIN) then
            begin
                X[I] := C_DivR(X[I],ABSXI);
            end
            else
            begin
                X[I] := C_Complex(1);
            end;
            Inc(I);
        end;
        KASE := 2;
        JUMP := 2;
        InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
        Exit;
    end;
    
    //
    // ENTRY   (JUMP = 2)
    // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
    //
    if JUMP=2 then
    begin
        J := InternalComplexRCondICMAX1(X, N);
        ITER := 2;
        
        //
        // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
        //
        I:=1;
        while I<=N do
        begin
            X[I] := C_Complex(0);
            Inc(I);
        end;
        X[J] := C_Complex(1);
        KASE := 1;
        JUMP := 3;
        InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
        Exit;
    end;
    
    //
    // ENTRY   (JUMP = 3)
    // X HAS BEEN OVERWRITTEN BY A*X.
    //
    if JUMP=3 then
    begin
        for i_ := 1 to N do
        begin
            V[i_] := X[i_];
        end;
        ESTOLD := EST;
        EST := InternalComplexRCondSCSUM1(V, N);
        
        //
        // TEST FOR CYCLING.
        //
        if AP_FP_Less_Eq(EST,ESTOLD) then
        begin
            
            //
            // ITERATION COMPLETE.  FINAL STAGE.
            //
            ALTSGN := 1;
            I:=1;
            while I<=N do
            begin
                X[I] := C_Complex(ALTSGN*(1+AP_Double((I-1))/(N-1)));
                ALTSGN := -ALTSGN;
                Inc(I);
            end;
            KASE := 1;
            JUMP := 5;
            InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
            Exit;
        end;
        I:=1;
        while I<=N do
        begin
            ABSXI := AbsComplex(X[I]);
            if AP_FP_Greater(ABSXI,SAFMIN) then
            begin
                X[I] := C_DivR(X[I],ABSXI);
            end
            else
            begin
                X[I] := C_Complex(1);
            end;
            Inc(I);
        end;
        KASE := 2;
        JUMP := 4;
        InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
        Exit;
    end;
    
    //
    // ENTRY   (JUMP = 4)
    // X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
    //
    if JUMP=4 then
    begin
        JLAST := J;
        J := InternalComplexRCondICMAX1(X, N);
        if AP_FP_Neq(AbsComplex(X[JLAST]),AbsComplex(X[J])) and (ITER<ITMAX) then
        begin
            ITER := ITER+1;
            
            //
            // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
            //
            I:=1;
            while I<=N do
            begin
                X[I] := C_Complex(0);
                Inc(I);
            end;
            X[J] := C_Complex(1);
            KASE := 1;
            JUMP := 3;
            InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
            Exit;
        end;
        
        //
        // ITERATION COMPLETE.  FINAL STAGE.
        //
        ALTSGN := 1;
        I:=1;
        while I<=N do
        begin
            X[I] := C_Complex(ALTSGN*(1+AP_Double((I-1))/(N-1)));
            ALTSGN := -ALTSGN;
            Inc(I);
        end;
        KASE := 1;
        JUMP := 5;
        InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
        Exit;
    end;
    
    //
    // ENTRY   (JUMP = 5)
    // X HAS BEEN OVERWRITTEN BY A*X.
    //
    if JUMP=5 then
    begin
        TEMP := 2*(InternalComplexRCondSCSUM1(X, N)/(3*N));
        if AP_FP_Greater(TEMP,EST) then
        begin
            for i_ := 1 to N do
            begin
                V[i_] := X[i_];
            end;
            EST := TEMP;
        end;
        KASE := 0;
        InternalComplexRCondSaveAll(ISAVE, RSAVE, I, ITER, J, JLAST, JUMP, ABSXI, ALTSGN, ESTOLD, TEMP);
        Exit;
    end;
end;


function InternalComplexRCondSCSUM1(const X : TComplex1DArray;
     N : AlglibInteger):Double;
var
    I : AlglibInteger;
begin
    Result := 0;
    I:=1;
    while I<=N do
    begin
        Result := Result+AbsComplex(X[I]);
        Inc(I);
    end;
end;


function InternalComplexRCondICMAX1(const X : TComplex1DArray;
     N : AlglibInteger):AlglibInteger;
var
    I : AlglibInteger;
    M : Double;
begin
    Result := 1;
    M := AbsComplex(X[1]);
    I:=2;
    while I<=N do
    begin
        if AP_FP_Greater(AbsComplex(X[I]),M) then
        begin
            Result := I;
            M := AbsComplex(X[I]);
        end;
        Inc(I);
    end;
end;


procedure InternalComplexRCondSaveAll(var ISAVE : TInteger1DArray;
     var RSAVE : TReal1DArray;
     var I : AlglibInteger;
     var ITER : AlglibInteger;
     var J : AlglibInteger;
     var JLAST : AlglibInteger;
     var JUMP : AlglibInteger;
     var ABSXI : Double;
     var ALTSGN : Double;
     var ESTOLD : Double;
     var TEMP : Double);
begin
    ISAVE[0] := I;
    ISAVE[1] := ITER;
    ISAVE[2] := J;
    ISAVE[3] := JLAST;
    ISAVE[4] := JUMP;
    RSAVE[0] := ABSXI;
    RSAVE[1] := ALTSGN;
    RSAVE[2] := ESTOLD;
    RSAVE[3] := TEMP;
end;


procedure InternalComplexRCondLoadAll(var ISAVE : TInteger1DArray;
     var RSAVE : TReal1DArray;
     var I : AlglibInteger;
     var ITER : AlglibInteger;
     var J : AlglibInteger;
     var JLAST : AlglibInteger;
     var JUMP : AlglibInteger;
     var ABSXI : Double;
     var ALTSGN : Double;
     var ESTOLD : Double;
     var TEMP : Double);
begin
    I := ISAVE[0];
    ITER := ISAVE[1];
    J := ISAVE[2];
    JLAST := ISAVE[3];
    JUMP := ISAVE[4];
    ABSXI := RSAVE[0];
    ALTSGN := RSAVE[1];
    ESTOLD := RSAVE[2];
    TEMP := RSAVE[3];
end;


end.