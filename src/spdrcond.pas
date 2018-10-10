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
unit spdrcond;
interface
uses Math, Sysutils, Ap, trlinsolve, cholesky, estnorm;

function SPDMatrixRCond(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
function SPDMatrixCholeskyRCond(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
function RCondSPD(A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
function RCondCholesky(const CD : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
procedure InternalCholeskyRCond(const ChFrm : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsNormProvided : Boolean;
     ANORM : Double;
     var RCOND : Double);

implementation

(*************************************************************************
Condition number estimate of a symmetric positive definite matrix.

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm of condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    A       -   symmetric positive definite matrix which is given by its
                upper or lower triangle depending on the value of
                IsUpper. Array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format.

Result:
    1/LowerBound(cond(A)), if matrix A is positive definite,
   -1, if matrix A is not positive definite, and its condition number
    could not be found by this algorithm.
*************************************************************************)
function SPDMatrixRCond(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    A1 : TReal2DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    IM : AlglibInteger;
    JM : AlglibInteger;
    V : Double;
    Nrm : Double;
    Pivots : TInteger1DArray;
begin
    SetLength(A1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        if IsUpper then
        begin
            J:=I;
            while J<=N do
            begin
                A1[I,J] := A[I-1,J-1];
                Inc(J);
            end;
        end
        else
        begin
            J:=1;
            while J<=I do
            begin
                A1[I,J] := A[I-1,J-1];
                Inc(J);
            end;
        end;
        Inc(I);
    end;
    Nrm := 0;
    J:=1;
    while J<=N do
    begin
        V := 0;
        I:=1;
        while I<=N do
        begin
            IM := I;
            JM := J;
            if IsUpper and (J<I) then
            begin
                IM := J;
                JM := I;
            end;
            if  not IsUpper and (J>I) then
            begin
                IM := J;
                JM := I;
            end;
            V := V+AbsReal(A1[IM,JM]);
            Inc(I);
        end;
        Nrm := Max(Nrm, V);
        Inc(J);
    end;
    if CholeskyDecomposition(A1, N, IsUpper) then
    begin
        InternalCholeskyRCond(A1, N, IsUpper, True, Nrm, V);
        Result := V;
    end
    else
    begin
        Result := -1;
    end;
end;


(*************************************************************************
Condition number estimate of a symmetric positive definite matrix given by
Cholesky decomposition.

The algorithm calculates a lower bound of the condition number. In this
case, the algorithm does not return a lower bound of the condition number,
but an inverse number (to avoid an overflow in case of a singular matrix).

It should be noted that 1-norm and inf-norm condition numbers of symmetric
matrices are equal, so the algorithm doesn't take into account the
differences between these types of norms.

Input parameters:
    CD  - Cholesky decomposition of matrix A,
          output of SMatrixCholesky subroutine.
    N   - size of matrix A.

Result: 1/LowerBound(cond(A))
*************************************************************************)
function SPDMatrixCholeskyRCond(const A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    A1 : TReal2DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
begin
    SetLength(A1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        if IsUpper then
        begin
            J:=I;
            while J<=N do
            begin
                A1[I,J] := A[I-1,J-1];
                Inc(J);
            end;
        end
        else
        begin
            J:=1;
            while J<=I do
            begin
                A1[I,J] := A[I-1,J-1];
                Inc(J);
            end;
        end;
        Inc(I);
    end;
    InternalCholeskyRCond(A1, N, IsUpper, False, 0, V);
    Result := V;
end;


function RCondSPD(A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    IM : AlglibInteger;
    JM : AlglibInteger;
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
            IM := I;
            JM := J;
            if IsUpper and (J<I) then
            begin
                IM := J;
                JM := I;
            end;
            if  not IsUpper and (J>I) then
            begin
                IM := J;
                JM := I;
            end;
            V := V+AbsReal(A[IM,JM]);
            Inc(I);
        end;
        Nrm := Max(Nrm, V);
        Inc(J);
    end;
    if CholeskyDecomposition(A, N, IsUpper) then
    begin
        InternalCholeskyRCond(A, N, IsUpper, True, Nrm, V);
        Result := V;
    end
    else
    begin
        Result := -1;
    end;
end;


function RCondCholesky(const CD : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Double;
var
    V : Double;
begin
    InternalCholeskyRCond(CD, N, IsUpper, False, 0, V);
    Result := V;
end;


procedure InternalCholeskyRCond(const ChFrm : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsNormProvided : Boolean;
     ANORM : Double;
     var RCOND : Double);
var
    NORMIN : Boolean;
    I : AlglibInteger;
    IX : AlglibInteger;
    KASE : AlglibInteger;
    AINVNM : Double;
    Scl : Double;
    SCALEL : Double;
    SCALEU : Double;
    SMLNUM : Double;
    WORK0 : TReal1DArray;
    WORK1 : TReal1DArray;
    WORK2 : TReal1DArray;
    IWORK : TInteger1DArray;
    V : Double;
    i_ : AlglibInteger;
begin
    Assert(N>=0);
    
    //
    // Estimate the norm of A.
    //
    if  not IsNormProvided then
    begin
        KASE := 0;
        ANORM := 0;
        while True do
        begin
            IterativeEstimate1Norm(N, WORK1, WORK0, IWORK, ANORM, KASE);
            if KASE=0 then
            begin
                Break;
            end;
            if IsUpper then
            begin
                
                //
                // Multiply by U
                //
                I:=1;
                while I<=N do
                begin
                    V := APVDotProduct(@ChFrm[I][0], I, N, @WORK0[0], I, N);
                    WORK0[I] := V;
                    Inc(I);
                end;
                
                //
                // Multiply by U'
                //
                I:=N;
                while I>=1 do
                begin
                    V := 0.0;
                    for i_ := 1 to I do
                    begin
                        V := V + ChFrm[i_,I]*WORK0[i_];
                    end;
                    WORK0[I] := V;
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
                    V := 0.0;
                    for i_ := I to N do
                    begin
                        V := V + ChFrm[i_,I]*WORK0[i_];
                    end;
                    WORK0[I] := V;
                    Inc(I);
                end;
                
                //
                // Multiply by L
                //
                I:=N;
                while I>=1 do
                begin
                    V := APVDotProduct(@ChFrm[I][0], 1, I, @WORK0[0], 1, I);
                    WORK0[I] := V;
                    Dec(I);
                end;
            end;
        end;
    end;
    
    //
    // Quick return if possible
    //
    RCOND := 0;
    if N=0 then
    begin
        RCOND := 1;
        Exit;
    end;
    if AP_FP_Eq(ANORM,0) then
    begin
        Exit;
    end;
    SMLNUM := MinRealNumber;
    
    //
    // Estimate the 1-norm of inv(A).
    //
    KASE := 0;
    NORMIN := False;
    while True do
    begin
        IterativeEstimate1Norm(N, WORK1, WORK0, IWORK, AINVNM, KASE);
        if KASE=0 then
        begin
            Break;
        end;
        if IsUpper then
        begin
            
            //
            // Multiply by inv(U').
            //
            SafeSolveTriangular(ChFrm, N, WORK0, SCALEL, IsUpper, True, False, NORMIN, WORK2);
            NORMIN := True;
            
            //
            // Multiply by inv(U).
            //
            SafeSolveTriangular(ChFrm, N, WORK0, SCALEU, IsUpper, False, False, NORMIN, WORK2);
        end
        else
        begin
            
            //
            // Multiply by inv(L).
            //
            SafeSolveTriangular(ChFrm, N, WORK0, SCALEL, IsUpper, False, False, NORMIN, WORK2);
            NORMIN := True;
            
            //
            // Multiply by inv(L').
            //
            SafeSolveTriangular(ChFrm, N, WORK0, SCALEU, IsUpper, True, False, NORMIN, WORK2);
        end;
        
        //
        // Multiply by 1/SCALE if doing so will not cause overflow.
        //
        Scl := SCALEL*SCALEU;
        if AP_FP_Neq(Scl,1) then
        begin
            IX := 1;
            I:=2;
            while I<=N do
            begin
                if AP_FP_Greater(AbsReal(WORK0[I]),AbsReal(WORK0[IX])) then
                begin
                    IX := I;
                end;
                Inc(I);
            end;
            if AP_FP_Less(Scl,AbsReal(WORK0[IX])*SMLNUM) or AP_FP_Eq(Scl,0) then
            begin
                Exit;
            end;
            I:=1;
            while I<=N do
            begin
                WORK0[I] := WORK0[I]/Scl;
                Inc(I);
            end;
        end;
    end;
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if AP_FP_Neq(AINVNM,0) then
    begin
        V := 1/AINVNM;
        RCOND := V/ANORM;
    end;
end;


end.