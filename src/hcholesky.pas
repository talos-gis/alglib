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
unit hcholesky;
interface
uses Math, Sysutils, Ap, cblas;

function HMatrixCholesky(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
function HermitianCholeskyDecomposition(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;

implementation

(*************************************************************************
Cholesky decomposition

The algorithm computes Cholesky decomposition  of  a  Hermitian  positive-
definite matrix.

The result of an algorithm is a representation of matrix A as A = U'*U  or
A = L*L' (here X' detones conj(X^T)).

Input parameters:
    A       -   upper or lower triangle of a factorized matrix.
                array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   if IsUpper=True, then A contains an upper triangle of
                a symmetric matrix, otherwise A contains a lower one.

Output parameters:
    A       -   the result of factorization. If IsUpper=True, then
                the upper triangle contains matrix U, so that A = U'*U,
                and the elements below the main diagonal are not modified.
                Similarly, if IsUpper = False.

Result:
    If the matrix is positive-definite, the function returns True.
    Otherwise, the function returns False. This means that the
    factorization could not be carried out.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************)
function HMatrixCholesky(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
var
    J : AlglibInteger;
    AJJ : Double;
    V : Complex;
    R : Double;
    T : TComplex1DArray;
    T2 : TComplex1DArray;
    T3 : TComplex1DArray;
    I : AlglibInteger;
    A1 : TComplex2DArray;
    i_ : AlglibInteger;
begin
    if  not IsUpper then
    begin
        SetLength(A1, N+1, N+1);
        I:=1;
        while I<=N do
        begin
            J:=1;
            while J<=N do
            begin
                A1[I,J] := A[I-1,J-1];
                Inc(J);
            end;
            Inc(I);
        end;
        Result := HermitianCholeskyDecomposition(A1, N, IsUpper);
        I:=1;
        while I<=N do
        begin
            J:=1;
            while J<=N do
            begin
                A[I-1,J-1] := A1[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
        Exit;
    end;
    SetLength(T, N-1+1);
    SetLength(T2, N-1+1);
    SetLength(T3, N-1+1);
    Result := True;
    if N<0 then
    begin
        Result := False;
        Exit;
    end;
    
    //
    // Quick return if possible
    //
    if N=0 then
    begin
        Exit;
    end;
    if IsUpper then
    begin
        
        //
        // Compute the Cholesky factorization A = U'*U.
        //
        J:=0;
        while J<=N-1 do
        begin
            
            //
            // Compute U(J,J) and test for non-positive-definiteness.
            //
            V := C_Complex(0.0);
            for i_ := 0 to J-1 do
            begin
                V := C_Add(V,C_Mul(Conj(A[i_,J]),A[i_,J]));
            end;
            AJJ := C_Sub(A[J,J],V).X;
            if AP_FP_Less_Eq(AJJ,0) then
            begin
                A[J,J] := C_Complex(AJJ);
                Result := False;
                Exit;
            end;
            AJJ := SQRT(AJJ);
            A[J,J] := C_Complex(AJJ);
            
            //
            // Compute elements J+1:N-1 of row J.
            //
            if J<N-1 then
            begin
                for i_ := 0 to J-1 do
                begin
                    T2[i_] := Conj(A[i_,J]);
                end;
                for i_ := J+1 to N-1 do
                begin
                    T3[i_] := A[J,i_];
                end;
                ComplexMatrixVectorMultiply(A, 0, J-1, J+1, N-1, True, False, T2, 0, J-1, C_Complex(-1.0), T3, J+1, N-1, C_Complex(1.0), T);
                for i_ := J+1 to N-1 do
                begin
                    A[J,i_] := T3[i_];
                end;
                R := 1/AJJ;
                for i_ := J+1 to N-1 do
                begin
                    A[J,i_] := C_MulR(A[J,i_],R);
                end;
            end;
            Inc(J);
        end;
    end
    else
    begin
        
        //
        // Compute the Cholesky factorization A = L*L'.
        //
        J:=0;
        while J<=N-1 do
        begin
            
            //
            // Compute L(J+1,J+1) and test for non-positive-definiteness.
            //
            V := C_Complex(0.0);
            for i_ := 0 to J-1 do
            begin
                V := C_Add(V,C_Mul(Conj(A[J,i_]),A[J,i_]));
            end;
            AJJ := C_Sub(A[J,J],V).X;
            if AP_FP_Less_Eq(AJJ,0) then
            begin
                A[J,J] := C_Complex(AJJ);
                Result := False;
                Exit;
            end;
            AJJ := SQRT(AJJ);
            A[J,J] := C_Complex(AJJ);
            
            //
            // Compute elements J+1:N of column J.
            //
            if J<N-1 then
            begin
                for i_ := 0 to J-1 do
                begin
                    T2[i_] := Conj(A[J,i_]);
                end;
                for i_ := J+1 to N-1 do
                begin
                    T3[i_] := A[i_,J];
                end;
                ComplexMatrixVectorMultiply(A, J+1, N-1, 0, J-1, False, False, T2, 0, J-1, C_Complex(-1.0), T3, J+1, N-1, C_Complex(1.0), T);
                for i_ := J+1 to N-1 do
                begin
                    A[i_,J] := T3[i_];
                end;
                R := 1/AJJ;
                for i_ := J+1 to N-1 do
                begin
                    A[i_,J] := C_MulR(A[i_,J],R);
                end;
            end;
            Inc(J);
        end;
    end;
end;


function HermitianCholeskyDecomposition(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
var
    J : AlglibInteger;
    AJJ : Double;
    V : Complex;
    R : Double;
    T : TComplex1DArray;
    T2 : TComplex1DArray;
    T3 : TComplex1DArray;
    i_ : AlglibInteger;
begin
    SetLength(T, N+1);
    SetLength(T2, N+1);
    SetLength(T3, N+1);
    Result := True;
    if N<0 then
    begin
        Result := False;
        Exit;
    end;
    
    //
    // Quick return if possible
    //
    if N=0 then
    begin
        Exit;
    end;
    if IsUpper then
    begin
        
        //
        // Compute the Cholesky factorization A = U'*U.
        //
        J:=1;
        while J<=N do
        begin
            
            //
            // Compute U(J,J) and test for non-positive-definiteness.
            //
            V := C_Complex(0.0);
            for i_ := 1 to J-1 do
            begin
                V := C_Add(V,C_Mul(Conj(A[i_,J]),A[i_,J]));
            end;
            AJJ := C_Sub(A[J,J],V).X;
            if AP_FP_Less_Eq(AJJ,0) then
            begin
                A[J,J] := C_Complex(AJJ);
                Result := False;
                Exit;
            end;
            AJJ := SQRT(AJJ);
            A[J,J] := C_Complex(AJJ);
            
            //
            // Compute elements J+1:N of row J.
            //
            if J<N then
            begin
                for i_ := 1 to J-1 do
                begin
                    A[i_,J] := Conj(A[i_,J]);
                end;
                for i_ := 1 to J-1 do
                begin
                    T2[i_] := A[i_,J];
                end;
                for i_ := J+1 to N do
                begin
                    T3[i_] := A[J,i_];
                end;
                ComplexMatrixVectorMultiply(A, 1, J-1, J+1, N, True, False, T2, 1, J-1, C_Complex(-1.0), T3, J+1, N, C_Complex(1.0), T);
                for i_ := J+1 to N do
                begin
                    A[J,i_] := T3[i_];
                end;
                for i_ := 1 to J-1 do
                begin
                    A[i_,J] := Conj(A[i_,J]);
                end;
                R := 1/AJJ;
                for i_ := J+1 to N do
                begin
                    A[J,i_] := C_MulR(A[J,i_],R);
                end;
            end;
            Inc(J);
        end;
    end
    else
    begin
        
        //
        // Compute the Cholesky factorization A = L*L'.
        //
        J:=1;
        while J<=N do
        begin
            
            //
            // Compute L(J,J) and test for non-positive-definiteness.
            //
            V := C_Complex(0.0);
            for i_ := 1 to J-1 do
            begin
                V := C_Add(V,C_Mul(Conj(A[J,i_]),A[J,i_]));
            end;
            AJJ := C_Sub(A[J,J],V).X;
            if AP_FP_Less_Eq(AJJ,0) then
            begin
                A[J,J] := C_Complex(AJJ);
                Result := False;
                Exit;
            end;
            AJJ := SQRT(AJJ);
            A[J,J] := C_Complex(AJJ);
            
            //
            // Compute elements J+1:N of column J.
            //
            if J<N then
            begin
                for i_ := 1 to J-1 do
                begin
                    A[J,i_] := Conj(A[J,i_]);
                end;
                for i_ := 1 to J-1 do
                begin
                    T2[i_] := A[J,i_];
                end;
                for i_ := J+1 to N do
                begin
                    T3[i_] := A[i_,J];
                end;
                ComplexMatrixVectorMultiply(A, J+1, N, 1, J-1, False, False, T2, 1, J-1, C_Complex(-1.0), T3, J+1, N, C_Complex(1.0), T);
                for i_ := J+1 to N do
                begin
                    A[i_,J] := T3[i_];
                end;
                for i_ := 1 to J-1 do
                begin
                    A[J,i_] := Conj(A[J,i_]);
                end;
                R := 1/AJJ;
                for i_ := J+1 to N do
                begin
                    A[i_,J] := C_MulR(A[i_,J],R);
                end;
            end;
            Inc(J);
        end;
    end;
end;


end.