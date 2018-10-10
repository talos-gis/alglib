(*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee. All rights reserved.

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
unit clu;
interface
uses Math, Sysutils, Ap;

procedure CMatrixLU(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
procedure ComplexLUDecomposition(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
procedure ComplexLUDecompositionUnpacked(A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var L : TComplex2DArray;
     var U : TComplex2DArray;
     var Pivots : TInteger1DArray);

implementation

(*************************************************************************
LU decomposition of a complex general matrix of size MxN

The subroutine calculates the LU decomposition of a rectangular general
matrix with partial pivoting (with row permutations).

Input parameters:
    A   -   matrix A whose indexes range within [0..M-1, 0..N-1].
    M   -   number of rows in matrix A.
    N   -   number of columns in matrix A.

Output parameters:
    A   -   matrices L and U in compact form (see below).
            Array whose indexes range within [0..M-1, 0..N-1].
    Pivots - permutation matrix in compact form (see below).
            Array whose index ranges within [0..Min(M-1,N-1)].

Matrix A is represented as A = P * L * U, where P is a permutation matrix,
matrix L - lower triangular (or lower trapezoid, if M>N) matrix,
U - upper triangular (or upper trapezoid, if M<N) matrix.

Let M be equal to 4 and N be equal to 3:

                   (  1          )    ( U11 U12 U13  )
A = P1 * P2 * P3 * ( L21  1      )  * (     U22 U23  )
                   ( L31 L32  1  )    (         U33  )
                   ( L41 L42 L43 )

Matrix L has size MxMin(M,N), matrix U has size Min(M,N)xN, matrix P(i) is
a permutation of the identity matrix of size MxM with numbers I and Pivots[I].

The algorithm returns array Pivots and the following matrix which replaces
matrix A and contains matrices L and U in compact form (the example applies
to M=4, N=3).

 ( U11 U12 U13 )
 ( L21 U22 U23 )
 ( L31 L32 U33 )
 ( L41 L42 L43 )

As we can see, the unit diagonal isn't stored.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1992
*************************************************************************)
procedure CMatrixLU(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    JP : AlglibInteger;
    T1 : TComplex1DArray;
    S : Complex;
    i_ : AlglibInteger;
begin
    SetLength(Pivots, Min(M-1, N-1)+1);
    SetLength(T1, Max(M-1, N-1)+1);
    Assert((M>=0) and (N>=0), 'Error in LUDecomposition: incorrect function arguments');
    
    //
    // Quick return if possible
    //
    if (M=0) or (N=0) then
    begin
        Exit;
    end;
    J:=0;
    while J<=Min(M-1, N-1) do
    begin
        
        //
        // Find pivot and test for singularity.
        //
        JP := J;
        I:=J+1;
        while I<=M-1 do
        begin
            if AP_FP_Greater(AbsComplex(A[I,J]),AbsComplex(A[JP,J])) then
            begin
                JP := I;
            end;
            Inc(I);
        end;
        Pivots[J] := JP;
        if C_NotEqualR(A[JP,J],0) then
        begin
            
            //
            //Apply the interchange to rows
            //
            if JP<>J then
            begin
                for i_ := 0 to N-1 do
                begin
                    T1[i_] := A[J,i_];
                end;
                for i_ := 0 to N-1 do
                begin
                    A[J,i_] := A[JP,i_];
                end;
                for i_ := 0 to N-1 do
                begin
                    A[JP,i_] := T1[i_];
                end;
            end;
            
            //
            //Compute elements J+1:M of J-th column.
            //
            if J<M then
            begin
                JP := J+1;
                S := C_RDiv(1,A[J,J]);
                for i_ := JP to M-1 do
                begin
                    A[i_,J] := C_Mul(S, A[i_,J]);
                end;
            end;
        end;
        if J<Min(M, N)-1 then
        begin
            
            //
            //Update trailing submatrix.
            //
            JP := J+1;
            I:=J+1;
            while I<=M-1 do
            begin
                S := A[I,J];
                for i_ := JP to N-1 do
                begin
                    A[I,i_] := C_Sub(A[I,i_], C_Mul(S, A[J,i_]));
                end;
                Inc(I);
            end;
        end;
        Inc(J);
    end;
end;


procedure ComplexLUDecomposition(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    JP : AlglibInteger;
    T1 : TComplex1DArray;
    s : Complex;
    i_ : AlglibInteger;
begin
    SetLength(Pivots, Min(M, N)+1);
    SetLength(T1, Max(M, N)+1);
    Assert((M>=0) and (N>=0), 'Error in ComplexLUDecomposition: incorrect function arguments');
    
    //
    // Quick return if possible
    //
    if (M=0) or (N=0) then
    begin
        Exit;
    end;
    J:=1;
    while J<=Min(M, N) do
    begin
        
        //
        // Find pivot and test for singularity.
        //
        JP := J;
        I:=J+1;
        while I<=M do
        begin
            if AP_FP_Greater(AbsComplex(A[I,J]),AbsComplex(A[JP,J])) then
            begin
                JP := I;
            end;
            Inc(I);
        end;
        Pivots[J] := JP;
        if C_NotEqualR(A[JP,J],0) then
        begin
            
            //
            //Apply the interchange to rows
            //
            if JP<>J then
            begin
                for i_ := 1 to N do
                begin
                    T1[i_] := A[J,i_];
                end;
                for i_ := 1 to N do
                begin
                    A[J,i_] := A[JP,i_];
                end;
                for i_ := 1 to N do
                begin
                    A[JP,i_] := T1[i_];
                end;
            end;
            
            //
            //Compute elements J+1:M of J-th column.
            //
            if J<M then
            begin
                JP := J+1;
                S := C_RDiv(1,A[J,J]);
                for i_ := JP to M do
                begin
                    A[i_,J] := C_Mul(S, A[i_,J]);
                end;
            end;
        end;
        if J<Min(M, N) then
        begin
            
            //
            //Update trailing submatrix.
            //
            JP := J+1;
            I:=J+1;
            while I<=M do
            begin
                S := A[I,J];
                for i_ := JP to N do
                begin
                    A[I,i_] := C_Sub(A[I,i_], C_Mul(S, A[J,i_]));
                end;
                Inc(I);
            end;
        end;
        Inc(J);
    end;
end;


procedure ComplexLUDecompositionUnpacked(A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var L : TComplex2DArray;
     var U : TComplex2DArray;
     var Pivots : TInteger1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    MinMN : AlglibInteger;
begin
    A := DynamicArrayCopy(A);
    if (M=0) or (N=0) then
    begin
        Exit;
    end;
    MinMN := Min(M, N);
    SetLength(L, M+1, MinMN+1);
    SetLength(U, MinMN+1, N+1);
    ComplexLUDecomposition(A, M, N, Pivots);
    I:=1;
    while I<=M do
    begin
        J:=1;
        while J<=MinMN do
        begin
            if J>I then
            begin
                L[I,J] := C_Complex(0);
            end;
            if J=I then
            begin
                L[I,J] := C_Complex(1);
            end;
            if J<I then
            begin
                L[I,J] := A[I,J];
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=1;
    while I<=MinMN do
    begin
        J:=1;
        while J<=N do
        begin
            if J<I then
            begin
                U[I,J] := C_Complex(0);
            end;
            if J>=I then
            begin
                U[I,J] := A[I,J];
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


end.