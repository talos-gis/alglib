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
unit inv;
interface
uses Math, Sysutils, Ap, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trinverse;

function RMatrixLUInverse(var A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger):Boolean;
function RMatrixInverse(var A : TReal2DArray; N : AlglibInteger):Boolean;

implementation

(*************************************************************************
Inversion of a matrix given by its LU decomposition.

Input parameters:
    A       -   LU decomposition of the matrix (output of RMatrixLU subroutine).
    Pivots  -   table of permutations which were made during the LU decomposition
                (the output of RMatrixLU subroutine).
    N       -   size of matrix A.

Output parameters:
    A       -   inverse of matrix A.
                Array whose indexes range within [0..N-1, 0..N-1].

Result:
    True, if the matrix is not singular.
    False, if the matrix is singular.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************)
function RMatrixLUInverse(var A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger):Boolean;
var
    WORK : TReal1DArray;
    I : AlglibInteger;
    IWS : AlglibInteger;
    J : AlglibInteger;
    JB : AlglibInteger;
    JJ : AlglibInteger;
    JP : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
begin
    Result := True;
    
    //
    // Quick return if possible
    //
    if N=0 then
    begin
        Exit;
    end;
    SetLength(WORK, N-1+1);
    
    //
    // Form inv(U)
    //
    if  not RMatrixTRInverse(A, N, True, False) then
    begin
        Result := False;
        Exit;
    end;
    
    //
    // Solve the equation inv(A)*L = inv(U) for inv(A).
    //
    J:=N-1;
    while J>=0 do
    begin
        
        //
        // Copy current column of L to WORK and replace with zeros.
        //
        I:=J+1;
        while I<=N-1 do
        begin
            WORK[I] := A[I,J];
            A[I,J] := 0;
            Inc(I);
        end;
        
        //
        // Compute current column of inv(A).
        //
        if J<N-1 then
        begin
            I:=0;
            while I<=N-1 do
            begin
                V := APVDotProduct(@A[I][0], J+1, N-1, @WORK[0], J+1, N-1);
                A[I,J] := A[I,J]-V;
                Inc(I);
            end;
        end;
        Dec(J);
    end;
    
    //
    // Apply column interchanges.
    //
    J:=N-2;
    while J>=0 do
    begin
        JP := Pivots[J];
        if JP<>J then
        begin
            for i_ := 0 to N-1 do
            begin
                WORK[i_] := A[i_,J];
            end;
            for i_ := 0 to N-1 do
            begin
                A[i_,J] := A[i_,JP];
            end;
            for i_ := 0 to N-1 do
            begin
                A[i_,JP] := WORK[i_];
            end;
        end;
        Dec(J);
    end;
end;


(*************************************************************************
Inversion of a general matrix.

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.

Output parameters:
    A   -   inverse of matrix A.
            Array whose indexes range within [0..N-1, 0..N-1].

Result:
    True, if the matrix is not singular.
    False, if the matrix is singular.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
function RMatrixInverse(var A : TReal2DArray; N : AlglibInteger):Boolean;
var
    Pivots : TInteger1DArray;
begin
    RMatrixLU(A, N, N, Pivots);
    Result := RMatrixLUInverse(A, Pivots, N);
end;


end.