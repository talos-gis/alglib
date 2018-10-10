(*************************************************************************
Copyright (c) 2005-2007, Sergey Bochkanov (ALGLIB project).

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
unit det;
interface
uses Math, Sysutils, Ap, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac;

function RMatrixLUDet(const A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger):Double;
function RMatrixDet(A : TReal2DArray; N : AlglibInteger):Double;

implementation

(*************************************************************************
Determinant calculation of the matrix given by its LU decomposition.

Input parameters:
    A       -   LU decomposition of the matrix (output of
                RMatrixLU subroutine).
    Pivots  -   table of permutations which were made during
                the LU decomposition.
                Output of RMatrixLU subroutine.
    N       -   size of matrix A.

Result: matrix determinant.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
function RMatrixLUDet(const A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger):Double;
var
    I : AlglibInteger;
    S : AlglibInteger;
begin
    Result := 1;
    S := 1;
    I:=0;
    while I<=N-1 do
    begin
        Result := Result*A[I,I];
        if Pivots[I]<>I then
        begin
            S := -S;
        end;
        Inc(I);
    end;
    Result := Result*S;
end;


(*************************************************************************
Calculation of the determinant of a general matrix

Input parameters:
    A       -   matrix, array[0..N-1, 0..N-1]
    N       -   size of matrix A.

Result: determinant of matrix A.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
function RMatrixDet(A : TReal2DArray; N : AlglibInteger):Double;
var
    Pivots : TInteger1DArray;
begin
    A := DynamicArrayCopy(A);
    RMatrixLU(A, N, N, Pivots);
    Result := RMatrixLUDet(A, Pivots, N);
end;


end.