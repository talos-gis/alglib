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
unit sbisinv;
interface
uses Math, Sysutils, Ap, sblas, reflections, tridiagonal, blas, tdbisinv;

function SMatrixEVDR(A : TReal2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     B1 : Double;
     B2 : Double;
     var M : AlglibInteger;
     var W : TReal1DArray;
     var Z : TReal2DArray):Boolean;
function SMatrixEVDI(A : TReal2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     var W : TReal1DArray;
     var Z : TReal2DArray):Boolean;
function SymmetricEigenValuesAndVectorsInInterval(A : TReal2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     B1 : Double;
     B2 : Double;
     var M : AlglibInteger;
     var W : TReal1DArray;
     var Z : TReal2DArray):Boolean;
function SymmetricEigenValuesAndVectorsByIndexes(A : TReal2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     var W : TReal1DArray;
     var Z : TReal2DArray):Boolean;

implementation

(*************************************************************************
Subroutine for finding the eigenvalues (and eigenvectors) of  a  symmetric
matrix  in  a  given half open interval (A, B] by using  a  bisection  and
inverse iteration

Input parameters:
    A       -   symmetric matrix which is given by its upper or lower
                triangular part. Array [0..N-1, 0..N-1].
    N       -   size of matrix A.
    ZNeeded -   flag controlling whether the eigenvectors are needed or not.
                If ZNeeded is equal to:
                 * 0, the eigenvectors are not returned;
                 * 1, the eigenvectors are returned.
    IsUpperA -  storage format of matrix A.
    B1, B2 -    half open interval (B1, B2] to search eigenvalues in.

Output parameters:
    M       -   number of eigenvalues found in a given half-interval (M>=0).
    W       -   array of the eigenvalues found.
                Array whose index ranges within [0..M-1].
    Z       -   if ZNeeded is equal to:
                 * 0, Z hasn’t changed;
                 * 1, Z contains eigenvectors.
                Array whose indexes range within [0..N-1, 0..M-1].
                The eigenvectors are stored in the matrix columns.

Result:
    True, if successful. M contains the number of eigenvalues in the given
    half-interval (could be equal to 0), W contains the eigenvalues,
    Z contains the eigenvectors (if needed).

    False, if the bisection method subroutine wasn't able to find the
    eigenvalues in the given interval or if the inverse iteration subroutine
    wasn't able to find all the corresponding eigenvectors.
    In that case, the eigenvalues and eigenvectors are not returned,
    M is equal to 0.

  -- ALGLIB --
     Copyright 07.01.2006 by Bochkanov Sergey
*************************************************************************)
function SMatrixEVDR(A : TReal2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     B1 : Double;
     B2 : Double;
     var M : AlglibInteger;
     var W : TReal1DArray;
     var Z : TReal2DArray):Boolean;
var
    Tau : TReal1DArray;
    E : TReal1DArray;
begin
    A := DynamicArrayCopy(A);
    Assert((ZNeeded=0) or (ZNeeded=1), 'SMatrixTDEVDR: incorrect ZNeeded');
    SMatrixTD(A, N, IsUpper, Tau, W, E);
    if ZNeeded=1 then
    begin
        SMatrixTDUnpackQ(A, N, IsUpper, Tau, Z);
    end;
    Result := SMatrixTDEVDR(W, E, N, ZNeeded, B1, B2, M, Z);
end;


(*************************************************************************
Subroutine for finding the eigenvalues and  eigenvectors  of  a  symmetric
matrix with given indexes by using bisection and inverse iteration methods.

Input parameters:
    A       -   symmetric matrix which is given by its upper or lower
                triangular part. Array whose indexes range within [0..N-1, 0..N-1].
    N       -   size of matrix A.
    ZNeeded -   flag controlling whether the eigenvectors are needed or not.
                If ZNeeded is equal to:
                 * 0, the eigenvectors are not returned;
                 * 1, the eigenvectors are returned.
    IsUpperA -  storage format of matrix A.
    I1, I2 -    index interval for searching (from I1 to I2).
                0 <= I1 <= I2 <= N-1.

Output parameters:
    W       -   array of the eigenvalues found.
                Array whose index ranges within [0..I2-I1].
    Z       -   if ZNeeded is equal to:
                 * 0, Z hasn’t changed;
                 * 1, Z contains eigenvectors.
                Array whose indexes range within [0..N-1, 0..I2-I1].
                In that case, the eigenvectors are stored in the matrix columns.

Result:
    True, if successful. W contains the eigenvalues, Z contains the
    eigenvectors (if needed).

    False, if the bisection method subroutine wasn't able to find the
    eigenvalues in the given interval or if the inverse iteration subroutine
    wasn't able to find all the corresponding eigenvectors.
    In that case, the eigenvalues and eigenvectors are not returned.

  -- ALGLIB --
     Copyright 07.01.2006 by Bochkanov Sergey
*************************************************************************)
function SMatrixEVDI(A : TReal2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     var W : TReal1DArray;
     var Z : TReal2DArray):Boolean;
var
    Tau : TReal1DArray;
    E : TReal1DArray;
begin
    A := DynamicArrayCopy(A);
    Assert((ZNeeded=0) or (ZNeeded=1), 'SMatrixEVDI: incorrect ZNeeded');
    SMatrixTD(A, N, IsUpper, Tau, W, E);
    if ZNeeded=1 then
    begin
        SMatrixTDUnpackQ(A, N, IsUpper, Tau, Z);
    end;
    Result := SMatrixTDEVDI(W, E, N, ZNeeded, I1, I2, Z);
end;


function SymmetricEigenValuesAndVectorsInInterval(A : TReal2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     B1 : Double;
     B2 : Double;
     var M : AlglibInteger;
     var W : TReal1DArray;
     var Z : TReal2DArray):Boolean;
var
    Tau : TReal1DArray;
    E : TReal1DArray;
begin
    A := DynamicArrayCopy(A);
    Assert((ZNeeded=0) or (ZNeeded=1), 'SymmetricEigenValuesAndVectorsInInterval: incorrect ZNeeded');
    ToTridiagonal(A, N, IsUpper, Tau, W, E);
    if ZNeeded=1 then
    begin
        UnpackQFromTridiagonal(A, N, IsUpper, Tau, Z);
    end;
    Result := TridiagonalEigenValuesAndVectorsInInterval(W, E, N, ZNeeded, B1, B2, M, Z);
end;


function SymmetricEigenValuesAndVectorsByIndexes(A : TReal2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     var W : TReal1DArray;
     var Z : TReal2DArray):Boolean;
var
    Tau : TReal1DArray;
    E : TReal1DArray;
begin
    A := DynamicArrayCopy(A);
    Assert((ZNeeded=0) or (ZNeeded=1), 'SymmetricEigenValuesAndVectorsInInterval: incorrect ZNeeded');
    ToTridiagonal(A, N, IsUpper, Tau, W, E);
    if ZNeeded=1 then
    begin
        UnpackQFromTridiagonal(A, N, IsUpper, Tau, Z);
    end;
    Result := TridiagonalEigenValuesAndVectorsByIndexes(W, E, N, ZNeeded, I1, I2, Z);
end;


end.