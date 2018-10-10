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
unit hevd;
interface
uses Math, Sysutils, Ap, blas, rotations, tdevd, cblas, creflections, hblas, htridiagonal;

function HMatrixEVD(A : TComplex2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     var D : TReal1DArray;
     var Z : TComplex2DArray):Boolean;
function HermitianEVD(A : TComplex2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     var D : TReal1DArray;
     var Z : TComplex2DArray):Boolean;

implementation

(*************************************************************************
Finding the eigenvalues and eigenvectors of a Hermitian matrix

The algorithm finds eigen pairs of a Hermitian matrix by  reducing  it  to
real tridiagonal form and using the QL/QR algorithm.

Input parameters:
    A       -   Hermitian matrix which is given  by  its  upper  or  lower
                triangular part.
                Array whose indexes range within [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format.
    ZNeeded -   flag controlling whether the eigenvectors  are  needed  or
                not. If ZNeeded is equal to:
                 * 0, the eigenvectors are not returned;
                 * 1, the eigenvectors are returned.

Output parameters:
    D       -   eigenvalues in ascending order.
                Array whose index ranges within [0..N-1].
    Z       -   if ZNeeded is equal to:
                 * 0, Z hasn’t changed;
                 * 1, Z contains the eigenvectors.
                Array whose indexes range within [0..N-1, 0..N-1].
                The eigenvectors are stored in the matrix columns.

Result:
    True, if the algorithm has converged.
    False, if the algorithm hasn't converged (rare case).

Note:
    eigen vectors of Hermitian matrix are defined up to multiplication  by
    a complex number L, such as |L|=1.

  -- ALGLIB --
     Copyright 2005, 23 March 2007 by Bochkanov Sergey
*************************************************************************)
function HMatrixEVD(A : TComplex2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     var D : TReal1DArray;
     var Z : TComplex2DArray):Boolean;
var
    Tau : TComplex1DArray;
    E : TReal1DArray;
    WORK : TReal1DArray;
    T : TReal2DArray;
    Q : TComplex2DArray;
    I : AlglibInteger;
    K : AlglibInteger;
    V : Double;
begin
    A := DynamicArrayCopy(A);
    Assert((ZNeeded=0) or (ZNeeded=1), 'HermitianEVD: incorrect ZNeeded');
    
    //
    // Reduce to tridiagonal form
    //
    HMatrixTD(A, N, IsUpper, Tau, D, E);
    if ZNeeded=1 then
    begin
        HMatrixTDUnpackQ(A, N, IsUpper, Tau, Q);
        ZNeeded := 2;
    end;
    
    //
    // TDEVD
    //
    Result := SMatrixTDEVD(D, E, N, ZNeeded, T);
    
    //
    // Eigenvectors are needed
    // Calculate Z = Q*T = Re(Q)*T + i*Im(Q)*T
    //
    if Result and (ZNeeded<>0) then
    begin
        SetLength(WORK, N-1+1);
        SetLength(Z, N-1+1, N-1+1);
        I:=0;
        while I<=N-1 do
        begin
            
            //
            // Calculate real part
            //
            K:=0;
            while K<=N-1 do
            begin
                WORK[K] := 0;
                Inc(K);
            end;
            K:=0;
            while K<=N-1 do
            begin
                V := Q[I,K].X;
                APVAdd(@WORK[0], 0, N-1, @T[K][0], 0, N-1, V);
                Inc(K);
            end;
            K:=0;
            while K<=N-1 do
            begin
                Z[I,K].X := WORK[K];
                Inc(K);
            end;
            
            //
            // Calculate imaginary part
            //
            K:=0;
            while K<=N-1 do
            begin
                WORK[K] := 0;
                Inc(K);
            end;
            K:=0;
            while K<=N-1 do
            begin
                V := Q[I,K].Y;
                APVAdd(@WORK[0], 0, N-1, @T[K][0], 0, N-1, V);
                Inc(K);
            end;
            K:=0;
            while K<=N-1 do
            begin
                Z[I,K].Y := WORK[K];
                Inc(K);
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************

  -- ALGLIB --
     Copyright 2005, 23 March 2007 by Bochkanov Sergey
*************************************************************************)
function HermitianEVD(A : TComplex2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     var D : TReal1DArray;
     var Z : TComplex2DArray):Boolean;
var
    Tau : TComplex1DArray;
    E : TReal1DArray;
    WORK : TReal1DArray;
    T : TReal2DArray;
    Q : TComplex2DArray;
    I : AlglibInteger;
    K : AlglibInteger;
    V : Double;
begin
    A := DynamicArrayCopy(A);
    Assert((ZNeeded=0) or (ZNeeded=1), 'HermitianEVD: incorrect ZNeeded');
    
    //
    // Reduce to tridiagonal form
    //
    HermitianToTridiagonal(A, N, IsUpper, Tau, D, E);
    if ZNeeded=1 then
    begin
        UnpackQFromHermitianTridiagonal(A, N, IsUpper, Tau, Q);
        ZNeeded := 2;
    end;
    
    //
    // TDEVD
    //
    Result := TridiagonalEVD(D, E, N, ZNeeded, T);
    
    //
    // Eigenvectors are needed
    // Calculate Z = Q*T = Re(Q)*T + i*Im(Q)*T
    //
    if Result and (ZNeeded<>0) then
    begin
        SetLength(WORK, N+1);
        SetLength(Z, N+1, N+1);
        I:=1;
        while I<=N do
        begin
            
            //
            // Calculate real part
            //
            K:=1;
            while K<=N do
            begin
                WORK[K] := 0;
                Inc(K);
            end;
            K:=1;
            while K<=N do
            begin
                V := Q[I,K].X;
                APVAdd(@WORK[0], 1, N, @T[K][0], 1, N, V);
                Inc(K);
            end;
            K:=1;
            while K<=N do
            begin
                Z[I,K].X := WORK[K];
                Inc(K);
            end;
            
            //
            // Calculate imaginary part
            //
            K:=1;
            while K<=N do
            begin
                WORK[K] := 0;
                Inc(K);
            end;
            K:=1;
            while K<=N do
            begin
                V := Q[I,K].Y;
                APVAdd(@WORK[0], 1, N, @T[K][0], 1, N, V);
                Inc(K);
            end;
            K:=1;
            while K<=N do
            begin
                Z[I,K].Y := WORK[K];
                Inc(K);
            end;
            Inc(I);
        end;
    end;
end;


end.