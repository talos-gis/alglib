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
unit hbisinv;
interface
uses Math, Sysutils, Ap, cblas, creflections, hblas, htridiagonal, blas, tdbisinv;

function HMatrixEVDR(A : TComplex2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     B1 : Double;
     B2 : Double;
     var M : AlglibInteger;
     var W : TReal1DArray;
     var Z : TComplex2DArray):Boolean;
function HMatrixEVDI(A : TComplex2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     var W : TReal1DArray;
     var Z : TComplex2DArray):Boolean;
function HermitianEigenValuesAndVectorsInInterval(A : TComplex2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     B1 : Double;
     B2 : Double;
     var M : AlglibInteger;
     var W : TReal1DArray;
     var Z : TComplex2DArray):Boolean;
function HermitianEigenValuesAndVectorsByIndexes(A : TComplex2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     var W : TReal1DArray;
     var Z : TComplex2DArray):Boolean;

implementation

(*************************************************************************
Subroutine for finding the eigenvalues (and eigenvectors) of  a  Hermitian
matrix  in  a  given half-interval (A, B] by using a bisection and inverse
iteration

Input parameters:
    A       -   Hermitian matrix which is given  by  its  upper  or  lower
                triangular  part.  Array  whose   indexes   range   within
                [0..N-1, 0..N-1].
    N       -   size of matrix A.
    ZNeeded -   flag controlling whether the eigenvectors  are  needed  or
                not. If ZNeeded is equal to:
                 * 0, the eigenvectors are not returned;
                 * 1, the eigenvectors are returned.
    IsUpperA -  storage format of matrix A.
    B1, B2 -    half-interval (B1, B2] to search eigenvalues in.

Output parameters:
    M       -   number of eigenvalues found in a given half-interval, M>=0
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

    False, if the bisection method subroutine  wasn't  able  to  find  the
    eigenvalues  in  the  given  interval  or  if  the  inverse  iteration
    subroutine  wasn't  able  to  find all the corresponding eigenvectors.
    In that case, the eigenvalues and eigenvectors are not returned, M  is
    equal to 0.

Note:
    eigen vectors of Hermitian matrix are defined up to multiplication  by
    a complex number L, such as |L|=1.

  -- ALGLIB --
     Copyright 07.01.2006, 24.03.2007 by Bochkanov Sergey.
*************************************************************************)
function HMatrixEVDR(A : TComplex2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     B1 : Double;
     B2 : Double;
     var M : AlglibInteger;
     var W : TReal1DArray;
     var Z : TComplex2DArray):Boolean;
var
    Q : TComplex2DArray;
    T : TReal2DArray;
    Tau : TComplex1DArray;
    E : TReal1DArray;
    WORK : TReal1DArray;
    I : AlglibInteger;
    K : AlglibInteger;
    V : Double;
begin
    A := DynamicArrayCopy(A);
    Assert((ZNeeded=0) or (ZNeeded=1), 'HermitianEigenValuesAndVectorsInInterval: incorrect ZNeeded');
    
    //
    // Reduce to tridiagonal form
    //
    HMatrixTD(A, N, IsUpper, Tau, W, E);
    if ZNeeded=1 then
    begin
        HMatrixTDUnpackQ(A, N, IsUpper, Tau, Q);
        ZNeeded := 2;
    end;
    
    //
    // Bisection and inverse iteration
    //
    Result := SMatrixTDEVDR(W, E, N, ZNeeded, B1, B2, M, T);
    
    //
    // Eigenvectors are needed
    // Calculate Z = Q*T = Re(Q)*T + i*Im(Q)*T
    //
    if Result and (ZNeeded<>0) and (M<>0) then
    begin
        SetLength(WORK, M-1+1);
        SetLength(Z, N-1+1, M-1+1);
        I:=0;
        while I<=N-1 do
        begin
            
            //
            // Calculate real part
            //
            K:=0;
            while K<=M-1 do
            begin
                WORK[K] := 0;
                Inc(K);
            end;
            K:=0;
            while K<=N-1 do
            begin
                V := Q[I,K].X;
                APVAdd(@WORK[0], 0, M-1, @T[K][0], 0, M-1, V);
                Inc(K);
            end;
            K:=0;
            while K<=M-1 do
            begin
                Z[I,K].X := WORK[K];
                Inc(K);
            end;
            
            //
            // Calculate imaginary part
            //
            K:=0;
            while K<=M-1 do
            begin
                WORK[K] := 0;
                Inc(K);
            end;
            K:=0;
            while K<=N-1 do
            begin
                V := Q[I,K].Y;
                APVAdd(@WORK[0], 0, M-1, @T[K][0], 0, M-1, V);
                Inc(K);
            end;
            K:=0;
            while K<=M-1 do
            begin
                Z[I,K].Y := WORK[K];
                Inc(K);
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Subroutine for finding the eigenvalues and  eigenvectors  of  a  Hermitian
matrix with given indexes by using bisection and inverse iteration methods

Input parameters:
    A       -   Hermitian matrix which is given  by  its  upper  or  lower
                triangular part.
                Array whose indexes range within [0..N-1, 0..N-1].
    N       -   size of matrix A.
    ZNeeded -   flag controlling whether the eigenvectors  are  needed  or
                not. If ZNeeded is equal to:
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
                In  that  case,  the eigenvectors are stored in the matrix
                columns.

Result:
    True, if successful. W contains the eigenvalues, Z contains the
    eigenvectors (if needed).

    False, if the bisection method subroutine  wasn't  able  to  find  the
    eigenvalues  in  the  given  interval  or  if  the  inverse  iteration
    subroutine wasn't able to find  all  the  corresponding  eigenvectors.
    In that case, the eigenvalues and eigenvectors are not returned.

Note:
    eigen vectors of Hermitian matrix are defined up to multiplication  by
    a complex number L, such as |L|=1.

  -- ALGLIB --
     Copyright 07.01.2006, 24.03.2007 by Bochkanov Sergey.
*************************************************************************)
function HMatrixEVDI(A : TComplex2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     var W : TReal1DArray;
     var Z : TComplex2DArray):Boolean;
var
    Q : TComplex2DArray;
    T : TReal2DArray;
    Tau : TComplex1DArray;
    E : TReal1DArray;
    WORK : TReal1DArray;
    I : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    M : AlglibInteger;
begin
    A := DynamicArrayCopy(A);
    Assert((ZNeeded=0) or (ZNeeded=1), 'HermitianEigenValuesAndVectorsByIndexes: incorrect ZNeeded');
    
    //
    // Reduce to tridiagonal form
    //
    HMatrixTD(A, N, IsUpper, Tau, W, E);
    if ZNeeded=1 then
    begin
        HMatrixTDUnpackQ(A, N, IsUpper, Tau, Q);
        ZNeeded := 2;
    end;
    
    //
    // Bisection and inverse iteration
    //
    Result := SMatrixTDEVDI(W, E, N, ZNeeded, I1, I2, T);
    
    //
    // Eigenvectors are needed
    // Calculate Z = Q*T = Re(Q)*T + i*Im(Q)*T
    //
    M := I2-I1+1;
    if Result and (ZNeeded<>0) then
    begin
        SetLength(WORK, M-1+1);
        SetLength(Z, N-1+1, M-1+1);
        I:=0;
        while I<=N-1 do
        begin
            
            //
            // Calculate real part
            //
            K:=0;
            while K<=M-1 do
            begin
                WORK[K] := 0;
                Inc(K);
            end;
            K:=0;
            while K<=N-1 do
            begin
                V := Q[I,K].X;
                APVAdd(@WORK[0], 0, M-1, @T[K][0], 0, M-1, V);
                Inc(K);
            end;
            K:=0;
            while K<=M-1 do
            begin
                Z[I,K].X := WORK[K];
                Inc(K);
            end;
            
            //
            // Calculate imaginary part
            //
            K:=0;
            while K<=M-1 do
            begin
                WORK[K] := 0;
                Inc(K);
            end;
            K:=0;
            while K<=N-1 do
            begin
                V := Q[I,K].Y;
                APVAdd(@WORK[0], 0, M-1, @T[K][0], 0, M-1, V);
                Inc(K);
            end;
            K:=0;
            while K<=M-1 do
            begin
                Z[I,K].Y := WORK[K];
                Inc(K);
            end;
            Inc(I);
        end;
    end;
end;


function HermitianEigenValuesAndVectorsInInterval(A : TComplex2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     B1 : Double;
     B2 : Double;
     var M : AlglibInteger;
     var W : TReal1DArray;
     var Z : TComplex2DArray):Boolean;
var
    Q : TComplex2DArray;
    T : TReal2DArray;
    Tau : TComplex1DArray;
    E : TReal1DArray;
    WORK : TReal1DArray;
    I : AlglibInteger;
    K : AlglibInteger;
    V : Double;
begin
    A := DynamicArrayCopy(A);
    Assert((ZNeeded=0) or (ZNeeded=1), 'HermitianEigenValuesAndVectorsInInterval: incorrect ZNeeded');
    
    //
    // Reduce to tridiagonal form
    //
    HermitianToTridiagonal(A, N, IsUpper, Tau, W, E);
    if ZNeeded=1 then
    begin
        UnpackQFromHermitianTridiagonal(A, N, IsUpper, Tau, Q);
        ZNeeded := 2;
    end;
    
    //
    // Bisection and inverse iteration
    //
    Result := TridiagonalEigenValuesAndVectorsInInterval(W, E, N, ZNeeded, B1, B2, M, T);
    
    //
    // Eigenvectors are needed
    // Calculate Z = Q*T = Re(Q)*T + i*Im(Q)*T
    //
    if Result and (ZNeeded<>0) and (M<>0) then
    begin
        SetLength(WORK, M+1);
        SetLength(Z, N+1, M+1);
        I:=1;
        while I<=N do
        begin
            
            //
            // Calculate real part
            //
            K:=1;
            while K<=M do
            begin
                WORK[K] := 0;
                Inc(K);
            end;
            K:=1;
            while K<=N do
            begin
                V := Q[I,K].X;
                APVAdd(@WORK[0], 1, M, @T[K][0], 1, M, V);
                Inc(K);
            end;
            K:=1;
            while K<=M do
            begin
                Z[I,K].X := WORK[K];
                Inc(K);
            end;
            
            //
            // Calculate imaginary part
            //
            K:=1;
            while K<=M do
            begin
                WORK[K] := 0;
                Inc(K);
            end;
            K:=1;
            while K<=N do
            begin
                V := Q[I,K].Y;
                APVAdd(@WORK[0], 1, M, @T[K][0], 1, M, V);
                Inc(K);
            end;
            K:=1;
            while K<=M do
            begin
                Z[I,K].Y := WORK[K];
                Inc(K);
            end;
            Inc(I);
        end;
    end;
end;


function HermitianEigenValuesAndVectorsByIndexes(A : TComplex2DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     IsUpper : Boolean;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     var W : TReal1DArray;
     var Z : TComplex2DArray):Boolean;
var
    Q : TComplex2DArray;
    T : TReal2DArray;
    Tau : TComplex1DArray;
    E : TReal1DArray;
    WORK : TReal1DArray;
    I : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    M : AlglibInteger;
begin
    A := DynamicArrayCopy(A);
    Assert((ZNeeded=0) or (ZNeeded=1), 'HermitianEigenValuesAndVectorsByIndexes: incorrect ZNeeded');
    
    //
    // Reduce to tridiagonal form
    //
    HermitianToTridiagonal(A, N, IsUpper, Tau, W, E);
    if ZNeeded=1 then
    begin
        UnpackQFromHermitianTridiagonal(A, N, IsUpper, Tau, Q);
        ZNeeded := 2;
    end;
    
    //
    // Bisection and inverse iteration
    //
    Result := TridiagonalEigenValuesAndVectorsByIndexes(W, E, N, ZNeeded, I1, I2, T);
    
    //
    // Eigenvectors are needed
    // Calculate Z = Q*T = Re(Q)*T + i*Im(Q)*T
    //
    M := I2-I1+1;
    if Result and (ZNeeded<>0) then
    begin
        SetLength(WORK, M+1);
        SetLength(Z, N+1, M+1);
        I:=1;
        while I<=N do
        begin
            
            //
            // Calculate real part
            //
            K:=1;
            while K<=M do
            begin
                WORK[K] := 0;
                Inc(K);
            end;
            K:=1;
            while K<=N do
            begin
                V := Q[I,K].X;
                APVAdd(@WORK[0], 1, M, @T[K][0], 1, M, V);
                Inc(K);
            end;
            K:=1;
            while K<=M do
            begin
                Z[I,K].X := WORK[K];
                Inc(K);
            end;
            
            //
            // Calculate imaginary part
            //
            K:=1;
            while K<=M do
            begin
                WORK[K] := 0;
                Inc(K);
            end;
            K:=1;
            while K<=N do
            begin
                V := Q[I,K].Y;
                APVAdd(@WORK[0], 1, M, @T[K][0], 1, M, V);
                Inc(K);
            end;
            K:=1;
            while K<=M do
            begin
                Z[I,K].Y := WORK[K];
                Inc(K);
            end;
            Inc(I);
        end;
    end;
end;


end.