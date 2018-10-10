(*************************************************************************
This file is a part of ALGLIB project.

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
unit csolve;
interface
uses Math, Sysutils, Ap, clu;

function CMatrixLUSolve(const A : TComplex2DArray;
     const Pivots : TInteger1DArray;
     B : TComplex1DArray;
     N : AlglibInteger;
     var X : TComplex1DArray):Boolean;
function CMatrixSolve(A : TComplex2DArray;
     B : TComplex1DArray;
     N : AlglibInteger;
     var X : TComplex1DArray):Boolean;
function ComplexSolveSystemLU(const A : TComplex2DArray;
     const Pivots : TInteger1DArray;
     B : TComplex1DArray;
     N : AlglibInteger;
     var X : TComplex1DArray):Boolean;
function ComplexSolveSystem(A : TComplex2DArray;
     B : TComplex1DArray;
     N : AlglibInteger;
     var X : TComplex1DArray):Boolean;

implementation

(*************************************************************************
Solving a system of linear equations with a system matrix given by its
LU decomposition.

The algorithm solves a system of linear equations whose matrix is given by
its LU decomposition. In case of a singular matrix, the algorithm  returns
False.

The algorithm solves systems with a square matrix only.

Input parameters:
    A       -   LU decomposition of a system matrix in compact  form  (the
                result of the RMatrixLU subroutine).
    Pivots  -   row permutation table (the result of a
                RMatrixLU subroutine).
    B       -   right side of a system.
                Array whose index ranges within [0..N-1].
    N       -   size of matrix A.

Output parameters:
    X       -   solution of a system.
                Array whose index ranges within [0..N-1].

Result:
    True, if the matrix is not singular.
    False, if the matrux is singular. In this case, X doesn't contain a
solution.

  -- ALGLIB --
     Copyright 2005-2008 by Bochkanov Sergey
*************************************************************************)
function CMatrixLUSolve(const A : TComplex2DArray;
     const Pivots : TInteger1DArray;
     B : TComplex1DArray;
     N : AlglibInteger;
     var X : TComplex1DArray):Boolean;
var
    Y : TComplex1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Complex;
    i_ : AlglibInteger;
begin
    B := DynamicArrayCopy(B);
    SetLength(Y, N-1+1);
    SetLength(X, N-1+1);
    Result := True;
    I:=0;
    while I<=N-1 do
    begin
        if C_EqualR(A[I,I],0) then
        begin
            Result := False;
            Exit;
        end;
        Inc(I);
    end;
    
    //
    // pivots
    //
    I:=0;
    while I<=N-1 do
    begin
        if Pivots[I]<>I then
        begin
            V := B[I];
            B[I] := B[Pivots[I]];
            B[Pivots[I]] := V;
        end;
        Inc(I);
    end;
    
    //
    // Ly = b
    //
    Y[0] := B[0];
    I:=1;
    while I<=N-1 do
    begin
        V := C_Complex(0.0);
        for i_ := 0 to I-1 do
        begin
            V := C_Add(V,C_Mul(A[I,i_],Y[i_]));
        end;
        Y[I] := C_Sub(B[I],V);
        Inc(I);
    end;
    
    //
    // Ux = y
    //
    X[N-1] := C_Div(Y[N-1],A[N-1,N-1]);
    I:=N-2;
    while I>=0 do
    begin
        V := C_Complex(0.0);
        for i_ := I+1 to N-1 do
        begin
            V := C_Add(V,C_Mul(A[I,i_],X[i_]));
        end;
        X[I] := C_Div(C_Sub(Y[I],V),A[I,I]);
        Dec(I);
    end;
end;


(*************************************************************************
Solving a system of linear equations.

The algorithm solves a system of linear equations by using the
LU decomposition. The algorithm solves systems with a square matrix only.

Input parameters:
    A   -   system matrix.
            Array whose indexes range within [0..N-1, 0..N-1].
    B   -   right side of a system.
            Array whose indexes range within [0..N-1].
    N   -   size of matrix A.

Output parameters:
    X   -   solution of a system.
            Array whose index ranges within [0..N-1].

Result:
    True, if the matrix is not singular.
    False, if the matrix is singular. In this case, X doesn't contain a
solution.

  -- ALGLIB --
     Copyright 2005-2008 by Bochkanov Sergey
*************************************************************************)
function CMatrixSolve(A : TComplex2DArray;
     B : TComplex1DArray;
     N : AlglibInteger;
     var X : TComplex1DArray):Boolean;
var
    Pivots : TInteger1DArray;
    I : AlglibInteger;
begin
    A := DynamicArrayCopy(A);
    B := DynamicArrayCopy(B);
    CMatrixLU(A, N, N, Pivots);
    Result := CMatrixLUSolve(A, Pivots, B, N, X);
end;


function ComplexSolveSystemLU(const A : TComplex2DArray;
     const Pivots : TInteger1DArray;
     B : TComplex1DArray;
     N : AlglibInteger;
     var X : TComplex1DArray):Boolean;
var
    Y : TComplex1DArray;
    I : AlglibInteger;
    V : Complex;
    IP1 : AlglibInteger;
    IM1 : AlglibInteger;
    i_ : AlglibInteger;
begin
    B := DynamicArrayCopy(B);
    SetLength(Y, N+1);
    SetLength(X, N+1);
    Result := True;
    I:=1;
    while I<=N do
    begin
        if C_EqualR(A[I,I],0) then
        begin
            Result := False;
            Exit;
        end;
        Inc(I);
    end;
    
    //
    // pivots
    //
    I:=1;
    while I<=N do
    begin
        if Pivots[I]<>I then
        begin
            V := B[I];
            B[I] := B[Pivots[I]];
            B[Pivots[I]] := V;
        end;
        Inc(I);
    end;
    
    //
    // Ly = b
    //
    Y[1] := B[1];
    I:=2;
    while I<=N do
    begin
        IM1 := I-1;
        V := C_Complex(0.0);
        for i_ := 1 to IM1 do
        begin
            V := C_Add(V,C_Mul(A[I,i_],Y[i_]));
        end;
        Y[I] := C_Sub(B[I],V);
        Inc(I);
    end;
    
    //
    // Ux = y
    //
    X[N] := C_Div(Y[N],A[N,N]);
    I:=N-1;
    while I>=1 do
    begin
        IP1 := I+1;
        V := C_Complex(0.0);
        for i_ := IP1 to N do
        begin
            V := C_Add(V,C_Mul(A[I,i_],X[i_]));
        end;
        X[I] := C_Div(C_Sub(Y[I],V),A[I,I]);
        Dec(I);
    end;
end;


function ComplexSolveSystem(A : TComplex2DArray;
     B : TComplex1DArray;
     N : AlglibInteger;
     var X : TComplex1DArray):Boolean;
var
    Pivots : TInteger1DArray;
begin
    A := DynamicArrayCopy(A);
    B := DynamicArrayCopy(B);
    ComplexLUDecomposition(A, N, N, Pivots);
    Result := ComplexSolveSystemLU(A, Pivots, B, N, X);
end;


end.