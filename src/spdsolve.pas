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
unit spdsolve;
interface
uses Math, Sysutils, Ap, cholesky;

function SPDMatrixCholeskySolve(const A : TReal2DArray;
     B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;
function SPDMatrixSolve(A : TReal2DArray;
     B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;
function SolveSystemCholesky(const A : TReal2DArray;
     B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;
function SolveSPDSystem(A : TReal2DArray;
     B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;

implementation

(*************************************************************************
Solving a system of linear equations with a system  matrix  given  by  its
Cholesky decomposition.

The algorithm solves systems with a square matrix only.

Input parameters:
    A       -   Cholesky decomposition of a system matrix (the result of
                the SMatrixCholesky subroutine).
    B       -   right side of a system.
                Array whose index ranges within [0..N-1].
    N       -   size of matrix A.
    IsUpper -   points to the triangle of matrix A in which the Cholesky
                decomposition is stored. If IsUpper=True,  the  Cholesky
                decomposition has the form of U'*U, and the upper triangle
                of matrix A stores matrix U (in  that  case,  the  lower
                triangle isn’t used and isn’t changed by the subroutine)
                Similarly, if IsUpper = False, the Cholesky decomposition
                has the form of L*L',  and  the  lower  triangle  stores
                matrix L.

Output parameters:
    X       -   solution of a system.
                Array whose index ranges within [0..N-1].

Result:
    True, if the system is not singular. X contains the solution.
    False, if the system is singular (there is a zero element on the main
diagonal). In this case, X doesn't contain a solution.

  -- ALGLIB --
     Copyright 2005-2008 by Bochkanov Sergey
*************************************************************************)
function SPDMatrixCholeskySolve(const A : TReal2DArray;
     B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;
var
    I : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
begin
    B := DynamicArrayCopy(B);
    Assert(N>0, 'Error: N<=0 in SolveSystemCholesky');
    
    //
    // det(A)=0?
    //
    Result := True;
    I:=0;
    while I<=N-1 do
    begin
        if AP_FP_Eq(A[I,I],0) then
        begin
            Result := False;
            Exit;
        end;
        Inc(I);
    end;
    
    //
    // det(A)<>0
    //
    SetLength(X, N-1+1);
    if IsUpper then
    begin
        
        //
        // A = U'*U, solve U'*y = b first
        //
        B[0] := B[0]/A[0,0];
        I:=1;
        while I<=N-1 do
        begin
            V := 0.0;
            for i_ := 0 to I-1 do
            begin
                V := V + A[i_,I]*B[i_];
            end;
            B[I] := (B[I]-V)/A[I,I];
            Inc(I);
        end;
        
        //
        // Solve U*x = y
        //
        B[N-1] := B[N-1]/A[N-1,N-1];
        I:=N-2;
        while I>=0 do
        begin
            V := APVDotProduct(@A[I][0], I+1, N-1, @B[0], I+1, N-1);
            B[I] := (B[I]-V)/A[I,I];
            Dec(I);
        end;
        APVMove(@X[0], 0, N-1, @B[0], 0, N-1);
    end
    else
    begin
        
        //
        // A = L*L', solve L'*y = b first
        //
        B[0] := B[0]/A[0,0];
        I:=1;
        while I<=N-1 do
        begin
            V := APVDotProduct(@A[I][0], 0, I-1, @B[0], 0, I-1);
            B[I] := (B[I]-V)/A[I,I];
            Inc(I);
        end;
        
        //
        // Solve L'*x = y
        //
        B[N-1] := B[N-1]/A[N-1,N-1];
        I:=N-2;
        while I>=0 do
        begin
            V := 0.0;
            for i_ := I+1 to N-1 do
            begin
                V := V + A[i_,I]*B[i_];
            end;
            B[I] := (B[I]-V)/A[I,I];
            Dec(I);
        end;
        APVMove(@X[0], 0, N-1, @B[0], 0, N-1);
    end;
end;


(*************************************************************************
Solving a system of linear equations with  a  symmetric  positive-definite
matrix by using the Cholesky decomposition.

The algorithm solves a system of linear equations whose matrix is symmetric
and positive-definite.

Input parameters:
    A       -   upper or lower triangle part of a symmetric system matrix.
                Array whose indexes range within [0..N-1, 0..N-1].
    B       -   right side of a system.
                Array whose index ranges within [0..N-1].
    N       -   size of matrix A.
    IsUpper -   points to the triangle of matrix A in which the matrix is stored.

Output parameters:
    X       -   solution of a system.
                Array whose index ranges within [0..N-1].

Result:
    True, if the system is not singular.
    False, if the system is singular. In this case, X doesn't contain a
solution.

  -- ALGLIB --
     Copyright 2005-2008 by Bochkanov Sergey
*************************************************************************)
function SPDMatrixSolve(A : TReal2DArray;
     B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;
begin
    A := DynamicArrayCopy(A);
    B := DynamicArrayCopy(B);
    Result := SPDMatrixCholesky(A, N, IsUpper);
    if  not Result then
    begin
        Exit;
    end;
    Result := SPDMatrixCholeskySolve(A, B, N, IsUpper, X);
end;


function SolveSystemCholesky(const A : TReal2DArray;
     B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;
var
    I : AlglibInteger;
    IM1 : AlglibInteger;
    IP1 : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
begin
    B := DynamicArrayCopy(B);
    Assert(N>0, 'Error: N<=0 in SolveSystemCholesky');
    
    //
    // det(A)=0?
    //
    Result := True;
    I:=1;
    while I<=N do
    begin
        if AP_FP_Eq(A[I,I],0) then
        begin
            Result := False;
            Exit;
        end;
        Inc(I);
    end;
    
    //
    // det(A)<>0
    //
    SetLength(X, N+1);
    if IsUpper then
    begin
        
        //
        // A = U'*U, solve U'*y = b first
        //
        B[1] := B[1]/A[1,1];
        I:=2;
        while I<=N do
        begin
            IM1 := I-1;
            V := 0.0;
            for i_ := 1 to IM1 do
            begin
                V := V + A[i_,I]*B[i_];
            end;
            B[I] := (B[I]-V)/A[I,I];
            Inc(I);
        end;
        
        //
        // Solve U*x = y
        //
        B[N] := B[N]/A[N,N];
        I:=N-1;
        while I>=1 do
        begin
            IP1 := I+1;
            V := APVDotProduct(@A[I][0], IP1, N, @B[0], IP1, N);
            B[I] := (B[I]-V)/A[I,I];
            Dec(I);
        end;
        APVMove(@X[0], 1, N, @B[0], 1, N);
    end
    else
    begin
        
        //
        // A = L*L', solve L'*y = b first
        //
        B[1] := B[1]/A[1,1];
        I:=2;
        while I<=N do
        begin
            IM1 := I-1;
            V := APVDotProduct(@A[I][0], 1, IM1, @B[0], 1, IM1);
            B[I] := (B[I]-V)/A[I,I];
            Inc(I);
        end;
        
        //
        // Solve L'*x = y
        //
        B[N] := B[N]/A[N,N];
        I:=N-1;
        while I>=1 do
        begin
            IP1 := I+1;
            V := 0.0;
            for i_ := IP1 to N do
            begin
                V := V + A[i_,I]*B[i_];
            end;
            B[I] := (B[I]-V)/A[I,I];
            Dec(I);
        end;
        APVMove(@X[0], 1, N, @B[0], 1, N);
    end;
end;


function SolveSPDSystem(A : TReal2DArray;
     B : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var X : TReal1DArray):Boolean;
begin
    A := DynamicArrayCopy(A);
    B := DynamicArrayCopy(B);
    Result := CholeskyDecomposition(A, N, IsUpper);
    if  not Result then
    begin
        Exit;
    end;
    Result := SolveSystemCholesky(A, B, N, IsUpper, X);
end;


end.