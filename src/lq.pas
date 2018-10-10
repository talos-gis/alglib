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
unit lq;
interface
uses Math, Sysutils, Ap, reflections;

procedure RMatrixLQ(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TReal1DArray);
procedure RMatrixLQUnpackQ(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     QRows : AlglibInteger;
     var Q : TReal2DArray);
procedure RMatrixLQUnpackL(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var L : TReal2DArray);
procedure LQDecomposition(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TReal1DArray);
procedure UnpackQFromLQ(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     QRows : AlglibInteger;
     var Q : TReal2DArray);
procedure LQDecompositionUnpacked(A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var L : TReal2DArray;
     var Q : TReal2DArray);

implementation

(*************************************************************************
LQ decomposition of a rectangular matrix of size MxN

Input parameters:
    A   -   matrix A whose indexes range within [0..M-1, 0..N-1].
    M   -   number of rows in matrix A.
    N   -   number of columns in matrix A.

Output parameters:
    A   -   matrices L and Q in compact form (see below)
    Tau -   array of scalar factors which are used to form
            matrix Q. Array whose index ranges within [0..Min(M,N)-1].

Matrix A is represented as A = LQ, where Q is an orthogonal matrix of size
MxM, L - lower triangular (or lower trapezoid) matrix of size M x N.

The elements of matrix L are located on and below  the  main  diagonal  of
matrix A. The elements which are located in Tau array and above  the  main
diagonal of matrix A are used to form matrix Q as follows:

Matrix Q is represented as a product of elementary reflections

Q = H(k-1)*H(k-2)*...*H(1)*H(0),

where k = min(m,n), and each H(i) is of the form

H(i) = 1 - tau * v * (v^T)

where tau is a scalar stored in Tau[I]; v - real vector, so that v(0:i-1)=0,
v(i) = 1, v(i+1:n-1) stored in A(i,i+1:n-1).

  -- ALGLIB --
     Copyright 2005-2007 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixLQ(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TReal1DArray);
var
    WORK : TReal1DArray;
    T : TReal1DArray;
    I : AlglibInteger;
    K : AlglibInteger;
    MinMN : AlglibInteger;
    MaxMN : AlglibInteger;
    Tmp : Double;
begin
    MinMN := Min(M, N);
    MaxMN := Max(M, N);
    SetLength(WORK, M+1);
    SetLength(T, N+1);
    SetLength(TAU, MinMN-1+1);
    K := Min(M, N);
    I:=0;
    while I<=K-1 do
    begin
        
        //
        // Generate elementary reflector H(i) to annihilate A(i,i+1:n-1)
        //
        APVMove(@T[0], 1, N-I, @A[I][0], I, N-1);
        GenerateReflection(T, N-I, Tmp);
        Tau[I] := Tmp;
        APVMove(@A[I][0], I, N-1, @T[0], 1, N-I);
        T[1] := 1;
        if I<N then
        begin
            
            //
            // Apply H(i) to A(i+1:m,i:n) from the right
            //
            ApplyReflectionFromTheRight(A, Tau[I], T, I+1, M-1, I, N-1, WORK);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Partial unpacking of matrix Q from the LQ decomposition of a matrix A

Input parameters:
    A       -   matrices L and Q in compact form.
                Output of RMatrixLQ subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.
    Tau     -   scalar factors which are used to form Q.
                Output of the RMatrixLQ subroutine.
    QRows   -   required number of rows in matrix Q. N>=QRows>=0.

Output parameters:
    Q       -   first QRows rows of matrix Q. Array whose indexes range
                within [0..QRows-1, 0..N-1]. If QRows=0, the array remains
                unchanged.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixLQUnpackQ(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     QRows : AlglibInteger;
     var Q : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    MinMN : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
begin
    Assert(QRows<=N, 'RMatrixLQUnpackQ: QRows>N!');
    if (M<=0) or (N<=0) or (QRows<=0) then
    begin
        Exit;
    end;
    
    //
    // init
    //
    MinMN := Min(M, N);
    K := Min(MinMN, QRows);
    SetLength(Q, QRows-1+1, N-1+1);
    SetLength(V, N+1);
    SetLength(WORK, QRows+1);
    I:=0;
    while I<=QRows-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if I=J then
            begin
                Q[I,J] := 1;
            end
            else
            begin
                Q[I,J] := 0;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // unpack Q
    //
    I:=K-1;
    while I>=0 do
    begin
        
        //
        // Apply H(i)
        //
        APVMove(@V[0], 1, N-I, @A[I][0], I, N-1);
        V[1] := 1;
        ApplyReflectionFromTheRight(Q, Tau[I], V, 0, QRows-1, I, N-1, WORK);
        Dec(I);
    end;
end;


(*************************************************************************
Unpacking of matrix L from the LQ decomposition of a matrix A

Input parameters:
    A       -   matrices Q and L in compact form.
                Output of RMatrixLQ subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.

Output parameters:
    L       -   matrix L, array[0..M-1, 0..N-1].

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixLQUnpackL(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var L : TReal2DArray);
var
    I : AlglibInteger;
    K : AlglibInteger;
begin
    if (M<=0) or (N<=0) then
    begin
        Exit;
    end;
    SetLength(L, M-1+1, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        L[0,I] := 0;
        Inc(I);
    end;
    I:=1;
    while I<=M-1 do
    begin
        APVMove(@L[I][0], 0, N-1, @L[0][0], 0, N-1);
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        K := Min(I, N-1);
        APVMove(@L[I][0], 0, K, @A[I][0], 0, K);
        Inc(I);
    end;
end;


procedure LQDecomposition(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TReal1DArray);
var
    WORK : TReal1DArray;
    T : TReal1DArray;
    I : AlglibInteger;
    K : AlglibInteger;
    NMIP1 : AlglibInteger;
    MinMN : AlglibInteger;
    MaxMN : AlglibInteger;
    Tmp : Double;
begin
    MinMN := Min(M, N);
    MaxMN := Max(M, N);
    SetLength(WORK, M+1);
    SetLength(T, N+1);
    SetLength(TAU, MinMN+1);
    
    //
    // Test the input arguments
    //
    K := Min(M, N);
    I:=1;
    while I<=K do
    begin
        
        //
        // Generate elementary reflector H(i) to annihilate A(i,i+1:n)
        //
        NMIP1 := N-I+1;
        APVMove(@T[0], 1, NMIP1, @A[I][0], I, N);
        GenerateReflection(T, NMIP1, Tmp);
        Tau[I] := Tmp;
        APVMove(@A[I][0], I, N, @T[0], 1, NMIP1);
        T[1] := 1;
        if I<N then
        begin
            
            //
            // Apply H(i) to A(i+1:m,i:n) from the right
            //
            ApplyReflectionFromTheRight(A, Tau[I], T, I+1, M, I, N, WORK);
        end;
        Inc(I);
    end;
end;


procedure UnpackQFromLQ(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     QRows : AlglibInteger;
     var Q : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    MinMN : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
    VM : AlglibInteger;
begin
    Assert(QRows<=N, 'UnpackQFromLQ: QRows>N!');
    if (M=0) or (N=0) or (QRows=0) then
    begin
        Exit;
    end;
    
    //
    // init
    //
    MinMN := Min(M, N);
    K := Min(MinMN, QRows);
    SetLength(Q, QRows+1, N+1);
    SetLength(V, N+1);
    SetLength(WORK, QRows+1);
    I:=1;
    while I<=QRows do
    begin
        J:=1;
        while J<=N do
        begin
            if I=J then
            begin
                Q[I,J] := 1;
            end
            else
            begin
                Q[I,J] := 0;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // unpack Q
    //
    I:=K;
    while I>=1 do
    begin
        
        //
        // Apply H(i)
        //
        VM := N-I+1;
        APVMove(@V[0], 1, VM, @A[I][0], I, N);
        V[1] := 1;
        ApplyReflectionFromTheRight(Q, Tau[I], V, 1, QRows, I, N, WORK);
        Dec(I);
    end;
end;


procedure LQDecompositionUnpacked(A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var L : TReal2DArray;
     var Q : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    Tau : TReal1DArray;
begin
    A := DynamicArrayCopy(A);
    if N<=0 then
    begin
        Exit;
    end;
    SetLength(Q, N+1, N+1);
    SetLength(L, M+1, N+1);
    
    //
    // LQDecomposition
    //
    LQDecomposition(A, M, N, Tau);
    
    //
    // L
    //
    I:=1;
    while I<=M do
    begin
        J:=1;
        while J<=N do
        begin
            if J>I then
            begin
                L[I,J] := 0;
            end
            else
            begin
                L[I,J] := A[I,J];
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Q
    //
    UnpackQFromLQ(A, M, N, Tau, N, Q);
end;


end.