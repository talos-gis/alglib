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
unit hessenberg;
interface
uses Math, Sysutils, Ap, reflections;

procedure RMatrixHessenberg(var A : TReal2DArray;
     N : AlglibInteger;
     var Tau : TReal1DArray);
procedure RMatrixHessenbergUnpackQ(const A : TReal2DArray;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     var Q : TReal2DArray);
procedure RMatrixHessenbergUnpackH(const A : TReal2DArray;
     N : AlglibInteger;
     var H : TReal2DArray);
procedure ToUpperHessenberg(var A : TReal2DArray;
     N : AlglibInteger;
     var Tau : TReal1DArray);
procedure UnpackQFromUpperHessenberg(const A : TReal2DArray;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     var Q : TReal2DArray);
procedure UnpackHFromUpperHessenberg(const A : TReal2DArray;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     var H : TReal2DArray);

implementation

(*************************************************************************
Reduction of a square matrix to  upper Hessenberg form: Q'*A*Q = H,
where Q is an orthogonal matrix, H - Hessenberg matrix.

Input parameters:
    A       -   matrix A with elements [0..N-1, 0..N-1]
    N       -   size of matrix A.

Output parameters:
    A       -   matrices Q and P in  compact form (see below).
    Tau     -   array of scalar factors which are used to form matrix Q.
                Array whose index ranges within [0..N-2]

Matrix H is located on the main diagonal, on the lower secondary  diagonal
and above the main diagonal of matrix A. The elements which are used to
form matrix Q are situated in array Tau and below the lower secondary
diagonal of matrix A as follows:

Matrix Q is represented as a product of elementary reflections

Q = H(0)*H(2)*...*H(n-2),

where each H(i) is given by

H(i) = 1 - tau * v * (v^T)

where tau is a scalar stored in Tau[I]; v - is a real vector,
so that v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) stored in A(i+2:n-1,i).

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1992
*************************************************************************)
procedure RMatrixHessenberg(var A : TReal2DArray;
     N : AlglibInteger;
     var Tau : TReal1DArray);
var
    I : AlglibInteger;
    V : Double;
    T : TReal1DArray;
    WORK : TReal1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert(N>=0, 'RMatrixHessenberg: incorrect N!');
    
    //
    // Quick return if possible
    //
    if N<=1 then
    begin
        Exit;
    end;
    SetLength(Tau, N-2+1);
    SetLength(T, N+1);
    SetLength(WORK, N-1+1);
    I:=0;
    while I<=N-2 do
    begin
        
        //
        // Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
        //
        i1_ := (I+1) - (1);
        for i_ := 1 to N-I-1 do
        begin
            T[i_] := A[i_+i1_,I];
        end;
        GenerateReflection(T, N-I-1, V);
        i1_ := (1) - (I+1);
        for i_ := I+1 to N-1 do
        begin
            A[i_,I] := T[i_+i1_];
        end;
        Tau[I] := V;
        T[1] := 1;
        
        //
        // Apply H(i) to A(1:ihi,i+1:ihi) from the right
        //
        ApplyReflectionFromTheRight(A, V, T, 0, N-1, I+1, N-1, WORK);
        
        //
        // Apply H(i) to A(i+1:ihi,i+1:n) from the left
        //
        ApplyReflectionFromTheLeft(A, V, T, I+1, N-1, I+1, N-1, WORK);
        Inc(I);
    end;
end;


(*************************************************************************
Unpacking matrix Q which reduces matrix A to upper Hessenberg form

Input parameters:
    A   -   output of RMatrixHessenberg subroutine.
    N   -   size of matrix A.
    Tau -   scalar factors which are used to form Q.
            Output of RMatrixHessenberg subroutine.

Output parameters:
    Q   -   matrix Q.
            Array whose indexes range within [0..N-1, 0..N-1].

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixHessenbergUnpackQ(const A : TReal2DArray;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     var Q : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if N=0 then
    begin
        Exit;
    end;
    
    //
    // init
    //
    SetLength(Q, N-1+1, N-1+1);
    SetLength(V, N-1+1);
    SetLength(WORK, N-1+1);
    I:=0;
    while I<=N-1 do
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
    I:=0;
    while I<=N-2 do
    begin
        
        //
        // Apply H(i)
        //
        i1_ := (I+1) - (1);
        for i_ := 1 to N-I-1 do
        begin
            V[i_] := A[i_+i1_,I];
        end;
        V[1] := 1;
        ApplyReflectionFromTheRight(Q, Tau[I], V, 0, N-1, I+1, N-1, WORK);
        Inc(I);
    end;
end;


(*************************************************************************
Unpacking matrix H (the result of matrix A reduction to upper Hessenberg form)

Input parameters:
    A   -   output of RMatrixHessenberg subroutine.
    N   -   size of matrix A.

Output parameters:
    H   -   matrix H. Array whose indexes range within [0..N-1, 0..N-1].

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixHessenbergUnpackH(const A : TReal2DArray;
     N : AlglibInteger;
     var H : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
begin
    if N=0 then
    begin
        Exit;
    end;
    SetLength(H, N-1+1, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=I-2 do
        begin
            H[I,J] := 0;
            Inc(J);
        end;
        J := Max(0, I-1);
        APVMove(@H[I][0], J, N-1, @A[I][0], J, N-1);
        Inc(I);
    end;
end;


procedure ToUpperHessenberg(var A : TReal2DArray;
     N : AlglibInteger;
     var Tau : TReal1DArray);
var
    I : AlglibInteger;
    IP1 : AlglibInteger;
    NMI : AlglibInteger;
    V : Double;
    T : TReal1DArray;
    WORK : TReal1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert(N>=0, 'ToUpperHessenberg: incorrect N!');
    
    //
    // Quick return if possible
    //
    if N<=1 then
    begin
        Exit;
    end;
    SetLength(Tau, N-1+1);
    SetLength(T, N+1);
    SetLength(WORK, N+1);
    I:=1;
    while I<=N-1 do
    begin
        
        //
        // Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
        //
        IP1 := I+1;
        NMI := N-I;
        i1_ := (IP1) - (1);
        for i_ := 1 to NMI do
        begin
            T[i_] := A[i_+i1_,I];
        end;
        GenerateReflection(T, NMI, V);
        i1_ := (1) - (IP1);
        for i_ := IP1 to N do
        begin
            A[i_,I] := T[i_+i1_];
        end;
        Tau[I] := V;
        T[1] := 1;
        
        //
        // Apply H(i) to A(1:ihi,i+1:ihi) from the right
        //
        ApplyReflectionFromTheRight(A, V, T, 1, N, I+1, N, WORK);
        
        //
        // Apply H(i) to A(i+1:ihi,i+1:n) from the left
        //
        ApplyReflectionFromTheLeft(A, V, T, I+1, N, I+1, N, WORK);
        Inc(I);
    end;
end;


procedure UnpackQFromUpperHessenberg(const A : TReal2DArray;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     var Q : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
    IP1 : AlglibInteger;
    NMI : AlglibInteger;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if N=0 then
    begin
        Exit;
    end;
    
    //
    // init
    //
    SetLength(Q, N+1, N+1);
    SetLength(V, N+1);
    SetLength(WORK, N+1);
    I:=1;
    while I<=N do
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
    I:=1;
    while I<=N-1 do
    begin
        
        //
        // Apply H(i)
        //
        IP1 := I+1;
        NMI := N-I;
        i1_ := (IP1) - (1);
        for i_ := 1 to NMI do
        begin
            V[i_] := A[i_+i1_,I];
        end;
        V[1] := 1;
        ApplyReflectionFromTheRight(Q, Tau[I], V, 1, N, I+1, N, WORK);
        Inc(I);
    end;
end;


procedure UnpackHFromUpperHessenberg(const A : TReal2DArray;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     var H : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
begin
    if N=0 then
    begin
        Exit;
    end;
    SetLength(H, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        J:=1;
        while J<=I-2 do
        begin
            H[I,J] := 0;
            Inc(J);
        end;
        J := Max(1, I-1);
        APVMove(@H[I][0], J, N, @A[I][0], J, N);
        Inc(I);
    end;
end;


end.