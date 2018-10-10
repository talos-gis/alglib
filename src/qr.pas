(*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

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
unit qr;
interface
uses Math, Sysutils, Ap, reflections;

procedure RMatrixQR(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TReal1DArray);
procedure RMatrixQRUnpackQ(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     QColumns : AlglibInteger;
     var Q : TReal2DArray);
procedure RMatrixQRUnpackR(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var R : TReal2DArray);
procedure QRDecomposition(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TReal1DArray);
procedure UnpackQFromQR(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     QColumns : AlglibInteger;
     var Q : TReal2DArray);
procedure QRDecompositionUnpacked(A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Q : TReal2DArray;
     var R : TReal2DArray);

implementation

(*************************************************************************
QR decomposition of a rectangular matrix of size MxN

Input parameters:
    A   -   matrix A whose indexes range within [0..M-1, 0..N-1].
    M   -   number of rows in matrix A.
    N   -   number of columns in matrix A.

Output parameters:
    A   -   matrices Q and R in compact form (see below).
    Tau -   array of scalar factors which are used to form
            matrix Q. Array whose index ranges within [0.. Min(M-1,N-1)].

Matrix A is represented as A = QR, where Q is an orthogonal matrix of size
MxM, R - upper triangular (or upper trapezoid) matrix of size M x N.

The elements of matrix R are located on and above the main diagonal of
matrix A. The elements which are located in Tau array and below the main
diagonal of matrix A are used to form matrix Q as follows:

Matrix Q is represented as a product of elementary reflections

Q = H(0)*H(2)*...*H(k-1),

where k = min(m,n), and each H(i) is in the form

H(i) = 1 - tau * v * (v^T)

where tau is a scalar stored in Tau[I]; v - real vector,
so that v(0:i-1) = 0, v(i) = 1, v(i+1:m-1) stored in A(i+1:m-1,i).

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992.
     Translation from FORTRAN to pseudocode (AlgoPascal)
     by Sergey Bochkanov, ALGLIB project, 2005-2007.
*************************************************************************)
procedure RMatrixQR(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TReal1DArray);
var
    WORK : TReal1DArray;
    T : TReal1DArray;
    I : AlglibInteger;
    K : AlglibInteger;
    MinMN : AlglibInteger;
    Tmp : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if (M<=0) or (N<=0) then
    begin
        Exit;
    end;
    MinMN := Min(M, N);
    SetLength(WORK, N-1+1);
    SetLength(T, M+1);
    SetLength(TAU, MinMN-1+1);
    
    //
    // Test the input arguments
    //
    K := MinMN;
    I:=0;
    while I<=K-1 do
    begin
        
        //
        // Generate elementary reflector H(i) to annihilate A(i+1:m,i)
        //
        i1_ := (I) - (1);
        for i_ := 1 to M-I do
        begin
            T[i_] := A[i_+i1_,I];
        end;
        GenerateReflection(T, M-I, Tmp);
        Tau[I] := Tmp;
        i1_ := (1) - (I);
        for i_ := I to M-1 do
        begin
            A[i_,I] := T[i_+i1_];
        end;
        T[1] := 1;
        if I<N then
        begin
            
            //
            // Apply H(i) to A(i:m-1,i+1:n-1) from the left
            //
            ApplyReflectionFromTheLeft(A, Tau[I], T, I, M-1, I+1, N-1, WORK);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Partial unpacking of matrix Q from the QR decomposition of a matrix A

Input parameters:
    A       -   matrices Q and R in compact form.
                Output of RMatrixQR subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.
    Tau     -   scalar factors which are used to form Q.
                Output of the RMatrixQR subroutine.
    QColumns -  required number of columns of matrix Q. M>=QColumns>=0.

Output parameters:
    Q       -   first QColumns columns of matrix Q.
                Array whose indexes range within [0..M-1, 0..QColumns-1].
                If QColumns=0, the array remains unchanged.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixQRUnpackQ(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     QColumns : AlglibInteger;
     var Q : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    MinMN : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert(QColumns<=M, 'UnpackQFromQR: QColumns>M!');
    if (M<=0) or (N<=0) or (QColumns<=0) then
    begin
        Exit;
    end;
    
    //
    // init
    //
    MinMN := Min(M, N);
    K := Min(MinMN, QColumns);
    SetLength(Q, M-1+1, QColumns-1+1);
    SetLength(V, M+1);
    SetLength(WORK, QColumns-1+1);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=QColumns-1 do
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
        i1_ := (I) - (1);
        for i_ := 1 to M-I do
        begin
            V[i_] := A[i_+i1_,I];
        end;
        V[1] := 1;
        ApplyReflectionFromTheLeft(Q, Tau[I], V, I, M-1, 0, QColumns-1, WORK);
        Dec(I);
    end;
end;


(*************************************************************************
Unpacking of matrix R from the QR decomposition of a matrix A

Input parameters:
    A       -   matrices Q and R in compact form.
                Output of RMatrixQR subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.

Output parameters:
    R       -   matrix R, array[0..M-1, 0..N-1].

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixQRUnpackR(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var R : TReal2DArray);
var
    I : AlglibInteger;
    K : AlglibInteger;
begin
    if (M<=0) or (N<=0) then
    begin
        Exit;
    end;
    K := Min(M, N);
    SetLength(R, M-1+1, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        R[0,I] := 0;
        Inc(I);
    end;
    I:=1;
    while I<=M-1 do
    begin
        APVMove(@R[I][0], 0, N-1, @R[0][0], 0, N-1);
        Inc(I);
    end;
    I:=0;
    while I<=K-1 do
    begin
        APVMove(@R[I][0], I, N-1, @A[I][0], I, N-1);
        Inc(I);
    end;
end;


procedure QRDecomposition(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TReal1DArray);
var
    WORK : TReal1DArray;
    T : TReal1DArray;
    I : AlglibInteger;
    K : AlglibInteger;
    MMIP1 : AlglibInteger;
    MinMN : AlglibInteger;
    Tmp : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    MinMN := Min(M, N);
    SetLength(WORK, N+1);
    SetLength(T, M+1);
    SetLength(TAU, MinMN+1);
    
    //
    // Test the input arguments
    //
    K := Min(M, N);
    I:=1;
    while I<=K do
    begin
        
        //
        // Generate elementary reflector H(i) to annihilate A(i+1:m,i)
        //
        MMIP1 := M-I+1;
        i1_ := (I) - (1);
        for i_ := 1 to MMIP1 do
        begin
            T[i_] := A[i_+i1_,I];
        end;
        GenerateReflection(T, MMIP1, Tmp);
        Tau[I] := Tmp;
        i1_ := (1) - (I);
        for i_ := I to M do
        begin
            A[i_,I] := T[i_+i1_];
        end;
        T[1] := 1;
        if I<N then
        begin
            
            //
            // Apply H(i) to A(i:m,i+1:n) from the left
            //
            ApplyReflectionFromTheLeft(A, Tau[I], T, I, M, I+1, N, WORK);
        end;
        Inc(I);
    end;
end;


procedure UnpackQFromQR(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TReal1DArray;
     QColumns : AlglibInteger;
     var Q : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    MinMN : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
    VM : AlglibInteger;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert(QColumns<=M, 'UnpackQFromQR: QColumns>M!');
    if (M=0) or (N=0) or (QColumns=0) then
    begin
        Exit;
    end;
    
    //
    // init
    //
    MinMN := Min(M, N);
    K := Min(MinMN, QColumns);
    SetLength(Q, M+1, QColumns+1);
    SetLength(V, M+1);
    SetLength(WORK, QColumns+1);
    I:=1;
    while I<=M do
    begin
        J:=1;
        while J<=QColumns do
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
        VM := M-I+1;
        i1_ := (I) - (1);
        for i_ := 1 to VM do
        begin
            V[i_] := A[i_+i1_,I];
        end;
        V[1] := 1;
        ApplyReflectionFromTheLeft(Q, Tau[I], V, I, M, 1, QColumns, WORK);
        Dec(I);
    end;
end;


procedure QRDecompositionUnpacked(A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Q : TReal2DArray;
     var R : TReal2DArray);
var
    I : AlglibInteger;
    K : AlglibInteger;
    Tau : TReal1DArray;
    WORK : TReal1DArray;
    V : TReal1DArray;
begin
    A := DynamicArrayCopy(A);
    K := Min(M, N);
    if N<=0 then
    begin
        Exit;
    end;
    SetLength(WORK, M+1);
    SetLength(V, M+1);
    SetLength(Q, M+1, M+1);
    SetLength(R, M+1, N+1);
    
    //
    // QRDecomposition
    //
    QRDecomposition(A, M, N, Tau);
    
    //
    // R
    //
    I:=1;
    while I<=N do
    begin
        R[1,I] := 0;
        Inc(I);
    end;
    I:=2;
    while I<=M do
    begin
        APVMove(@R[I][0], 1, N, @R[1][0], 1, N);
        Inc(I);
    end;
    I:=1;
    while I<=K do
    begin
        APVMove(@R[I][0], I, N, @A[I][0], I, N);
        Inc(I);
    end;
    
    //
    // Q
    //
    UnpackQFromQR(A, M, N, Tau, M, Q);
end;


end.