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
unit cqr;
interface
uses Math, Sysutils, Ap, creflections;

procedure CMatrixQR(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TComplex1DArray);
procedure CMatrixQRUnpackQ(const QR : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TComplex1DArray;
     QColumns : AlglibInteger;
     var Q : TComplex2DArray);
procedure CMatrixQRUnpackR(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var R : TComplex2DArray);
procedure ComplexQRDecomposition(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TComplex1DArray);
procedure ComplexUnpackQFromQR(const QR : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TComplex1DArray;
     QColumns : AlglibInteger;
     var Q : TComplex2DArray);
procedure ComplexQRDecompositionUnpacked(A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Q : TComplex2DArray;
     var R : TComplex2DArray);

implementation

(*************************************************************************
QR decomposition of a rectangular complex matrix of size MxN

Input parameters:
    A   -   matrix A whose indexes range within [0..M-1, 0..N-1]
    M   -   number of rows in matrix A.
    N   -   number of columns in matrix A.

Output parameters:
    A   -   matrices Q and R in compact form
    Tau -   array of scalar factors which are used to form matrix Q. Array
            whose indexes range within [0.. Min(M,N)-1]

Matrix A is represented as A = QR, where Q is an orthogonal matrix of size
MxM, R - upper triangular (or upper trapezoid) matrix of size MxN.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************)
procedure CMatrixQR(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TComplex1DArray);
var
    WORK : TComplex1DArray;
    T : TComplex1DArray;
    I : AlglibInteger;
    K : AlglibInteger;
    MMI : AlglibInteger;
    MinMN : AlglibInteger;
    Tmp : Complex;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    MinMN := Min(M, N);
    if MinMN<=0 then
    begin
        Exit;
    end;
    SetLength(WORK, N-1+1);
    SetLength(T, M+1);
    SetLength(TAU, MinMN-1+1);
    
    //
    // Test the input arguments
    //
    K := Min(M, N);
    I:=0;
    while I<=K-1 do
    begin
        
        //
        // Generate elementary reflector H(i) to annihilate A(i+1:m,i)
        //
        MMI := M-I;
        i1_ := (I) - (1);
        for i_ := 1 to MMI do
        begin
            T[i_] := A[i_+i1_,I];
        end;
        ComplexGenerateReflection(T, MMI, Tmp);
        Tau[I] := Tmp;
        i1_ := (1) - (I);
        for i_ := I to M-1 do
        begin
            A[i_,I] := T[i_+i1_];
        end;
        T[1] := C_Complex(1);
        if I<N-1 then
        begin
            
            //
            // Apply H'(i) to A(i:m,i+1:n) from the left
            //
            ComplexApplyReflectionFromTheLeft(A, Conj(Tau[I]), T, I, M-1, I+1, N-1, WORK);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Partial unpacking of matrix Q from QR decomposition of a complex matrix A.

Input parameters:
    QR          -   matrices Q and R in compact form.
                    Output of CMatrixQR subroutine .
    M           -   number of rows in matrix A. M>=0.
    N           -   number of rows in matrix A. N>=0.
    Tau         -   scalar factors which are used to form Q.
                    Output of CMatrixQR subroutine .
    QColumns    -   required number of columns in matrix Q. M>=QColumns>=0.

Output parameters:
    Q           -   first QColumns columns of matrix Q.
                    Array whose index ranges within [0..M-1, 0..QColumns-1].
                    If QColumns=0, array isn't changed.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
procedure CMatrixQRUnpackQ(const QR : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TComplex1DArray;
     QColumns : AlglibInteger;
     var Q : TComplex2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    MinMN : AlglibInteger;
    V : TComplex1DArray;
    WORK : TComplex1DArray;
    VM : AlglibInteger;
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
                Q[I,J] := C_Complex(1);
            end
            else
            begin
                Q[I,J] := C_Complex(0);
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
        VM := M-I;
        i1_ := (I) - (1);
        for i_ := 1 to VM do
        begin
            V[i_] := QR[i_+i1_,I];
        end;
        V[1] := C_Complex(1);
        ComplexApplyReflectionFromTheLeft(Q, Tau[I], V, I, M-1, 0, QColumns-1, WORK);
        Dec(I);
    end;
end;


(*************************************************************************
Unpacking of matrix R from the QR decomposition of a matrix A

Input parameters:
    A       -   matrices Q and R in compact form.
                Output of CMatrixQR subroutine.
    M       -   number of rows in given matrix A. M>=0.
    N       -   number of columns in given matrix A. N>=0.

Output parameters:
    R       -   matrix R, array[0..M-1, 0..N-1].

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
procedure CMatrixQRUnpackR(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var R : TComplex2DArray);
var
    I : AlglibInteger;
    K : AlglibInteger;
    i_ : AlglibInteger;
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
        R[0,I] := C_Complex(0);
        Inc(I);
    end;
    I:=1;
    while I<=M-1 do
    begin
        for i_ := 0 to N-1 do
        begin
            R[I,i_] := R[0,i_];
        end;
        Inc(I);
    end;
    I:=0;
    while I<=K-1 do
    begin
        for i_ := I to N-1 do
        begin
            R[I,i_] := A[I,i_];
        end;
        Inc(I);
    end;
end;


procedure ComplexQRDecomposition(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Tau : TComplex1DArray);
var
    WORK : TComplex1DArray;
    T : TComplex1DArray;
    I : AlglibInteger;
    K : AlglibInteger;
    MMIP1 : AlglibInteger;
    MinMN : AlglibInteger;
    Tmp : Complex;
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
        ComplexGenerateReflection(T, MMIP1, Tmp);
        Tau[I] := Tmp;
        i1_ := (1) - (I);
        for i_ := I to M do
        begin
            A[i_,I] := T[i_+i1_];
        end;
        T[1] := C_Complex(1);
        if I<N then
        begin
            
            //
            // Apply H'(i) to A(i:m,i+1:n) from the left
            //
            ComplexApplyReflectionFromTheLeft(A, Conj(Tau[I]), T, I, M, I+1, N, WORK);
        end;
        Inc(I);
    end;
end;


procedure ComplexUnpackQFromQR(const QR : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const Tau : TComplex1DArray;
     QColumns : AlglibInteger;
     var Q : TComplex2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    MinMN : AlglibInteger;
    V : TComplex1DArray;
    WORK : TComplex1DArray;
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
                Q[I,J] := C_Complex(1);
            end
            else
            begin
                Q[I,J] := C_Complex(0);
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
            V[i_] := QR[i_+i1_,I];
        end;
        V[1] := C_Complex(1);
        ComplexApplyReflectionFromTheLeft(Q, Tau[I], V, I, M, 1, QColumns, WORK);
        Dec(I);
    end;
end;


procedure ComplexQRDecompositionUnpacked(A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Q : TComplex2DArray;
     var R : TComplex2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    L : AlglibInteger;
    VM : AlglibInteger;
    Tau : TComplex1DArray;
    WORK : TComplex1DArray;
    V : TComplex1DArray;
    Tmp : Double;
    i_ : AlglibInteger;
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
    ComplexQRDecomposition(A, M, N, Tau);
    
    //
    // R
    //
    I:=1;
    while I<=N do
    begin
        R[1,I] := C_Complex(0);
        Inc(I);
    end;
    I:=2;
    while I<=M do
    begin
        for i_ := 1 to N do
        begin
            R[I,i_] := R[1,i_];
        end;
        Inc(I);
    end;
    I:=1;
    while I<=K do
    begin
        for i_ := I to N do
        begin
            R[I,i_] := A[I,i_];
        end;
        Inc(I);
    end;
    
    //
    // Q
    //
    ComplexUnpackQFromQR(A, M, N, Tau, M, Q);
end;


end.