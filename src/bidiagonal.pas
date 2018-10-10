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
unit bidiagonal;
interface
uses Math, Sysutils, Ap, reflections;

procedure RMatrixBD(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var TauQ : TReal1DArray;
     var TauP : TReal1DArray);
procedure RMatrixBDUnpackQ(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauQ : TReal1DArray;
     QColumns : AlglibInteger;
     var Q : TReal2DArray);
procedure RMatrixBDMultiplyByQ(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauQ : TReal1DArray;
     var Z : TReal2DArray;
     ZRows : AlglibInteger;
     ZColumns : AlglibInteger;
     FromTheRight : Boolean;
     DoTranspose : Boolean);
procedure RMatrixBDUnpackPT(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauP : TReal1DArray;
     PTRows : AlglibInteger;
     var PT : TReal2DArray);
procedure RMatrixBDMultiplyByP(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauP : TReal1DArray;
     var Z : TReal2DArray;
     ZRows : AlglibInteger;
     ZColumns : AlglibInteger;
     FromTheRight : Boolean;
     DoTranspose : Boolean);
procedure RMatrixBDUnpackDiagonals(const B : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var IsUpper : Boolean;
     var D : TReal1DArray;
     var E : TReal1DArray);
procedure ToBidiagonal(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var TauQ : TReal1DArray;
     var TauP : TReal1DArray);
procedure UnpackQFromBidiagonal(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauQ : TReal1DArray;
     QColumns : AlglibInteger;
     var Q : TReal2DArray);
procedure MultiplyByQFromBidiagonal(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauQ : TReal1DArray;
     var Z : TReal2DArray;
     ZRows : AlglibInteger;
     ZColumns : AlglibInteger;
     FromTheRight : Boolean;
     DoTranspose : Boolean);
procedure UnpackPTFromBidiagonal(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauP : TReal1DArray;
     PTRows : AlglibInteger;
     var PT : TReal2DArray);
procedure MultiplyByPFromBidiagonal(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauP : TReal1DArray;
     var Z : TReal2DArray;
     ZRows : AlglibInteger;
     ZColumns : AlglibInteger;
     FromTheRight : Boolean;
     DoTranspose : Boolean);
procedure UnpackDiagonalsFromBidiagonal(const B : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var IsUpper : Boolean;
     var D : TReal1DArray;
     var E : TReal1DArray);

implementation

(*************************************************************************
Reduction of a rectangular matrix to  bidiagonal form

The algorithm reduces the rectangular matrix A to  bidiagonal form by
orthogonal transformations P and Q: A = Q*B*P.

Input parameters:
    A       -   source matrix. array[0..M-1, 0..N-1]
    M       -   number of rows in matrix A.
    N       -   number of columns in matrix A.

Output parameters:
    A       -   matrices Q, B, P in compact form (see below).
    TauQ    -   scalar factors which are used to form matrix Q.
    TauP    -   scalar factors which are used to form matrix P.

The main diagonal and one of the  secondary  diagonals  of  matrix  A  are
replaced with bidiagonal  matrix  B.  Other  elements  contain  elementary
reflections which form MxM matrix Q and NxN matrix P, respectively.

If M>=N, B is the upper  bidiagonal  MxN  matrix  and  is  stored  in  the
corresponding  elements  of  matrix  A.  Matrix  Q  is  represented  as  a
product   of   elementary   reflections   Q = H(0)*H(1)*...*H(n-1),  where
H(i) = 1-tau*v*v'. Here tau is a scalar which is stored  in  TauQ[i],  and
vector v has the following  structure:  v(0:i-1)=0, v(i)=1, v(i+1:m-1)  is
stored   in   elements   A(i+1:m-1,i).   Matrix   P  is  as  follows:  P =
G(0)*G(1)*...*G(n-2), where G(i) = 1 - tau*u*u'. Tau is stored in TauP[i],
u(0:i)=0, u(i+1)=1, u(i+2:n-1) is stored in elements A(i,i+2:n-1).

If M<N, B is the  lower  bidiagonal  MxN  matrix  and  is  stored  in  the
corresponding   elements  of  matrix  A.  Q = H(0)*H(1)*...*H(m-2),  where
H(i) = 1 - tau*v*v', tau is stored in TauQ, v(0:i)=0, v(i+1)=1, v(i+2:m-1)
is    stored    in   elements   A(i+2:m-1,i).    P = G(0)*G(1)*...*G(m-1),
G(i) = 1-tau*u*u', tau is stored in  TauP,  u(0:i-1)=0, u(i)=1, u(i+1:n-1)
is stored in A(i,i+1:n-1).

EXAMPLE:

m=6, n=5 (m > n):               m=5, n=6 (m < n):

(  d   e   u1  u1  u1 )         (  d   u1  u1  u1  u1  u1 )
(  v1  d   e   u2  u2 )         (  e   d   u2  u2  u2  u2 )
(  v1  v2  d   e   u3 )         (  v1  e   d   u3  u3  u3 )
(  v1  v2  v3  d   e  )         (  v1  v2  e   d   u4  u4 )
(  v1  v2  v3  v4  d  )         (  v1  v2  v3  e   d   u5 )
(  v1  v2  v3  v4  v5 )

Here vi and ui are vectors which form H(i) and G(i), and d and e -
are the diagonal and off-diagonal elements of matrix B.
*************************************************************************)
procedure RMatrixBD(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var TauQ : TReal1DArray;
     var TauP : TReal1DArray);
var
    WORK : TReal1DArray;
    T : TReal1DArray;
    MinMN : AlglibInteger;
    MaxMN : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    LTau : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // Prepare
    //
    if (N<=0) or (M<=0) then
    begin
        Exit;
    end;
    MinMN := Min(M, N);
    MaxMN := Max(M, N);
    SetLength(Work, MaxMN+1);
    SetLength(T, MaxMN+1);
    if M>=N then
    begin
        SetLength(TauQ, N-1+1);
        SetLength(TauP, N-1+1);
    end
    else
    begin
        SetLength(TauQ, M-1+1);
        SetLength(TauP, M-1+1);
    end;
    if M>=N then
    begin
        
        //
        // Reduce to upper bidiagonal form
        //
        I:=0;
        while I<=N-1 do
        begin
            
            //
            // Generate elementary reflector H(i) to annihilate A(i+1:m-1,i)
            //
            i1_ := (I) - (1);
            for i_ := 1 to M-I do
            begin
                T[i_] := A[i_+i1_,I];
            end;
            GenerateReflection(T, M-I, LTau);
            TauQ[I] := LTau;
            i1_ := (1) - (I);
            for i_ := I to M-1 do
            begin
                A[i_,I] := T[i_+i1_];
            end;
            T[1] := 1;
            
            //
            // Apply H(i) to A(i:m-1,i+1:n-1) from the left
            //
            ApplyReflectionFromTheLeft(A, LTau, T, I, M-1, I+1, N-1, WORK);
            if I<N-1 then
            begin
                
                //
                // Generate elementary reflector G(i) to annihilate
                // A(i,i+2:n-1)
                //
                APVMove(@T[0], 1, N-I-1, @A[I][0], I+1, N-1);
                GenerateReflection(T, N-1-I, LTau);
                TauP[I] := LTau;
                APVMove(@A[I][0], I+1, N-1, @T[0], 1, N-1-I);
                T[1] := 1;
                
                //
                // Apply G(i) to A(i+1:m-1,i+1:n-1) from the right
                //
                ApplyReflectionFromTheRight(A, LTau, T, I+1, M-1, I+1, N-1, WORK);
            end
            else
            begin
                TAUP[I] := 0;
            end;
            Inc(I);
        end;
    end
    else
    begin
        
        //
        // Reduce to lower bidiagonal form
        //
        I:=0;
        while I<=M-1 do
        begin
            
            //
            // Generate elementary reflector G(i) to annihilate A(i,i+1:n-1)
            //
            APVMove(@T[0], 1, N-I, @A[I][0], I, N-1);
            GenerateReflection(T, N-I, LTau);
            TauP[I] := LTau;
            APVMove(@A[I][0], I, N-1, @T[0], 1, N-I);
            T[1] := 1;
            
            //
            // Apply G(i) to A(i+1:m-1,i:n-1) from the right
            //
            ApplyReflectionFromTheRight(A, LTau, T, I+1, M-1, I, N-1, WORK);
            if I<M-1 then
            begin
                
                //
                // Generate elementary reflector H(i) to annihilate
                // A(i+2:m-1,i)
                //
                i1_ := (I+1) - (1);
                for i_ := 1 to M-1-I do
                begin
                    T[i_] := A[i_+i1_,I];
                end;
                GenerateReflection(T, M-1-I, LTau);
                TauQ[I] := LTau;
                i1_ := (1) - (I+1);
                for i_ := I+1 to M-1 do
                begin
                    A[i_,I] := T[i_+i1_];
                end;
                T[1] := 1;
                
                //
                // Apply H(i) to A(i+1:m-1,i+1:n-1) from the left
                //
                ApplyReflectionFromTheLeft(A, LTau, T, I+1, M-1, I+1, N-1, WORK);
            end
            else
            begin
                TAUQ[I] := 0;
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Unpacking matrix Q which reduces a matrix to bidiagonal form.

Input parameters:
    QP          -   matrices Q and P in compact form.
                    Output of ToBidiagonal subroutine.
    M           -   number of rows in matrix A.
    N           -   number of columns in matrix A.
    TAUQ        -   scalar factors which are used to form Q.
                    Output of ToBidiagonal subroutine.
    QColumns    -   required number of columns in matrix Q.
                    M>=QColumns>=0.

Output parameters:
    Q           -   first QColumns columns of matrix Q.
                    Array[0..M-1, 0..QColumns-1]
                    If QColumns=0, the array is not modified.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixBDUnpackQ(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauQ : TReal1DArray;
     QColumns : AlglibInteger;
     var Q : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Assert(QColumns<=M, 'RMatrixBDUnpackQ: QColumns>M!');
    Assert(QColumns>=0, 'RMatrixBDUnpackQ: QColumns<0!');
    if (M=0) or (N=0) or (QColumns=0) then
    begin
        Exit;
    end;
    
    //
    // prepare Q
    //
    SetLength(Q, M-1+1, QColumns-1+1);
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
    // Calculate
    //
    RMatrixBDMultiplyByQ(QP, M, N, TauQ, Q, M, QColumns, False, False);
end;


(*************************************************************************
Multiplication by matrix Q which reduces matrix A to  bidiagonal form.

The algorithm allows pre- or post-multiply by Q or Q'.

Input parameters:
    QP          -   matrices Q and P in compact form.
                    Output of ToBidiagonal subroutine.
    M           -   number of rows in matrix A.
    N           -   number of columns in matrix A.
    TAUQ        -   scalar factors which are used to form Q.
                    Output of ToBidiagonal subroutine.
    Z           -   multiplied matrix.
                    array[0..ZRows-1,0..ZColumns-1]
    ZRows       -   number of rows in matrix Z. If FromTheRight=False,
                    ZRows=M, otherwise ZRows can be arbitrary.
    ZColumns    -   number of columns in matrix Z. If FromTheRight=True,
                    ZColumns=M, otherwise ZColumns can be arbitrary.
    FromTheRight -  pre- or post-multiply.
    DoTranspose -   multiply by Q or Q'.

Output parameters:
    Z           -   product of Z and Q.
                    Array[0..ZRows-1,0..ZColumns-1]
                    If ZRows=0 or ZColumns=0, the array is not modified.

  -- ALGLIB --
     Copyright 2005 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixBDMultiplyByQ(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauQ : TReal1DArray;
     var Z : TReal2DArray;
     ZRows : AlglibInteger;
     ZColumns : AlglibInteger;
     FromTheRight : Boolean;
     DoTranspose : Boolean);
var
    I : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    IStep : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
    Mx : AlglibInteger;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if (M<=0) or (N<=0) or (ZRows<=0) or (ZColumns<=0) then
    begin
        Exit;
    end;
    Assert(FromTheRight and (ZColumns=M) or  not FromTheRight and (ZRows=M), 'RMatrixBDMultiplyByQ: incorrect Z size!');
    
    //
    // init
    //
    Mx := Max(M, N);
    Mx := Max(Mx, ZRows);
    Mx := Max(Mx, ZColumns);
    SetLength(V, Mx+1);
    SetLength(WORK, Mx+1);
    if M>=N then
    begin
        
        //
        // setup
        //
        if FromTheRight then
        begin
            I1 := 0;
            I2 := N-1;
            IStep := +1;
        end
        else
        begin
            I1 := N-1;
            I2 := 0;
            IStep := -1;
        end;
        if DoTranspose then
        begin
            I := I1;
            I1 := I2;
            I2 := I;
            IStep := -IStep;
        end;
        
        //
        // Process
        //
        I := I1;
        repeat
            i1_ := (I) - (1);
            for i_ := 1 to M-I do
            begin
                V[i_] := QP[i_+i1_,I];
            end;
            V[1] := 1;
            if FromTheRight then
            begin
                ApplyReflectionFromTheRight(Z, TAUQ[I], V, 0, ZRows-1, I, M-1, WORK);
            end
            else
            begin
                ApplyReflectionFromTheLeft(Z, TAUQ[I], V, I, M-1, 0, ZColumns-1, WORK);
            end;
            I := I+IStep;
        until I=I2+IStep;
    end
    else
    begin
        
        //
        // setup
        //
        if FromTheRight then
        begin
            I1 := 0;
            I2 := M-2;
            IStep := +1;
        end
        else
        begin
            I1 := M-2;
            I2 := 0;
            IStep := -1;
        end;
        if DoTranspose then
        begin
            I := I1;
            I1 := I2;
            I2 := I;
            IStep := -IStep;
        end;
        
        //
        // Process
        //
        if M-1>0 then
        begin
            I := I1;
            repeat
                i1_ := (I+1) - (1);
                for i_ := 1 to M-I-1 do
                begin
                    V[i_] := QP[i_+i1_,I];
                end;
                V[1] := 1;
                if FromTheRight then
                begin
                    ApplyReflectionFromTheRight(Z, TAUQ[I], V, 0, ZRows-1, I+1, M-1, WORK);
                end
                else
                begin
                    ApplyReflectionFromTheLeft(Z, TAUQ[I], V, I+1, M-1, 0, ZColumns-1, WORK);
                end;
                I := I+IStep;
            until I=I2+IStep;
        end;
    end;
end;


(*************************************************************************
Unpacking matrix P which reduces matrix A to bidiagonal form.
The subroutine returns transposed matrix P.

Input parameters:
    QP      -   matrices Q and P in compact form.
                Output of ToBidiagonal subroutine.
    M       -   number of rows in matrix A.
    N       -   number of columns in matrix A.
    TAUP    -   scalar factors which are used to form P.
                Output of ToBidiagonal subroutine.
    PTRows  -   required number of rows of matrix P^T. N >= PTRows >= 0.

Output parameters:
    PT      -   first PTRows columns of matrix P^T
                Array[0..PTRows-1, 0..N-1]
                If PTRows=0, the array is not modified.

  -- ALGLIB --
     Copyright 2005-2007 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixBDUnpackPT(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauP : TReal1DArray;
     PTRows : AlglibInteger;
     var PT : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Assert(PTRows<=N, 'RMatrixBDUnpackPT: PTRows>N!');
    Assert(PTRows>=0, 'RMatrixBDUnpackPT: PTRows<0!');
    if (M=0) or (N=0) or (PTRows=0) then
    begin
        Exit;
    end;
    
    //
    // prepare PT
    //
    SetLength(PT, PTRows-1+1, N-1+1);
    I:=0;
    while I<=PTRows-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if I=J then
            begin
                PT[I,J] := 1;
            end
            else
            begin
                PT[I,J] := 0;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Calculate
    //
    RMatrixBDMultiplyByP(QP, M, N, TauP, PT, PTRows, N, True, True);
end;


(*************************************************************************
Multiplication by matrix P which reduces matrix A to  bidiagonal form.

The algorithm allows pre- or post-multiply by P or P'.

Input parameters:
    QP          -   matrices Q and P in compact form.
                    Output of RMatrixBD subroutine.
    M           -   number of rows in matrix A.
    N           -   number of columns in matrix A.
    TAUP        -   scalar factors which are used to form P.
                    Output of RMatrixBD subroutine.
    Z           -   multiplied matrix.
                    Array whose indexes range within [0..ZRows-1,0..ZColumns-1].
    ZRows       -   number of rows in matrix Z. If FromTheRight=False,
                    ZRows=N, otherwise ZRows can be arbitrary.
    ZColumns    -   number of columns in matrix Z. If FromTheRight=True,
                    ZColumns=N, otherwise ZColumns can be arbitrary.
    FromTheRight -  pre- or post-multiply.
    DoTranspose -   multiply by P or P'.

Output parameters:
    Z - product of Z and P.
                Array whose indexes range within [0..ZRows-1,0..ZColumns-1].
                If ZRows=0 or ZColumns=0, the array is not modified.

  -- ALGLIB --
     Copyright 2005-2007 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixBDMultiplyByP(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauP : TReal1DArray;
     var Z : TReal2DArray;
     ZRows : AlglibInteger;
     ZColumns : AlglibInteger;
     FromTheRight : Boolean;
     DoTranspose : Boolean);
var
    I : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
    Mx : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    IStep : AlglibInteger;
begin
    if (M<=0) or (N<=0) or (ZRows<=0) or (ZColumns<=0) then
    begin
        Exit;
    end;
    Assert(FromTheRight and (ZColumns=N) or  not FromTheRight and (ZRows=N), 'RMatrixBDMultiplyByP: incorrect Z size!');
    
    //
    // init
    //
    Mx := Max(M, N);
    Mx := Max(Mx, ZRows);
    Mx := Max(Mx, ZColumns);
    SetLength(V, Mx+1);
    SetLength(WORK, Mx+1);
    SetLength(V, Mx+1);
    SetLength(WORK, Mx+1);
    if M>=N then
    begin
        
        //
        // setup
        //
        if FromTheRight then
        begin
            I1 := N-2;
            I2 := 0;
            IStep := -1;
        end
        else
        begin
            I1 := 0;
            I2 := N-2;
            IStep := +1;
        end;
        if  not DoTranspose then
        begin
            I := I1;
            I1 := I2;
            I2 := I;
            IStep := -IStep;
        end;
        
        //
        // Process
        //
        if N-1>0 then
        begin
            I := I1;
            repeat
                APVMove(@V[0], 1, N-1-I, @QP[I][0], I+1, N-1);
                V[1] := 1;
                if FromTheRight then
                begin
                    ApplyReflectionFromTheRight(Z, TAUP[I], V, 0, ZRows-1, I+1, N-1, WORK);
                end
                else
                begin
                    ApplyReflectionFromTheLeft(Z, TAUP[I], V, I+1, N-1, 0, ZColumns-1, WORK);
                end;
                I := I+IStep;
            until I=I2+IStep;
        end;
    end
    else
    begin
        
        //
        // setup
        //
        if FromTheRight then
        begin
            I1 := M-1;
            I2 := 0;
            IStep := -1;
        end
        else
        begin
            I1 := 0;
            I2 := M-1;
            IStep := +1;
        end;
        if  not DoTranspose then
        begin
            I := I1;
            I1 := I2;
            I2 := I;
            IStep := -IStep;
        end;
        
        //
        // Process
        //
        I := I1;
        repeat
            APVMove(@V[0], 1, N-I, @QP[I][0], I, N-1);
            V[1] := 1;
            if FromTheRight then
            begin
                ApplyReflectionFromTheRight(Z, TauP[I], V, 0, ZRows-1, I, N-1, WORK);
            end
            else
            begin
                ApplyReflectionFromTheLeft(Z, TauP[I], V, I, N-1, 0, ZColumns-1, WORK);
            end;
            I := I+IStep;
        until I=I2+IStep;
    end;
end;


(*************************************************************************
Unpacking of the main and secondary diagonals of bidiagonal decomposition
of matrix A.

Input parameters:
    B   -   output of RMatrixBD subroutine.
    M   -   number of rows in matrix B.
    N   -   number of columns in matrix B.

Output parameters:
    IsUpper -   True, if the matrix is upper bidiagonal.
                otherwise IsUpper is False.
    D       -   the main diagonal.
                Array whose index ranges within [0..Min(M,N)-1].
    E       -   the secondary diagonal (upper or lower, depending on
                the value of IsUpper).
                Array index ranges within [0..Min(M,N)-1], the last
                element is not used.

  -- ALGLIB --
     Copyright 2005-2007 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixBDUnpackDiagonals(const B : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var IsUpper : Boolean;
     var D : TReal1DArray;
     var E : TReal1DArray);
var
    I : AlglibInteger;
begin
    IsUpper := M>=N;
    if (M<=0) or (N<=0) then
    begin
        Exit;
    end;
    if IsUpper then
    begin
        SetLength(D, N-1+1);
        SetLength(E, N-1+1);
        I:=0;
        while I<=N-2 do
        begin
            D[I] := B[I,I];
            E[I] := B[I,I+1];
            Inc(I);
        end;
        D[N-1] := B[N-1,N-1];
    end
    else
    begin
        SetLength(D, M-1+1);
        SetLength(E, M-1+1);
        I:=0;
        while I<=M-2 do
        begin
            D[I] := B[I,I];
            E[I] := B[I+1,I];
            Inc(I);
        end;
        D[M-1] := B[M-1,M-1];
    end;
end;


procedure ToBidiagonal(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var TauQ : TReal1DArray;
     var TauP : TReal1DArray);
var
    WORK : TReal1DArray;
    T : TReal1DArray;
    MinMN : AlglibInteger;
    MaxMN : AlglibInteger;
    I : AlglibInteger;
    LTau : Double;
    MMIP1 : AlglibInteger;
    NMI : AlglibInteger;
    IP1 : AlglibInteger;
    NMIP1 : AlglibInteger;
    MMI : AlglibInteger;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    MinMN := Min(M, N);
    MaxMN := Max(M, N);
    SetLength(Work, MaxMN+1);
    SetLength(T, MaxMN+1);
    SetLength(TauP, MinMN+1);
    SetLength(TauQ, MinMN+1);
    if M>=N then
    begin
        
        //
        // Reduce to upper bidiagonal form
        //
        I:=1;
        while I<=N do
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
            GenerateReflection(T, MMIP1, LTau);
            TauQ[I] := LTau;
            i1_ := (1) - (I);
            for i_ := I to M do
            begin
                A[i_,I] := T[i_+i1_];
            end;
            T[1] := 1;
            
            //
            // Apply H(i) to A(i:m,i+1:n) from the left
            //
            ApplyReflectionFromTheLeft(A, LTau, T, I, M, I+1, N, WORK);
            if I<N then
            begin
                
                //
                // Generate elementary reflector G(i) to annihilate
                // A(i,i+2:n)
                //
                NMI := N-I;
                IP1 := I+1;
                APVMove(@T[0], 1, NMI, @A[I][0], IP1, N);
                GenerateReflection(T, NMI, LTau);
                TauP[I] := LTau;
                APVMove(@A[I][0], IP1, N, @T[0], 1, NMI);
                T[1] := 1;
                
                //
                // Apply G(i) to A(i+1:m,i+1:n) from the right
                //
                ApplyReflectionFromTheRight(A, LTau, T, I+1, M, I+1, N, WORK);
            end
            else
            begin
                TAUP[I] := 0;
            end;
            Inc(I);
        end;
    end
    else
    begin
        
        //
        // Reduce to lower bidiagonal form
        //
        I:=1;
        while I<=M do
        begin
            
            //
            // Generate elementary reflector G(i) to annihilate A(i,i+1:n)
            //
            NMIP1 := N-I+1;
            APVMove(@T[0], 1, NMIP1, @A[I][0], I, N);
            GenerateReflection(T, NMIP1, LTau);
            TauP[I] := LTau;
            APVMove(@A[I][0], I, N, @T[0], 1, NMIP1);
            T[1] := 1;
            
            //
            // Apply G(i) to A(i+1:m,i:n) from the right
            //
            ApplyReflectionFromTheRight(A, LTau, T, I+1, M, I, N, WORK);
            if I<M then
            begin
                
                //
                // Generate elementary reflector H(i) to annihilate
                // A(i+2:m,i)
                //
                MMI := M-I;
                IP1 := I+1;
                i1_ := (IP1) - (1);
                for i_ := 1 to MMI do
                begin
                    T[i_] := A[i_+i1_,I];
                end;
                GenerateReflection(T, MMI, LTau);
                TauQ[I] := LTau;
                i1_ := (1) - (IP1);
                for i_ := IP1 to M do
                begin
                    A[i_,I] := T[i_+i1_];
                end;
                T[1] := 1;
                
                //
                // Apply H(i) to A(i+1:m,i+1:n) from the left
                //
                ApplyReflectionFromTheLeft(A, LTau, T, I+1, M, I+1, N, WORK);
            end
            else
            begin
                TAUQ[I] := 0;
            end;
            Inc(I);
        end;
    end;
end;


procedure UnpackQFromBidiagonal(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauQ : TReal1DArray;
     QColumns : AlglibInteger;
     var Q : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    IP1 : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
    VM : AlglibInteger;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert(QColumns<=M, 'UnpackQFromBidiagonal: QColumns>M!');
    if (M=0) or (N=0) or (QColumns=0) then
    begin
        Exit;
    end;
    
    //
    // init
    //
    SetLength(Q, M+1, QColumns+1);
    SetLength(V, M+1);
    SetLength(WORK, QColumns+1);
    
    //
    // prepare Q
    //
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
    if M>=N then
    begin
        I:=Min(N, QColumns);
        while I>=1 do
        begin
            VM := M-I+1;
            i1_ := (I) - (1);
            for i_ := 1 to VM do
            begin
                V[i_] := QP[i_+i1_,I];
            end;
            V[1] := 1;
            ApplyReflectionFromTheLeft(Q, TAUQ[I], V, I, M, 1, QColumns, WORK);
            Dec(I);
        end;
    end
    else
    begin
        I:=Min(M-1, QColumns-1);
        while I>=1 do
        begin
            VM := M-I;
            IP1 := I+1;
            i1_ := (IP1) - (1);
            for i_ := 1 to VM do
            begin
                V[i_] := QP[i_+i1_,I];
            end;
            V[1] := 1;
            ApplyReflectionFromTheLeft(Q, TAUQ[I], V, I+1, M, 1, QColumns, WORK);
            Dec(I);
        end;
    end;
end;


procedure MultiplyByQFromBidiagonal(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauQ : TReal1DArray;
     var Z : TReal2DArray;
     ZRows : AlglibInteger;
     ZColumns : AlglibInteger;
     FromTheRight : Boolean;
     DoTranspose : Boolean);
var
    I : AlglibInteger;
    IP1 : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    IStep : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
    VM : AlglibInteger;
    Mx : AlglibInteger;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if (M<=0) or (N<=0) or (ZRows<=0) or (ZColumns<=0) then
    begin
        Exit;
    end;
    Assert(FromTheRight and (ZColumns=M) or  not FromTheRight and (ZRows=M), 'MultiplyByQFromBidiagonal: incorrect Z size!');
    
    //
    // init
    //
    Mx := Max(M, N);
    Mx := Max(Mx, ZRows);
    Mx := Max(Mx, ZColumns);
    SetLength(V, Mx+1);
    SetLength(WORK, Mx+1);
    if M>=N then
    begin
        
        //
        // setup
        //
        if FromTheRight then
        begin
            I1 := 1;
            I2 := N;
            IStep := +1;
        end
        else
        begin
            I1 := N;
            I2 := 1;
            IStep := -1;
        end;
        if DoTranspose then
        begin
            I := I1;
            I1 := I2;
            I2 := I;
            IStep := -IStep;
        end;
        
        //
        // Process
        //
        I := I1;
        repeat
            VM := M-I+1;
            i1_ := (I) - (1);
            for i_ := 1 to VM do
            begin
                V[i_] := QP[i_+i1_,I];
            end;
            V[1] := 1;
            if FromTheRight then
            begin
                ApplyReflectionFromTheRight(Z, TAUQ[I], V, 1, ZRows, I, M, WORK);
            end
            else
            begin
                ApplyReflectionFromTheLeft(Z, TAUQ[I], V, I, M, 1, ZColumns, WORK);
            end;
            I := I+IStep;
        until I=I2+IStep;
    end
    else
    begin
        
        //
        // setup
        //
        if FromTheRight then
        begin
            I1 := 1;
            I2 := M-1;
            IStep := +1;
        end
        else
        begin
            I1 := M-1;
            I2 := 1;
            IStep := -1;
        end;
        if DoTranspose then
        begin
            I := I1;
            I1 := I2;
            I2 := I;
            IStep := -IStep;
        end;
        
        //
        // Process
        //
        if M-1>0 then
        begin
            I := I1;
            repeat
                VM := M-I;
                IP1 := I+1;
                i1_ := (IP1) - (1);
                for i_ := 1 to VM do
                begin
                    V[i_] := QP[i_+i1_,I];
                end;
                V[1] := 1;
                if FromTheRight then
                begin
                    ApplyReflectionFromTheRight(Z, TAUQ[I], V, 1, ZRows, I+1, M, WORK);
                end
                else
                begin
                    ApplyReflectionFromTheLeft(Z, TAUQ[I], V, I+1, M, 1, ZColumns, WORK);
                end;
                I := I+IStep;
            until I=I2+IStep;
        end;
    end;
end;


procedure UnpackPTFromBidiagonal(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauP : TReal1DArray;
     PTRows : AlglibInteger;
     var PT : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    IP1 : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
    VM : AlglibInteger;
begin
    Assert(PTRows<=N, 'UnpackPTFromBidiagonal: PTRows>N!');
    if (M=0) or (N=0) or (PTRows=0) then
    begin
        Exit;
    end;
    
    //
    // init
    //
    SetLength(PT, PTRows+1, N+1);
    SetLength(V, N+1);
    SetLength(WORK, PTRows+1);
    
    //
    // prepare PT
    //
    I:=1;
    while I<=PTRows do
    begin
        J:=1;
        while J<=N do
        begin
            if I=J then
            begin
                PT[I,J] := 1;
            end
            else
            begin
                PT[I,J] := 0;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    if M>=N then
    begin
        I:=Min(N-1, PTRows-1);
        while I>=1 do
        begin
            VM := N-I;
            IP1 := I+1;
            APVMove(@V[0], 1, VM, @QP[I][0], IP1, N);
            V[1] := 1;
            ApplyReflectionFromTheRight(PT, TAUP[I], V, 1, PTRows, I+1, N, WORK);
            Dec(I);
        end;
    end
    else
    begin
        I:=Min(M, PTRows);
        while I>=1 do
        begin
            VM := N-I+1;
            APVMove(@V[0], 1, VM, @QP[I][0], I, N);
            V[1] := 1;
            ApplyReflectionFromTheRight(PT, TAUP[I], V, 1, PTRows, I, N, WORK);
            Dec(I);
        end;
    end;
end;


procedure MultiplyByPFromBidiagonal(const QP : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const TauP : TReal1DArray;
     var Z : TReal2DArray;
     ZRows : AlglibInteger;
     ZColumns : AlglibInteger;
     FromTheRight : Boolean;
     DoTranspose : Boolean);
var
    I : AlglibInteger;
    IP1 : AlglibInteger;
    V : TReal1DArray;
    WORK : TReal1DArray;
    VM : AlglibInteger;
    Mx : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    IStep : AlglibInteger;
begin
    if (M<=0) or (N<=0) or (ZRows<=0) or (ZColumns<=0) then
    begin
        Exit;
    end;
    Assert(FromTheRight and (ZColumns=N) or  not FromTheRight and (ZRows=N), 'MultiplyByQFromBidiagonal: incorrect Z size!');
    
    //
    // init
    //
    Mx := Max(M, N);
    Mx := Max(Mx, ZRows);
    Mx := Max(Mx, ZColumns);
    SetLength(V, Mx+1);
    SetLength(WORK, Mx+1);
    SetLength(V, Mx+1);
    SetLength(WORK, Mx+1);
    if M>=N then
    begin
        
        //
        // setup
        //
        if FromTheRight then
        begin
            I1 := N-1;
            I2 := 1;
            IStep := -1;
        end
        else
        begin
            I1 := 1;
            I2 := N-1;
            IStep := +1;
        end;
        if  not DoTranspose then
        begin
            I := I1;
            I1 := I2;
            I2 := I;
            IStep := -IStep;
        end;
        
        //
        // Process
        //
        if N-1>0 then
        begin
            I := I1;
            repeat
                VM := N-I;
                IP1 := I+1;
                APVMove(@V[0], 1, VM, @QP[I][0], IP1, N);
                V[1] := 1;
                if FromTheRight then
                begin
                    ApplyReflectionFromTheRight(Z, TAUP[I], V, 1, ZRows, I+1, N, WORK);
                end
                else
                begin
                    ApplyReflectionFromTheLeft(Z, TAUP[I], V, I+1, N, 1, ZColumns, WORK);
                end;
                I := I+IStep;
            until I=I2+IStep;
        end;
    end
    else
    begin
        
        //
        // setup
        //
        if FromTheRight then
        begin
            I1 := M;
            I2 := 1;
            IStep := -1;
        end
        else
        begin
            I1 := 1;
            I2 := M;
            IStep := +1;
        end;
        if  not DoTranspose then
        begin
            I := I1;
            I1 := I2;
            I2 := I;
            IStep := -IStep;
        end;
        
        //
        // Process
        //
        I := I1;
        repeat
            VM := N-I+1;
            APVMove(@V[0], 1, VM, @QP[I][0], I, N);
            V[1] := 1;
            if FromTheRight then
            begin
                ApplyReflectionFromTheRight(Z, TauP[I], V, 1, ZRows, I, N, WORK);
            end
            else
            begin
                ApplyReflectionFromTheLeft(Z, TauP[I], V, I, N, 1, ZColumns, WORK);
            end;
            I := I+IStep;
        until I=I2+IStep;
    end;
end;


procedure UnpackDiagonalsFromBidiagonal(const B : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var IsUpper : Boolean;
     var D : TReal1DArray;
     var E : TReal1DArray);
var
    I : AlglibInteger;
begin
    IsUpper := M>=N;
    if (M=0) or (N=0) then
    begin
        Exit;
    end;
    if IsUpper then
    begin
        SetLength(D, N+1);
        SetLength(E, N+1);
        I:=1;
        while I<=N-1 do
        begin
            D[I] := B[I,I];
            E[I] := B[I,I+1];
            Inc(I);
        end;
        D[N] := B[N,N];
    end
    else
    begin
        SetLength(D, M+1);
        SetLength(E, M+1);
        I:=1;
        while I<=M-1 do
        begin
            D[I] := B[I,I];
            E[I] := B[I+1,I];
            Inc(I);
        end;
        D[M] := B[M,M];
    end;
end;


end.