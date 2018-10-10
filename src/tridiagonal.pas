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
unit tridiagonal;
interface
uses Math, Sysutils, Ap, sblas, reflections;

procedure SMatrixTD(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tau : TReal1DArray;
     var D : TReal1DArray;
     var E : TReal1DArray);
procedure SMatrixTDUnpackQ(const A : TReal2DArray;
     const N : AlglibInteger;
     const IsUpper : Boolean;
     const Tau : TReal1DArray;
     var Q : TReal2DArray);
procedure ToTridiagonal(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tau : TReal1DArray;
     var D : TReal1DArray;
     var E : TReal1DArray);
procedure UnpackQFromTridiagonal(const A : TReal2DArray;
     const N : AlglibInteger;
     const IsUpper : Boolean;
     const Tau : TReal1DArray;
     var Q : TReal2DArray);

implementation

(*************************************************************************
Reduction of a symmetric matrix which is given by its higher or lower
triangular part to a tridiagonal matrix using orthogonal similarity
transformation: Q'*A*Q=T.

Input parameters:
    A       -   matrix to be transformed
                array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format. If IsUpper = True, then matrix A is given
                by its upper triangle, and the lower triangle is not used
                and not modified by the algorithm, and vice versa
                if IsUpper = False.

Output parameters:
    A       -   matrices T and Q in  compact form (see lower)
    Tau     -   array of factors which are forming matrices H(i)
                array with elements [0..N-2].
    D       -   main diagonal of symmetric matrix T.
                array with elements [0..N-1].
    E       -   secondary diagonal of symmetric matrix T.
                array with elements [0..N-2].


  If IsUpper=True, the matrix Q is represented as a product of elementary
  reflectors

     Q = H(n-2) . . . H(2) H(0).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a real scalar, and v is a real vector with
  v(i+1:n-1) = 0, v(i) = 1, v(0:i-1) is stored on exit in
  A(0:i-1,i+1), and tau in TAU(i).

  If IsUpper=False, the matrix Q is represented as a product of elementary
  reflectors

     Q = H(0) H(2) . . . H(n-2).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a real scalar, and v is a real vector with
  v(0:i) = 0, v(i+1) = 1, v(i+2:n-1) is stored on exit in A(i+2:n-1,i),
  and tau in TAU(i).

  The contents of A on exit are illustrated by the following examples
  with n = 5:

  if UPLO = 'U':                       if UPLO = 'L':

    (  d   e   v1  v2  v3 )              (  d                  )
    (      d   e   v2  v3 )              (  e   d              )
    (          d   e   v3 )              (  v0  e   d          )
    (              d   e  )              (  v0  v1  e   d      )
    (                  d  )              (  v0  v1  v2  e   d  )

  where d and e denote diagonal and off-diagonal elements of T, and vi
  denotes an element of the vector defining H(i).

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1992
*************************************************************************)
procedure SMatrixTD(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tau : TReal1DArray;
     var D : TReal1DArray;
     var E : TReal1DArray);
var
    I : AlglibInteger;
    ALPHA : Double;
    TAUI : Double;
    V : Double;
    T : TReal1DArray;
    T2 : TReal1DArray;
    T3 : TReal1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Exit;
    end;
    SetLength(T, N+1);
    SetLength(T2, N+1);
    SetLength(T3, N+1);
    if N>1 then
    begin
        SetLength(Tau, N-2+1);
    end;
    SetLength(D, N-1+1);
    if N>1 then
    begin
        SetLength(E, N-2+1);
    end;
    if IsUpper then
    begin
        
        //
        // Reduce the upper triangle of A
        //
        I:=N-2;
        while I>=0 do
        begin
            
            //
            // Generate elementary reflector H() = E - tau * v * v'
            //
            if I>=1 then
            begin
                i1_ := (0) - (2);
                for i_ := 2 to I+1 do
                begin
                    T[i_] := A[i_+i1_,I+1];
                end;
            end;
            T[1] := A[I,I+1];
            GenerateReflection(T, I+1, TauI);
            if I>=1 then
            begin
                i1_ := (2) - (0);
                for i_ := 0 to I-1 do
                begin
                    A[i_,I+1] := T[i_+i1_];
                end;
            end;
            A[I,I+1] := T[1];
            E[I] := A[I,I+1];
            if AP_FP_Neq(TAUI,0) then
            begin
                
                //
                // Apply H from both sides to A
                //
                A[I,I+1] := 1;
                
                //
                // Compute  x := tau * A * v  storing x in TAU
                //
                i1_ := (0) - (1);
                for i_ := 1 to I+1 do
                begin
                    T[i_] := A[i_+i1_,I+1];
                end;
                SymmetricMatrixVectorMultiply(A, IsUpper, 0, I, T, TauI, T3);
                APVMove(@Tau[0], 0, I, @T3[0], 1, I+1);
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                V := 0.0;
                for i_ := 0 to I do
                begin
                    V := V + Tau[i_]*A[i_,I+1];
                end;
                ALPHA := -0.5*TAUI*V;
                for i_ := 0 to I do
                begin
                    Tau[i_] := Tau[i_] + Alpha*A[i_,I+1];
                end;
                
                //
                // Apply the transformation as a rank-2 update:
                //    A := A - v * w' - w * v'
                //
                i1_ := (0) - (1);
                for i_ := 1 to I+1 do
                begin
                    T[i_] := A[i_+i1_,I+1];
                end;
                APVMove(@T3[0], 1, I+1, @Tau[0], 0, I);
                SymmetricRank2Update(A, IsUpper, 0, I, T, T3, T2, -1);
                A[I,I+1] := E[I];
            end;
            D[I+1] := A[I+1,I+1];
            TAU[I] := TAUI;
            Dec(I);
        end;
        D[0] := A[0,0];
    end
    else
    begin
        
        //
        // Reduce the lower triangle of A
        //
        I:=0;
        while I<=N-2 do
        begin
            
            //
            // Generate elementary reflector H = E - tau * v * v'
            //
            i1_ := (I+1) - (1);
            for i_ := 1 to N-I-1 do
            begin
                T[i_] := A[i_+i1_,I];
            end;
            GenerateReflection(T, N-I-1, TauI);
            i1_ := (1) - (I+1);
            for i_ := I+1 to N-1 do
            begin
                A[i_,I] := T[i_+i1_];
            end;
            E[I] := A[I+1,I];
            if AP_FP_Neq(TAUI,0) then
            begin
                
                //
                // Apply H from both sides to A
                //
                A[I+1,I] := 1;
                
                //
                // Compute  x := tau * A * v  storing y in TAU
                //
                i1_ := (I+1) - (1);
                for i_ := 1 to N-I-1 do
                begin
                    T[i_] := A[i_+i1_,I];
                end;
                SymmetricMatrixVectorMultiply(A, IsUpper, I+1, N-1, T, TauI, T2);
                APVMove(@Tau[0], I, N-2, @T2[0], 1, N-I-1);
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                i1_ := (I+1)-(I);
                V := 0.0;
                for i_ := I to N-2 do
                begin
                    V := V + Tau[i_]*A[i_+i1_,I];
                end;
                ALPHA := -0.5*TAUI*V;
                i1_ := (I+1) - (I);
                for i_ := I to N-2 do
                begin
                    Tau[i_] := Tau[i_] + Alpha*A[i_+i1_,I];
                end;
                
                //
                // Apply the transformation as a rank-2 update:
                //     A := A - v * w' - w * v'
                //
                //
                i1_ := (I+1) - (1);
                for i_ := 1 to N-I-1 do
                begin
                    T[i_] := A[i_+i1_,I];
                end;
                APVMove(@T2[0], 1, N-I-1, @Tau[0], I, N-2);
                SymmetricRank2Update(A, IsUpper, I+1, N-1, T, T2, T3, -1);
                A[I+1,I] := E[I];
            end;
            D[I] := A[I,I];
            TAU[I] := TAUI;
            Inc(I);
        end;
        D[N-1] := A[N-1,N-1];
    end;
end;


(*************************************************************************
Unpacking matrix Q which reduces symmetric matrix to a tridiagonal
form.

Input parameters:
    A       -   the result of a SMatrixTD subroutine
    N       -   size of matrix A.
    IsUpper -   storage format (a parameter of SMatrixTD subroutine)
    Tau     -   the result of a SMatrixTD subroutine

Output parameters:
    Q       -   transformation matrix.
                array with elements [0..N-1, 0..N-1].

  -- ALGLIB --
     Copyright 2005-2008 by Bochkanov Sergey
*************************************************************************)
procedure SMatrixTDUnpackQ(const A : TReal2DArray;
     const N : AlglibInteger;
     const IsUpper : Boolean;
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
    SetLength(V, N+1);
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
    if IsUpper then
    begin
        I:=0;
        while I<=N-2 do
        begin
            
            //
            // Apply H(i)
            //
            i1_ := (0) - (1);
            for i_ := 1 to I+1 do
            begin
                V[i_] := A[i_+i1_,I+1];
            end;
            V[I+1] := 1;
            ApplyReflectionFromTheLeft(Q, Tau[I], V, 0, I, 0, N-1, WORK);
            Inc(I);
        end;
    end
    else
    begin
        I:=N-2;
        while I>=0 do
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
            ApplyReflectionFromTheLeft(Q, Tau[I], V, I+1, N-1, 0, N-1, WORK);
            Dec(I);
        end;
    end;
end;


procedure ToTridiagonal(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tau : TReal1DArray;
     var D : TReal1DArray;
     var E : TReal1DArray);
var
    I : AlglibInteger;
    IP1 : AlglibInteger;
    IM1 : AlglibInteger;
    NMI : AlglibInteger;
    NM1 : AlglibInteger;
    ALPHA : Double;
    TAUI : Double;
    V : Double;
    T : TReal1DArray;
    T2 : TReal1DArray;
    T3 : TReal1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Exit;
    end;
    SetLength(T, N+1);
    SetLength(T2, N+1);
    SetLength(T3, N+1);
    SetLength(Tau, Max(1, N-1)+1);
    SetLength(D, N+1);
    SetLength(E, Max(1, N-1)+1);
    if IsUpper then
    begin
        
        //
        // Reduce the upper triangle of A
        //
        I:=N-1;
        while I>=1 do
        begin
            
            //
            // Generate elementary reflector H(i) = I - tau * v * v'
            // to annihilate A(1:i-1,i+1)
            //
            // DLARFG( I, A( I, I+1 ), A( 1, I+1 ), 1, TAUI );
            //
            IP1 := I+1;
            IM1 := I-1;
            if I>=2 then
            begin
                i1_ := (1) - (2);
                for i_ := 2 to I do
                begin
                    T[i_] := A[i_+i1_,IP1];
                end;
            end;
            T[1] := A[I,IP1];
            GenerateReflection(T, I, TauI);
            if I>=2 then
            begin
                i1_ := (2) - (1);
                for i_ := 1 to IM1 do
                begin
                    A[i_,IP1] := T[i_+i1_];
                end;
            end;
            A[I,IP1] := T[1];
            E[I] := A[I,I+1];
            if AP_FP_Neq(TAUI,0) then
            begin
                
                //
                // Apply H(i) from both sides to A(1:i,1:i)
                //
                A[I,I+1] := 1;
                
                //
                // Compute  x := tau * A * v  storing x in TAU(1:i)
                //
                // DSYMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO, TAU, 1 );
                //
                IP1 := I+1;
                for i_ := 1 to I do
                begin
                    T[i_] := A[i_,IP1];
                end;
                SymmetricMatrixVectorMultiply(A, IsUpper, 1, I, T, TauI, Tau);
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                IP1 := I+1;
                V := 0.0;
                for i_ := 1 to I do
                begin
                    V := V + Tau[i_]*A[i_,IP1];
                end;
                ALPHA := -0.5*TAUI*V;
                for i_ := 1 to I do
                begin
                    Tau[i_] := Tau[i_] + Alpha*A[i_,IP1];
                end;
                
                //
                // Apply the transformation as a rank-2 update:
                //    A := A - v * w' - w * v'
                //
                // DSYR2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A, LDA );
                //
                for i_ := 1 to I do
                begin
                    T[i_] := A[i_,IP1];
                end;
                SymmetricRank2Update(A, IsUpper, 1, I, T, Tau, T2, -1);
                A[I,I+1] := E[I];
            end;
            D[I+1] := A[I+1,I+1];
            TAU[I] := TAUI;
            Dec(I);
        end;
        D[1] := A[1,1];
    end
    else
    begin
        
        //
        // Reduce the lower triangle of A
        //
        I:=1;
        while I<=N-1 do
        begin
            
            //
            // Generate elementary reflector H(i) = I - tau * v * v'
            // to annihilate A(i+2:n,i)
            //
            //DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), 1, TAUI );
            //
            NMI := N-I;
            IP1 := I+1;
            i1_ := (IP1) - (1);
            for i_ := 1 to NMI do
            begin
                T[i_] := A[i_+i1_,I];
            end;
            GenerateReflection(T, NMI, TauI);
            i1_ := (1) - (IP1);
            for i_ := IP1 to N do
            begin
                A[i_,I] := T[i_+i1_];
            end;
            E[I] := A[I+1,I];
            if AP_FP_Neq(TAUI,0) then
            begin
                
                //
                // Apply H(i) from both sides to A(i+1:n,i+1:n)
                //
                A[I+1,I] := 1;
                
                //
                // Compute  x := tau * A * v  storing y in TAU(i:n-1)
                //
                //DSYMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO, TAU( I ), 1 );
                //
                IP1 := I+1;
                NMI := N-I;
                NM1 := N-1;
                i1_ := (IP1) - (1);
                for i_ := 1 to NMI do
                begin
                    T[i_] := A[i_+i1_,I];
                end;
                SymmetricMatrixVectorMultiply(A, IsUpper, I+1, N, T, TauI, T2);
                APVMove(@Tau[0], I, NM1, @T2[0], 1, NMI);
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                NM1 := N-1;
                IP1 := I+1;
                i1_ := (IP1)-(I);
                V := 0.0;
                for i_ := I to NM1 do
                begin
                    V := V + Tau[i_]*A[i_+i1_,I];
                end;
                ALPHA := -0.5*TAUI*V;
                i1_ := (IP1) - (I);
                for i_ := I to NM1 do
                begin
                    Tau[i_] := Tau[i_] + Alpha*A[i_+i1_,I];
                end;
                
                //
                // Apply the transformation as a rank-2 update:
                //     A := A - v * w' - w * v'
                //
                //DSYR2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1, A( I+1, I+1 ), LDA );
                //
                NM1 := N-1;
                NMI := N-I;
                IP1 := I+1;
                i1_ := (IP1) - (1);
                for i_ := 1 to NMI do
                begin
                    T[i_] := A[i_+i1_,I];
                end;
                APVMove(@T2[0], 1, NMI, @Tau[0], I, NM1);
                SymmetricRank2Update(A, IsUpper, I+1, N, T, T2, T3, -1);
                A[I+1,I] := E[I];
            end;
            D[I] := A[I,I];
            TAU[I] := TAUI;
            Inc(I);
        end;
        D[N] := A[N,N];
    end;
end;


procedure UnpackQFromTridiagonal(const A : TReal2DArray;
     const N : AlglibInteger;
     const IsUpper : Boolean;
     const Tau : TReal1DArray;
     var Q : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    IP1 : AlglibInteger;
    NMI : AlglibInteger;
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
    if IsUpper then
    begin
        I:=1;
        while I<=N-1 do
        begin
            
            //
            // Apply H(i)
            //
            IP1 := I+1;
            for i_ := 1 to I do
            begin
                V[i_] := A[i_,IP1];
            end;
            V[I] := 1;
            ApplyReflectionFromTheLeft(Q, Tau[I], V, 1, I, 1, N, WORK);
            Inc(I);
        end;
    end
    else
    begin
        I:=N-1;
        while I>=1 do
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
            ApplyReflectionFromTheLeft(Q, Tau[I], V, I+1, N, 1, N, WORK);
            Dec(I);
        end;
    end;
end;


end.