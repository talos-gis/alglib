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
unit htridiagonal;
interface
uses Math, Sysutils, Ap, cblas, creflections, hblas;

procedure HMatrixTD(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tau : TComplex1DArray;
     var D : TReal1DArray;
     var E : TReal1DArray);
procedure HMatrixTDUnpackQ(const A : TComplex2DArray;
     const N : AlglibInteger;
     const IsUpper : Boolean;
     const Tau : TComplex1DArray;
     var Q : TComplex2DArray);
procedure HermitianToTridiagonal(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tau : TComplex1DArray;
     var D : TReal1DArray;
     var E : TReal1DArray);
procedure UnpackQFromHermitianTridiagonal(const A : TComplex2DArray;
     const N : AlglibInteger;
     const IsUpper : Boolean;
     const Tau : TComplex1DArray;
     var Q : TComplex2DArray);

implementation

(*************************************************************************
Reduction of a Hermitian matrix which is given  by  its  higher  or  lower
triangular part to a real  tridiagonal  matrix  using  unitary  similarity
transformation: Q'*A*Q = T.

Input parameters:
    A       -   matrix to be transformed
                array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format. If IsUpper = True, then matrix A is  given
                by its upper triangle, and the lower triangle is not  used
                and not modified by the algorithm, and vice versa
                if IsUpper = False.

Output parameters:
    A       -   matrices T and Q in  compact form (see lower)
    Tau     -   array of factors which are forming matrices H(i)
                array with elements [0..N-2].
    D       -   main diagonal of real symmetric matrix T.
                array with elements [0..N-1].
    E       -   secondary diagonal of real symmetric matrix T.
                array with elements [0..N-2].


  If IsUpper=True, the matrix Q is represented as a product of elementary
  reflectors

     Q = H(n-2) . . . H(2) H(0).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a complex scalar, and v is a complex vector with
  v(i+1:n-1) = 0, v(i) = 1, v(0:i-1) is stored on exit in
  A(0:i-1,i+1), and tau in TAU(i).

  If IsUpper=False, the matrix Q is represented as a product of elementary
  reflectors

     Q = H(0) H(2) . . . H(n-2).

  Each H(i) has the form

     H(i) = I - tau * v * v'

  where tau is a complex scalar, and v is a complex vector with
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
procedure HMatrixTD(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tau : TComplex1DArray;
     var D : TReal1DArray;
     var E : TReal1DArray);
var
    I : AlglibInteger;
    Alpha : Complex;
    TauI : Complex;
    V : Complex;
    T : TComplex1DArray;
    T2 : TComplex1DArray;
    T3 : TComplex1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        Assert(AP_FP_Eq(A[I,I].Y,0));
        Inc(I);
    end;
    if N>1 then
    begin
        SetLength(Tau, N-2+1);
        SetLength(E, N-2+1);
    end;
    SetLength(D, N-1+1);
    SetLength(T, N-1+1);
    SetLength(T2, N-1+1);
    SetLength(T3, N-1+1);
    if IsUpper then
    begin
        
        //
        // Reduce the upper triangle of A
        //
        A[N-1,N-1] := C_Complex(A[N-1,N-1].X);
        I:=N-2;
        while I>=0 do
        begin
            
            //
            // Generate elementary reflector H = I+1 - tau * v * v'
            //
            ALPHA := A[I,I+1];
            T[1] := ALPHA;
            if I>=1 then
            begin
                i1_ := (0) - (2);
                for i_ := 2 to I+1 do
                begin
                    T[i_] := A[i_+i1_,I+1];
                end;
            end;
            ComplexGenerateReflection(T, I+1, TauI);
            if I>=1 then
            begin
                i1_ := (2) - (0);
                for i_ := 0 to I-1 do
                begin
                    A[i_,I+1] := T[i_+i1_];
                end;
            end;
            Alpha := T[1];
            E[I] := ALPHA.X;
            if C_NotEqualR(TAUI,0) then
            begin
                
                //
                // Apply H(I+1) from both sides to A
                //
                A[I,I+1] := C_Complex(1);
                
                //
                // Compute  x := tau * A * v  storing x in TAU
                //
                i1_ := (0) - (1);
                for i_ := 1 to I+1 do
                begin
                    T[i_] := A[i_+i1_,I+1];
                end;
                HermitianMatrixVectorMultiply(A, IsUpper, 0, I, T, TauI, T2);
                i1_ := (1) - (0);
                for i_ := 0 to I do
                begin
                    Tau[i_] := T2[i_+i1_];
                end;
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                V := C_Complex(0.0);
                for i_ := 0 to I do
                begin
                    V := C_Add(V,C_Mul(Conj(Tau[i_]),A[i_,I+1]));
                end;
                ALPHA := C_Opposite(C_Mul(C_MulR(TauI,0.5),V));
                for i_ := 0 to I do
                begin
                    Tau[i_] := C_Add(Tau[i_], C_Mul(Alpha, A[i_,I+1]));
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
                i1_ := (0) - (1);
                for i_ := 1 to I+1 do
                begin
                    T3[i_] := Tau[i_+i1_];
                end;
                HermitianRank2Update(A, IsUpper, 0, I, T, T3, T2, C_Complex(-1));
            end
            else
            begin
                A[I,I] := C_Complex(A[I,I].X);
            end;
            A[I,I+1] := C_Complex(E[I]);
            D[I+1] := A[I+1,I+1].X;
            TAU[I] := TAUI;
            Dec(I);
        end;
        D[0] := A[0,0].X;
    end
    else
    begin
        
        //
        // Reduce the lower triangle of A
        //
        A[0,0] := C_Complex(A[0,0].X);
        I:=0;
        while I<=N-2 do
        begin
            
            //
            // Generate elementary reflector H = I - tau * v * v'
            //
            i1_ := (I+1) - (1);
            for i_ := 1 to N-I-1 do
            begin
                T[i_] := A[i_+i1_,I];
            end;
            ComplexGenerateReflection(T, N-I-1, TauI);
            i1_ := (1) - (I+1);
            for i_ := I+1 to N-1 do
            begin
                A[i_,I] := T[i_+i1_];
            end;
            E[I] := A[I+1,I].X;
            if C_NotEqualR(TauI,0) then
            begin
                
                //
                // Apply H(i) from both sides to A(i+1:n,i+1:n)
                //
                A[I+1,I] := C_Complex(1);
                
                //
                // Compute  x := tau * A * v  storing y in TAU
                //
                i1_ := (I+1) - (1);
                for i_ := 1 to N-I-1 do
                begin
                    T[i_] := A[i_+i1_,I];
                end;
                HermitianMatrixVectorMultiply(A, IsUpper, I+1, N-1, T, TauI, T2);
                i1_ := (1) - (I);
                for i_ := I to N-2 do
                begin
                    Tau[i_] := T2[i_+i1_];
                end;
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                i1_ := (I+1)-(I);
                V := C_Complex(0.0);
                for i_ := I to N-2 do
                begin
                    V := C_Add(V,C_Mul(Conj(Tau[i_]),A[i_+i1_,I]));
                end;
                ALPHA := C_Opposite(C_Mul(C_MulR(TauI,0.5),V));
                i1_ := (I+1) - (I);
                for i_ := I to N-2 do
                begin
                    Tau[i_] := C_Add(Tau[i_], C_Mul(Alpha, A[i_+i1_,I]));
                end;
                
                //
                // Apply the transformation as a rank-2 update:
                // A := A - v * w' - w * v'
                //
                i1_ := (I+1) - (1);
                for i_ := 1 to N-I-1 do
                begin
                    T[i_] := A[i_+i1_,I];
                end;
                i1_ := (I) - (1);
                for i_ := 1 to N-I-1 do
                begin
                    T2[i_] := Tau[i_+i1_];
                end;
                HermitianRank2Update(A, IsUpper, I+1, N-1, T, T2, T3, C_Complex(-1));
            end
            else
            begin
                A[I+1,I+1] := C_Complex(A[I+1,I+1].X);
            end;
            A[I+1,I] := C_Complex(E[I]);
            D[I] := A[I,I].X;
            TAU[I] := TauI;
            Inc(I);
        end;
        D[N-1] := A[N-1,N-1].X;
    end;
end;


(*************************************************************************
Unpacking matrix Q which reduces a Hermitian matrix to a real  tridiagonal
form.

Input parameters:
    A       -   the result of a HMatrixTD subroutine
    N       -   size of matrix A.
    IsUpper -   storage format (a parameter of HMatrixTD subroutine)
    Tau     -   the result of a HMatrixTD subroutine

Output parameters:
    Q       -   transformation matrix.
                array with elements [0..N-1, 0..N-1].

  -- ALGLIB --
     Copyright 2005, 2007, 2008 by Bochkanov Sergey
*************************************************************************)
procedure HMatrixTDUnpackQ(const A : TComplex2DArray;
     const N : AlglibInteger;
     const IsUpper : Boolean;
     const Tau : TComplex1DArray;
     var Q : TComplex2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : TComplex1DArray;
    WORK : TComplex1DArray;
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
            V[I+1] := C_Complex(1);
            ComplexApplyReflectionFromTheLeft(Q, Tau[I], V, 0, I, 0, N-1, WORK);
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
            V[1] := C_Complex(1);
            ComplexApplyReflectionFromTheLeft(Q, Tau[I], V, I+1, N-1, 0, N-1, WORK);
            Dec(I);
        end;
    end;
end;


(*************************************************************************

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1992
*************************************************************************)
procedure HermitianToTridiagonal(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     var Tau : TComplex1DArray;
     var D : TReal1DArray;
     var E : TReal1DArray);
var
    I : AlglibInteger;
    Alpha : Complex;
    TauI : Complex;
    V : Complex;
    T : TComplex1DArray;
    T2 : TComplex1DArray;
    T3 : TComplex1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Exit;
    end;
    I:=1;
    while I<=N do
    begin
        Assert(AP_FP_Eq(A[I,I].Y,0));
        Inc(I);
    end;
    SetLength(Tau, Max(1, N-1)+1);
    SetLength(D, N+1);
    SetLength(E, Max(1, N-1)+1);
    SetLength(T, N+1);
    SetLength(T2, N+1);
    SetLength(T3, N+1);
    if IsUpper then
    begin
        
        //
        // Reduce the upper triangle of A
        //
        A[N,N] := C_Complex(A[N,N].X);
        I:=N-1;
        while I>=1 do
        begin
            
            //
            // Generate elementary reflector H(i) = I - tau * v * v'
            // to annihilate A(1:i-1,i+1)
            //
            ALPHA := A[I,I+1];
            T[1] := ALPHA;
            if I>=2 then
            begin
                i1_ := (1) - (2);
                for i_ := 2 to I do
                begin
                    T[i_] := A[i_+i1_,I+1];
                end;
            end;
            ComplexGenerateReflection(T, I, TauI);
            if I>=2 then
            begin
                i1_ := (2) - (1);
                for i_ := 1 to I-1 do
                begin
                    A[i_,I+1] := T[i_+i1_];
                end;
            end;
            Alpha := T[1];
            E[I] := ALPHA.X;
            if C_NotEqualR(TAUI,0) then
            begin
                
                //
                // Apply H(i) from both sides to A(1:i,1:i)
                //
                A[I,I+1] := C_Complex(1);
                
                //
                // Compute  x := tau * A * v  storing x in TAU(1:i)
                //
                for i_ := 1 to I do
                begin
                    T[i_] := A[i_,I+1];
                end;
                HermitianMatrixVectorMultiply(A, IsUpper, 1, I, T, TauI, Tau);
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                V := C_Complex(0.0);
                for i_ := 1 to I do
                begin
                    V := C_Add(V,C_Mul(Conj(Tau[i_]),A[i_,I+1]));
                end;
                ALPHA := C_Opposite(C_Mul(C_MulR(TauI,0.5),V));
                for i_ := 1 to I do
                begin
                    Tau[i_] := C_Add(Tau[i_], C_Mul(Alpha, A[i_,I+1]));
                end;
                
                //
                // Apply the transformation as a rank-2 update:
                //    A := A - v * w' - w * v'
                //
                for i_ := 1 to I do
                begin
                    T[i_] := A[i_,I+1];
                end;
                HermitianRank2Update(A, IsUpper, 1, I, T, Tau, T2, C_Complex(-1));
            end
            else
            begin
                A[I,I] := C_Complex(A[I,I].X);
            end;
            A[I,I+1] := C_Complex(E[I]);
            D[I+1] := A[I+1,I+1].X;
            TAU[I] := TAUI;
            Dec(I);
        end;
        D[1] := A[1,1].X;
    end
    else
    begin
        
        //
        // Reduce the lower triangle of A
        //
        A[1,1] := C_Complex(A[1,1].X);
        I:=1;
        while I<=N-1 do
        begin
            
            //
            // Generate elementary reflector H(i) = I - tau * v * v'
            // to annihilate A(i+2:n,i)
            //
            i1_ := (I+1) - (1);
            for i_ := 1 to N-I do
            begin
                T[i_] := A[i_+i1_,I];
            end;
            ComplexGenerateReflection(T, N-I, TauI);
            i1_ := (1) - (I+1);
            for i_ := I+1 to N do
            begin
                A[i_,I] := T[i_+i1_];
            end;
            E[I] := A[I+1,I].X;
            if C_NotEqualR(TauI,0) then
            begin
                
                //
                // Apply H(i) from both sides to A(i+1:n,i+1:n)
                //
                A[I+1,I] := C_Complex(1);
                
                //
                // Compute  x := tau * A * v  storing y in TAU(i:n-1)
                //
                i1_ := (I+1) - (1);
                for i_ := 1 to N-I do
                begin
                    T[i_] := A[i_+i1_,I];
                end;
                HermitianMatrixVectorMultiply(A, IsUpper, I+1, N, T, TauI, T2);
                i1_ := (1) - (I);
                for i_ := I to N-1 do
                begin
                    Tau[i_] := T2[i_+i1_];
                end;
                
                //
                // Compute  w := x - 1/2 * tau * (x'*v) * v
                //
                i1_ := (I+1)-(I);
                V := C_Complex(0.0);
                for i_ := I to N-1 do
                begin
                    V := C_Add(V,C_Mul(Conj(Tau[i_]),A[i_+i1_,I]));
                end;
                ALPHA := C_Opposite(C_Mul(C_MulR(TauI,0.5),V));
                i1_ := (I+1) - (I);
                for i_ := I to N-1 do
                begin
                    Tau[i_] := C_Add(Tau[i_], C_Mul(Alpha, A[i_+i1_,I]));
                end;
                
                //
                // Apply the transformation as a rank-2 update:
                // A := A - v * w' - w * v'
                //
                i1_ := (I+1) - (1);
                for i_ := 1 to N-I do
                begin
                    T[i_] := A[i_+i1_,I];
                end;
                i1_ := (I) - (1);
                for i_ := 1 to N-I do
                begin
                    T2[i_] := Tau[i_+i1_];
                end;
                HermitianRank2Update(A, IsUpper, I+1, N, T, T2, T3, C_Complex(-1));
            end
            else
            begin
                A[I+1,I+1] := C_Complex(A[I+1,I+1].X);
            end;
            A[I+1,I] := C_Complex(E[I]);
            D[I] := A[I,I].X;
            TAU[I] := TauI;
            Inc(I);
        end;
        D[N] := A[N,N].X;
    end;
end;


(*************************************************************************

  -- ALGLIB --
     Copyright 2005, 2007 by Bochkanov Sergey
*************************************************************************)
procedure UnpackQFromHermitianTridiagonal(const A : TComplex2DArray;
     const N : AlglibInteger;
     const IsUpper : Boolean;
     const Tau : TComplex1DArray;
     var Q : TComplex2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : TComplex1DArray;
    WORK : TComplex1DArray;
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
    if IsUpper then
    begin
        I:=1;
        while I<=N-1 do
        begin
            
            //
            // Apply H(i)
            //
            for i_ := 1 to I do
            begin
                V[i_] := A[i_,I+1];
            end;
            V[I] := C_Complex(1);
            ComplexApplyReflectionFromTheLeft(Q, Tau[I], V, 1, I, 1, N, WORK);
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
            i1_ := (I+1) - (1);
            for i_ := 1 to N-I do
            begin
                V[i_] := A[i_+i1_,I];
            end;
            V[1] := C_Complex(1);
            ComplexApplyReflectionFromTheLeft(Q, Tau[I], V, I+1, N, 1, N, WORK);
            Dec(I);
        end;
    end;
end;


end.