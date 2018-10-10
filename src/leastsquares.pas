(*************************************************************************
Copyright (c) 2006-2007, Sergey Bochkanov (ALGLIB project).

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
unit leastsquares;
interface
uses Math, Sysutils, Ap, spline3, reflections, lq, bidiagonal, rotations, bdsvd, qr, blas, svd;

procedure BuildGeneralLeastSquares(const Y : TReal1DArray;
     const W : TReal1DArray;
     const FMatrix : TReal2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var C : TReal1DArray);
procedure BuildLinearLeastSquares(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     var a : Double;
     var b : Double);
procedure BuildSplineLeastSquares(const X : TReal1DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     A : Double;
     B : Double;
     N : AlglibInteger;
     M : AlglibInteger;
     var CTbl : TReal1DArray);
procedure BuildPolynomialLeastSquares(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var C : TReal1DArray);
procedure BuildChebyshevLeastSquares(const X : TReal1DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     A : Double;
     B : Double;
     N : AlglibInteger;
     M : AlglibInteger;
     var CTbl : TReal1DArray);
function BuildChebyshevLeastSquaresConstrained(const X : TReal1DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     A : Double;
     B : Double;
     N : AlglibInteger;
     const XC : TReal1DArray;
     const YC : TReal1DArray;
     const DC : TInteger1DArray;
     NC : AlglibInteger;
     M : AlglibInteger;
     var CTbl : TReal1DArray):Boolean;
function CalculateChebyshevLeastSquares(const M : AlglibInteger;
     const A : TReal1DArray;
     X : Double):Double;

implementation

(*************************************************************************
Weighted approximation by arbitrary function basis in a space of arbitrary
dimension using linear least squares method.

Input parameters:
    Y   -   array[0..N-1]
            It contains a set  of  function  values  in  N  points.  Space
            dimension  and  points  don't  matter.  Procedure  works  with
            function values in these points and values of basis  functions
            only.

    W   -   array[0..N-1]
            It contains weights corresponding  to  function  values.  Each
            summand in square sum of approximation deviations  from  given
            values is multiplied by the square of corresponding weight.

    FMatrix-a table of basis functions values, array[0..N-1, 0..M-1].
            FMatrix[I, J] - value of J-th basis function in I-th point.

    N   -   number of points used. N>=1.
    M   -   number of basis functions, M>=1.

Output parameters:
    C   -   decomposition coefficients.
            Array of real numbers whose index goes from 0 to M-1.
            C[j] - j-th basis function coefficient.

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************)
procedure BuildGeneralLeastSquares(const Y : TReal1DArray;
     const W : TReal1DArray;
     const FMatrix : TReal2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var C : TReal1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    A : TReal2DArray;
    Q : TReal2DArray;
    VT : TReal2DArray;
    B : TReal1DArray;
    Tau : TReal1DArray;
    B2 : TReal2DArray;
    TauQ : TReal1DArray;
    TauP : TReal1DArray;
    D : TReal1DArray;
    E : TReal1DArray;
    IsUpperA : Boolean;
    MI : AlglibInteger;
    NI : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    MI := N;
    NI := M;
    SetLength(C, NI-1+1);
    
    //
    // Initialize design matrix.
    // Here we are making MI>=NI.
    //
    SetLength(A, NI+1, Max(MI, NI)+1);
    SetLength(B, Max(MI, NI)+1);
    I:=1;
    while I<=MI do
    begin
        B[I] := W[I-1]*Y[I-1];
        Inc(I);
    end;
    I:=MI+1;
    while I<=NI do
    begin
        B[I] := 0;
        Inc(I);
    end;
    J:=1;
    while J<=NI do
    begin
        i1_ := (0) - (1);
        for i_ := 1 to MI do
        begin
            A[J,i_] := FMatrix[i_+i1_,J-1];
        end;
        Inc(J);
    end;
    J:=1;
    while J<=NI do
    begin
        I:=MI+1;
        while I<=NI do
        begin
            A[J,I] := 0;
            Inc(I);
        end;
        Inc(J);
    end;
    J:=1;
    while J<=NI do
    begin
        I:=1;
        while I<=MI do
        begin
            A[J,I] := A[J,I]*W[I-1];
            Inc(I);
        end;
        Inc(J);
    end;
    MI := Max(MI, NI);
    
    //
    // LQ-decomposition of A'
    // B2 := Q*B
    //
    LQDecomposition(A, NI, MI, Tau);
    UnpackQFromLQ(A, NI, MI, Tau, NI, Q);
    SetLength(B2, 1+1, NI+1);
    J:=1;
    while J<=NI do
    begin
        B2[1,J] := 0;
        Inc(J);
    end;
    I:=1;
    while I<=NI do
    begin
        V := APVDotProduct(@B[0], 1, MI, @Q[I][0], 1, MI);
        B2[1,I] := V;
        Inc(I);
    end;
    
    //
    // Back from A' to A
    // Making cols(A)=rows(A)
    //
    I:=1;
    while I<=NI-1 do
    begin
        for i_ := I+1 to NI do
        begin
            A[I,i_] := A[i_,I];
        end;
        Inc(I);
    end;
    I:=2;
    while I<=NI do
    begin
        J:=1;
        while J<=I-1 do
        begin
            A[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Bidiagonal decomposition of A
    // A = Q * d2 * P'
    // B2 := (Q'*B2')'
    //
    ToBidiagonal(A, NI, NI, TauQ, TauP);
    MultiplyByQFromBidiagonal(A, NI, NI, TauQ, B2, 1, NI, True, False);
    UnpackPTFromBidiagonal(A, NI, NI, TauP, NI, VT);
    UnpackDiagonalsFromBidiagonal(A, NI, NI, IsUpperA, D, E);
    
    //
    // Singular value decomposition of A
    // A = U * d * V'
    // B2 := (U'*B2')'
    //
    if  not BidiagonalSVDDecomposition(D, E, NI, IsUpperA, False, B2, 1, Q, 0, VT, NI) then
    begin
        I:=0;
        while I<=NI-1 do
        begin
            C[I] := 0;
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // B2 := (d^(-1) * B2')'
    //
    if AP_FP_Neq(D[1],0) then
    begin
        I:=1;
        while I<=NI do
        begin
            if AP_FP_Greater(D[I],MachineEpsilon*10*Sqrt(NI)*D[1]) then
            begin
                B2[1,I] := B2[1,I]/D[I];
            end
            else
            begin
                B2[1,I] := 0;
            end;
            Inc(I);
        end;
    end;
    
    //
    // B := (V * B2')'
    //
    I:=1;
    while I<=NI do
    begin
        B[I] := 0;
        Inc(I);
    end;
    I:=1;
    while I<=NI do
    begin
        V := B2[1,I];
        APVAdd(@B[0], 1, NI, @VT[I][0], 1, NI, V);
        Inc(I);
    end;
    
    //
    // Out
    //
    I:=1;
    while I<=NI do
    begin
        C[I-1] := B[I];
        Inc(I);
    end;
end;


(*************************************************************************
Linear approximation using least squares method

The subroutine calculates coefficients of  the  line  approximating  given
function.

Input parameters:
    X   -   array[0..N-1], it contains a set of abscissas.
    Y   -   array[0..N-1], function values.
    N   -   number of points, N>=1

Output parameters:
    a, b-   coefficients of linear approximation a+b*t

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************)
procedure BuildLinearLeastSquares(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     var a : Double;
     var b : Double);
var
    PP : Double;
    QQ : Double;
    PQ : Double;
    b1 : Double;
    b2 : Double;
    d1 : Double;
    d2 : Double;
    t1 : Double;
    t2 : Double;
    Phi : Double;
    c : Double;
    s : Double;
    m : Double;
    I : AlglibInteger;
begin
    PP := N;
    QQ := 0;
    PQ := 0;
    b1 := 0;
    b2 := 0;
    I:=0;
    while I<=N-1 do
    begin
        PQ := PQ+X[I];
        QQ := QQ+AP_Sqr(X[I]);
        b1 := B1+Y[I];
        b2 := B2+X[I]*Y[I];
        Inc(I);
    end;
    Phi := ArcTan2(2*PQ, QQ-PP)/2;
    c := Cos(Phi);
    s := Sin(Phi);
    d1 := AP_Sqr(c)*pp+AP_Sqr(s)*qq-2*s*c*pq;
    d2 := AP_Sqr(s)*pp+AP_Sqr(c)*qq+2*s*c*pq;
    if AP_FP_Greater(AbsReal(d1),AbsReal(d2)) then
    begin
        m := AbsReal(d1);
    end
    else
    begin
        m := AbsReal(d2);
    end;
    t1 := c*b1-s*b2;
    t2 := s*b1+c*b2;
    if AP_FP_Greater(AbsReal(d1),m*MachineEpsilon*1000) then
    begin
        t1 := t1/d1;
    end
    else
    begin
        t1 := 0;
    end;
    if AP_FP_Greater(AbsReal(d2),m*MachineEpsilon*1000) then
    begin
        t2 := t2/d2;
    end
    else
    begin
        t2 := 0;
    end;
    a := c*t1+s*t2;
    b := -s*t1+c*t2;
end;


(*************************************************************************
Weighted cubic spline approximation using linear least squares

Input parameters:
    X   -   array[0..N-1], abscissas
    Y   -   array[0..N-1], function values
    W   -   array[0..N-1], weights.
    A, B-   interval to build splines in.
    N   -   number of points used. N>=1.
    M   -   number of basic splines, M>=2.

Output parameters:
    CTbl-   coefficients table to be used by SplineInterpolation function.
  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************)
procedure BuildSplineLeastSquares(const X : TReal1DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     A : Double;
     B : Double;
     N : AlglibInteger;
     M : AlglibInteger;
     var CTbl : TReal1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    MA : TReal2DArray;
    Q : TReal2DArray;
    VT : TReal2DArray;
    MB : TReal1DArray;
    Tau : TReal1DArray;
    B2 : TReal2DArray;
    TauQ : TReal1DArray;
    TauP : TReal1DArray;
    D : TReal1DArray;
    E : TReal1DArray;
    IsUpperA : Boolean;
    MI : AlglibInteger;
    NI : AlglibInteger;
    V : Double;
    SX : TReal1DArray;
    SY : TReal1DArray;
    i_ : AlglibInteger;
begin
    Assert(M>=2, 'BuildSplineLeastSquares: M is too small!');
    MI := N;
    NI := M;
    SetLength(SX, NI-1+1);
    SetLength(SY, NI-1+1);
    
    //
    // Initializing design matrix
    // Here we are making MI>=NI
    //
    SetLength(MA, NI+1, Max(MI, NI)+1);
    SetLength(MB, Max(MI, NI)+1);
    J:=0;
    while J<=NI-1 do
    begin
        SX[J] := A+(B-A)*J/(NI-1);
        Inc(J);
    end;
    J:=0;
    while J<=NI-1 do
    begin
        I:=0;
        while I<=NI-1 do
        begin
            SY[I] := 0;
            Inc(I);
        end;
        SY[J] := 1;
        BuildCubicSpline(SX, SY, NI, 0, 0.0, 0, 0.0, CTbl);
        I:=0;
        while I<=MI-1 do
        begin
            MA[J+1,I+1] := W[I]*SplineInterpolation(CTbl, X[I]);
            Inc(I);
        end;
        Inc(J);
    end;
    J:=1;
    while J<=NI do
    begin
        I:=MI+1;
        while I<=NI do
        begin
            MA[J,I] := 0;
            Inc(I);
        end;
        Inc(J);
    end;
    
    //
    // Initializing right part
    //
    I:=0;
    while I<=MI-1 do
    begin
        MB[I+1] := W[I]*Y[I];
        Inc(I);
    end;
    I:=MI+1;
    while I<=NI do
    begin
        MB[I] := 0;
        Inc(I);
    end;
    MI := Max(MI, NI);
    
    //
    // LQ-decomposition of A'
    // B2 := Q*B
    //
    LQDecomposition(MA, NI, MI, Tau);
    UnpackQFromLQ(MA, NI, MI, Tau, NI, Q);
    SetLength(B2, 1+1, NI+1);
    J:=1;
    while J<=NI do
    begin
        B2[1,J] := 0;
        Inc(J);
    end;
    I:=1;
    while I<=NI do
    begin
        V := APVDotProduct(@MB[0], 1, MI, @Q[I][0], 1, MI);
        B2[1,I] := V;
        Inc(I);
    end;
    
    //
    // Back from A' to A
    // Making cols(A)=rows(A)
    //
    I:=1;
    while I<=NI-1 do
    begin
        for i_ := I+1 to NI do
        begin
            MA[I,i_] := MA[i_,I];
        end;
        Inc(I);
    end;
    I:=2;
    while I<=NI do
    begin
        J:=1;
        while J<=I-1 do
        begin
            MA[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Bidiagonal decomposition of A
    // A = Q * d2 * P'
    // B2 := (Q'*B2')'
    //
    ToBidiagonal(MA, NI, NI, TauQ, TauP);
    MultiplyByQFromBidiagonal(MA, NI, NI, TauQ, B2, 1, NI, True, False);
    UnpackPTFromBidiagonal(MA, NI, NI, TauP, NI, VT);
    UnpackDiagonalsFromBidiagonal(MA, NI, NI, IsUpperA, D, E);
    
    //
    // Singular value decomposition of A
    // A = U * d * V'
    // B2 := (U'*B2')'
    //
    if  not BidiagonalSVDDecomposition(D, E, NI, IsUpperA, False, B2, 1, Q, 0, VT, NI) then
    begin
        I:=1;
        while I<=NI do
        begin
            D[I] := 0;
            B2[1,I] := 0;
            J:=1;
            while J<=NI do
            begin
                if I=J then
                begin
                    VT[I,J] := 1;
                end
                else
                begin
                    VT[I,J] := 0;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        B2[1,1] := 1;
    end;
    
    //
    // B2 := (d^(-1) * B2')'
    //
    I:=1;
    while I<=NI do
    begin
        if AP_FP_Greater(D[I],MachineEpsilon*10*Sqrt(NI)*D[1]) then
        begin
            B2[1,I] := B2[1,I]/D[I];
        end
        else
        begin
            B2[1,I] := 0;
        end;
        Inc(I);
    end;
    
    //
    // B := (V * B2')'
    //
    I:=1;
    while I<=NI do
    begin
        MB[I] := 0;
        Inc(I);
    end;
    I:=1;
    while I<=NI do
    begin
        V := B2[1,I];
        APVAdd(@MB[0], 1, NI, @VT[I][0], 1, NI, V);
        Inc(I);
    end;
    
    //
    // Forming result spline
    //
    I:=0;
    while I<=NI-1 do
    begin
        SY[I] := MB[I+1];
        Inc(I);
    end;
    BuildCubicSpline(SX, SY, NI, 0, 0.0, 0, 0.0, CTbl);
end;


(*************************************************************************
Polynomial approximation using least squares method

The subroutine calculates coefficients  of  the  polynomial  approximating
given function. It is recommended to use this function only if you need to
obtain coefficients of approximation polynomial. If you have to build  and
calculate polynomial approximation (NOT coefficients), it's better to  use
BuildChebyshevLeastSquares      subroutine     in     combination     with
CalculateChebyshevLeastSquares   subroutine.   The   result  of  Chebyshev
polynomial approximation is equivalent to the result obtained using powers
of X, but has higher  accuracy  due  to  better  numerical  properties  of
Chebyshev polynomials.

Input parameters:
    X   -   array[0..N-1], abscissas
    Y   -   array[0..N-1], function values
    N   -   number of points, N>=1
    M   -   order of polynomial required, M>=0

Output parameters:
    C   -   approximating polynomial coefficients, array[0..M],
            C[i] - coefficient at X^i.

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************)
procedure BuildPolynomialLeastSquares(const X : TReal1DArray;
     const Y : TReal1DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     var C : TReal1DArray);
var
    CTbl : TReal1DArray;
    W : TReal1DArray;
    C1 : TReal1DArray;
    MaxX : Double;
    MinX : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    E : Double;
    D : Double;
    L1 : Double;
    L2 : Double;
    Z2 : TReal1DArray;
    Z1 : TReal1DArray;
begin
    
    //
    // Initialize
    //
    MaxX := X[0];
    MinX := X[0];
    I:=1;
    while I<=N-1 do
    begin
        if AP_FP_Greater(X[I],MaxX) then
        begin
            MaxX := X[I];
        end;
        if AP_FP_Less(X[I],MinX) then
        begin
            MinX := X[I];
        end;
        Inc(I);
    end;
    if AP_FP_Eq(MinX,MaxX) then
    begin
        MinX := MinX-0.5;
        MaxX := MaxX+0.5;
    end;
    SetLength(W, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        W[I] := 1;
        Inc(I);
    end;
    
    //
    // Build Chebyshev approximation
    //
    BuildChebyshevLeastSquares(X, Y, W, MinX, MaxX, N, M, CTbl);
    
    //
    // From Chebyshev to powers of X
    //
    SetLength(C1, M+1);
    I:=0;
    while I<=M do
    begin
        C1[i] := 0;
        Inc(I);
    end;
    d := 0;
    i:=0;
    while i<=M do
    begin
        k:=i;
        while k<=M do
        begin
            e := C1[k];
            C1[k] := 0;
            if (i<=1) and (k=i) then
            begin
                C1[k] := 1;
            end
            else
            begin
                if i<>0 then
                begin
                    C1[k] := 2*d;
                end;
                if k>i+1 then
                begin
                    C1[k] := C1[k]-C1[k-2];
                end;
            end;
            d := e;
            Inc(k);
        end;
        d := C1[i];
        e := 0;
        k := i;
        while k<=M do
        begin
            e := e+C1[k]*CTbl[k];
            k := k+2;
        end;
        C1[i] := e;
        Inc(i);
    end;
    
    //
    // Linear translation
    //
    L1 := 2/(CTbl[M+2]-CTbl[M+1]);
    L2 := -2*CTbl[M+1]/(CTbl[M+2]-CTbl[M+1])-1;
    SetLength(C, M+1);
    SetLength(Z2, M+1);
    SetLength(Z1, M+1);
    C[0] := C1[0];
    Z1[0] := 1;
    Z2[0] := 1;
    i:=1;
    while i<=M do
    begin
        Z2[i] := 1;
        Z1[i] := L2*Z1[i-1];
        C[0] := C[0]+C1[i]*Z1[i];
        Inc(i);
    end;
    j:=1;
    while j<=M do
    begin
        Z2[0] := L1*Z2[0];
        C[j] := C1[j]*Z2[0];
        i:=j+1;
        while i<=M do
        begin
            k := i-j;
            Z2[k] := L1*Z2[k]+Z2[k-1];
            C[j] := C[j]+C1[i]*Z2[k]*Z1[k];
            Inc(i);
        end;
        Inc(j);
    end;
end;


(*************************************************************************
Chebyshev polynomial approximation using least squares method.

The algorithm reduces interval [A, B] to the interval [-1,1], then  builds
least squares approximation using Chebyshev polynomials.

Input parameters:
    X   -   array[0..N-1], abscissas
    Y   -   array[0..N-1], function values
    W   -   array[0..N-1], weights
    A, B-   interval to build approximating polynomials in.
    N   -   number of points used. N>=1.
    M   -   order of polynomial, M>=0. This parameter is passed into
            CalculateChebyshevLeastSquares function.

Output parameters:
    CTbl - coefficient table. This parameter is passed into
            CalculateChebyshevLeastSquares function.
  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************)
procedure BuildChebyshevLeastSquares(const X : TReal1DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     A : Double;
     B : Double;
     N : AlglibInteger;
     M : AlglibInteger;
     var CTbl : TReal1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    MA : TReal2DArray;
    Q : TReal2DArray;
    VT : TReal2DArray;
    MB : TReal1DArray;
    Tau : TReal1DArray;
    B2 : TReal2DArray;
    TauQ : TReal1DArray;
    TauP : TReal1DArray;
    D : TReal1DArray;
    E : TReal1DArray;
    IsUpperA : Boolean;
    MI : AlglibInteger;
    NI : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
begin
    MI := N;
    NI := M+1;
    
    //
    // Initializing design matrix
    // Here we are making MI>=NI
    //
    SetLength(MA, NI+1, Max(MI, NI)+1);
    SetLength(MB, Max(MI, NI)+1);
    J:=1;
    while J<=NI do
    begin
        I:=1;
        while I<=MI do
        begin
            V := 2*(X[I-1]-A)/(B-A)-1;
            if J=1 then
            begin
                MA[J,I] := 1.0;
            end;
            if J=2 then
            begin
                MA[J,I] := V;
            end;
            if J>2 then
            begin
                MA[J,I] := 2.0*V*MA[J-1,I]-MA[J-2,I];
            end;
            Inc(I);
        end;
        Inc(J);
    end;
    J:=1;
    while J<=NI do
    begin
        I:=1;
        while I<=MI do
        begin
            MA[J,I] := W[I-1]*MA[J,I];
            Inc(I);
        end;
        Inc(J);
    end;
    J:=1;
    while J<=NI do
    begin
        I:=MI+1;
        while I<=NI do
        begin
            MA[J,I] := 0;
            Inc(I);
        end;
        Inc(J);
    end;
    
    //
    // Initializing right part
    //
    I:=0;
    while I<=MI-1 do
    begin
        MB[I+1] := W[I]*Y[I];
        Inc(I);
    end;
    I:=MI+1;
    while I<=NI do
    begin
        MB[I] := 0;
        Inc(I);
    end;
    MI := Max(MI, NI);
    
    //
    // LQ-decomposition of A'
    // B2 := Q*B
    //
    LQDecomposition(MA, NI, MI, Tau);
    UnpackQFromLQ(MA, NI, MI, Tau, NI, Q);
    SetLength(B2, 1+1, NI+1);
    J:=1;
    while J<=NI do
    begin
        B2[1,J] := 0;
        Inc(J);
    end;
    I:=1;
    while I<=NI do
    begin
        V := APVDotProduct(@MB[0], 1, MI, @Q[I][0], 1, MI);
        B2[1,I] := V;
        Inc(I);
    end;
    
    //
    // Back from A' to A
    // Making cols(A)=rows(A)
    //
    I:=1;
    while I<=NI-1 do
    begin
        for i_ := I+1 to NI do
        begin
            MA[I,i_] := MA[i_,I];
        end;
        Inc(I);
    end;
    I:=2;
    while I<=NI do
    begin
        J:=1;
        while J<=I-1 do
        begin
            MA[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Bidiagonal decomposition of A
    // A = Q * d2 * P'
    // B2 := (Q'*B2')'
    //
    ToBidiagonal(MA, NI, NI, TauQ, TauP);
    MultiplyByQFromBidiagonal(MA, NI, NI, TauQ, B2, 1, NI, True, False);
    UnpackPTFromBidiagonal(MA, NI, NI, TauP, NI, VT);
    UnpackDiagonalsFromBidiagonal(MA, NI, NI, IsUpperA, D, E);
    
    //
    // Singular value decomposition of A
    // A = U * d * V'
    // B2 := (U'*B2')'
    //
    if  not BidiagonalSVDDecomposition(D, E, NI, IsUpperA, False, B2, 1, Q, 0, VT, NI) then
    begin
        I:=1;
        while I<=NI do
        begin
            D[I] := 0;
            B2[1,I] := 0;
            J:=1;
            while J<=NI do
            begin
                if I=J then
                begin
                    VT[I,J] := 1;
                end
                else
                begin
                    VT[I,J] := 0;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        B2[1,1] := 1;
    end;
    
    //
    // B2 := (d^(-1) * B2')'
    //
    I:=1;
    while I<=NI do
    begin
        if AP_FP_Greater(D[I],MachineEpsilon*10*Sqrt(NI)*D[1]) then
        begin
            B2[1,I] := B2[1,I]/D[I];
        end
        else
        begin
            B2[1,I] := 0;
        end;
        Inc(I);
    end;
    
    //
    // B := (V * B2')'
    //
    I:=1;
    while I<=NI do
    begin
        MB[I] := 0;
        Inc(I);
    end;
    I:=1;
    while I<=NI do
    begin
        V := B2[1,I];
        APVAdd(@MB[0], 1, NI, @VT[I][0], 1, NI, V);
        Inc(I);
    end;
    
    //
    // Forming result
    //
    SetLength(CTbl, NI+1+1);
    I:=1;
    while I<=NI do
    begin
        CTbl[I-1] := MB[I];
        Inc(I);
    end;
    CTbl[NI] := A;
    CTbl[NI+1] := B;
end;


(*************************************************************************
Weighted Chebyshev polynomial constrained least squares approximation.

The algorithm reduces [A,B] to [-1,1] and builds the Chebyshev polynomials
series by approximating a given function using the least squares method.

Input parameters:
    X   -   abscissas, array[0..N-1]
    Y   -   function values, array[0..N-1]
    W   -   weights, array[0..N-1].  Each  item  in  the  squared  sum  of
            deviations from given values is  multiplied  by  a  square  of
            corresponding weight.
    A, B-   interval in which the approximating polynomials are built.
    N   -   number of points, N>0.
    XC, YC, DC-
            constraints (see description below)., array[0..NC-1]
    NC  -   number of constraints. 0 <= NC < M+1.
    M   -   degree of polynomial, M>=0. This parameter is passed into  the
            CalculateChebyshevLeastSquares subroutine.

Output parameters:
    CTbl-   coefficient  table.  This  parameter  is   passed   into   the
            CalculateChebyshevLeastSquares subroutine.

Result:
    True, if the algorithm succeeded.
    False, if the internal singular value decomposition subroutine  hasn't
converged or the given constraints could not be met  simultaneously  (e.g.
P(0)=0 è P(0)=1).

Specifying constraints:
    This subroutine can solve  the  problem  having  constrained  function
values or its derivatives in several points. NC specifies  the  number  of
constraints, DC - the type of constraints, XC and YC - constraints as such.
Thus, for each i from 0 to NC-1 the following constraint is given:
    P(xc[i]) = yc[i],       if DC[i]=0
or
    d/dx(P(xc[i])) = yc[i], if DC[i]=1
(here P(x) is approximating polynomial).
    This version of the subroutine supports only either polynomial or  its
derivative value constraints.  If  DC[i]  is  not  equal  to  0 and 1, the
subroutine will be aborted. The number of constraints should be less  than
the number of degrees of freedom of approximating  polynomial  -  M+1  (at
that, it could be equal to 0).

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************)
function BuildChebyshevLeastSquaresConstrained(const X : TReal1DArray;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     A : Double;
     B : Double;
     N : AlglibInteger;
     const XC : TReal1DArray;
     const YC : TReal1DArray;
     const DC : TInteger1DArray;
     NC : AlglibInteger;
     M : AlglibInteger;
     var CTbl : TReal1DArray):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    ReducedSize : AlglibInteger;
    DesignMatrix : TReal2DArray;
    RightPart : TReal1DArray;
    CMatrix : TReal2DArray;
    C : TReal2DArray;
    U : TReal2DArray;
    VT : TReal2DArray;
    D : TReal1DArray;
    CR : TReal1DArray;
    WS : TReal1DArray;
    TJ : TReal1DArray;
    UJ : TReal1DArray;
    DTJ : TReal1DArray;
    Tmp : TReal1DArray;
    Tmp2 : TReal1DArray;
    TmpMatrix : TReal2DArray;
    V : Double;
    i_ : AlglibInteger;
begin
    Assert(N>0);
    Assert(M>=0);
    Assert((NC>=0) and (NC<M+1));
    Result := True;
    
    //
    // Initialize design matrix and right part.
    // Add fictional rows if needed to ensure that N>=M+1.
    //
    SetLength(DesignMatrix, Max(N, M+1)+1, M+1+1);
    SetLength(RightPart, Max(N, M+1)+1);
    I:=1;
    while I<=N do
    begin
        J:=1;
        while J<=M+1 do
        begin
            V := 2*(X[I-1]-A)/(B-A)-1;
            if J=1 then
            begin
                DesignMatrix[I,J] := 1.0;
            end;
            if J=2 then
            begin
                DesignMatrix[I,J] := V;
            end;
            if J>2 then
            begin
                DesignMatrix[I,J] := 2.0*V*DesignMatrix[I,J-1]-DesignMatrix[I,J-2];
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=1;
    while I<=N do
    begin
        J:=1;
        while J<=M+1 do
        begin
            DesignMatrix[I,J] := W[I-1]*DesignMatrix[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
    I:=N+1;
    while I<=M+1 do
    begin
        J:=1;
        while J<=M+1 do
        begin
            DesignMatrix[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        RightPart[I+1] := W[I]*Y[I];
        Inc(I);
    end;
    I:=N+1;
    while I<=M+1 do
    begin
        RightPart[I] := 0;
        Inc(I);
    end;
    N := Max(N, M+1);
    
    //
    // Now N>=M+1 and we are ready to the next stage.
    // Handle constraints.
    // Represent feasible set of coefficients as x = C*t + d
    //
    SetLength(C, M+1+1, M+1+1);
    SetLength(D, M+1+1);
    if NC=0 then
    begin
        
        //
        // No constraints
        //
        I:=1;
        while I<=M+1 do
        begin
            J:=1;
            while J<=M+1 do
            begin
                C[I,J] := 0;
                Inc(J);
            end;
            D[I] := 0;
            Inc(I);
        end;
        I:=1;
        while I<=M+1 do
        begin
            C[I,I] := 1;
            Inc(I);
        end;
        ReducedSize := M+1;
    end
    else
    begin
        
        //
        // Constraints are present.
        // Fill constraints matrix CMatrix and solve CMatrix*x = cr.
        //
        SetLength(CMatrix, NC+1, M+1+1);
        SetLength(CR, NC+1);
        SetLength(TJ, M+1);
        SetLength(UJ, M+1);
        SetLength(DTJ, M+1);
        I:=0;
        while I<=NC-1 do
        begin
            V := 2*(XC[I]-A)/(B-A)-1;
            J:=0;
            while J<=M do
            begin
                if J=0 then
                begin
                    TJ[J] := 1;
                    UJ[J] := 1;
                    DTJ[J] := 0;
                end;
                if J=1 then
                begin
                    TJ[J] := V;
                    UJ[J] := 2*V;
                    DTJ[J] := 1;
                end;
                if J>1 then
                begin
                    TJ[J] := 2*V*TJ[J-1]-TJ[J-2];
                    UJ[J] := 2*V*UJ[J-1]-UJ[J-2];
                    DTJ[J] := J*UJ[J-1];
                end;
                Assert((DC[I]=0) or (DC[I]=1));
                if DC[I]=0 then
                begin
                    CMatrix[I+1,J+1] := TJ[J];
                end;
                if DC[I]=1 then
                begin
                    CMatrix[I+1,J+1] := DTJ[J];
                end;
                Inc(J);
            end;
            CR[I+1] := YC[I];
            Inc(I);
        end;
        
        //
        // Solve CMatrix*x = cr.
        // Fill C and d:
        // 1. SVD: CMatrix = U * WS * V^T
        // 2. C := V[1:M+1,NC+1:M+1]
        // 3. tmp := WS^-1 * U^T * cr
        // 4. d := V[1:M+1,1:NC] * tmp
        //
        if  not SVDDecomposition(CMatrix, NC, M+1, 2, 2, 2, WS, U, VT) then
        begin
            Result := False;
            Exit;
        end;
        if AP_FP_Eq(WS[1],0) or AP_FP_Less_Eq(WS[NC],MachineEpsilon*10*Sqrt(NC)*WS[1]) then
        begin
            Result := False;
            Exit;
        end;
        SetLength(C, M+1+1, M+1-NC+1);
        SetLength(D, M+1+1);
        I:=1;
        while I<=M+1-NC do
        begin
            for i_ := 1 to M+1 do
            begin
                C[i_,I] := VT[NC+I,i_];
            end;
            Inc(I);
        end;
        SetLength(Tmp, NC+1);
        I:=1;
        while I<=NC do
        begin
            V := 0.0;
            for i_ := 1 to NC do
            begin
                V := V + U[i_,I]*CR[i_];
            end;
            Tmp[I] := V/WS[I];
            Inc(I);
        end;
        I:=1;
        while I<=M+1 do
        begin
            D[I] := 0;
            Inc(I);
        end;
        I:=1;
        while I<=NC do
        begin
            V := Tmp[I];
            APVAdd(@D[0], 1, M+1, @VT[I][0], 1, M+1, V);
            Inc(I);
        end;
        
        //
        // Reduce problem:
        // 1. RightPart := RightPart - DesignMatrix*d
        // 2. DesignMatrix := DesignMatrix*C
        //
        I:=1;
        while I<=N do
        begin
            V := APVDotProduct(@DesignMatrix[I][0], 1, M+1, @D[0], 1, M+1);
            RightPart[I] := RightPart[I]-V;
            Inc(I);
        end;
        ReducedSize := M+1-NC;
        SetLength(TmpMatrix, N+1, ReducedSize+1);
        SetLength(Tmp, N+1);
        MatrixMatrixMultiply(DesignMatrix, 1, N, 1, M+1, False, C, 1, M+1, 1, ReducedSize, False, 1.0, TmpMatrix, 1, N, 1, ReducedSize, 0.0, Tmp);
        CopyMatrix(TmpMatrix, 1, N, 1, ReducedSize, DesignMatrix, 1, N, 1, ReducedSize);
    end;
    
    //
    // Solve reduced problem DesignMatrix*t = RightPart.
    //
    if  not SVDDecomposition(DesignMatrix, N, ReducedSize, 1, 1, 2, WS, U, VT) then
    begin
        Result := False;
        Exit;
    end;
    SetLength(Tmp, ReducedSize+1);
    SetLength(Tmp2, ReducedSize+1);
    I:=1;
    while I<=ReducedSize do
    begin
        Tmp[I] := 0;
        Inc(I);
    end;
    I:=1;
    while I<=N do
    begin
        V := RightPart[I];
        APVAdd(@Tmp[0], 1, ReducedSize, @U[I][0], 1, ReducedSize, V);
        Inc(I);
    end;
    I:=1;
    while I<=ReducedSize do
    begin
        if AP_FP_Neq(WS[I],0) and AP_FP_Greater(WS[I],MachineEpsilon*10*Sqrt(NC)*WS[1]) then
        begin
            Tmp[I] := Tmp[I]/WS[I];
        end
        else
        begin
            Tmp[I] := 0;
        end;
        Inc(I);
    end;
    I:=1;
    while I<=ReducedSize do
    begin
        Tmp2[I] := 0;
        Inc(I);
    end;
    I:=1;
    while I<=ReducedSize do
    begin
        V := Tmp[I];
        APVAdd(@Tmp2[0], 1, ReducedSize, @VT[I][0], 1, ReducedSize, V);
        Inc(I);
    end;
    
    //
    // Solution is in the tmp2.
    // Transform it from t to x.
    //
    SetLength(CTbl, M+2+1);
    I:=1;
    while I<=M+1 do
    begin
        V := APVDotProduct(@C[I][0], 1, ReducedSize, @Tmp2[0], 1, ReducedSize);
        CTbl[I-1] := V+D[I];
        Inc(I);
    end;
    CTbl[M+1] := A;
    CTbl[M+2] := B;
end;


(*************************************************************************
Calculation of a Chebyshev  polynomial  obtained   during  least  squares
approximaion at the given point.

Input parameters:
    M   -   order of polynomial (parameter of the
            BuildChebyshevLeastSquares function).
    A   -   coefficient table.
            A[0..M] contains coefficients of the i-th Chebyshev polynomial.
            A[M+1] contains left boundary of approximation interval.
            A[M+2] contains right boundary of approximation interval.
    X   -   point to perform calculations in.

The result is the value at the given point.

It should be noted that array A contains coefficients  of  the  Chebyshev
polynomials defined on interval [-1,1].   Argument  is  reduced  to  this
interval before calculating polynomial value.
  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************)
function CalculateChebyshevLeastSquares(const M : AlglibInteger;
     const A : TReal1DArray;
     X : Double):Double;
var
    b1 : Double;
    b2 : Double;
    i : AlglibInteger;
begin
    X := 2*(X-A[M+1])/(A[M+2]-A[M+1])-1;
    b1 := 0;
    b2 := 0;
    i := M;
    repeat
        Result := 2*x*b1-b2+a[i];
        b2 := b1;
        b1 := Result;
        i := i-1;
    until  not (i>=0);
    Result := Result-x*b2;
end;


end.