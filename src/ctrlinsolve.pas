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
unit ctrlinsolve;
interface
uses Math, Sysutils, Ap;

function ComplexSafeSolveTriangular(const A : TComplex2DArray;
     N : AlglibInteger;
     var X : TComplex1DArray;
     IsUpper : Boolean;
     Trans : AlglibInteger;
     Isunit : Boolean;
     var WORKA : TComplex1DArray;
     var WORKX : TComplex1DArray):Boolean;

implementation

(*************************************************************************
Utility subroutine performing the "safe" solution of a  system  of  linear
equations with triangular complex coefficient matrices.

The feature of an algorithm is that it could not cause an  overflow  or  a
division by zero regardless of the matrix used as the input. If an overflow
is possible, an error code is returned.

The algorithm can solve systems of equations with upper/lower triangular
matrices,  with/without unit diagonal, and systems of types A*x=b, A^T*x=b,
A^H*x=b.

Input parameters:
    A       -   system matrix.
                Array whose indexes range within [1..N, 1..N].
    N       -   size of matrix A.
    X       -   right-hand member of a system.
                Array whose index ranges within [1..N].
    IsUpper -   matrix type. If it is True, the system matrix is the upper
                triangular matrix and is located in the corresponding part
                of matrix A.
    Trans   -   problem type.
                If Trans is:
                    * 0, A*x=b
                    * 1, A^T*x=b
                    * 2, A^H*x=b
    Isunit  -   matrix type. If it is True, the system matrix has  a  unit
                diagonal (the elements on the main diagonal are  not  used
                in the calculation process), otherwise the matrix is
                considered to be a general triangular matrix.
    CNORM   -   array which is stored in norms of rows and columns of  the
                matrix. If the array hasn't been filled up during previous
                executions  of  an  algorithm  with the same matrix as the
                input,  it  will  be  filled  up by the subroutine. If the
                array is filled up, the subroutine uses it without filling
                it up again.
    NORMIN  -   flag defining the state of column norms array. If True, the
                array is filled up.
    WORKA   -   working array whose index ranges within [1..N].
    WORKX   -   working array whose index ranges within [1..N].

Output parameters (if the result is True):
    X       -   solution. Array whose index ranges within [1..N].
    CNORM   -   array of column norms whose index ranges within [1..N].

Result:
    True, if the matrix is not singular  and  the  algorithm  finished its
        work correctly without causing an overflow.
    False, if  the  matrix  is  singular  or  the  algorithm was cancelled
        because of an overflow possibility.

Note:
    The disadvantage of an algorithm is that  sometimes  it  overestimates
    an overflow probability. This is not a problem when  solving  ordinary
    systems. If the elements of the matrix used as the input are close  to
    MaxRealNumber, a false overflow detection is possible, but in practice
    such matrices can rarely be found.
    You can find more reliable subroutines in the LAPACK library
    (xLATRS subroutine ).

  -- ALGLIB --
     Copyright 31.03.2006 by Bochkanov Sergey
*************************************************************************)
function ComplexSafeSolveTriangular(const A : TComplex2DArray;
     N : AlglibInteger;
     var X : TComplex1DArray;
     IsUpper : Boolean;
     Trans : AlglibInteger;
     Isunit : Boolean;
     var WORKA : TComplex1DArray;
     var WORKX : TComplex1DArray):Boolean;
var
    I : AlglibInteger;
    L : AlglibInteger;
    J : AlglibInteger;
    DoLSwp : Boolean;
    MA : Double;
    MX : Double;
    V : Double;
    T : Complex;
    R : Complex;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert((Trans>=0) and (Trans<=2), 'ComplexSafeSolveTriangular: incorrect parameters!');
    Result := True;
    
    //
    // Quick return if possible
    //
    if N<=0 then
    begin
        Exit;
    end;
    
    //
    // Main cycle
    //
    L:=1;
    while L<=N do
    begin
        
        //
        // Prepare subtask L
        //
        DoLSwp := False;
        if Trans=0 then
        begin
            if IsUpper then
            begin
                I := N+1-L;
                i1_ := (I) - (1);
                for i_ := 1 to L do
                begin
                    WORKA[i_] := A[I,i_+i1_];
                end;
                i1_ := (I) - (1);
                for i_ := 1 to L do
                begin
                    WORKX[i_] := X[i_+i1_];
                end;
                DoLSwp := True;
            end;
            if  not IsUpper then
            begin
                I := L;
                for i_ := 1 to L do
                begin
                    WORKA[i_] := A[I,i_];
                end;
                for i_ := 1 to L do
                begin
                    WORKX[i_] := X[i_];
                end;
            end;
        end;
        if Trans=1 then
        begin
            if IsUpper then
            begin
                I := L;
                for i_ := 1 to L do
                begin
                    WORKA[i_] := A[i_,I];
                end;
                for i_ := 1 to L do
                begin
                    WORKX[i_] := X[i_];
                end;
            end;
            if  not IsUpper then
            begin
                I := N+1-L;
                i1_ := (I) - (1);
                for i_ := 1 to L do
                begin
                    WORKA[i_] := A[i_+i1_,I];
                end;
                i1_ := (I) - (1);
                for i_ := 1 to L do
                begin
                    WORKX[i_] := X[i_+i1_];
                end;
                DoLSwp := True;
            end;
        end;
        if Trans=2 then
        begin
            if IsUpper then
            begin
                I := L;
                for i_ := 1 to L do
                begin
                    WORKA[i_] := Conj(A[i_,I]);
                end;
                for i_ := 1 to L do
                begin
                    WORKX[i_] := X[i_];
                end;
            end;
            if  not IsUpper then
            begin
                I := N+1-L;
                i1_ := (I) - (1);
                for i_ := 1 to L do
                begin
                    WORKA[i_] := Conj(A[i_+i1_,I]);
                end;
                i1_ := (I) - (1);
                for i_ := 1 to L do
                begin
                    WORKX[i_] := X[i_+i1_];
                end;
                DoLSwp := True;
            end;
        end;
        if DoLSwp then
        begin
            T := WORKX[L];
            WORKX[L] := WORKX[1];
            WORKX[1] := T;
            T := WORKA[L];
            WORKA[L] := WORKA[1];
            WORKA[1] := T;
        end;
        if Isunit then
        begin
            WORKA[L] := C_Complex(1);
        end;
        
        //
        // Test if workA[L]=0
        //
        if C_EqualR(WORKA[L],0) then
        begin
            Result := False;
            Exit;
        end;
        
        //
        // Now we have:
        //
        //  workA[1:L]*workX[1:L] = b[I]
        //
        // with known workA[1:L] and workX[1:L-1]
        // and unknown workX[L]
        //
        T := C_Complex(0);
        if L>=2 then
        begin
            MA := 0;
            J:=1;
            while J<=L-1 do
            begin
                MA := Max(MA, AbsComplex(WORKA[J]));
                Inc(J);
            end;
            MX := 0;
            J:=1;
            while J<=L-1 do
            begin
                MX := Max(MX, AbsComplex(WORKX[J]));
                Inc(J);
            end;
            if AP_FP_Greater(Max(MA, MX),1) then
            begin
                V := MaxRealNumber/Max(MA, MX);
                V := V/(L-1);
                if AP_FP_Less(V,Min(MA, MX)) then
                begin
                    Result := False;
                    Exit;
                end;
            end;
            T := C_Complex(0.0);
            for i_ := 1 to L-1 do
            begin
                T := C_Add(T,C_Mul(WORKA[i_],WORKX[i_]));
            end;
        end;
        
        //
        // Now we have:
        //
        //  workA[L]*workX[L] + T = b[I]
        //
        if AP_FP_Greater_Eq(Max(AbsComplex(T), AbsComplex(X[I])),0.5*MaxRealNumber) then
        begin
            Result := False;
            Exit;
        end;
        R := C_Sub(X[I],T);
        
        //
        // Now we have:
        //
        //  workA[L]*workX[L] = R
        //
        if C_NotEqualR(R,0) then
        begin
            if AP_FP_Greater_Eq(Ln(AbsComplex(R))-Ln(AbsComplex(WORKA[L])),Ln(MaxRealNumber)) then
            begin
                Result := False;
                Exit;
            end;
        end;
        
        //
        // X[I]
        //
        X[I] := C_Div(R,WORKA[L]);
        Inc(L);
    end;
end;


end.