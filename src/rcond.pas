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
unit rcond;
interface
uses Math, Sysutils, Ap, lu, trlinsolve;

function RMatrixRCond1(const A : TReal2DArray; N : AlglibInteger):Double;
function RMatrixLURCond1(const LUDcmp : TReal2DArray;
     N : AlglibInteger):Double;
function RMatrixRCondInf(const A : TReal2DArray; N : AlglibInteger):Double;
function RMatrixLURCondInf(const LUDcmp : TReal2DArray;
     N : AlglibInteger):Double;
function RCond1(A : TReal2DArray; N : AlglibInteger):Double;
function RCond1LU(const LUDcmp : TReal2DArray; N : AlglibInteger):Double;
function RCondInf(A : TReal2DArray; N : AlglibInteger):Double;
function RCondInfLU(const LUDcmp : TReal2DArray; N : AlglibInteger):Double;

implementation

procedure InternalEstimateRCondLU(const LUDcmp : TReal2DArray;
     N : AlglibInteger;
     OneNorm : Boolean;
     IsANormProvided : Boolean;
     ANORM : Double;
     var RC : Double);forward;
procedure InternalEstimateNorm(N : AlglibInteger;
     var V : TReal1DArray;
     var X : TReal1DArray;
     var ISGN : TInteger1DArray;
     var EST : Double;
     var KASE : AlglibInteger);forward;


(*************************************************************************
Estimate of a matrix condition number (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.

Result: 1/LowerBound(cond(A))
*************************************************************************)
function RMatrixRCond1(const A : TReal2DArray; N : AlglibInteger):Double;
var
    I : AlglibInteger;
    A1 : TReal2DArray;
begin
    Assert(N>=1, 'RMatrixRCond1: N<1!');
    SetLength(A1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        APVMove(@A1[I][0], 1, N, @A[I-1][0], 0, N-1);
        Inc(I);
    end;
    Result := RCond1(A1, N);
end;


(*************************************************************************
Estimate of the condition number of a matrix given by its LU decomposition (1-norm)

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    LUDcmp      -   LU decomposition of a matrix in compact form. Output of
                    the RMatrixLU subroutine.
    N           -   size of matrix A.

Result: 1/LowerBound(cond(A))
*************************************************************************)
function RMatrixLURCond1(const LUDcmp : TReal2DArray;
     N : AlglibInteger):Double;
var
    I : AlglibInteger;
    A1 : TReal2DArray;
begin
    Assert(N>=1, 'RMatrixLURCond1: N<1!');
    SetLength(A1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        APVMove(@A1[I][0], 1, N, @LUDcmp[I-1][0], 0, N-1);
        Inc(I);
    end;
    Result := RCond1LU(A1, N);
end;


(*************************************************************************
Estimate of a matrix condition number (infinity-norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N   -   size of matrix A.

Result: 1/LowerBound(cond(A))
*************************************************************************)
function RMatrixRCondInf(const A : TReal2DArray; N : AlglibInteger):Double;
var
    I : AlglibInteger;
    A1 : TReal2DArray;
begin
    Assert(N>=1, 'RMatrixRCondInf: N<1!');
    SetLength(A1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        APVMove(@A1[I][0], 1, N, @A[I-1][0], 0, N-1);
        Inc(I);
    end;
    Result := RCondInf(A1, N);
end;


(*************************************************************************
Estimate of the condition number of a matrix given by its LU decomposition
(infinity norm).

The algorithm calculates a lower bound of the condition number. In this case,
the algorithm does not return a lower bound of the condition number, but an
inverse number (to avoid an overflow in case of a singular matrix).

Input parameters:
    LUDcmp  -   LU decomposition of a matrix in compact form. Output of
                the RMatrixLU subroutine.
    N       -   size of matrix A.

Result: 1/LowerBound(cond(A))
*************************************************************************)
function RMatrixLURCondInf(const LUDcmp : TReal2DArray;
     N : AlglibInteger):Double;
var
    I : AlglibInteger;
    A1 : TReal2DArray;
begin
    Assert(N>=1, 'RMatrixLURCondInf: N<1!');
    SetLength(A1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        APVMove(@A1[I][0], 1, N, @LUDcmp[I-1][0], 0, N-1);
        Inc(I);
    end;
    Result := RCondInfLU(A1, N);
end;


function RCond1(A : TReal2DArray; N : AlglibInteger):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Nrm : Double;
    Pivots : TInteger1DArray;
begin
    A := DynamicArrayCopy(A);
    Nrm := 0;
    J:=1;
    while J<=N do
    begin
        V := 0;
        I:=1;
        while I<=N do
        begin
            V := V+AbsReal(A[I,J]);
            Inc(I);
        end;
        Nrm := Max(Nrm, V);
        Inc(J);
    end;
    LUDecomposition(A, N, N, Pivots);
    InternalEstimateRCondLU(A, N, True, True, Nrm, V);
    Result := V;
end;


function RCond1LU(const LUDcmp : TReal2DArray; N : AlglibInteger):Double;
var
    V : Double;
begin
    InternalEstimateRCondLU(LUDcmp, N, True, False, 0, V);
    Result := V;
end;


function RCondInf(A : TReal2DArray; N : AlglibInteger):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    Nrm : Double;
    Pivots : TInteger1DArray;
begin
    A := DynamicArrayCopy(A);
    Nrm := 0;
    I:=1;
    while I<=N do
    begin
        V := 0;
        J:=1;
        while J<=N do
        begin
            V := V+AbsReal(A[I,J]);
            Inc(J);
        end;
        Nrm := Max(Nrm, V);
        Inc(I);
    end;
    LUDecomposition(A, N, N, Pivots);
    InternalEstimateRCondLU(A, N, False, True, Nrm, V);
    Result := V;
end;


function RCondInfLU(const LUDcmp : TReal2DArray; N : AlglibInteger):Double;
var
    V : Double;
begin
    InternalEstimateRCondLU(LUDcmp, N, False, False, 0, V);
    Result := V;
end;


procedure InternalEstimateRCondLU(const LUDcmp : TReal2DArray;
     N : AlglibInteger;
     OneNorm : Boolean;
     IsANormProvided : Boolean;
     ANORM : Double;
     var RC : Double);
var
    WORK0 : TReal1DArray;
    WORK1 : TReal1DArray;
    WORK2 : TReal1DArray;
    WORK3 : TReal1DArray;
    IWORK : TInteger1DArray;
    V : Double;
    NORMIN : Boolean;
    I : AlglibInteger;
    IM1 : AlglibInteger;
    IP1 : AlglibInteger;
    IX : AlglibInteger;
    KASE : AlglibInteger;
    KASE1 : AlglibInteger;
    AINVNM : Double;
    ASCALE : Double;
    SL : Double;
    SMLNUM : Double;
    SU : Double;
    MUpper : Boolean;
    MTrans : Boolean;
    Munit : Boolean;
    i_ : AlglibInteger;
begin
    
    //
    // Quick return if possible
    //
    if N=0 then
    begin
        RC := 1;
        Exit;
    end;
    
    //
    // init
    //
    if OneNorm then
    begin
        KASE1 := 1;
    end
    else
    begin
        KASE1 := 2;
    end;
    MUpper := True;
    MTrans := True;
    Munit := True;
    SetLength(WORK0, N+1);
    SetLength(WORK1, N+1);
    SetLength(WORK2, N+1);
    SetLength(WORK3, N+1);
    SetLength(IWORK, N+1);
    
    //
    // Estimate the norm of A.
    //
    if  not IsANormProvided then
    begin
        KASE := 0;
        ANORM := 0;
        while True do
        begin
            InternalEstimateNorm(N, WORK1, WORK0, IWORK, ANORM, KASE);
            if KASE=0 then
            begin
                Break;
            end;
            if KASE=KASE1 then
            begin
                
                //
                // Multiply by U
                //
                I:=1;
                while I<=N do
                begin
                    V := APVDotProduct(@LUDcmp[I][0], I, N, @WORK0[0], I, N);
                    WORK0[I] := V;
                    Inc(I);
                end;
                
                //
                // Multiply by L
                //
                I:=N;
                while I>=1 do
                begin
                    IM1 := I-1;
                    if I>1 then
                    begin
                        V := APVDotProduct(@LUDcmp[I][0], 1, IM1, @WORK0[0], 1, IM1);
                    end
                    else
                    begin
                        V := 0;
                    end;
                    WORK0[I] := WORK0[I]+V;
                    Dec(I);
                end;
            end
            else
            begin
                
                //
                // Multiply by L'
                //
                I:=1;
                while I<=N do
                begin
                    IP1 := I+1;
                    V := 0.0;
                    for i_ := IP1 to N do
                    begin
                        V := V + LUDcmp[i_,I]*WORK0[i_];
                    end;
                    WORK0[I] := WORK0[I]+V;
                    Inc(I);
                end;
                
                //
                // Multiply by U'
                //
                I:=N;
                while I>=1 do
                begin
                    V := 0.0;
                    for i_ := 1 to I do
                    begin
                        V := V + LUDcmp[i_,I]*WORK0[i_];
                    end;
                    WORK0[I] := V;
                    Dec(I);
                end;
            end;
        end;
    end;
    
    //
    // Quick return if possible
    //
    RC := 0;
    if AP_FP_Eq(ANORM,0) then
    begin
        Exit;
    end;
    
    //
    // Estimate the norm of inv(A).
    //
    SMLNUM := MinRealNumber;
    AINVNM := 0;
    NORMIN := False;
    KASE := 0;
    while True do
    begin
        InternalEstimateNorm(N, WORK1, WORK0, IWORK, AINVNM, KASE);
        if KASE=0 then
        begin
            Break;
        end;
        if KASE=KASE1 then
        begin
            
            //
            // Multiply by inv(L).
            //
            SafeSolveTriangular(LUDcmp, N, WORK0, SL,  not MUpper,  not MTrans, Munit, NORMIN, WORK2);
            
            //
            // Multiply by inv(U).
            //
            SafeSolveTriangular(LUDcmp, N, WORK0, SU, MUpper,  not MTrans,  not Munit, NORMIN, WORK3);
        end
        else
        begin
            
            //
            // Multiply by inv(U').
            //
            SafeSolveTriangular(LUDcmp, N, WORK0, SU, MUpper, MTrans,  not Munit, NORMIN, WORK3);
            
            //
            // Multiply by inv(L').
            //
            SafeSolveTriangular(LUDcmp, N, WORK0, SL,  not MUpper, MTrans, Munit, NORMIN, WORK2);
        end;
        
        //
        // Divide X by 1/(SL*SU) if doing so will not cause overflow.
        //
        ASCALE := SL*SU;
        NORMIN := True;
        if AP_FP_Neq(ASCALE,1) then
        begin
            IX := 1;
            I:=2;
            while I<=N do
            begin
                if AP_FP_Greater(AbsReal(WORK0[I]),AbsReal(WORK0[IX])) then
                begin
                    IX := I;
                end;
                Inc(I);
            end;
            if AP_FP_Less(ASCALE,ABSReal(WORK0[IX])*SMLNUM) or AP_FP_Eq(ASCALE,0) then
            begin
                Exit;
            end;
            I:=1;
            while I<=N do
            begin
                WORK0[I] := WORK0[I]/ASCALE;
                Inc(I);
            end;
        end;
    end;
    
    //
    // Compute the estimate of the reciprocal condition number.
    //
    if AP_FP_Neq(AINVNM,0) then
    begin
        RC := 1/AINVNM;
        RC := RC/ANORM;
    end;
end;


procedure InternalEstimateNorm(N : AlglibInteger;
     var V : TReal1DArray;
     var X : TReal1DArray;
     var ISGN : TInteger1DArray;
     var EST : Double;
     var KASE : AlglibInteger);
var
    ITMAX : AlglibInteger;
    I : AlglibInteger;
    T : Double;
    Flg : Boolean;
    PosITER : AlglibInteger;
    PosJ : AlglibInteger;
    PosJLAST : AlglibInteger;
    PosJUMP : AlglibInteger;
    PosALTSGN : AlglibInteger;
    PosESTOLD : AlglibInteger;
    PosTEMP : AlglibInteger;
begin
    ITMAX := 5;
    PosALTSGN := N+1;
    PosESTOLD := N+2;
    PosTEMP := N+3;
    PosITER := N+1;
    PosJ := N+2;
    PosJLAST := N+3;
    PosJUMP := N+4;
    if KASE=0 then
    begin
        SetLength(V, N+3+1);
        SetLength(X, N+1);
        SetLength(ISGN, N+4+1);
        T := AP_Double(1)/N;
        I:=1;
        while I<=N do
        begin
            X[I] := T;
            Inc(I);
        end;
        KASE := 1;
        ISGN[PosJUMP] := 1;
        Exit;
    end;
    
    //
    //     ................ ENTRY   (JUMP = 1)
    //     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
    //
    if ISGN[PosJUMP]=1 then
    begin
        if N=1 then
        begin
            V[1] := X[1];
            EST := ABSReal(V[1]);
            KASE := 0;
            Exit;
        end;
        EST := 0;
        I:=1;
        while I<=N do
        begin
            EST := EST+AbsReal(X[I]);
            Inc(I);
        end;
        I:=1;
        while I<=N do
        begin
            if AP_FP_Greater_Eq(X[I],0) then
            begin
                X[I] := 1;
            end
            else
            begin
                X[I] := -1;
            end;
            ISGN[I] := Sign(X[I]);
            Inc(I);
        end;
        KASE := 2;
        ISGN[PosJUMP] := 2;
        Exit;
    end;
    
    //
    //     ................ ENTRY   (JUMP = 2)
    //     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
    //
    if ISGN[PosJUMP]=2 then
    begin
        ISGN[PosJ] := 1;
        I:=2;
        while I<=N do
        begin
            if AP_FP_Greater(AbsReal(X[I]),AbsReal(X[ISGN[PosJ]])) then
            begin
                ISGN[PosJ] := I;
            end;
            Inc(I);
        end;
        ISGN[PosITER] := 2;
        
        //
        // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
        //
        I:=1;
        while I<=N do
        begin
            X[I] := 0;
            Inc(I);
        end;
        X[ISGN[PosJ]] := 1;
        KASE := 1;
        ISGN[PosJUMP] := 3;
        Exit;
    end;
    
    //
    //     ................ ENTRY   (JUMP = 3)
    //     X HAS BEEN OVERWRITTEN BY A*X.
    //
    if ISGN[PosJUMP]=3 then
    begin
        APVMove(@V[0], 1, N, @X[0], 1, N);
        V[PosESTOLD] := EST;
        EST := 0;
        I:=1;
        while I<=N do
        begin
            EST := EST+AbsReal(V[I]);
            Inc(I);
        end;
        Flg := False;
        I:=1;
        while I<=N do
        begin
            if AP_FP_Greater_Eq(X[I],0) and (ISGN[I]<0) or AP_FP_Less(X[I],0) and (ISGN[I]>=0) then
            begin
                Flg := True;
            end;
            Inc(I);
        end;
        
        //
        // REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
        // OR MAY BE CYCLING.
        //
        if  not Flg or AP_FP_Less_Eq(EST,V[PosESTOLD]) then
        begin
            V[PosALTSGN] := 1;
            I:=1;
            while I<=N do
            begin
                X[I] := V[PosALTSGN]*(1+AP_Double((I-1))/(N-1));
                V[PosALTSGN] := -V[PosALTSGN];
                Inc(I);
            end;
            KASE := 1;
            ISGN[PosJUMP] := 5;
            Exit;
        end;
        I:=1;
        while I<=N do
        begin
            if AP_FP_Greater_Eq(X[I],0) then
            begin
                X[I] := 1;
                ISGN[I] := 1;
            end
            else
            begin
                X[I] := -1;
                ISGN[I] := -1;
            end;
            Inc(I);
        end;
        KASE := 2;
        ISGN[PosJUMP] := 4;
        Exit;
    end;
    
    //
    //     ................ ENTRY   (JUMP = 4)
    //     X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
    //
    if ISGN[PosJUMP]=4 then
    begin
        ISGN[PosJLAST] := ISGN[PosJ];
        ISGN[PosJ] := 1;
        I:=2;
        while I<=N do
        begin
            if AP_FP_Greater(AbsReal(X[I]),AbsReal(X[ISGN[PosJ]])) then
            begin
                ISGN[PosJ] := I;
            end;
            Inc(I);
        end;
        if AP_FP_Neq(X[ISGN[PosJLAST]],ABSReal(X[ISGN[PosJ]])) and (ISGN[PosITER]<ITMAX) then
        begin
            ISGN[PosITER] := ISGN[PosITER]+1;
            I:=1;
            while I<=N do
            begin
                X[I] := 0;
                Inc(I);
            end;
            X[ISGN[PosJ]] := 1;
            KASE := 1;
            ISGN[PosJUMP] := 3;
            Exit;
        end;
        
        //
        // ITERATION COMPLETE.  FINAL STAGE.
        //
        V[PosALTSGN] := 1;
        I:=1;
        while I<=N do
        begin
            X[I] := V[PosALTSGN]*(1+AP_Double((I-1))/(N-1));
            V[PosALTSGN] := -V[PosALTSGN];
            Inc(I);
        end;
        KASE := 1;
        ISGN[PosJUMP] := 5;
        Exit;
    end;
    
    //
    //     ................ ENTRY   (JUMP = 5)
    //     X HAS BEEN OVERWRITTEN BY A*X.
    //
    if ISGN[PosJUMP]=5 then
    begin
        V[PosTEMP] := 0;
        I:=1;
        while I<=N do
        begin
            V[PosTEMP] := V[PosTEMP]+AbsReal(X[I]);
            Inc(I);
        end;
        V[PosTEMP] := 2*V[PosTEMP]/(3*N);
        if AP_FP_Greater(V[PosTEMP],EST) then
        begin
            APVMove(@V[0], 1, N, @X[0], 1, N);
            EST := V[PosTEMP];
        end;
        KASE := 0;
        Exit;
    end;
end;


end.