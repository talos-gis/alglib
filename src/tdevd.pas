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
unit tdevd;
interface
uses Math, Sysutils, Ap, blas, rotations;

function SMatrixTDEVD(var D : TReal1DArray;
     E : TReal1DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     var Z : TReal2DArray):Boolean;
function TridiagonalEVD(var D : TReal1DArray;
     E : TReal1DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     var Z : TReal2DArray):Boolean;

implementation

procedure TdEVDE2(const A : Double;
     const B : Double;
     const C : Double;
     var RT1 : Double;
     var RT2 : Double);forward;
procedure TdEVDEV2(const A : Double;
     const B : Double;
     const C : Double;
     var RT1 : Double;
     var RT2 : Double;
     var CS1 : Double;
     var SN1 : Double);forward;
function TdEVDPythag(A : Double; B : Double):Double;forward;
function TdEVDExtSign(a : Double; b : Double):Double;forward;


(*************************************************************************
Finding the eigenvalues and eigenvectors of a tridiagonal symmetric matrix

The algorithm finds the eigen pairs of a tridiagonal symmetric matrix by
using an QL/QR algorithm with implicit shifts.

Input parameters:
    D       -   the main diagonal of a tridiagonal matrix.
                Array whose index ranges within [0..N-1].
    E       -   the secondary diagonal of a tridiagonal matrix.
                Array whose index ranges within [0..N-2].
    N       -   size of matrix A.
    ZNeeded -   flag controlling whether the eigenvectors are needed or not.
                If ZNeeded is equal to:
                 * 0, the eigenvectors are not needed;
                 * 1, the eigenvectors of a tridiagonal matrix
                   are multiplied by the square matrix Z. It is used if the
                   tridiagonal matrix is obtained by the similarity
                   transformation of a symmetric matrix;
                 * 2, the eigenvectors of a tridiagonal matrix replace the
                   square matrix Z;
                 * 3, matrix Z contains the first row of the eigenvectors
                   matrix.
    Z       -   if ZNeeded=1, Z contains the square matrix by which the
                eigenvectors are multiplied.
                Array whose indexes range within [0..N-1, 0..N-1].

Output parameters:
    D       -   eigenvalues in ascending order.
                Array whose index ranges within [0..N-1].
    Z       -   if ZNeeded is equal to:
                 * 0, Z hasn’t changed;
                 * 1, Z contains the product of a given matrix (from the left)
                   and the eigenvectors matrix (from the right);
                 * 2, Z contains the eigenvectors.
                 * 3, Z contains the first row of the eigenvectors matrix.
                If ZNeeded<3, Z is the array whose indexes range within [0..N-1, 0..N-1].
                In that case, the eigenvectors are stored in the matrix columns.
                If ZNeeded=3, Z is the array whose indexes range within [0..0, 0..N-1].

Result:
    True, if the algorithm has converged.
    False, if the algorithm hasn't converged.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************)
function SMatrixTDEVD(var D : TReal1DArray;
     E : TReal1DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     var Z : TReal2DArray):Boolean;
var
    D1 : TReal1DArray;
    E1 : TReal1DArray;
    Z1 : TReal2DArray;
    I : AlglibInteger;
begin
    E := DynamicArrayCopy(E);
    
    //
    // Prepare 1-based task
    //
    SetLength(D1, N+1);
    SetLength(E1, N+1);
    APVMove(@D1[0], 1, N, @D[0], 0, N-1);
    if N>1 then
    begin
        APVMove(@E1[0], 1, N-1, @E[0], 0, N-2);
    end;
    if ZNeeded=1 then
    begin
        SetLength(Z1, N+1, N+1);
        I:=1;
        while I<=N do
        begin
            APVMove(@Z1[I][0], 1, N, @Z[I-1][0], 0, N-1);
            Inc(I);
        end;
    end;
    
    //
    // Solve 1-based task
    //
    Result := TridiagonalEVD(D1, E1, N, ZNeeded, Z1);
    if  not Result then
    begin
        Exit;
    end;
    
    //
    // Convert back to 0-based result
    //
    APVMove(@D[0], 0, N-1, @D1[0], 1, N);
    if ZNeeded<>0 then
    begin
        if ZNeeded=1 then
        begin
            I:=1;
            while I<=N do
            begin
                APVMove(@Z[I-1][0], 0, N-1, @Z1[I][0], 1, N);
                Inc(I);
            end;
            Exit;
        end;
        if ZNeeded=2 then
        begin
            SetLength(Z, N-1+1, N-1+1);
            I:=1;
            while I<=N do
            begin
                APVMove(@Z[I-1][0], 0, N-1, @Z1[I][0], 1, N);
                Inc(I);
            end;
            Exit;
        end;
        if ZNeeded=3 then
        begin
            SetLength(Z, 0+1, N-1+1);
            APVMove(@Z[0][0], 0, N-1, @Z1[1][0], 1, N);
            Exit;
        end;
        Assert(False, 'SMatrixTDEVD: Incorrect ZNeeded!');
    end;
end;


function TridiagonalEVD(var D : TReal1DArray;
     E : TReal1DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     var Z : TReal2DArray):Boolean;
var
    MAXIT : AlglibInteger;
    I : AlglibInteger;
    II : AlglibInteger;
    ISCALE : AlglibInteger;
    J : AlglibInteger;
    JTOT : AlglibInteger;
    K : AlglibInteger;
    T : AlglibInteger;
    L : AlglibInteger;
    L1 : AlglibInteger;
    LEND : AlglibInteger;
    LENDM1 : AlglibInteger;
    LENDP1 : AlglibInteger;
    LENDSV : AlglibInteger;
    LM1 : AlglibInteger;
    LSV : AlglibInteger;
    M : AlglibInteger;
    MM : AlglibInteger;
    MM1 : AlglibInteger;
    NM1 : AlglibInteger;
    NMAXIT : AlglibInteger;
    TmpInt : AlglibInteger;
    ANORM : Double;
    B : Double;
    C : Double;
    EPS : Double;
    EPS2 : Double;
    F : Double;
    G : Double;
    P : Double;
    R : Double;
    RT1 : Double;
    RT2 : Double;
    S : Double;
    SAFMAX : Double;
    SAFMIN : Double;
    SSFMAX : Double;
    SSFMIN : Double;
    TST : Double;
    Tmp : Double;
    WORK1 : TReal1DArray;
    WORK2 : TReal1DArray;
    WORKC : TReal1DArray;
    WORKS : TReal1DArray;
    WTEMP : TReal1DArray;
    GotoFlag : Boolean;
    ZRows : AlglibInteger;
    WasTranspose : Boolean;
    i_ : AlglibInteger;
begin
    E := DynamicArrayCopy(E);
    Assert((ZNeeded>=0) and (ZNeeded<=3), 'TridiagonalEVD: Incorrent ZNeeded');
    
    //
    // Quick return if possible
    //
    if (ZNeeded<0) or (ZNeeded>3) then
    begin
        Result := False;
        Exit;
    end;
    Result := True;
    if N=0 then
    begin
        Exit;
    end;
    if N=1 then
    begin
        if (ZNeeded=2) or (ZNeeded=3) then
        begin
            SetLength(Z, 1+1, 1+1);
            Z[1,1] := 1;
        end;
        Exit;
    end;
    MAXIT := 30;
    
    //
    // Initialize arrays
    //
    SetLength(WTEMP, N+1);
    SetLength(WORK1, N-1+1);
    SetLength(WORK2, N-1+1);
    SetLength(WORKC, N+1);
    SetLength(WORKS, N+1);
    
    //
    // Determine the unit roundoff and over/underflow thresholds.
    //
    EPS := MachineEpsilon;
    EPS2 := AP_Sqr(EPS);
    SAFMIN := MinRealNumber;
    SAFMAX := MaxRealNumber;
    SSFMAX := Sqrt(SAFMAX)/3;
    SSFMIN := Sqrt(SAFMIN)/EPS2;
    
    //
    // Prepare Z
    //
    // Here we are using transposition to get rid of column operations
    //
    //
    WasTranspose := False;
    if ZNeeded=0 then
    begin
        ZRows := 0;
    end;
    if ZNeeded=1 then
    begin
        ZRows := N;
    end;
    if ZNeeded=2 then
    begin
        ZRows := N;
    end;
    if ZNeeded=3 then
    begin
        ZRows := 1;
    end;
    if ZNeeded=1 then
    begin
        WasTranspose := True;
        InplaceTranspose(Z, 1, N, 1, N, WTEMP);
    end;
    if ZNeeded=2 then
    begin
        WasTranspose := True;
        SetLength(Z, N+1, N+1);
        I:=1;
        while I<=N do
        begin
            J:=1;
            while J<=N do
            begin
                if I=J then
                begin
                    Z[I,J] := 1;
                end
                else
                begin
                    Z[I,J] := 0;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    if ZNeeded=3 then
    begin
        WasTranspose := False;
        SetLength(Z, 1+1, N+1);
        J:=1;
        while J<=N do
        begin
            if J=1 then
            begin
                Z[1,J] := 1;
            end
            else
            begin
                Z[1,J] := 0;
            end;
            Inc(J);
        end;
    end;
    NMAXIT := N*MAXIT;
    JTOT := 0;
    
    //
    // Determine where the matrix splits and choose QL or QR iteration
    // for each block, according to whether top or bottom diagonal
    // element is smaller.
    //
    L1 := 1;
    NM1 := N-1;
    while True do
    begin
        if L1>N then
        begin
            Break;
        end;
        if L1>1 then
        begin
            E[L1-1] := 0;
        end;
        GotoFlag := False;
        if L1<=NM1 then
        begin
            M:=L1;
            while M<=NM1 do
            begin
                TST := ABSReal(E[M]);
                if AP_FP_Eq(TST,0) then
                begin
                    GotoFlag := True;
                    Break;
                end;
                if AP_FP_Less_Eq(TST,Sqrt(AbsReal(D[M]))*Sqrt(AbsReal(D[M+1]))*EPS) then
                begin
                    E[M] := 0;
                    GotoFlag := True;
                    Break;
                end;
                Inc(M);
            end;
        end;
        if  not GotoFlag then
        begin
            M := N;
        end;
        
        //
        // label 30:
        //
        L := L1;
        LSV := L;
        LEND := M;
        LENDSV := LEND;
        L1 := M+1;
        if LEND=L then
        begin
            Continue;
        end;
        
        //
        // Scale submatrix in rows and columns L to LEND
        //
        if L=LEND then
        begin
            ANORM := AbsReal(D[L]);
        end
        else
        begin
            ANORM := Max(AbsReal(D[L])+AbsReal(E[L]), AbsReal(E[LEND-1])+AbsReal(D[LEND]));
            I:=L+1;
            while I<=LEND-1 do
            begin
                ANORM := Max(ANORM, AbsReal(D[I])+AbsReal(E[I])+AbsReal(E[I-1]));
                Inc(I);
            end;
        end;
        ISCALE := 0;
        if AP_FP_Eq(ANORM,0) then
        begin
            Continue;
        end;
        if AP_FP_Greater(ANORM,SSFMAX) then
        begin
            ISCALE := 1;
            Tmp := SSFMAX/ANORM;
            TmpInt := LEND-1;
            APVMul(@D[0], L, LEND, Tmp);
            APVMul(@E[0], L, TmpInt, Tmp);
        end;
        if AP_FP_Less(ANORM,SSFMIN) then
        begin
            ISCALE := 2;
            Tmp := SSFMIN/ANORM;
            TmpInt := LEND-1;
            APVMul(@D[0], L, LEND, Tmp);
            APVMul(@E[0], L, TmpInt, Tmp);
        end;
        
        //
        // Choose between QL and QR iteration
        //
        if AP_FP_Less(AbsReal(D[LEND]),AbsReal(D[L])) then
        begin
            LEND := LSV;
            L := LENDSV;
        end;
        if LEND>L then
        begin
            
            //
            // QL Iteration
            //
            // Look for small subdiagonal element.
            //
            while True do
            begin
                GotoFlag := False;
                if L<>LEND then
                begin
                    LENDM1 := LEND-1;
                    M:=L;
                    while M<=LENDM1 do
                    begin
                        TST := AP_Sqr(AbsReal(E[M]));
                        if AP_FP_Less_Eq(TST,EPS2*AbsReal(D[M])*AbsReal(D[M+1])+SAFMIN) then
                        begin
                            GotoFlag := True;
                            Break;
                        end;
                        Inc(M);
                    end;
                end;
                if  not GotoFlag then
                begin
                    M := LEND;
                end;
                if M<LEND then
                begin
                    E[M] := 0;
                end;
                P := D[L];
                if M<>L then
                begin
                    
                    //
                    // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
                    // to compute its eigensystem.
                    //
                    if M=L+1 then
                    begin
                        if ZNeeded>0 then
                        begin
                            TdEVDEV2(D[L], E[L], D[L+1], RT1, RT2, C, S);
                            WORK1[L] := C;
                            WORK2[L] := S;
                            WORKC[1] := WORK1[L];
                            WORKS[1] := WORK2[L];
                            if  not WasTranspose then
                            begin
                                ApplyRotationsFromTheRight(False, 1, ZRows, L, L+1, WORKC, WORKS, Z, WTEMP);
                            end
                            else
                            begin
                                ApplyRotationsFromTheLeft(False, L, L+1, 1, ZRows, WORKC, WORKS, Z, WTEMP);
                            end;
                        end
                        else
                        begin
                            TdEVDE2(D[L], E[L], D[L+1], RT1, RT2);
                        end;
                        D[L] := RT1;
                        D[L+1] := RT2;
                        E[L] := 0;
                        L := L+2;
                        if L<=LEND then
                        begin
                            Continue;
                        end;
                        
                        //
                        // GOTO 140
                        //
                        Break;
                    end;
                    if JTOT=NMAXIT then
                    begin
                        
                        //
                        // GOTO 140
                        //
                        Break;
                    end;
                    JTOT := JTOT+1;
                    
                    //
                    // Form shift.
                    //
                    G := (D[L+1]-P)/(2*E[L]);
                    R := TdEVDPythag(G, 1);
                    G := D[M]-P+E[L]/(G+TdEVDExtSign(R, G));
                    S := 1;
                    C := 1;
                    P := 0;
                    
                    //
                    // Inner loop
                    //
                    MM1 := M-1;
                    I:=MM1;
                    while I>=L do
                    begin
                        F := S*E[I];
                        B := C*E[I];
                        GenerateRotation(G, F, C, S, R);
                        if I<>M-1 then
                        begin
                            E[I+1] := R;
                        end;
                        G := D[I+1]-P;
                        R := (D[I]-G)*S+2*C*B;
                        P := S*R;
                        D[I+1] := G+P;
                        G := C*R-B;
                        
                        //
                        // If eigenvectors are desired, then save rotations.
                        //
                        if ZNeeded>0 then
                        begin
                            WORK1[I] := C;
                            WORK2[I] := -S;
                        end;
                        Dec(I);
                    end;
                    
                    //
                    // If eigenvectors are desired, then apply saved rotations.
                    //
                    if ZNeeded>0 then
                    begin
                        I:=L;
                        while I<=M-1 do
                        begin
                            WORKC[I-L+1] := WORK1[I];
                            WORKS[I-L+1] := WORK2[I];
                            Inc(I);
                        end;
                        if  not WasTranspose then
                        begin
                            ApplyRotationsFromTheRight(False, 1, ZRows, L, M, WORKC, WORKS, Z, WTEMP);
                        end
                        else
                        begin
                            ApplyRotationsFromTheLeft(False, L, M, 1, ZRows, WORKC, WORKS, Z, WTEMP);
                        end;
                    end;
                    D[L] := D[L]-P;
                    E[L] := G;
                    Continue;
                end;
                
                //
                // Eigenvalue found.
                //
                D[L] := P;
                L := L+1;
                if L<=LEND then
                begin
                    Continue;
                end;
                Break;
            end;
        end
        else
        begin
            
            //
            // QR Iteration
            //
            // Look for small superdiagonal element.
            //
            while True do
            begin
                GotoFlag := False;
                if L<>LEND then
                begin
                    LENDP1 := LEND+1;
                    M:=L;
                    while M>=LENDP1 do
                    begin
                        TST := AP_Sqr(ABSReal(E[M-1]));
                        if AP_FP_Less_Eq(TST,EPS2*ABSReal(D[M])*ABSReal(D[M-1])+SAFMIN) then
                        begin
                            GotoFlag := True;
                            Break;
                        end;
                        Dec(M);
                    end;
                end;
                if  not GotoFlag then
                begin
                    M := LEND;
                end;
                if M>LEND then
                begin
                    E[M-1] := 0;
                end;
                P := D[L];
                if M<>L then
                begin
                    
                    //
                    // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
                    // to compute its eigensystem.
                    //
                    if M=L-1 then
                    begin
                        if ZNeeded>0 then
                        begin
                            TdEVDEV2(D[L-1], E[L-1], D[L], RT1, RT2, C, S);
                            WORK1[M] := C;
                            WORK2[M] := S;
                            WORKC[1] := C;
                            WORKS[1] := S;
                            if  not WasTranspose then
                            begin
                                ApplyRotationsFromTheRight(True, 1, ZRows, L-1, L, WORKC, WORKS, Z, WTEMP);
                            end
                            else
                            begin
                                ApplyRotationsFromTheLeft(True, L-1, L, 1, ZRows, WORKC, WORKS, Z, WTEMP);
                            end;
                        end
                        else
                        begin
                            TdEVDE2(D[L-1], E[L-1], D[L], RT1, RT2);
                        end;
                        D[L-1] := RT1;
                        D[L] := RT2;
                        E[L-1] := 0;
                        L := L-2;
                        if L>=LEND then
                        begin
                            Continue;
                        end;
                        Break;
                    end;
                    if JTOT=NMAXIT then
                    begin
                        Break;
                    end;
                    JTOT := JTOT+1;
                    
                    //
                    // Form shift.
                    //
                    G := (D[L-1]-P)/(2*E[L-1]);
                    R := TdEVDPythag(G, 1);
                    G := D[M]-P+E[L-1]/(G+TdEVDExtSign(R, G));
                    S := 1;
                    C := 1;
                    P := 0;
                    
                    //
                    // Inner loop
                    //
                    LM1 := L-1;
                    I:=M;
                    while I<=LM1 do
                    begin
                        F := S*E[I];
                        B := C*E[I];
                        GenerateRotation(G, F, C, S, R);
                        if I<>M then
                        begin
                            E[I-1] := R;
                        end;
                        G := D[I]-P;
                        R := (D[I+1]-G)*S+2*C*B;
                        P := S*R;
                        D[I] := G+P;
                        G := C*R-B;
                        
                        //
                        // If eigenvectors are desired, then save rotations.
                        //
                        if ZNeeded>0 then
                        begin
                            WORK1[I] := C;
                            WORK2[I] := S;
                        end;
                        Inc(I);
                    end;
                    
                    //
                    // If eigenvectors are desired, then apply saved rotations.
                    //
                    if ZNeeded>0 then
                    begin
                        MM := L-M+1;
                        I:=M;
                        while I<=L-1 do
                        begin
                            WORKC[I-M+1] := WORK1[I];
                            WORKS[I-M+1] := WORK2[I];
                            Inc(I);
                        end;
                        if  not WasTranspose then
                        begin
                            ApplyRotationsFromTheRight(True, 1, ZRows, M, L, WORKC, WORKS, Z, WTEMP);
                        end
                        else
                        begin
                            ApplyRotationsFromTheLeft(True, M, L, 1, ZRows, WORKC, WORKS, Z, WTEMP);
                        end;
                    end;
                    D[L] := D[L]-P;
                    E[LM1] := G;
                    Continue;
                end;
                
                //
                // Eigenvalue found.
                //
                D[L] := P;
                L := L-1;
                if L>=LEND then
                begin
                    Continue;
                end;
                Break;
            end;
        end;
        
        //
        // Undo scaling if necessary
        //
        if ISCALE=1 then
        begin
            Tmp := ANORM/SSFMAX;
            TmpInt := LENDSV-1;
            APVMul(@D[0], LSV, LENDSV, Tmp);
            APVMul(@E[0], LSV, TmpInt, Tmp);
        end;
        if ISCALE=2 then
        begin
            Tmp := ANORM/SSFMIN;
            TmpInt := LENDSV-1;
            APVMul(@D[0], LSV, LENDSV, Tmp);
            APVMul(@E[0], LSV, TmpInt, Tmp);
        end;
        
        //
        // Check for no convergence to an eigenvalue after a total
        // of N*MAXIT iterations.
        //
        if JTOT>=NMAXIT then
        begin
            Result := False;
            if WasTranspose then
            begin
                InplaceTranspose(Z, 1, N, 1, N, WTEMP);
            end;
            Exit;
        end;
    end;
    
    //
    // Order eigenvalues and eigenvectors.
    //
    if ZNeeded=0 then
    begin
        
        //
        // Sort
        //
        if N=1 then
        begin
            Exit;
        end;
        if N=2 then
        begin
            if AP_FP_Greater(D[1],D[2]) then
            begin
                Tmp := D[1];
                D[1] := D[2];
                D[2] := Tmp;
            end;
            Exit;
        end;
        i := 2;
        repeat
            t := i;
            while t<>1 do
            begin
                k := t div 2;
                if AP_FP_Greater_Eq(D[k],D[t]) then
                begin
                    t := 1;
                end
                else
                begin
                    Tmp := D[k];
                    D[k] := D[t];
                    D[t] := Tmp;
                    t := k;
                end;
            end;
            i := i+1;
        until  not (i<=n);
        i := n-1;
        repeat
            Tmp := D[i+1];
            D[i+1] := D[1];
            D[+1] := Tmp;
            t := 1;
            while t<>0 do
            begin
                k := 2*t;
                if k>i then
                begin
                    t := 0;
                end
                else
                begin
                    if k<i then
                    begin
                        if AP_FP_Greater(D[k+1],D[k]) then
                        begin
                            k := k+1;
                        end;
                    end;
                    if AP_FP_Greater_Eq(D[t],D[k]) then
                    begin
                        t := 0;
                    end
                    else
                    begin
                        Tmp := D[k];
                        D[k] := D[t];
                        D[t] := Tmp;
                        t := k;
                    end;
                end;
            end;
            i := i-1;
        until  not (i>=1);
    end
    else
    begin
        
        //
        // Use Selection Sort to minimize swaps of eigenvectors
        //
        II:=2;
        while II<=N do
        begin
            I := II-1;
            K := I;
            P := D[I];
            J:=II;
            while J<=N do
            begin
                if AP_FP_Less(D[J],P) then
                begin
                    K := J;
                    P := D[J];
                end;
                Inc(J);
            end;
            if K<>I then
            begin
                D[K] := D[I];
                D[I] := P;
                if WasTranspose then
                begin
                    APVMove(@WTEMP[0], 1, N, @Z[I][0], 1, N);
                    APVMove(@Z[I][0], 1, N, @Z[K][0], 1, N);
                    APVMove(@Z[K][0], 1, N, @WTEMP[0], 1, N);
                end
                else
                begin
                    for i_ := 1 to ZRows do
                    begin
                        WTEMP[i_] := Z[i_,I];
                    end;
                    for i_ := 1 to ZRows do
                    begin
                        Z[i_,I] := Z[i_,K];
                    end;
                    for i_ := 1 to ZRows do
                    begin
                        Z[i_,K] := WTEMP[i_];
                    end;
                end;
            end;
            Inc(II);
        end;
        if WasTranspose then
        begin
            InplaceTranspose(Z, 1, N, 1, N, WTEMP);
        end;
    end;
end;


(*************************************************************************
DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
   [  A   B  ]
   [  B   C  ].
On return, RT1 is the eigenvalue of larger absolute value, and RT2
is the eigenvalue of smaller absolute value.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1992
*************************************************************************)
procedure TdEVDE2(const A : Double;
     const B : Double;
     const C : Double;
     var RT1 : Double;
     var RT2 : Double);
var
    AB : Double;
    ACMN : Double;
    ACMX : Double;
    ADF : Double;
    DF : Double;
    RT : Double;
    SM : Double;
    TB : Double;
begin
    SM := A+C;
    DF := A-C;
    ADF := AbsReal(DF);
    TB := B+B;
    AB := AbsReal(TB);
    if AP_FP_Greater(AbsReal(A),AbsReal(C)) then
    begin
        ACMX := A;
        ACMN := C;
    end
    else
    begin
        ACMX := C;
        ACMN := A;
    end;
    if AP_FP_Greater(ADF,AB) then
    begin
        RT := ADF*Sqrt(1+AP_Sqr(AB/ADF));
    end
    else
    begin
        if AP_FP_Less(ADF,AB) then
        begin
            RT := AB*Sqrt(1+AP_Sqr(ADF/AB));
        end
        else
        begin
            
            //
            // Includes case AB=ADF=0
            //
            RT := AB*Sqrt(2);
        end;
    end;
    if AP_FP_Less(SM,0) then
    begin
        RT1 := 0.5*(SM-RT);
        
        //
        // Order of execution important.
        // To get fully accurate smaller eigenvalue,
        // next line needs to be executed in higher precision.
        //
        RT2 := ACMX/RT1*ACMN-B/RT1*B;
    end
    else
    begin
        if AP_FP_Greater(SM,0) then
        begin
            RT1 := 0.5*(SM+RT);
            
            //
            // Order of execution important.
            // To get fully accurate smaller eigenvalue,
            // next line needs to be executed in higher precision.
            //
            RT2 := ACMX/RT1*ACMN-B/RT1*B;
        end
        else
        begin
            
            //
            // Includes case RT1 = RT2 = 0
            //
            RT1 := 0.5*RT;
            RT2 := -0.5*RT;
        end;
    end;
end;


(*************************************************************************
DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix

   [  A   B  ]
   [  B   C  ].

On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
eigenvector for RT1, giving the decomposition

   [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
   [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].


  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     October 31, 1992
*************************************************************************)
procedure TdEVDEV2(const A : Double;
     const B : Double;
     const C : Double;
     var RT1 : Double;
     var RT2 : Double;
     var CS1 : Double;
     var SN1 : Double);
var
    SGN1 : AlglibInteger;
    SGN2 : AlglibInteger;
    AB : Double;
    ACMN : Double;
    ACMX : Double;
    ACS : Double;
    ADF : Double;
    CS : Double;
    CT : Double;
    DF : Double;
    RT : Double;
    SM : Double;
    TB : Double;
    TN : Double;
begin
    
    //
    // Compute the eigenvalues
    //
    SM := A+C;
    DF := A-C;
    ADF := AbsReal(DF);
    TB := B+B;
    AB := AbsReal(TB);
    if AP_FP_Greater(AbsReal(A),AbsReal(C)) then
    begin
        ACMX := A;
        ACMN := C;
    end
    else
    begin
        ACMX := C;
        ACMN := A;
    end;
    if AP_FP_Greater(ADF,AB) then
    begin
        RT := ADF*Sqrt(1+AP_Sqr(AB/ADF));
    end
    else
    begin
        if AP_FP_Less(ADF,AB) then
        begin
            RT := AB*Sqrt(1+AP_Sqr(ADF/AB));
        end
        else
        begin
            
            //
            // Includes case AB=ADF=0
            //
            RT := AB*Sqrt(2);
        end;
    end;
    if AP_FP_Less(SM,0) then
    begin
        RT1 := 0.5*(SM-RT);
        SGN1 := -1;
        
        //
        // Order of execution important.
        // To get fully accurate smaller eigenvalue,
        // next line needs to be executed in higher precision.
        //
        RT2 := ACMX/RT1*ACMN-B/RT1*B;
    end
    else
    begin
        if AP_FP_Greater(SM,0) then
        begin
            RT1 := 0.5*(SM+RT);
            SGN1 := 1;
            
            //
            // Order of execution important.
            // To get fully accurate smaller eigenvalue,
            // next line needs to be executed in higher precision.
            //
            RT2 := ACMX/RT1*ACMN-B/RT1*B;
        end
        else
        begin
            
            //
            // Includes case RT1 = RT2 = 0
            //
            RT1 := 0.5*RT;
            RT2 := -0.5*RT;
            SGN1 := 1;
        end;
    end;
    
    //
    // Compute the eigenvector
    //
    if AP_FP_Greater_Eq(DF,0) then
    begin
        CS := DF+RT;
        SGN2 := 1;
    end
    else
    begin
        CS := DF-RT;
        SGN2 := -1;
    end;
    ACS := AbsReal(CS);
    if AP_FP_Greater(ACS,AB) then
    begin
        CT := -TB/CS;
        SN1 := 1/SQRT(1+CT*CT);
        CS1 := CT*SN1;
    end
    else
    begin
        if AP_FP_Eq(AB,0) then
        begin
            CS1 := 1;
            SN1 := 0;
        end
        else
        begin
            TN := -CS/TB;
            CS1 := 1/SQRT(1+TN*TN);
            SN1 := TN*CS1;
        end;
    end;
    if SGN1=SGN2 then
    begin
        TN := CS1;
        CS1 := -SN1;
        SN1 := TN;
    end;
end;


(*************************************************************************
Internal routine
*************************************************************************)
function TdEVDPythag(A : Double; B : Double):Double;
begin
    if AP_FP_Less(AbsReal(A),AbsReal(B)) then
    begin
        Result := AbsReal(B)*Sqrt(1+AP_Sqr(A/B));
    end
    else
    begin
        Result := AbsReal(A)*Sqrt(1+AP_Sqr(B/A));
    end;
end;


(*************************************************************************
Internal routine
*************************************************************************)
function TdEVDExtSign(a : Double; b : Double):Double;
begin
    if AP_FP_Greater_Eq(b,0) then
    begin
        Result := AbsReal(a);
    end
    else
    begin
        Result := -AbsReal(a);
    end;
end;


end.