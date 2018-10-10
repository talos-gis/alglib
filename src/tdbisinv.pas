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
unit tdbisinv;
interface
uses Math, Sysutils, Ap, blas;

function SMatrixTDEVDR(var D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     A : Double;
     B : Double;
     var M : AlglibInteger;
     var Z : TReal2DArray):Boolean;
function SMatrixTDEVDI(var D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     var Z : TReal2DArray):Boolean;
function TridiagonalEigenValuesAndVectorsInInterval(var D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     A : Double;
     B : Double;
     var M : AlglibInteger;
     var Z : TReal2DArray):Boolean;
function TridiagonalEigenValuesAndVectorsByIndexes(var D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     var Z : TReal2DArray):Boolean;
function InternalBisectionEigenValues(D : TReal1DArray;
     E : TReal1DArray;
     N : AlglibInteger;
     IRANGE : AlglibInteger;
     IORDER : AlglibInteger;
     VL : Double;
     VU : Double;
     IL : AlglibInteger;
     IU : AlglibInteger;
     ABSTOL : Double;
     var W : TReal1DArray;
     var M : AlglibInteger;
     var NSPLIT : AlglibInteger;
     var IBLOCK : TInteger1DArray;
     var ISPLIT : TInteger1DArray;
     var ErrorCode : AlglibInteger):Boolean;
procedure InternalDSTEIN(const N : AlglibInteger;
     const D : TReal1DArray;
     E : TReal1DArray;
     const M : AlglibInteger;
     W : TReal1DArray;
     const IBLOCK : TInteger1DArray;
     const ISPLIT : TInteger1DArray;
     var Z : TReal2DArray;
     var IFAIL : TInteger1DArray;
     var INFO : AlglibInteger);

implementation

procedure TDINInternalDLAGTF(const N : AlglibInteger;
     var A : TReal1DArray;
     const LAMBDA : Double;
     var B : TReal1DArray;
     var C : TReal1DArray;
     const TOL : Double;
     var D : TReal1DArray;
     var IIN : TInteger1DArray;
     var INFO : AlglibInteger);forward;
procedure TDINInternalDLAGTS(const N : AlglibInteger;
     const A : TReal1DArray;
     const B : TReal1DArray;
     const C : TReal1DArray;
     const D : TReal1DArray;
     const IIN : TInteger1DArray;
     var Y : TReal1DArray;
     var TOL : Double;
     var INFO : AlglibInteger);forward;
procedure InternalDLAEBZ(const IJOB : AlglibInteger;
     const NITMAX : AlglibInteger;
     const N : AlglibInteger;
     const MMAX : AlglibInteger;
     const MINP : AlglibInteger;
     const ABSTOL : Double;
     const RELTOL : Double;
     const PIVMIN : Double;
     const D : TReal1DArray;
     const E : TReal1DArray;
     const E2 : TReal1DArray;
     var NVAL : TInteger1DArray;
     var AB : TReal2DArray;
     var C : TReal1DArray;
     var MOUT : AlglibInteger;
     var NAB : TInteger2DArray;
     var WORK : TReal1DArray;
     var IWORK : TInteger1DArray;
     var INFO : AlglibInteger);forward;


(*************************************************************************
Subroutine for finding the tridiagonal matrix eigenvalues/vectors in a
given half-interval (A, B] by using bisection and inverse iteration.

Input parameters:
    D       -   the main diagonal of a tridiagonal matrix.
                Array whose index ranges within [0..N-1].
    E       -   the secondary diagonal of a tridiagonal matrix.
                Array whose index ranges within [0..N-2].
    N       -   size of matrix, N>=0.
    ZNeeded -   flag controlling whether the eigenvectors are needed or not.
                If ZNeeded is equal to:
                 * 0, the eigenvectors are not needed;
                 * 1, the eigenvectors of a tridiagonal matrix are multiplied
                   by the square matrix Z. It is used if the tridiagonal
                   matrix is obtained by the similarity transformation
                   of a symmetric matrix.
                 * 2, the eigenvectors of a tridiagonal matrix replace matrix Z.
    A, B    -   half-interval (A, B] to search eigenvalues in.
    Z       -   if ZNeeded is equal to:
                 * 0, Z isn't used and remains unchanged;
                 * 1, Z contains the square matrix (array whose indexes range
                   within [0..N-1, 0..N-1]) which reduces the given symmetric
                   matrix to tridiagonal form;
                 * 2, Z isn't used (but changed on the exit).

Output parameters:
    D       -   array of the eigenvalues found.
                Array whose index ranges within [0..M-1].
    M       -   number of eigenvalues found in the given half-interval (M>=0).
    Z       -   if ZNeeded is equal to:
                 * 0, doesn't contain any information;
                 * 1, contains the product of a given NxN matrix Z (from the
                   left) and NxM matrix of the eigenvectors found (from the
                   right). Array whose indexes range within [0..N-1, 0..M-1].
                 * 2, contains the matrix of the eigenvectors found.
                   Array whose indexes range within [0..N-1, 0..M-1].

Result:

    True, if successful. In that case, M contains the number of eigenvalues
    in the given half-interval (could be equal to 0), D contains the eigenvalues,
    Z contains the eigenvectors (if needed).
    It should be noted that the subroutine changes the size of arrays D and Z.

    False, if the bisection method subroutine wasn't able to find the
    eigenvalues in the given interval or if the inverse iteration subroutine
    wasn't able to find all the corresponding eigenvectors. In that case,
    the eigenvalues and eigenvectors are not returned, M is equal to 0.

  -- ALGLIB --
     Copyright 31.03.2008 by Bochkanov Sergey
*************************************************************************)
function SMatrixTDEVDR(var D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     A : Double;
     B : Double;
     var M : AlglibInteger;
     var Z : TReal2DArray):Boolean;
var
    ErrorCode : AlglibInteger;
    NSPLIT : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    CR : AlglibInteger;
    IBLOCK : TInteger1DArray;
    ISPLIT : TInteger1DArray;
    IFAIL : TInteger1DArray;
    D1 : TReal1DArray;
    E1 : TReal1DArray;
    W : TReal1DArray;
    Z2 : TReal2DArray;
    Z3 : TReal2DArray;
    V : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert((ZNeeded>=0) and (ZNeeded<=2), 'SMatrixTDEVDR: incorrect ZNeeded!');
    
    //
    // Special cases
    //
    if AP_FP_Less_Eq(B,A) then
    begin
        M := 0;
        Result := True;
        Exit;
    end;
    if N<=0 then
    begin
        M := 0;
        Result := True;
        Exit;
    end;
    
    //
    // Copy D,E to D1, E1
    //
    SetLength(D1, N+1);
    APVMove(@D1[0], 1, N, @D[0], 0, N-1);
    if N>1 then
    begin
        SetLength(E1, N-1+1);
        APVMove(@E1[0], 1, N-1, @E[0], 0, N-2);
    end;
    
    //
    // No eigen vectors
    //
    if ZNeeded=0 then
    begin
        Result := InternalBisectionEigenValues(D1, E1, N, 2, 1, A, B, 0, 0, -1, W, M, NSPLIT, IBLOCK, ISPLIT, ErrorCode);
        if  not Result or (M=0) then
        begin
            M := 0;
            Exit;
        end;
        SetLength(D, M-1+1);
        APVMove(@D[0], 0, M-1, @W[0], 1, M);
        Exit;
    end;
    
    //
    // Eigen vectors are multiplied by Z
    //
    if ZNeeded=1 then
    begin
        
        //
        // Find eigen pairs
        //
        Result := InternalBisectionEigenValues(D1, E1, N, 2, 2, A, B, 0, 0, -1, W, M, NSPLIT, IBLOCK, ISPLIT, ErrorCode);
        if  not Result or (M=0) then
        begin
            M := 0;
            Exit;
        end;
        InternalDSTEIN(N, D1, E1, M, W, IBLOCK, ISPLIT, Z2, IFAIL, CR);
        if CR<>0 then
        begin
            M := 0;
            Result := False;
            Exit;
        end;
        
        //
        // Sort eigen values and vectors
        //
        I:=1;
        while I<=M do
        begin
            K := I;
            J:=I;
            while J<=M do
            begin
                if AP_FP_Less(W[J],W[K]) then
                begin
                    K := J;
                end;
                Inc(J);
            end;
            V := W[I];
            W[I] := W[K];
            W[K] := V;
            J:=1;
            while J<=N do
            begin
                V := Z2[J,I];
                Z2[J,I] := Z2[J,K];
                Z2[J,K] := V;
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Transform Z2 and overwrite Z
        //
        SetLength(Z3, M+1, N+1);
        I:=1;
        while I<=M do
        begin
            for i_ := 1 to N do
            begin
                Z3[I,i_] := Z2[i_,I];
            end;
            Inc(I);
        end;
        I:=1;
        while I<=N do
        begin
            J:=1;
            while J<=M do
            begin
                V := APVDotProduct(@Z[I-1][0], 0, N-1, @Z3[J][0], 1, N);
                Z2[I,J] := V;
                Inc(J);
            end;
            Inc(I);
        end;
        SetLength(Z, N-1+1, M-1+1);
        I:=1;
        while I<=M do
        begin
            i1_ := (1) - (0);
            for i_ := 0 to N-1 do
            begin
                Z[i_,I-1] := Z2[i_+i1_,I];
            end;
            Inc(I);
        end;
        
        //
        // Store W
        //
        SetLength(D, M-1+1);
        I:=1;
        while I<=M do
        begin
            D[I-1] := W[I];
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // Eigen vectors are stored in Z
    //
    if ZNeeded=2 then
    begin
        
        //
        // Find eigen pairs
        //
        Result := InternalBisectionEigenValues(D1, E1, N, 2, 2, A, B, 0, 0, -1, W, M, NSPLIT, IBLOCK, ISPLIT, ErrorCode);
        if  not Result or (M=0) then
        begin
            M := 0;
            Exit;
        end;
        InternalDSTEIN(N, D1, E1, M, W, IBLOCK, ISPLIT, Z2, IFAIL, CR);
        if CR<>0 then
        begin
            M := 0;
            Result := False;
            Exit;
        end;
        
        //
        // Sort eigen values and vectors
        //
        I:=1;
        while I<=M do
        begin
            K := I;
            J:=I;
            while J<=M do
            begin
                if AP_FP_Less(W[J],W[K]) then
                begin
                    K := J;
                end;
                Inc(J);
            end;
            V := W[I];
            W[I] := W[K];
            W[K] := V;
            J:=1;
            while J<=N do
            begin
                V := Z2[J,I];
                Z2[J,I] := Z2[J,K];
                Z2[J,K] := V;
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Store W
        //
        SetLength(D, M-1+1);
        I:=1;
        while I<=M do
        begin
            D[I-1] := W[I];
            Inc(I);
        end;
        SetLength(Z, N-1+1, M-1+1);
        I:=1;
        while I<=M do
        begin
            i1_ := (1) - (0);
            for i_ := 0 to N-1 do
            begin
                Z[i_,I-1] := Z2[i_+i1_,I];
            end;
            Inc(I);
        end;
        Exit;
    end;
    Result := False;
end;


(*************************************************************************
Subroutine for finding tridiagonal matrix eigenvalues/vectors with given
indexes (in ascending order) by using the bisection and inverse iteraion.

Input parameters:
    D       -   the main diagonal of a tridiagonal matrix.
                Array whose index ranges within [0..N-1].
    E       -   the secondary diagonal of a tridiagonal matrix.
                Array whose index ranges within [0..N-2].
    N       -   size of matrix. N>=0.
    ZNeeded -   flag controlling whether the eigenvectors are needed or not.
                If ZNeeded is equal to:
                 * 0, the eigenvectors are not needed;
                 * 1, the eigenvectors of a tridiagonal matrix are multiplied
                   by the square matrix Z. It is used if the
                   tridiagonal matrix is obtained by the similarity transformation
                   of a symmetric matrix.
                 * 2, the eigenvectors of a tridiagonal matrix replace
                   matrix Z.
    I1, I2  -   index interval for searching (from I1 to I2).
                0 <= I1 <= I2 <= N-1.
    Z       -   if ZNeeded is equal to:
                 * 0, Z isn't used and remains unchanged;
                 * 1, Z contains the square matrix (array whose indexes range within [0..N-1, 0..N-1])
                   which reduces the given symmetric matrix to  tridiagonal form;
                 * 2, Z isn't used (but changed on the exit).

Output parameters:
    D       -   array of the eigenvalues found.
                Array whose index ranges within [0..I2-I1].
    Z       -   if ZNeeded is equal to:
                 * 0, doesn't contain any information;
                 * 1, contains the product of a given NxN matrix Z (from the left) and
                   Nx(I2-I1) matrix of the eigenvectors found (from the right).
                   Array whose indexes range within [0..N-1, 0..I2-I1].
                 * 2, contains the matrix of the eigenvalues found.
                   Array whose indexes range within [0..N-1, 0..I2-I1].


Result:

    True, if successful. In that case, D contains the eigenvalues,
    Z contains the eigenvectors (if needed).
    It should be noted that the subroutine changes the size of arrays D and Z.

    False, if the bisection method subroutine wasn't able to find the eigenvalues
    in the given interval or if the inverse iteration subroutine wasn't able
    to find all the corresponding eigenvectors. In that case, the eigenvalues
    and eigenvectors are not returned.

  -- ALGLIB --
     Copyright 25.12.2005 by Bochkanov Sergey
*************************************************************************)
function SMatrixTDEVDI(var D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     var Z : TReal2DArray):Boolean;
var
    ErrorCode : AlglibInteger;
    NSPLIT : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    M : AlglibInteger;
    CR : AlglibInteger;
    IBLOCK : TInteger1DArray;
    ISPLIT : TInteger1DArray;
    IFAIL : TInteger1DArray;
    W : TReal1DArray;
    D1 : TReal1DArray;
    E1 : TReal1DArray;
    Z2 : TReal2DArray;
    Z3 : TReal2DArray;
    V : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Assert((0<=I1) and (I1<=I2) and (I2<N), 'SMatrixTDEVDI: incorrect I1/I2!');
    
    //
    // Copy D,E to D1, E1
    //
    SetLength(D1, N+1);
    APVMove(@D1[0], 1, N, @D[0], 0, N-1);
    if N>1 then
    begin
        SetLength(E1, N-1+1);
        APVMove(@E1[0], 1, N-1, @E[0], 0, N-2);
    end;
    
    //
    // No eigen vectors
    //
    if ZNeeded=0 then
    begin
        Result := InternalBisectionEigenValues(D1, E1, N, 3, 1, 0, 0, I1+1, I2+1, -1, W, M, NSPLIT, IBLOCK, ISPLIT, ErrorCode);
        if  not Result then
        begin
            Exit;
        end;
        if M<>I2-I1+1 then
        begin
            Result := False;
            Exit;
        end;
        SetLength(D, M-1+1);
        I:=1;
        while I<=M do
        begin
            D[I-1] := W[I];
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // Eigen vectors are multiplied by Z
    //
    if ZNeeded=1 then
    begin
        
        //
        // Find eigen pairs
        //
        Result := InternalBisectionEigenValues(D1, E1, N, 3, 2, 0, 0, I1+1, I2+1, -1, W, M, NSPLIT, IBLOCK, ISPLIT, ErrorCode);
        if  not Result then
        begin
            Exit;
        end;
        if M<>I2-I1+1 then
        begin
            Result := False;
            Exit;
        end;
        InternalDSTEIN(N, D1, E1, M, W, IBLOCK, ISPLIT, Z2, IFAIL, CR);
        if CR<>0 then
        begin
            Result := False;
            Exit;
        end;
        
        //
        // Sort eigen values and vectors
        //
        I:=1;
        while I<=M do
        begin
            K := I;
            J:=I;
            while J<=M do
            begin
                if AP_FP_Less(W[J],W[K]) then
                begin
                    K := J;
                end;
                Inc(J);
            end;
            V := W[I];
            W[I] := W[K];
            W[K] := V;
            J:=1;
            while J<=N do
            begin
                V := Z2[J,I];
                Z2[J,I] := Z2[J,K];
                Z2[J,K] := V;
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Transform Z2 and overwrite Z
        //
        SetLength(Z3, M+1, N+1);
        I:=1;
        while I<=M do
        begin
            for i_ := 1 to N do
            begin
                Z3[I,i_] := Z2[i_,I];
            end;
            Inc(I);
        end;
        I:=1;
        while I<=N do
        begin
            J:=1;
            while J<=M do
            begin
                V := APVDotProduct(@Z[I-1][0], 0, N-1, @Z3[J][0], 1, N);
                Z2[I,J] := V;
                Inc(J);
            end;
            Inc(I);
        end;
        SetLength(Z, N-1+1, M-1+1);
        I:=1;
        while I<=M do
        begin
            i1_ := (1) - (0);
            for i_ := 0 to N-1 do
            begin
                Z[i_,I-1] := Z2[i_+i1_,I];
            end;
            Inc(I);
        end;
        
        //
        // Store W
        //
        SetLength(D, M-1+1);
        I:=1;
        while I<=M do
        begin
            D[I-1] := W[I];
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // Eigen vectors are stored in Z
    //
    if ZNeeded=2 then
    begin
        
        //
        // Find eigen pairs
        //
        Result := InternalBisectionEigenValues(D1, E1, N, 3, 2, 0, 0, I1+1, I2+1, -1, W, M, NSPLIT, IBLOCK, ISPLIT, ErrorCode);
        if  not Result then
        begin
            Exit;
        end;
        if M<>I2-I1+1 then
        begin
            Result := False;
            Exit;
        end;
        InternalDSTEIN(N, D1, E1, M, W, IBLOCK, ISPLIT, Z2, IFAIL, CR);
        if CR<>0 then
        begin
            Result := False;
            Exit;
        end;
        
        //
        // Sort eigen values and vectors
        //
        I:=1;
        while I<=M do
        begin
            K := I;
            J:=I;
            while J<=M do
            begin
                if AP_FP_Less(W[J],W[K]) then
                begin
                    K := J;
                end;
                Inc(J);
            end;
            V := W[I];
            W[I] := W[K];
            W[K] := V;
            J:=1;
            while J<=N do
            begin
                V := Z2[J,I];
                Z2[J,I] := Z2[J,K];
                Z2[J,K] := V;
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Store Z
        //
        SetLength(Z, N-1+1, M-1+1);
        I:=1;
        while I<=M do
        begin
            i1_ := (1) - (0);
            for i_ := 0 to N-1 do
            begin
                Z[i_,I-1] := Z2[i_+i1_,I];
            end;
            Inc(I);
        end;
        
        //
        // Store W
        //
        SetLength(D, M-1+1);
        I:=1;
        while I<=M do
        begin
            D[I-1] := W[I];
            Inc(I);
        end;
        Exit;
    end;
    Result := False;
end;


function TridiagonalEigenValuesAndVectorsInInterval(var D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     A : Double;
     B : Double;
     var M : AlglibInteger;
     var Z : TReal2DArray):Boolean;
var
    ErrorCode : AlglibInteger;
    NSPLIT : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    CR : AlglibInteger;
    IBLOCK : TInteger1DArray;
    ISPLIT : TInteger1DArray;
    IFAIL : TInteger1DArray;
    W : TReal1DArray;
    Z2 : TReal2DArray;
    Z3 : TReal2DArray;
    V : Double;
    i_ : AlglibInteger;
begin
    
    //
    // No eigen vectors
    //
    if ZNeeded=0 then
    begin
        Result := InternalBisectionEigenValues(D, E, N, 2, 1, A, B, 0, 0, -1, W, M, NSPLIT, IBLOCK, ISPLIT, ErrorCode);
        if  not Result or (M=0) then
        begin
            M := 0;
            Exit;
        end;
        SetLength(D, M+1);
        I:=1;
        while I<=M do
        begin
            D[I] := W[I];
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // Eigen vectors are multiplied by Z
    //
    if ZNeeded=1 then
    begin
        
        //
        // Find eigen pairs
        //
        Result := InternalBisectionEigenValues(D, E, N, 2, 2, A, B, 0, 0, -1, W, M, NSPLIT, IBLOCK, ISPLIT, ErrorCode);
        if  not Result or (M=0) then
        begin
            M := 0;
            Exit;
        end;
        InternalDSTEIN(N, D, E, M, W, IBLOCK, ISPLIT, Z2, IFAIL, CR);
        if CR<>0 then
        begin
            M := 0;
            Result := False;
            Exit;
        end;
        
        //
        // Sort eigen values and vectors
        //
        I:=1;
        while I<=M do
        begin
            K := I;
            J:=I;
            while J<=M do
            begin
                if AP_FP_Less(W[J],W[K]) then
                begin
                    K := J;
                end;
                Inc(J);
            end;
            V := W[I];
            W[I] := W[K];
            W[K] := V;
            J:=1;
            while J<=N do
            begin
                V := Z2[J,I];
                Z2[J,I] := Z2[J,K];
                Z2[J,K] := V;
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Transform Z2 and overwrite Z
        //
        SetLength(Z3, M+1, N+1);
        I:=1;
        while I<=M do
        begin
            for i_ := 1 to N do
            begin
                Z3[I,i_] := Z2[i_,I];
            end;
            Inc(I);
        end;
        I:=1;
        while I<=N do
        begin
            J:=1;
            while J<=M do
            begin
                V := APVDotProduct(@Z[I][0], 1, N, @Z3[J][0], 1, N);
                Z2[I,J] := V;
                Inc(J);
            end;
            Inc(I);
        end;
        SetLength(Z, N+1, M+1);
        I:=1;
        while I<=M do
        begin
            for i_ := 1 to N do
            begin
                Z[i_,I] := Z2[i_,I];
            end;
            Inc(I);
        end;
        
        //
        // Store W
        //
        SetLength(D, M+1);
        I:=1;
        while I<=M do
        begin
            D[I] := W[I];
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // Eigen vectors are stored in Z
    //
    if ZNeeded=2 then
    begin
        
        //
        // Find eigen pairs
        //
        Result := InternalBisectionEigenValues(D, E, N, 2, 2, A, B, 0, 0, -1, W, M, NSPLIT, IBLOCK, ISPLIT, ErrorCode);
        if  not Result or (M=0) then
        begin
            M := 0;
            Exit;
        end;
        InternalDSTEIN(N, D, E, M, W, IBLOCK, ISPLIT, Z, IFAIL, CR);
        if CR<>0 then
        begin
            M := 0;
            Result := False;
            Exit;
        end;
        
        //
        // Sort eigen values and vectors
        //
        I:=1;
        while I<=M do
        begin
            K := I;
            J:=I;
            while J<=M do
            begin
                if AP_FP_Less(W[J],W[K]) then
                begin
                    K := J;
                end;
                Inc(J);
            end;
            V := W[I];
            W[I] := W[K];
            W[K] := V;
            J:=1;
            while J<=N do
            begin
                V := Z[J,I];
                Z[J,I] := Z[J,K];
                Z[J,K] := V;
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Store W
        //
        SetLength(D, M+1);
        I:=1;
        while I<=M do
        begin
            D[I] := W[I];
            Inc(I);
        end;
        Exit;
    end;
    Result := False;
end;


function TridiagonalEigenValuesAndVectorsByIndexes(var D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     ZNeeded : AlglibInteger;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     var Z : TReal2DArray):Boolean;
var
    ErrorCode : AlglibInteger;
    NSPLIT : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    M : AlglibInteger;
    CR : AlglibInteger;
    IBLOCK : TInteger1DArray;
    ISPLIT : TInteger1DArray;
    IFAIL : TInteger1DArray;
    W : TReal1DArray;
    Z2 : TReal2DArray;
    Z3 : TReal2DArray;
    V : Double;
    i_ : AlglibInteger;
begin
    
    //
    // No eigen vectors
    //
    if ZNeeded=0 then
    begin
        Result := InternalBisectionEigenValues(D, E, N, 3, 1, 0, 0, I1, I2, -1, W, M, NSPLIT, IBLOCK, ISPLIT, ErrorCode);
        if  not Result then
        begin
            Exit;
        end;
        SetLength(D, M+1);
        I:=1;
        while I<=M do
        begin
            D[I] := W[I];
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // Eigen vectors are multiplied by Z
    //
    if ZNeeded=1 then
    begin
        
        //
        // Find eigen pairs
        //
        Result := InternalBisectionEigenValues(D, E, N, 3, 2, 0, 0, I1, I2, -1, W, M, NSPLIT, IBLOCK, ISPLIT, ErrorCode);
        if  not Result then
        begin
            Exit;
        end;
        InternalDSTEIN(N, D, E, M, W, IBLOCK, ISPLIT, Z2, IFAIL, CR);
        if CR<>0 then
        begin
            Result := False;
            Exit;
        end;
        
        //
        // Sort eigen values and vectors
        //
        I:=1;
        while I<=M do
        begin
            K := I;
            J:=I;
            while J<=M do
            begin
                if AP_FP_Less(W[J],W[K]) then
                begin
                    K := J;
                end;
                Inc(J);
            end;
            V := W[I];
            W[I] := W[K];
            W[K] := V;
            J:=1;
            while J<=N do
            begin
                V := Z2[J,I];
                Z2[J,I] := Z2[J,K];
                Z2[J,K] := V;
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Transform Z2 and overwrite Z
        //
        SetLength(Z3, M+1, N+1);
        I:=1;
        while I<=M do
        begin
            for i_ := 1 to N do
            begin
                Z3[I,i_] := Z2[i_,I];
            end;
            Inc(I);
        end;
        I:=1;
        while I<=N do
        begin
            J:=1;
            while J<=M do
            begin
                V := APVDotProduct(@Z[I][0], 1, N, @Z3[J][0], 1, N);
                Z2[I,J] := V;
                Inc(J);
            end;
            Inc(I);
        end;
        SetLength(Z, N+1, M+1);
        I:=1;
        while I<=M do
        begin
            for i_ := 1 to N do
            begin
                Z[i_,I] := Z2[i_,I];
            end;
            Inc(I);
        end;
        
        //
        // Store W
        //
        SetLength(D, M+1);
        I:=1;
        while I<=M do
        begin
            D[I] := W[I];
            Inc(I);
        end;
        Exit;
    end;
    
    //
    // Eigen vectors are stored in Z
    //
    if ZNeeded=2 then
    begin
        
        //
        // Find eigen pairs
        //
        Result := InternalBisectionEigenValues(D, E, N, 3, 2, 0, 0, I1, I2, -1, W, M, NSPLIT, IBLOCK, ISPLIT, ErrorCode);
        if  not Result then
        begin
            Exit;
        end;
        InternalDSTEIN(N, D, E, M, W, IBLOCK, ISPLIT, Z, IFAIL, CR);
        if CR<>0 then
        begin
            Result := False;
            Exit;
        end;
        
        //
        // Sort eigen values and vectors
        //
        I:=1;
        while I<=M do
        begin
            K := I;
            J:=I;
            while J<=M do
            begin
                if AP_FP_Less(W[J],W[K]) then
                begin
                    K := J;
                end;
                Inc(J);
            end;
            V := W[I];
            W[I] := W[K];
            W[K] := V;
            J:=1;
            while J<=N do
            begin
                V := Z[J,I];
                Z[J,I] := Z[J,K];
                Z[J,K] := V;
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Store W
        //
        SetLength(D, M+1);
        I:=1;
        while I<=M do
        begin
            D[I] := W[I];
            Inc(I);
        end;
        Exit;
    end;
    Result := False;
end;


function InternalBisectionEigenValues(D : TReal1DArray;
     E : TReal1DArray;
     N : AlglibInteger;
     IRANGE : AlglibInteger;
     IORDER : AlglibInteger;
     VL : Double;
     VU : Double;
     IL : AlglibInteger;
     IU : AlglibInteger;
     ABSTOL : Double;
     var W : TReal1DArray;
     var M : AlglibInteger;
     var NSPLIT : AlglibInteger;
     var IBLOCK : TInteger1DArray;
     var ISPLIT : TInteger1DArray;
     var ErrorCode : AlglibInteger):Boolean;
var
    FUDGE : Double;
    RELFAC : Double;
    NCNVRG : Boolean;
    TOOFEW : Boolean;
    IB : AlglibInteger;
    IBEGIN : AlglibInteger;
    IDISCL : AlglibInteger;
    IDISCU : AlglibInteger;
    IE : AlglibInteger;
    IEND : AlglibInteger;
    IINFO : AlglibInteger;
    IM : AlglibInteger;
    IIN : AlglibInteger;
    IOFF : AlglibInteger;
    IOUT : AlglibInteger;
    ITMAX : AlglibInteger;
    IW : AlglibInteger;
    IWOFF : AlglibInteger;
    J : AlglibInteger;
    ITMP1 : AlglibInteger;
    JB : AlglibInteger;
    JDISC : AlglibInteger;
    JE : AlglibInteger;
    NWL : AlglibInteger;
    NWU : AlglibInteger;
    ATOLI : Double;
    BNORM : Double;
    GL : Double;
    GU : Double;
    PIVMIN : Double;
    RTOLI : Double;
    SAFEMN : Double;
    TMP1 : Double;
    TMP2 : Double;
    TNORM : Double;
    ULP : Double;
    WKILL : Double;
    WL : Double;
    WLU : Double;
    WU : Double;
    WUL : Double;
    ScaleFactor : Double;
    T : Double;
    IDUMMA : TInteger1DArray;
    WORK : TReal1DArray;
    IWORK : TInteger1DArray;
    IA1S2 : TInteger1DArray;
    RA1S2 : TReal1DArray;
    RA1S2X2 : TReal2DArray;
    IA1S2X2 : TInteger2DArray;
    RA1SIIN : TReal1DArray;
    RA2SIIN : TReal1DArray;
    RA3SIIN : TReal1DArray;
    RA4SIIN : TReal1DArray;
    RA1SIINX2 : TReal2DArray;
    IA1SIINX2 : TInteger2DArray;
    IWORKSPACE : TInteger1DArray;
    RWORKSPACE : TReal1DArray;
    TmpI : AlglibInteger;
begin
    D := DynamicArrayCopy(D);
    E := DynamicArrayCopy(E);
    
    //
    // Quick return if possible
    //
    M := 0;
    if N=0 then
    begin
        Result := True;
        Exit;
    end;
    
    //
    // Get machine constants
    // NB is the minimum vector length for vector bisection, or 0
    // if only scalar is to be done.
    //
    FUDGE := 2;
    RELFAC := 2;
    SAFEMN := MinRealNumber;
    ULP := 2*MachineEpsilon;
    RTOLI := ULP*RELFAC;
    SetLength(IDUMMA, 1+1);
    SetLength(WORK, 4*N+1);
    SetLength(IWORK, 3*N+1);
    SetLength(W, N+1);
    SetLength(IBLOCK, N+1);
    SetLength(ISPLIT, N+1);
    SetLength(IA1S2, 2+1);
    SetLength(RA1S2, 2+1);
    SetLength(RA1S2X2, 2+1, 2+1);
    SetLength(IA1S2X2, 2+1, 2+1);
    SetLength(RA1SIIN, N+1);
    SetLength(RA2SIIN, N+1);
    SetLength(RA3SIIN, N+1);
    SetLength(RA4SIIN, N+1);
    SetLength(RA1SIINX2, N+1, 2+1);
    SetLength(IA1SIINX2, N+1, 2+1);
    SetLength(IWORKSPACE, N+1);
    SetLength(RWORKSPACE, N+1);
    
    //
    // Check for Errors
    //
    Result := False;
    ErrorCode := 0;
    if (IRANGE<=0) or (IRANGE>=4) then
    begin
        ErrorCode := -4;
    end;
    if (IORDER<=0) or (IORDER>=3) then
    begin
        ErrorCode := -5;
    end;
    if N<0 then
    begin
        ErrorCode := -3;
    end;
    if (IRANGE=2) and AP_FP_Greater_Eq(VL,VU) then
    begin
        ErrorCode := -6;
    end;
    if (IRANGE=3) and ((IL<1) or (IL>Max(1, N))) then
    begin
        ErrorCode := -8;
    end;
    if (IRANGE=3) and ((IU<Min(N, IL)) or (IU>N)) then
    begin
        ErrorCode := -9;
    end;
    if ErrorCode<>0 then
    begin
        Exit;
    end;
    
    //
    // Initialize error flags
    //
    NCNVRG := False;
    TOOFEW := False;
    
    //
    // Simplifications:
    //
    if (IRANGE=3) and (IL=1) and (IU=N) then
    begin
        IRANGE := 1;
    end;
    
    //
    // Special Case when N=1
    //
    if N=1 then
    begin
        NSPLIT := 1;
        ISPLIT[1] := 1;
        if (IRANGE=2) and (AP_FP_Greater_Eq(VL,D[1]) or AP_FP_Less(VU,D[1])) then
        begin
            M := 0;
        end
        else
        begin
            W[1] := D[1];
            IBLOCK[1] := 1;
            M := 1;
        end;
        Result := True;
        Exit;
    end;
    
    //
    // Scaling
    //
    T := AbsReal(D[N]);
    J:=1;
    while J<=N-1 do
    begin
        T := Max(T, AbsReal(D[J]));
        T := Max(T, AbsReal(E[J]));
        Inc(J);
    end;
    ScaleFactor := 1;
    if AP_FP_Neq(T,0) then
    begin
        if AP_FP_Greater(T,Sqrt(Sqrt(MinRealNumber))*Sqrt(MaxRealNumber)) then
        begin
            ScaleFactor := T;
        end;
        if AP_FP_Less(T,Sqrt(Sqrt(MaxRealNumber))*Sqrt(MinRealNumber)) then
        begin
            ScaleFactor := T;
        end;
        J:=1;
        while J<=N-1 do
        begin
            D[J] := D[J]/ScaleFactor;
            E[J] := E[J]/ScaleFactor;
            Inc(J);
        end;
        D[N] := D[N]/ScaleFactor;
    end;
    
    //
    // Compute Splitting Points
    //
    NSPLIT := 1;
    WORK[N] := 0;
    PIVMIN := 1;
    J:=2;
    while J<=N do
    begin
        TMP1 := AP_Sqr(E[J-1]);
        if AP_FP_Greater(AbsReal(D[J]*D[J-1])*AP_Sqr(ULP)+SAFEMN,TMP1) then
        begin
            ISPLIT[NSPLIT] := J-1;
            NSPLIT := NSPLIT+1;
            WORK[J-1] := 0;
        end
        else
        begin
            WORK[J-1] := TMP1;
            PIVMIN := Max(PIVMIN, TMP1);
        end;
        Inc(J);
    end;
    ISPLIT[NSPLIT] := N;
    PIVMIN := PIVMIN*SAFEMN;
    
    //
    // Compute Interval and ATOLI
    //
    if IRANGE=3 then
    begin
        
        //
        // RANGE='I': Compute the interval containing eigenvalues
        //     IL through IU.
        //
        // Compute Gershgorin interval for entire (split) matrix
        // and use it as the initial interval
        //
        GU := D[1];
        GL := D[1];
        TMP1 := 0;
        J:=1;
        while J<=N-1 do
        begin
            TMP2 := Sqrt(WORK[J]);
            GU := Max(GU, D[J]+TMP1+TMP2);
            GL := Min(GL, D[J]-TMP1-TMP2);
            TMP1 := TMP2;
            Inc(J);
        end;
        GU := Max(GU, D[N]+TMP1);
        GL := Min(GL, D[N]-TMP1);
        TNORM := Max(AbsReal(GL), AbsReal(GU));
        GL := GL-FUDGE*TNORM*ULP*N-FUDGE*2*PIVMIN;
        GU := GU+FUDGE*TNORM*ULP*N+FUDGE*PIVMIN;
        
        //
        // Compute Iteration parameters
        //
        ITMAX := Ceil((Ln(TNORM+PIVMIN)-Ln(PIVMIN))/Ln(2))+2;
        if AP_FP_Less_Eq(ABSTOL,0) then
        begin
            ATOLI := ULP*TNORM;
        end
        else
        begin
            ATOLI := ABSTOL;
        end;
        WORK[N+1] := GL;
        WORK[N+2] := GL;
        WORK[N+3] := GU;
        WORK[N+4] := GU;
        WORK[N+5] := GL;
        WORK[N+6] := GU;
        IWORK[1] := -1;
        IWORK[2] := -1;
        IWORK[3] := N+1;
        IWORK[4] := N+1;
        IWORK[5] := IL-1;
        IWORK[6] := IU;
        
        //
        // Calling DLAEBZ
        //
        // DLAEBZ( 3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D, E,
        //    WORK, IWORK( 5 ), WORK( N+1 ), WORK( N+5 ), IOUT,
        //    IWORK, W, IBLOCK, IINFO )
        //
        IA1S2[1] := IWORK[5];
        IA1S2[2] := IWORK[6];
        RA1S2[1] := WORK[N+5];
        RA1S2[2] := WORK[N+6];
        RA1S2X2[1,1] := WORK[N+1];
        RA1S2X2[2,1] := WORK[N+2];
        RA1S2X2[1,2] := WORK[N+3];
        RA1S2X2[2,2] := WORK[N+4];
        IA1S2X2[1,1] := IWORK[1];
        IA1S2X2[2,1] := IWORK[2];
        IA1S2X2[1,2] := IWORK[3];
        IA1S2X2[2,2] := IWORK[4];
        InternalDLAEBZ(3, ITMAX, N, 2, 2, ATOLI, RTOLI, PIVMIN, D, E, WORK, IA1S2, RA1S2X2, RA1S2, IOUT, IA1S2X2, W, IBLOCK, IINFO);
        IWORK[5] := IA1S2[1];
        IWORK[6] := IA1S2[2];
        WORK[N+5] := RA1S2[1];
        WORK[N+6] := RA1S2[2];
        WORK[N+1] := RA1S2X2[1,1];
        WORK[N+2] := RA1S2X2[2,1];
        WORK[N+3] := RA1S2X2[1,2];
        WORK[N+4] := RA1S2X2[2,2];
        IWORK[1] := IA1S2X2[1,1];
        IWORK[2] := IA1S2X2[2,1];
        IWORK[3] := IA1S2X2[1,2];
        IWORK[4] := IA1S2X2[2,2];
        if IWORK[6]=IU then
        begin
            WL := WORK[N+1];
            WLU := WORK[N+3];
            NWL := IWORK[1];
            WU := WORK[N+4];
            WUL := WORK[N+2];
            NWU := IWORK[4];
        end
        else
        begin
            WL := WORK[N+2];
            WLU := WORK[N+4];
            NWL := IWORK[2];
            WU := WORK[N+3];
            WUL := WORK[N+1];
            NWU := IWORK[3];
        end;
        if (NWL<0) or (NWL>=N) or (NWU<1) or (NWU>N) then
        begin
            ErrorCode := 4;
            Result := False;
            Exit;
        end;
    end
    else
    begin
        
        //
        // RANGE='A' or 'V' -- Set ATOLI
        //
        TNORM := Max(ABSReal(D[1])+ABSReal(E[1]), ABSReal(D[N])+ABSReal(E[N-1]));
        J:=2;
        while J<=N-1 do
        begin
            TNORM := Max(TNORM, ABSReal(D[J])+ABSReal(E[J-1])+ABSReal(E[J]));
            Inc(J);
        end;
        if AP_FP_Less_Eq(ABSTOL,0) then
        begin
            ATOLI := ULP*TNORM;
        end
        else
        begin
            ATOLI := ABSTOL;
        end;
        if IRANGE=2 then
        begin
            WL := VL;
            WU := VU;
        end
        else
        begin
            WL := 0;
            WU := 0;
        end;
    end;
    
    //
    // Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
    // NWL accumulates the number of eigenvalues .le. WL,
    // NWU accumulates the number of eigenvalues .le. WU
    //
    M := 0;
    IEND := 0;
    ErrorCode := 0;
    NWL := 0;
    NWU := 0;
    JB:=1;
    while JB<=NSPLIT do
    begin
        IOFF := IEND;
        IBEGIN := IOFF+1;
        IEND := ISPLIT[JB];
        IIN := IEND-IOFF;
        if IIN=1 then
        begin
            
            //
            // Special Case -- IIN=1
            //
            if (IRANGE=1) or AP_FP_Greater_Eq(WL,D[IBEGIN]-PIVMIN) then
            begin
                NWL := NWL+1;
            end;
            if (IRANGE=1) or AP_FP_Greater_Eq(WU,D[IBEGIN]-PIVMIN) then
            begin
                NWU := NWU+1;
            end;
            if (IRANGE=1) or AP_FP_Less(WL,D[IBEGIN]-PIVMIN) and AP_FP_Greater_Eq(WU,D[IBEGIN]-PIVMIN) then
            begin
                M := M+1;
                W[M] := D[IBEGIN];
                IBLOCK[M] := JB;
            end;
        end
        else
        begin
            
            //
            // General Case -- IIN > 1
            //
            // Compute Gershgorin Interval
            // and use it as the initial interval
            //
            GU := D[IBEGIN];
            GL := D[IBEGIN];
            TMP1 := 0;
            J:=IBEGIN;
            while J<=IEND-1 do
            begin
                TMP2 := ABSReal(E[J]);
                GU := Max(GU, D[J]+TMP1+TMP2);
                GL := Min(GL, D[J]-TMP1-TMP2);
                TMP1 := TMP2;
                Inc(J);
            end;
            GU := Max(GU, D[IEND]+TMP1);
            GL := Min(GL, D[IEND]-TMP1);
            BNORM := Max(ABSReal(GL), ABSReal(GU));
            GL := GL-FUDGE*BNORM*ULP*IIN-FUDGE*PIVMIN;
            GU := GU+FUDGE*BNORM*ULP*IIN+FUDGE*PIVMIN;
            
            //
            // Compute ATOLI for the current submatrix
            //
            if AP_FP_Less_Eq(ABSTOL,0) then
            begin
                ATOLI := ULP*Max(ABSReal(GL), ABSReal(GU));
            end
            else
            begin
                ATOLI := ABSTOL;
            end;
            if IRANGE>1 then
            begin
                if AP_FP_Less(GU,WL) then
                begin
                    NWL := NWL+IIN;
                    NWU := NWU+IIN;
                    Inc(JB);
                    Continue;
                end;
                GL := Max(GL, WL);
                GU := Min(GU, WU);
                if AP_FP_Greater_Eq(GL,GU) then
                begin
                    Inc(JB);
                    Continue;
                end;
            end;
            
            //
            // Set Up Initial Interval
            //
            WORK[N+1] := GL;
            WORK[N+IIN+1] := GU;
            
            //
            // Calling DLAEBZ
            //
            // CALL DLAEBZ( 1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,
            //    D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ),
            //    IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IM,
            //    IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
            //
            TmpI:=1;
            while TmpI<=IIN do
            begin
                RA1SIIN[TmpI] := D[IBEGIN-1+TmpI];
                if IBEGIN-1+TmpI<N then
                begin
                    RA2SIIN[TmpI] := E[IBEGIN-1+TmpI];
                end;
                RA3SIIN[TmpI] := WORK[IBEGIN-1+TmpI];
                RA1SIINX2[TmpI,1] := WORK[N+TmpI];
                RA1SIINX2[TmpI,2] := WORK[N+TmpI+IIN];
                RA4SIIN[TmpI] := WORK[N+2*IIN+TmpI];
                RWORKSPACE[TmpI] := W[M+TmpI];
                IWORKSPACE[TmpI] := IBLOCK[M+TmpI];
                IA1SIINX2[TmpI,1] := IWORK[TmpI];
                IA1SIINX2[TmpI,2] := IWORK[TmpI+IIN];
                Inc(TmpI);
            end;
            InternalDLAEBZ(1, 0, IIN, IIN, 1, ATOLI, RTOLI, PIVMIN, RA1SIIN, RA2SIIN, RA3SIIN, IDUMMA, RA1SIINX2, RA4SIIN, IM, IA1SIINX2, RWORKSPACE, IWORKSPACE, IINFO);
            TmpI:=1;
            while TmpI<=IIN do
            begin
                WORK[N+TmpI] := RA1SIINX2[TmpI,1];
                WORK[N+TmpI+IIN] := RA1SIINX2[TmpI,2];
                WORK[N+2*IIN+TmpI] := RA4SIIN[TmpI];
                W[M+TmpI] := RWORKSPACE[TmpI];
                IBLOCK[M+TmpI] := IWORKSPACE[TmpI];
                IWORK[TmpI] := IA1SIINX2[TmpI,1];
                IWORK[TmpI+IIN] := IA1SIINX2[TmpI,2];
                Inc(TmpI);
            end;
            NWL := NWL+IWORK[1];
            NWU := NWU+IWORK[IIN+1];
            IWOFF := M-IWORK[1];
            
            //
            // Compute Eigenvalues
            //
            ITMAX := Ceil((Ln(GU-GL+PIVMIN)-Ln(PIVMIN))/Ln(2))+2;
            
            //
            // Calling DLAEBZ
            //
            //CALL DLAEBZ( 2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,
            //    D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ),
            //    IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IOUT,
            //    IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
            //
            TmpI:=1;
            while TmpI<=IIN do
            begin
                RA1SIIN[TmpI] := D[IBEGIN-1+TmpI];
                if IBEGIN-1+TmpI<N then
                begin
                    RA2SIIN[TmpI] := E[IBEGIN-1+TmpI];
                end;
                RA3SIIN[TmpI] := WORK[IBEGIN-1+TmpI];
                RA1SIINX2[TmpI,1] := WORK[N+TmpI];
                RA1SIINX2[TmpI,2] := WORK[N+TmpI+IIN];
                RA4SIIN[TmpI] := WORK[N+2*IIN+TmpI];
                RWORKSPACE[TmpI] := W[M+TmpI];
                IWORKSPACE[TmpI] := IBLOCK[M+TmpI];
                IA1SIINX2[TmpI,1] := IWORK[TmpI];
                IA1SIINX2[TmpI,2] := IWORK[TmpI+IIN];
                Inc(TmpI);
            end;
            InternalDLAEBZ(2, ITMAX, IIN, IIN, 1, ATOLI, RTOLI, PIVMIN, RA1SIIN, RA2SIIN, RA3SIIN, IDUMMA, RA1SIINX2, RA4SIIN, IOUT, IA1SIINX2, RWORKSPACE, IWORKSPACE, IINFO);
            TmpI:=1;
            while TmpI<=IIN do
            begin
                WORK[N+TmpI] := RA1SIINX2[TmpI,1];
                WORK[N+TmpI+IIN] := RA1SIINX2[TmpI,2];
                WORK[N+2*IIN+TmpI] := RA4SIIN[TmpI];
                W[M+TmpI] := RWORKSPACE[TmpI];
                IBLOCK[M+TmpI] := IWORKSPACE[TmpI];
                IWORK[TmpI] := IA1SIINX2[TmpI,1];
                IWORK[TmpI+IIN] := IA1SIINX2[TmpI,2];
                Inc(TmpI);
            end;
            
            //
            // Copy Eigenvalues Into W and IBLOCK
            // Use -JB for block number for unconverged eigenvalues.
            //
            J:=1;
            while J<=IOUT do
            begin
                TMP1 := 0.5*(WORK[J+N]+WORK[J+IIN+N]);
                
                //
                // Flag non-convergence.
                //
                if J>IOUT-IINFO then
                begin
                    NCNVRG := True;
                    IB := -JB;
                end
                else
                begin
                    IB := JB;
                end;
                JE:=IWORK[J]+1+IWOFF;
                while JE<=IWORK[J+IIN]+IWOFF do
                begin
                    W[JE] := TMP1;
                    IBLOCK[JE] := IB;
                    Inc(JE);
                end;
                Inc(J);
            end;
            M := M+IM;
        end;
        Inc(JB);
    end;
    
    //
    // If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
    // If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
    //
    if IRANGE=3 then
    begin
        IM := 0;
        IDISCL := IL-1-NWL;
        IDISCU := NWU-IU;
        if (IDISCL>0) or (IDISCU>0) then
        begin
            JE:=1;
            while JE<=M do
            begin
                if AP_FP_Less_Eq(W[JE],WLU) and (IDISCL>0) then
                begin
                    IDISCL := IDISCL-1;
                end
                else
                begin
                    if AP_FP_Greater_Eq(W[JE],WUL) and (IDISCU>0) then
                    begin
                        IDISCU := IDISCU-1;
                    end
                    else
                    begin
                        IM := IM+1;
                        W[IM] := W[JE];
                        IBLOCK[IM] := IBLOCK[JE];
                    end;
                end;
                Inc(JE);
            end;
            M := IM;
        end;
        if (IDISCL>0) or (IDISCU>0) then
        begin
            
            //
            // Code to deal with effects of bad arithmetic:
            // Some low eigenvalues to be discarded are not in (WL,WLU],
            // or high eigenvalues to be discarded are not in (WUL,WU]
            // so just kill off the smallest IDISCL/largest IDISCU
            // eigenvalues, by simply finding the smallest/largest
            // eigenvalue(s).
            //
            // (If N(w) is monotone non-decreasing, this should never
            //  happen.)
            //
            if IDISCL>0 then
            begin
                WKILL := WU;
                JDISC:=1;
                while JDISC<=IDISCL do
                begin
                    IW := 0;
                    JE:=1;
                    while JE<=M do
                    begin
                        if (IBLOCK[JE]<>0) and (AP_FP_Less(W[JE],WKILL) or (IW=0)) then
                        begin
                            IW := JE;
                            WKILL := W[JE];
                        end;
                        Inc(JE);
                    end;
                    IBLOCK[IW] := 0;
                    Inc(JDISC);
                end;
            end;
            if IDISCU>0 then
            begin
                WKILL := WL;
                JDISC:=1;
                while JDISC<=IDISCU do
                begin
                    IW := 0;
                    JE:=1;
                    while JE<=M do
                    begin
                        if (IBLOCK[JE]<>0) and (AP_FP_Greater(W[JE],WKILL) or (IW=0)) then
                        begin
                            IW := JE;
                            WKILL := W[JE];
                        end;
                        Inc(JE);
                    end;
                    IBLOCK[IW] := 0;
                    Inc(JDISC);
                end;
            end;
            IM := 0;
            JE:=1;
            while JE<=M do
            begin
                if IBLOCK[JE]<>0 then
                begin
                    IM := IM+1;
                    W[IM] := W[JE];
                    IBLOCK[IM] := IBLOCK[JE];
                end;
                Inc(JE);
            end;
            M := IM;
        end;
        if (IDISCL<0) or (IDISCU<0) then
        begin
            TOOFEW := True;
        end;
    end;
    
    //
    // If ORDER='B', do nothing -- the eigenvalues are already sorted
    //    by block.
    // If ORDER='E', sort the eigenvalues from smallest to largest
    //
    if (IORDER=1) and (NSPLIT>1) then
    begin
        JE:=1;
        while JE<=M-1 do
        begin
            IE := 0;
            TMP1 := W[JE];
            J:=JE+1;
            while J<=M do
            begin
                if AP_FP_Less(W[J],TMP1) then
                begin
                    IE := J;
                    TMP1 := W[J];
                end;
                Inc(J);
            end;
            if IE<>0 then
            begin
                ITMP1 := IBLOCK[IE];
                W[IE] := W[JE];
                IBLOCK[IE] := IBLOCK[JE];
                W[JE] := TMP1;
                IBLOCK[JE] := ITMP1;
            end;
            Inc(JE);
        end;
    end;
    J:=1;
    while J<=M do
    begin
        W[J] := W[J]*ScaleFactor;
        Inc(J);
    end;
    ErrorCode := 0;
    if NCNVRG then
    begin
        ErrorCode := ErrorCode+1;
    end;
    if TOOFEW then
    begin
        ErrorCode := ErrorCode+2;
    end;
    Result := ErrorCode=0;
end;


procedure InternalDSTEIN(const N : AlglibInteger;
     const D : TReal1DArray;
     E : TReal1DArray;
     const M : AlglibInteger;
     W : TReal1DArray;
     const IBLOCK : TInteger1DArray;
     const ISPLIT : TInteger1DArray;
     var Z : TReal2DArray;
     var IFAIL : TInteger1DArray;
     var INFO : AlglibInteger);
var
    MAXITS : AlglibInteger;
    EXTRA : AlglibInteger;
    B1 : AlglibInteger;
    BLKSIZ : AlglibInteger;
    BN : AlglibInteger;
    GPIND : AlglibInteger;
    I : AlglibInteger;
    IINFO : AlglibInteger;
    ITS : AlglibInteger;
    J : AlglibInteger;
    J1 : AlglibInteger;
    JBLK : AlglibInteger;
    JMAX : AlglibInteger;
    NBLK : AlglibInteger;
    NRMCHK : AlglibInteger;
    DTPCRT : Double;
    EPS : Double;
    EPS1 : Double;
    NRM : Double;
    ONENRM : Double;
    ORTOL : Double;
    PERTOL : Double;
    SCL : Double;
    SEP : Double;
    TOL : Double;
    XJ : Double;
    XJM : Double;
    ZTR : Double;
    WORK1 : TReal1DArray;
    WORK2 : TReal1DArray;
    WORK3 : TReal1DArray;
    WORK4 : TReal1DArray;
    WORK5 : TReal1DArray;
    IWORK : TInteger1DArray;
    TmpCriterion : Boolean;
    TI : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    E := DynamicArrayCopy(E);
    W := DynamicArrayCopy(W);
    MAXITS := 5;
    EXTRA := 2;
    SetLength(WORK1, Max(N, 1)+1);
    SetLength(WORK2, Max(N-1, 1)+1);
    SetLength(WORK3, Max(N, 1)+1);
    SetLength(WORK4, Max(N, 1)+1);
    SetLength(WORK5, Max(N, 1)+1);
    SetLength(IWORK, Max(N, 1)+1);
    SetLength(IFAIL, Max(M, 1)+1);
    SetLength(Z, Max(N, 1)+1, Max(M, 1)+1);
    
    //
    // Test the input parameters.
    //
    INFO := 0;
    I:=1;
    while I<=M do
    begin
        IFAIL[I] := 0;
        Inc(I);
    end;
    if N<0 then
    begin
        INFO := -1;
        Exit;
    end;
    if (M<0) or (M>N) then
    begin
        INFO := -4;
        Exit;
    end;
    J:=2;
    while J<=M do
    begin
        if IBLOCK[J]<IBLOCK[J-1] then
        begin
            INFO := -6;
            Break;
        end;
        if (IBLOCK[J]=IBLOCK[J-1]) and AP_FP_Less(W[J],W[J-1]) then
        begin
            INFO := -5;
            Break;
        end;
        Inc(J);
    end;
    if INFO<>0 then
    begin
        Exit;
    end;
    
    //
    // Quick return if possible
    //
    if (N=0) or (M=0) then
    begin
        Exit;
    end;
    if N=1 then
    begin
        Z[1,1] := 1;
        Exit;
    end;
    
    //
    // Some preparations
    //
    TI := N-1;
    APVMove(@WORK1[0], 1, TI, @E[0], 1, TI);
    SetLength(E, N+1);
    APVMove(@E[0], 1, TI, @WORK1[0], 1, TI);
    APVMove(@WORK1[0], 1, M, @W[0], 1, M);
    SetLength(W, N+1);
    APVMove(@W[0], 1, M, @WORK1[0], 1, M);
    
    //
    // Get machine constants.
    //
    EPS := MachineEpsilon;
    
    //
    // Compute eigenvectors of matrix blocks.
    //
    J1 := 1;
    NBLK:=1;
    while NBLK<=IBLOCK[M] do
    begin
        
        //
        // Find starting and ending indices of block nblk.
        //
        if NBLK=1 then
        begin
            B1 := 1;
        end
        else
        begin
            B1 := ISPLIT[NBLK-1]+1;
        end;
        BN := ISPLIT[NBLK];
        BLKSIZ := BN-B1+1;
        if BLKSIZ<>1 then
        begin
            
            //
            // Compute reorthogonalization criterion and stopping criterion.
            //
            GPIND := B1;
            ONENRM := ABSReal(D[B1])+ABSReal(E[B1]);
            ONENRM := Max(ONENRM, ABSReal(D[BN])+ABSReal(E[BN-1]));
            I:=B1+1;
            while I<=BN-1 do
            begin
                ONENRM := Max(ONENRM, ABSReal(D[I])+ABSReal(E[I-1])+ABSReal(E[I]));
                Inc(I);
            end;
            ORTOL := 0.001*ONENRM;
            DTPCRT := SQRT(0.1/BLKSIZ);
        end;
        
        //
        // Loop through eigenvalues of block nblk.
        //
        JBLK := 0;
        J:=J1;
        while J<=M do
        begin
            if IBLOCK[J]<>NBLK then
            begin
                J1 := J;
                Break;
            end;
            JBLK := JBLK+1;
            XJ := W[J];
            if BLKSIZ=1 then
            begin
                
                //
                // Skip all the work if the block size is one.
                //
                WORK1[1] := 1;
            end
            else
            begin
                
                //
                // If eigenvalues j and j-1 are too close, add a relatively
                // small perturbation.
                //
                if JBLK>1 then
                begin
                    EPS1 := ABSReal(EPS*XJ);
                    PERTOL := 10*EPS1;
                    SEP := XJ-XJM;
                    if AP_FP_Less(SEP,PERTOL) then
                    begin
                        XJ := XJM+PERTOL;
                    end;
                end;
                ITS := 0;
                NRMCHK := 0;
                
                //
                // Get random starting vector.
                //
                TI:=1;
                while TI<=BLKSIZ do
                begin
                    WORK1[TI] := 2*RandomReal-1;
                    Inc(TI);
                end;
                
                //
                // Copy the matrix T so it won't be destroyed in factorization.
                //
                TI:=1;
                while TI<=BLKSIZ-1 do
                begin
                    WORK2[TI] := E[B1+TI-1];
                    WORK3[TI] := E[B1+TI-1];
                    WORK4[TI] := D[B1+TI-1];
                    Inc(TI);
                end;
                WORK4[BLKSIZ] := D[B1+BLKSIZ-1];
                
                //
                // Compute LU factors with partial pivoting  ( PT = LU )
                //
                TOL := 0;
                TDINInternalDLAGTF(BLKSIZ, WORK4, XJ, WORK2, WORK3, TOL, WORK5, IWORK, IINFO);
                
                //
                // Update iteration count.
                //
                repeat
                    ITS := ITS+1;
                    if ITS>MAXITS then
                    begin
                        
                        //
                        // If stopping criterion was not satisfied, update info and
                        // store eigenvector number in array ifail.
                        //
                        INFO := INFO+1;
                        IFAIL[INFO] := J;
                        Break;
                    end;
                    
                    //
                    // Normalize and scale the righthand side vector Pb.
                    //
                    V := 0;
                    TI:=1;
                    while TI<=BLKSIZ do
                    begin
                        V := V+AbsReal(WORK1[TI]);
                        Inc(TI);
                    end;
                    SCL := BLKSIZ*ONENRM*Max(EPS, ABSReal(WORK4[BLKSIZ]))/V;
                    APVMul(@WORK1[0], 1, BLKSIZ, SCL);
                    
                    //
                    // Solve the system LU = Pb.
                    //
                    TDINInternalDLAGTS(BLKSIZ, WORK4, WORK2, WORK3, WORK5, IWORK, WORK1, TOL, IINFO);
                    
                    //
                    // Reorthogonalize by modified Gram-Schmidt if eigenvalues are
                    // close enough.
                    //
                    if JBLK<>1 then
                    begin
                        if AP_FP_Greater(ABSReal(XJ-XJM),ORTOL) then
                        begin
                            GPIND := J;
                        end;
                        if GPIND<>J then
                        begin
                            I:=GPIND;
                            while I<=J-1 do
                            begin
                                I1 := B1;
                                I2 := B1+BLKSIZ-1;
                                i1_ := (I1)-(1);
                                ZTR := 0.0;
                                for i_ := 1 to BLKSIZ do
                                begin
                                    ZTR := ZTR + WORK1[i_]*Z[i_+i1_,I];
                                end;
                                i1_ := (I1) - (1);
                                for i_ := 1 to BLKSIZ do
                                begin
                                    WORK1[i_] := WORK1[i_] - ZTR*Z[i_+i1_,I];
                                end;
                                Inc(I);
                            end;
                        end;
                    end;
                    
                    //
                    // Check the infinity norm of the iterate.
                    //
                    JMAX := VectorIdxAbsMax(WORK1, 1, BLKSIZ);
                    NRM := AbsReal(WORK1[JMAX]);
                    
                    //
                    // Continue for additional iterations after norm reaches
                    // stopping criterion.
                    //
                    TmpCriterion := False;
                    if AP_FP_Less(NRM,DTPCRT) then
                    begin
                        TmpCriterion := True;
                    end
                    else
                    begin
                        NRMCHK := NRMCHK+1;
                        if NRMCHK<EXTRA+1 then
                        begin
                            TmpCriterion := True;
                        end;
                    end;
                until  not TmpCriterion;
                
                //
                // Accept iterate as jth eigenvector.
                //
                SCL := 1/VectorNorm2(WORK1, 1, BLKSIZ);
                JMAX := VectorIdxAbsMax(WORK1, 1, BLKSIZ);
                if AP_FP_Less(WORK1[JMAX],0) then
                begin
                    SCL := -SCL;
                end;
                APVMul(@WORK1[0], 1, BLKSIZ, SCL);
            end;
            I:=1;
            while I<=N do
            begin
                Z[I,J] := 0;
                Inc(I);
            end;
            I:=1;
            while I<=BLKSIZ do
            begin
                Z[B1+I-1,J] := WORK1[I];
                Inc(I);
            end;
            
            //
            // Save the shift to check eigenvalue spacing at next
            // iteration.
            //
            XJM := XJ;
            Inc(J);
        end;
        Inc(NBLK);
    end;
end;


procedure TDINInternalDLAGTF(const N : AlglibInteger;
     var A : TReal1DArray;
     const LAMBDA : Double;
     var B : TReal1DArray;
     var C : TReal1DArray;
     const TOL : Double;
     var D : TReal1DArray;
     var IIN : TInteger1DArray;
     var INFO : AlglibInteger);
var
    K : AlglibInteger;
    EPS : Double;
    MULT : Double;
    PIV1 : Double;
    PIV2 : Double;
    SCALE1 : Double;
    SCALE2 : Double;
    TEMP : Double;
    TL : Double;
begin
    INFO := 0;
    if N<0 then
    begin
        INFO := -1;
        Exit;
    end;
    if N=0 then
    begin
        Exit;
    end;
    A[1] := A[1]-LAMBDA;
    IIN[N] := 0;
    if N=1 then
    begin
        if AP_FP_Eq(A[1],0) then
        begin
            IIN[1] := 1;
        end;
        Exit;
    end;
    EPS := MachineEpsilon;
    TL := Max(TOL, EPS);
    SCALE1 := ABSReal(A[1])+ABSReal(B[1]);
    K:=1;
    while K<=N-1 do
    begin
        A[K+1] := A[K+1]-LAMBDA;
        SCALE2 := ABSReal(C[K])+ABSReal(A[K+1]);
        if K<N-1 then
        begin
            SCALE2 := SCALE2+ABSReal(B[K+1]);
        end;
        if AP_FP_Eq(A[K],0) then
        begin
            PIV1 := 0;
        end
        else
        begin
            PIV1 := ABSReal(A[K])/SCALE1;
        end;
        if AP_FP_Eq(C[K],0) then
        begin
            IIN[K] := 0;
            PIV2 := 0;
            SCALE1 := SCALE2;
            if K<N-1 then
            begin
                D[K] := 0;
            end;
        end
        else
        begin
            PIV2 := ABSReal(C[K])/SCALE2;
            if AP_FP_Less_Eq(PIV2,PIV1) then
            begin
                IIN[K] := 0;
                SCALE1 := SCALE2;
                C[K] := C[K]/A[K];
                A[K+1] := A[K+1]-C[K]*B[K];
                if K<N-1 then
                begin
                    D[K] := 0;
                end;
            end
            else
            begin
                IIN[K] := 1;
                MULT := A[K]/C[K];
                A[K] := C[K];
                TEMP := A[K+1];
                A[K+1] := B[K]-MULT*TEMP;
                if K<N-1 then
                begin
                    D[K] := B[K+1];
                    B[K+1] := -MULT*D[K];
                end;
                B[K] := TEMP;
                C[K] := MULT;
            end;
        end;
        if AP_FP_Less_Eq(Max(PIV1, PIV2),TL) and (IIN[N]=0) then
        begin
            IIN[N] := K;
        end;
        Inc(K);
    end;
    if AP_FP_Less_Eq(ABSReal(A[N]),SCALE1*TL) and (IIN[N]=0) then
    begin
        IIN[N] := N;
    end;
end;


procedure TDINInternalDLAGTS(const N : AlglibInteger;
     const A : TReal1DArray;
     const B : TReal1DArray;
     const C : TReal1DArray;
     const D : TReal1DArray;
     const IIN : TInteger1DArray;
     var Y : TReal1DArray;
     var TOL : Double;
     var INFO : AlglibInteger);
var
    K : AlglibInteger;
    ABSAK : Double;
    AK : Double;
    BIGNUM : Double;
    EPS : Double;
    PERT : Double;
    SFMIN : Double;
    TEMP : Double;
begin
    INFO := 0;
    if N<0 then
    begin
        INFO := -1;
        Exit;
    end;
    if N=0 then
    begin
        Exit;
    end;
    EPS := MachineEpsilon;
    SFMIN := MinRealNumber;
    BIGNUM := 1/SFMIN;
    if AP_FP_Less_Eq(TOL,0) then
    begin
        TOL := ABSReal(A[1]);
        if N>1 then
        begin
            TOL := Max(TOL, Max(ABSReal(A[2]), ABSReal(B[1])));
        end;
        K:=3;
        while K<=N do
        begin
            TOL := Max(TOL, Max(ABSReal(A[K]), Max(ABSReal(B[K-1]), ABSReal(D[K-2]))));
            Inc(K);
        end;
        TOL := TOL*EPS;
        if AP_FP_Eq(TOL,0) then
        begin
            TOL := EPS;
        end;
    end;
    K:=2;
    while K<=N do
    begin
        if IIN[K-1]=0 then
        begin
            Y[K] := Y[K]-C[K-1]*Y[K-1];
        end
        else
        begin
            TEMP := Y[K-1];
            Y[K-1] := Y[K];
            Y[K] := TEMP-C[K-1]*Y[K];
        end;
        Inc(K);
    end;
    K:=N;
    while K>=1 do
    begin
        if K<=N-2 then
        begin
            TEMP := Y[K]-B[K]*Y[K+1]-D[K]*Y[K+2];
        end
        else
        begin
            if K=N-1 then
            begin
                TEMP := Y[K]-B[K]*Y[K+1];
            end
            else
            begin
                TEMP := Y[K];
            end;
        end;
        AK := A[K];
        PERT := AbsReal(TOL);
        if AP_FP_Less(AK,0) then
        begin
            PERT := -PERT;
        end;
        while True do
        begin
            ABSAK := ABSReal(AK);
            if AP_FP_Less(ABSAK,1) then
            begin
                if AP_FP_Less(ABSAK,SFMIN) then
                begin
                    if AP_FP_Eq(ABSAK,0) or AP_FP_Greater(ABSReal(TEMP)*SFMIN,ABSAK) then
                    begin
                        AK := AK+PERT;
                        PERT := 2*PERT;
                        Continue;
                    end
                    else
                    begin
                        TEMP := TEMP*BIGNUM;
                        AK := AK*BIGNUM;
                    end;
                end
                else
                begin
                    if AP_FP_Greater(ABSReal(TEMP),ABSAK*BIGNUM) then
                    begin
                        AK := AK+PERT;
                        PERT := 2*PERT;
                        Continue;
                    end;
                end;
            end;
            Break;
        end;
        Y[K] := TEMP/AK;
        Dec(K);
    end;
end;


procedure InternalDLAEBZ(const IJOB : AlglibInteger;
     const NITMAX : AlglibInteger;
     const N : AlglibInteger;
     const MMAX : AlglibInteger;
     const MINP : AlglibInteger;
     const ABSTOL : Double;
     const RELTOL : Double;
     const PIVMIN : Double;
     const D : TReal1DArray;
     const E : TReal1DArray;
     const E2 : TReal1DArray;
     var NVAL : TInteger1DArray;
     var AB : TReal2DArray;
     var C : TReal1DArray;
     var MOUT : AlglibInteger;
     var NAB : TInteger2DArray;
     var WORK : TReal1DArray;
     var IWORK : TInteger1DArray;
     var INFO : AlglibInteger);
var
    ITMP1 : AlglibInteger;
    ITMP2 : AlglibInteger;
    J : AlglibInteger;
    JI : AlglibInteger;
    JIT : AlglibInteger;
    JP : AlglibInteger;
    KF : AlglibInteger;
    KFNEW : AlglibInteger;
    KL : AlglibInteger;
    KLNEW : AlglibInteger;
    TMP1 : Double;
    TMP2 : Double;
begin
    INFO := 0;
    if (IJOB<1) or (IJOB>3) then
    begin
        INFO := -1;
        Exit;
    end;
    
    //
    // Initialize NAB
    //
    if IJOB=1 then
    begin
        
        //
        // Compute the number of eigenvalues in the initial intervals.
        //
        MOUT := 0;
        
        //
        //DIR$ NOVECTOR
        //
        JI:=1;
        while JI<=MINP do
        begin
            JP:=1;
            while JP<=2 do
            begin
                TMP1 := D[1]-AB[JI,JP];
                if AP_FP_Less(ABSReal(TMP1),PIVMIN) then
                begin
                    TMP1 := -PIVMIN;
                end;
                NAB[JI,JP] := 0;
                if AP_FP_Less_Eq(TMP1,0) then
                begin
                    NAB[JI,JP] := 1;
                end;
                J:=2;
                while J<=N do
                begin
                    TMP1 := D[J]-E2[J-1]/TMP1-AB[JI,JP];
                    if AP_FP_Less(ABSReal(TMP1),PIVMIN) then
                    begin
                        TMP1 := -PIVMIN;
                    end;
                    if AP_FP_Less_Eq(TMP1,0) then
                    begin
                        NAB[JI,JP] := NAB[JI,JP]+1;
                    end;
                    Inc(J);
                end;
                Inc(JP);
            end;
            MOUT := MOUT+NAB[JI,2]-NAB[JI,1];
            Inc(JI);
        end;
        Exit;
    end;
    
    //
    // Initialize for loop
    //
    // KF and KL have the following meaning:
    //   Intervals 1,...,KF-1 have converged.
    //   Intervals KF,...,KL  still need to be refined.
    //
    KF := 1;
    KL := MINP;
    
    //
    // If IJOB=2, initialize C.
    // If IJOB=3, use the user-supplied starting point.
    //
    if IJOB=2 then
    begin
        JI:=1;
        while JI<=MINP do
        begin
            C[JI] := 0.5*(AB[JI,1]+AB[JI,2]);
            Inc(JI);
        end;
    end;
    
    //
    // Iteration loop
    //
    JIT:=1;
    while JIT<=NITMAX do
    begin
        
        //
        // Loop over intervals
        //
        //
        // Serial Version of the loop
        //
        KLNEW := KL;
        JI:=KF;
        while JI<=KL do
        begin
            
            //
            // Compute N(w), the number of eigenvalues less than w
            //
            TMP1 := C[JI];
            TMP2 := D[1]-TMP1;
            ITMP1 := 0;
            if AP_FP_Less_Eq(TMP2,PIVMIN) then
            begin
                ITMP1 := 1;
                TMP2 := Min(TMP2, -PIVMIN);
            end;
            
            //
            // A series of compiler directives to defeat vectorization
            // for the next loop
            //
            //*$PL$ CMCHAR=' '
            //CDIR$          NEXTSCALAR
            //C$DIR          SCALAR
            //CDIR$          NEXT SCALAR
            //CVD$L          NOVECTOR
            //CDEC$          NOVECTOR
            //CVD$           NOVECTOR
            //*VDIR          NOVECTOR
            //*VOCL          LOOP,SCALAR
            //CIBM           PREFER SCALAR
            //*$PL$ CMCHAR='*'
            //
            J:=2;
            while J<=N do
            begin
                TMP2 := D[J]-E2[J-1]/TMP2-TMP1;
                if AP_FP_Less_Eq(TMP2,PIVMIN) then
                begin
                    ITMP1 := ITMP1+1;
                    TMP2 := Min(TMP2, -PIVMIN);
                end;
                Inc(J);
            end;
            if IJOB<=2 then
            begin
                
                //
                // IJOB=2: Choose all intervals containing eigenvalues.
                //
                // Insure that N(w) is monotone
                //
                ITMP1 := Min(NAB[JI,2], Max(NAB[JI,1], ITMP1));
                
                //
                // Update the Queue -- add intervals if both halves
                // contain eigenvalues.
                //
                if ITMP1=NAB[JI,2] then
                begin
                    
                    //
                    // No eigenvalue in the upper interval:
                    // just use the lower interval.
                    //
                    AB[JI,2] := TMP1;
                end
                else
                begin
                    if ITMP1=NAB[JI,1] then
                    begin
                        
                        //
                        // No eigenvalue in the lower interval:
                        // just use the upper interval.
                        //
                        AB[JI,1] := TMP1;
                    end
                    else
                    begin
                        if KLNEW<MMAX then
                        begin
                            
                            //
                            // Eigenvalue in both intervals -- add upper to queue.
                            //
                            KLNEW := KLNEW+1;
                            AB[KLNEW,2] := AB[JI,2];
                            NAB[KLNEW,2] := NAB[JI,2];
                            AB[KLNEW,1] := TMP1;
                            NAB[KLNEW,1] := ITMP1;
                            AB[JI,2] := TMP1;
                            NAB[JI,2] := ITMP1;
                        end
                        else
                        begin
                            INFO := MMAX+1;
                            Exit;
                        end;
                    end;
                end;
            end
            else
            begin
                
                //
                // IJOB=3: Binary search.  Keep only the interval
                // containing  w  s.t. N(w) = NVAL
                //
                if ITMP1<=NVAL[JI] then
                begin
                    AB[JI,1] := TMP1;
                    NAB[JI,1] := ITMP1;
                end;
                if ITMP1>=NVAL[JI] then
                begin
                    AB[JI,2] := TMP1;
                    NAB[JI,2] := ITMP1;
                end;
            end;
            Inc(JI);
        end;
        KL := KLNEW;
        
        //
        // Check for convergence
        //
        KFNEW := KF;
        JI:=KF;
        while JI<=KL do
        begin
            TMP1 := ABSReal(AB[JI,2]-AB[JI,1]);
            TMP2 := Max(ABSReal(AB[JI,2]), ABSReal(AB[JI,1]));
            if AP_FP_Less(TMP1,Max(ABSTOL, Max(PIVMIN, RELTOL*TMP2))) or (NAB[JI,1]>=NAB[JI,2]) then
            begin
                
                //
                // Converged -- Swap with position KFNEW,
                // then increment KFNEW
                //
                if JI>KFNEW then
                begin
                    TMP1 := AB[JI,1];
                    TMP2 := AB[JI,2];
                    ITMP1 := NAB[JI,1];
                    ITMP2 := NAB[JI,2];
                    AB[JI,1] := AB[KFNEW,1];
                    AB[JI,2] := AB[KFNEW,2];
                    NAB[JI,1] := NAB[KFNEW,1];
                    NAB[JI,2] := NAB[KFNEW,2];
                    AB[KFNEW,1] := TMP1;
                    AB[KFNEW,2] := TMP2;
                    NAB[KFNEW,1] := ITMP1;
                    NAB[KFNEW,2] := ITMP2;
                    if IJOB=3 then
                    begin
                        ITMP1 := NVAL[JI];
                        NVAL[JI] := NVAL[KFNEW];
                        NVAL[KFNEW] := ITMP1;
                    end;
                end;
                KFNEW := KFNEW+1;
            end;
            Inc(JI);
        end;
        KF := KFNEW;
        
        //
        // Choose Midpoints
        //
        JI:=KF;
        while JI<=KL do
        begin
            C[JI] := 0.5*(AB[JI,1]+AB[JI,2]);
            Inc(JI);
        end;
        
        //
        // If no more intervals to refine, quit.
        //
        if KF>KL then
        begin
            Break;
        end;
        Inc(JIT);
    end;
    
    //
    // Converged
    //
    INFO := Max(KL+1-KF, 0);
    MOUT := KL;
end;


end.