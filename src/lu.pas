(*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee. All rights reserved.

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
unit lu;
interface
uses Math, Sysutils, Ap;

procedure RMatrixLU(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
procedure LUDecomposition(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
procedure LUDecompositionUnpacked(A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var L : TReal2DArray;
     var U : TReal2DArray;
     var Pivots : TInteger1DArray);

implementation

const
    LUNB = 8;

procedure RMatrixLU2(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);forward;


(*************************************************************************
LU decomposition of a general matrix of size MxN

The subroutine calculates the LU decomposition of a rectangular general
matrix with partial pivoting (with row permutations).

Input parameters:
    A   -   matrix A whose indexes range within [0..M-1, 0..N-1].
    M   -   number of rows in matrix A.
    N   -   number of columns in matrix A.

Output parameters:
    A   -   matrices L and U in compact form (see below).
            Array whose indexes range within [0..M-1, 0..N-1].
    Pivots - permutation matrix in compact form (see below).
            Array whose index ranges within [0..Min(M-1,N-1)].

Matrix A is represented as A = P * L * U, where P is a permutation matrix,
matrix L - lower triangular (or lower trapezoid, if M>N) matrix,
U - upper triangular (or upper trapezoid, if M<N) matrix.

Let M be equal to 4 and N be equal to 3:

                   (  1          )    ( U11 U12 U13  )
A = P1 * P2 * P3 * ( L21  1      )  * (     U22 U23  )
                   ( L31 L32  1  )    (         U33  )
                   ( L41 L42 L43 )

Matrix L has size MxMin(M,N), matrix U has size Min(M,N)xN, matrix P(i) is
a permutation of the identity matrix of size MxM with numbers I and Pivots[I].

The algorithm returns array Pivots and the following matrix which replaces
matrix A and contains matrices L and U in compact form (the example applies
to M=4, N=3).

 ( U11 U12 U13 )
 ( L21 U22 U23 )
 ( L31 L32 U33 )
 ( L41 L42 L43 )

As we can see, the unit diagonal isn't stored.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1992
*************************************************************************)
procedure RMatrixLU(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
var
    B : TReal2DArray;
    T : TReal1DArray;
    BP : TInteger1DArray;
    MinMN : AlglibInteger;
    I : AlglibInteger;
    IP : AlglibInteger;
    J : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
    CB : AlglibInteger;
    NB : AlglibInteger;
    V : Double;
begin
    Assert(LUNB>=1, 'RMatrixLU internal error');
    NB := LUNB;
    
    //
    // Decide what to use - blocked or unblocked code
    //
    if (N<=1) or (Min(M, N)<=NB) or (NB=1) then
    begin
        
        //
        // Unblocked code
        //
        RMatrixLU2(A, M, N, Pivots);
    end
    else
    begin
        
        //
        // Blocked code.
        // First, prepare temporary matrix and indices
        //
        SetLength(B, M-1+1, NB-1+1);
        SetLength(T, N-1+1);
        SetLength(Pivots, Min(M, N)-1+1);
        MinMN := Min(M, N);
        J1 := 0;
        J2 := Min(MinMN, NB)-1;
        
        //
        // Main cycle
        //
        while J1<MinMN do
        begin
            CB := J2-J1+1;
            
            //
            // LU factorization of diagonal and subdiagonal blocks:
            // 1. Copy columns J1..J2 of A to B
            // 2. LU(B)
            // 3. Copy result back to A
            // 4. Copy pivots, apply pivots
            //
            I:=J1;
            while I<=M-1 do
            begin
                APVMove(@B[I-J1][0], 0, CB-1, @A[I][0], J1, J2);
                Inc(I);
            end;
            RMatrixLU2(B, M-J1, CB, BP);
            I:=J1;
            while I<=M-1 do
            begin
                APVMove(@A[I][0], J1, J2, @B[I-J1][0], 0, CB-1);
                Inc(I);
            end;
            I:=0;
            while I<=CB-1 do
            begin
                IP := BP[I];
                Pivots[J1+I] := J1+IP;
                if BP[I]<>I then
                begin
                    if J1<>0 then
                    begin
                        
                        //
                        // Interchange columns 0:J1-1
                        //
                        APVMove(@T[0], 0, J1-1, @A[J1+I][0], 0, J1-1);
                        APVMove(@A[J1+I][0], 0, J1-1, @A[J1+IP][0], 0, J1-1);
                        APVMove(@A[J1+IP][0], 0, J1-1, @T[0], 0, J1-1);
                    end;
                    if J2<N-1 then
                    begin
                        
                        //
                        // Interchange the rest of the matrix, if needed
                        //
                        APVMove(@T[0], J2+1, N-1, @A[J1+I][0], J2+1, N-1);
                        APVMove(@A[J1+I][0], J2+1, N-1, @A[J1+IP][0], J2+1, N-1);
                        APVMove(@A[J1+IP][0], J2+1, N-1, @T[0], J2+1, N-1);
                    end;
                end;
                Inc(I);
            end;
            
            //
            // Compute block row of U
            //
            if J2<N-1 then
            begin
                I:=J1+1;
                while I<=J2 do
                begin
                    J:=J1;
                    while J<=I-1 do
                    begin
                        V := A[I,J];
                        APVSub(@A[I][0], J2+1, N-1, @A[J][0], J2+1, N-1, V);
                        Inc(J);
                    end;
                    Inc(I);
                end;
            end;
            
            //
            // Update trailing submatrix
            //
            if J2<N-1 then
            begin
                I:=J2+1;
                while I<=M-1 do
                begin
                    J:=J1;
                    while J<=J2 do
                    begin
                        V := A[I,J];
                        APVSub(@A[I][0], J2+1, N-1, @A[J][0], J2+1, N-1, V);
                        Inc(J);
                    end;
                    Inc(I);
                end;
            end;
            
            //
            // Next step
            //
            J1 := J2+1;
            J2 := Min(MinMN, J1+NB)-1;
        end;
    end;
end;


procedure LUDecomposition(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    JP : AlglibInteger;
    T1 : TReal1DArray;
    s : Double;
    i_ : AlglibInteger;
begin
    SetLength(Pivots, Min(M, N)+1);
    SetLength(T1, Max(M, N)+1);
    Assert((M>=0) and (N>=0), 'Error in LUDecomposition: incorrect function arguments');
    
    //
    // Quick return if possible
    //
    if (M=0) or (N=0) then
    begin
        Exit;
    end;
    J:=1;
    while J<=Min(M, N) do
    begin
        
        //
        // Find pivot and test for singularity.
        //
        JP := J;
        I:=J+1;
        while I<=M do
        begin
            if AP_FP_Greater(AbsReal(A[I,J]),AbsReal(A[JP,J])) then
            begin
                JP := I;
            end;
            Inc(I);
        end;
        Pivots[J] := JP;
        if AP_FP_Neq(A[JP,J],0) then
        begin
            
            //
            //Apply the interchange to rows
            //
            if JP<>J then
            begin
                APVMove(@T1[0], 1, N, @A[J][0], 1, N);
                APVMove(@A[J][0], 1, N, @A[JP][0], 1, N);
                APVMove(@A[JP][0], 1, N, @T1[0], 1, N);
            end;
            
            //
            //Compute elements J+1:M of J-th column.
            //
            if J<M then
            begin
                
                //
                // CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
                //
                JP := J+1;
                S := 1/A[J,J];
                for i_ := JP to M do
                begin
                    A[i_,J] := S*A[i_,J];
                end;
            end;
        end;
        if J<Min(M, N) then
        begin
            
            //
            //Update trailing submatrix.
            //CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,A( J+1, J+1 ), LDA )
            //
            JP := J+1;
            I:=J+1;
            while I<=M do
            begin
                S := A[I,J];
                APVSub(@A[I][0], JP, N, @A[J][0], JP, N, S);
                Inc(I);
            end;
        end;
        Inc(J);
    end;
end;


procedure LUDecompositionUnpacked(A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var L : TReal2DArray;
     var U : TReal2DArray;
     var Pivots : TInteger1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    MinMN : AlglibInteger;
begin
    A := DynamicArrayCopy(A);
    if (M=0) or (N=0) then
    begin
        Exit;
    end;
    MinMN := Min(M, N);
    SetLength(L, M+1, MinMN+1);
    SetLength(U, MinMN+1, N+1);
    LUDecomposition(A, M, N, Pivots);
    I:=1;
    while I<=M do
    begin
        J:=1;
        while J<=MinMN do
        begin
            if J>I then
            begin
                L[I,J] := 0;
            end;
            if J=I then
            begin
                L[I,J] := 1;
            end;
            if J<I then
            begin
                L[I,J] := A[I,J];
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=1;
    while I<=MinMN do
    begin
        J:=1;
        while J<=N do
        begin
            if J<I then
            begin
                U[I,J] := 0;
            end;
            if J>=I then
            begin
                U[I,J] := A[I,J];
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Level 2 BLAS version of RMatrixLU

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1992
*************************************************************************)
procedure RMatrixLU2(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    JP : AlglibInteger;
    T1 : TReal1DArray;
    s : Double;
    i_ : AlglibInteger;
begin
    SetLength(Pivots, Min(M-1, N-1)+1);
    SetLength(T1, Max(M-1, N-1)+1);
    Assert((M>=0) and (N>=0), 'Error in LUDecomposition: incorrect function arguments');
    
    //
    // Quick return if possible
    //
    if (M=0) or (N=0) then
    begin
        Exit;
    end;
    J:=0;
    while J<=Min(M-1, N-1) do
    begin
        
        //
        // Find pivot and test for singularity.
        //
        JP := J;
        I:=J+1;
        while I<=M-1 do
        begin
            if AP_FP_Greater(AbsReal(A[I,J]),AbsReal(A[JP,J])) then
            begin
                JP := I;
            end;
            Inc(I);
        end;
        Pivots[J] := JP;
        if AP_FP_Neq(A[JP,J],0) then
        begin
            
            //
            //Apply the interchange to rows
            //
            if JP<>J then
            begin
                APVMove(@T1[0], 0, N-1, @A[J][0], 0, N-1);
                APVMove(@A[J][0], 0, N-1, @A[JP][0], 0, N-1);
                APVMove(@A[JP][0], 0, N-1, @T1[0], 0, N-1);
            end;
            
            //
            //Compute elements J+1:M of J-th column.
            //
            if J<M then
            begin
                JP := J+1;
                S := 1/A[J,J];
                for i_ := JP to M-1 do
                begin
                    A[i_,J] := S*A[i_,J];
                end;
            end;
        end;
        if J<Min(M, N)-1 then
        begin
            
            //
            //Update trailing submatrix.
            //
            JP := J+1;
            I:=J+1;
            while I<=M-1 do
            begin
                S := A[I,J];
                APVSub(@A[I][0], JP, N-1, @A[J][0], JP, N-1, S);
                Inc(I);
            end;
        end;
        Inc(J);
    end;
end;


end.