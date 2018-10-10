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
unit spdinverse;
interface
uses Math, Sysutils, Ap, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac;

function SPDMatrixCholeskyInverse(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
function SPDMatrixInverse(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;

implementation

(*************************************************************************
Inversion of a symmetric positive definite matrix which is given
by Cholesky decomposition.

Input parameters:
    A       -   Cholesky decomposition of the matrix to be inverted:
                A=U’*U or A = L*L'.
                Output of  CholeskyDecomposition subroutine.
                Array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper –   storage format.
                If IsUpper = True, then matrix A is given as A = U'*U
                (matrix contains upper triangle).
                Similarly, if IsUpper = False, then A = L*L'.

Output parameters:
    A       -   upper or lower triangle of symmetric matrix A^-1, depending
                on the value of IsUpper.

Result:
    True, if the inversion succeeded.
    False, if matrix A contains zero elements on its main diagonal.
    Matrix A could not be inverted.

The algorithm is the modification of DPOTRI and DLAUU2 subroutines from
LAPACK library.
*************************************************************************)
function SPDMatrixCholeskyInverse(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    AJJ : Double;
    AII : Double;
    T : TReal1DArray;
    A1 : TReal2DArray;
    i_ : AlglibInteger;
begin
    Result := True;
    
    //
    // Test the input parameters.
    //
    SetLength(T, N-1+1);
    if IsUpper then
    begin
        
        //
        // Compute inverse of upper triangular matrix.
        //
        J:=0;
        while J<=N-1 do
        begin
            if AP_FP_Eq(A[J,J],0) then
            begin
                Result := False;
                Exit;
            end;
            A[J,J] := 1/A[J,J];
            AJJ := -A[J,J];
            
            //
            // Compute elements 1:j-1 of j-th column.
            //
            for i_ := 0 to J-1 do
            begin
                T[i_] := A[i_,J];
            end;
            I:=0;
            while I<=J-1 do
            begin
                V := APVDotProduct(@A[I][0], I, J-1, @T[0], I, J-1);
                A[I,J] := V;
                Inc(I);
            end;
            for i_ := 0 to J-1 do
            begin
                A[i_,J] := AJJ*A[i_,J];
            end;
            Inc(J);
        end;
        
        //
        // InvA = InvU * InvU'
        //
        I:=0;
        while I<=N-1 do
        begin
            AII := A[I,I];
            if I<N-1 then
            begin
                V := APVDotProduct(@A[I][0], I, N-1, @A[I][0], I, N-1);
                A[I,I] := V;
                K:=0;
                while K<=I-1 do
                begin
                    V := APVDotProduct(@A[K][0], I+1, N-1, @A[I][0], I+1, N-1);
                    A[K,I] := A[K,I]*AII+V;
                    Inc(K);
                end;
            end
            else
            begin
                for i_ := 0 to I do
                begin
                    A[i_,I] := AII*A[i_,I];
                end;
            end;
            Inc(I);
        end;
    end
    else
    begin
        
        //
        // Compute inverse of lower triangular matrix.
        //
        J:=N-1;
        while J>=0 do
        begin
            if AP_FP_Eq(A[J,J],0) then
            begin
                Result := False;
                Exit;
            end;
            A[J,J] := 1/A[J,J];
            AJJ := -A[J,J];
            if J<N-1 then
            begin
                
                //
                // Compute elements j+1:n of j-th column.
                //
                for i_ := J+1 to N-1 do
                begin
                    T[i_] := A[i_,J];
                end;
                I:=J+1+1;
                while I<=N do
                begin
                    V := APVDotProduct(@A[I-1][0], J+1, I-1, @T[0], J+1, I-1);
                    A[I-1,J] := V;
                    Inc(I);
                end;
                for i_ := J+1 to N-1 do
                begin
                    A[i_,J] := AJJ*A[i_,J];
                end;
            end;
            Dec(J);
        end;
        
        //
        // InvA = InvL' * InvL
        //
        I:=0;
        while I<=N-1 do
        begin
            AII := A[I,I];
            if I<N-1 then
            begin
                V := 0.0;
                for i_ := I to N-1 do
                begin
                    V := V + A[i_,I]*A[i_,I];
                end;
                A[I,I] := V;
                K:=0;
                while K<=I-1 do
                begin
                    V := 0.0;
                    for i_ := I+1 to N-1 do
                    begin
                        V := V + A[i_,K]*A[i_,I];
                    end;
                    A[I,K] := AII*A[I,K]+V;
                    Inc(K);
                end;
            end
            else
            begin
                APVMul(@A[I][0], 0, I, AII);
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Inversion of a symmetric positive definite matrix.

Given an upper or lower triangle of a symmetric positive definite matrix,
the algorithm generates matrix A^-1 and saves the upper or lower triangle
depending on the input.

Input parameters:
    A       -   matrix to be inverted (upper or lower triangle).
                Array with elements [0..N-1,0..N-1].
    N       -   size of matrix A.
    IsUpper -   storage format.
                If IsUpper = True, then the upper triangle of matrix A is
                given, otherwise the lower triangle is given.

Output parameters:
    A       -   inverse of matrix A.
                Array with elements [0..N-1,0..N-1].
                If IsUpper = True, then the upper triangle of matrix A^-1
                is used, and the elements below the main diagonal are not
                used nor changed. The same applies if IsUpper = False.

Result:
    True, if the matrix is positive definite.
    False, if the matrix is not positive definite (and it could not be
    inverted by this algorithm).
*************************************************************************)
function SPDMatrixInverse(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
begin
    Result := False;
    if SPDMatrixCholesky(A, N, IsUpper) then
    begin
        if SPDMatrixCholeskyInverse(A, N, IsUpper) then
        begin
            Result := True;
        end;
    end;
end;


end.