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
unit cholesky;
interface
uses Math, Sysutils, Ap;

function SPDMatrixCholesky(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
function CholeskyDecomposition(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;

implementation

(*************************************************************************
Cholesky decomposition

The algorithm computes Cholesky decomposition of a symmetric
positive-definite matrix.

The result of an algorithm is a representation of matrix A as A = U'*U or
A = L*L'.

Input parameters:
    A       -   upper or lower triangle of a factorized matrix.
                array with elements [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   if IsUpper=True, then A contains an upper triangle of
                a symmetric matrix, otherwise A contains a lower one.

Output parameters:
    A       -   the result of factorization. If IsUpper=True, then
                the upper triangle contains matrix U, so that A = U'*U,
                and the elements below the main diagonal are not modified.
                Similarly, if IsUpper = False.

Result:
    If the matrix is positive-definite, the function returns True.
    Otherwise, the function returns False. This means that the
    factorization could not be carried out.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************)
function SPDMatrixCholesky(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    AJJ : Double;
    v : Double;
    i_ : AlglibInteger;
begin
    
    //
    //     Test the input parameters.
    //
    Assert(N>=0, 'Error in SMatrixCholesky: incorrect function arguments');
    
    //
    //     Quick return if possible
    //
    Result := True;
    if N<=0 then
    begin
        Exit;
    end;
    if IsUpper then
    begin
        
        //
        // Compute the Cholesky factorization A = U'*U.
        //
        J:=0;
        while J<=N-1 do
        begin
            
            //
            // Compute U(J,J) and test for non-positive-definiteness.
            //
            V := 0.0;
            for i_ := 0 to J-1 do
            begin
                V := V + A[i_,J]*A[i_,J];
            end;
            AJJ := A[J,J]-V;
            if AP_FP_Less_Eq(AJJ,0) then
            begin
                Result := False;
                Exit;
            end;
            AJJ := SQRT(AJJ);
            A[J,J] := AJJ;
            
            //
            // Compute elements J+1:N of row J.
            //
            if J<N-1 then
            begin
                I:=0;
                while I<=J-1 do
                begin
                    V := A[I,J];
                    APVSub(@A[J][0], J+1, N-1, @A[I][0], J+1, N-1, V);
                    Inc(I);
                end;
                V := 1/AJJ;
                APVMul(@A[J][0], J+1, N-1, V);
            end;
            Inc(J);
        end;
    end
    else
    begin
        
        //
        // Compute the Cholesky factorization A = L*L'.
        //
        J:=0;
        while J<=N-1 do
        begin
            
            //
            // Compute L(J,J) and test for non-positive-definiteness.
            //
            V := APVDotProduct(@A[J][0], 0, J-1, @A[J][0], 0, J-1);
            AJJ := A[J,J]-V;
            if AP_FP_Less_Eq(AJJ,0) then
            begin
                Result := False;
                Exit;
            end;
            AJJ := SQRT(AJJ);
            A[J,J] := AJJ;
            
            //
            // Compute elements J+1:N of column J.
            //
            if J<N-1 then
            begin
                I:=J+1;
                while I<=N-1 do
                begin
                    V := APVDotProduct(@A[I][0], 0, J-1, @A[J][0], 0, J-1);
                    A[I,J] := A[I,J]-V;
                    Inc(I);
                end;
                V := 1/AJJ;
                for i_ := J+1 to N-1 do
                begin
                    A[i_,J] := V*A[i_,J];
                end;
            end;
            Inc(J);
        end;
    end;
end;


function CholeskyDecomposition(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    AJJ : Double;
    v : Double;
    JM1 : AlglibInteger;
    JP1 : AlglibInteger;
    i_ : AlglibInteger;
begin
    
    //
    //     Test the input parameters.
    //
    Assert(N>=0, 'Error in CholeskyDecomposition: incorrect function arguments');
    
    //
    //     Quick return if possible
    //
    Result := True;
    if N=0 then
    begin
        Exit;
    end;
    if IsUpper then
    begin
        
        //
        // Compute the Cholesky factorization A = U'*U.
        //
        J:=1;
        while J<=N do
        begin
            
            //
            // Compute U(J,J) and test for non-positive-definiteness.
            //
            JM1 := J-1;
            v := 0.0;
            for i_ := 1 to JM1 do
            begin
                v := v + A[i_,J]*A[i_,J];
            end;
            AJJ := A[J,J]-v;
            if AP_FP_Less_Eq(AJJ,0) then
            begin
                Result := False;
                Exit;
            end;
            AJJ := SQRT(AJJ);
            A[J,J] := AJJ;
            
            //
            // Compute elements J+1:N of row J.
            //
            if J<N then
            begin
                I:=J+1;
                while I<=N do
                begin
                    JM1 := J-1;
                    V := 0.0;
                    for i_ := 1 to JM1 do
                    begin
                        V := V + A[i_,I]*A[i_,J];
                    end;
                    A[J,I] := A[J,I]-V;
                    Inc(I);
                end;
                V := 1/AJJ;
                JP1 := J+1;
                APVMul(@A[J][0], JP1, N, V);
            end;
            Inc(J);
        end;
    end
    else
    begin
        
        //
        // Compute the Cholesky factorization A = L*L'.
        //
        J:=1;
        while J<=N do
        begin
            
            //
            // Compute L(J,J) and test for non-positive-definiteness.
            //
            JM1 := J-1;
            V := APVDotProduct(@A[J][0], 1, JM1, @A[J][0], 1, JM1);
            AJJ := A[J,J]-V;
            if AP_FP_Less_Eq(AJJ,0) then
            begin
                Result := False;
                Exit;
            end;
            AJJ := SQRT(AJJ);
            A[J,J] := AJJ;
            
            //
            // Compute elements J+1:N of column J.
            //
            if J<N then
            begin
                I:=J+1;
                while I<=N do
                begin
                    JM1 := J-1;
                    V := APVDotProduct(@A[I][0], 1, JM1, @A[J][0], 1, JM1);
                    A[I,J] := A[I,J]-V;
                    Inc(I);
                end;
                V := 1/AJJ;
                JP1 := J+1;
                for i_ := JP1 to N do
                begin
                    A[i_,J] := V*A[i_,J];
                end;
            end;
            Inc(J);
        end;
    end;
end;


end.