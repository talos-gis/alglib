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
unit trinverse;
interface
uses Math, Sysutils, Ap;

function RMatrixTRInverse(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsunitTriangular : Boolean):Boolean;
function InvTriangular(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsunitTriangular : Boolean):Boolean;

implementation

(*************************************************************************
Triangular matrix inversion

The subroutine inverts the following types of matrices:
    * upper triangular
    * upper triangular with unit diagonal
    * lower triangular
    * lower triangular with unit diagonal

In case of an upper (lower) triangular matrix,  the  inverse  matrix  will
also be upper (lower) triangular, and after the end of the algorithm,  the
inverse matrix replaces the source matrix. The elements  below (above) the
main diagonal are not changed by the algorithm.

If  the matrix  has a unit diagonal, the inverse matrix also  has  a  unit
diagonal, and the diagonal elements are not passed to the algorithm.

Input parameters:
    A       -   matrix.
                Array whose indexes range within [0..N-1, 0..N-1].
    N       -   size of matrix A.
    IsUpper -   True, if the matrix is upper triangular.
    IsunitTriangular
            -   True, if the matrix has a unit diagonal.

Output parameters:
    A       -   inverse matrix (if the problem is not degenerate).

Result:
    True, if the matrix is not singular.
    False, if the matrix is singular.

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************)
function RMatrixTRInverse(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsunitTriangular : Boolean):Boolean;
var
    NOunit : Boolean;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    AJJ : Double;
    T : TReal1DArray;
    i_ : AlglibInteger;
begin
    Result := True;
    SetLength(T, N-1+1);
    
    //
    // Test the input parameters.
    //
    NOunit :=  not IsunitTriangular;
    if IsUpper then
    begin
        
        //
        // Compute inverse of upper triangular matrix.
        //
        J:=0;
        while J<=N-1 do
        begin
            if NOunit then
            begin
                if AP_FP_Eq(A[J,J],0) then
                begin
                    Result := False;
                    Exit;
                end;
                A[J,J] := 1/A[J,J];
                AJJ := -A[J,J];
            end
            else
            begin
                AJJ := -1;
            end;
            
            //
            // Compute elements 1:j-1 of j-th column.
            //
            if J>0 then
            begin
                for i_ := 0 to J-1 do
                begin
                    T[i_] := A[i_,J];
                end;
                I:=0;
                while I<=J-1 do
                begin
                    if I<J-1 then
                    begin
                        V := APVDotProduct(@A[I][0], I+1, J-1, @T[0], I+1, J-1);
                    end
                    else
                    begin
                        V := 0;
                    end;
                    if NOunit then
                    begin
                        A[I,J] := V+A[I,I]*T[I];
                    end
                    else
                    begin
                        A[I,J] := V+T[I];
                    end;
                    Inc(I);
                end;
                for i_ := 0 to J-1 do
                begin
                    A[i_,J] := AJJ*A[i_,J];
                end;
            end;
            Inc(J);
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
            if NOunit then
            begin
                if AP_FP_Eq(A[J,J],0) then
                begin
                    Result := False;
                    Exit;
                end;
                A[J,J] := 1/A[J,J];
                AJJ := -A[J,J];
            end
            else
            begin
                AJJ := -1;
            end;
            if J<N-1 then
            begin
                
                //
                // Compute elements j+1:n of j-th column.
                //
                for i_ := J+1 to N-1 do
                begin
                    T[i_] := A[i_,J];
                end;
                I:=J+1;
                while I<=N-1 do
                begin
                    if I>J+1 then
                    begin
                        V := APVDotProduct(@A[I][0], J+1, I-1, @T[0], J+1, I-1);
                    end
                    else
                    begin
                        V := 0;
                    end;
                    if NOunit then
                    begin
                        A[I,J] := V+A[I,I]*T[I];
                    end
                    else
                    begin
                        A[I,J] := V+T[I];
                    end;
                    Inc(I);
                end;
                for i_ := J+1 to N-1 do
                begin
                    A[i_,J] := AJJ*A[i_,J];
                end;
            end;
            Dec(J);
        end;
    end;
end;


function InvTriangular(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsunitTriangular : Boolean):Boolean;
var
    NOunit : Boolean;
    I : AlglibInteger;
    J : AlglibInteger;
    NMJ : AlglibInteger;
    JM1 : AlglibInteger;
    JP1 : AlglibInteger;
    V : Double;
    AJJ : Double;
    T : TReal1DArray;
    i_ : AlglibInteger;
begin
    Result := True;
    SetLength(T, N+1);
    
    //
    // Test the input parameters.
    //
    NOunit :=  not IsunitTriangular;
    if IsUpper then
    begin
        
        //
        // Compute inverse of upper triangular matrix.
        //
        J:=1;
        while J<=N do
        begin
            if NOunit then
            begin
                if AP_FP_Eq(A[J,J],0) then
                begin
                    Result := False;
                    Exit;
                end;
                A[J,J] := 1/A[J,J];
                AJJ := -A[J,J];
            end
            else
            begin
                AJJ := -1;
            end;
            
            //
            // Compute elements 1:j-1 of j-th column.
            //
            if J>1 then
            begin
                JM1 := J-1;
                for i_ := 1 to JM1 do
                begin
                    T[i_] := A[i_,J];
                end;
                I:=1;
                while I<=J-1 do
                begin
                    if I<J-1 then
                    begin
                        V := APVDotProduct(@A[I][0], I+1, JM1, @T[0], I+1, JM1);
                    end
                    else
                    begin
                        V := 0;
                    end;
                    if NOunit then
                    begin
                        A[I,J] := V+A[I,I]*T[I];
                    end
                    else
                    begin
                        A[I,J] := V+T[I];
                    end;
                    Inc(I);
                end;
                for i_ := 1 to JM1 do
                begin
                    A[i_,J] := AJJ*A[i_,J];
                end;
            end;
            Inc(J);
        end;
    end
    else
    begin
        
        //
        // Compute inverse of lower triangular matrix.
        //
        J:=N;
        while J>=1 do
        begin
            if NOunit then
            begin
                if AP_FP_Eq(A[J,J],0) then
                begin
                    Result := False;
                    Exit;
                end;
                A[J,J] := 1/A[J,J];
                AJJ := -A[J,J];
            end
            else
            begin
                AJJ := -1;
            end;
            if J<N then
            begin
                
                //
                // Compute elements j+1:n of j-th column.
                //
                NMJ := N-J;
                JP1 := J+1;
                for i_ := JP1 to N do
                begin
                    T[i_] := A[i_,J];
                end;
                I:=J+1;
                while I<=N do
                begin
                    if I>J+1 then
                    begin
                        V := APVDotProduct(@A[I][0], JP1, I-1, @T[0], JP1, I-1);
                    end
                    else
                    begin
                        V := 0;
                    end;
                    if NOunit then
                    begin
                        A[I,J] := V+A[I,I]*T[I];
                    end
                    else
                    begin
                        A[I,J] := V+T[I];
                    end;
                    Inc(I);
                end;
                for i_ := JP1 to N do
                begin
                    A[i_,J] := AJJ*A[i_,J];
                end;
            end;
            Dec(J);
        end;
    end;
end;


end.