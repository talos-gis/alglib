(*************************************************************************
Copyright (c) 2005-2007, Sergey Bochkanov (ALGLIB project).

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
unit cblas;
interface
uses Math, Sysutils, Ap;

procedure ComplexMatrixVectorMultiply(const A : TComplex2DArray;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     J1 : AlglibInteger;
     J2 : AlglibInteger;
     TransA : Boolean;
     ConjA : Boolean;
     const X : TComplex1DArray;
     IX1 : AlglibInteger;
     IX2 : AlglibInteger;
     Alpha : Complex;
     var Y : TComplex1DArray;
     IY1 : AlglibInteger;
     IY2 : AlglibInteger;
     Beta : Complex;
     var T : TComplex1DArray);

implementation

procedure ComplexMatrixVectorMultiply(const A : TComplex2DArray;
     I1 : AlglibInteger;
     I2 : AlglibInteger;
     J1 : AlglibInteger;
     J2 : AlglibInteger;
     TransA : Boolean;
     ConjA : Boolean;
     const X : TComplex1DArray;
     IX1 : AlglibInteger;
     IX2 : AlglibInteger;
     Alpha : Complex;
     var Y : TComplex1DArray;
     IY1 : AlglibInteger;
     IY2 : AlglibInteger;
     Beta : Complex;
     var T : TComplex1DArray);
var
    I : AlglibInteger;
    V : Complex;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if  not TransA then
    begin
        
        //
        // y := alpha*A*x + beta*y
        //
        // or
        //
        // y := alpha*conj(A)*x + beta*y
        //
        if (I1>I2) or (J1>J2) then
        begin
            Exit;
        end;
        Assert(J2-J1=IX2-IX1, 'ComplexMatrixVectorMultiply: A and X dont match!');
        Assert(I2-I1=IY2-IY1, 'ComplexMatrixVectorMultiply: A and Y dont match!');
        
        //
        // beta*y
        //
        if C_EqualR(Beta,0) then
        begin
            I:=IY1;
            while I<=IY2 do
            begin
                Y[I] := C_Complex(0);
                Inc(I);
            end;
        end
        else
        begin
            for i_ := IY1 to IY2 do
            begin
                Y[i_] := C_Mul(Beta, Y[i_]);
            end;
        end;
        
        //
        // conj?
        //
        if ConjA then
        begin
            for i_ := IX1 to IX2 do
            begin
                T[i_] := Conj(X[i_]);
            end;
            Alpha := Conj(Alpha);
            for i_ := IY1 to IY2 do
            begin
                Y[i_] := Conj(Y[i_]);
            end;
        end
        else
        begin
            for i_ := IX1 to IX2 do
            begin
                T[i_] := X[i_];
            end;
        end;
        
        //
        // alpha*A*x
        //
        I:=I1;
        while I<=I2 do
        begin
            i1_ := (IX1)-(J1);
            V := C_Complex(0.0);
            for i_ := J1 to J2 do
            begin
                V := C_Add(V,C_Mul(A[I,i_],X[i_+i1_]));
            end;
            Y[IY1+I-I1] := C_Add(Y[IY1+I-I1],C_Mul(Alpha,V));
            Inc(I);
        end;
        
        //
        // conj?
        //
        if ConjA then
        begin
            for i_ := IY1 to IY2 do
            begin
                Y[i_] := Conj(Y[i_]);
            end;
        end;
    end
    else
    begin
        
        //
        // y := alpha*A'*x + beta*y;
        //
        // or
        //
        // y := alpha*conj(A')*x + beta*y;
        //
        if (I1>I2) or (J1>J2) then
        begin
            Exit;
        end;
        Assert(I2-I1=IX2-IX1, 'ComplexMatrixVectorMultiply: A and X dont match!');
        Assert(J2-J1=IY2-IY1, 'ComplexMatrixVectorMultiply: A and Y dont match!');
        
        //
        // beta*y
        //
        if C_EqualR(Beta,0) then
        begin
            I:=IY1;
            while I<=IY2 do
            begin
                Y[I] := C_Complex(0);
                Inc(I);
            end;
        end
        else
        begin
            for i_ := IY1 to IY2 do
            begin
                Y[i_] := C_Mul(Beta, Y[i_]);
            end;
        end;
        
        //
        // conj?
        //
        if ConjA then
        begin
            for i_ := IX1 to IX2 do
            begin
                T[i_] := Conj(X[i_]);
            end;
            Alpha := Conj(Alpha);
            for i_ := IY1 to IY2 do
            begin
                Y[i_] := Conj(Y[i_]);
            end;
        end
        else
        begin
            for i_ := IX1 to IX2 do
            begin
                T[i_] := X[i_];
            end;
        end;
        
        //
        // alpha*A'*x
        //
        I:=I1;
        while I<=I2 do
        begin
            V := C_Mul(Alpha,X[IX1+I-I1]);
            i1_ := (J1) - (IY1);
            for i_ := IY1 to IY2 do
            begin
                Y[i_] := C_Add(Y[i_], C_Mul(V, A[I,i_+i1_]));
            end;
            Inc(I);
        end;
        
        //
        // conj?
        //
        if ConjA then
        begin
            for i_ := IY1 to IY2 do
            begin
                Y[i_] := Conj(Y[i_]);
            end;
        end;
    end;
end;


end.