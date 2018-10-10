(*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

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
unit taskgen;
interface
uses Math, Sysutils, Ap;

procedure TaskGenInt1D(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);
procedure TaskGenInt1DEquidist(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);
procedure TaskGenInt1DCheb1(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);
procedure TaskGenInt1DCheb2(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);

implementation

(*************************************************************************
This  function  generates  1-dimensional  general  interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure TaskGenInt1D(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);
var
    I : AlglibInteger;
    H : Double;
begin
    Assert(N>=1, 'TaskGenInterpolationEqdist1D: N<1!');
    SetLength(X, N);
    SetLength(Y, N);
    if N>1 then
    begin
        X[0] := A;
        Y[0] := 2*RandomReal-1;
        H := (B-A)/(N-1);
        I:=1;
        while I<=N-1 do
        begin
            if I<>N-1 then
            begin
                X[I] := A+(I+0.2*(2*RandomReal-1))*H;
            end
            else
            begin
                X[I] := B;
            end;
            Y[I] := Y[I-1]+(2*RandomReal-1)*(X[I]-X[I-1]);
            Inc(I);
        end;
    end
    else
    begin
        X[0] := 0.5*(A+B);
        Y[0] := 2*RandomReal-1;
    end;
end;


(*************************************************************************
This function generates  1-dimensional equidistant interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure TaskGenInt1DEquidist(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);
var
    I : AlglibInteger;
    H : Double;
begin
    Assert(N>=1, 'TaskGenInterpolationEqdist1D: N<1!');
    SetLength(X, N);
    SetLength(Y, N);
    if N>1 then
    begin
        X[0] := A;
        Y[0] := 2*RandomReal-1;
        H := (B-A)/(N-1);
        I:=1;
        while I<=N-1 do
        begin
            X[I] := A+I*H;
            Y[I] := Y[I-1]+(2*RandomReal-1)*H;
            Inc(I);
        end;
    end
    else
    begin
        X[0] := 0.5*(A+B);
        Y[0] := 2*RandomReal-1;
    end;
end;


(*************************************************************************
This function generates  1-dimensional Chebyshev-1 interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure TaskGenInt1DCheb1(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);
var
    I : AlglibInteger;
    H : Double;
begin
    Assert(N>=1, 'TaskGenInterpolation1DCheb1: N<1!');
    SetLength(X, N);
    SetLength(Y, N);
    if N>1 then
    begin
        I:=0;
        while I<=N-1 do
        begin
            X[I] := 0.5*(B+A)+0.5*(B-A)*Cos(Pi*(2*i+1)/(2*N));
            if I=0 then
            begin
                Y[I] := 2*RandomReal-1;
            end
            else
            begin
                Y[I] := Y[I-1]+(2*RandomReal-1)*(X[I]-X[I-1]);
            end;
            Inc(I);
        end;
    end
    else
    begin
        X[0] := 0.5*(A+B);
        Y[0] := 2*RandomReal-1;
    end;
end;


(*************************************************************************
This function generates  1-dimensional Chebyshev-2 interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************)
procedure TaskGenInt1DCheb2(A : Double;
     B : Double;
     N : AlglibInteger;
     var X : TReal1DArray;
     var Y : TReal1DArray);
var
    I : AlglibInteger;
    H : Double;
begin
    Assert(N>=1, 'TaskGenInterpolation1DCheb2: N<1!');
    SetLength(X, N);
    SetLength(Y, N);
    if N>1 then
    begin
        I:=0;
        while I<=N-1 do
        begin
            X[I] := 0.5*(B+A)+0.5*(B-A)*Cos(Pi*i/(n-1));
            if I=0 then
            begin
                Y[I] := 2*RandomReal-1;
            end
            else
            begin
                Y[I] := Y[I-1]+(2*RandomReal-1)*(X[I]-X[I-1]);
            end;
            Inc(I);
        end;
    end
    else
    begin
        X[0] := 0.5*(A+B);
        Y[0] := 2*RandomReal-1;
    end;
end;


end.