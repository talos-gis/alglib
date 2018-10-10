(*************************************************************************
Copyright 2008 by Sergey Bochkanov (ALGLIB project).

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
unit tsort;
interface
uses Math, Sysutils, Ap;

procedure TagSort(var A : TReal1DArray;
     N : AlglibInteger;
     var P1 : TInteger1DArray;
     var P2 : TInteger1DArray);
procedure TagSortFastI(var A : TReal1DArray;
     var B : TInteger1DArray;
     N : AlglibInteger);
procedure TagSortFastR(var A : TReal1DArray;
     var B : TReal1DArray;
     N : AlglibInteger);
procedure TagSortFast(var A : TReal1DArray; N : AlglibInteger);

implementation

procedure TagSort(var A : TReal1DArray;
     N : AlglibInteger;
     var P1 : TInteger1DArray;
     var P2 : TInteger1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    T : AlglibInteger;
    Tmp : Double;
    TmpI : AlglibInteger;
    PV : TInteger1DArray;
    VP : TInteger1DArray;
    LV : AlglibInteger;
    LP : AlglibInteger;
    RV : AlglibInteger;
    RP : AlglibInteger;
begin
    
    //
    // Special cases
    //
    if N<=0 then
    begin
        Exit;
    end;
    if N=1 then
    begin
        SetLength(P1, 0+1);
        SetLength(P2, 0+1);
        P1[0] := 0;
        P2[0] := 0;
        Exit;
    end;
    
    //
    // General case, N>1: prepare permutations table P1
    //
    SetLength(P1, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        P1[I] := I;
        Inc(I);
    end;
    
    //
    // General case, N>1: sort, update P1
    //
    TagSortFastI(A, P1, N);
    
    //
    // General case, N>1: fill permutations table P2
    //
    // To fill P2 we maintain two arrays:
    // * PV, Position(Value). PV[i] contains position of I-th key at the moment
    // * VP, Value(Position). VP[i] contains key which has position I at the moment
    //
    // At each step we making permutation of two items:
    //   Left, which is given by position/value pair LP/LV
    //   and Right, which is given by RP/RV
    // and updating PV[] and VP[] correspondingly.
    //
    SetLength(PV, N-1+1);
    SetLength(VP, N-1+1);
    SetLength(P2, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        PV[I] := I;
        VP[I] := I;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        
        //
        // calculate LP, LV, RP, RV
        //
        LP := I;
        LV := VP[LP];
        RV := P1[I];
        RP := PV[RV];
        
        //
        // Fill P2
        //
        P2[I] := RP;
        
        //
        // update PV and VP
        //
        VP[LP] := RV;
        VP[RP] := LV;
        PV[LV] := RP;
        PV[RV] := LP;
        Inc(I);
    end;
end;


procedure TagSortFastI(var A : TReal1DArray;
     var B : TInteger1DArray;
     N : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    T : AlglibInteger;
    Tmp : Double;
    TmpI : AlglibInteger;
begin
    
    //
    // Special cases
    //
    if N<=1 then
    begin
        Exit;
    end;
    
    //
    // General case, N>1: sort, update B
    //
    i := 2;
    repeat
        t := i;
        while t<>1 do
        begin
            k := t div 2;
            if AP_FP_Greater_Eq(A[k-1],A[t-1]) then
            begin
                t := 1;
            end
            else
            begin
                Tmp := A[k-1];
                A[k-1] := A[t-1];
                A[t-1] := Tmp;
                TmpI := B[k-1];
                B[k-1] := B[t-1];
                B[t-1] := TmpI;
                t := k;
            end;
        end;
        i := i+1;
    until  not (i<=n);
    i := n-1;
    repeat
        Tmp := A[i];
        A[i] := A[0];
        A[0] := Tmp;
        TmpI := B[i];
        B[i] := B[0];
        B[0] := TmpI;
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
                    if AP_FP_Greater(A[k],A[k-1]) then
                    begin
                        k := k+1;
                    end;
                end;
                if AP_FP_Greater_Eq(A[t-1],A[k-1]) then
                begin
                    t := 0;
                end
                else
                begin
                    Tmp := A[k-1];
                    A[k-1] := A[t-1];
                    A[t-1] := Tmp;
                    TmpI := B[k-1];
                    B[k-1] := B[t-1];
                    B[t-1] := TmpI;
                    t := k;
                end;
            end;
        end;
        i := i-1;
    until  not (i>=1);
end;


procedure TagSortFastR(var A : TReal1DArray;
     var B : TReal1DArray;
     N : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    T : AlglibInteger;
    Tmp : Double;
    TmpR : Double;
begin
    
    //
    // Special cases
    //
    if N<=1 then
    begin
        Exit;
    end;
    
    //
    // General case, N>1: sort, update B
    //
    i := 2;
    repeat
        t := i;
        while t<>1 do
        begin
            k := t div 2;
            if AP_FP_Greater_Eq(A[k-1],A[t-1]) then
            begin
                t := 1;
            end
            else
            begin
                Tmp := A[k-1];
                A[k-1] := A[t-1];
                A[t-1] := Tmp;
                TmpR := B[k-1];
                B[k-1] := B[t-1];
                B[t-1] := TmpR;
                t := k;
            end;
        end;
        i := i+1;
    until  not (i<=n);
    i := n-1;
    repeat
        Tmp := A[i];
        A[i] := A[0];
        A[0] := Tmp;
        TmpR := B[i];
        B[i] := B[0];
        B[0] := TmpR;
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
                    if AP_FP_Greater(A[k],A[k-1]) then
                    begin
                        k := k+1;
                    end;
                end;
                if AP_FP_Greater_Eq(A[t-1],A[k-1]) then
                begin
                    t := 0;
                end
                else
                begin
                    Tmp := A[k-1];
                    A[k-1] := A[t-1];
                    A[t-1] := Tmp;
                    TmpR := B[k-1];
                    B[k-1] := B[t-1];
                    B[t-1] := TmpR;
                    t := k;
                end;
            end;
        end;
        i := i-1;
    until  not (i>=1);
end;


procedure TagSortFast(var A : TReal1DArray; N : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    T : AlglibInteger;
    Tmp : Double;
    TmpR : Double;
begin
    
    //
    // Special cases
    //
    if N<=1 then
    begin
        Exit;
    end;
    
    //
    // General case, N>1: sort, update B
    //
    i := 2;
    repeat
        t := i;
        while t<>1 do
        begin
            k := t div 2;
            if AP_FP_Greater_Eq(A[k-1],A[t-1]) then
            begin
                t := 1;
            end
            else
            begin
                Tmp := A[k-1];
                A[k-1] := A[t-1];
                A[t-1] := Tmp;
                t := k;
            end;
        end;
        i := i+1;
    until  not (i<=n);
    i := n-1;
    repeat
        Tmp := A[i];
        A[i] := A[0];
        A[0] := Tmp;
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
                    if AP_FP_Greater(A[k],A[k-1]) then
                    begin
                        k := k+1;
                    end;
                end;
                if AP_FP_Greater_Eq(A[t-1],A[k-1]) then
                begin
                    t := 0;
                end
                else
                begin
                    Tmp := A[k-1];
                    A[k-1] := A[t-1];
                    A[t-1] := Tmp;
                    t := k;
                end;
            end;
        end;
        i := i-1;
    until  not (i>=1);
end;


end.