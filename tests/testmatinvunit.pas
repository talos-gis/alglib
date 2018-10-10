unit testmatinvunit;
interface
uses Math, Sysutils, Ap, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv;

function TestMatInv(Silent : Boolean):Boolean;
function testmatinvunit_test_silent():Boolean;
function testmatinvunit_test():Boolean;

implementation

procedure RMatrixMakeACopy(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);forward;
procedure CMatrixMakeACopy(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TComplex2DArray);forward;
function RMatrixCheckInverse(const A : TReal2DArray;
     const InvA : TReal2DArray;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : MatInvReport):Boolean;forward;
function SPDMatrixCheckInverse(A : TReal2DArray;
     InvA : TReal2DArray;
     IsUpper : Boolean;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : MatInvReport):Boolean;forward;
function HPDMatrixCheckInverse(A : TComplex2DArray;
     InvA : TComplex2DArray;
     IsUpper : Boolean;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : MatInvReport):Boolean;forward;
function RMatrixCheckInverseSingular(const InvA : TReal2DArray;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : MatInvReport):Boolean;forward;
function CMatrixCheckInverse(const A : TComplex2DArray;
     const InvA : TComplex2DArray;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : MatInvReport):Boolean;forward;
function CMatrixCheckInverseSingular(const InvA : TComplex2DArray;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : MatInvReport):Boolean;forward;
procedure RMatrixDropHalf(var A : TReal2DArray;
     N : AlglibInteger;
     DropLower : Boolean);forward;
procedure CMatrixDropHalf(var A : TComplex2DArray;
     N : AlglibInteger;
     DropLower : Boolean);forward;
procedure TestRTRInv(MaxN : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var RTRErrors : Boolean);forward;
procedure TestCTRInv(MaxN : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var CTRErrors : Boolean);forward;
procedure TestRInv(MaxN : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var RErrors : Boolean);forward;
procedure TestCInv(MaxN : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var CErrors : Boolean);forward;
procedure TestSPDInv(MaxN : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var SPDErrors : Boolean);forward;
procedure TestHPDInv(MaxN : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var HPDErrors : Boolean);forward;
procedure Unset2D(var X : TReal2DArray);forward;
procedure Unset1D(var X : TReal1DArray);forward;
procedure CUnset2D(var X : TComplex2DArray);forward;
procedure CUnset1D(var X : TComplex1DArray);forward;
procedure UnsetRep(var R : MatInvReport);forward;


(*************************************************************************
Test
*************************************************************************)
function TestMatInv(Silent : Boolean):Boolean;
var
    MaxRN : AlglibInteger;
    MaxCN : AlglibInteger;
    PassCount : AlglibInteger;
    Threshold : Double;
    RCondTol : Double;
    RTRErrors : Boolean;
    CTRErrors : Boolean;
    RErrors : Boolean;
    CErrors : Boolean;
    SPDErrors : Boolean;
    HPDErrors : Boolean;
    WasErrors : Boolean;
    EmptyRA : TReal2DArray;
    EmptyCA : TReal2DArray;
begin
    MaxRN := 3*ABLASBlockSize(EmptyRA)+1;
    MaxCN := 3*ABLASBlockSize(EmptyCA)+1;
    PassCount := 1;
    Threshold := 10000*MachineEpsilon;
    RCondTol := 0.01;
    RTRErrors := False;
    CTRErrors := False;
    RErrors := False;
    CErrors := False;
    SPDErrors := False;
    HPDErrors := False;
    TestRTRInv(MaxRN, PassCount, Threshold, RTRErrors);
    TestCTRInv(MaxCN, PassCount, Threshold, CTRErrors);
    TestRInv(MaxRN, PassCount, Threshold, RErrors);
    TestSPDInv(MaxRN, PassCount, Threshold, SPDErrors);
    TestCInv(MaxCN, PassCount, Threshold, CErrors);
    TestHPDInv(MaxCN, PassCount, Threshold, HPDErrors);
    WasErrors := RTRErrors or CTRErrors or RErrors or CErrors or SPDErrors or HPDErrors;
    if  not Silent then
    begin
        Write(Format('TESTING MATINV'#13#10'',[]));
        Write(Format('* REAL TRIANGULAR:                        ',[]));
        if RTRErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* COMPLEX TRIANGULAR:                     ',[]));
        if CTRErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* REAL:                                   ',[]));
        if RErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* COMPLEX:                                ',[]));
        if CErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* SPD:                                    ',[]));
        if SPDErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* HPD:                                    ',[]));
        if HPDErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        if WasErrors then
        begin
            Write(Format('TEST FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('TEST PASSED'#13#10'',[]));
        end;
    end;
    Result :=  not WasErrors;
end;


(*************************************************************************
Copy
*************************************************************************)
procedure RMatrixMakeACopy(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    SetLength(B, M-1+1, N-1+1);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            B[I,J] := A[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Copy
*************************************************************************)
procedure CMatrixMakeACopy(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TComplex2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    SetLength(B, M-1+1, N-1+1);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            B[I,J] := A[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Checks whether inverse is correct
Returns True on success.
*************************************************************************)
function RMatrixCheckInverse(const A : TReal2DArray;
     const InvA : TReal2DArray;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : MatInvReport):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
begin
    Result := True;
    if Info<=0 then
    begin
        Result := False;
    end
    else
    begin
        Result := Result and  not (AP_FP_Less(Rep.R1,100*MachineEpsilon) or AP_FP_Greater(Rep.R1,1+1000*MachineEpsilon));
        Result := Result and  not (AP_FP_Less(Rep.RInf,100*MachineEpsilon) or AP_FP_Greater(Rep.RInf,1+1000*MachineEpsilon));
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                V := 0.0;
                for i_ := 0 to N-1 do
                begin
                    V := V + A[I,i_]*InvA[i_,J];
                end;
                if I=J then
                begin
                    V := V-1;
                end;
                Result := Result and AP_FP_Less_Eq(AbsReal(V),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Checks whether inverse is correct
Returns True on success.
*************************************************************************)
function SPDMatrixCheckInverse(A : TReal2DArray;
     InvA : TReal2DArray;
     IsUpper : Boolean;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : MatInvReport):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
begin
    A := DynamicArrayCopy(A);
    InvA := DynamicArrayCopy(InvA);
    I:=0;
    while I<=N-2 do
    begin
        if IsUpper then
        begin
            for i_ := I+1 to N-1 do
            begin
                A[i_,I] := A[I,i_];
            end;
            for i_ := I+1 to N-1 do
            begin
                InvA[i_,I] := InvA[I,i_];
            end;
        end
        else
        begin
            for i_ := I+1 to N-1 do
            begin
                A[I,i_] := A[i_,I];
            end;
            for i_ := I+1 to N-1 do
            begin
                InvA[I,i_] := InvA[i_,I];
            end;
        end;
        Inc(I);
    end;
    Result := True;
    if Info<=0 then
    begin
        Result := False;
    end
    else
    begin
        Result := Result and  not (AP_FP_Less(Rep.R1,100*MachineEpsilon) or AP_FP_Greater(Rep.R1,1+1000*MachineEpsilon));
        Result := Result and  not (AP_FP_Less(Rep.RInf,100*MachineEpsilon) or AP_FP_Greater(Rep.RInf,1+1000*MachineEpsilon));
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                V := 0.0;
                for i_ := 0 to N-1 do
                begin
                    V := V + A[I,i_]*InvA[i_,J];
                end;
                if I=J then
                begin
                    V := V-1;
                end;
                Result := Result and AP_FP_Less_Eq(AbsReal(V),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Checks whether inverse is correct
Returns True on success.
*************************************************************************)
function HPDMatrixCheckInverse(A : TComplex2DArray;
     InvA : TComplex2DArray;
     IsUpper : Boolean;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : MatInvReport):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Complex;
    i_ : AlglibInteger;
begin
    A := DynamicArrayCopy(A);
    InvA := DynamicArrayCopy(InvA);
    I:=0;
    while I<=N-2 do
    begin
        if IsUpper then
        begin
            for i_ := I+1 to N-1 do
            begin
                A[i_,I] := Conj(A[I,i_]);
            end;
            for i_ := I+1 to N-1 do
            begin
                InvA[i_,I] := Conj(InvA[I,i_]);
            end;
        end
        else
        begin
            for i_ := I+1 to N-1 do
            begin
                A[I,i_] := Conj(A[i_,I]);
            end;
            for i_ := I+1 to N-1 do
            begin
                InvA[I,i_] := Conj(InvA[i_,I]);
            end;
        end;
        Inc(I);
    end;
    Result := True;
    if Info<=0 then
    begin
        Result := False;
    end
    else
    begin
        Result := Result and  not (AP_FP_Less(Rep.R1,100*MachineEpsilon) or AP_FP_Greater(Rep.R1,1+1000*MachineEpsilon));
        Result := Result and  not (AP_FP_Less(Rep.RInf,100*MachineEpsilon) or AP_FP_Greater(Rep.RInf,1+1000*MachineEpsilon));
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                V := C_Complex(0.0);
                for i_ := 0 to N-1 do
                begin
                    V := C_Add(V,C_Mul(A[I,i_],InvA[i_,J]));
                end;
                if I=J then
                begin
                    V := C_SubR(V,1);
                end;
                Result := Result and AP_FP_Less_Eq(AbsComplex(V),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Checks whether inversion result indicate singular matrix
Returns True on success.
*************************************************************************)
function RMatrixCheckInverseSingular(const InvA : TReal2DArray;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : MatInvReport):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Result := True;
    if (Info<>-3) and (Info<>1) then
    begin
        Result := False;
    end
    else
    begin
        Result := Result and  not (AP_FP_Less(Rep.R1,0) or AP_FP_Greater(Rep.R1,1000*MachineEpsilon));
        Result := Result and  not (AP_FP_Less(Rep.RInf,0) or AP_FP_Greater(Rep.RInf,1000*MachineEpsilon));
        if Info=-3 then
        begin
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    Result := Result and AP_FP_Eq(InvA[I,J],0);
                    Inc(J);
                end;
                Inc(I);
            end;
        end;
    end;
end;


(*************************************************************************
Checks whether inverse is correct
Returns True on success.
*************************************************************************)
function CMatrixCheckInverse(const A : TComplex2DArray;
     const InvA : TComplex2DArray;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : MatInvReport):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Complex;
    i_ : AlglibInteger;
begin
    Result := True;
    if Info<=0 then
    begin
        Result := False;
    end
    else
    begin
        Result := Result and  not (AP_FP_Less(Rep.R1,100*MachineEpsilon) or AP_FP_Greater(Rep.R1,1+1000*MachineEpsilon));
        Result := Result and  not (AP_FP_Less(Rep.RInf,100*MachineEpsilon) or AP_FP_Greater(Rep.RInf,1+1000*MachineEpsilon));
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                V := C_Complex(0.0);
                for i_ := 0 to N-1 do
                begin
                    V := C_Add(V,C_Mul(A[I,i_],InvA[i_,J]));
                end;
                if I=J then
                begin
                    V := C_SubR(V,1);
                end;
                Result := Result and AP_FP_Less_Eq(AbsComplex(V),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Checks whether inversion result indicate singular matrix
Returns True on success.
*************************************************************************)
function CMatrixCheckInverseSingular(const InvA : TComplex2DArray;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : MatInvReport):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Result := True;
    if (Info<>-3) and (Info<>1) then
    begin
        Result := False;
    end
    else
    begin
        Result := Result and  not (AP_FP_Less(Rep.R1,0) or AP_FP_Greater(Rep.R1,1000*MachineEpsilon));
        Result := Result and  not (AP_FP_Less(Rep.RInf,0) or AP_FP_Greater(Rep.RInf,1000*MachineEpsilon));
        if Info=-3 then
        begin
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    Result := Result and C_EqualR(InvA[I,J],0);
                    Inc(J);
                end;
                Inc(I);
            end;
        end;
    end;
end;


(*************************************************************************
Drops upper or lower half of the matrix - fills it by special pattern
which may be used later to ensure that this part wasn't changed
*************************************************************************)
procedure RMatrixDropHalf(var A : TReal2DArray;
     N : AlglibInteger;
     DropLower : Boolean);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if DropLower and (I>J) or  not DropLower and (I<J) then
            begin
                A[I,J] := 1+2*I+3*J;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Drops upper or lower half of the matrix - fills it by special pattern
which may be used later to ensure that this part wasn't changed
*************************************************************************)
procedure CMatrixDropHalf(var A : TComplex2DArray;
     N : AlglibInteger;
     DropLower : Boolean);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if DropLower and (I>J) or  not DropLower and (I<J) then
            begin
                A[I,J] := C_Complex(1+2*I+3*J);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Real TR inverse
*************************************************************************)
procedure TestRTRInv(MaxN : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var RTRErrors : Boolean);
var
    A : TReal2DArray;
    B : TReal2DArray;
    N : AlglibInteger;
    Pass : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    Task : AlglibInteger;
    IsUpper : Boolean;
    IsUnit : Boolean;
    V : Double;
    WasErrors : Boolean;
    Info : AlglibInteger;
    Rep : MatInvReport;
    i_ : AlglibInteger;
begin
    WasErrors := False;
    
    //
    // Test
    //
    N:=1;
    while N<=MaxN do
    begin
        SetLength(A, N, N);
        SetLength(B, N, N);
        Task:=0;
        while Task<=3 do
        begin
            Pass:=1;
            while Pass<=PassCount do
            begin
                
                //
                // Determine task
                //
                IsUpper := Task mod 2=0;
                IsUnit := Task div 2 mod 2=0;
                
                //
                // Generate matrix
                //
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        if I=J then
                        begin
                            A[I,I] := 1+RandomReal;
                        end
                        else
                        begin
                            A[I,J] := 0.2*RandomReal-0.1;
                        end;
                        B[I,J] := A[I,J];
                        Inc(J);
                    end;
                    Inc(I);
                end;
                
                //
                // Inverse
                //
                RMatrixTRInverse(B, N, IsUpper, IsUnit, Info, Rep);
                if Info<=0 then
                begin
                    RTRErrors := True;
                    Exit;
                end;
                
                //
                // Structural test
                //
                if IsUnit then
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        RTRErrors := RTRErrors or AP_FP_Neq(A[I,I],B[I,I]);
                        Inc(I);
                    end;
                end;
                if IsUpper then
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=I-1 do
                        begin
                            RTRErrors := RTRErrors or AP_FP_Neq(A[I,J],B[I,J]);
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                end
                else
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=I+1;
                        while J<=N-1 do
                        begin
                            RTRErrors := RTRErrors or AP_FP_Neq(A[I,J],B[I,J]);
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                end;
                
                //
                // Inverse test
                //
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        if (J<I) and IsUpper or (J>I) and  not IsUpper then
                        begin
                            A[I,J] := 0;
                            B[I,J] := 0;
                        end;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                if IsUnit then
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        A[I,I] := 1;
                        B[I,I] := 1;
                        Inc(I);
                    end;
                end;
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        V := 0.0;
                        for i_ := 0 to N-1 do
                        begin
                            V := V + A[I,i_]*B[i_,J];
                        end;
                        if J<>I then
                        begin
                            RTRErrors := RTRErrors or AP_FP_Greater(AbsReal(V),Threshold);
                        end
                        else
                        begin
                            RTRErrors := RTRErrors or AP_FP_Greater(AbsReal(V-1),Threshold);
                        end;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                Inc(Pass);
            end;
            Inc(Task);
        end;
        Inc(N);
    end;
end;


(*************************************************************************
Complex TR inverse
*************************************************************************)
procedure TestCTRInv(MaxN : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var CTRErrors : Boolean);
var
    A : TComplex2DArray;
    B : TComplex2DArray;
    N : AlglibInteger;
    Pass : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    Task : AlglibInteger;
    IsUpper : Boolean;
    IsUnit : Boolean;
    V : Complex;
    WasErrors : Boolean;
    Info : AlglibInteger;
    Rep : MatInvReport;
    i_ : AlglibInteger;
begin
    WasErrors := False;
    
    //
    // Test
    //
    N:=1;
    while N<=MaxN do
    begin
        SetLength(A, N, N);
        SetLength(B, N, N);
        Task:=0;
        while Task<=3 do
        begin
            Pass:=1;
            while Pass<=PassCount do
            begin
                
                //
                // Determine task
                //
                IsUpper := Task mod 2=0;
                IsUnit := Task div 2 mod 2=0;
                
                //
                // Generate matrix
                //
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        if I=J then
                        begin
                            A[I,I].X := 1+RandomReal;
                            A[I,I].Y := 1+RandomReal;
                        end
                        else
                        begin
                            A[I,J].X := 0.2*RandomReal-0.1;
                            A[I,J].Y := 0.2*RandomReal-0.1;
                        end;
                        B[I,J] := A[I,J];
                        Inc(J);
                    end;
                    Inc(I);
                end;
                
                //
                // Inverse
                //
                CMatrixTRInverse(B, N, IsUpper, IsUnit, Info, Rep);
                if Info<=0 then
                begin
                    CTRErrors := True;
                    Exit;
                end;
                
                //
                // Structural test
                //
                if IsUnit then
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        CTRErrors := CTRErrors or C_NotEqual(A[I,I],B[I,I]);
                        Inc(I);
                    end;
                end;
                if IsUpper then
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=I-1 do
                        begin
                            CTRErrors := CTRErrors or C_NotEqual(A[I,J],B[I,J]);
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                end
                else
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=I+1;
                        while J<=N-1 do
                        begin
                            CTRErrors := CTRErrors or C_NotEqual(A[I,J],B[I,J]);
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                end;
                
                //
                // Inverse test
                //
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        if (J<I) and IsUpper or (J>I) and  not IsUpper then
                        begin
                            A[I,J] := C_Complex(0);
                            B[I,J] := C_Complex(0);
                        end;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                if IsUnit then
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        A[I,I] := C_Complex(1);
                        B[I,I] := C_Complex(1);
                        Inc(I);
                    end;
                end;
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        V := C_Complex(0.0);
                        for i_ := 0 to N-1 do
                        begin
                            V := C_Add(V,C_Mul(A[I,i_],B[i_,J]));
                        end;
                        if J<>I then
                        begin
                            CTRErrors := CTRErrors or AP_FP_Greater(AbsComplex(V),Threshold);
                        end
                        else
                        begin
                            CTRErrors := CTRErrors or AP_FP_Greater(AbsComplex(C_SubR(V,1)),Threshold);
                        end;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                Inc(Pass);
            end;
            Inc(Task);
        end;
        Inc(N);
    end;
end;


(*************************************************************************
Real test
*************************************************************************)
procedure TestRInv(MaxN : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var RErrors : Boolean);
var
    A : TReal2DArray;
    LUA : TReal2DArray;
    InvA : TReal2DArray;
    InvLUA : TReal2DArray;
    P : TInteger1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    N : AlglibInteger;
    Pass : AlglibInteger;
    TaskKind : AlglibInteger;
    V : Double;
    Info : AlglibInteger;
    Rep : MatInvReport;
    i_ : AlglibInteger;
begin
    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        N:=1;
        while N<=MaxN do
        begin
            
            //
            // ********************************************************
            // WELL CONDITIONED TASKS
            // ability to find correct solution is tested
            // ********************************************************
            //
            // 1. generate random well conditioned matrix A.
            // 2. generate random solution vector xe
            // 3. generate right part b=A*xe
            // 4. test different methods on original A
            //
            RMatrixRndCond(N, 1000, A);
            RMatrixMakeACopy(A, N, N, LUA);
            RMatrixLU(LUA, N, N, P);
            RMatrixMakeACopy(A, N, N, InvA);
            RMatrixMakeACopy(LUA, N, N, InvLUA);
            Info := 0;
            UnsetRep(Rep);
            RMatrixInverse(InvA, N, Info, Rep);
            RErrors := RErrors or  not RMatrixCheckInverse(A, InvA, N, Threshold, Info, Rep);
            Info := 0;
            UnsetRep(Rep);
            RMatrixLUInverse(InvLUA, P, N, Info, Rep);
            RErrors := RErrors or  not RMatrixCheckInverse(A, InvLUA, N, Threshold, Info, Rep);
            
            //
            // ********************************************************
            // EXACTLY SINGULAR MATRICES
            // ability to detect singularity is tested
            // ********************************************************
            //
            // 1. generate different types of singular matrices:
            //    * zero
            //    * with zero columns
            //    * with zero rows
            //    * with equal rows/columns
            // 2. test different methods
            //
            TaskKind:=0;
            while TaskKind<=4 do
            begin
                Unset2D(A);
                if TaskKind=0 then
                begin
                    
                    //
                    // all zeros
                    //
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J] := 0;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                end;
                if TaskKind=1 then
                begin
                    
                    //
                    // there is zero column
                    //
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J] := 2*RandomReal-1;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    K := RandomInteger(N);
                    for i_ := 0 to N-1 do
                    begin
                        A[i_,K] := 0*A[i_,K];
                    end;
                end;
                if TaskKind=2 then
                begin
                    
                    //
                    // there is zero row
                    //
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J] := 2*RandomReal-1;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    K := RandomInteger(N);
                    APVMul(@A[K][0], 0, N-1, 0);
                end;
                if TaskKind=3 then
                begin
                    
                    //
                    // equal columns
                    //
                    if N<2 then
                    begin
                        Inc(TaskKind);
                        Continue;
                    end;
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J] := 2*RandomReal-1;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    K := 1+RandomInteger(N-1);
                    for i_ := 0 to N-1 do
                    begin
                        A[i_,0] := A[i_,K];
                    end;
                end;
                if TaskKind=4 then
                begin
                    
                    //
                    // equal rows
                    //
                    if N<2 then
                    begin
                        Inc(TaskKind);
                        Continue;
                    end;
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J] := 2*RandomReal-1;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    K := 1+RandomInteger(N-1);
                    APVMove(@A[0][0], 0, N-1, @A[K][0], 0, N-1);
                end;
                RMatrixMakeACopy(A, N, N, LUA);
                RMatrixLU(LUA, N, N, P);
                Info := 0;
                UnsetRep(Rep);
                RMatrixInverse(A, N, Info, Rep);
                RErrors := RErrors or  not RMatrixCheckInverseSingular(A, N, Threshold, Info, Rep);
                Info := 0;
                UnsetRep(Rep);
                RMatrixLUInverse(LUA, P, N, Info, Rep);
                RErrors := RErrors or  not RMatrixCheckInverseSingular(LUA, N, Threshold, Info, Rep);
                Inc(TaskKind);
            end;
            Inc(N);
        end;
        Inc(Pass);
    end;
end;


(*************************************************************************
Complex test
*************************************************************************)
procedure TestCInv(MaxN : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var CErrors : Boolean);
var
    A : TComplex2DArray;
    LUA : TComplex2DArray;
    InvA : TComplex2DArray;
    InvLUA : TComplex2DArray;
    P : TInteger1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    N : AlglibInteger;
    Pass : AlglibInteger;
    TaskKind : AlglibInteger;
    V : Double;
    Info : AlglibInteger;
    Rep : MatInvReport;
    i_ : AlglibInteger;
begin
    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        N:=1;
        while N<=MaxN do
        begin
            
            //
            // ********************************************************
            // WELL CONDITIONED TASKS
            // ability to find correct solution is tested
            // ********************************************************
            //
            // 1. generate random well conditioned matrix A.
            // 2. generate random solution vector xe
            // 3. generate right part b=A*xe
            // 4. test different methods on original A
            //
            CMatrixRndCond(N, 1000, A);
            CMatrixMakeACopy(A, N, N, LUA);
            CMatrixLU(LUA, N, N, P);
            CMatrixMakeACopy(A, N, N, InvA);
            CMatrixMakeACopy(LUA, N, N, InvLUA);
            Info := 0;
            UnsetRep(Rep);
            CMatrixInverse(InvA, N, Info, Rep);
            CErrors := CErrors or  not CMatrixCheckInverse(A, InvA, N, Threshold, Info, Rep);
            Info := 0;
            UnsetRep(Rep);
            CMatrixLUInverse(InvLUA, P, N, Info, Rep);
            CErrors := CErrors or  not CMatrixCheckInverse(A, InvLUA, N, Threshold, Info, Rep);
            
            //
            // ********************************************************
            // EXACTLY SINGULAR MATRICES
            // ability to detect singularity is tested
            // ********************************************************
            //
            // 1. generate different types of singular matrices:
            //    * zero
            //    * with zero columns
            //    * with zero rows
            //    * with equal rows/columns
            // 2. test different methods
            //
            TaskKind:=0;
            while TaskKind<=4 do
            begin
                CUnset2D(A);
                if TaskKind=0 then
                begin
                    
                    //
                    // all zeros
                    //
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J] := C_Complex(0);
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                end;
                if TaskKind=1 then
                begin
                    
                    //
                    // there is zero column
                    //
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J].X := 2*RandomReal-1;
                            A[I,J].Y := 2*RandomReal-1;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    K := RandomInteger(N);
                    for i_ := 0 to N-1 do
                    begin
                        A[i_,K] := C_MulR(A[i_,K],0);
                    end;
                end;
                if TaskKind=2 then
                begin
                    
                    //
                    // there is zero row
                    //
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J].X := 2*RandomReal-1;
                            A[I,J].Y := 2*RandomReal-1;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    K := RandomInteger(N);
                    for i_ := 0 to N-1 do
                    begin
                        A[K,i_] := C_MulR(A[K,i_],0);
                    end;
                end;
                if TaskKind=3 then
                begin
                    
                    //
                    // equal columns
                    //
                    if N<2 then
                    begin
                        Inc(TaskKind);
                        Continue;
                    end;
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J].X := 2*RandomReal-1;
                            A[I,J].Y := 2*RandomReal-1;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    K := 1+RandomInteger(N-1);
                    for i_ := 0 to N-1 do
                    begin
                        A[i_,0] := A[i_,K];
                    end;
                end;
                if TaskKind=4 then
                begin
                    
                    //
                    // equal rows
                    //
                    if N<2 then
                    begin
                        Inc(TaskKind);
                        Continue;
                    end;
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J].X := 2*RandomReal-1;
                            A[I,J].Y := 2*RandomReal-1;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    K := 1+RandomInteger(N-1);
                    for i_ := 0 to N-1 do
                    begin
                        A[0,i_] := A[K,i_];
                    end;
                end;
                CMatrixMakeACopy(A, N, N, LUA);
                CMatrixLU(LUA, N, N, P);
                Info := 0;
                UnsetRep(Rep);
                CMatrixInverse(A, N, Info, Rep);
                CErrors := CErrors or  not CMatrixCheckInverseSingular(A, N, Threshold, Info, Rep);
                Info := 0;
                UnsetRep(Rep);
                CMatrixLUInverse(LUA, P, N, Info, Rep);
                CErrors := CErrors or  not CMatrixCheckInverseSingular(LUA, N, Threshold, Info, Rep);
                Inc(TaskKind);
            end;
            Inc(N);
        end;
        Inc(Pass);
    end;
end;


(*************************************************************************
SPD test
*************************************************************************)
procedure TestSPDInv(MaxN : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var SPDErrors : Boolean);
var
    A : TReal2DArray;
    CHA : TReal2DArray;
    InvA : TReal2DArray;
    InvCHA : TReal2DArray;
    IsUpper : Boolean;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    N : AlglibInteger;
    Pass : AlglibInteger;
    TaskKind : AlglibInteger;
    V : Double;
    Info : AlglibInteger;
    Rep : MatInvReport;
    i_ : AlglibInteger;
begin
    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        N:=1;
        while N<=MaxN do
        begin
            IsUpper := AP_FP_Greater(RandomReal,0.5);
            
            //
            // ********************************************************
            // WELL CONDITIONED TASKS
            // ability to find correct solution is tested
            // ********************************************************
            //
            // 1. generate random well conditioned matrix A.
            // 2. generate random solution vector xe
            // 3. generate right part b=A*xe
            // 4. test different methods on original A
            //
            SPDMatrixRndCond(N, 1000, A);
            RMatrixDropHalf(A, N, IsUpper);
            RMatrixMakeACopy(A, N, N, CHA);
            if  not SPDMatrixCholesky(CHA, N, IsUpper) then
            begin
                Inc(N);
                Continue;
            end;
            RMatrixMakeACopy(A, N, N, InvA);
            RMatrixMakeACopy(CHA, N, N, InvCHA);
            Info := 0;
            UnsetRep(Rep);
            SPDMatrixInverse(InvA, N, IsUpper, Info, Rep);
            SPDErrors := SPDErrors or  not SPDMatrixCheckInverse(A, InvA, IsUpper, N, Threshold, Info, Rep);
            Info := 0;
            UnsetRep(Rep);
            SPDMatrixCholeskyInverse(InvCHA, N, IsUpper, Info, Rep);
            SPDErrors := SPDErrors or  not SPDMatrixCheckInverse(A, InvCHA, IsUpper, N, Threshold, Info, Rep);
            
            //
            // ********************************************************
            // EXACTLY SINGULAR MATRICES
            // ability to detect singularity is tested
            // ********************************************************
            //
            // 1. generate different types of singular matrices:
            //    * zero
            //    * with zero columns
            //    * with zero rows
            // 2. test different methods
            //
            TaskKind:=0;
            while TaskKind<=2 do
            begin
                Unset2D(A);
                if TaskKind=0 then
                begin
                    
                    //
                    // all zeros
                    //
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J] := 0;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                end;
                if TaskKind=1 then
                begin
                    
                    //
                    // there is zero column
                    //
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J] := 2*RandomReal-1;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    K := RandomInteger(N);
                    for i_ := 0 to N-1 do
                    begin
                        A[i_,K] := 0*A[i_,K];
                    end;
                end;
                if TaskKind=2 then
                begin
                    
                    //
                    // there is zero row
                    //
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J] := 2*RandomReal-1;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    K := RandomInteger(N);
                    APVMul(@A[K][0], 0, N-1, 0);
                end;
                Info := 0;
                UnsetRep(Rep);
                SPDMatrixCholeskyInverse(A, N, IsUpper, Info, Rep);
                if (Info<>-3) and (Info<>1) then
                begin
                    SPDErrors := True;
                end
                else
                begin
                    SPDErrors := SPDErrors or AP_FP_Less(Rep.R1,0) or AP_FP_Greater(Rep.R1,1000*MachineEpsilon);
                    SPDErrors := SPDErrors or AP_FP_Less(Rep.RInf,0) or AP_FP_Greater(Rep.RInf,1000*MachineEpsilon);
                end;
                Inc(TaskKind);
            end;
            Inc(N);
        end;
        Inc(Pass);
    end;
end;


(*************************************************************************
HPD test
*************************************************************************)
procedure TestHPDInv(MaxN : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var HPDErrors : Boolean);
var
    A : TComplex2DArray;
    CHA : TComplex2DArray;
    InvA : TComplex2DArray;
    InvCHA : TComplex2DArray;
    IsUpper : Boolean;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    N : AlglibInteger;
    Pass : AlglibInteger;
    TaskKind : AlglibInteger;
    V : Complex;
    Info : AlglibInteger;
    Rep : MatInvReport;
    i_ : AlglibInteger;
begin
    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        N:=1;
        while N<=MaxN do
        begin
            IsUpper := AP_FP_Greater(RandomReal,0.5);
            
            //
            // ********************************************************
            // WELL CONDITIONED TASKS
            // ability to find correct solution is tested
            // ********************************************************
            //
            // 1. generate random well conditioned matrix A.
            // 2. generate random solution vector xe
            // 3. generate right part b=A*xe
            // 4. test different methods on original A
            //
            HPDMatrixRndCond(N, 1000, A);
            CMatrixDropHalf(A, N, IsUpper);
            CMatrixMakeACopy(A, N, N, CHA);
            if  not HPDMatrixCholesky(CHA, N, IsUpper) then
            begin
                Inc(N);
                Continue;
            end;
            CMatrixMakeACopy(A, N, N, InvA);
            CMatrixMakeACopy(CHA, N, N, InvCHA);
            Info := 0;
            UnsetRep(Rep);
            HPDMatrixInverse(InvA, N, IsUpper, Info, Rep);
            HPDErrors := HPDErrors or  not HPDMatrixCheckInverse(A, InvA, IsUpper, N, Threshold, Info, Rep);
            Info := 0;
            UnsetRep(Rep);
            HPDMatrixCholeskyInverse(InvCHA, N, IsUpper, Info, Rep);
            HPDErrors := HPDErrors or  not HPDMatrixCheckInverse(A, InvCHA, IsUpper, N, Threshold, Info, Rep);
            
            //
            // ********************************************************
            // EXACTLY SINGULAR MATRICES
            // ability to detect singularity is tested
            // ********************************************************
            //
            // 1. generate different types of singular matrices:
            //    * zero
            //    * with zero columns
            //    * with zero rows
            // 2. test different methods
            //
            TaskKind:=0;
            while TaskKind<=2 do
            begin
                CUnset2D(A);
                if TaskKind=0 then
                begin
                    
                    //
                    // all zeros
                    //
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J] := C_Complex(0);
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                end;
                if TaskKind=1 then
                begin
                    
                    //
                    // there is zero column
                    //
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J].X := 2*RandomReal-1;
                            A[I,J].Y := 2*RandomReal-1;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    K := RandomInteger(N);
                    for i_ := 0 to N-1 do
                    begin
                        A[i_,K] := C_MulR(A[i_,K],0);
                    end;
                    for i_ := 0 to N-1 do
                    begin
                        A[K,i_] := C_MulR(A[K,i_],0);
                    end;
                end;
                if TaskKind=2 then
                begin
                    
                    //
                    // there is zero row
                    //
                    SetLength(A, N, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J].X := 2*RandomReal-1;
                            A[I,J].Y := 2*RandomReal-1;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    K := RandomInteger(N);
                    for i_ := 0 to N-1 do
                    begin
                        A[K,i_] := C_MulR(A[K,i_],0);
                    end;
                    for i_ := 0 to N-1 do
                    begin
                        A[i_,K] := C_MulR(A[i_,K],0);
                    end;
                end;
                Info := 0;
                UnsetRep(Rep);
                HPDMatrixCholeskyInverse(A, N, IsUpper, Info, Rep);
                if (Info<>-3) and (Info<>1) then
                begin
                    HPDErrors := True;
                end
                else
                begin
                    HPDErrors := HPDErrors or AP_FP_Less(Rep.R1,0) or AP_FP_Greater(Rep.R1,1000*MachineEpsilon);
                    HPDErrors := HPDErrors or AP_FP_Less(Rep.RInf,0) or AP_FP_Greater(Rep.RInf,1000*MachineEpsilon);
                end;
                Inc(TaskKind);
            end;
            Inc(N);
        end;
        Inc(Pass);
    end;
end;


(*************************************************************************
Unsets real matrix
*************************************************************************)
procedure Unset2D(var X : TReal2DArray);
begin
    SetLength(X, 1, 1);
    X[0,0] := 2*RandomReal-1;
end;


(*************************************************************************
Unsets real vector
*************************************************************************)
procedure Unset1D(var X : TReal1DArray);
begin
    SetLength(X, 1);
    X[0] := 2*RandomReal-1;
end;


(*************************************************************************
Unsets real matrix
*************************************************************************)
procedure CUnset2D(var X : TComplex2DArray);
begin
    SetLength(X, 1, 1);
    X[0,0] := C_Complex(2*RandomReal-1);
end;


(*************************************************************************
Unsets real vector
*************************************************************************)
procedure CUnset1D(var X : TComplex1DArray);
begin
    SetLength(X, 1);
    X[0] := C_Complex(2*RandomReal-1);
end;


(*************************************************************************
Unsets report
*************************************************************************)
procedure UnsetRep(var R : MatInvReport);
begin
    R.R1 := -1;
    R.RInf := -1;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testmatinvunit_test_silent():Boolean;
begin
    Result := TestMatInv(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testmatinvunit_test():Boolean;
begin
    Result := TestMatInv(False);
end;


end.