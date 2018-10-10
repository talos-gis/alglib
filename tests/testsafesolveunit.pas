unit testsafesolveunit;
interface
uses Math, Sysutils, Ap, safesolve;

function TestSafeSolve(Silent : Boolean):Boolean;
function testsafesolveunit_test_silent():Boolean;
function testsafesolveunit_test():Boolean;

implementation

procedure RMatrixMakeACopy(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);forward;
procedure CMatrixMakeACopy(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TComplex2DArray);forward;


(*************************************************************************
Main unittest subroutine
*************************************************************************)
function TestSafeSolve(Silent : Boolean):Boolean;
var
    MaxMN : AlglibInteger;
    Threshold : Double;
    RErrors : Boolean;
    CErrors : Boolean;
    WasErrors : Boolean;
    IsUpper : Boolean;
    Trans : AlglibInteger;
    IsUnit : Boolean;
    ScaleA : Double;
    Growth : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    N : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
    CV : Complex;
    CA : TComplex2DArray;
    CEA : TComplex2DArray;
    CTmpA : TComplex2DArray;
    CXS : TComplex1DArray;
    CXE : TComplex1DArray;
    RV : Double;
    RA : TReal2DArray;
    REA : TReal2DArray;
    RTmpA : TReal2DArray;
    RXS : TReal1DArray;
    RXE : TReal1DArray;
    i_ : AlglibInteger;
begin
    MaxMN := 30;
    Threshold := 100000*MachineEpsilon;
    RErrors := False;
    CErrors := False;
    WasErrors := False;
    
    //
    // Different problems: general tests
    //
    N:=1;
    while N<=MaxMN do
    begin
        
        //
        // test complex solver with well-conditioned matrix:
        // 1. generate A: fill off-diagonal elements with small values,
        //    diagonal elements are filled with larger values
        // 2. generate 'effective' A
        // 3. prepare task (exact X is stored in CXE, right part - in CXS),
        //    solve and compare CXS and CXE
        //
        IsUpper := AP_FP_Greater(RandomReal,0.5);
        Trans := RandomInteger(3);
        IsUnit := AP_FP_Greater(RandomReal,0.5);
        ScaleA := RandomReal+0.5;
        SetLength(CA, N, N);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                if I=J then
                begin
                    CA[I,J].X := (2*RandomInteger(2)-1)*(5+RandomReal);
                    CA[I,J].Y := (2*RandomInteger(2)-1)*(5+RandomReal);
                end
                else
                begin
                    CA[I,J].X := 0.2*RandomReal-0.1;
                    CA[I,J].Y := 0.2*RandomReal-0.1;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        CMatrixMakeACopy(CA, N, N, CTmpA);
        I:=0;
        while I<=N-1 do
        begin
            if IsUpper then
            begin
                J1 := 0;
                J2 := I-1;
            end
            else
            begin
                J1 := I+1;
                J2 := N-1;
            end;
            J:=J1;
            while J<=J2 do
            begin
                CTmpA[I,J] := C_Complex(0);
                Inc(J);
            end;
            if IsUnit then
            begin
                CTmpA[I,I] := C_Complex(1);
            end;
            Inc(I);
        end;
        SetLength(CEA, N, N);
        I:=0;
        while I<=N-1 do
        begin
            if Trans=0 then
            begin
                for i_ := 0 to N-1 do
                begin
                    CEA[I,i_] := C_MulR(CTmpA[I,i_],ScaleA);
                end;
            end;
            if Trans=1 then
            begin
                for i_ := 0 to N-1 do
                begin
                    CEA[i_,I] := C_MulR(CTmpA[I,i_],ScaleA);
                end;
            end;
            if Trans=2 then
            begin
                for i_ := 0 to N-1 do
                begin
                    CEA[i_,I] := C_MulR(Conj(CTmpA[I,i_]),ScaleA);
                end;
            end;
            Inc(I);
        end;
        SetLength(CXE, N);
        I:=0;
        while I<=N-1 do
        begin
            CXE[I].X := 2*RandomReal-1;
            CXE[I].Y := 2*RandomReal-1;
            Inc(I);
        end;
        SetLength(CXS, N);
        I:=0;
        while I<=N-1 do
        begin
            CV := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                CV := C_Add(CV,C_Mul(CEA[I,i_],CXE[i_]));
            end;
            CXS[I] := CV;
            Inc(I);
        end;
        if CMatrixScaledTRSafeSolve(CA, ScaleA, N, CXS, IsUpper, Trans, IsUnit, Sqrt(MaxRealNumber)) then
        begin
            I:=0;
            while I<=N-1 do
            begin
                CErrors := CErrors or AP_FP_Greater(AbsComplex(C_Sub(CXS[I],CXE[I])),Threshold);
                Inc(I);
            end;
        end
        else
        begin
            CErrors := True;
        end;
        
        //
        // same with real
        //
        IsUpper := AP_FP_Greater(RandomReal,0.5);
        Trans := RandomInteger(2);
        IsUnit := AP_FP_Greater(RandomReal,0.5);
        ScaleA := RandomReal+0.5;
        SetLength(RA, N, N);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                if I=J then
                begin
                    RA[I,J] := (2*RandomInteger(2)-1)*(5+RandomReal);
                end
                else
                begin
                    RA[I,J] := 0.2*RandomReal-0.1;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        RMatrixMakeACopy(RA, N, N, RTmpA);
        I:=0;
        while I<=N-1 do
        begin
            if IsUpper then
            begin
                J1 := 0;
                J2 := I-1;
            end
            else
            begin
                J1 := I+1;
                J2 := N-1;
            end;
            J:=J1;
            while J<=J2 do
            begin
                RTmpA[I,J] := 0;
                Inc(J);
            end;
            if IsUnit then
            begin
                RTmpA[I,I] := 1;
            end;
            Inc(I);
        end;
        SetLength(REA, N, N);
        I:=0;
        while I<=N-1 do
        begin
            if Trans=0 then
            begin
                APVMove(@REA[I][0], 0, N-1, @RTmpA[I][0], 0, N-1, ScaleA);
            end;
            if Trans=1 then
            begin
                for i_ := 0 to N-1 do
                begin
                    REA[i_,I] := ScaleA*RTmpA[I,i_];
                end;
            end;
            Inc(I);
        end;
        SetLength(RXE, N);
        I:=0;
        while I<=N-1 do
        begin
            RXE[I] := 2*RandomReal-1;
            Inc(I);
        end;
        SetLength(RXS, N);
        I:=0;
        while I<=N-1 do
        begin
            RV := APVDotProduct(@REA[I][0], 0, N-1, @RXE[0], 0, N-1);
            RXS[I] := RV;
            Inc(I);
        end;
        if RMatrixScaledTRSafeSolve(RA, ScaleA, N, RXS, IsUpper, Trans, IsUnit, Sqrt(MaxRealNumber)) then
        begin
            I:=0;
            while I<=N-1 do
            begin
                RErrors := RErrors or AP_FP_Greater(AbsReal(RXS[I]-RXE[I]),Threshold);
                Inc(I);
            end;
        end
        else
        begin
            RErrors := True;
        end;
        Inc(N);
    end;
    
    //
    // Special test with diagonal ill-conditioned matrix:
    // * ability to solve it when resulting growth is less than threshold
    // * ability to stop solve when resulting growth is greater than threshold
    //
    // A = diag(1, 1/growth)
    // b = (1, 0.5)
    //
    N := 2;
    Growth := 10;
    SetLength(CA, N, N);
    CA[0,0] := C_Complex(1);
    CA[0,1] := C_Complex(0);
    CA[1,0] := C_Complex(0);
    CA[1,1] := C_Complex(1/Growth);
    SetLength(CXS, N);
    CXS[0] := C_Complex(1.0);
    CXS[1] := C_Complex(0.5);
    CErrors := CErrors or  not CMatrixScaledTRSafeSolve(CA, 1.0, N, CXS, AP_FP_Greater(RandomReal,0.5), RandomInteger(3), False, 1.05*Max(AbsComplex(CXS[1])*Growth, 1.0));
    CErrors := CErrors or  not CMatrixScaledTRSafeSolve(CA, 1.0, N, CXS, AP_FP_Greater(RandomReal,0.5), RandomInteger(3), False, 0.95*Max(AbsComplex(CXS[1])*Growth, 1.0));
    SetLength(RA, N, N);
    RA[0,0] := 1;
    RA[0,1] := 0;
    RA[1,0] := 0;
    RA[1,1] := 1/Growth;
    SetLength(RXS, N);
    RXS[0] := 1.0;
    RXS[1] := 0.5;
    RErrors := RErrors or  not RMatrixScaledTRSafeSolve(RA, 1.0, N, RXS, AP_FP_Greater(RandomReal,0.5), RandomInteger(2), False, 1.05*Max(AbsReal(RXS[1])*Growth, 1.0));
    RErrors := RErrors or  not RMatrixScaledTRSafeSolve(RA, 1.0, N, RXS, AP_FP_Greater(RandomReal,0.5), RandomInteger(2), False, 0.95*Max(AbsReal(RXS[1])*Growth, 1.0));
    
    //
    // Special test with diagonal degenerate matrix:
    // * ability to solve it when resulting growth is less than threshold
    // * ability to stop solve when resulting growth is greater than threshold
    //
    // A = diag(1, 0)
    // b = (1, 0.5)
    //
    N := 2;
    SetLength(CA, N, N);
    CA[0,0] := C_Complex(1);
    CA[0,1] := C_Complex(0);
    CA[1,0] := C_Complex(0);
    CA[1,1] := C_Complex(0);
    SetLength(CXS, N);
    CXS[0] := C_Complex(1.0);
    CXS[1] := C_Complex(0.5);
    CErrors := CErrors or CMatrixScaledTRSafeSolve(CA, 1.0, N, CXS, AP_FP_Greater(RandomReal,0.5), RandomInteger(3), False, Sqrt(MaxRealNumber));
    SetLength(RA, N, N);
    RA[0,0] := 1;
    RA[0,1] := 0;
    RA[1,0] := 0;
    RA[1,1] := 0;
    SetLength(RXS, N);
    RXS[0] := 1.0;
    RXS[1] := 0.5;
    RErrors := RErrors or RMatrixScaledTRSafeSolve(RA, 1.0, N, RXS, AP_FP_Greater(RandomReal,0.5), RandomInteger(2), False, Sqrt(MaxRealNumber));
    
    //
    // report
    //
    WasErrors := RErrors or CErrors;
    if  not Silent then
    begin
        Write(Format('TESTING SAFE TR SOLVER'#13#10'',[]));
        Write(Format('REAL:                                    ',[]));
        if  not RErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('COMPLEX:                                 ',[]));
        if  not CErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        if WasErrors then
        begin
            Write(Format('TEST FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('TEST PASSED'#13#10'',[]));
        end;
        Write(Format(''#13#10''#13#10'',[]));
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
Silent unit test
*************************************************************************)
function testsafesolveunit_test_silent():Boolean;
begin
    Result := TestSafeSolve(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testsafesolveunit_test():Boolean;
begin
    Result := TestSafeSolve(False);
end;


end.