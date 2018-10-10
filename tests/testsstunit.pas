unit testsstunit;
interface
uses Math, Sysutils, Ap, trlinsolve;

function TestSST(Silent : Boolean):Boolean;
function testsstunit_test_silent():Boolean;
function testsstunit_test():Boolean;

implementation

procedure MakeACopy(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);forward;


(*************************************************************************
Main unittest subroutine
*************************************************************************)
function TestSST(Silent : Boolean):Boolean;
var
    MaxMN : AlglibInteger;
    PassCount : AlglibInteger;
    Threshold : Double;
    AEffective : TReal2DArray;
    AParam : TReal2DArray;
    XE : TReal1DArray;
    B : TReal1DArray;
    N : AlglibInteger;
    Pass : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    CntS : AlglibInteger;
    CntU : AlglibInteger;
    CntT : AlglibInteger;
    CntM : AlglibInteger;
    WasErrors : Boolean;
    IsUpper : Boolean;
    IsTrans : Boolean;
    Isunit : Boolean;
    V : Double;
    S : Double;
begin
    WasErrors := False;
    MaxMN := 15;
    PassCount := 15;
    Threshold := 1000*MachineEpsilon;
    
    //
    // Different problems
    //
    N:=1;
    while N<=MaxMN do
    begin
        SetLength(AEffective, N-1+1, N-1+1);
        SetLength(AParam, N-1+1, N-1+1);
        SetLength(XE, N-1+1);
        SetLength(B, N-1+1);
        Pass:=1;
        while Pass<=PassCount do
        begin
            CntS:=0;
            while CntS<=1 do
            begin
                CntU:=0;
                while CntU<=1 do
                begin
                    CntT:=0;
                    while CntT<=1 do
                    begin
                        CntM:=0;
                        while CntM<=2 do
                        begin
                            IsUpper := CntS=0;
                            Isunit := CntU=0;
                            IsTrans := CntT=0;
                            
                            //
                            // Skip meaningless combinations of parameters:
                            // (matrix is singular) AND (matrix is unit diagonal)
                            //
                            if (CntM=2) and Isunit then
                            begin
                                Inc(CntM);
                                Continue;
                            end;
                            
                            //
                            // Clear matrices
                            //
                            I:=0;
                            while I<=N-1 do
                            begin
                                J:=0;
                                while J<=N-1 do
                                begin
                                    AEffective[I,J] := 0;
                                    AParam[I,J] := 0;
                                    Inc(J);
                                end;
                                Inc(I);
                            end;
                            
                            //
                            // Prepare matrices
                            //
                            if IsUpper then
                            begin
                                I:=0;
                                while I<=N-1 do
                                begin
                                    J:=I;
                                    while J<=N-1 do
                                    begin
                                        AEffective[I,J] := 0.9*(2*RandomReal-1);
                                        AParam[I,J] := AEffective[I,J];
                                        Inc(J);
                                    end;
                                    AEffective[I,I] := (2*RandomInteger(2)-1)*(0.8+RandomReal);
                                    AParam[I,I] := AEffective[I,I];
                                    Inc(I);
                                end;
                            end
                            else
                            begin
                                I:=0;
                                while I<=N-1 do
                                begin
                                    J:=0;
                                    while J<=I do
                                    begin
                                        AEffective[I,J] := 0.9*(2*RandomReal-1);
                                        AParam[I,J] := AEffective[I,J];
                                        Inc(J);
                                    end;
                                    AEffective[I,I] := (2*RandomInteger(2)-1)*(0.8+RandomReal);
                                    AParam[I,I] := AEffective[I,I];
                                    Inc(I);
                                end;
                            end;
                            if Isunit then
                            begin
                                I:=0;
                                while I<=N-1 do
                                begin
                                    AEffective[I,I] := 1;
                                    AParam[I,I] := 0;
                                    Inc(I);
                                end;
                            end;
                            if IsTrans then
                            begin
                                if IsUpper then
                                begin
                                    I:=0;
                                    while I<=N-1 do
                                    begin
                                        J:=I+1;
                                        while J<=N-1 do
                                        begin
                                            AEffective[J,I] := AEffective[I,J];
                                            AEffective[I,J] := 0;
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
                                            AEffective[I,J] := AEffective[J,I];
                                            AEffective[J,I] := 0;
                                            Inc(J);
                                        end;
                                        Inc(I);
                                    end;
                                end;
                            end;
                            
                            //
                            // Prepare task, solve, compare
                            //
                            I:=0;
                            while I<=N-1 do
                            begin
                                XE[I] := 2*RandomReal-1;
                                Inc(I);
                            end;
                            I:=0;
                            while I<=N-1 do
                            begin
                                V := APVDotProduct(@AEffective[I][0], 0, N-1, @XE[0], 0, N-1);
                                B[I] := V;
                                Inc(I);
                            end;
                            RMatrixTRSafeSolve(AParam, N, B, S, IsUpper, IsTrans, Isunit);
                            APVMul(@XE[0], 0, N-1, S);
                            APVSub(@XE[0], 0, N-1, @B[0], 0, N-1);
                            V := APVDotProduct(@XE[0], 0, N-1, @XE[0], 0, N-1);
                            V := Sqrt(V);
                            WasErrors := WasErrors or AP_FP_Greater(V,Threshold);
                            Inc(CntM);
                        end;
                        Inc(CntT);
                    end;
                    Inc(CntU);
                end;
                Inc(CntS);
            end;
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    if  not Silent then
    begin
        Write(Format('TESTING RMatrixTRSafeSolve'#13#10'',[]));
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
procedure MakeACopy(const A : TReal2DArray;
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
Silent unit test
*************************************************************************)
function testsstunit_test_silent():Boolean;
begin
    Result := TestSST(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testsstunit_test():Boolean;
begin
    Result := TestSST(False);
end;


end.