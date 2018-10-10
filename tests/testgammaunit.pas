unit testgammaunit;
interface
uses Math, Sysutils, Ap, gammafunc;

function TestGamma(Silent : Boolean):Boolean;
function testgammaunit_test_silent():Boolean;
function testgammaunit_test():Boolean;

implementation

function TestGamma(Silent : Boolean):Boolean;
var
    Threshold : Double;
    V : Double;
    S : Double;
    WasErrors : Boolean;
    GammaErrors : Boolean;
    LnGammaErrors : Boolean;
begin
    GammaErrors := False;
    LnGammaErrors := False;
    WasErrors := False;
    Threshold := 100*MachineEpsilon;
    
    //
    //
    //
    GammaErrors := GammaErrors or AP_FP_Greater(AbsReal(Gamma(0.5)-Sqrt(Pi)),Threshold);
    GammaErrors := GammaErrors or AP_FP_Greater(AbsReal(Gamma(1.5)-0.5*Sqrt(Pi)),Threshold);
    V := LnGamma(0.5, S);
    LnGammaErrors := LnGammaErrors or AP_FP_Greater(AbsReal(V-Ln(Sqrt(Pi))),Threshold) or AP_FP_Neq(S,1);
    V := LnGamma(1.5, S);
    LnGammaErrors := LnGammaErrors or AP_FP_Greater(AbsReal(V-Ln(0.5*Sqrt(Pi))),Threshold) or AP_FP_Neq(S,1);
    
    //
    // report
    //
    WasErrors := GammaErrors or LnGammaErrors;
    if  not Silent then
    begin
        Write(Format('TESTING GAMMA FUNCTION'#13#10'',[]));
        Write(Format('GAMMA:                                   ',[]));
        if GammaErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('LN GAMMA:                                ',[]));
        if LnGammaErrors then
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
        Write(Format(''#13#10''#13#10'',[]));
    end;
    
    //
    // end
    //
    Result :=  not WasErrors;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testgammaunit_test_silent():Boolean;
begin
    Result := TestGamma(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testgammaunit_test():Boolean;
begin
    Result := TestGamma(False);
end;


end.