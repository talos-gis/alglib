unit testautogk;
interface
uses Math, Sysutils, Ap, tsort, blas, rotations, tdevd, gammafunc, gq, gkq, autogk;

function TestAutoGKunit(Silent : Boolean):Boolean;
function testautogk_test_silent():Boolean;
function testautogk_test():Boolean;

implementation

(*************************************************************************
Test
*************************************************************************)
function TestAutoGKunit(Silent : Boolean):Boolean;
var
    A : Double;
    B : Double;
    State : AutoGKState;
    Rep : AutoGKReport;
    V : Double;
    Exact : Double;
    EAbs : Double;
    Alpha : Double;
    PKind : AlglibInteger;
    ErrTol : Double;
    SimpleErrors : Boolean;
    SngEndErrors : Boolean;
    WasErrors : Boolean;
begin
    SimpleErrors := False;
    SngEndErrors := False;
    WasErrors := False;
    ErrTol := 10000*MachineEpsilon;
    
    //
    // Simple test: integral(exp(x),+-1,+-2), no maximum width requirements
    //
    A := (2*RandomInteger(2)-1)*1.0;
    B := (2*RandomInteger(2)-1)*2.0;
    AutoGKSmooth(A, B, State);
    while AutoGKIteration(State) do
    begin
        State.F := Exp(State.X);
    end;
    AutoGKResults(State, V, Rep);
    Exact := Exp(B)-Exp(A);
    EAbs := AbsReal(Exp(B)-Exp(A));
    if Rep.TerminationType<=0 then
    begin
        SimpleErrors := True;
    end
    else
    begin
        SimpleErrors := SimpleErrors or AP_FP_Greater(AbsReal(Exact-V),ErrTol*EAbs);
    end;
    
    //
    // Simple test: integral(exp(x),+-1,+-2), XWidth=0.1
    //
    A := (2*RandomInteger(2)-1)*1.0;
    B := (2*RandomInteger(2)-1)*2.0;
    AutoGKSmoothW(A, B, 0.1, State);
    while AutoGKIteration(State) do
    begin
        State.F := Exp(State.X);
    end;
    AutoGKResults(State, V, Rep);
    Exact := Exp(B)-Exp(A);
    EAbs := AbsReal(Exp(B)-Exp(A));
    if Rep.TerminationType<=0 then
    begin
        SimpleErrors := True;
    end
    else
    begin
        SimpleErrors := SimpleErrors or AP_FP_Greater(AbsReal(Exact-V),ErrTol*EAbs);
    end;
    
    //
    // Simple test: integral(cos(100*x),0,2*pi), no maximum width requirements
    //
    A := 0;
    B := 2*Pi;
    AutoGKSmooth(A, B, State);
    while AutoGKIteration(State) do
    begin
        State.F := Cos(100*State.X);
    end;
    AutoGKResults(State, V, Rep);
    Exact := 0;
    EAbs := 4;
    if Rep.TerminationType<=0 then
    begin
        SimpleErrors := True;
    end
    else
    begin
        SimpleErrors := SimpleErrors or AP_FP_Greater(AbsReal(Exact-V),ErrTol*EAbs);
    end;
    
    //
    // Simple test: integral(cos(100*x),0,2*pi), XWidth=0.3
    //
    A := 0;
    B := 2*Pi;
    AutoGKSmoothW(A, B, 0.3, State);
    while AutoGKIteration(State) do
    begin
        State.F := Cos(100*State.X);
    end;
    AutoGKResults(State, V, Rep);
    Exact := 0;
    EAbs := 4;
    if Rep.TerminationType<=0 then
    begin
        SimpleErrors := True;
    end
    else
    begin
        SimpleErrors := SimpleErrors or AP_FP_Greater(AbsReal(Exact-V),ErrTol*EAbs);
    end;
    
    //
    // singular problem on [a,b] = [0.1, 0.5]
    //     f2(x) = (1+x)*(b-x)^alpha, -1 < alpha < 1
    //
    PKind:=0;
    while PKind<=6 do
    begin
        A := 0.1;
        B := 0.5;
        if PKind=0 then
        begin
            Alpha := -0.9;
        end;
        if PKind=1 then
        begin
            Alpha := -0.5;
        end;
        if PKind=2 then
        begin
            Alpha := -0.1;
        end;
        if PKind=3 then
        begin
            Alpha := 0.0;
        end;
        if PKind=4 then
        begin
            Alpha := 0.1;
        end;
        if PKind=5 then
        begin
            Alpha := 0.5;
        end;
        if PKind=6 then
        begin
            Alpha := 0.9;
        end;
        
        //
        // f1(x) = (1+x)*(x-a)^alpha, -1 < alpha < 1
        // 1. use singular integrator for [a,b]
        // 2. use singular integrator for [b,a]
        //
        Exact := Power(B-A, Alpha+2)/(Alpha+2)+(1+A)*Power(B-A, Alpha+1)/(Alpha+1);
        EAbs := AbsReal(Exact);
        AutoGKSingular(A, B, Alpha, 0.0, State);
        while AutoGKIteration(State) do
        begin
            if AP_FP_Less(State.XMinusA,0.01) then
            begin
                State.F := Power(State.XMinusA, Alpha)*(1+State.X);
            end
            else
            begin
                State.F := Power(State.X-A, Alpha)*(1+State.X);
            end;
        end;
        AutoGKResults(State, V, Rep);
        if Rep.TerminationType<=0 then
        begin
            SngEndErrors := True;
        end
        else
        begin
            SngEndErrors := SngEndErrors or AP_FP_Greater(AbsReal(V-Exact),ErrTol*EAbs);
        end;
        AutoGKSingular(B, A, 0.0, Alpha, State);
        while AutoGKIteration(State) do
        begin
            if AP_FP_Greater(State.BMinusX,-0.01) then
            begin
                State.F := Power(-State.BMinusX, Alpha)*(1+State.X);
            end
            else
            begin
                State.F := Power(State.X-A, Alpha)*(1+State.X);
            end;
        end;
        AutoGKResults(State, V, Rep);
        if Rep.TerminationType<=0 then
        begin
            SngEndErrors := True;
        end
        else
        begin
            SngEndErrors := SngEndErrors or AP_FP_Greater(AbsReal(-V-Exact),ErrTol*EAbs);
        end;
        
        //
        // f1(x) = (1+x)*(b-x)^alpha, -1 < alpha < 1
        // 1. use singular integrator for [a,b]
        // 2. use singular integrator for [b,a]
        //
        Exact := (1+B)*Power(B-A, Alpha+1)/(Alpha+1)-Power(B-A, Alpha+2)/(Alpha+2);
        EAbs := AbsReal(Exact);
        AutoGKSingular(A, B, 0.0, Alpha, State);
        while AutoGKIteration(State) do
        begin
            if AP_FP_Less(State.BMinusX,0.01) then
            begin
                State.F := Power(State.BMinusX, Alpha)*(1+State.X);
            end
            else
            begin
                State.F := Power(B-State.X, Alpha)*(1+State.X);
            end;
        end;
        AutoGKResults(State, V, Rep);
        if Rep.TerminationType<=0 then
        begin
            SngEndErrors := True;
        end
        else
        begin
            SngEndErrors := SngEndErrors or AP_FP_Greater(AbsReal(V-Exact),ErrTol*EAbs);
        end;
        AutoGKSingular(B, A, Alpha, 0.0, State);
        while AutoGKIteration(State) do
        begin
            if AP_FP_Greater(State.XMinusA,-0.01) then
            begin
                State.F := Power(-State.XMinusA, Alpha)*(1+State.X);
            end
            else
            begin
                State.F := Power(B-State.X, Alpha)*(1+State.X);
            end;
        end;
        AutoGKResults(State, V, Rep);
        if Rep.TerminationType<=0 then
        begin
            SngEndErrors := True;
        end
        else
        begin
            SngEndErrors := SngEndErrors or AP_FP_Greater(AbsReal(-V-Exact),ErrTol*EAbs);
        end;
        Inc(PKind);
    end;
    
    //
    // end
    //
    WasErrors := SimpleErrors or SngEndErrors;
    if  not Silent then
    begin
        Write(Format('TESTING AUTOGK'#13#10'',[]));
        Write(Format('INTEGRATION WITH GIVEN ACCURACY:          ',[]));
        if SimpleErrors or SngEndErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* SIMPLE PROBLEMS:                        ',[]));
        if SimpleErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* SINGULAR PROBLEMS (ENDS OF INTERVAL):   ',[]));
        if SngEndErrors then
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
    Result :=  not WasErrors;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testautogk_test_silent():Boolean;
begin
    Result := TestAutoGKunit(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testautogk_test():Boolean;
begin
    Result := TestAutoGKunit(False);
end;


end.