
program _demo;
Array[0]
var
    State : LMState;
    Rep : LMReport;
    S : TReal1DArray;
    X : Double;
    Y : Double;
begin
    
    //
    // Example of solving simple task using FGH scheme.
    //
    // Function minimized:
    //     F = (x-2*y)^2 + (x-2)^2 + (y-1)^2
    // exact solution is (2,1).
    //
    SetLength(S, 2);
    S[0] := RandomReal-0.5;
    S[1] := RandomReal-0.5;
    MinLMFGH(2, S, 0.0, 0.001, 0, State);
    while MinLMIteration(State) do
    begin
        X := State.X[0];
        Y := State.X[1];
        if State.NeedF then
        begin
            State.F := AP_Sqr(X-2*Y)+AP_Sqr(X-2)+AP_Sqr(Y-1);
        end;
        if State.NeedFG then
        begin
            State.F := AP_Sqr(X-2*Y)+AP_Sqr(X-2)+AP_Sqr(Y-1);
            State.G[0] := 2*(X-2*Y)+2*(X-2)+0;
            State.G[1] := -4*(X-2*Y)+0+2*(Y-1);
        end;
        if State.NeedFGH then
        begin
            State.F := AP_Sqr(X-2*Y)+AP_Sqr(X-2)+AP_Sqr(Y-1);
            State.G[0] := 2*(X-2*Y)+2*(X-2)+0;
            State.G[1] := -4*(X-2*Y)+0+2*(Y-1);
            State.H[0,0] := 4;
            State.H[1,0] := -4;
            State.H[0,1] := -4;
            State.H[1,1] := 10;
        end;
    end;
    MinLMResults(State, S, Rep);
    
    //
    // output results
    //
    Write(Format('X = %4.2f (correct value - 2.00)'#13#10'',[
        S[0]]));
    Write(Format('Y = %4.2f (correct value - 1.00)'#13#10'',[
        S[1]]));
    Write(Format('TerminationType = %0d (should be 2 - stopping when step is small enough)'#13#10'',[
        Rep.TerminationType]));
    Write(Format('NFunc = %0d'#13#10'',[
        Rep.NFunc]));
    Write(Format('NJac  = %0d'#13#10'',[
        Rep.NJac]));
    Write(Format('NGrad = %0d'#13#10'',[
        Rep.NGrad]));
    Write(Format('NHess = %0d (should be 1 - task is very simple)'#13#10'',[
        Rep.NHess]));
end.