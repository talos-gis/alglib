
program _demo;
Array[0]
var
    N : AlglibInteger;
    State : MinCGState;
    Rep : MinCGReport;
    S : TReal1DArray;
    X : Double;
    Y : Double;
begin
    
    //
    // Function minimized:
    //     F = exp(x-1) + exp(1-x) + (y-x)^2
    // N = 2 - task dimension.
    //
    // Take a look at MinCGSetStpMax() call - it prevents us
    // from overflow (which may be result of too large step).
    // Try to comment it and see what will happen.
    //
    N := 2;
    SetLength(S, 2);
    S[0] := 10;
    S[1] := RandomReal-0.5;
    MinCGCreate(N, S, State);
    MinCGSetCond(State, 0.0, 0.0, 0.0001, 0);
    MinCGSetXRep(State, True);
    MinCGSetStpMax(State, 1.0);
    Write(Format(''#13#10''#13#10'F = exp(x-1) + exp(1-x) + (y-x)^2'#13#10'',[]));
    Write(Format('OPTIMIZATION STARTED'#13#10'',[]));
    while MinCGIteration(State) do
    begin
        if State.NeedFG then
        begin
            X := State.X[0];
            Y := State.X[1];
            State.F := Exp(X-1)+Exp(1-X)+AP_Sqr(Y-X);
            State.G[0] := Exp(X-1)-Exp(1-X)+2*(X-Y);
            State.G[1] := 2*(Y-X);
        end;
        if State.XUpdated then
        begin
            Write(Format('    F(%8.5f,%8.5f)=%0.5f'#13#10'',[
                State.X[0],
                State.X[1],
                State.F]));
        end;
    end;
    Write(Format('OPTIMIZATION STOPPED'#13#10'',[]));
    MinCGResults(State, S, Rep);
    
    //
    // output results
    //
    Write(Format('X = %4.2f (should be 1.00)'#13#10'',[
        S[0]]));
    Write(Format('Y = %4.2f (should be 1.00)'#13#10''#13#10''#13#10'',[
        S[1]]));
end.