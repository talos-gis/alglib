
program _demo;
Array[0]
var
    N : AlglibInteger;
    M : AlglibInteger;
    State : LBFGSState;
    Rep : LBFGSReport;
    S : TReal1DArray;
    X : Double;
    Y : Double;
begin
    
    //
    // Function minimized:
    //     F = exp(x-1) + exp(1-x) + (y-x)^2
    // N = 2 - task dimension
    // M = 1 - build tank-1 model
    //
    N := 2;
    M := 1;
    SetLength(S, 2);
    S[0] := RandomReal-0.5;
    S[1] := RandomReal-0.5;
    MinLBFGS(N, M, S, 0.0, 0.0, 0.0001, 0, 0, State);
    while MinLBFGSIteration(State) do
    begin
        X := State.X[0];
        Y := State.X[1];
        State.F := Exp(X-1)+Exp(1-X)+AP_Sqr(Y-X);
        State.G[0] := Exp(X-1)-Exp(1-X)+2*(X-Y);
        State.G[1] := 2*(Y-X);
    end;
    MinLBFGSResults(State, S, Rep);
    
    //
    // output results
    //
    Write(Format('F = exp(x-1) + exp(1-x) + (y-x)^2'#13#10'',[]));
    Write(Format('X = %4.2f (should be 1.00)'#13#10'',[
        S[0]]));
    Write(Format('Y = %4.2f (should be 1.00)'#13#10'',[
        S[1]]));
end.