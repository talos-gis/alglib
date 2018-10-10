
program _demo;
Array[0]
var
    X : TReal1DArray;
    Y : TReal1DArray;
    YTbl : TReal2DArray;
    Eps : Double;
    H : Double;
    M : AlglibInteger;
    I : AlglibInteger;
    State : ODESolverState;
    Rep : ODESolverReport;
begin
    
    //
    // ODESolver unit is used to solve simple ODE:
    // y' = y, y(0) = 1.
    //
    // Its solution is well known in academic circles :)
    //
    // No intermediate values are calculated,
    // just starting and final points.
    //
    SetLength(Y, 1);
    Y[0] := 1;
    SetLength(X, 2);
    X[0] := 0;
    X[1] := 1;
    Eps := 1.0E-4;
    H := 0.01;
    ODESolverRKCK(Y, 1, X, 2, Eps, H, State);
    while ODESolverIteration(State) do
    begin
        State.DY[0] := State.Y[0];
    end;
    ODESolverResults(State, M, X, YTbl, Rep);
    Write(Format('    X  Y(X)'#13#10'',[]));
    I:=0;
    while I<=M-1 do
    begin
        Write(Format('%5.3f %5.3f'#13#10'',[
            X[I],
            YTbl[I,0]]));
        Inc(I);
    end;
end.