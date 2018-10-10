
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
    // y'' = -y, y(0) = 0, y'(0)=1.
    //
    // This ODE may be written as first-order system:
    // y' =  z
    // z' = -y
    //
    // Its solution is well known in academic circles :)
    //
    // Three intermediate values are calculated,
    // plus starting and final points.
    //
    SetLength(Y, 2);
    Y[0] := 0;
    Y[1] := 1;
    SetLength(X, 5);
    X[0] := Pi*0/4;
    X[1] := Pi*1/4;
    X[2] := Pi*2/4;
    X[3] := Pi*3/4;
    X[4] := Pi*4/4;
    Eps := 1.0E-8;
    H := 0.01;
    ODESolverRKCK(Y, 2, X, 5, Eps, H, State);
    while ODESolverIteration(State) do
    begin
        State.DY[0] := State.Y[1];
        State.DY[1] := -State.Y[0];
    end;
    ODESolverResults(State, M, X, YTbl, Rep);
    Write(Format('     X   Y(X)     Error'#13#10'',[]));
    I:=0;
    while I<=M-1 do
    begin
        Write(Format('%6.3f %6.3f  %8.2e'#13#10'',[
            X[I],
            YTbl[I,0],
            AbsReal(YTbl[I,0]-Sin(X[I]))]));
        Inc(I);
    end;
end.