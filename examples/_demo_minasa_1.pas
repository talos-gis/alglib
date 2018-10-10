
program _demo;
Array[0]
var
    N : AlglibInteger;
    I : AlglibInteger;
    State : MinASAState;
    Rep : MinASAReport;
    S : TReal1DArray;
    BndL : TReal1DArray;
    BndU : TReal1DArray;
    X : Double;
    Y : Double;
    Z : Double;
begin
    
    //
    // Function being minimized:
    //     F = x+2y+3z subject to 0<=x<=1, 0<=y<=1, 0<=z<=1.
    //
    N := 3;
    SetLength(S, N);
    SetLength(BndL, N);
    SetLength(BndU, N);
    I:=0;
    while I<=N-1 do
    begin
        S[I] := 1;
        BndL[I] := 0;
        BndU[I] := 1;
        Inc(I);
    end;
    MinASACreate(N, S, BndL, BndU, State);
    MinASASetCond(State, 0.0, 0.0, 0.00001, 0);
    MinASASetXRep(State, True);
    Write(Format(''#13#10''#13#10'F = x+2y+3z subject to 0<=x<=1, 0<=y<=1, 0<=z<=1'#13#10'',[]));
    Write(Format('OPTIMIZATION STARTED'#13#10'',[]));
    while MinASAIteration(State) do
    begin
        if State.NeedFG then
        begin
            X := State.X[0];
            Y := State.X[1];
            Z := State.X[2];
            State.F := X+2*Y+3*Z;
            State.G[0] := 1;
            State.G[1] := 2;
            State.G[2] := 3;
        end;
        if State.XUpdated then
        begin
            Write(Format('    F(%4.2f,%4.2f,%4.2f)=%0.3f'#13#10'',[
                State.X[0],
                State.X[1],
                State.X[2],
                State.F]));
        end;
    end;
    Write(Format('OPTIMIZATION STOPPED'#13#10'',[]));
    MinASAResults(State, S, Rep);
    
    //
    // output results
    //
    Write(Format('X = %4.2f (should be 0.00)'#13#10'',[
        S[0]]));
    Write(Format('Y = %4.2f (should be 0.00)'#13#10'',[
        S[1]]));
    Write(Format('Z = %4.2f (should be 0.00)'#13#10''#13#10''#13#10'',[
        S[2]]));
end.