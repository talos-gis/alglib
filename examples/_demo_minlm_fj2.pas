
program _demo;
Array[0]
var
    State : MinLMState;
    Rep : MinLMReport;
    I : AlglibInteger;
    S : TReal1DArray;
    X : TReal1DArray;
    Y : TReal1DArray;
    FI : Double;
    N : AlglibInteger;
    M : AlglibInteger;
begin
    
    //
    // Example of solving polynomial approximation task using FJ scheme.
    //
    // Data points:
    //     xi are random numbers from [-1,+1],
    //
    // Function being fitted:
    //     yi = exp(xi) - sin(xi) - x^3/3
    //
    // Function being minimized:
    //     F(a,b,c) =
    //         (a + b*x0 + c*x0^2 - y0)^2 +
    //         (a + b*x1 + c*x1^2 - y1)^2 + ...
    //
    N := 3;
    SetLength(S, N);
    I:=0;
    while I<=N-1 do
    begin
        S[I] := RandomReal-0.5;
        Inc(I);
    end;
    M := 100;
    SetLength(X, M);
    SetLength(Y, M);
    I:=0;
    while I<=M-1 do
    begin
        X[I] := AP_Double(2*I)/(M-1)-1;
        Y[I] := Exp(X[I])-Sin(X[I])-X[I]*X[I]*X[I]/3;
        Inc(I);
    end;
    
    //
    // Now S stores starting point, X and Y store points being fitted.
    //
    MinLMCreateFJ(N, M, S, State);
    MinLMSetCond(State, 0.0, 0.0, 0.001, 0);
    while MinLMIteration(State) do
    begin
        if State.NeedF then
        begin
            State.F := 0;
        end;
        I:=0;
        while I<=M-1 do
        begin
            
            //
            // "a" is stored in State.X[0]
            // "b" - State.X[1]
            // "c" - State.X[2]
            //
            FI := State.X[0]+State.X[1]*X[I]+State.X[2]*AP_Sqr(X[I])-Y[I];
            if State.NeedF then
            begin
                
                //
                // F is equal to sum of fi squared.
                //
                State.F := State.F+AP_Sqr(FI);
            end;
            if State.NeedFiJ then
            begin
                
                //
                // Fi
                //
                State.Fi[I] := FI;
                
                //
                // dFi/da
                //
                State.J[I,0] := 1;
                
                //
                // dFi/db
                //
                State.J[I,1] := X[I];
                
                //
                // dFi/dc
                //
                State.J[I,2] := AP_Sqr(X[I]);
            end;
            Inc(I);
        end;
    end;
    MinLMResults(State, S, Rep);
    
    //
    // output results
    //
    Write(Format('A = %4.2f'#13#10'',[
        S[0]]));
    Write(Format('B = %4.2f'#13#10'',[
        S[1]]));
    Write(Format('C = %4.2f'#13#10'',[
        S[2]]));
    Write(Format('TerminationType = %0d (should be 2 - stopping when step is small enough)'#13#10'',[
        Rep.TerminationType]));
end.