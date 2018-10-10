
program _demo;
Array[0]
var
    M : AlglibInteger;
    N : AlglibInteger;
    K : AlglibInteger;
    Y : TReal1DArray;
    X : TReal2DArray;
    C : TReal1DArray;
    Rep : LSFitReport;
    State : LSFitState;
    Info : AlglibInteger;
    EpsF : Double;
    EpsX : Double;
    MaxIts : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    A : Double;
    B : Double;
begin
    Write(Format('Fitting 0.5(1+cos(x)) on [-pi,+pi] with exp(-alpha*x^2)'#13#10'',[]));
    
    //
    // Fitting 0.5(1+cos(x)) on [-pi,+pi] with Gaussian exp(-alpha*x^2):
    // * without Hessian (gradient only)
    // * using alpha=1 as initial value
    // * using 1000 uniformly distributed points to fit to
    //
    // Notes:
    // * N - number of points
    // * M - dimension of space where points reside
    // * K - number of parameters being fitted
    //
    N := 1000;
    M := 1;
    K := 1;
    A := -Pi;
    B := +Pi;
    
    //
    // Prepare task matrix
    //
    SetLength(Y, N);
    SetLength(X, N, M);
    SetLength(C, K);
    I:=0;
    while I<=N-1 do
    begin
        X[I,0] := A+(B-A)*I/(N-1);
        Y[I] := 0.5*(1+Cos(X[I,0]));
        Inc(I);
    end;
    C[0] := 1.0;
    EpsF := 0.0;
    EpsX := 0.0001;
    MaxIts := 0;
    
    //
    // Solve
    //
    LSFitNonlinearFG(X, Y, C, N, M, K, EpsF, EpsX, MaxIts, True, State);
    while LSFitNonlinearIteration(State) do
    begin
        if State.NeedF then
        begin
            
            //
            // F(x) = Exp(-alpha*x^2)
            //
            State.F := Exp(-State.C[0]*AP_Sqr(State.X[0]));
        end;
        if State.NeedFG then
        begin
            
            //
            // F(x)      = Exp(-alpha*x^2)
            // dF/dAlpha = (-x^2)*Exp(-alpha*x^2)
            //
            State.F := Exp(-State.C[0]*AP_Sqr(State.X[0]));
            State.G[0] := -AP_Sqr(State.X[0])*State.F;
        end;
    end;
    LSFitNonlinearResults(State, Info, C, Rep);
    Write(Format('alpha:   %0.3f'#13#10'',[
        C[0]]));
    Write(Format('rms.err: %0.3f'#13#10'',[
        Rep.RMSError]));
    Write(Format('max.err: %0.3f'#13#10'',[
        Rep.MaxError]));
    Write(Format('Termination type: %0d'#13#10'',[
        Info]));
    Write(Format(''#13#10''#13#10'',[]));
end.