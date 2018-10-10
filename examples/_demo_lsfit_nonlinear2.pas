
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
    Write(Format('Fitting 1-x^2 on [-1,+1] with cos(alpha*pi*x)+beta'#13#10'',[]));
    
    //
    // Fitting 1-x^2 on [-1,+1] with cos(alpha*pi*x)+beta:
    // * using Hessian
    // * using alpha=1 and beta=0 as initial values
    // * using 1000 uniformly distributed points to fit to
    //
    // Notes:
    // * N - number of points
    // * M - dimension of space where points reside
    // * K - number of parameters being fitted
    //
    N := 1000;
    M := 1;
    K := 2;
    A := -1;
    B := +1;
    
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
        Y[I] := 1-AP_Sqr(X[I,0]);
        Inc(I);
    end;
    C[0] := 1.0;
    C[1] := 0.0;
    EpsF := 0.0;
    EpsX := 0.0001;
    MaxIts := 0;
    
    //
    // Solve
    //
    LSFitNonlinearFGH(X, Y, C, N, M, K, EpsF, EpsX, MaxIts, State);
    while LSFitNonlinearIteration(State) do
    begin
        
        //
        // F(x) = Cos(alpha*pi*x)+beta
        //
        State.F := Cos(State.C[0]*Pi*State.X[0])+State.C[1];
        
        //
        // F(x)      = Cos(alpha*pi*x)+beta
        // dF/dAlpha = -pi*x*Sin(alpha*pi*x)
        // dF/dBeta  = 1.0
        //
        if State.NeedFG or State.NeedFGH then
        begin
            State.G[0] := -Pi*State.X[0]*Sin(State.C[0]*Pi*State.X[0]);
            State.G[1] := 1.0;
        end;
        
        //
        // F(x)            = Cos(alpha*pi*x)+beta
        // d2F/dAlpha2     = -(pi*x)^2*Cos(alpha*pi*x)
        // d2F/dAlphadBeta = 0
        // d2F/dBeta2     =  0
        //
        if State.NeedFGH then
        begin
            State.H[0,0] := -AP_Sqr(Pi*State.X[0])*Cos(State.C[0]*Pi*State.X[0]);
            State.H[0,1] := 0.0;
            State.H[1,0] := 0.0;
            State.H[1,1] := 0.0;
        end;
    end;
    LSFitNonlinearResults(State, Info, C, Rep);
    Write(Format('alpha:   %0.3f'#13#10'',[
        C[0]]));
    Write(Format('beta:    %0.3f'#13#10'',[
        C[1]]));
    Write(Format('rms.err: %0.3f'#13#10'',[
        Rep.RMSError]));
    Write(Format('max.err: %0.3f'#13#10'',[
        Rep.MaxError]));
    Write(Format('Termination type: %0d'#13#10'',[
        Info]));
    Write(Format(''#13#10''#13#10'',[]));
end.