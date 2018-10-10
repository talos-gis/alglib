
program _demo;
Array[0]
var
    M : AlglibInteger;
    N : AlglibInteger;
    Y : TReal1DArray;
    FMatrix : TReal2DArray;
    CMatrix : TReal2DArray;
    Rep : LSFitReport;
    Info : AlglibInteger;
    C : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    X : Double;
    A : Double;
    B : Double;
begin
    Write(Format(''#13#10''#13#10'Fitting tan(x) by third degree polynomial'#13#10''#13#10'',[]));
    Write(Format('Fit type             rms.err max.err    p(0)   dp(0)'#13#10'',[]));
    
    //
    // Fitting tan(x) at [0, 0.4*pi] by third degree polynomial:
    // a) without constraints
    // b) constrained at x=0: p(0)=0
    // c) constrained at x=0: p'(0)=1
    // c) constrained at x=0: p(0)=0, p'(0)=1
    //
    M := 4;
    N := 100;
    A := 0;
    B := 0.4*Pi;
    
    //
    // Prepare task matrix
    //
    SetLength(Y, N);
    SetLength(FMatrix, N, M);
    I:=0;
    while I<=N-1 do
    begin
        X := A+(B-A)*I/(N-1);
        Y[I] := Tan(X);
        FMatrix[I,0] := 1.0;
        J:=1;
        while J<=M-1 do
        begin
            FMatrix[I,J] := X*FMatrix[I,J-1];
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Solve unconstrained task
    //
    LSFitLinear(Y, FMatrix, N, M, Info, C, Rep);
    Write(Format('Unconstrained        %7.4f %7.4f %7.4f %7.4f'#13#10'',[
        Rep.RMSError,
        Rep.MaxError,
        C[0],
        C[1]]));
    
    //
    // Solve constrained task, p(0)=0
    // Prepare constraints matrix:
    // * first M columns store values of basis functions at X=0
    // * last column stores zero (desired value at X=0)
    //
    SetLength(CMatrix, 1, M+1);
    CMatrix[0,0] := 1;
    I:=1;
    while I<=M-1 do
    begin
        CMatrix[0,I] := 0;
        Inc(I);
    end;
    CMatrix[0,M] := 0;
    LSFitLinearC(Y, FMatrix, CMatrix, N, M, 1, Info, C, Rep);
    Write(Format('Constrained, p(0)=0  %7.4f %7.4f %7.4f %7.4f'#13#10'',[
        Rep.RMSError,
        Rep.MaxError,
        C[0],
        C[1]]));
    
    //
    // Solve constrained task, p'(0)=0
    // Prepare constraints matrix:
    // * first M columns store derivatives of basis functions at X=0
    // * last column stores 1.0 (desired derivative at X=0)
    //
    SetLength(CMatrix, 1, M+1);
    I:=0;
    while I<=M-1 do
    begin
        CMatrix[0,I] := 0;
        Inc(I);
    end;
    CMatrix[0,1] := 1;
    CMatrix[0,M] := 1;
    LSFitLinearC(Y, FMatrix, CMatrix, N, M, 1, Info, C, Rep);
    Write(Format('Constrained, dp(0)=1 %7.4f %7.4f %7.4f %7.4f'#13#10'',[
        Rep.RMSError,
        Rep.MaxError,
        C[0],
        C[1]]));
    
    //
    // Solve constrained task, p(0)=0, p'(0)=0
    // Prepare constraints matrix:
    // * first M columns store values/derivatives of basis functions at X=0
    // * last column stores desired values/derivative at X=0
    //
    SetLength(CMatrix, 2, M+1);
    CMatrix[0,0] := 1;
    I:=1;
    while I<=M-1 do
    begin
        CMatrix[0,I] := 0;
        Inc(I);
    end;
    CMatrix[0,M] := 0;
    I:=0;
    while I<=M-1 do
    begin
        CMatrix[1,I] := 0;
        Inc(I);
    end;
    CMatrix[1,1] := 1;
    CMatrix[1,M] := 1;
    LSFitLinearC(Y, FMatrix, CMatrix, N, M, 2, Info, C, Rep);
    Write(Format('Constrained, both    %7.4f %7.4f %7.4f %7.4f'#13#10'',[
        Rep.RMSError,
        Rep.MaxError,
        C[0],
        C[1]]));
    Write(Format(''#13#10''#13#10'',[]));
end.