
program _demo;
Array[0]
var
    M : AlglibInteger;
    N : AlglibInteger;
    X : TReal1DArray;
    Y : TReal1DArray;
    W : TReal1DArray;
    XC : TReal1DArray;
    YC : TReal1DArray;
    DC : TInteger1DArray;
    Rep : PolynomialFitReport;
    Info : AlglibInteger;
    P : BarycentricInterpolant;
    I : AlglibInteger;
    J : AlglibInteger;
    A : Double;
    B : Double;
    V : Double;
    DV : Double;
begin
    Write(Format(''#13#10''#13#10'Fitting exp(2*x) at [-1,+1] by polinomial'#13#10''#13#10'',[]));
    Write(Format('Fit type             rms.err max.err    p(0)   dp(0)'#13#10'',[]));
    
    //
    // Prepare points
    //
    M := 5;
    A := -1;
    B := +1;
    N := 1000;
    SetLength(X, N);
    SetLength(Y, N);
    SetLength(W, N);
    I:=0;
    while I<=N-1 do
    begin
        X[I] := A+(B-A)*I/(N-1);
        Y[I] := Exp(2*X[I]);
        W[I] := 1.0;
        Inc(I);
    end;
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5th degree polynomial
    // c) without constraints
    //
    PolynomialFit(X, Y, N, M, Info, P, Rep);
    BarycentricDiff1(P, 0.0, V, DV);
    Write(Format('Unconstrained        %7.4f %7.4f %7.4f %7.4f'#13#10'',[
        Rep.RMSError,
        Rep.MaxError,
        V,
        DV]));
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5th degree polynomial
    // c) constrained: p(0)=1
    //
    SetLength(XC, 1);
    SetLength(YC, 1);
    SetLength(DC, 1);
    XC[0] := 0;
    YC[0] := 1;
    DC[0] := 0;
    PolynomialFitWC(X, Y, W, N, XC, YC, DC, 1, M, Info, P, Rep);
    BarycentricDiff1(P, 0.0, V, DV);
    Write(Format('Constrained, p(0)=1  %7.4f %7.4f %7.4f %7.4f'#13#10'',[
        Rep.RMSError,
        Rep.MaxError,
        V,
        DV]));
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5th degree polynomial
    // c) constrained: dp(0)=2
    //
    SetLength(XC, 1);
    SetLength(YC, 1);
    SetLength(DC, 1);
    XC[0] := 0;
    YC[0] := 2;
    DC[0] := 1;
    PolynomialFitWC(X, Y, W, N, XC, YC, DC, 1, M, Info, P, Rep);
    BarycentricDiff1(P, 0.0, V, DV);
    Write(Format('Constrained, dp(0)=2 %7.4f %7.4f %7.4f %7.4f'#13#10'',[
        Rep.RMSError,
        Rep.MaxError,
        V,
        DV]));
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5th degree polynomial
    // c) constrained: p(0)=1, dp(0)=2
    //
    SetLength(XC, 2);
    SetLength(YC, 2);
    SetLength(DC, 2);
    XC[0] := 0;
    YC[0] := 1;
    DC[0] := 0;
    XC[1] := 0;
    YC[1] := 2;
    DC[1] := 1;
    PolynomialFitWC(X, Y, W, N, XC, YC, DC, 2, M, Info, P, Rep);
    BarycentricDiff1(P, 0.0, V, DV);
    Write(Format('Constrained, both    %7.4f %7.4f %7.4f %7.4f'#13#10'',[
        Rep.RMSError,
        Rep.MaxError,
        V,
        DV]));
    Write(Format(''#13#10''#13#10'',[]));
end.