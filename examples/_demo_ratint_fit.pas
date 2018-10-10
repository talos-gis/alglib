
program _demo;
Array[0]
var
    M : AlglibInteger;
    N : AlglibInteger;
    D : AlglibInteger;
    X : TReal1DArray;
    Y : TReal1DArray;
    W : TReal1DArray;
    XC : TReal1DArray;
    YC : TReal1DArray;
    DC : TInteger1DArray;
    Rep : BarycentricFitReport;
    Info : AlglibInteger;
    R : BarycentricInterpolant;
    I : AlglibInteger;
    J : AlglibInteger;
    A : Double;
    B : Double;
    V : Double;
    DV : Double;
begin
    Write(Format(''#13#10''#13#10'Fitting exp(2*x) at [-1,+1] by:'#13#10'1. constrained/unconstrained Floater-Hormann functions'#13#10'',[]));
    Write(Format(''#13#10'',[]));
    Write(Format('Fit type                rms.err max.err    p(0)   dp(0)  DBest'#13#10'',[]));
    
    //
    // Prepare points
    //
    M := 5;
    A := -1;
    B := +1;
    N := 10000;
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
    // b) by 5 Floater-Hormann functions
    // c) without constraints
    //
    BarycentricFitFloaterHormann(X, Y, N, M, Info, R, Rep);
    BarycentricDiff1(R, 0.0, V, DV);
    Write(Format('Unconstrained FH        %7.4f %7.4f %7.4f %7.4f      %0d'#13#10'',[
        Rep.RMSError,
        Rep.MaxError,
        V,
        DV,
        Rep.DBest]));
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5 Floater-Hormann functions
    // c) constrained: p(0)=1
    //
    SetLength(XC, 1);
    SetLength(YC, 1);
    SetLength(DC, 1);
    XC[0] := 0;
    YC[0] := 1;
    DC[0] := 0;
    BarycentricFitFloaterHormannWC(X, Y, W, N, XC, YC, DC, 1, M, Info, R, Rep);
    BarycentricDiff1(R, 0.0, V, DV);
    Write(Format('Constrained FH, p(0)=1  %7.4f %7.4f %7.4f %7.4f      %0d'#13#10'',[
        Rep.RMSError,
        Rep.MaxError,
        V,
        DV,
        Rep.DBest]));
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5 Floater-Hormann functions
    // c) constrained: dp(0)=2
    //
    SetLength(XC, 1);
    SetLength(YC, 1);
    SetLength(DC, 1);
    XC[0] := 0;
    YC[0] := 2;
    DC[0] := 1;
    BarycentricFitFloaterHormannWC(X, Y, W, N, XC, YC, DC, 1, M, Info, R, Rep);
    BarycentricDiff1(R, 0.0, V, DV);
    Write(Format('Constrained FH, dp(0)=2 %7.4f %7.4f %7.4f %7.4f      %0d'#13#10'',[
        Rep.RMSError,
        Rep.MaxError,
        V,
        DV,
        Rep.DBest]));
    
    //
    // Fitting:
    // a) f(x)=exp(2*x) at [-1,+1]
    // b) by 5 Floater-Hormann functions
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
    BarycentricFitFloaterHormannWC(X, Y, W, N, XC, YC, DC, 2, M, Info, R, Rep);
    BarycentricDiff1(R, 0.0, V, DV);
    Write(Format('Constrained FH, both    %7.4f %7.4f %7.4f %7.4f      %0d'#13#10'',[
        Rep.RMSError,
        Rep.MaxError,
        V,
        DV,
        Rep.DBest]));
    Write(Format(''#13#10''#13#10'',[]));
end.