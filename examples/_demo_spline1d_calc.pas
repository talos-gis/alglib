
program _demo;
Array[0]
var
    X : TReal1DArray;
    Y : TReal1DArray;
    N : AlglibInteger;
    I : AlglibInteger;
    T : Double;
    S : Spline1DInterpolant;
    V : Double;
    DV : Double;
    D2V : Double;
    Err : Double;
    MaxErr : Double;
begin
    
    //
    // Demonstration of Spline1DCalc(), Spline1DDiff(), Spline1DIntegrate()
    //
    Write(Format('DEMONSTRATION OF Spline1DCalc(), Spline1DDiff(), Spline1DIntegrate()'#13#10''#13#10'',[]));
    Write(Format('F(x)=sin(x), [0, pi]'#13#10'',[]));
    Write(Format('Natural cubic spline with 3 nodes is used'#13#10''#13#10'',[]));
    
    //
    // Create spline
    //
    N := 3;
    SetLength(X, N);
    SetLength(Y, N);
    I:=0;
    while I<=N-1 do
    begin
        X[I] := Pi*I/(N-1);
        Y[I] := Sin(X[I]);
        Inc(I);
    end;
    Spline1DBuildCubic(X, Y, N, 2, 0.0, 2, 0.0, S);
    
    //
    // Output results
    //
    Spline1DDiff(S, 0, V, DV, D2V);
    Write(Format('                 S(x)    F(x) '#13#10'',[]));
    Write(Format('function       %6.3f  %6.3f '#13#10'',[
        Spline1DCalc(S, 0),
        Double(Variant(0))]));
    Write(Format('d/dx(0)        %6.3f  %6.3f '#13#10'',[
        DV,
        Double(Variant(1))]));
    Write(Format('d2/dx2(0)      %6.3f  %6.3f '#13#10'',[
        D2V,
        Double(Variant(0))]));
    Write(Format('integral(0,pi) %6.3f  %6.3f '#13#10'',[
        Spline1DIntegrate(S, Pi),
        Double(Variant(2))]));
    Write(Format(''#13#10''#13#10'',[]));
end.