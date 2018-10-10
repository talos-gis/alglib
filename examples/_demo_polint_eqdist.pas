
program _demo;
Array[0]
var
    Y : TReal1DArray;
    N : AlglibInteger;
    I : AlglibInteger;
    T : Double;
    P : BarycentricInterpolant;
    V : Double;
    DV : Double;
    D2V : Double;
    Err : Double;
    MaxErr : Double;
begin
    
    //
    // Demonstration
    //
    Write(Format('POLYNOMIAL INTERPOLATION'#13#10''#13#10'',[]));
    Write(Format('F(x)=sin(x), [0, pi]'#13#10'',[]));
    Write(Format('Second degree polynomial is used'#13#10''#13#10'',[]));
    
    //
    // Create polynomial interpolant
    //
    N := 3;
    SetLength(Y, N);
    I:=0;
    while I<=N-1 do
    begin
        Y[I] := Sin(Pi*I/(N-1));
        Inc(I);
    end;
    PolynomialBuildEqDist(0, Pi, Y, N, P);
    
    //
    // Output results
    //
    BarycentricDiff2(P, 0, V, DV, D2V);
    Write(Format('                 P(x)    F(x) '#13#10'',[]));
    Write(Format('function       %6.3f  %6.3f '#13#10'',[
        BarycentricCalc(P, 0),
        Double(Variant(0))]));
    Write(Format('d/dx(0)        %6.3f  %6.3f '#13#10'',[
        DV,
        Double(Variant(1))]));
    Write(Format('d2/dx2(0)      %6.3f  %6.3f '#13#10'',[
        D2V,
        Double(Variant(0))]));
    Write(Format(''#13#10''#13#10'',[]));
end.