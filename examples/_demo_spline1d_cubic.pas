
program _demo;
Array[0]
var
    X : TReal1DArray;
    Y : TReal1DArray;
    N : AlglibInteger;
    I : AlglibInteger;
    T : Double;
    S : Spline1DInterpolant;
    Err : Double;
    MaxErr : Double;
begin
    
    //
    // Interpolation by natural Cubic spline.
    //
    Write(Format('INTERPOLATION BY NATURAL CUBIC SPLINE'#13#10''#13#10'',[]));
    Write(Format('F(x)=sin(x), [0, pi], 3 nodes'#13#10''#13#10'',[]));
    Write(Format('     x   F(x)   S(x)  Error'#13#10'',[]));
    
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
    Spline1DBuildCubic(X, Y, N, 1, +1, 1, -1, S);
    
    //
    // Output results
    //
    T := 0;
    MaxErr := 0;
    while AP_FP_Less(T,0.999999*Pi) do
    begin
        Err := AbsReal(Spline1DCalc(S, T)-Sin(T));
        MaxErr := Max(Err, MaxErr);
        Write(Format('%6.3f %6.3f %6.3f %6.3f'#13#10'',[
            T,
            Sin(T),
            Spline1DCalc(S, T),
            Err]));
        T := Min(Pi, T+0.25);
    end;
    Err := AbsReal(Spline1DCalc(S, Pi)-Sin(Pi));
    MaxErr := Max(Err, MaxErr);
    Write(Format('%6.3f %6.3f %6.3f %6.3f'#13#10''#13#10'',[
        Pi,
        Sin(Pi),
        Spline1DCalc(S, Pi),
        Err]));
    Write(Format('max|error| = %0.3f'#13#10'',[
        MaxErr]));
    Write(Format('Try other demos (spline1d_linear, spline1d_hermite) and compare errors...'#13#10''#13#10''#13#10'',[]));
end.