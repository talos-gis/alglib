
program _demo;
Array[0]
var
    X : TReal1DArray;
    Y : TReal1DArray;
    N : AlglibInteger;
    I : AlglibInteger;
    Info : AlglibInteger;
    S : Spline1DInterpolant;
    T : Double;
    Rep : Spline1DFitReport;
begin
    
    //
    // Fitting by unconstrained natural cubic spline
    //
    Write(Format('FITTING BY UNCONSTRAINED NATURAL CUBIC SPLINE'#13#10''#13#10'',[]));
    Write(Format('F(x)=sin(x)      function being fitted'#13#10'',[]));
    Write(Format('[0, pi]          interval'#13#10'',[]));
    Write(Format('M=4              number of basis functions to use'#13#10'',[]));
    Write(Format('N=100            number of points to fit'#13#10'',[]));
    
    //
    // Create and fit
    //
    N := 100;
    SetLength(X, N);
    SetLength(Y, N);
    I:=0;
    while I<=N-1 do
    begin
        X[I] := Pi*I/(N-1);
        Y[I] := Sin(X[I]);
        Inc(I);
    end;
    Spline1DFitCubic(X, Y, N, 4, Info, S, Rep);
    
    //
    // Output results
    //
    if Info>0 then
    begin
        Write(Format(''#13#10'OK, we have finished'#13#10''#13#10'',[]));
        Write(Format('     x   F(x)   S(x)  Error'#13#10'',[]));
        T := 0;
        while AP_FP_Less(T,0.999999*Pi) do
        begin
            Write(Format('%6.3f %6.3f %6.3f %6.3f'#13#10'',[
                T,
                Sin(T),
                Spline1DCalc(S, T),
                AbsReal(Spline1DCalc(S, T)-Sin(T))]));
            T := Min(Pi, T+0.25);
        end;
        Write(Format('%6.3f %6.3f %6.3f %6.3f'#13#10''#13#10'',[
            T,
            Sin(T),
            Spline1DCalc(S, T),
            AbsReal(Spline1DCalc(S, T)-Sin(T))]));
        Write(Format('rms error is %6.3f'#13#10'',[
            Rep.RMSError]));
        Write(Format('max error is %6.3f'#13#10'',[
            Rep.MaxError]));
    end
    else
    begin
        Write(Format(''#13#10'Something wrong, Info=%0d',[
            Info]));
    end;
end.