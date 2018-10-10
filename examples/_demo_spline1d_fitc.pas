
program _demo;
Array[0]
var
    X : TReal1DArray;
    Y : TReal1DArray;
    W : TReal1DArray;
    XC : TReal1DArray;
    YC : TReal1DArray;
    DC : TInteger1DArray;
    N : AlglibInteger;
    I : AlglibInteger;
    Info : AlglibInteger;
    S : Spline1DInterpolant;
    T : Double;
    Rep : Spline1DFitReport;
begin
    
    //
    // Fitting by constrained Hermite spline
    //
    Write(Format('FITTING BY CONSTRAINED HERMITE SPLINE'#13#10''#13#10'',[]));
    Write(Format('F(x)=sin(x)      function being fitted'#13#10'',[]));
    Write(Format('[0, pi]          interval'#13#10'',[]));
    Write(Format('M=6              number of basis functions to use'#13#10'',[]));
    Write(Format('S(0)=0           first constraint'#13#10'',[]));
    Write(Format('S(pi)=0          second constraint'#13#10'',[]));
    Write(Format('N=100            number of points to fit'#13#10'',[]));
    
    //
    // Create and fit:
    // * X  contains points
    // * Y  contains values
    // * W  contains weights
    // * XC contains constraints locations
    // * YC contains constraints values
    // * DC contains derivative indexes (0 = constrained function value)
    //
    N := 100;
    SetLength(X, N);
    SetLength(Y, N);
    SetLength(W, N);
    I:=0;
    while I<=N-1 do
    begin
        X[I] := Pi*I/(N-1);
        Y[I] := Sin(X[I]);
        W[I] := 1;
        Inc(I);
    end;
    SetLength(XC, 2);
    SetLength(YC, 2);
    SetLength(DC, 2);
    XC[0] := 0;
    YC[0] := 0;
    DC[0] := 0;
    XC[0] := Pi;
    YC[0] := 0;
    DC[0] := 0;
    Spline1DFitHermiteWC(X, Y, W, N, XC, YC, DC, 2, 6, Info, S, Rep);
    
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
        Write(Format('S(0) = S(pi) = 0 (exactly)'#13#10''#13#10'',[]));
    end
    else
    begin
        Write(Format(''#13#10'Something wrong, Info=%0d',[
            Info]));
    end;
end.