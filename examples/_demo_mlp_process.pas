
program _demo;
Array[0]
var
    Net : MultiLayerPerceptron;
    X : TReal1DArray;
    Y : TReal1DArray;
begin
    
    //
    // regression task with 2 inputs (independent variables)
    // and 2 outputs (dependent variables).
    //
    // network weights are initialized with small random values.
    //
    MLPCreate0(2, 2, Net);
    SetLength(X, 2);
    SetLength(Y, 2);
    X[0] := RandomReal-0.5;
    X[1] := RandomReal-0.5;
    MLPProcess(Net, X, Y);
    Write(Format('Regression task'#13#10'',[]));
    Write(Format('IN[0]  = %5.2f'#13#10'',[
        X[0]]));
    Write(Format('IN[1]  = %5.2f'#13#10'',[
        X[1]]));
    Write(Format('OUT[0] = %5.2f'#13#10'',[
        Y[0]]));
    Write(Format('OUT[1] = %5.2f'#13#10'',[
        Y[1]]));
end.