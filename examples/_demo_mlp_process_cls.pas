
program _demo;
Array[0]
var
    Net : MultiLayerPerceptron;
    X : TReal1DArray;
    Y : TReal1DArray;
begin
    
    //
    // classification task with 2 inputs and 3 classes.
    //
    // network weights are initialized with small random values.
    //
    MLPCreateC0(2, 3, Net);
    SetLength(X, 2);
    SetLength(Y, 3);
    X[0] := RandomReal-0.5;
    X[1] := RandomReal-0.5;
    MLPProcess(Net, X, Y);
    
    //
    // output results
    //
    Write(Format('Classification task'#13#10'',[]));
    Write(Format('IN[0]  = %5.2f'#13#10'',[
        X[0]]));
    Write(Format('IN[1]  = %5.2f'#13#10'',[
        X[1]]));
    Write(Format('Prob(Class=0|IN) = %5.2f'#13#10'',[
        Y[0]]));
    Write(Format('Prob(Class=1|IN) = %5.2f'#13#10'',[
        Y[1]]));
    Write(Format('Prob(Class=2|IN) = %5.2f'#13#10'',[
        Y[2]]));
end.