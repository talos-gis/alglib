
program _demo;
Array[0]
var
    Network1 : MultiLayerPerceptron;
    Network2 : MultiLayerPerceptron;
    Network3 : MultiLayerPerceptron;
    X : TReal1DArray;
    Y : TReal1DArray;
    R : TReal1DArray;
    RLen : AlglibInteger;
    V1 : Double;
    V2 : Double;
begin
    
    //
    // Generate two networks filled with small random values.
    // Use MLPSerialize/MLPUnserialize to make network copy.
    //
    MLPCreate0(1, 1, Network1);
    MLPCreate0(1, 1, Network2);
    MLPSerialize(Network1, R, RLen);
    MLPUnserialize(R, Network2);
    
    //
    // Now Network1 and Network2 should be identical.
    // Let's demonstrate it.
    //
    Write(Format('Test serialization/unserialization'#13#10'',[]));
    SetLength(X, 1);
    SetLength(Y, 1);
    X[0] := 2*RandomReal-1;
    MLPProcess(Network1, X, Y);
    V1 := Y[0];
    Write(Format('Network1(X) = %0.2f'#13#10'',[
        Y[0]]));
    MLPProcess(Network2, X, Y);
    V2 := Y[0];
    Write(Format('Network2(X) = %0.2f'#13#10'',[
        Y[0]]));
    if AP_FP_Eq(V1,V2) then
    begin
        Write(Format('Results are equal, OK.'#13#10'',[]));
    end
    else
    begin
        Write(Format('Results are not equal... Strange...',[]));
    end;
end.