
program _demo;
Array[0]
var
    State : AutoGKState;
    V : Double;
    Rep : AutoGKReport;
    A : Double;
    B : Double;
    Alpha : Double;
begin
    
    //
    // f1(x) = (1+x)*(x-a)^alpha, alpha=-0.3
    // Exact answer is (B-A)^(Alpha+2)/(Alpha+2) + (1+A)*(B-A)^(Alpha+1)/(Alpha+1)
    //
    // This code demonstrates use of the State.XMinusA (State.BMinusX) field.
    //
    // If we try to use State.X instead of State.XMinusA,
    // we will end up dividing by zero! (in 64-bit precision)
    //
    A := 1.0;
    B := 5.0;
    Alpha := -0.9;
    AutoGKSingular(A, B, Alpha, 0.0, State);
    while AutoGKIteration(State) do
    begin
        State.F := Power(State.XMinusA, Alpha)*(1+State.X);
    end;
    AutoGKResults(State, V, Rep);
    Write(Format('integral((1+x)*(x-a)^alpha) on [%0.1f; %0.1f] = %0.2f'#13#10'Exact answer is %0.2f'#13#10'',[
        A,
        B,
        V,
        Power(B-A, Alpha+2)/(Alpha+2)+(1+A)*Power(B-A, Alpha+1)/(Alpha+1)]));
end.