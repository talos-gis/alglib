
program _demo;
Array[0]
var
    State : AutoGKState;
    V : Double;
    Rep : AutoGKReport;
begin
    
    //
    // f(x) = x*sin(x), integrated at [-pi, pi].
    // Exact answer is 2*pi
    //
    AutoGKSmooth(-Pi, +Pi, State);
    while AutoGKIteration(State) do
    begin
        State.F := State.X*Sin(State.X);
    end;
    AutoGKResults(State, V, Rep);
    Write(Format('integral(x*sin(x),-pi,+pi) = %0.2f'#13#10'Exact answer is %0.2f'#13#10'',[
        V,
        2*Pi]));
end.