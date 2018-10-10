unit testlbfgs;
interface
uses Math, Sysutils, Ap, lbfgs;

function TestMinLBFGS(Silent : Boolean):Boolean;
function testlbfgs_test_silent():Boolean;
function testlbfgs_test():Boolean;

implementation

function TestMinLBFGS(Silent : Boolean):Boolean;
var
    WasErrors : Boolean;
    RefError : Boolean;
    EqError : Boolean;
    ConvError : Boolean;
    N : AlglibInteger;
    M : AlglibInteger;
    X : TReal1DArray;
    XE : TReal1DArray;
    B : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    A : TReal2DArray;
    State : LBFGSState;
    Rep : LBFGSReport;
begin
    WasErrors := False;
    
    //
    // Reference problem
    //
    SetLength(X, 2+1);
    N := 3;
    M := 2;
    x[0] := 100*RandomReal-50;
    x[1] := 100*RandomReal-50;
    x[2] := 100*RandomReal-50;
    MinLBFGS(N, M, X, 0.0, 0.0, 0.0, 0, 0, State);
    while MinLBFGSIteration(State) do
    begin
        State.F := AP_Sqr(State.X[0]-2)+AP_Sqr(State.X[1])+AP_Sqr(State.X[2]-State.X[0]);
        State.G[0] := 2*(State.X[0]-2)+2*(State.X[0]-State.X[2]);
        State.G[1] := 2*State.X[1];
        State.G[2] := 2*(State.X[2]-State.X[0]);
    end;
    MinLBFGSResults(State, X, Rep);
    RefError := (Rep.TerminationType<=0) or AP_FP_Greater(AbsReal(X[0]-2),0.001) or AP_FP_Greater(AbsReal(X[1]),0.001) or AP_FP_Greater(AbsReal(X[2]-2),0.001);
    
    //
    // Linear equations
    //
    EqError := False;
    N:=1;
    while N<=10 do
    begin
        
        //
        // Prepare task
        //
        SetLength(A, N-1+1, N-1+1);
        SetLength(X, N-1+1);
        SetLength(XE, N-1+1);
        SetLength(B, N-1+1);
        I:=0;
        while I<=N-1 do
        begin
            XE[I] := 2*RandomReal-1;
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                A[I,J] := 2*RandomReal-1;
                Inc(J);
            end;
            A[I,I] := A[I,I]+3*Sign(A[I,I]);
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            V := APVDotProduct(@A[I][0], 0, N-1, @XE[0], 0, N-1);
            B[I] := V;
            Inc(I);
        end;
        
        //
        // Test different M
        //
        M:=1;
        while M<=N do
        begin
            
            //
            // Solve task
            //
            I:=0;
            while I<=N-1 do
            begin
                X[I] := 2*RandomReal-1;
                Inc(I);
            end;
            MinLBFGS(N, M, X, 0.0, 0.0, 0.0, 0, 0, State);
            while MinLBFGSIteration(State) do
            begin
                State.F := 0;
                I:=0;
                while I<=N-1 do
                begin
                    State.G[I] := 0;
                    Inc(I);
                end;
                I:=0;
                while I<=N-1 do
                begin
                    V := APVDotProduct(@A[I][0], 0, N-1, @State.X[0], 0, N-1);
                    State.F := State.F+AP_Sqr(V-B[I]);
                    J:=0;
                    while J<=N-1 do
                    begin
                        State.G[J] := State.G[J]+2*(V-B[I])*A[I,J];
                        Inc(J);
                    end;
                    Inc(I);
                end;
            end;
            MinLBFGSResults(State, X, Rep);
            EqError := EqError or (Rep.TerminationType<=0);
            I:=0;
            while I<=N-1 do
            begin
                EqError := EqError or AP_FP_Greater(AbsReal(X[I]-XE[I]),0.001);
                Inc(I);
            end;
            Inc(M);
        end;
        Inc(N);
    end;
    
    //
    // Testing convergence properties
    //
    ConvError := False;
    SetLength(X, 2+1);
    N := 3;
    M := 2;
    I:=0;
    while I<=2 do
    begin
        X[I] := 6*RandomReal-3;
        Inc(I);
    end;
    MinLBFGS(N, M, X, 0.0001, 0.0, 0.0, 0, 0, State);
    while MinLBFGSIteration(State) do
    begin
        State.F := AP_Sqr(Exp(State.X[0])-2)+AP_Sqr(State.X[1])+AP_Sqr(State.X[2]-State.X[0]);
        State.G[0] := 2*(Exp(State.X[0])-2)*Exp(State.X[0])+2*(State.X[0]-State.X[2]);
        State.G[1] := 2*State.X[1];
        State.G[2] := 2*(State.X[2]-State.X[0]);
    end;
    MinLBFGSResults(State, X, Rep);
    ConvError := ConvError or AP_FP_Greater(AbsReal(X[0]-Ln(2)),0.05);
    ConvError := ConvError or AP_FP_Greater(AbsReal(X[1]),0.05);
    ConvError := ConvError or AP_FP_Greater(AbsReal(X[2]-Ln(2)),0.05);
    ConvError := ConvError or (Rep.TerminationType<>4);
    I:=0;
    while I<=2 do
    begin
        X[I] := 6*RandomReal-3;
        Inc(I);
    end;
    MinLBFGS(N, M, X, 0.0, 0.0001, 0.0, 0, 0, State);
    while MinLBFGSIteration(State) do
    begin
        State.F := AP_Sqr(Exp(State.X[0])-2)+AP_Sqr(State.X[1])+AP_Sqr(State.X[2]-State.X[0]);
        State.G[0] := 2*(Exp(State.X[0])-2)*Exp(State.X[0])+2*(State.X[0]-State.X[2]);
        State.G[1] := 2*State.X[1];
        State.G[2] := 2*(State.X[2]-State.X[0]);
    end;
    MinLBFGSResults(State, X, Rep);
    ConvError := ConvError or AP_FP_Greater(AbsReal(X[0]-Ln(2)),0.05);
    ConvError := ConvError or AP_FP_Greater(AbsReal(X[1]),0.05);
    ConvError := ConvError or AP_FP_Greater(AbsReal(X[2]-Ln(2)),0.05);
    ConvError := ConvError or (Rep.TerminationType<>1);
    I:=0;
    while I<=2 do
    begin
        X[I] := 6*RandomReal-3;
        Inc(I);
    end;
    MinLBFGS(N, M, X, 0.0, 0.0, 0.0001, 0, 0, State);
    while MinLBFGSIteration(State) do
    begin
        State.F := AP_Sqr(Exp(State.X[0])-2)+AP_Sqr(State.X[1])+AP_Sqr(State.X[2]-State.X[0]);
        State.G[0] := 2*(Exp(State.X[0])-2)*Exp(State.X[0])+2*(State.X[0]-State.X[2]);
        State.G[1] := 2*State.X[1];
        State.G[2] := 2*(State.X[2]-State.X[0]);
    end;
    MinLBFGSResults(State, X, Rep);
    ConvError := ConvError or AP_FP_Greater(AbsReal(X[0]-Ln(2)),0.05);
    ConvError := ConvError or AP_FP_Greater(AbsReal(X[1]),0.05);
    ConvError := ConvError or AP_FP_Greater(AbsReal(X[2]-Ln(2)),0.05);
    ConvError := ConvError or (Rep.TerminationType<>2);
    I:=0;
    while I<=2 do
    begin
        X[I] := 2*RandomReal-1;
        Inc(I);
    end;
    MinLBFGS(N, M, X, 0.0, 0.0, 0.0, 10, 0, State);
    while MinLBFGSIteration(State) do
    begin
        State.F := AP_Sqr(Exp(State.X[0])-2)+AP_Sqr(State.X[1])+AP_Sqr(State.X[2]-State.X[0]);
        State.G[0] := 2*(Exp(State.X[0])-2)*Exp(State.X[0])+2*(State.X[0]-State.X[2]);
        State.G[1] := 2*State.X[1];
        State.G[2] := 2*(State.X[2]-State.X[0]);
    end;
    MinLBFGSResults(State, X, Rep);
    ConvError := ConvError or (Rep.TerminationType<>5) or (Rep.IterationsCount<>10);
    
    //
    // end
    //
    WasErrors := RefError or EqError or ConvError;
    if  not Silent then
    begin
        Write(Format('TESTING L-BFGS OPTIMIZATION'#13#10'',[]));
        Write(Format('REFERENCE PROBLEM:                        ',[]));
        if RefError then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('LINEAR EQUATIONS:                         ',[]));
        if EqError then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('CONVERGENCE PROPERTIES:                   ',[]));
        if ConvError then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        if WasErrors then
        begin
            Write(Format('TEST FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('TEST PASSED'#13#10'',[]));
        end;
        Write(Format(''#13#10''#13#10'',[]));
    end;
    Result :=  not WasErrors;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testlbfgs_test_silent():Boolean;
begin
    Result := TestMinLBFGS(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testlbfgs_test():Boolean;
begin
    Result := TestMinLBFGS(False);
end;


end.