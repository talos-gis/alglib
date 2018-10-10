unit testminlbfgsunit;
interface
uses Math, Sysutils, Ap, linmin, minlbfgs;

function TestMinLBFGS(Silent : Boolean):Boolean;
function testminlbfgsunit_test_silent():Boolean;
function testminlbfgsunit_test():Boolean;

implementation

procedure TestFunc1(var State : MinLBFGSState);forward;
procedure TestFunc2(var State : MinLBFGSState);forward;
procedure TestFunc3(var State : MinLBFGSState);forward;


function TestMinLBFGS(Silent : Boolean):Boolean;
var
    WasErrors : Boolean;
    RefError : Boolean;
    NonConvError : Boolean;
    EqError : Boolean;
    ConvError : Boolean;
    CrashTest : Boolean;
    OtherErrors : Boolean;
    N : AlglibInteger;
    M : AlglibInteger;
    X : TReal1DArray;
    XE : TReal1DArray;
    B : TReal1DArray;
    XLast : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    A : TReal2DArray;
    MaxIts : AlglibInteger;
    State : MinLBFGSState;
    Rep : MinLBFGSReport;
    FPrev : Double;
    XPrev : Double;
    StpMax : Double;
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
    MinLBFGSCreate(N, M, X, State);
    MinLBFGSSetCond(State, 0, 0, 0, 0);
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
    // nonconvex problems with hard relief: we start from point with very small
    // gradient, but we need ever smaller gradient in the next step due to
    // Wolfe conditions.
    //
    NonConvError := False;
    SetLength(X, 1);
    N := 1;
    M := 1;
    V := -100;
    while AP_FP_Less(V,0.1) do
    begin
        X[0] := V;
        MinLBFGSCreate(N, M, X, State);
        MinLBFGSSetCond(State, 1.0E-9, 0, 0, 0);
        while MinLBFGSIteration(State) do
        begin
            State.F := AP_Sqr(State.X[0])/(1+AP_Sqr(State.X[0]));
            State.G[0] := (2*State.X[0]*(1+AP_Sqr(State.X[0]))-AP_Sqr(State.X[0])*2*State.X[0])/AP_Sqr(1+AP_Sqr(State.X[0]));
        end;
        MinLBFGSResults(State, X, Rep);
        NonConvError := NonConvError or (Rep.TerminationType<=0) or AP_FP_Greater(AbsReal(X[0]),0.001);
        V := V+0.1;
    end;
    
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
            MinLBFGSCreate(N, M, X, State);
            MinLBFGSSetCond(State, 0, 0, 0, 0);
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
    MinLBFGSCreate(N, M, X, State);
    MinLBFGSSetCond(State, 0.001, 0, 0, 0);
    while MinLBFGSIteration(State) do
    begin
        TestFunc3(State);
    end;
    MinLBFGSResults(State, X, Rep);
    ConvError := ConvError or (Rep.TerminationType<>4);
    I:=0;
    while I<=2 do
    begin
        X[I] := 6*RandomReal-3;
        Inc(I);
    end;
    MinLBFGSCreate(N, M, X, State);
    MinLBFGSSetCond(State, 0, 0.001, 0, 0);
    while MinLBFGSIteration(State) do
    begin
        TestFunc3(State);
    end;
    MinLBFGSResults(State, X, Rep);
    ConvError := ConvError or (Rep.TerminationType<>1);
    I:=0;
    while I<=2 do
    begin
        X[I] := 6*RandomReal-3;
        Inc(I);
    end;
    MinLBFGSCreate(N, M, X, State);
    MinLBFGSSetCond(State, 0, 0, 0.001, 0);
    while MinLBFGSIteration(State) do
    begin
        TestFunc3(State);
    end;
    MinLBFGSResults(State, X, Rep);
    ConvError := ConvError or (Rep.TerminationType<>2);
    I:=0;
    while I<=2 do
    begin
        X[I] := 2*RandomReal-1;
        Inc(I);
    end;
    MinLBFGSCreate(N, M, X, State);
    MinLBFGSSetCond(State, 0, 0, 0, 10);
    while MinLBFGSIteration(State) do
    begin
        TestFunc3(State);
    end;
    MinLBFGSResults(State, X, Rep);
    ConvError := ConvError or (Rep.TerminationType<>5) or (Rep.IterationsCount<>10);
    
    //
    // Crash test: too many iterations on a simple tasks
    // May fail when encounter zero step, underflow or something like that
    //
    CrashTest := False;
    SetLength(X, 2+1);
    N := 3;
    M := 2;
    MaxIts := 10000;
    I:=0;
    while I<=2 do
    begin
        X[I] := 6*RandomReal-3;
        Inc(I);
    end;
    MinLBFGSCreate(N, M, X, State);
    MinLBFGSSetCond(State, 0, 0, 0, MaxIts);
    while MinLBFGSIteration(State) do
    begin
        State.F := AP_Sqr(Exp(State.X[0])-2)+AP_Sqr(State.X[1])+AP_Sqr(State.X[2]-State.X[0]);
        State.G[0] := 2*(Exp(State.X[0])-2)*Exp(State.X[0])+2*(State.X[0]-State.X[2]);
        State.G[1] := 2*State.X[1];
        State.G[2] := 2*(State.X[2]-State.X[0]);
    end;
    MinLBFGSResults(State, X, Rep);
    CrashTest := CrashTest or (Rep.TerminationType<=0);
    
    //
    // Other properties:
    // 1. test reports (F should form monotone sequence)
    // 2. test maximum step
    //
    OtherErrors := False;
    N := 50;
    M := 2;
    SetLength(X, N);
    SetLength(XLast, N);
    I:=0;
    while I<=N-1 do
    begin
        X[I] := 1;
        Inc(I);
    end;
    MinLBFGSCreate(N, M, X, State);
    MinLBFGSSetCond(State, 0, 0, 0, 100);
    MinLBFGSSetXRep(State, True);
    FPrev := MaxRealNumber;
    while MinLBFGSIteration(State) do
    begin
        if State.NeedFG then
        begin
            State.F := 0;
            I:=0;
            while I<=N-1 do
            begin
                State.F := State.F+AP_Sqr((1+I)*State.X[I]);
                State.G[I] := 2*(1+I)*State.X[I];
                Inc(I);
            end;
        end;
        if State.XUpdated then
        begin
            OtherErrors := OtherErrors or AP_FP_Greater(State.F,FPrev);
            if AP_FP_Eq(FPrev,MaxRealNumber) then
            begin
                I:=0;
                while I<=N-1 do
                begin
                    OtherErrors := OtherErrors or AP_FP_Neq(State.X[I],X[I]);
                    Inc(I);
                end;
            end;
            FPrev := State.F;
            APVMove(@XLast[0], 0, N-1, @State.X[0], 0, N-1);
        end;
    end;
    MinLBFGSResults(State, X, Rep);
    I:=0;
    while I<=N-1 do
    begin
        OtherErrors := OtherErrors or AP_FP_Neq(X[I],XLast[I]);
        Inc(I);
    end;
    N := 1;
    M := 1;
    SetLength(X, N);
    X[0] := 100;
    StpMax := 0.05+0.05*RandomReal;
    MinLBFGSCreate(N, M, X, State);
    MinLBFGSSetCond(State, 1.0E-9, 0, 0, 0);
    MinLBFGSSetStpMax(State, StpMax);
    MinLBFGSSetXRep(State, True);
    XPrev := X[0];
    while MinLBFGSIteration(State) do
    begin
        if State.NeedFG then
        begin
            State.F := Exp(State.X[0])+Exp(-State.X[0]);
            State.G[0] := Exp(State.X[0])-Exp(-State.X[0]);
            OtherErrors := OtherErrors or AP_FP_Greater(AbsReal(State.X[0]-XPrev),(1+Sqrt(MachineEpsilon))*StpMax);
        end;
        if State.XUpdated then
        begin
            OtherErrors := OtherErrors or AP_FP_Greater(AbsReal(State.X[0]-XPrev),(1+Sqrt(MachineEpsilon))*StpMax);
            XPrev := State.X[0];
        end;
    end;
    
    //
    // end
    //
    WasErrors := RefError or NonConvError or EqError or ConvError or CrashTest or OtherErrors;
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
        Write(Format('NON-CONVEX PROBLEM:                       ',[]));
        if NonConvError then
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
        Write(Format('CRASH TEST:                               ',[]));
        if CrashTest then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('OTHER PROPERTIES:                         ',[]));
        if OtherErrors then
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
Calculate test function #1

It may show very interesting behavior when optimized with 'x[0]>=ln(2)'
constraint.
*************************************************************************)
procedure TestFunc1(var State : MinLBFGSState);
begin
    if AP_FP_Less(State.X[0],100) then
    begin
        State.F := AP_Sqr(Exp(State.X[0])-2)+AP_Sqr(State.X[1])+AP_Sqr(State.X[2]-State.X[0]);
        State.G[0] := 2*(Exp(State.X[0])-2)*Exp(State.X[0])+2*(State.X[0]-State.X[2]);
        State.G[1] := 2*State.X[1];
        State.G[2] := 2*(State.X[2]-State.X[0]);
    end
    else
    begin
        State.F := Sqrt(MaxRealNumber);
        State.G[0] := Sqrt(MaxRealNumber);
        State.G[1] := 0;
        State.G[2] := 0;
    end;
end;


(*************************************************************************
Calculate test function #2

Simple variation of #1, much more nonlinear, which makes unlikely premature
convergence of algorithm .
*************************************************************************)
procedure TestFunc2(var State : MinLBFGSState);
begin
    if AP_FP_Less(State.X[0],100) then
    begin
        State.F := AP_Sqr(Exp(State.X[0])-2)+AP_Sqr(AP_Sqr(State.X[1]))+AP_Sqr(State.X[2]-State.X[0]);
        State.G[0] := 2*(Exp(State.X[0])-2)*Exp(State.X[0])+2*(State.X[0]-State.X[2]);
        State.G[1] := 4*State.X[1]*AP_Sqr(State.X[1]);
        State.G[2] := 2*(State.X[2]-State.X[0]);
    end
    else
    begin
        State.F := Sqrt(MaxRealNumber);
        State.G[0] := Sqrt(MaxRealNumber);
        State.G[1] := 0;
        State.G[2] := 0;
    end;
end;


(*************************************************************************
Calculate test function #3

Simple variation of #1, much more nonlinear, with non-zero value at minimum.
It achieve two goals:
* makes unlikely premature convergence of algorithm .
* solves some issues with EpsF stopping condition which arise when
  F(minimum) is zero

*************************************************************************)
procedure TestFunc3(var State : MinLBFGSState);
var
    S : Double;
begin
    S := 0.001;
    if AP_FP_Less(State.X[0],100) then
    begin
        State.F := AP_Sqr(Exp(State.X[0])-2)+AP_Sqr(AP_Sqr(State.X[1])+S)+AP_Sqr(State.X[2]-State.X[0]);
        State.G[0] := 2*(Exp(State.X[0])-2)*Exp(State.X[0])+2*(State.X[0]-State.X[2]);
        State.G[1] := 2*(AP_Sqr(State.X[1])+S)*2*State.X[1];
        State.G[2] := 2*(State.X[2]-State.X[0]);
    end
    else
    begin
        State.F := Sqrt(MaxRealNumber);
        State.G[0] := Sqrt(MaxRealNumber);
        State.G[1] := 0;
        State.G[2] := 0;
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testminlbfgsunit_test_silent():Boolean;
begin
    Result := TestMinLBFGS(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testminlbfgsunit_test():Boolean;
begin
    Result := TestMinLBFGS(False);
end;


end.