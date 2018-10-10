unit testasa;
interface
uses Math, Sysutils, Ap, linmin, minasa;

function TestMinASA(Silent : Boolean):Boolean;
function testasa_test_silent():Boolean;
function testasa_test():Boolean;

implementation

procedure TestFunc1(var State : MinASAState);forward;
procedure TestFunc2(var State : MinASAState);forward;
procedure TestFunc3(var State : MinASAState);forward;
procedure CheckBounds(const X : TReal1DArray;
     const BndL : TReal1DArray;
     const BndU : TReal1DArray;
     N : AlglibInteger;
     var Err : Boolean);forward;
function ASABoundVal(X : Double; B1 : Double; B2 : Double):Double;forward;


function TestMinASA(Silent : Boolean):Boolean;
var
    WasErrors : Boolean;
    RefError : Boolean;
    ConvError : Boolean;
    OtherErrors : Boolean;
    N : AlglibInteger;
    X : TReal1DArray;
    XE : TReal1DArray;
    C : TReal1DArray;
    BndL : TReal1DArray;
    BndU : TReal1DArray;
    XLast : TReal1DArray;
    FPrev : Double;
    XPrev : Double;
    StpMax : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    S : Double;
    Tol : Double;
    AlgoType : AlglibInteger;
    A : TReal2DArray;
    State : MinASAState;
    Rep : MinASAReport;
begin
    WasErrors := False;
    RefError := False;
    ConvError := False;
    OtherErrors := False;
    
    //
    // Different algorithms
    //
    AlgoType:=-1;
    while AlgoType<=1 do
    begin
        
        //
        // reference problem, simple convex optimization
        //
        N:=1;
        while N<=5 do
        begin
            
            //
            // min(x'*diag(c)*x) on a random box
            //
            SetLength(X, N);
            SetLength(XE, N);
            SetLength(C, N);
            SetLength(BndL, N);
            SetLength(BndU, N);
            I:=0;
            while I<=N-1 do
            begin
                C[I] := 1+RandomReal;
                XE[I] := 4*RandomReal-2;
                BndL[I] := -Max(RandomReal, 0.2);
                BndU[I] := +Max(RandomReal, 0.2);
                X[I] := 0.5*(BndL[I]+BndU[I]);
                Inc(I);
            end;
            Tol := 0.001;
            MinASACreate(N, X, BndL, BndU, State);
            MinASASetCond(State, Tol, 0.0, 0.0, 0);
            MinASASetAlgorithm(State, AlgoType);
            while MinASAIteration(State) do
            begin
                CheckBounds(State.X, BndL, BndU, N, OtherErrors);
                State.F := 0;
                I:=0;
                while I<=N-1 do
                begin
                    State.F := State.F+C[I]*AP_Sqr(State.X[I]-XE[I]);
                    State.G[I] := 2*C[I]*(State.X[I]-XE[I]);
                    Inc(I);
                end;
            end;
            MinASAResults(State, X, Rep);
            RefError := RefError or (Rep.TerminationType<=0);
            I:=0;
            while I<=N-1 do
            begin
                RefError := RefError or AP_FP_Greater(AbsReal(ASABoundVal(XE[I], BndL[I], BndU[I])-X[I]),0.01);
                Inc(I);
            end;
            Inc(N);
        end;
        
        //
        // reference problem 2: non-convex optimization on [-2,2] x [1,2]
        //
        // A saddle function is minimized:
        // * stationary point [0,0] (non-feasible)
        // * constrained minimum [-2,2].
        // * starting point [+2,2]
        //
        // Path from start to end may be very complex, with multiple changes
        // in active constraints, so it is interesting task for our method.
        //
        // Scale parameter is used to make optimization more interesting
        // during GPA runs.
        //
        SetLength(X, 2);
        SetLength(BndL, 2);
        SetLength(BndU, 2);
        BndL[0] := -2;
        BndU[0] := 2;
        X[0] := 2;
        BndL[1] := 1;
        BndU[1] := 2;
        X[1] := 2;
        Tol := 0.001;
        S := 0.01;
        MinASACreate(2, X, BndL, BndU, State);
        MinASASetCond(State, Tol, 0.0, 0.0, 0);
        MinASASetAlgorithm(State, AlgoType);
        while MinASAIteration(State) do
        begin
            CheckBounds(State.X, BndL, BndU, 2, OtherErrors);
            State.F := S*(AP_Sqr(State.X[0]+State.X[1])-AP_Sqr(State.X[0]-State.X[1]));
            State.G[0] := S*(2*(State.X[0]+State.X[1])-2*(State.X[0]-State.X[1]));
            State.G[1] := S*(2*(State.X[0]+State.X[1])+2*(State.X[0]-State.X[1]));
        end;
        MinASAResults(State, X, Rep);
        RefError := RefError or (Rep.TerminationType<=0) or AP_FP_Greater(AbsReal(State.X[0]+2),0.01) or AP_FP_Greater(AbsReal(State.X[1]-2),0.01);
        
        //
        // function #1 with 'x[0]>=ln(2)' constraint.
        // may show very interesting behavior.
        //
        SetLength(X, 3);
        SetLength(BndL, 3);
        SetLength(BndU, 3);
        N := 3;
        I:=0;
        while I<=2 do
        begin
            BndL[I] := -10000;
            BndU[I] := +10000;
            Inc(I);
        end;
        BndL[0] := Ln(2);
        I:=0;
        while I<=2 do
        begin
            X[I] := 3*RandomReal+3;
            Inc(I);
        end;
        MinASACreate(N, X, BndL, BndU, State);
        MinASASetCond(State, 0.0000001, 0.0, 0.0, 0);
        MinASASetAlgorithm(State, AlgoType);
        while MinASAIteration(State) do
        begin
            CheckBounds(State.X, BndL, BndU, N, OtherErrors);
            TestFunc1(State);
        end;
        MinASAResults(State, X, Rep);
        RefError := RefError or (Rep.TerminationType<=0);
        RefError := RefError or AP_FP_Greater(AbsReal(X[0]-Ln(2)),0.05);
        RefError := RefError or AP_FP_Greater(AbsReal(X[1]),0.05);
        RefError := RefError or AP_FP_Greater(AbsReal(X[2]-Ln(2)),0.05);
        
        //
        // Testing convergence properties
        //
        SetLength(X, 3);
        SetLength(BndL, 3);
        SetLength(BndU, 3);
        N := 3;
        I:=0;
        while I<=2 do
        begin
            BndL[I] := -10000;
            BndU[I] := +10000;
            Inc(I);
        end;
        BndL[0] := Ln(2);
        I:=0;
        while I<=2 do
        begin
            X[I] := 3*RandomReal+3;
            Inc(I);
        end;
        MinASACreate(N, X, BndL, BndU, State);
        MinASASetCond(State, 0.001, 0.0, 0.0, 0);
        MinASASetAlgorithm(State, AlgoType);
        while MinASAIteration(State) do
        begin
            CheckBounds(State.X, BndL, BndU, N, OtherErrors);
            TestFunc3(State);
        end;
        MinASAResults(State, X, Rep);
        ConvError := ConvError or (Rep.TerminationType<>4);
        I:=0;
        while I<=2 do
        begin
            X[I] := 3*RandomReal+3;
            Inc(I);
        end;
        MinASACreate(N, X, BndL, BndU, State);
        MinASASetCond(State, 0.0, 0.001, 0.0, 0);
        MinASASetAlgorithm(State, AlgoType);
        while MinASAIteration(State) do
        begin
            CheckBounds(State.X, BndL, BndU, N, OtherErrors);
            TestFunc3(State);
        end;
        MinASAResults(State, X, Rep);
        ConvError := ConvError or (Rep.TerminationType<>1);
        I:=0;
        while I<=2 do
        begin
            X[I] := 3*RandomReal+3;
            Inc(I);
        end;
        MinASACreate(N, X, BndL, BndU, State);
        MinASASetCond(State, 0.0, 0.0, 0.001, 0);
        MinASASetAlgorithm(State, AlgoType);
        while MinASAIteration(State) do
        begin
            CheckBounds(State.X, BndL, BndU, N, OtherErrors);
            TestFunc3(State);
        end;
        MinASAResults(State, X, Rep);
        ConvError := ConvError or (Rep.TerminationType<>2);
        I:=0;
        while I<=2 do
        begin
            X[I] := 3*RandomReal+3;
            Inc(I);
        end;
        MinASACreate(N, X, BndL, BndU, State);
        MinASASetCond(State, 0.0, 0.0, 0.0, 3);
        MinASASetAlgorithm(State, AlgoType);
        while MinASAIteration(State) do
        begin
            CheckBounds(State.X, BndL, BndU, N, OtherErrors);
            TestFunc3(State);
        end;
        MinASAResults(State, X, Rep);
        ConvError := ConvError or  not ((Rep.TerminationType=5) and (Rep.IterationsCount=3) or (Rep.TerminationType=7));
        
        //
        // Other properties
        //
        //
        // Other properties:
        // 1. test reports (F should form monotone sequence)
        // 2. test maximum step
        //
        N := 50;
        SetLength(X, N);
        SetLength(XLast, N);
        SetLength(BndL, N);
        SetLength(BndU, N);
        I:=0;
        while I<=N-1 do
        begin
            X[I] := 1;
            XLast[I] := RandomReal;
            BndL[I] := -100000;
            BndU[I] := +100000;
            Inc(I);
        end;
        MinASACreate(N, X, BndL, BndU, State);
        MinASASetCond(State, 0, 0, 0, 100);
        MinASASetXRep(State, True);
        FPrev := MaxRealNumber;
        while MinASAIteration(State) do
        begin
            CheckBounds(State.X, BndL, BndU, N, OtherErrors);
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
        MinASAResults(State, X, Rep);
        I:=0;
        while I<=N-1 do
        begin
            OtherErrors := OtherErrors or AP_FP_Neq(X[I],XLast[I]);
            Inc(I);
        end;
        N := 1;
        SetLength(X, N);
        SetLength(BndL, N);
        SetLength(BndU, N);
        X[0] := 100;
        BndL[0] := -1000000;
        BndU[0] := +1000000;
        StpMax := 0.05+0.05*RandomReal;
        MinASACreate(N, X, BndL, BndU, State);
        MinASASetCond(State, 1.0E-9, 0, 0, 0);
        MinASASetStpMax(State, StpMax);
        MinASASetXRep(State, True);
        XPrev := X[0];
        while MinASAIteration(State) do
        begin
            CheckBounds(State.X, BndL, BndU, N, OtherErrors);
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
        Inc(AlgoType);
    end;
    
    //
    // end
    //
    WasErrors := RefError or ConvError or OtherErrors;
    if  not Silent then
    begin
        Write(Format('TESTING ASA OPTIMIZATION'#13#10'',[]));
        Write(Format('REFERENCE PROBLEMS:                       ',[]));
        if RefError then
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
procedure TestFunc1(var State : MinASAState);
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
procedure TestFunc2(var State : MinASAState);
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
procedure TestFunc3(var State : MinASAState);
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
Checks that X is bounded with respect to BndL/BndU.

If it is not, True is assigned to the Err variable (which is not changed
otherwise).
*************************************************************************)
procedure CheckBounds(const X : TReal1DArray;
     const BndL : TReal1DArray;
     const BndU : TReal1DArray;
     N : AlglibInteger;
     var Err : Boolean);
var
    I : AlglibInteger;
begin
    I:=0;
    while I<=N-1 do
    begin
        if AP_FP_Less(X[I],BndL[I]) or AP_FP_Greater(X[I],BndU[I]) then
        begin
            Err := True;
        end;
        Inc(I);
    end;
end;


(*************************************************************************
'bound' value: map X to [B1,B2]
*************************************************************************)
function ASABoundVal(X : Double; B1 : Double; B2 : Double):Double;
begin
    if AP_FP_Less_Eq(X,B1) then
    begin
        Result := B1;
        Exit;
    end;
    if AP_FP_Greater_Eq(X,B2) then
    begin
        Result := B2;
        Exit;
    end;
    Result := X;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testasa_test_silent():Boolean;
begin
    Result := TestMinASA(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testasa_test():Boolean;
begin
    Result := TestMinASA(False);
end;


end.