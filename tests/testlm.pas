unit testlm;
interface
uses Math, Sysutils, Ap, blas, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, hblas, sblas, ortfac, rotations, bdsvd, svd, xblas, densesolver, lbfgs, minlm;

function TestMinLM(Silent : Boolean):Boolean;
function testlm_test_silent():Boolean;
function testlm_test():Boolean;

implementation

function RKindVsStateCheck(RKind : AlglibInteger;
     const State : LMState):Boolean;forward;


function TestMinLM(Silent : Boolean):Boolean;
var
    WasErrors : Boolean;
    RefError : Boolean;
    Lin1Error : Boolean;
    Lin2Error : Boolean;
    EqError : Boolean;
    ConvError : Boolean;
    SCError : Boolean;
    RKind : AlglibInteger;
    CKind : AlglibInteger;
    EpsF : Double;
    EpsX : Double;
    MaxIts : AlglibInteger;
    N : AlglibInteger;
    M : AlglibInteger;
    X : TReal1DArray;
    XE : TReal1DArray;
    B : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    A : TReal2DArray;
    State : LMState;
    Rep : LMReport;
begin
    WasErrors := False;
    RefError := False;
    Lin1Error := False;
    Lin2Error := False;
    EqError := False;
    ConvError := False;
    SCError := False;
    
    //
    // Reference problem.
    // RKind is a algorithm selector:
    // * 0 = FJ
    // * 1 = FGJ
    // * 2 = FGH
    //
    SetLength(X, 2+1);
    N := 3;
    M := 3;
    RKind:=0;
    while RKind<=2 do
    begin
        X[0] := 100*RandomReal-50;
        X[1] := 100*RandomReal-50;
        X[2] := 100*RandomReal-50;
        if RKind=0 then
        begin
            MinLMFJ(N, M, X, 0.0, 0.0, 0, State);
        end;
        if RKind=1 then
        begin
            MinLMFGJ(N, M, X, 0.0, 0.0, 0, State);
        end;
        if RKind=2 then
        begin
            MinLMFGH(N, X, 0.0, 0.0, 0, State);
        end;
        while MinLMIteration(State) do
        begin
            
            //
            // (x-2)^2 + y^2 + (z-x)^2
            //
            State.F := AP_Sqr(State.X[0]-2)+AP_Sqr(State.X[1])+AP_Sqr(State.X[2]-State.X[0]);
            if State.NeedFG or State.NeedFGH then
            begin
                State.G[0] := 2*(State.X[0]-2)+2*(State.X[0]-State.X[2]);
                State.G[1] := 2*State.X[1];
                State.G[2] := 2*(State.X[2]-State.X[0]);
            end;
            if State.NeedFiJ then
            begin
                State.Fi[0] := State.X[0]-2;
                State.Fi[1] := State.X[1];
                State.Fi[2] := State.X[2]-State.X[0];
                State.J[0,0] := 1;
                State.J[0,1] := 0;
                State.J[0,2] := 0;
                State.J[1,0] := 0;
                State.J[1,1] := 1;
                State.J[1,2] := 0;
                State.J[2,0] := -1;
                State.J[2,1] := 0;
                State.J[2,2] := 1;
            end;
            if State.NeedFGH then
            begin
                State.H[0,0] := 4;
                State.H[0,1] := 0;
                State.H[0,2] := -2;
                State.H[1,0] := 0;
                State.H[1,1] := 2;
                State.H[1,2] := 0;
                State.H[2,0] := -2;
                State.H[2,1] := 0;
                State.H[2,2] := 2;
            end;
            SCError := SCError or  not RKindVsStateCheck(RKind, State);
        end;
        MinLMResults(State, X, Rep);
        RefError := RefError or (Rep.TerminationType<=0) or AP_FP_Greater(AbsReal(X[0]-2),0.001) or AP_FP_Greater(AbsReal(X[1]),0.001) or AP_FP_Greater(AbsReal(X[2]-2),0.001);
        Inc(RKind);
    end;
    
    //
    // 1D problem #1
    //
    RKind:=0;
    while RKind<=2 do
    begin
        SetLength(X, 1);
        N := 1;
        M := 1;
        x[0] := 100*RandomReal-50;
        if RKind=0 then
        begin
            MinLMFJ(N, M, X, 0.0, 0.0, 0, State);
        end;
        if RKind=1 then
        begin
            MinLMFGJ(N, M, X, 0.0, 0.0, 0, State);
        end;
        if RKind=2 then
        begin
            MinLMFGH(N, X, 0.0, 0.0, 0, State);
        end;
        while MinLMIteration(State) do
        begin
            State.F := AP_Sqr(Sin(State.X[0]));
            if State.NeedFG or State.NeedFGH then
            begin
                State.G[0] := 2*Sin(State.X[0])*Cos(State.X[0]);
            end;
            if State.NeedFiJ then
            begin
                State.Fi[0] := Sin(State.X[0]);
                State.J[0,0] := Cos(State.X[0]);
            end;
            if State.NeedFGH then
            begin
                State.H[0,0] := 2*(Cos(State.X[0])*Cos(State.X[0])-Sin(State.X[0])*Sin(State.X[0]));
            end;
            SCError := SCError or  not RKindVsStateCheck(RKind, State);
        end;
        MinLMResults(State, X, Rep);
        Lin1Error := (Rep.TerminationType<=0) or AP_FP_Greater(AbsReal(X[0]/Pi-Round(X[0]/Pi)),0.001);
        Inc(RKind);
    end;
    
    //
    // Linear equations
    //
    N:=1;
    while N<=10 do
    begin
        
        //
        // Prepare task
        //
        RMatrixRndCond(N, 100, A);
        SetLength(X, N);
        SetLength(XE, N);
        SetLength(B, N);
        I:=0;
        while I<=N-1 do
        begin
            XE[I] := 2*RandomReal-1;
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
        // Test different RKind
        //
        RKind:=0;
        while RKind<=2 do
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
            if RKind=0 then
            begin
                MinLMFJ(N, N, X, 0.0, 0.0, 0, State);
            end;
            if RKind=1 then
            begin
                MinLMFGJ(N, N, X, 0.0, 0.0, 0, State);
            end;
            if RKind=2 then
            begin
                MinLMFGH(N, X, 0.0, 0.0, 0, State);
            end;
            while MinLMIteration(State) do
            begin
                if State.NeedF or State.NeedFG or State.NeedFGH then
                begin
                    State.F := 0;
                end;
                if State.NeedFG or State.NeedFGH then
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        State.G[I] := 0;
                        Inc(I);
                    end;
                end;
                if State.NeedFGH then
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            State.H[I,J] := 0;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                end;
                I:=0;
                while I<=N-1 do
                begin
                    V := APVDotProduct(@A[I][0], 0, N-1, @State.X[0], 0, N-1);
                    if State.NeedF or State.NeedFG or State.NeedFGH then
                    begin
                        State.F := State.F+AP_Sqr(V-B[I]);
                    end;
                    if State.NeedFG or State.NeedFGH then
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            State.G[J] := State.G[J]+2*(V-B[I])*A[I,J];
                            Inc(J);
                        end;
                    end;
                    if State.NeedFGH then
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            K:=0;
                            while K<=N-1 do
                            begin
                                State.H[J,K] := State.H[J,K]+2*A[I,J]*A[I,K];
                                Inc(K);
                            end;
                            Inc(J);
                        end;
                    end;
                    if State.NeedFiJ then
                    begin
                        State.Fi[I] := V-B[I];
                        APVMove(@State.J[I][0], 0, N-1, @A[I][0], 0, N-1);
                    end;
                    Inc(I);
                end;
                SCError := SCError or  not RKindVsStateCheck(RKind, State);
            end;
            MinLMResults(State, X, Rep);
            EqError := EqError or (Rep.TerminationType<=0);
            I:=0;
            while I<=N-1 do
            begin
                EqError := EqError or AP_FP_Greater(AbsReal(X[I]-XE[I]),0.001);
                Inc(I);
            end;
            Inc(RKind);
        end;
        Inc(N);
    end;
    
    //
    // Testing convergence properties using
    // different optimizer types and different conditions
    //
    RKind:=0;
    while RKind<=2 do
    begin
        CKind:=0;
        while CKind<=2 do
        begin
            EpsF := 0;
            EpsX := 0;
            MaxIts := 0;
            if CKind=0 then
            begin
                EpsF := 0.0001;
            end;
            if CKind=1 then
            begin
                EpsX := 0.0001;
            end;
            if CKind=2 then
            begin
                MaxIts := 2;
            end;
            SetLength(X, 3);
            N := 3;
            M := 3;
            I:=0;
            while I<=2 do
            begin
                X[I] := 6;
                Inc(I);
            end;
            if RKind=0 then
            begin
                MinLMFJ(N, M, X, EpsF, EpsX, MaxIts, State);
            end;
            if RKind=1 then
            begin
                MinLMFGJ(N, M, X, EpsF, EpsX, MaxIts, State);
            end;
            if RKind=2 then
            begin
                MinLMFGH(N, X, EpsF, EpsX, MaxIts, State);
            end;
            while MinLMiteration(State) do
            begin
                if State.NeedF or State.NeedFG or State.NeedFGH then
                begin
                    State.F := AP_Sqr(Exp(State.X[0])-2)+AP_Sqr(State.X[1])+AP_Sqr(State.X[2]-State.X[0]);
                end;
                if State.NeedFG or State.NeedFGH then
                begin
                    State.G[0] := 2*(Exp(State.X[0])-2)*Exp(State.X[0])+2*(State.X[0]-State.X[2]);
                    State.G[1] := 2*State.X[1];
                    State.G[2] := 2*(State.X[2]-State.X[0]);
                end;
                if State.NeedFGH then
                begin
                    State.H[0,0] := 4*AP_Sqr(Exp(State.X[0]))-4*Exp(State.X[0])+2;
                    State.H[0,1] := 0;
                    State.H[0,2] := -2;
                    State.H[1,0] := 0;
                    State.H[1,1] := 2;
                    State.H[1,2] := 0;
                    State.H[2,0] := -2;
                    State.H[2,1] := 0;
                    State.H[2,2] := 2;
                end;
                if State.NeedFiJ then
                begin
                    State.Fi[0] := Exp(State.X[0])-2;
                    State.J[0,0] := Exp(State.X[0]);
                    State.J[0,1] := 0;
                    State.J[0,2] := 0;
                    State.Fi[1] := State.X[1];
                    State.J[1,0] := 0;
                    State.J[1,1] := 1;
                    State.J[1,2] := 0;
                    State.Fi[2] := State.X[2]-State.X[0];
                    State.J[2,0] := -1;
                    State.J[2,1] := 0;
                    State.J[2,2] := 1;
                end;
                SCError := SCError or  not RKindVsStateCheck(RKind, State);
            end;
            MinLMResults(State, X, Rep);
            if CKind=0 then
            begin
                ConvError := ConvError or AP_FP_Greater(AbsReal(X[0]-Ln(2)),0.05);
                ConvError := ConvError or AP_FP_Greater(AbsReal(X[1]),0.05);
                ConvError := ConvError or AP_FP_Greater(AbsReal(X[2]-Ln(2)),0.05);
                ConvError := ConvError or (Rep.TerminationType<>1);
            end;
            if CKind=1 then
            begin
                ConvError := ConvError or AP_FP_Greater(AbsReal(X[0]-Ln(2)),0.05);
                ConvError := ConvError or AP_FP_Greater(AbsReal(X[1]),0.05);
                ConvError := ConvError or AP_FP_Greater(AbsReal(X[2]-Ln(2)),0.05);
                ConvError := ConvError or (Rep.TerminationType<>2);
            end;
            if CKind=2 then
            begin
                ConvError := ConvError or (Rep.TerminationType<>5) or (Rep.IterationsCount<>MaxIts);
            end;
            Inc(CKind);
        end;
        Inc(RKind);
    end;
    
    //
    // end
    //
    WasErrors := RefError or Lin1Error or Lin2Error or EqError or ConvError or SCError;
    if  not Silent then
    begin
        Write(Format('TESTING LEVENBERG-MARQUARDT OPTIMIZATION'#13#10'',[]));
        Write(Format('REFERENCE PROBLEM:                        ',[]));
        if RefError then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('1-D PROBLEM #1:                           ',[]));
        if Lin1Error then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('1-D PROBLEM #2:                           ',[]));
        if Lin2Error then
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
        Write(Format('STATE FIELDS CONSISTENCY:                 ',[]));
        if SCError then
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
Asserts that State fields are consistent with RKind.
Returns False otherwise.
*************************************************************************)
function RKindVsStateCheck(RKind : AlglibInteger;
     const State : LMState):Boolean;
var
    NSet : AlglibInteger;
begin
    NSet := 0;
    if State.NeedF then
    begin
        NSet := NSet+1;
    end;
    if State.NeedFG then
    begin
        NSet := NSet+1;
    end;
    if State.NeedFiJ then
    begin
        NSet := NSet+1;
    end;
    if State.NeedFGH then
    begin
        NSet := NSet+1;
    end;
    if NSet<>1 then
    begin
        Result := False;
        Exit;
    end;
    if (RKind=0) and (State.NeedFG or State.NeedFGH) then
    begin
        Result := False;
        Exit;
    end;
    if (RKind=1) and State.NeedFGH then
    begin
        Result := False;
        Exit;
    end;
    if (RKind=2) and State.NeedFiJ then
    begin
        Result := False;
        Exit;
    end;
    Result := True;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testlm_test_silent():Boolean;
begin
    Result := TestMinLM(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testlm_test():Boolean;
begin
    Result := TestMinLM(False);
end;


end.