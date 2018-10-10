unit llstestunit;
interface
uses Math, Sysutils, Ap, blas, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, hblas, sblas, ortfac, rotations, bdsvd, svd, xblas, densesolver, lbfgs, minlm, lsfit;

function TestLLS(Silent : Boolean):Boolean;
function llstestunit_test_silent():Boolean;
function llstestunit_test():Boolean;

implementation

function IsGLSSolution(N : AlglibInteger;
     M : AlglibInteger;
     K : AlglibInteger;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     const FMatrix : TReal2DArray;
     const CMatrix : TReal2DArray;
     C : TReal1DArray):Boolean;forward;
function GetGLSError(N : AlglibInteger;
     M : AlglibInteger;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     const FMatrix : TReal2DArray;
     const C : TReal1DArray):Double;forward;
procedure FitLinearNonlinear(M : AlglibInteger;
     GradOnly : Boolean;
     const XY : TReal2DArray;
     var State : LSFitState;
     var NLSErrors : Boolean);forward;


function TestLLS(Silent : Boolean):Boolean;
var
    WasErrors : Boolean;
    LLSErrors : Boolean;
    NLSErrors : Boolean;
    Threshold : Double;
    NLThreshold : Double;
    MaxN : AlglibInteger;
    MaxM : AlglibInteger;
    PassCount : AlglibInteger;
    N : AlglibInteger;
    M : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    Pass : AlglibInteger;
    XScale : Double;
    X : TReal1DArray;
    Y : TReal1DArray;
    W : TReal1DArray;
    W2 : TReal1DArray;
    C : TReal1DArray;
    C2 : TReal1DArray;
    A : TReal2DArray;
    A2 : TReal2DArray;
    CM : TReal2DArray;
    V : Double;
    V1 : Double;
    V2 : Double;
    Rep : LSFitReport;
    Rep2 : LSFitReport;
    Info : AlglibInteger;
    Info2 : AlglibInteger;
    RefRMS : Double;
    RefAvg : Double;
    RefAvgRel : Double;
    RefMax : Double;
    State : LSFitState;
begin
    WasErrors := False;
    LLSErrors := False;
    NLSErrors := False;
    Threshold := 10000*MachineEpsilon;
    NLThreshold := 0.00001;
    MaxN := 6;
    MaxM := 6;
    PassCount := 4;
    
    //
    // Testing unconstrained least squares (linear/nonlinear)
    //
    N:=1;
    while N<=MaxN do
    begin
        M:=1;
        while M<=MaxM do
        begin
            Pass:=1;
            while Pass<=PassCount do
            begin
                
                //
                // Solve non-degenerate linear least squares task
                // Use Chebyshev basis. Its condition number is very good.
                //
                SetLength(A, N, M);
                SetLength(X, N);
                SetLength(Y, N);
                SetLength(W, N);
                XScale := 0.9+0.1*RandomReal;
                I:=0;
                while I<=N-1 do
                begin
                    if N=1 then
                    begin
                        X[I] := 2*RandomReal-1;
                    end
                    else
                    begin
                        X[I] := XScale*(AP_Double(2*I)/(N-1)-1);
                    end;
                    Y[I] := 3*X[I]+Exp(X[I]);
                    W[I] := 1+RandomReal;
                    A[I,0] := 1;
                    if M>1 then
                    begin
                        A[I,1] := X[I];
                    end;
                    J:=2;
                    while J<=M-1 do
                    begin
                        A[I,J] := 2*X[I]*A[I,J-1]-A[I,J-2];
                        Inc(J);
                    end;
                    Inc(I);
                end;
                
                //
                // 1. test weighted fitting (optimality)
                // 2. Solve degenerate least squares task built on the basis
                //    of previous task
                //
                LSFitLinearW(Y, W, A, N, M, Info, C, Rep);
                if Info<=0 then
                begin
                    LLSErrors := True;
                end
                else
                begin
                    LLSErrors := LLSErrors or  not IsGLSSolution(N, M, 0, Y, W, A, CM, C);
                end;
                SetLength(A2, N, 2*M);
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=M-1 do
                    begin
                        A2[I,2*J+0] := A[I,J];
                        A2[I,2*J+1] := A[I,J];
                        Inc(J);
                    end;
                    Inc(I);
                end;
                LSFitLinearW(Y, W, A2, N, 2*M, Info, C2, Rep);
                if Info<=0 then
                begin
                    LLSErrors := True;
                end
                else
                begin
                    
                    //
                    // test answer correctness using design matrix properties
                    // and previous task solution
                    //
                    J:=0;
                    while J<=M-1 do
                    begin
                        LLSErrors := LLSErrors or AP_FP_Greater(AbsReal(C2[2*J+0]+C2[2*J+1]-C[J]),Threshold);
                        Inc(J);
                    end;
                end;
                
                //
                // test non-weighted fitting
                //
                SetLength(W2, N);
                I:=0;
                while I<=N-1 do
                begin
                    W2[I] := 1;
                    Inc(I);
                end;
                LSFitLinearW(Y, W2, A, N, M, Info, C, Rep);
                LSFitLinear(Y, A, N, M, Info2, C2, Rep2);
                if (Info<=0) or (Info2<=0) then
                begin
                    LLSErrors := True;
                end
                else
                begin
                    
                    //
                    // test answer correctness
                    //
                    J:=0;
                    while J<=M-1 do
                    begin
                        LLSErrors := LLSErrors or AP_FP_Greater(AbsReal(C[J]-C2[J]),Threshold);
                        Inc(J);
                    end;
                    LLSErrors := LLSErrors or AP_FP_Greater(AbsReal(Rep.TaskRCond-Rep2.TaskRCond),Threshold);
                end;
                
                //
                // test nonlinear fitting on the linear task
                // (only non-degenerate task are tested)
                // and compare with answer from linear fitting subroutine
                //
                if N>=M then
                begin
                    SetLength(C2, M);
                    
                    //
                    // test gradient-only or Hessian-based weighted fitting
                    //
                    LSFitLinearW(Y, W, A, N, M, Info, C, Rep);
                    I:=0;
                    while I<=M-1 do
                    begin
                        C2[I] := 2*RandomReal-1;
                        Inc(I);
                    end;
                    LSFitNonlinearWFG(A, Y, W, C2, N, M, M, 0.0, NLThreshold, 0, AP_FP_Greater(RandomReal,0.5), State);
                    FitLinearNonlinear(M, True, A, State, NLSErrors);
                    LSFitNonlinearResults(State, Info, C2, Rep2);
                    if Info<=0 then
                    begin
                        NLSErrors := True;
                    end
                    else
                    begin
                        I:=0;
                        while I<=M-1 do
                        begin
                            NLSErrors := NLSErrors or AP_FP_Greater(AbsReal(C[I]-C2[I]),100*NLThreshold);
                            Inc(I);
                        end;
                    end;
                    I:=0;
                    while I<=M-1 do
                    begin
                        C2[I] := 2*RandomReal-1;
                        Inc(I);
                    end;
                    LSFitNonlinearWFGH(A, Y, W, C2, N, M, M, 0.0, NLThreshold, 0, State);
                    FitLinearNonlinear(M, False, A, State, NLSErrors);
                    LSFitNonlinearResults(State, Info, C2, Rep2);
                    if Info<=0 then
                    begin
                        NLSErrors := True;
                    end
                    else
                    begin
                        I:=0;
                        while I<=M-1 do
                        begin
                            NLSErrors := NLSErrors or AP_FP_Greater(AbsReal(C[I]-C2[I]),100*NLThreshold);
                            Inc(I);
                        end;
                    end;
                    
                    //
                    // test gradient-only or Hessian-based fitting without weights
                    //
                    LSFitLinear(Y, A, N, M, Info, C, Rep);
                    I:=0;
                    while I<=M-1 do
                    begin
                        C2[I] := 2*RandomReal-1;
                        Inc(I);
                    end;
                    LSFitNonlinearFG(A, Y, C2, N, M, M, 0.0, NLThreshold, 0, AP_FP_Greater(RandomReal,0.5), State);
                    FitLinearNonlinear(M, True, A, State, NLSErrors);
                    LSFitNonlinearResults(State, Info, C2, Rep2);
                    if Info<=0 then
                    begin
                        NLSErrors := True;
                    end
                    else
                    begin
                        I:=0;
                        while I<=M-1 do
                        begin
                            NLSErrors := NLSErrors or AP_FP_Greater(AbsReal(C[I]-C2[I]),100*NLThreshold);
                            Inc(I);
                        end;
                    end;
                    I:=0;
                    while I<=M-1 do
                    begin
                        C2[I] := 2*RandomReal-1;
                        Inc(I);
                    end;
                    LSFitNonlinearFGH(A, Y, C2, N, M, M, 0.0, NLThreshold, 0, State);
                    FitLinearNonlinear(M, False, A, State, NLSErrors);
                    LSFitNonlinearResults(State, Info, C2, Rep2);
                    if Info<=0 then
                    begin
                        NLSErrors := True;
                    end
                    else
                    begin
                        I:=0;
                        while I<=M-1 do
                        begin
                            NLSErrors := NLSErrors or AP_FP_Greater(AbsReal(C[I]-C2[I]),100*NLThreshold);
                            Inc(I);
                        end;
                    end;
                end;
                Inc(Pass);
            end;
            Inc(M);
        end;
        
        //
        // test correctness of the RCond field
        //
        SetLength(A, N-1+1, N-1+1);
        SetLength(X, N-1+1);
        SetLength(Y, N-1+1);
        SetLength(W, N-1+1);
        V1 := MaxRealNumber;
        V2 := MinRealNumber;
        I:=0;
        while I<=N-1 do
        begin
            X[I] := 0.1+0.9*RandomReal;
            Y[I] := 0.1+0.9*RandomReal;
            W[I] := 1;
            J:=0;
            while J<=N-1 do
            begin
                if I=J then
                begin
                    A[I,I] := 0.1+0.9*RandomReal;
                    V1 := Min(V1, A[I,I]);
                    V2 := Max(V2, A[I,I]);
                end
                else
                begin
                    A[I,J] := 0;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        LSFitLinearW(Y, W, A, N, N, Info, C, Rep);
        if Info<=0 then
        begin
            LLSErrors := True;
        end
        else
        begin
            LLSErrors := LLSErrors or AP_FP_Greater(AbsReal(Rep.TaskRCond-V1/V2),Threshold);
        end;
        Inc(N);
    end;
    
    //
    // Test constrained least squares
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        N:=1;
        while N<=MaxN do
        begin
            M:=1;
            while M<=MaxM do
            begin
                
                //
                // test for K<>0
                //
                K:=1;
                while K<=M-1 do
                begin
                    
                    //
                    // Prepare Chebyshev basis. Its condition number is very good.
                    // Prepare constraints (random numbers)
                    //
                    SetLength(A, N, M);
                    SetLength(X, N);
                    SetLength(Y, N);
                    SetLength(W, N);
                    XScale := 0.9+0.1*RandomReal;
                    I:=0;
                    while I<=N-1 do
                    begin
                        if N=1 then
                        begin
                            X[I] := 2*RandomReal-1;
                        end
                        else
                        begin
                            X[I] := XScale*(AP_Double(2*I)/(N-1)-1);
                        end;
                        Y[I] := 3*X[I]+Exp(X[I]);
                        W[I] := 1+RandomReal;
                        A[I,0] := 1;
                        if M>1 then
                        begin
                            A[I,1] := X[I];
                        end;
                        J:=2;
                        while J<=M-1 do
                        begin
                            A[I,J] := 2*X[I]*A[I,J-1]-A[I,J-2];
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    SetLength(CM, K, M+1);
                    I:=0;
                    while I<=K-1 do
                    begin
                        J:=0;
                        while J<=M do
                        begin
                            CM[I,J] := 2*RandomReal-1;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    
                    //
                    // Solve constrained task
                    //
                    LSFitLinearWC(Y, W, A, CM, N, M, K, Info, C, Rep);
                    if Info<=0 then
                    begin
                        LLSErrors := True;
                    end
                    else
                    begin
                        LLSErrors := LLSErrors or  not IsGLSSolution(N, M, K, Y, W, A, CM, C);
                    end;
                    
                    //
                    // test non-weighted fitting
                    //
                    SetLength(W2, N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        W2[I] := 1;
                        Inc(I);
                    end;
                    LSFitLinearWC(Y, W2, A, CM, N, M, K, Info, C, Rep);
                    LSFitLinearC(Y, A, CM, N, M, K, Info2, C2, Rep2);
                    if (Info<=0) or (Info2<=0) then
                    begin
                        LLSErrors := True;
                    end
                    else
                    begin
                        
                        //
                        // test answer correctness
                        //
                        J:=0;
                        while J<=M-1 do
                        begin
                            LLSErrors := LLSErrors or AP_FP_Greater(AbsReal(C[J]-C2[J]),Threshold);
                            Inc(J);
                        end;
                        LLSErrors := LLSErrors or AP_FP_Greater(AbsReal(Rep.TaskRCond-Rep2.TaskRCond),Threshold);
                    end;
                    Inc(K);
                end;
                Inc(M);
            end;
            Inc(N);
        end;
        Inc(Pass);
    end;
    
    //
    // nonlinear task for nonlinear fitting:
    //
    //     f(X,C) = 1/(1+C*X^2),
    //     C(true) = 2.
    //
    N := 100;
    SetLength(C, 1);
    C[0] := 1+2*RandomReal;
    SetLength(A, N, 1);
    SetLength(Y, N);
    I:=0;
    while I<=N-1 do
    begin
        A[I,0] := 4*RandomReal-2;
        Y[I] := 1/(1+2*AP_Sqr(A[I,0]));
        Inc(I);
    end;
    LSFitNonlinearFG(A, Y, C, N, 1, 1, 0.0, NLThreshold, 0, True, State);
    while LSFitNonlinearIteration(State) do
    begin
        if State.NeedF then
        begin
            State.F := 1/(1+State.C[0]*AP_Sqr(State.X[0]));
        end;
        if State.NeedFG then
        begin
            State.F := 1/(1+State.C[0]*AP_Sqr(State.X[0]));
            State.G[0] := -AP_Sqr(State.X[0])/AP_Sqr(1+State.C[0]*AP_Sqr(State.X[0]));
        end;
    end;
    LSFitNonlinearResults(State, Info, C, Rep);
    if Info<=0 then
    begin
        NLSErrors := True;
    end
    else
    begin
        NLSErrors := NLSErrors or AP_FP_Greater(AbsReal(C[0]-2),100*NLThreshold);
    end;
    
    //
    // solve simple task (fitting by constant function) and check
    // correctness of the errors calculated by subroutines
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        
        //
        // test on task with non-zero Yi
        //
        N := 4;
        V1 := RandomReal;
        V2 := RandomReal;
        V := 1+RandomReal;
        SetLength(C, 1);
        C[0] := 1+2*RandomReal;
        SetLength(A, 4, 1);
        SetLength(Y, 4);
        A[0,0] := 1;
        Y[0] := V-V2;
        A[1,0] := 1;
        Y[1] := V-V1;
        A[2,0] := 1;
        Y[2] := V+V1;
        A[3,0] := 1;
        Y[3] := V+V2;
        RefRms := Sqrt((AP_Sqr(V1)+AP_Sqr(V2))/2);
        RefAvg := (AbsReal(V1)+AbsReal(V2))/2;
        RefAvgRel := 0.25*(AbsReal(V2)/AbsReal(V-V2)+AbsReal(V1)/AbsReal(V-V1)+AbsReal(V1)/AbsReal(V+V1)+AbsReal(V2)/AbsReal(V+V2));
        RefMax := Max(V1, V2);
        
        //
        // Test LLS
        //
        LSFitLinear(Y, A, 4, 1, Info, C, Rep);
        if Info<=0 then
        begin
            LLSErrors := True;
        end
        else
        begin
            LLSErrors := LLSErrors or AP_FP_Greater(AbsReal(C[0]-V),Threshold);
            LLSErrors := LLSErrors or AP_FP_Greater(AbsReal(Rep.RMSError-RefRMS),Threshold);
            LLSErrors := LLSErrors or AP_FP_Greater(AbsReal(Rep.AvgError-RefAvg),Threshold);
            LLSErrors := LLSErrors or AP_FP_Greater(AbsReal(Rep.AvgRelError-RefAvgRel),Threshold);
            LLSErrors := LLSErrors or AP_FP_Greater(AbsReal(Rep.MaxError-RefMax),Threshold);
        end;
        
        //
        // Test NLS
        //
        LSFitNonlinearFG(A, Y, C, 4, 1, 1, 0.0, NLThreshold, 0, True, State);
        while LSFitNonlinearIteration(State) do
        begin
            if State.NeedF then
            begin
                State.F := State.C[0];
            end;
            if State.NeedFG then
            begin
                State.F := State.C[0];
                State.G[0] := 1;
            end;
        end;
        LSFitNonlinearResults(State, Info, C, Rep);
        if Info<=0 then
        begin
            NLSErrors := True;
        end
        else
        begin
            NLSErrors := NLSErrors or AP_FP_Greater(AbsReal(C[0]-V),Threshold);
            NLSErrors := NLSErrors or AP_FP_Greater(AbsReal(Rep.RMSError-RefRMS),Threshold);
            NLSErrors := NLSErrors or AP_FP_Greater(AbsReal(Rep.AvgError-RefAvg),Threshold);
            NLSErrors := NLSErrors or AP_FP_Greater(AbsReal(Rep.AvgRelError-RefAvgRel),Threshold);
            NLSErrors := NLSErrors or AP_FP_Greater(AbsReal(Rep.MaxError-RefMax),Threshold);
        end;
        Inc(Pass);
    end;
    
    //
    // report
    //
    WasErrors := LLSErrors or NLSErrors;
    if  not Silent then
    begin
        Write(Format('TESTING LEAST SQUARES'#13#10'',[]));
        Write(Format('LINEAR LEAST SQUARES:                    ',[]));
        if LLSErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('NON-LINEAR LEAST SQUARES:                ',[]));
        if NLSErrors then
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
    
    //
    // end
    //
    Result :=  not WasErrors;
end;


(*************************************************************************
Tests whether C is solution of (possibly) constrained LLS problem
*************************************************************************)
function IsGLSSolution(N : AlglibInteger;
     M : AlglibInteger;
     K : AlglibInteger;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     const FMatrix : TReal2DArray;
     const CMatrix : TReal2DArray;
     C : TReal1DArray):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    C2 : TReal1DArray;
    SV : TReal1DArray;
    DeltaC : TReal1DArray;
    DeltaProj : TReal1DArray;
    U : TReal2DArray;
    VT : TReal2DArray;
    V : Double;
    S1 : Double;
    S2 : Double;
    S3 : Double;
    Delta : Double;
    Threshold : Double;
begin
    C := DynamicArrayCopy(C);
    
    //
    // Setup.
    // Threshold is small because CMatrix may be ill-conditioned
    //
    Delta := 0.001;
    Threshold := Sqrt(MachineEpsilon);
    SetLength(C2, M);
    SetLength(DeltaC, M);
    SetLength(DeltaProj, M);
    
    //
    // test whether C is feasible point or not (projC must be close to C)
    //
    I:=0;
    while I<=K-1 do
    begin
        V := APVDotProduct(@CMatrix[I][0], 0, M-1, @C[0], 0, M-1);
        if AP_FP_Greater(AbsReal(V-CMatrix[I,M]),Threshold) then
        begin
            Result := False;
            Exit;
        end;
        Inc(I);
    end;
    
    //
    // find orthogonal basis of Null(CMatrix) (stored in rows from K to M-1)
    //
    if K>0 then
    begin
        RMatrixSVD(CMatrix, K, M, 0, 2, 2, SV, U, VT);
    end;
    
    //
    // Test result
    //
    Result := True;
    S1 := GetGLSError(N, M, Y, W, FMatrix, C);
    J:=0;
    while J<=M-1 do
    begin
        
        //
        // prepare modification of C which leave us in the feasible set.
        //
        // let deltaC be increment on Jth coordinate, then project
        // deltaC in the Null(CMatrix) and store result in DeltaProj
        //
        APVMove(@C2[0], 0, M-1, @C[0], 0, M-1);
        I:=0;
        while I<=M-1 do
        begin
            if I=J then
            begin
                DeltaC[I] := Delta;
            end
            else
            begin
                DeltaC[I] := 0;
            end;
            Inc(I);
        end;
        if K=0 then
        begin
            APVMove(@DeltaProj[0], 0, M-1, @DeltaC[0], 0, M-1);
        end
        else
        begin
            I:=0;
            while I<=M-1 do
            begin
                DeltaProj[I] := 0;
                Inc(I);
            end;
            I:=K;
            while I<=M-1 do
            begin
                V := APVDotProduct(@VT[I][0], 0, M-1, @DeltaC[0], 0, M-1);
                APVAdd(@DeltaProj[0], 0, M-1, @VT[I][0], 0, M-1, V);
                Inc(I);
            end;
        end;
        
        //
        // now we have DeltaProj such that if C is feasible,
        // then C+DeltaProj is feasible too
        //
        APVMove(@C2[0], 0, M-1, @C[0], 0, M-1);
        APVAdd(@C2[0], 0, M-1, @DeltaProj[0], 0, M-1);
        S2 := GetGLSError(N, M, Y, W, FMatrix, C2);
        APVMove(@C2[0], 0, M-1, @C[0], 0, M-1);
        APVSub(@C2[0], 0, M-1, @DeltaProj[0], 0, M-1);
        S3 := GetGLSError(N, M, Y, W, FMatrix, C2);
        Result := Result and AP_FP_Greater_Eq(S2,S1/(1+Threshold)) and AP_FP_Greater_Eq(S3,S1/(1+Threshold));
        Inc(J);
    end;
end;


(*************************************************************************
Tests whether C is solution of LLS problem
*************************************************************************)
function GetGLSError(N : AlglibInteger;
     M : AlglibInteger;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     const FMatrix : TReal2DArray;
     const C : TReal1DArray):Double;
var
    I : AlglibInteger;
    V : Double;
begin
    Result := 0;
    I:=0;
    while I<=N-1 do
    begin
        V := APVDotProduct(@FMatrix[I][0], 0, M-1, @C[0], 0, M-1);
        Result := Result+AP_Sqr(W[I]*(V-Y[I]));
        Inc(I);
    end;
end;


(*************************************************************************
Subroutine for nonlinear fitting of linear problem
*************************************************************************)
procedure FitLinearNonlinear(M : AlglibInteger;
     GradOnly : Boolean;
     const XY : TReal2DArray;
     var State : LSFitState;
     var NLSErrors : Boolean);
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
begin
    while LSFitNonlinearIteration(State) do
    begin
        
        //
        // assume that one and only one of flags is set
        // test that we didn't request hessian in hessian-free setting
        //
        if GradOnly and State.NeedFGH then
        begin
            NLSErrors := True;
        end;
        I := 0;
        if State.NeedF then
        begin
            I := I+1;
        end;
        if State.NeedFG then
        begin
            I := I+1;
        end;
        if State.NeedFGH then
        begin
            I := I+1;
        end;
        if I<>1 then
        begin
            NLSErrors := True;
        end;
        
        //
        // test that PointIndex is consistent with actual point passed
        //
        I:=0;
        while I<=M-1 do
        begin
            NLSErrors := NLSErrors or AP_FP_Neq(XY[State.PointIndex,I],State.X[I]);
            Inc(I);
        end;
        
        //
        // calculate
        //
        if State.NeedF then
        begin
            V := APVDotProduct(@State.X[0], 0, M-1, @State.C[0], 0, M-1);
            State.F := V;
            Continue;
        end;
        if State.NeedFG then
        begin
            V := APVDotProduct(@State.X[0], 0, M-1, @State.C[0], 0, M-1);
            State.F := V;
            APVMove(@State.G[0], 0, M-1, @State.X[0], 0, M-1);
            Continue;
        end;
        if State.NeedFGH then
        begin
            V := APVDotProduct(@State.X[0], 0, M-1, @State.C[0], 0, M-1);
            State.F := V;
            APVMove(@State.G[0], 0, M-1, @State.X[0], 0, M-1);
            I:=0;
            while I<=M-1 do
            begin
                J:=0;
                while J<=M-1 do
                begin
                    State.H[I,J] := 0;
                    Inc(J);
                end;
                Inc(I);
            end;
            Continue;
        end;
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function llstestunit_test_silent():Boolean;
begin
    Result := TestLLS(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function llstestunit_test():Boolean;
begin
    Result := TestLLS(False);
end;


end.