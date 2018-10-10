unit testdensesolverunit;
interface
uses Math, Sysutils, Ap, reflections, bidiagonal, qr, lq, blas, rotations, bdsvd, svd, lu, trlinsolve, rcond, tsort, xblas, densesolver;

function TestDenseSolver(Silent : Boolean):Boolean;
function testdensesolverunit_test_silent():Boolean;
function testdensesolverunit_test():Boolean;

implementation

procedure Unset2D(var X : TReal2DArray);forward;
procedure Unset1D(var X : TReal1DArray);forward;
procedure UnsetRep(var R : DenseSolverReport);forward;
procedure UnsetLSRep(var R : DenseSolverLSReport);forward;


(*************************************************************************
Test
*************************************************************************)
function TestDenseSolver(Silent : Boolean):Boolean;
var
    A : TReal2DArray;
    ATmp : TReal2DArray;
    XE : TReal2DArray;
    B : TReal2DArray;
    BV : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    N : AlglibInteger;
    M : AlglibInteger;
    Pass : AlglibInteger;
    TaskKind : AlglibInteger;
    MX : Double;
    V : Double;
    VErr : Double;
    Info : AlglibInteger;
    Rep : DenseSolverReport;
    RepLS : DenseSolverLSReport;
    X : TReal2DArray;
    XV : TReal1DArray;
    Y : TReal1DArray;
    TX : TReal1DArray;
    MaxN : AlglibInteger;
    MaxM : AlglibInteger;
    PassCount : AlglibInteger;
    Threshold : Double;
    WeakThreshold : Double;
    GenRealErrors : Boolean;
    LSRealErrors : Boolean;
    RfsErrors : Boolean;
    WasErrors : Boolean;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    MaxN := 10;
    MaxM := 5;
    PassCount := 5;
    Threshold := 10000*MachineEpsilon;
    WeakThreshold := 0.001;
    RfsErrors := False;
    GenRealErrors := False;
    LSRealErrors := False;
    WasErrors := False;
    
    //
    // General square matrices:
    // * test general properties
    // * test iterative improvement
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
                // ********************************************************
                // WELL CONDITIONED TASKS
                // ability to find correct solution is tested
                // ********************************************************
                //
                // 1. generate random well conditioned matrix A:
                //        A = R + W,
                //    where R has small random elements, and W is
                //    permutation of diagonal matrix with large
                //    random entries. So A should be well conditioned.
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods on original A
                //
                SetLength(A, N, N);
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        A[I,J] := 2*RandomReal-1;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                K := RandomInteger(N);
                I:=0;
                while I<=N-1 do
                begin
                    A[I,(I+K) mod N] := (RandomInteger(2)-0.5)*(5+5*RandomReal);
                    Inc(I);
                end;
                SetLength(XE, N, M);
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=M-1 do
                    begin
                        XE[I,J] := 2*RandomReal-1;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                SetLength(B, N, M);
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=M-1 do
                    begin
                        V := 0.0;
                        for i_ := 0 to N-1 do
                        begin
                            V := V + A[I,i_]*XE[i_,J];
                        end;
                        B[I,J] := V;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                
                //
                // Test RMatrixSolveM()
                //
                Info := 0;
                UnsetRep(Rep);
                Unset2D(X);
                RMatrixSolveM(A, N, B, M, Info, Rep, X);
                if Info<=0 then
                begin
                    GenRealErrors := True;
                end
                else
                begin
                    GenRealErrors := GenRealErrors or AP_FP_Less(Rep.R1,100*MachineEpsilon) or AP_FP_Greater(Rep.R1,1+1000*MachineEpsilon);
                    GenRealErrors := GenRealErrors or AP_FP_Less(Rep.RInf,100*MachineEpsilon) or AP_FP_Greater(Rep.RInf,1+1000*MachineEpsilon);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=M-1 do
                        begin
                            GenRealErrors := GenRealErrors or AP_FP_Greater(AbsReal(XE[I,J]-X[I,J]),Threshold);
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                end;
                
                //
                // Test RMatrixSolve()
                //
                Info := 0;
                UnsetRep(Rep);
                Unset2D(X);
                SetLength(BV, N);
                for i_ := 0 to N-1 do
                begin
                    BV[i_] := B[i_,0];
                end;
                RMatrixSolve(A, N, BV, Info, Rep, XV);
                if Info<=0 then
                begin
                    GenRealErrors := True;
                end
                else
                begin
                    GenRealErrors := GenRealErrors or AP_FP_Less(Rep.R1,100*MachineEpsilon) or AP_FP_Greater(Rep.R1,1+1000*MachineEpsilon);
                    GenRealErrors := GenRealErrors or AP_FP_Less(Rep.RInf,100*MachineEpsilon) or AP_FP_Greater(Rep.RInf,1+1000*MachineEpsilon);
                    I:=0;
                    while I<=N-1 do
                    begin
                        GenRealErrors := GenRealErrors or AP_FP_Greater(AbsReal(XE[I,0]-XV[I]),Threshold);
                        Inc(I);
                    end;
                end;
                
                //
                // Test DenseSolverRLS():
                // * test on original system A*x = b
                // * test on overdetermined system with the same solution: (A' A')'*x = b
                // * test on underdetermined system with the same solution: (A 0 0 0 ) * z = b
                //
                Info := 0;
                UnsetLSRep(RepLS);
                Unset1D(XV);
                SetLength(BV, N);
                for i_ := 0 to N-1 do
                begin
                    BV[i_] := B[i_,0];
                end;
                RMatrixSolveLS(A, N, N, BV, 0.0, Info, RepLS, XV);
                if Info<=0 then
                begin
                    LSRealErrors := True;
                end
                else
                begin
                    LSRealErrors := LSRealErrors or AP_FP_Less(RepLS.R2,100*MachineEpsilon) or AP_FP_Greater(RepLS.R2,1+1000*MachineEpsilon);
                    LSRealErrors := LSRealErrors or (RepLS.N<>N) or (RepLS.K<>0);
                    I:=0;
                    while I<=N-1 do
                    begin
                        LSRealErrors := LSRealErrors or AP_FP_Greater(AbsReal(XE[I,0]-XV[I]),Threshold);
                        Inc(I);
                    end;
                end;
                Info := 0;
                UnsetLSRep(RepLS);
                Unset1D(XV);
                SetLength(BV, 2*N);
                for i_ := 0 to N-1 do
                begin
                    BV[i_] := B[i_,0];
                end;
                i1_ := (0) - (N);
                for i_ := N to 2*N-1 do
                begin
                    BV[i_] := B[i_+i1_,0];
                end;
                SetLength(ATmp, 2*N, N);
                CopyMatrix(A, 0, N-1, 0, N-1, ATmp, 0, N-1, 0, N-1);
                CopyMatrix(A, 0, N-1, 0, N-1, ATmp, N, 2*N-1, 0, N-1);
                RMatrixSolveLS(ATmp, 2*N, N, BV, 0.0, Info, RepLS, XV);
                if Info<=0 then
                begin
                    LSRealErrors := True;
                end
                else
                begin
                    LSRealErrors := LSRealErrors or AP_FP_Less(RepLS.R2,100*MachineEpsilon) or AP_FP_Greater(RepLS.R2,1+1000*MachineEpsilon);
                    LSRealErrors := LSRealErrors or (RepLS.N<>N) or (RepLS.K<>0);
                    I:=0;
                    while I<=N-1 do
                    begin
                        LSRealErrors := LSRealErrors or AP_FP_Greater(AbsReal(XE[I,0]-XV[I]),Threshold);
                        Inc(I);
                    end;
                end;
                Info := 0;
                UnsetLSRep(RepLS);
                Unset1D(XV);
                SetLength(BV, N);
                for i_ := 0 to N-1 do
                begin
                    BV[i_] := B[i_,0];
                end;
                SetLength(ATmp, N, 2*N);
                CopyMatrix(A, 0, N-1, 0, N-1, ATmp, 0, N-1, 0, N-1);
                I:=0;
                while I<=N-1 do
                begin
                    J:=N;
                    while J<=2*N-1 do
                    begin
                        ATmp[I,J] := 0;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                RMatrixSolveLS(ATmp, N, 2*N, BV, 0.0, Info, RepLS, XV);
                if Info<=0 then
                begin
                    LSRealErrors := True;
                end
                else
                begin
                    LSRealErrors := LSRealErrors or AP_FP_Neq(RepLS.R2,0);
                    LSRealErrors := LSRealErrors or (RepLS.N<>2*N) or (RepLS.K<>N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        LSRealErrors := LSRealErrors or AP_FP_Greater(AbsReal(XE[I,0]-XV[I]),Threshold);
                        Inc(I);
                    end;
                    I:=N;
                    while I<=2*N-1 do
                    begin
                        LSRealErrors := LSRealErrors or AP_FP_Greater(AbsReal(XV[I]),Threshold);
                        Inc(I);
                    end;
                end;
                
                //
                // ********************************************************
                // EXACTLY SINGULAR MATRICES
                // ability to detect singularity is tested
                // ********************************************************
                //
                // 1. generate different types of singular matrices:
                //    * zero
                //    * with zero columns
                //    * with zero rows
                //    * with equal rows/columns
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods
                //
                TaskKind:=0;
                while TaskKind<=4 do
                begin
                    Unset2D(A);
                    if TaskKind=0 then
                    begin
                        
                        //
                        // all zeros
                        //
                        SetLength(A, N, N);
                        I:=0;
                        while I<=N-1 do
                        begin
                            J:=0;
                            while J<=N-1 do
                            begin
                                A[I,J] := 0;
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                    end;
                    if TaskKind=1 then
                    begin
                        
                        //
                        // there is zero column
                        //
                        SetLength(A, N, N);
                        I:=0;
                        while I<=N-1 do
                        begin
                            J:=0;
                            while J<=N-1 do
                            begin
                                A[I,J] := 2*RandomReal-1;
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        K := RandomInteger(N);
                        for i_ := 0 to N-1 do
                        begin
                            A[i_,K] := 0*A[i_,K];
                        end;
                    end;
                    if TaskKind=2 then
                    begin
                        
                        //
                        // there is zero row
                        //
                        SetLength(A, N, N);
                        I:=0;
                        while I<=N-1 do
                        begin
                            J:=0;
                            while J<=N-1 do
                            begin
                                A[I,J] := 2*RandomReal-1;
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        K := RandomInteger(N);
                        APVMul(@A[K][0], 0, N-1, 0);
                    end;
                    if TaskKind=3 then
                    begin
                        
                        //
                        // equal columns
                        //
                        if N<2 then
                        begin
                            Inc(TaskKind);
                            Continue;
                        end;
                        SetLength(A, N, N);
                        I:=0;
                        while I<=N-1 do
                        begin
                            J:=0;
                            while J<=N-1 do
                            begin
                                A[I,J] := 2*RandomReal-1;
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        K := 1+RandomInteger(N-1);
                        for i_ := 0 to N-1 do
                        begin
                            A[i_,0] := A[i_,K];
                        end;
                    end;
                    if TaskKind=4 then
                    begin
                        
                        //
                        // equal rows
                        //
                        if N<2 then
                        begin
                            Inc(TaskKind);
                            Continue;
                        end;
                        SetLength(A, N, N);
                        I:=0;
                        while I<=N-1 do
                        begin
                            J:=0;
                            while J<=N-1 do
                            begin
                                A[I,J] := 2*RandomReal-1;
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        K := 1+RandomInteger(N-1);
                        APVMove(@A[0][0], 0, N-1, @A[K][0], 0, N-1);
                    end;
                    SetLength(XE, N, M);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=M-1 do
                        begin
                            XE[I,J] := 2*RandomReal-1;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    SetLength(B, N, M);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=M-1 do
                        begin
                            V := 0.0;
                            for i_ := 0 to N-1 do
                            begin
                                V := V + A[I,i_]*XE[i_,J];
                            end;
                            B[I,J] := V;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    
                    //
                    // Test RMatrixSolveM()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    Unset2D(X);
                    RMatrixSolveM(A, N, B, M, Info, Rep, X);
                    GenRealErrors := GenRealErrors or (Info<>-3);
                    GenRealErrors := GenRealErrors or AP_FP_Less(Rep.R1,0) or AP_FP_Greater(Rep.R1,1000*MachineEpsilon);
                    GenRealErrors := GenRealErrors or AP_FP_Less(Rep.RInf,0) or AP_FP_Greater(Rep.RInf,1000*MachineEpsilon);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=M-1 do
                        begin
                            GenRealErrors := GenRealErrors or AP_FP_Neq(X[I,J],0);
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    
                    //
                    // Test RMatrixSolve()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    Unset2D(X);
                    SetLength(BV, N);
                    for i_ := 0 to N-1 do
                    begin
                        BV[i_] := B[i_,0];
                    end;
                    RMatrixSolve(A, N, BV, Info, Rep, XV);
                    GenRealErrors := GenRealErrors or (Info<>-3);
                    GenRealErrors := GenRealErrors or AP_FP_Less(Rep.R1,0) or AP_FP_Greater(Rep.R1,1000*MachineEpsilon);
                    GenRealErrors := GenRealErrors or AP_FP_Less(Rep.RInf,0) or AP_FP_Greater(Rep.RInf,1000*MachineEpsilon);
                    I:=0;
                    while I<=N-1 do
                    begin
                        GenRealErrors := GenRealErrors or AP_FP_Neq(XV[I],0);
                        Inc(I);
                    end;
                    Inc(TaskKind);
                end;
                Inc(M);
            end;
            Inc(N);
        end;
        Inc(Pass);
    end;
    Pass:=1;
    while Pass<=PassCount do
    begin
        
        //
        // Test iterative improvement for real matrices
        //
        // A matrix/right part are constructed such that both matrix
        // and solution components are within (0,1). Such matrix/right part
        // have nice properties - system can be solved using iterative
        // improvement with |A*x-b| about several ulps of max(1,|b|).
        //
        N := 100;
        SetLength(A, N, N);
        SetLength(B, N, 1);
        SetLength(BV, N);
        SetLength(TX, N);
        SetLength(XV, N);
        SetLength(Y, N);
        I:=0;
        while I<=N-1 do
        begin
            XV[I] := RandomReal;
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                A[I,J] := RandomReal;
                Inc(J);
            end;
            APVMove(@Y[0], 0, N-1, @A[I][0], 0, N-1);
            XDot(Y, XV, N, TX, V, VErr);
            BV[I] := V;
            Inc(I);
        end;
        for i_ := 0 to N-1 do
        begin
            B[i_,0] := BV[i_];
        end;
        
        //
        // Test RMatrixSolveM()
        //
        Unset2D(X);
        RMatrixSolveM(A, N, B, 1, Info, Rep, X);
        if Info<=0 then
        begin
            RfsErrors := True;
        end
        else
        begin
            SetLength(XV, N);
            for i_ := 0 to N-1 do
            begin
                XV[i_] := X[i_,0];
            end;
            I:=0;
            while I<=N-1 do
            begin
                APVMove(@Y[0], 0, N-1, @A[I][0], 0, N-1);
                XDot(Y, XV, N, TX, V, VErr);
                RfsErrors := RfsErrors or AP_FP_Greater(AbsReal(V-B[I,0]),4*MachineEpsilon*Max(1, AbsReal(B[I,0])));
                Inc(I);
            end;
        end;
        
        //
        // Test RMatrixSolve()
        //
        Unset1D(XV);
        RMatrixSolve(A, N, BV, Info, Rep, XV);
        if Info<=0 then
        begin
            RfsErrors := True;
        end
        else
        begin
            I:=0;
            while I<=N-1 do
            begin
                APVMove(@Y[0], 0, N-1, @A[I][0], 0, N-1);
                XDot(Y, XV, N, TX, V, VErr);
                RfsErrors := RfsErrors or AP_FP_Greater(AbsReal(V-BV[I]),4*MachineEpsilon*Max(1, AbsReal(BV[I])));
                Inc(I);
            end;
        end;
        
        //
        // Test LS-solver on the same matrix
        //
        RMatrixSolveLS(A, N, N, BV, 0.0, Info, RepLS, XV);
        if Info<=0 then
        begin
            RfsErrors := True;
        end
        else
        begin
            I:=0;
            while I<=N-1 do
            begin
                APVMove(@Y[0], 0, N-1, @A[I][0], 0, N-1);
                XDot(Y, XV, N, TX, V, VErr);
                RfsErrors := RfsErrors or AP_FP_Greater(AbsReal(V-BV[I]),4*MachineEpsilon*Max(1, AbsReal(BV[I])));
                Inc(I);
            end;
        end;
        Inc(Pass);
    end;
    
    //
    // end
    //
    WasErrors := GenRealErrors or LSRealErrors or RfsErrors;
    if  not Silent then
    begin
        Write(Format('TESTING DENSE SOLVER'#13#10'',[]));
        Write(Format('* GENERAL REAL:                           ',[]));
        if GenRealErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* LEAST SQUARES REAL:                     ',[]));
        if LSRealErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* ITERATIVE IMPROVEMENT:                  ',[]));
        if RfsErrors then
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
    end;
    Result :=  not WasErrors;
end;


(*************************************************************************
Unsets real matrix
*************************************************************************)
procedure Unset2D(var X : TReal2DArray);
begin
    SetLength(X, 1, 1);
    X[0,0] := 2*RandomReal-1;
end;


(*************************************************************************
Unsets real vector
*************************************************************************)
procedure Unset1D(var X : TReal1DArray);
begin
    SetLength(X, 1);
    X[0] := 2*RandomReal-1;
end;


(*************************************************************************
Unsets report
*************************************************************************)
procedure UnsetRep(var R : DenseSolverReport);
begin
    R.R1 := -1;
    R.RInf := -1;
end;


(*************************************************************************
Unsets report
*************************************************************************)
procedure UnsetLSRep(var R : DenseSolverLSReport);
begin
    R.R2 := -1;
    R.N := -1;
    R.K := -1;
    Unset2D(R.CX);
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testdensesolverunit_test_silent():Boolean;
begin
    Result := TestDenseSolver(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testdensesolverunit_test():Boolean;
begin
    Result := TestDenseSolver(False);
end;


end.