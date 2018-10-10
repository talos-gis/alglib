unit testdensesolverunit;
interface
uses Math, Sysutils, Ap, reflections, bidiagonal, qr, lq, blas, rotations, bdsvd, svd, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, tsort, xblas, densesolver;

function TestDenseSolver(Silent : Boolean):Boolean;
function testdensesolverunit_test_silent():Boolean;
function testdensesolverunit_test():Boolean;

implementation

function RMatrixCheckSolutionM(const XE : TReal2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TReal2DArray):Boolean;forward;
function RMatrixCheckSolution(const XE : TReal2DArray;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TReal1DArray):Boolean;forward;
function RMatrixCheckSingularM(N : AlglibInteger;
     M : AlglibInteger;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TReal2DArray):Boolean;forward;
function RMatrixCheckSingular(N : AlglibInteger;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TReal1DArray):Boolean;forward;
function CMatrixCheckSolutionM(const XE : TComplex2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TComplex2DArray):Boolean;forward;
function CMatrixCheckSolution(const XE : TComplex2DArray;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TComplex1DArray):Boolean;forward;
function CMatrixCheckSingularM(N : AlglibInteger;
     M : AlglibInteger;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TComplex2DArray):Boolean;forward;
function CMatrixCheckSingular(N : AlglibInteger;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TComplex1DArray):Boolean;forward;
procedure RMatrixMakeACopy(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);forward;
procedure CMatrixMakeACopy(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TComplex2DArray);forward;
procedure RMatrixDropHalf(var A : TReal2DArray;
     N : AlglibInteger;
     DropLower : Boolean);forward;
procedure CMatrixDropHalf(var A : TComplex2DArray;
     N : AlglibInteger;
     DropLower : Boolean);forward;
procedure TestRSolver(MaxN : AlglibInteger;
     MaxM : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var RErrors : Boolean;
     var RfsErrors : Boolean);forward;
procedure TestSPDSolver(MaxN : AlglibInteger;
     MaxM : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var SPDErrors : Boolean;
     var RfsErrors : Boolean);forward;
procedure TestCSolver(MaxN : AlglibInteger;
     MaxM : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var CErrors : Boolean;
     var RfsErrors : Boolean);forward;
procedure TestHPDSolver(MaxN : AlglibInteger;
     MaxM : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var HPDErrors : Boolean;
     var RfsErrors : Boolean);forward;
procedure Unset2D(var X : TReal2DArray);forward;
procedure Unset1D(var X : TReal1DArray);forward;
procedure CUnset2D(var X : TComplex2DArray);forward;
procedure CUnset1D(var X : TComplex1DArray);forward;
procedure UnsetRep(var R : DenseSolverReport);forward;
procedure UnsetLSRep(var R : DenseSolverLSReport);forward;


(*************************************************************************
Test
*************************************************************************)
function TestDenseSolver(Silent : Boolean):Boolean;
var
    A : TReal2DArray;
    LUA : TReal2DArray;
    ATmp : TReal2DArray;
    P : TInteger1DArray;
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
    RErrors : Boolean;
    CErrors : Boolean;
    SPDErrors : Boolean;
    HPDErrors : Boolean;
    RfsErrors : Boolean;
    WasErrors : Boolean;
begin
    MaxN := 10;
    MaxM := 5;
    PassCount := 5;
    Threshold := 10000*MachineEpsilon;
    RfsErrors := False;
    RErrors := False;
    CErrors := False;
    SPDErrors := False;
    HPDErrors := False;
    TestRSolver(MaxN, MaxM, PassCount, Threshold, RErrors, RfsErrors);
    TestSPDSolver(MaxN, MaxM, PassCount, Threshold, SPDErrors, RfsErrors);
    TestCSolver(MaxN, MaxM, PassCount, Threshold, CErrors, RfsErrors);
    TestHPDSolver(MaxN, MaxM, PassCount, Threshold, HPDErrors, RfsErrors);
    WasErrors := RErrors or CErrors or SPDErrors or HPDErrors or RfsErrors;
    if  not Silent then
    begin
        Write(Format('TESTING DENSE SOLVER'#13#10'',[]));
        Write(Format('* REAL:                                   ',[]));
        if RErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* COMPLEX:                                ',[]));
        if CErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* SPD:                                    ',[]));
        if SPDErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* HPD:                                    ',[]));
        if HPDErrors then
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
Checks whether solver results are correct solution.
Returns True on success.
*************************************************************************)
function RMatrixCheckSolutionM(const XE : TReal2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TReal2DArray):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Result := True;
    if Info<=0 then
    begin
        Result := False;
    end
    else
    begin
        Result := Result and  not (AP_FP_Less(Rep.R1,100*MachineEpsilon) or AP_FP_Greater(Rep.R1,1+1000*MachineEpsilon));
        Result := Result and  not (AP_FP_Less(Rep.RInf,100*MachineEpsilon) or AP_FP_Greater(Rep.RInf,1+1000*MachineEpsilon));
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                Result := Result and AP_FP_Less_Eq(AbsReal(XE[I,J]-XS[I,J]),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Checks whether solver results are correct solution.
Returns True on success.
*************************************************************************)
function RMatrixCheckSolution(const XE : TReal2DArray;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TReal1DArray):Boolean;
var
    XSM : TReal2DArray;
    i_ : AlglibInteger;
begin
    SetLength(XSM, N, 1);
    for i_ := 0 to N-1 do
    begin
        XSM[i_,0] := XS[i_];
    end;
    Result := RMatrixCheckSolutionM(XE, N, 1, Threshold, Info, Rep, XSM);
end;


(*************************************************************************
Checks whether solver results indicate singular matrix.
Returns True on success.
*************************************************************************)
function RMatrixCheckSingularM(N : AlglibInteger;
     M : AlglibInteger;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TReal2DArray):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Result := True;
    if (Info<>-3) and (Info<>1) then
    begin
        Result := False;
    end
    else
    begin
        Result := Result and  not (AP_FP_Less(Rep.R1,0) or AP_FP_Greater(Rep.R1,1000*MachineEpsilon));
        Result := Result and  not (AP_FP_Less(Rep.RInf,0) or AP_FP_Greater(Rep.RInf,1000*MachineEpsilon));
        if Info=-3 then
        begin
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=M-1 do
                begin
                    Result := Result and AP_FP_Eq(XS[I,J],0);
                    Inc(J);
                end;
                Inc(I);
            end;
        end;
    end;
end;


(*************************************************************************
Checks whether solver results indicate singular matrix.
Returns True on success.
*************************************************************************)
function RMatrixCheckSingular(N : AlglibInteger;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TReal1DArray):Boolean;
var
    XSM : TReal2DArray;
    i_ : AlglibInteger;
begin
    SetLength(XSM, N, 1);
    for i_ := 0 to N-1 do
    begin
        XSM[i_,0] := XS[i_];
    end;
    Result := RMatrixCheckSingularM(N, 1, Info, Rep, XSM);
end;


(*************************************************************************
Checks whether solver results are correct solution.
Returns True on success.
*************************************************************************)
function CMatrixCheckSolutionM(const XE : TComplex2DArray;
     N : AlglibInteger;
     M : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TComplex2DArray):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Result := True;
    if Info<=0 then
    begin
        Result := False;
    end
    else
    begin
        Result := Result and  not (AP_FP_Less(Rep.R1,100*MachineEpsilon) or AP_FP_Greater(Rep.R1,1+1000*MachineEpsilon));
        Result := Result and  not (AP_FP_Less(Rep.RInf,100*MachineEpsilon) or AP_FP_Greater(Rep.RInf,1+1000*MachineEpsilon));
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                Result := Result and AP_FP_Less_Eq(AbsComplex(C_Sub(XE[I,J],XS[I,J])),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Checks whether solver results are correct solution.
Returns True on success.
*************************************************************************)
function CMatrixCheckSolution(const XE : TComplex2DArray;
     N : AlglibInteger;
     Threshold : Double;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TComplex1DArray):Boolean;
var
    XSM : TComplex2DArray;
    i_ : AlglibInteger;
begin
    SetLength(XSM, N, 1);
    for i_ := 0 to N-1 do
    begin
        XSM[i_,0] := XS[i_];
    end;
    Result := CMatrixCheckSolutionM(XE, N, 1, Threshold, Info, Rep, XSM);
end;


(*************************************************************************
Checks whether solver results indicate singular matrix.
Returns True on success.
*************************************************************************)
function CMatrixCheckSingularM(N : AlglibInteger;
     M : AlglibInteger;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TComplex2DArray):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Result := True;
    if (Info<>-3) and (Info<>1) then
    begin
        Result := False;
    end
    else
    begin
        Result := Result and  not (AP_FP_Less(Rep.R1,0) or AP_FP_Greater(Rep.R1,1000*MachineEpsilon));
        Result := Result and  not (AP_FP_Less(Rep.RInf,0) or AP_FP_Greater(Rep.RInf,1000*MachineEpsilon));
        if Info=-3 then
        begin
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=M-1 do
                begin
                    Result := Result and C_EqualR(XS[I,J],0);
                    Inc(J);
                end;
                Inc(I);
            end;
        end;
    end;
end;


(*************************************************************************
Checks whether solver results indicate singular matrix.
Returns True on success.
*************************************************************************)
function CMatrixCheckSingular(N : AlglibInteger;
     Info : AlglibInteger;
     const Rep : DenseSolverReport;
     const XS : TComplex1DArray):Boolean;
var
    XSM : TComplex2DArray;
    i_ : AlglibInteger;
begin
    SetLength(XSM, N, 1);
    for i_ := 0 to N-1 do
    begin
        XSM[i_,0] := XS[i_];
    end;
    Result := CMatrixCheckSingularM(N, 1, Info, Rep, XSM);
end;


(*************************************************************************
Copy
*************************************************************************)
procedure RMatrixMakeACopy(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    SetLength(B, M-1+1, N-1+1);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            B[I,J] := A[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Copy
*************************************************************************)
procedure CMatrixMakeACopy(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TComplex2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    SetLength(B, M-1+1, N-1+1);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            B[I,J] := A[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Drops upper or lower half of the matrix - fills it by special pattern
which may be used later to ensure that this part wasn't changed
*************************************************************************)
procedure RMatrixDropHalf(var A : TReal2DArray;
     N : AlglibInteger;
     DropLower : Boolean);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if DropLower and (I>J) or  not DropLower and (I<J) then
            begin
                A[I,J] := 1+2*I+3*J;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Drops upper or lower half of the matrix - fills it by special pattern
which may be used later to ensure that this part wasn't changed
*************************************************************************)
procedure CMatrixDropHalf(var A : TComplex2DArray;
     N : AlglibInteger;
     DropLower : Boolean);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if DropLower and (I>J) or  not DropLower and (I<J) then
            begin
                A[I,J] := C_Complex(1+2*I+3*J);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Real test
*************************************************************************)
procedure TestRSolver(MaxN : AlglibInteger;
     MaxM : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var RErrors : Boolean;
     var RfsErrors : Boolean);
var
    A : TReal2DArray;
    LUA : TReal2DArray;
    ATmp : TReal2DArray;
    P : TInteger1DArray;
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
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
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
                // 1. generate random well conditioned matrix A.
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods on original A
                //
                RMatrixRndCond(N, 1000, A);
                RMatrixMakeACopy(A, N, N, LUA);
                RMatrixLU(LUA, N, N, P);
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
                // Test solvers
                //
                Info := 0;
                UnsetRep(Rep);
                Unset2D(X);
                RMatrixSolveM(A, N, B, M, AP_FP_Greater(RandomReal,0.5), Info, Rep, X);
                RErrors := RErrors or  not RMatrixCheckSolutionM(XE, N, M, Threshold, Info, Rep, X);
                Info := 0;
                UnsetRep(Rep);
                Unset1D(XV);
                SetLength(BV, N);
                for i_ := 0 to N-1 do
                begin
                    BV[i_] := B[i_,0];
                end;
                RMatrixSolve(A, N, BV, Info, Rep, XV);
                RErrors := RErrors or  not RMatrixCheckSolution(XE, N, Threshold, Info, Rep, XV);
                Info := 0;
                UnsetRep(Rep);
                Unset2D(X);
                RMatrixLUSolveM(LUA, P, N, B, M, Info, Rep, X);
                RErrors := RErrors or  not RMatrixCheckSolutionM(XE, N, M, Threshold, Info, Rep, X);
                Info := 0;
                UnsetRep(Rep);
                Unset1D(XV);
                SetLength(BV, N);
                for i_ := 0 to N-1 do
                begin
                    BV[i_] := B[i_,0];
                end;
                RMatrixLUSolve(LUA, P, N, BV, Info, Rep, XV);
                RErrors := RErrors or  not RMatrixCheckSolution(XE, N, Threshold, Info, Rep, XV);
                Info := 0;
                UnsetRep(Rep);
                Unset2D(X);
                RMatrixMixedSolveM(A, LUA, P, N, B, M, Info, Rep, X);
                RErrors := RErrors or  not RMatrixCheckSolutionM(XE, N, M, Threshold, Info, Rep, X);
                Info := 0;
                UnsetRep(Rep);
                Unset1D(XV);
                SetLength(BV, N);
                for i_ := 0 to N-1 do
                begin
                    BV[i_] := B[i_,0];
                end;
                RMatrixMixedSolve(A, LUA, P, N, BV, Info, Rep, XV);
                RErrors := RErrors or  not RMatrixCheckSolution(XE, N, Threshold, Info, Rep, XV);
                
                //
                // Test DenseSolverRLS():
                // * test on original system A*x = b
                // * test on overdetermined system with the same solution: (A' A')'*x = (b' b')'
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
                    RErrors := True;
                end
                else
                begin
                    RErrors := RErrors or AP_FP_Less(RepLS.R2,100*MachineEpsilon) or AP_FP_Greater(RepLS.R2,1+1000*MachineEpsilon);
                    RErrors := RErrors or (RepLS.N<>N) or (RepLS.K<>0);
                    I:=0;
                    while I<=N-1 do
                    begin
                        RErrors := RErrors or AP_FP_Greater(AbsReal(XE[I,0]-XV[I]),Threshold);
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
                    RErrors := True;
                end
                else
                begin
                    RErrors := RErrors or AP_FP_Less(RepLS.R2,100*MachineEpsilon) or AP_FP_Greater(RepLS.R2,1+1000*MachineEpsilon);
                    RErrors := RErrors or (RepLS.N<>N) or (RepLS.K<>0);
                    I:=0;
                    while I<=N-1 do
                    begin
                        RErrors := RErrors or AP_FP_Greater(AbsReal(XE[I,0]-XV[I]),Threshold);
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
                    RErrors := True;
                end
                else
                begin
                    RErrors := RErrors or AP_FP_Neq(RepLS.R2,0);
                    RErrors := RErrors or (RepLS.N<>2*N) or (RepLS.K<>N);
                    I:=0;
                    while I<=N-1 do
                    begin
                        RErrors := RErrors or AP_FP_Greater(AbsReal(XE[I,0]-XV[I]),Threshold);
                        Inc(I);
                    end;
                    I:=N;
                    while I<=2*N-1 do
                    begin
                        RErrors := RErrors or AP_FP_Greater(AbsReal(XV[I]),Threshold);
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
                    RMatrixMakeACopy(A, N, N, LUA);
                    RMatrixLU(LUA, N, N, P);
                    
                    //
                    // Test RMatrixSolveM()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    Unset2D(X);
                    RMatrixSolveM(A, N, B, M, AP_FP_Greater(RandomReal,0.5), Info, Rep, X);
                    RErrors := RErrors or  not RMatrixCheckSingularM(N, M, Info, Rep, X);
                    
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
                    RErrors := RErrors or  not RMatrixCheckSingular(N, Info, Rep, XV);
                    
                    //
                    // Test RMatrixLUSolveM()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    Unset2D(X);
                    RMatrixLUSolveM(LUA, P, N, B, M, Info, Rep, X);
                    RErrors := RErrors or  not RMatrixCheckSingularM(N, M, Info, Rep, X);
                    
                    //
                    // Test RMatrixLUSolve()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    Unset2D(X);
                    SetLength(BV, N);
                    for i_ := 0 to N-1 do
                    begin
                        BV[i_] := B[i_,0];
                    end;
                    RMatrixLUSolve(LUA, P, N, BV, Info, Rep, XV);
                    RErrors := RErrors or  not RMatrixCheckSingular(N, Info, Rep, XV);
                    
                    //
                    // Test RMatrixMixedSolveM()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    Unset2D(X);
                    RMatrixMixedSolveM(A, LUA, P, N, B, M, Info, Rep, X);
                    RErrors := RErrors or  not RMatrixCheckSingularM(N, M, Info, Rep, X);
                    
                    //
                    // Test RMatrixMixedSolve()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    Unset2D(X);
                    SetLength(BV, N);
                    for i_ := 0 to N-1 do
                    begin
                        BV[i_] := B[i_,0];
                    end;
                    RMatrixMixedSolve(A, LUA, P, N, BV, Info, Rep, XV);
                    RErrors := RErrors or  not RMatrixCheckSingular(N, Info, Rep, XV);
                    Inc(TaskKind);
                end;
                Inc(M);
            end;
            Inc(N);
        end;
        Inc(Pass);
    end;
    
    //
    // test iterative improvement
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        
        //
        // Test iterative improvement matrices
        //
        // A matrix/right part are constructed such that both matrix
        // and solution components are within (-1,+1). Such matrix/right part
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
            XV[I] := 2*RandomReal-1;
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
        RMatrixSolveM(A, N, B, 1, True, Info, Rep, X);
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
                RfsErrors := RfsErrors or AP_FP_Greater(AbsReal(V-B[I,0]),8*MachineEpsilon*Max(1, AbsReal(B[I,0])));
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
                RfsErrors := RfsErrors or AP_FP_Greater(AbsReal(V-BV[I]),8*MachineEpsilon*Max(1, AbsReal(BV[I])));
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
                RfsErrors := RfsErrors or AP_FP_Greater(AbsReal(V-BV[I]),8*MachineEpsilon*Max(1, AbsReal(BV[I])));
                Inc(I);
            end;
        end;
        Inc(Pass);
    end;
end;


(*************************************************************************
SPD test
*************************************************************************)
procedure TestSPDSolver(MaxN : AlglibInteger;
     MaxM : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var SPDErrors : Boolean;
     var RfsErrors : Boolean);
var
    A : TReal2DArray;
    CHA : TReal2DArray;
    ATmp : TReal2DArray;
    P : TInteger1DArray;
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
    IsUpper : Boolean;
    Info : AlglibInteger;
    Rep : DenseSolverReport;
    RepLS : DenseSolverLSReport;
    X : TReal2DArray;
    XV : TReal1DArray;
    Y : TReal1DArray;
    TX : TReal1DArray;
    i_ : AlglibInteger;
begin
    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
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
                // 1. generate random well conditioned matrix A.
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods on original A
                //
                IsUpper := AP_FP_Greater(RandomReal,0.5);
                SPDMatrixRndCond(N, 1000, A);
                RMatrixMakeACopy(A, N, N, CHA);
                if  not SPDMatrixCholesky(CHA, N, IsUpper) then
                begin
                    SPDErrors := True;
                    Exit;
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
                RMatrixDropHalf(A, N, IsUpper);
                RMatrixDropHalf(CHA, N, IsUpper);
                
                //
                // Test solvers
                //
                Info := 0;
                UnsetRep(Rep);
                Unset2D(X);
                SPDMatrixSolveM(A, N, IsUpper, B, M, Info, Rep, X);
                SPDErrors := SPDErrors or  not RMatrixCheckSolutionM(XE, N, M, Threshold, Info, Rep, X);
                Info := 0;
                UnsetRep(Rep);
                Unset1D(XV);
                SetLength(BV, N);
                for i_ := 0 to N-1 do
                begin
                    BV[i_] := B[i_,0];
                end;
                SPDMatrixSolve(A, N, IsUpper, BV, Info, Rep, XV);
                SPDErrors := SPDErrors or  not RMatrixCheckSolution(XE, N, Threshold, Info, Rep, XV);
                Info := 0;
                UnsetRep(Rep);
                Unset2D(X);
                SPDMatrixCholeskySolveM(CHA, N, IsUpper, B, M, Info, Rep, X);
                SPDErrors := SPDErrors or  not RMatrixCheckSolutionM(XE, N, M, Threshold, Info, Rep, X);
                Info := 0;
                UnsetRep(Rep);
                Unset1D(XV);
                SetLength(BV, N);
                for i_ := 0 to N-1 do
                begin
                    BV[i_] := B[i_,0];
                end;
                SPDMatrixCholeskySolve(CHA, N, IsUpper, BV, Info, Rep, XV);
                SPDErrors := SPDErrors or  not RMatrixCheckSolution(XE, N, Threshold, Info, Rep, XV);
                
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
                while TaskKind<=3 do
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
                            J:=I;
                            while J<=N-1 do
                            begin
                                A[I,J] := 2*RandomReal-1;
                                A[J,I] := A[I,J];
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        K := RandomInteger(N);
                        for i_ := 0 to N-1 do
                        begin
                            A[i_,K] := 0*A[i_,K];
                        end;
                        APVMul(@A[K][0], 0, N-1, 0);
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
                            J:=I;
                            while J<=N-1 do
                            begin
                                A[I,J] := 2*RandomReal-1;
                                A[J,I] := A[I,J];
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        K := RandomInteger(N);
                        APVMul(@A[K][0], 0, N-1, 0);
                        for i_ := 0 to N-1 do
                        begin
                            A[i_,K] := 0*A[i_,K];
                        end;
                    end;
                    if TaskKind=3 then
                    begin
                        
                        //
                        // equal columns/rows
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
                            J:=I;
                            while J<=N-1 do
                            begin
                                A[I,J] := 2*RandomReal-1;
                                A[J,I] := A[I,J];
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        K := 1+RandomInteger(N-1);
                        for i_ := 0 to N-1 do
                        begin
                            A[i_,0] := A[i_,K];
                        end;
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
                    RMatrixMakeACopy(A, N, N, CHA);
                    RMatrixDropHalf(A, N, IsUpper);
                    RMatrixDropHalf(CHA, N, IsUpper);
                    
                    //
                    // Test SPDMatrixSolveM()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    Unset2D(X);
                    SPDMatrixSolveM(A, N, IsUpper, B, M, Info, Rep, X);
                    SPDErrors := SPDErrors or  not RMatrixCheckSingularM(N, M, Info, Rep, X);
                    
                    //
                    // Test SPDMatrixSolve()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    Unset2D(X);
                    SetLength(BV, N);
                    for i_ := 0 to N-1 do
                    begin
                        BV[i_] := B[i_,0];
                    end;
                    SPDMatrixSolve(A, N, IsUpper, BV, Info, Rep, XV);
                    SPDErrors := SPDErrors or  not RMatrixCheckSingular(N, Info, Rep, XV);
                    
                    //
                    // 'equal columns/rows' are degenerate, but
                    // Cholesky matrix with equal columns/rows IS NOT degenerate,
                    // so it is not used for testing purposes.
                    //
                    if TaskKind<>3 then
                    begin
                        
                        //
                        // Test SPDMatrixLUSolveM()
                        //
                        Info := 0;
                        UnsetRep(Rep);
                        Unset2D(X);
                        SPDMatrixCholeskySolveM(CHA, N, IsUpper, B, M, Info, Rep, X);
                        SPDErrors := SPDErrors or  not RMatrixCheckSingularM(N, M, Info, Rep, X);
                        
                        //
                        // Test SPDMatrixLUSolve()
                        //
                        Info := 0;
                        UnsetRep(Rep);
                        Unset2D(X);
                        SetLength(BV, N);
                        for i_ := 0 to N-1 do
                        begin
                            BV[i_] := B[i_,0];
                        end;
                        SPDMatrixCholeskySolve(CHA, N, IsUpper, BV, Info, Rep, XV);
                        SPDErrors := SPDErrors or  not RMatrixCheckSingular(N, Info, Rep, XV);
                    end;
                    Inc(TaskKind);
                end;
                Inc(M);
            end;
            Inc(N);
        end;
        Inc(Pass);
    end;
end;


(*************************************************************************
Real test
*************************************************************************)
procedure TestCSolver(MaxN : AlglibInteger;
     MaxM : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var CErrors : Boolean;
     var RfsErrors : Boolean);
var
    A : TComplex2DArray;
    LUA : TComplex2DArray;
    ATmp : TComplex2DArray;
    P : TInteger1DArray;
    XE : TComplex2DArray;
    B : TComplex2DArray;
    BV : TComplex1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    N : AlglibInteger;
    M : AlglibInteger;
    Pass : AlglibInteger;
    TaskKind : AlglibInteger;
    MX : Double;
    VErr : Double;
    V : Complex;
    Info : AlglibInteger;
    Rep : DenseSolverReport;
    RepLS : DenseSolverLSReport;
    X : TComplex2DArray;
    XV : TComplex1DArray;
    Y : TComplex1DArray;
    TX : TReal1DArray;
    i_ : AlglibInteger;
begin
    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
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
                // 1. generate random well conditioned matrix A.
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods on original A
                //
                CMatrixRndCond(N, 1000, A);
                CMatrixMakeACopy(A, N, N, LUA);
                CMatrixLU(LUA, N, N, P);
                SetLength(XE, N, M);
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=M-1 do
                    begin
                        XE[I,J].X := 2*RandomReal-1;
                        XE[I,J].Y := 2*RandomReal-1;
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
                        V := C_Complex(0.0);
                        for i_ := 0 to N-1 do
                        begin
                            V := C_Add(V,C_Mul(A[I,i_],XE[i_,J]));
                        end;
                        B[I,J] := V;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                
                //
                // Test solvers
                //
                Info := 0;
                UnsetRep(Rep);
                CUnset2D(X);
                CMatrixSolveM(A, N, B, M, AP_FP_Greater(RandomReal,0.5), Info, Rep, X);
                CErrors := CErrors or  not CMatrixCheckSolutionM(XE, N, M, Threshold, Info, Rep, X);
                Info := 0;
                UnsetRep(Rep);
                CUnset1D(XV);
                SetLength(BV, N);
                for i_ := 0 to N-1 do
                begin
                    BV[i_] := B[i_,0];
                end;
                CMatrixSolve(A, N, BV, Info, Rep, XV);
                CErrors := CErrors or  not CMatrixCheckSolution(XE, N, Threshold, Info, Rep, XV);
                Info := 0;
                UnsetRep(Rep);
                CUnset2D(X);
                CMatrixLUSolveM(LUA, P, N, B, M, Info, Rep, X);
                CErrors := CErrors or  not CMatrixCheckSolutionM(XE, N, M, Threshold, Info, Rep, X);
                Info := 0;
                UnsetRep(Rep);
                CUnset1D(XV);
                SetLength(BV, N);
                for i_ := 0 to N-1 do
                begin
                    BV[i_] := B[i_,0];
                end;
                CMatrixLUSolve(LUA, P, N, BV, Info, Rep, XV);
                CErrors := CErrors or  not CMatrixCheckSolution(XE, N, Threshold, Info, Rep, XV);
                Info := 0;
                UnsetRep(Rep);
                CUnset2D(X);
                CMatrixMixedSolveM(A, LUA, P, N, B, M, Info, Rep, X);
                CErrors := CErrors or  not CMatrixCheckSolutionM(XE, N, M, Threshold, Info, Rep, X);
                Info := 0;
                UnsetRep(Rep);
                CUnset1D(XV);
                SetLength(BV, N);
                for i_ := 0 to N-1 do
                begin
                    BV[i_] := B[i_,0];
                end;
                CMatrixMixedSolve(A, LUA, P, N, BV, Info, Rep, XV);
                CErrors := CErrors or  not CMatrixCheckSolution(XE, N, Threshold, Info, Rep, XV);
                
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
                    CUnset2D(A);
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
                                A[I,J] := C_Complex(0);
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
                                A[I,J].X := 2*RandomReal-1;
                                A[I,J].Y := 2*RandomReal-1;
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        K := RandomInteger(N);
                        for i_ := 0 to N-1 do
                        begin
                            A[i_,K] := C_MulR(A[i_,K],0);
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
                                A[I,J].X := 2*RandomReal-1;
                                A[I,J].Y := 2*RandomReal-1;
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        K := RandomInteger(N);
                        for i_ := 0 to N-1 do
                        begin
                            A[K,i_] := C_MulR(A[K,i_],0);
                        end;
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
                                A[I,J].X := 2*RandomReal-1;
                                A[I,J].Y := 2*RandomReal-1;
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
                                A[I,J].X := 2*RandomReal-1;
                                A[I,J].Y := 2*RandomReal-1;
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        K := 1+RandomInteger(N-1);
                        for i_ := 0 to N-1 do
                        begin
                            A[0,i_] := A[K,i_];
                        end;
                    end;
                    SetLength(XE, N, M);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=M-1 do
                        begin
                            XE[I,J] := C_Complex(2*RandomReal-1);
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
                            V := C_Complex(0.0);
                            for i_ := 0 to N-1 do
                            begin
                                V := C_Add(V,C_Mul(A[I,i_],XE[i_,J]));
                            end;
                            B[I,J] := V;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    CMatrixMakeACopy(A, N, N, LUA);
                    CMatrixLU(LUA, N, N, P);
                    
                    //
                    // Test CMatrixSolveM()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    CUnset2D(X);
                    CMatrixSolveM(A, N, B, M, AP_FP_Greater(RandomReal,0.5), Info, Rep, X);
                    CErrors := CErrors or  not CMatrixCheckSingularM(N, M, Info, Rep, X);
                    
                    //
                    // Test CMatrixSolve()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    CUnset2D(X);
                    SetLength(BV, N);
                    for i_ := 0 to N-1 do
                    begin
                        BV[i_] := B[i_,0];
                    end;
                    CMatrixSolve(A, N, BV, Info, Rep, XV);
                    CErrors := CErrors or  not CMatrixCheckSingular(N, Info, Rep, XV);
                    
                    //
                    // Test CMatrixLUSolveM()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    CUnset2D(X);
                    CMatrixLUSolveM(LUA, P, N, B, M, Info, Rep, X);
                    CErrors := CErrors or  not CMatrixCheckSingularM(N, M, Info, Rep, X);
                    
                    //
                    // Test CMatrixLUSolve()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    CUnset2D(X);
                    SetLength(BV, N);
                    for i_ := 0 to N-1 do
                    begin
                        BV[i_] := B[i_,0];
                    end;
                    CMatrixLUSolve(LUA, P, N, BV, Info, Rep, XV);
                    CErrors := CErrors or  not CMatrixCheckSingular(N, Info, Rep, XV);
                    
                    //
                    // Test CMatrixMixedSolveM()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    CUnset2D(X);
                    CMatrixMixedSolveM(A, LUA, P, N, B, M, Info, Rep, X);
                    CErrors := CErrors or  not CMatrixCheckSingularM(N, M, Info, Rep, X);
                    
                    //
                    // Test CMatrixMixedSolve()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    CUnset2D(X);
                    SetLength(BV, N);
                    for i_ := 0 to N-1 do
                    begin
                        BV[i_] := B[i_,0];
                    end;
                    CMatrixMixedSolve(A, LUA, P, N, BV, Info, Rep, XV);
                    CErrors := CErrors or  not CMatrixCheckSingular(N, Info, Rep, XV);
                    Inc(TaskKind);
                end;
                Inc(M);
            end;
            Inc(N);
        end;
        Inc(Pass);
    end;
    
    //
    // test iterative improvement
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        
        //
        // Test iterative improvement matrices
        //
        // A matrix/right part are constructed such that both matrix
        // and solution components magnitudes are within (-1,+1).
        // Such matrix/right part have nice properties - system can
        // be solved using iterative improvement with |A*x-b| about
        // several ulps of max(1,|b|).
        //
        N := 100;
        SetLength(A, N, N);
        SetLength(B, N, 1);
        SetLength(BV, N);
        SetLength(TX, 2*N);
        SetLength(XV, N);
        SetLength(Y, N);
        I:=0;
        while I<=N-1 do
        begin
            XV[I].X := 2*RandomReal-1;
            XV[I].Y := 2*RandomReal-1;
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                A[I,J].X := 2*RandomReal-1;
                A[I,J].Y := 2*RandomReal-1;
                Inc(J);
            end;
            for i_ := 0 to N-1 do
            begin
                Y[i_] := A[I,i_];
            end;
            XCDot(Y, XV, N, TX, V, VErr);
            BV[I] := V;
            Inc(I);
        end;
        for i_ := 0 to N-1 do
        begin
            B[i_,0] := BV[i_];
        end;
        
        //
        // Test CMatrixSolveM()
        //
        CUnset2D(X);
        CMatrixSolveM(A, N, B, 1, True, Info, Rep, X);
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
                for i_ := 0 to N-1 do
                begin
                    Y[i_] := A[I,i_];
                end;
                XCDot(Y, XV, N, TX, V, VErr);
                RfsErrors := RfsErrors or AP_FP_Greater(AbsComplex(C_Sub(V,B[I,0])),8*MachineEpsilon*Max(1, AbsComplex(B[I,0])));
                Inc(I);
            end;
        end;
        
        //
        // Test CMatrixSolve()
        //
        CUnset1D(XV);
        CMatrixSolve(A, N, BV, Info, Rep, XV);
        if Info<=0 then
        begin
            RfsErrors := True;
        end
        else
        begin
            I:=0;
            while I<=N-1 do
            begin
                for i_ := 0 to N-1 do
                begin
                    Y[i_] := A[I,i_];
                end;
                XCDot(Y, XV, N, TX, V, VErr);
                RfsErrors := RfsErrors or AP_FP_Greater(AbsComplex(C_Sub(V,BV[I])),8*MachineEpsilon*Max(1, AbsComplex(BV[I])));
                Inc(I);
            end;
        end;
        
        //
        // TODO: Test LS-solver on the same matrix
        //
        Inc(Pass);
    end;
end;


(*************************************************************************
HPD test
*************************************************************************)
procedure TestHPDSolver(MaxN : AlglibInteger;
     MaxM : AlglibInteger;
     PassCount : AlglibInteger;
     Threshold : Double;
     var HPDErrors : Boolean;
     var RfsErrors : Boolean);
var
    A : TComplex2DArray;
    CHA : TComplex2DArray;
    ATmp : TComplex2DArray;
    P : TInteger1DArray;
    XE : TComplex2DArray;
    B : TComplex2DArray;
    BV : TComplex1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    N : AlglibInteger;
    M : AlglibInteger;
    Pass : AlglibInteger;
    TaskKind : AlglibInteger;
    MX : Double;
    V : Complex;
    IsUpper : Boolean;
    Info : AlglibInteger;
    Rep : DenseSolverReport;
    RepLS : DenseSolverLSReport;
    X : TComplex2DArray;
    XV : TComplex1DArray;
    Y : TComplex1DArray;
    TX : TComplex1DArray;
    i_ : AlglibInteger;
begin
    
    //
    // General square matrices:
    // * test general solvers
    // * test least squares solver
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
                // 1. generate random well conditioned matrix A.
                // 2. generate random solution vector xe
                // 3. generate right part b=A*xe
                // 4. test different methods on original A
                //
                IsUpper := AP_FP_Greater(RandomReal,0.5);
                HPDMatrixRndCond(N, 1000, A);
                CMatrixMakeACopy(A, N, N, CHA);
                if  not HPDMatrixCholesky(CHA, N, IsUpper) then
                begin
                    HPDErrors := True;
                    Exit;
                end;
                SetLength(XE, N, M);
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=M-1 do
                    begin
                        XE[I,J].X := 2*RandomReal-1;
                        XE[I,J].Y := 2*RandomReal-1;
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
                        V := C_Complex(0.0);
                        for i_ := 0 to N-1 do
                        begin
                            V := C_Add(V,C_Mul(A[I,i_],XE[i_,J]));
                        end;
                        B[I,J] := V;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                CMatrixDropHalf(A, N, IsUpper);
                CMatrixDropHalf(CHA, N, IsUpper);
                
                //
                // Test solvers
                //
                Info := 0;
                UnsetRep(Rep);
                CUnset2D(X);
                HPDMatrixSolveM(A, N, IsUpper, B, M, Info, Rep, X);
                HPDErrors := HPDErrors or  not CMatrixCheckSolutionM(XE, N, M, Threshold, Info, Rep, X);
                Info := 0;
                UnsetRep(Rep);
                CUnset1D(XV);
                SetLength(BV, N);
                for i_ := 0 to N-1 do
                begin
                    BV[i_] := B[i_,0];
                end;
                HPDMatrixSolve(A, N, IsUpper, BV, Info, Rep, XV);
                HPDErrors := HPDErrors or  not CMatrixCheckSolution(XE, N, Threshold, Info, Rep, XV);
                Info := 0;
                UnsetRep(Rep);
                CUnset2D(X);
                HPDMatrixCholeskySolveM(CHA, N, IsUpper, B, M, Info, Rep, X);
                HPDErrors := HPDErrors or  not CMatrixCheckSolutionM(XE, N, M, Threshold, Info, Rep, X);
                Info := 0;
                UnsetRep(Rep);
                CUnset1D(XV);
                SetLength(BV, N);
                for i_ := 0 to N-1 do
                begin
                    BV[i_] := B[i_,0];
                end;
                HPDMatrixCholeskySolve(CHA, N, IsUpper, BV, Info, Rep, XV);
                HPDErrors := HPDErrors or  not CMatrixCheckSolution(XE, N, Threshold, Info, Rep, XV);
                
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
                while TaskKind<=3 do
                begin
                    CUnset2D(A);
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
                                A[I,J] := C_Complex(0);
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
                            J:=I;
                            while J<=N-1 do
                            begin
                                A[I,J].X := 2*RandomReal-1;
                                A[I,J].Y := 2*RandomReal-1;
                                if I=J then
                                begin
                                    A[I,J].Y := 0;
                                end;
                                A[J,I] := A[I,J];
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        K := RandomInteger(N);
                        for i_ := 0 to N-1 do
                        begin
                            A[i_,K] := C_MulR(A[i_,K],0);
                        end;
                        for i_ := 0 to N-1 do
                        begin
                            A[K,i_] := C_MulR(A[K,i_],0);
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
                            J:=I;
                            while J<=N-1 do
                            begin
                                A[I,J].X := 2*RandomReal-1;
                                A[I,J].Y := 2*RandomReal-1;
                                if I=J then
                                begin
                                    A[I,J].Y := 0;
                                end;
                                A[J,I] := A[I,J];
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        K := RandomInteger(N);
                        for i_ := 0 to N-1 do
                        begin
                            A[K,i_] := C_MulR(A[K,i_],0);
                        end;
                        for i_ := 0 to N-1 do
                        begin
                            A[i_,K] := C_MulR(A[i_,K],0);
                        end;
                    end;
                    if TaskKind=3 then
                    begin
                        
                        //
                        // equal columns/rows
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
                            J:=I;
                            while J<=N-1 do
                            begin
                                A[I,J].X := 2*RandomReal-1;
                                A[I,J].Y := 2*RandomReal-1;
                                if I=J then
                                begin
                                    A[I,J].Y := 0;
                                end;
                                A[J,I] := A[I,J];
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        K := 1+RandomInteger(N-1);
                        for i_ := 0 to N-1 do
                        begin
                            A[i_,0] := A[i_,K];
                        end;
                        for i_ := 0 to N-1 do
                        begin
                            A[0,i_] := A[K,i_];
                        end;
                    end;
                    SetLength(XE, N, M);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=M-1 do
                        begin
                            XE[I,J] := C_Complex(2*RandomReal-1);
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
                            V := C_Complex(0.0);
                            for i_ := 0 to N-1 do
                            begin
                                V := C_Add(V,C_Mul(A[I,i_],XE[i_,J]));
                            end;
                            B[I,J] := V;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    CMatrixMakeACopy(A, N, N, CHA);
                    CMatrixDropHalf(A, N, IsUpper);
                    CMatrixDropHalf(CHA, N, IsUpper);
                    
                    //
                    // Test SPDMatrixSolveM()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    CUnset2D(X);
                    HPDMatrixSolveM(A, N, IsUpper, B, M, Info, Rep, X);
                    HPDErrors := HPDErrors or  not CMatrixCheckSingularM(N, M, Info, Rep, X);
                    
                    //
                    // Test SPDMatrixSolve()
                    //
                    Info := 0;
                    UnsetRep(Rep);
                    CUnset2D(X);
                    SetLength(BV, N);
                    for i_ := 0 to N-1 do
                    begin
                        BV[i_] := B[i_,0];
                    end;
                    HPDMatrixSolve(A, N, IsUpper, BV, Info, Rep, XV);
                    HPDErrors := HPDErrors or  not CMatrixCheckSingular(N, Info, Rep, XV);
                    
                    //
                    // 'equal columns/rows' are degenerate, but
                    // Cholesky matrix with equal columns/rows IS NOT degenerate,
                    // so it is not used for testing purposes.
                    //
                    if TaskKind<>3 then
                    begin
                        
                        //
                        // Test SPDMatrixLUSolveM()
                        //
                        Info := 0;
                        UnsetRep(Rep);
                        CUnset2D(X);
                        HPDMatrixCholeskySolveM(CHA, N, IsUpper, B, M, Info, Rep, X);
                        HPDErrors := HPDErrors or  not CMatrixCheckSingularM(N, M, Info, Rep, X);
                        
                        //
                        // Test SPDMatrixLUSolve()
                        //
                        Info := 0;
                        UnsetRep(Rep);
                        CUnset2D(X);
                        SetLength(BV, N);
                        for i_ := 0 to N-1 do
                        begin
                            BV[i_] := B[i_,0];
                        end;
                        HPDMatrixCholeskySolve(CHA, N, IsUpper, BV, Info, Rep, XV);
                        HPDErrors := HPDErrors or  not CMatrixCheckSingular(N, Info, Rep, XV);
                    end;
                    Inc(TaskKind);
                end;
                Inc(M);
            end;
            Inc(N);
        end;
        Inc(Pass);
    end;
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
Unsets real matrix
*************************************************************************)
procedure CUnset2D(var X : TComplex2DArray);
begin
    SetLength(X, 1, 1);
    X[0,0] := C_Complex(2*RandomReal-1);
end;


(*************************************************************************
Unsets real vector
*************************************************************************)
procedure CUnset1D(var X : TComplex1DArray);
begin
    SetLength(X, 1);
    X[0] := C_Complex(2*RandomReal-1);
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