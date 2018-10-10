(*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************)
unit minlm;
interface
uses Math, Sysutils, Ap, blas, reflections, creflections, hqrnd, matgen, trinverse, ablasf, ablas, trfac, bidiagonal, qr, lq, rotations, bdsvd, svd, trlinsolve, safesolve, rcond, tsort, xblas, densesolver, lbfgs;

type
LMState = record
    WrongParams : Boolean;
    N : AlglibInteger;
    M : AlglibInteger;
    EpsF : Double;
    EpsX : Double;
    MaxIts : AlglibInteger;
    Flags : AlglibInteger;
    UserMode : AlglibInteger;
    X : TReal1DArray;
    F : Double;
    FI : TReal1DArray;
    J : TReal2DArray;
    H : TReal2DArray;
    G : TReal1DArray;
    NeedF : Boolean;
    NeedFG : Boolean;
    NeedFGH : Boolean;
    NeedFiJ : Boolean;
    XUpdated : Boolean;
    InternalState : LBFGSState;
    InternalRep : LBFGSReport;
    XPrec : TReal1DArray;
    XBase : TReal1DArray;
    XDir : TReal1DArray;
    GBase : TReal1DArray;
    RawModel : TReal2DArray;
    Model : TReal2DArray;
    WORK : TReal1DArray;
    RState : RCommState;
    RepIterationsCount : AlglibInteger;
    RepTerminationType : AlglibInteger;
    RepNFunc : AlglibInteger;
    RepNJac : AlglibInteger;
    RepNGrad : AlglibInteger;
    RepNHess : AlglibInteger;
    RepNCholesky : AlglibInteger;
    SolverInfo : AlglibInteger;
    SolverRep : DenseSolverReport;
end;


LMReport = record
    IterationsCount : AlglibInteger;
    TerminationType : AlglibInteger;
    NFunc : AlglibInteger;
    NJac : AlglibInteger;
    NGrad : AlglibInteger;
    NHess : AlglibInteger;
    NCholesky : AlglibInteger;
end;



procedure MinLMFGH(const N : AlglibInteger;
     const X : TReal1DArray;
     const EpsF : Double;
     const EpsX : Double;
     const MaxIts : AlglibInteger;
     var State : LMState);
procedure MinLMFGJ(const N : AlglibInteger;
     const M : AlglibInteger;
     const X : TReal1DArray;
     const EpsF : Double;
     const EpsX : Double;
     const MaxIts : AlglibInteger;
     var State : LMState);
procedure MinLMFJ(const N : AlglibInteger;
     const M : AlglibInteger;
     const X : TReal1DArray;
     const EpsF : Double;
     const EpsX : Double;
     const MaxIts : AlglibInteger;
     var State : LMState);
function MinLMIteration(var State : LMState):Boolean;
procedure MinLMResults(const State : LMState;
     var X : TReal1DArray;
     var Rep : LMReport);

implementation

const
    LMModeFJ = 0;
    LMModeFGJ = 1;
    LMModeFGH = 2;
    LMFlagNoPreLBFGS = 1;
    LMFlagNoIntLBFGS = 2;
    LMPreLBFGSM = 5;
    LMIntLBFGSIts = 5;
    LBFGSNoRealloc = 1;

procedure LMPrepare(N : AlglibInteger;
     M : AlglibInteger;
     HaveGrad : Boolean;
     var State : LMState);forward;
procedure LMClearRequestFields(var State : LMState);forward;


(*************************************************************************
    LEVENBERG-MARQUARDT-LIKE METHOD FOR NON-LINEAR OPTIMIZATION

Optimization using function gradient and Hessian.  Algorithm -  Levenberg-
Marquardt   modification   with   L-BFGS   pre-optimization  and  internal
pre-conditioned L-BFGS optimization after each Levenberg-Marquardt step.

Function F has general form (not "sum-of-squares"):

    F = F(x[0], ..., x[n-1])

EXAMPLE

See HTML-documentation.

INPUT PARAMETERS:
    N       -   dimension, N>1
    X       -   initial solution, array[0..N-1]
    EpsF    -   stopping criterion. Algorithm stops if
                |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}
    EpsX    -   stopping criterion. Algorithm stops if
                |X(k+1)-X(k)| <= EpsX*(1+|X(k)|)
    MaxIts  -   stopping criterion. Algorithm stops after MaxIts iterations.
                MaxIts=0 means no stopping criterion.

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state between subsequent
                calls of MinLMIteration. Used for reverse communication.
                This structure should be passed to MinLMIteration subroutine.

See also MinLMIteration, MinLMResults.

NOTE

Passing EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to automatic
stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure MinLMFGH(const N : AlglibInteger;
     const X : TReal1DArray;
     const EpsF : Double;
     const EpsX : Double;
     const MaxIts : AlglibInteger;
     var State : LMState);
begin
    
    //
    // Prepare RComm
    //
    SetLength(State.RState.IA, 3+1);
    SetLength(State.RState.BA, 0+1);
    SetLength(State.RState.RA, 8+1);
    State.RState.Stage := -1;
    
    //
    // prepare internal structures
    //
    LMPrepare(N, 0, True, State);
    
    //
    // initialize, check parameters
    //
    State.XUpdated := False;
    State.N := N;
    State.M := 0;
    State.EpsF := EpsF;
    State.EpsX := EpsX;
    State.MaxIts := MaxIts;
    State.Flags := 0;
    if AP_FP_Eq(State.EpsF,0) and AP_FP_Eq(State.EpsX,0) and (State.MaxIts=0) then
    begin
        State.EpsX := 1.0E-6;
    end;
    State.UserMode := LMModeFGH;
    State.WrongParams := False;
    if (N<1) or AP_FP_Less(EpsF,0) or AP_FP_Less(EpsX,0) or (MaxIts<0) then
    begin
        State.WrongParams := True;
        Exit;
    end;
    APVMove(@State.X[0], 0, N-1, @X[0], 0, N-1);
end;


(*************************************************************************
    LEVENBERG-MARQUARDT-LIKE METHOD FOR NON-LINEAR OPTIMIZATION

Optimization using function gradient and Jacobian.  Algorithm -  Levenberg-
Marquardt   modification   with   L-BFGS   pre-optimization  and  internal
pre-conditioned L-BFGS optimization after each Levenberg-Marquardt step.

Function F is represented as sum of squares:

    F = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])

EXAMPLE

See HTML-documentation.

INPUT PARAMETERS:
    N       -   dimension, N>1
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]
    EpsF    -   stopping criterion. Algorithm stops if
                |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}
    EpsX    -   stopping criterion. Algorithm stops if
                |X(k+1)-X(k)| <= EpsX*(1+|X(k)|)
    MaxIts  -   stopping criterion. Algorithm stops after MaxIts iterations.
                MaxIts=0 means no stopping criterion.

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state between subsequent
                calls of MinLMIteration. Used for reverse communication.
                This structure should be passed to MinLMIteration subroutine.

See also MinLMIteration, MinLMResults.

NOTE

Passing EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to automatic
stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure MinLMFGJ(const N : AlglibInteger;
     const M : AlglibInteger;
     const X : TReal1DArray;
     const EpsF : Double;
     const EpsX : Double;
     const MaxIts : AlglibInteger;
     var State : LMState);
begin
    
    //
    // Prepare RComm
    //
    SetLength(State.RState.IA, 3+1);
    SetLength(State.RState.BA, 0+1);
    SetLength(State.RState.RA, 8+1);
    State.RState.Stage := -1;
    
    //
    // prepare internal structures
    //
    LMPrepare(N, M, True, State);
    
    //
    // initialize, check parameters
    //
    State.XUpdated := False;
    State.N := N;
    State.M := M;
    State.EpsF := EpsF;
    State.EpsX := EpsX;
    State.MaxIts := MaxIts;
    State.Flags := 0;
    if AP_FP_Eq(State.EpsF,0) and AP_FP_Eq(State.EpsX,0) and (State.MaxIts=0) then
    begin
        State.EpsX := 1.0E-6;
    end;
    State.UserMode := LMModeFGJ;
    State.WrongParams := False;
    if (N<1) or (M<1) or AP_FP_Less(EpsF,0) or AP_FP_Less(EpsX,0) or (MaxIts<0) then
    begin
        State.WrongParams := True;
        Exit;
    end;
    APVMove(@State.X[0], 0, N-1, @X[0], 0, N-1);
end;


(*************************************************************************
    CLASSIC LEVENBERG-MARQUARDT METHOD FOR NON-LINEAR OPTIMIZATION

Optimization using Jacobi matrix. Algorithm  -  classic Levenberg-Marquardt
method.

Function F is represented as sum of squares:

    F = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])

EXAMPLE

See HTML-documentation.

INPUT PARAMETERS:
    N       -   dimension, N>1
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]
    EpsF    -   stopping criterion. Algorithm stops if
                |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}
    EpsX    -   stopping criterion. Algorithm stops if
                |X(k+1)-X(k)| <= EpsX*(1+|X(k)|)
    MaxIts  -   stopping criterion. Algorithm stops after MaxIts iterations.
                MaxIts=0 means no stopping criterion.

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state between subsequent
                calls of MinLMIteration. Used for reverse communication.
                This structure should be passed to MinLMIteration subroutine.

See also MinLMIteration, MinLMResults.

NOTE

Passing EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to automatic
stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure MinLMFJ(const N : AlglibInteger;
     const M : AlglibInteger;
     const X : TReal1DArray;
     const EpsF : Double;
     const EpsX : Double;
     const MaxIts : AlglibInteger;
     var State : LMState);
begin
    
    //
    // Prepare RComm
    //
    SetLength(State.RState.IA, 3+1);
    SetLength(State.RState.BA, 0+1);
    SetLength(State.RState.RA, 8+1);
    State.RState.Stage := -1;
    
    //
    // prepare internal structures
    //
    LMPrepare(N, M, True, State);
    
    //
    // initialize, check parameters
    //
    State.XUpdated := False;
    State.N := N;
    State.M := M;
    State.EpsF := EpsF;
    State.EpsX := EpsX;
    State.MaxIts := MaxIts;
    State.Flags := 0;
    if AP_FP_Eq(State.EpsF,0) and AP_FP_Eq(State.EpsX,0) and (State.MaxIts=0) then
    begin
        State.EpsX := 1.0E-6;
    end;
    State.UserMode := LMModeFJ;
    State.WrongParams := False;
    if (N<1) or (M<1) or AP_FP_Less(EpsF,0) or AP_FP_Less(EpsX,0) or (MaxIts<0) then
    begin
        State.WrongParams := True;
        Exit;
    end;
    APVMove(@State.X[0], 0, N-1, @X[0], 0, N-1);
end;


(*************************************************************************
One Levenberg-Marquardt iteration.

Called after inialization of State structure with MinLMXXX subroutine.
See HTML docs for examples.

Input parameters:
    State   -   structure which stores algorithm state between subsequent
                calls and which is used for reverse communication. Must be
                initialized with MinLMXXX call first.

If subroutine returned False, iterative algorithm has converged.

If subroutine returned True, then:
* if State.NeedF=True,      -   function value F at State.X[0..N-1]
                                is required
* if State.NeedFG=True      -   function value F and gradient G
                                are required
* if State.NeedFiJ=True     -   function vector f[i] and Jacobi matrix J
                                are required
* if State.NeedFGH=True     -   function value F, gradient G and Hesian H
                                are required

One and only one of this fields can be set at time.

Results are stored:
* function value            -   in LMState.F
* gradient                  -   in LMState.G[0..N-1]
* Jacobi matrix             -   in LMState.J[0..M-1,0..N-1]
* Hessian                   -   in LMState.H[0..N-1,0..N-1]

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************)
function MinLMIteration(var State : LMState):Boolean;
var
    N : AlglibInteger;
    M : AlglibInteger;
    I : AlglibInteger;
    XNorm : Double;
    StepNorm : Double;
    SPD : Boolean;
    FBase : Double;
    FNew : Double;
    Lambda : Double;
    Nu : Double;
    LambdaUp : Double;
    LambdaDown : Double;
    LBFGSFlags : AlglibInteger;
    V : Double;
label
lbl_9, lbl_0, lbl_10, lbl_7, lbl_1, lbl_11, lbl_2, lbl_13, lbl_15, lbl_3, lbl_21, lbl_4, lbl_22, lbl_19, lbl_17, lbl_5, lbl_23, lbl_6, lbl_25, lbl_16, lbl_rcomm;
begin
    
    //
    // Reverse communication preparations
    // I know it looks ugly, but it works the same way
    // anywhere from C++ to Python.
    //
    // This code initializes locals by:
    // * random values determined during code
    //   generation - on first subroutine call
    // * values from previous call - on subsequent calls
    //
    if State.RState.Stage>=0 then
    begin
        N := State.RState.IA[0];
        M := State.RState.IA[1];
        I := State.RState.IA[2];
        LBFGSFlags := State.RState.IA[3];
        SPD := State.RState.BA[0];
        XNorm := State.RState.RA[0];
        StepNorm := State.RState.RA[1];
        FBase := State.RState.RA[2];
        FNew := State.RState.RA[3];
        Lambda := State.RState.RA[4];
        Nu := State.RState.RA[5];
        LambdaUp := State.RState.RA[6];
        LambdaDown := State.RState.RA[7];
        V := State.RState.RA[8];
    end
    else
    begin
        N := -983;
        M := -989;
        I := -834;
        LBFGSFlags := 900;
        SPD := True;
        XNorm := 364;
        StepNorm := 214;
        FBase := -338;
        FNew := -686;
        Lambda := 912;
        Nu := 585;
        LambdaUp := 497;
        LambdaDown := -271;
        V := -581;
    end;
    if State.RState.Stage=0 then
    begin
        goto lbl_0;
    end;
    if State.RState.Stage=1 then
    begin
        goto lbl_1;
    end;
    if State.RState.Stage=2 then
    begin
        goto lbl_2;
    end;
    if State.RState.Stage=3 then
    begin
        goto lbl_3;
    end;
    if State.RState.Stage=4 then
    begin
        goto lbl_4;
    end;
    if State.RState.Stage=5 then
    begin
        goto lbl_5;
    end;
    if State.RState.Stage=6 then
    begin
        goto lbl_6;
    end;
    
    //
    // Routine body
    //
    Assert((State.UserMode=LMModeFJ) or (State.UserMode=LMModeFGJ) or (State.UserMode=LMModeFGH), 'LM: internal error');
    if State.WrongParams then
    begin
        State.RepTerminationType := -1;
        Result := False;
        Exit;
    end;
    
    //
    // prepare params
    //
    N := State.N;
    M := State.M;
    LambdaUp := 10;
    LambdaDown := 0.3;
    Nu := 2;
    LBFGSFlags := 0;
    
    //
    // if we have F and G
    //
    if not(((State.UserMode=LMModeFGJ) or (State.UserMode=LMModeFGH)) and (State.Flags div LMFlagNoPreLBFGS mod 2=0)) then
    begin
        goto lbl_7;
    end;
    
    //
    // First stage of the hybrid algorithm: LBFGS
    //
    MinLBFGS(N, Min(N, LMPreLBFGSM), State.X, 0.0, 0.0, 0.0, Max(5, N), 0, State.InternalState);
lbl_9:
    if not MinLBFGSIteration(State.InternalState) then
    begin
        goto lbl_10;
    end;
    
    //
    // RComm
    //
    APVMove(@State.X[0], 0, N-1, @State.InternalState.X[0], 0, N-1);
    LMClearRequestFields(State);
    State.NeedFG := True;
    State.RState.Stage := 0;
    goto lbl_rcomm;
lbl_0:
    State.RepNFunc := State.RepNFunc+1;
    State.RepNGrad := State.RepNGrad+1;
    
    //
    // Call LBFGS
    //
    State.InternalState.F := State.F;
    APVMove(@State.InternalState.G[0], 0, N-1, @State.G[0], 0, N-1);
    goto lbl_9;
lbl_10:
    MinLBFGSResults(State.InternalState, State.X, State.InternalRep);
lbl_7:
    
    //
    // Second stage of the hybrid algorithm: LM
    // Initialize quadratic model.
    //
    if State.UserMode<>LMModeFGH then
    begin
        goto lbl_11;
    end;
    
    //
    // RComm
    //
    LMClearRequestFields(State);
    State.NeedFGH := True;
    State.RState.Stage := 1;
    goto lbl_rcomm;
lbl_1:
    State.RepNFunc := State.RepNFunc+1;
    State.RepNGrad := State.RepNGrad+1;
    State.RepNHess := State.RepNHess+1;
    
    //
    // generate raw quadratic model
    //
    I:=0;
    while I<=N-1 do
    begin
        APVMove(@State.RawModel[I][0], 0, N-1, @State.H[I][0], 0, N-1);
        Inc(I);
    end;
    APVMove(@State.GBase[0], 0, N-1, @State.G[0], 0, N-1);
    FBase := State.F;
lbl_11:
    if not((State.UserMode=LMModeFGJ) or (State.UserMode=LMModeFJ)) then
    begin
        goto lbl_13;
    end;
    
    //
    // RComm
    //
    LMClearRequestFields(State);
    State.NeedFiJ := True;
    State.RState.Stage := 2;
    goto lbl_rcomm;
lbl_2:
    State.RepNFunc := State.RepNFunc+1;
    State.RepNJac := State.RepNJac+1;
    
    //
    // generate raw quadratic model
    //
    MatrixMatrixMultiply(State.J, 0, M-1, 0, N-1, True, State.J, 0, M-1, 0, N-1, False, 1.0, State.RawModel, 0, N-1, 0, N-1, 0.0, State.WORK);
    MatrixVectorMultiply(State.J, 0, M-1, 0, N-1, True, State.FI, 0, M-1, 1.0, State.GBase, 0, N-1, 0.0);
    FBase := APVDotProduct(@State.FI[0], 0, M-1, @State.FI[0], 0, M-1);
lbl_13:
    Lambda := 0.001;
lbl_15:
    if False then
    begin
        goto lbl_16;
    end;
    
    //
    // 1. Model = RawModel+lambda*I
    // 2. Try to solve (RawModel+Lambda*I)*dx = -g.
    //    Increase lambda if left part is not positive definite.
    //
    I:=0;
    while I<=N-1 do
    begin
        APVMove(@State.Model[I][0], 0, N-1, @State.RawModel[I][0], 0, N-1);
        State.Model[I,I] := State.Model[I,I]+Lambda;
        Inc(I);
    end;
    SPD := SPDMatrixCholesky(State.Model, N, True);
    State.RepNCholesky := State.RepNCholesky+1;
    if  not SPD then
    begin
        Lambda := Lambda*LambdaUp*Nu;
        Nu := Nu*2;
        goto lbl_15;
    end;
    SPDMatrixCholeskySolve(State.Model, N, True, State.GBase, State.SolverInfo, State.SolverRep, State.XDir);
    if State.SolverInfo<0 then
    begin
        Lambda := Lambda*LambdaUp*Nu;
        Nu := Nu*2;
        goto lbl_15;
    end;
    APVMul(@State.XDir[0], 0, N-1, -1);
    
    //
    // Candidate lambda found.
    // 1. Save old w in WBase
    // 1. Test some stopping criterions
    // 2. If error(w+wdir)>error(w), increase lambda
    //
    APVMove(@State.XBase[0], 0, N-1, @State.X[0], 0, N-1);
    APVAdd(@State.X[0], 0, N-1, @State.XDir[0], 0, N-1);
    XNorm := APVDotProduct(@State.XBase[0], 0, N-1, @State.XBase[0], 0, N-1);
    StepNorm := APVDotProduct(@State.XDir[0], 0, N-1, @State.XDir[0], 0, N-1);
    XNorm := Sqrt(XNorm);
    StepNorm := Sqrt(StepNorm);
    if AP_FP_Less_Eq(StepNorm,State.EpsX*(1+XNorm)) then
    begin
        
        //
        // step size if small enough
        //
        State.RepTerminationType := 2;
        goto lbl_16;
    end;
    LMClearRequestFields(State);
    State.NeedF := True;
    State.RState.Stage := 3;
    goto lbl_rcomm;
lbl_3:
    State.RepNFunc := State.RepNFunc+1;
    FNew := State.F;
    if AP_FP_Less_Eq(AbsReal(FNew-FBase),State.EpsF*Max(1, Max(AbsReal(FBase), AbsReal(FNew)))) then
    begin
        
        //
        // function change is small enough
        //
        State.RepTerminationType := 1;
        goto lbl_16;
    end;
    if AP_FP_Greater(FNew,FBase) then
    begin
        
        //
        // restore state and continue out search for lambda
        //
        APVMove(@State.X[0], 0, N-1, @State.XBase[0], 0, N-1);
        Lambda := Lambda*LambdaUp*Nu;
        Nu := Nu*2;
        goto lbl_15;
    end;
    if not(((State.UserMode=LMModeFGJ) or (State.UserMode=LMModeFGH)) and (State.Flags div LMFlagNoIntLBFGS mod 2=0)) then
    begin
        goto lbl_17;
    end;
    
    //
    // Optimize using inv(cholesky(H)) as preconditioner
    //
    if not RMatrixTRInverse(State.Model, N, True, False) then
    begin
        goto lbl_19;
    end;
    
    //
    // if matrix can be inverted use it.
    // just silently move to next iteration otherwise.
    // (will be very, very rare, mostly for specially
    // designed near-degenerate tasks)
    //
    APVMove(@State.XBase[0], 0, N-1, @State.X[0], 0, N-1);
    I:=0;
    while I<=N-1 do
    begin
        State.XPrec[I] := 0;
        Inc(I);
    end;
    MinLBFGS(N, Min(N, LMIntLBFGSIts), State.XPrec, 0.0, 0.0, 0.0, LMIntLBFGSIts, LBFGSFlags, State.InternalState);
lbl_21:
    if not MinLBFGSIteration(State.InternalState) then
    begin
        goto lbl_22;
    end;
    
    //
    // convert XPrec to unpreconditioned form, then call RComm.
    //
    I:=0;
    while I<=N-1 do
    begin
        V := APVDotProduct(@State.InternalState.X[0], I, N-1, @State.Model[I][0], I, N-1);
        State.X[I] := State.XBase[I]+V;
        Inc(I);
    end;
    LMClearRequestFields(State);
    State.NeedFG := True;
    State.RState.Stage := 4;
    goto lbl_rcomm;
lbl_4:
    State.RepNFunc := State.RepNFunc+1;
    State.RepNGrad := State.RepNGrad+1;
    
    //
    // 1. pass State.F to State.InternalState.F
    // 2. convert gradient back to preconditioned form
    //
    State.InternalState.F := State.F;
    I:=0;
    while I<=N-1 do
    begin
        State.InternalState.G[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        V := State.G[I];
        APVAdd(@State.InternalState.G[0], I, N-1, @State.Model[I][0], I, N-1, V);
        Inc(I);
    end;
    
    //
    // next iteration
    //
    goto lbl_21;
lbl_22:
    
    //
    // change LBFGS flags to NoRealloc.
    // L-BFGS subroutine will use memory allocated from previous run.
    // it is possible since all subsequent calls will be with same N/M.
    //
    LBFGSFlags := LBFGSNoRealloc;
    
    //
    // back to unpreconditioned X
    //
    MinLBFGSResults(State.InternalState, State.XPrec, State.InternalRep);
    I:=0;
    while I<=N-1 do
    begin
        V := APVDotProduct(@State.XPrec[0], I, N-1, @State.Model[I][0], I, N-1);
        State.X[I] := State.XBase[I]+V;
        Inc(I);
    end;
lbl_19:
lbl_17:
    
    //
    // Accept new position.
    // Calculate Hessian
    //
    if State.UserMode<>LMModeFGH then
    begin
        goto lbl_23;
    end;
    
    //
    // RComm
    //
    LMClearRequestFields(State);
    State.NeedFGH := True;
    State.RState.Stage := 5;
    goto lbl_rcomm;
lbl_5:
    State.RepNFunc := State.RepNFunc+1;
    State.RepNGrad := State.RepNGrad+1;
    State.RepNHess := State.RepNHess+1;
    
    //
    // Update raw quadratic model
    //
    I:=0;
    while I<=N-1 do
    begin
        APVMove(@State.RawModel[I][0], 0, N-1, @State.H[I][0], 0, N-1);
        Inc(I);
    end;
    APVMove(@State.GBase[0], 0, N-1, @State.G[0], 0, N-1);
    FBase := State.F;
lbl_23:
    if not((State.UserMode=LMModeFGJ) or (State.UserMode=LMModeFJ)) then
    begin
        goto lbl_25;
    end;
    
    //
    // RComm
    //
    LMClearRequestFields(State);
    State.NeedFiJ := True;
    State.RState.Stage := 6;
    goto lbl_rcomm;
lbl_6:
    State.RepNFunc := State.RepNFunc+1;
    State.RepNJac := State.RepNJac+1;
    
    //
    // generate raw quadratic model
    //
    MatrixMatrixMultiply(State.J, 0, M-1, 0, N-1, True, State.J, 0, M-1, 0, N-1, False, 1.0, State.RawModel, 0, N-1, 0, N-1, 0.0, State.WORK);
    MatrixVectorMultiply(State.J, 0, M-1, 0, N-1, True, State.FI, 0, M-1, 1.0, State.GBase, 0, N-1, 0.0);
    FBase := APVDotProduct(@State.FI[0], 0, M-1, @State.FI[0], 0, M-1);
lbl_25:
    State.RepIterationsCount := State.RepIterationsCount+1;
    if (State.RepIterationsCount>=State.MaxIts) and (State.MaxIts>0) then
    begin
        State.RepTerminationType := 5;
        goto lbl_16;
    end;
    
    //
    // Update lambda
    //
    Lambda := Lambda*LambdaDown;
    Nu := 2;
    goto lbl_15;
lbl_16:
    Result := False;
    Exit;
    
    //
    // Saving state
    //
lbl_rcomm:
    Result := True;
    State.RState.IA[0] := N;
    State.RState.IA[1] := M;
    State.RState.IA[2] := I;
    State.RState.IA[3] := LBFGSFlags;
    State.RState.BA[0] := SPD;
    State.RState.RA[0] := XNorm;
    State.RState.RA[1] := StepNorm;
    State.RState.RA[2] := FBase;
    State.RState.RA[3] := FNew;
    State.RState.RA[4] := Lambda;
    State.RState.RA[5] := Nu;
    State.RState.RA[6] := LambdaUp;
    State.RState.RA[7] := LambdaDown;
    State.RState.RA[8] := V;
end;


(*************************************************************************
Levenberg-Marquardt algorithm results

Called after MinLMIteration returned False.

Input parameters:
    State   -   algorithm state (used by MinLMIteration).

Output parameters:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -1    incorrect parameters were specified
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  5    MaxIts steps was taken
                * Rep.IterationsCount contains iterations count
                * Rep.NFunc     - number of function calculations
                * Rep.NJac      - number of Jacobi matrix calculations
                * Rep.NGrad     - number of gradient calculations
                * Rep.NHess     - number of Hessian calculations
                * Rep.NCholesky - number of Cholesky decomposition calculations

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************)
procedure MinLMResults(const State : LMState;
     var X : TReal1DArray;
     var Rep : LMReport);
begin
    SetLength(X, State.N-1+1);
    APVMove(@X[0], 0, State.N-1, @State.X[0], 0, State.N-1);
    Rep.IterationsCount := State.RepIterationsCount;
    Rep.TerminationType := State.RepTerminationType;
    Rep.NFunc := State.RepNFunc;
    Rep.NJac := State.RepNJac;
    Rep.NGrad := State.RepNGrad;
    Rep.NHess := State.RepNHess;
    Rep.NCholesky := State.RepNCholesky;
end;


(*************************************************************************
Prepare internal structures (except for RComm).

Note: M must be zero for FGH mode, non-zero for FJ/FGJ mode.
*************************************************************************)
procedure LMPrepare(N : AlglibInteger;
     M : AlglibInteger;
     HaveGrad : Boolean;
     var State : LMState);
begin
    State.RepIterationsCount := 0;
    State.RepTerminationType := 0;
    State.RepNFunc := 0;
    State.RepNJac := 0;
    State.RepNGrad := 0;
    State.RepNHess := 0;
    State.RepNCholesky := 0;
    if (N<0) or (M<0) then
    begin
        Exit;
    end;
    if HaveGrad then
    begin
        SetLength(State.G, N-1+1);
    end;
    if M<>0 then
    begin
        SetLength(State.J, M-1+1, N-1+1);
        SetLength(State.FI, M-1+1);
        SetLength(State.H, 0+1, 0+1);
    end
    else
    begin
        SetLength(State.J, 0+1, 0+1);
        SetLength(State.FI, 0+1);
        SetLength(State.H, N-1+1, N-1+1);
    end;
    SetLength(State.X, N-1+1);
    SetLength(State.RawModel, N-1+1, N-1+1);
    SetLength(State.Model, N-1+1, N-1+1);
    SetLength(State.XBase, N-1+1);
    SetLength(State.XPrec, N-1+1);
    SetLength(State.GBase, N-1+1);
    SetLength(State.XDir, N-1+1);
    SetLength(State.Work, Max(N, M)+1);
end;


(*************************************************************************
Clears request fileds (to be sure that we don't forgot to clear something)
*************************************************************************)
procedure LMClearRequestFields(var State : LMState);
begin
    State.NeedF := False;
    State.NeedFG := False;
    State.NeedFGH := False;
    State.NeedFiJ := False;
end;


end.