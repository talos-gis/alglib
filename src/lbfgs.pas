(*************************************************************************
Copyright (c) 2007-2008, Sergey Bochkanov (ALGLIB project).

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
unit lbfgs;
interface
uses Math, Sysutils, Ap;

type
LBFGSState = record
    N : AlglibInteger;
    M : AlglibInteger;
    EpsG : Double;
    EpsF : Double;
    EpsX : Double;
    MaxIts : AlglibInteger;
    Flags : AlglibInteger;
    NFEV : AlglibInteger;
    MCStage : AlglibInteger;
    K : AlglibInteger;
    Q : AlglibInteger;
    P : AlglibInteger;
    Rho : TReal1DArray;
    Y : TReal2DArray;
    S : TReal2DArray;
    Theta : TReal1DArray;
    D : TReal1DArray;
    Stp : Double;
    WORK : TReal1DArray;
    FOld : Double;
    GammaK : Double;
    X : TReal1DArray;
    F : Double;
    G : TReal1DArray;
    XUpdated : Boolean;
    RState : RCommState;
    RepIterationsCount : AlglibInteger;
    RepNFEV : AlglibInteger;
    RepTerminationType : AlglibInteger;
    BRACKT : Boolean;
    STAGE1 : Boolean;
    INFOC : AlglibInteger;
    DG : Double;
    DGM : Double;
    DGINIT : Double;
    DGTEST : Double;
    DGX : Double;
    DGXM : Double;
    DGY : Double;
    DGYM : Double;
    FINIT : Double;
    FTEST1 : Double;
    FM : Double;
    FX : Double;
    FXM : Double;
    FY : Double;
    FYM : Double;
    STX : Double;
    STY : Double;
    STMIN : Double;
    STMAX : Double;
    WIDTH : Double;
    WIDTH1 : Double;
    XTRAPF : Double;
end;


LBFGSReport = record
    IterationsCount : AlglibInteger;
    NFEV : AlglibInteger;
    TerminationType : AlglibInteger;
end;



procedure MinLBFGS(N : AlglibInteger;
     M : AlglibInteger;
     const X : TReal1DArray;
     EpsG : Double;
     EpsF : Double;
     EpsX : Double;
     MaxIts : AlglibInteger;
     Flags : AlglibInteger;
     var State : LBFGSState);
function MinLBFGSIteration(var State : LBFGSState):Boolean;
procedure MinLBFGSResults(const State : LBFGSState;
     var X : TReal1DArray;
     var Rep : LBFGSReport);

implementation

const
    FTOL = 0.0001;
    XTOL = 100*MachineEpsilon;
    GTOL = 0.9;
    MAXFEV = 20;
    STPMIN = 1.0E-20;
    STPMAX = 1.0E20;

procedure MCSRCH(const N : AlglibInteger;
     var X : TReal1DArray;
     var F : Double;
     var G : TReal1DArray;
     const S : TReal1DArray;
     var STP : Double;
     var INFO : AlglibInteger;
     var NFEV : AlglibInteger;
     var WA : TReal1DArray;
     var State : LBFGSState;
     var Stage : AlglibInteger);forward;
procedure MCSTEP(var STX : Double;
     var FX : Double;
     var DX : Double;
     var STY : Double;
     var FY : Double;
     var DY : Double;
     var STP : Double;
     const FP : Double;
     const DP : Double;
     var BRACKT : Boolean;
     const STMIN : Double;
     const STMAX : Double;
     var INFO : AlglibInteger);forward;


(*************************************************************************
        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION

The subroutine minimizes function F(x) of N arguments by  using  a  quasi-
Newton method (LBFGS scheme) which is optimized to use  a  minimum  amount
of memory.

The subroutine generates the approximation of an inverse Hessian matrix by
using information about the last M steps of the algorithm  (instead of N).
It lessens a required amount of memory from a value  of  order  N^2  to  a
value of order 2*N*M.

Input parameters:
    N   -   problem dimension. N>0
    M   -   number of corrections in the BFGS scheme of Hessian
            approximation update. Recommended value:  3<=M<=7. The smaller
            value causes worse convergence, the bigger will  not  cause  a
            considerably better convergence, but will cause a fall in  the
            performance. M<=N.
    X   -   initial solution approximation, array[0..N-1].
    EpsG -  positive number which  defines  a  precision  of  search.  The
            subroutine finishes its work if the condition ||G|| < EpsG  is
            satisfied, where ||.|| means Euclidian norm, G - gradient, X -
            current approximation.
    EpsF -  positive number which  defines  a  precision  of  search.  The
            subroutine finishes its work if on iteration  number  k+1  the
            condition |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}    is
            satisfied.
    EpsX -  positive number which  defines  a  precision  of  search.  The
            subroutine finishes its work if on iteration number k+1    the
            condition |X(k+1)-X(k)| <= EpsX is fulfilled.
    MaxIts- maximum number of iterations. If MaxIts=0, the number of
            iterations is unlimited.
    Flags - additional settings:
            * Flags = 0     means no additional settings
            * Flags = 1     "do not allocate memory". used when solving
                            a many subsequent tasks with  same N/M  values.
                            First  call MUST  be without this flag bit set,
                            subsequent calls of MinLBFGS with same LBFGSState
                            structure can set Flags to 1.

Output parameters:
    State - structure used for reverse communication.
    
See also MinLBFGSIteration, MinLBFGSResults

NOTE:

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 14.11.2007 by Bochkanov Sergey
*************************************************************************)
procedure MinLBFGS(N : AlglibInteger;
     M : AlglibInteger;
     const X : TReal1DArray;
     EpsG : Double;
     EpsF : Double;
     EpsX : Double;
     MaxIts : AlglibInteger;
     Flags : AlglibInteger;
     var State : LBFGSState);
var
    AllocateMem : Boolean;
begin
    Assert(N>=1, 'MinLBFGS: N too small!');
    Assert(M>=1, 'MinLBFGS: M too small!');
    Assert(M<=N, 'MinLBFGS: M too large!');
    Assert(AP_FP_Greater_Eq(EpsG,0), 'MinLBFGS: negative EpsG!');
    Assert(AP_FP_Greater_Eq(EpsF,0), 'MinLBFGS: negative EpsF!');
    Assert(AP_FP_Greater_Eq(EpsX,0), 'MinLBFGS: negative EpsX!');
    Assert(MaxIts>=0, 'MinLBFGS: negative MaxIts!');
    
    //
    // Initialize
    //
    if AP_FP_Eq(EpsG,0) and AP_FP_Eq(EpsF,0) and AP_FP_Eq(EpsX,0) and (MaxIts=0) then
    begin
        EpsX := 1.0E-6;
    end;
    State.N := N;
    State.M := M;
    State.EpsG := EpsG;
    State.EpsF := EpsF;
    State.EpsX := EpsX;
    State.MaxIts := MaxIts;
    State.Flags := Flags;
    AllocateMem := Flags mod 2=0;
    Flags := Flags div 2;
    if AllocateMem then
    begin
        SetLength(State.Rho, M-1+1);
        SetLength(State.Theta, M-1+1);
        SetLength(State.Y, M-1+1, N-1+1);
        SetLength(State.S, M-1+1, N-1+1);
        SetLength(State.D, N-1+1);
        SetLength(State.X, N-1+1);
        SetLength(State.G, N-1+1);
        SetLength(State.WORK, N-1+1);
    end;
    
    //
    // Initialize Rep structure
    //
    State.XUpdated := False;
    
    //
    // Prepare first run
    //
    State.K := 0;
    APVMove(@State.X[0], 0, N-1, @X[0], 0, N-1);
    SetLength(State.RState.IA, 6+1);
    SetLength(State.RState.RA, 4+1);
    State.RState.Stage := -1;
end;


(*************************************************************************
One L-BFGS iteration

Called after initialization with MinLBFGS.
See HTML documentation for examples.

Input parameters:
    State   -   structure which stores algorithm state between calls and
                which is used for reverse communication. Must be initialized
                with MinLBFGS.

If suborutine returned False, iterative proces has converged.

If subroutine returned True, caller should calculate function value
State.F an gradient State.G[0..N-1] at State.X[0..N-1] and call
MinLBFGSIteration again.

  -- ALGLIB --
     Copyright 20.04.2009 by Bochkanov Sergey
*************************************************************************)
function MinLBFGSIteration(var State : LBFGSState):Boolean;
var
    N : AlglibInteger;
    M : AlglibInteger;
    MaxIts : AlglibInteger;
    EpsF : Double;
    EpsG : Double;
    EpsX : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    IC : AlglibInteger;
    MCINFO : AlglibInteger;
    V : Double;
    VV : Double;
label
lbl_0, lbl_2, lbl_4, lbl_1, lbl_5, lbl_3, lbl_rcomm;
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
        MaxIts := State.RState.IA[2];
        I := State.RState.IA[3];
        J := State.RState.IA[4];
        IC := State.RState.IA[5];
        MCINFO := State.RState.IA[6];
        EpsF := State.RState.RA[0];
        EpsG := State.RState.RA[1];
        EpsX := State.RState.RA[2];
        V := State.RState.RA[3];
        VV := State.RState.RA[4];
    end
    else
    begin
        N := -983;
        M := -989;
        MaxIts := -834;
        I := 900;
        J := -287;
        IC := 364;
        MCINFO := 214;
        EpsF := -338;
        EpsG := -686;
        EpsX := 912;
        V := 585;
        VV := 497;
    end;
    if State.RState.Stage=0 then
    begin
        goto lbl_0;
    end;
    if State.RState.Stage=1 then
    begin
        goto lbl_1;
    end;
    
    //
    // Routine body
    //
    
    //
    // Unload frequently used variables from State structure
    // (just for typing convinience)
    //
    N := State.N;
    M := State.M;
    EpsG := State.EpsG;
    EpsF := State.EpsF;
    EpsX := State.EpsX;
    MaxIts := State.MaxIts;
    State.RepTerminationType := 0;
    State.RepIterationsCount := 0;
    State.RepNFEV := 0;
    
    //
    // Update info
    //
    State.XUpdated := False;
    
    //
    // Calculate F/G
    //
    State.RState.Stage := 0;
    goto lbl_rcomm;
lbl_0:
    State.RepNFEV := 1;
    
    //
    // Preparations
    //
    State.FOld := State.F;
    V := APVDotProduct(@State.G[0], 0, N-1, @State.G[0], 0, N-1);
    V := Sqrt(V);
    if AP_FP_Eq(V,0) then
    begin
        State.RepTerminationType := 4;
        Result := False;
        Exit;
    end;
    State.Stp := Min(1.0/V, 1);
    APVMoveNeg(@State.D[0], 0, N-1, @State.G[0], 0, N-1);
    
    //
    // Main cycle
    //
lbl_2:
    if False then
    begin
        goto lbl_3;
    end;
    
    //
    // Main cycle: prepare to 1-D line search
    //
    State.P := State.K mod M;
    State.Q := Min(State.K, M-1);
    
    //
    // Store X[k], G[k]
    //
    APVMoveNeg(@State.S[State.P][0], 0, N-1, @State.X[0], 0, N-1);
    APVMoveNeg(@State.Y[State.P][0], 0, N-1, @State.G[0], 0, N-1);
    
    //
    // Minimize F(x+alpha*d)
    //
    State.MCStage := 0;
    if State.K<>0 then
    begin
        State.Stp := 1.0;
    end;
    MCSRCH(N, State.X, State.F, State.G, State.D, State.Stp, MCINFO, State.NFEV, State.WORK, State, State.MCStage);
lbl_4:
    if State.MCStage=0 then
    begin
        goto lbl_5;
    end;
    State.RState.Stage := 1;
    goto lbl_rcomm;
lbl_1:
    MCSRCH(N, State.X, State.F, State.G, State.D, State.Stp, MCINFO, State.NFEV, State.WORK, State, State.MCStage);
    goto lbl_4;
lbl_5:
    
    //
    // Main cycle: update information and Hessian.
    // Check stopping conditions.
    //
    State.RepNFEV := State.RepNFEV+State.NFEV;
    State.RepIterationsCount := State.RepIterationsCount+1;
    
    //
    // Calculate S[k], Y[k], Rho[k], GammaK
    //
    APVAdd(@State.S[State.P][0], 0, N-1, @State.X[0], 0, N-1);
    APVAdd(@State.Y[State.P][0], 0, N-1, @State.G[0], 0, N-1);
    
    //
    // Stopping conditions
    //
    if (State.RepIterationsCount>=MaxIts) and (MaxIts>0) then
    begin
        
        //
        // Too many iterations
        //
        State.RepTerminationType := 5;
        Result := False;
        Exit;
    end;
    V := APVDotProduct(@State.G[0], 0, N-1, @State.G[0], 0, N-1);
    if AP_FP_Less_Eq(Sqrt(V),EpsG) then
    begin
        
        //
        // Gradient is small enough
        //
        State.RepTerminationType := 4;
        Result := False;
        Exit;
    end;
    if AP_FP_Less_Eq(State.FOld-State.F,EpsF*Max(AbsReal(State.FOld), Max(AbsReal(State.F), 1.0))) then
    begin
        
        //
        // F(k+1)-F(k) is small enough
        //
        State.RepTerminationType := 1;
        Result := False;
        Exit;
    end;
    V := APVDotProduct(@State.S[State.P][0], 0, N-1, @State.S[State.P][0], 0, N-1);
    if AP_FP_Less_Eq(Sqrt(V),EpsX) then
    begin
        
        //
        // X(k+1)-X(k) is small enough
        //
        State.RepTerminationType := 2;
        Result := False;
        Exit;
    end;
    
    //
    // Calculate Rho[k], GammaK
    //
    V := APVDotProduct(@State.Y[State.P][0], 0, N-1, @State.S[State.P][0], 0, N-1);
    VV := APVDotProduct(@State.Y[State.P][0], 0, N-1, @State.Y[State.P][0], 0, N-1);
    if AP_FP_Eq(V,0) or AP_FP_Eq(VV,0) then
    begin
        
        //
        // Rounding errors make further iterations impossible.
        //
        State.RepTerminationType := -2;
        Result := False;
        Exit;
    end;
    State.Rho[State.P] := 1/V;
    State.GammaK := V/VV;
    
    //
    //  Calculate d(k+1) = -H(k+1)*g(k+1)
    //
    //  for I:=K downto K-Q do
    //      V = s(i)^T * work(iteration:I)
    //      theta(i) = V
    //      work(iteration:I+1) = work(iteration:I) - V*Rho(i)*y(i)
    //  work(last iteration) = H0*work(last iteration)
    //  for I:=K-Q to K do
    //      V = y(i)^T*work(iteration:I)
    //      work(iteration:I+1) = work(iteration:I) +(-V+theta(i))*Rho(i)*s(i)
    //
    //  NOW WORK CONTAINS d(k+1)
    //
    APVMove(@State.WORK[0], 0, N-1, @State.G[0], 0, N-1);
    I:=State.K;
    while I>=State.K-State.Q do
    begin
        IC := I mod M;
        V := APVDotProduct(@State.S[IC][0], 0, N-1, @State.Work[0], 0, N-1);
        State.Theta[IC] := V;
        VV := V*State.Rho[IC];
        APVSub(@State.Work[0], 0, N-1, @State.Y[IC][0], 0, N-1, VV);
        Dec(I);
    end;
    V := State.GammaK;
    APVMul(@State.Work[0], 0, N-1, V);
    I:=State.K-State.Q;
    while I<=State.K do
    begin
        IC := I mod M;
        V := APVDotProduct(@State.Y[IC][0], 0, N-1, @State.Work[0], 0, N-1);
        VV := State.Rho[IC]*(-V+State.Theta[IC]);
        APVAdd(@State.Work[0], 0, N-1, @State.S[IC][0], 0, N-1, VV);
        Inc(I);
    end;
    APVMoveNeg(@State.D[0], 0, N-1, @State.Work[0], 0, N-1);
    
    //
    // Next step
    //
    State.FOld := State.F;
    State.K := State.K+1;
    State.XUpdated := True;
    goto lbl_2;
lbl_3:
    Result := False;
    Exit;
    
    //
    // Saving state
    //
lbl_rcomm:
    Result := True;
    State.RState.IA[0] := N;
    State.RState.IA[1] := M;
    State.RState.IA[2] := MaxIts;
    State.RState.IA[3] := I;
    State.RState.IA[4] := J;
    State.RState.IA[5] := IC;
    State.RState.IA[6] := MCINFO;
    State.RState.RA[0] := EpsF;
    State.RState.RA[1] := EpsG;
    State.RState.RA[2] := EpsX;
    State.RState.RA[3] := V;
    State.RState.RA[4] := VV;
end;


(*************************************************************************
L-BFGS algorithm results

Called after MinLBFGSIteration returned False.

Input parameters:
    State   -   algorithm state (used by MinLBFGSIteration).

Output parameters:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -2    rounding errors prevent further improvement.
                            X contains best point found.
                    * -1    incorrect parameters were specified
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient norm is no more than EpsG
                    *  5    MaxIts steps was taken
                * Rep.IterationsCount contains iterations count
                * NFEV countains number of function calculations

  -- ALGLIB --
     Copyright 14.11.2007 by Bochkanov Sergey
*************************************************************************)
procedure MinLBFGSResults(const State : LBFGSState;
     var X : TReal1DArray;
     var Rep : LBFGSReport);
begin
    SetLength(X, State.N-1+1);
    APVMove(@X[0], 0, State.N-1, @State.X[0], 0, State.N-1);
    Rep.IterationsCount := State.RepIterationsCount;
    Rep.NFEV := State.RepNFEV;
    Rep.TerminationType := State.RepTerminationType;
end;


(*************************************************************************
THE  PURPOSE  OF  MCSRCH  IS  TO  FIND A STEP WHICH SATISFIES A SUFFICIENT
DECREASE CONDITION AND A CURVATURE CONDITION.

AT EACH STAGE THE SUBROUTINE  UPDATES  AN  INTERVAL  OF  UNCERTAINTY  WITH
ENDPOINTS  STX  AND  STY.  THE INTERVAL OF UNCERTAINTY IS INITIALLY CHOSEN
SO THAT IT CONTAINS A MINIMIZER OF THE MODIFIED FUNCTION

    F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).

IF  A STEP  IS OBTAINED FOR  WHICH THE MODIFIED FUNCTION HAS A NONPOSITIVE
FUNCTION  VALUE  AND  NONNEGATIVE  DERIVATIVE,   THEN   THE   INTERVAL  OF
UNCERTAINTY IS CHOSEN SO THAT IT CONTAINS A MINIMIZER OF F(X+STP*S).

THE  ALGORITHM  IS  DESIGNED TO FIND A STEP WHICH SATISFIES THE SUFFICIENT
DECREASE CONDITION

    F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),

AND THE CURVATURE CONDITION

    ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).

IF  FTOL  IS  LESS  THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION IS BOUNDED
BELOW,  THEN  THERE  IS  ALWAYS  A  STEP  WHICH SATISFIES BOTH CONDITIONS.
IF  NO  STEP  CAN BE FOUND  WHICH  SATISFIES  BOTH  CONDITIONS,  THEN  THE
ALGORITHM  USUALLY STOPS  WHEN  ROUNDING ERRORS  PREVENT FURTHER PROGRESS.
IN THIS CASE STP ONLY SATISFIES THE SUFFICIENT DECREASE CONDITION.

PARAMETERS DESCRIPRION

N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF VARIABLES.

X IS  AN  ARRAY  OF  LENGTH N. ON INPUT IT MUST CONTAIN THE BASE POINT FOR
THE LINE SEARCH. ON OUTPUT IT CONTAINS X+STP*S.

F IS  A  VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F AT X. ON OUTPUT
IT CONTAINS THE VALUE OF F AT X + STP*S.

G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE GRADIENT OF F AT X.
ON OUTPUT IT CONTAINS THE GRADIENT OF F AT X + STP*S.

S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE SEARCH DIRECTION.

STP  IS  A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN INITIAL ESTIMATE
OF A SATISFACTORY STEP. ON OUTPUT STP CONTAINS THE FINAL ESTIMATE.

FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. TERMINATION OCCURS WHEN THE
SUFFICIENT DECREASE CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
SATISFIED.

XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS WHEN THE RELATIVE
WIDTH OF THE INTERVAL OF UNCERTAINTY IS AT MOST XTOL.

STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH SPECIFY LOWER  AND
UPPER BOUNDS FOR THE STEP.

MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION OCCURS WHEN THE
NUMBER OF CALLS TO FCN IS AT LEAST MAXFEV BY THE END OF AN ITERATION.

INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
    INFO = 0  IMPROPER INPUT PARAMETERS.

    INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
              DIRECTIONAL DERIVATIVE CONDITION HOLD.

    INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
              IS AT MOST XTOL.

    INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.

    INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.

    INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.

    INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
              THERE MAY NOT BE A STEP WHICH SATISFIES THE
              SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
              TOLERANCES MAY BE TOO SMALL.

NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF CALLS TO FCN.

WA IS A WORK ARRAY OF LENGTH N.

ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
JORGE J. MORE', DAVID J. THUENTE
*************************************************************************)
procedure MCSRCH(const N : AlglibInteger;
     var X : TReal1DArray;
     var F : Double;
     var G : TReal1DArray;
     const S : TReal1DArray;
     var STP : Double;
     var INFO : AlglibInteger;
     var NFEV : AlglibInteger;
     var WA : TReal1DArray;
     var State : LBFGSState;
     var Stage : AlglibInteger);
var
    V : Double;
    P5 : Double;
    P66 : Double;
    ZERO : Double;
begin
    
    //
    // init
    //
    P5 := 0.5;
    P66 := 0.66;
    State.XTRAPF := 4.0;
    ZERO := 0;
    
    //
    // Main cycle
    //
    while True do
    begin
        if Stage=0 then
        begin
            
            //
            // NEXT
            //
            Stage := 2;
            Continue;
        end;
        if Stage=2 then
        begin
            State.INFOC := 1;
            INFO := 0;
            
            //
            //     CHECK THE INPUT PARAMETERS FOR ERRORS.
            //
            if (N<=0) or AP_FP_Less_Eq(STP,0) or AP_FP_Less(FTOL,0) or AP_FP_Less(GTOL,ZERO) or AP_FP_Less(XTOL,ZERO) or AP_FP_Less(STPMIN,ZERO) or AP_FP_Less(STPMAX,STPMIN) or (MAXFEV<=0) then
            begin
                Stage := 0;
                Exit;
            end;
            
            //
            //     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
            //     AND CHECK THAT S IS A DESCENT DIRECTION.
            //
            V := APVDotProduct(@G[0], 0, N-1, @S[0], 0, N-1);
            State.DGINIT := V;
            if AP_FP_Greater_Eq(State.DGINIT,0) then
            begin
                Stage := 0;
                Exit;
            end;
            
            //
            //     INITIALIZE LOCAL VARIABLES.
            //
            State.BRACKT := False;
            State.STAGE1 := True;
            NFEV := 0;
            State.FINIT := F;
            State.DGTEST := FTOL*State.DGINIT;
            State.WIDTH := STPMAX-STPMIN;
            State.WIDTH1 := State.WIDTH/P5;
            APVMove(@WA[0], 0, N-1, @X[0], 0, N-1);
            
            //
            //     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
            //     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
            //     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
            //     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
            //     THE INTERVAL OF UNCERTAINTY.
            //     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
            //     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
            //
            State.STX := 0;
            State.FX := State.FINIT;
            State.DGX := State.DGINIT;
            State.STY := 0;
            State.FY := State.FINIT;
            State.DGY := State.DGINIT;
            
            //
            // NEXT
            //
            Stage := 3;
            Continue;
        end;
        if Stage=3 then
        begin
            
            //
            //     START OF ITERATION.
            //
            //     SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
            //     TO THE PRESENT INTERVAL OF UNCERTAINTY.
            //
            if State.BRACKT then
            begin
                if AP_FP_Less(State.STX,State.STY) then
                begin
                    State.STMIN := State.STX;
                    State.STMAX := State.STY;
                end
                else
                begin
                    State.STMIN := State.STY;
                    State.STMAX := State.STX;
                end;
            end
            else
            begin
                State.STMIN := State.STX;
                State.STMAX := STP+State.XTRAPF*(STP-State.STX);
            end;
            
            //
            //        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
            //
            if AP_FP_Greater(STP,STPMAX) then
            begin
                STP := STPMAX;
            end;
            if AP_FP_Less(STP,STPMIN) then
            begin
                STP := STPMIN;
            end;
            
            //
            //        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
            //        STP BE THE LOWEST POINT OBTAINED SO FAR.
            //
            if State.BRACKT and (AP_FP_Less_Eq(STP,State.STMIN) or AP_FP_Greater_Eq(STP,State.STMAX)) or (NFEV>=MAXFEV-1) or (State.INFOC=0) or State.BRACKT and AP_FP_Less_Eq(State.STMAX-State.STMIN,XTOL*State.STMAX) then
            begin
                STP := State.STX;
            end;
            
            //
            //        EVALUATE THE FUNCTION AND GRADIENT AT STP
            //        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
            //
            APVMove(@X[0], 0, N-1, @WA[0], 0, N-1);
            APVAdd(@X[0], 0, N-1, @S[0], 0, N-1, STP);
            
            //
            // NEXT
            //
            Stage := 4;
            Exit;
        end;
        if Stage=4 then
        begin
            INFO := 0;
            NFEV := NFEV+1;
            V := APVDotProduct(@G[0], 0, N-1, @S[0], 0, N-1);
            State.DG := V;
            State.FTEST1 := State.FINIT+STP*State.DGTEST;
            
            //
            //        TEST FOR CONVERGENCE.
            //
            if State.BRACKT and (AP_FP_Less_Eq(STP,State.STMIN) or AP_FP_Greater_Eq(STP,State.STMAX)) or (State.INFOC=0) then
            begin
                INFO := 6;
            end;
            if AP_FP_Eq(STP,STPMAX) and AP_FP_Less_Eq(F,State.FTEST1) and AP_FP_Less_Eq(State.DG,State.DGTEST) then
            begin
                INFO := 5;
            end;
            if AP_FP_Eq(STP,STPMIN) and (AP_FP_Greater(F,State.FTEST1) or AP_FP_Greater_Eq(State.DG,State.DGTEST)) then
            begin
                INFO := 4;
            end;
            if NFEV>=MAXFEV then
            begin
                INFO := 3;
            end;
            if State.BRACKT and AP_FP_Less_Eq(State.STMAX-State.STMIN,XTOL*State.STMAX) then
            begin
                INFO := 2;
            end;
            if AP_FP_Less_Eq(F,State.FTEST1) and AP_FP_Less_Eq(AbsReal(State.DG),-GTOL*State.DGINIT) then
            begin
                INFO := 1;
            end;
            
            //
            //        CHECK FOR TERMINATION.
            //
            if INFO<>0 then
            begin
                Stage := 0;
                Exit;
            end;
            
            //
            //        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
            //        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
            //
            if State.STAGE1 and AP_FP_Less_Eq(F,State.FTEST1) and AP_FP_Greater_Eq(State.DG,Min(FTOL, GTOL)*State.DGINIT) then
            begin
                State.STAGE1 := False;
            end;
            
            //
            //        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
            //        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
            //        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
            //        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
            //        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
            //
            if State.STAGE1 and AP_FP_Less_Eq(F,State.FX) and AP_FP_Greater(F,State.FTEST1) then
            begin
                
                //
                //           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
                //
                State.FM := F-STP*State.DGTEST;
                State.FXM := State.FX-State.STX*State.DGTEST;
                State.FYM := State.FY-State.STY*State.DGTEST;
                State.DGM := State.DG-State.DGTEST;
                State.DGXM := State.DGX-State.DGTEST;
                State.DGYM := State.DGY-State.DGTEST;
                
                //
                //           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
                //           AND TO COMPUTE THE NEW STEP.
                //
                MCSTEP(State.STX, State.FXM, State.DGXM, State.STY, State.FYM, State.DGYM, STP, State.FM, State.DGM, State.BRACKT, State.STMIN, State.STMAX, State.INFOC);
                
                //
                //           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
                //
                State.FX := State.FXM+State.STX*State.DGTEST;
                State.FY := State.FYM+State.STY*State.DGTEST;
                State.DGX := State.DGXM+State.DGTEST;
                State.DGY := State.DGYM+State.DGTEST;
            end
            else
            begin
                
                //
                //           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
                //           AND TO COMPUTE THE NEW STEP.
                //
                MCSTEP(State.STX, State.FX, State.DGX, State.STY, State.FY, State.DGY, STP, F, State.DG, State.BRACKT, State.STMIN, State.STMAX, State.INFOC);
            end;
            
            //
            //        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
            //        INTERVAL OF UNCERTAINTY.
            //
            if State.BRACKT then
            begin
                if AP_FP_Greater_Eq(AbsReal(State.STY-State.STX),P66*State.WIDTH1) then
                begin
                    STP := State.STX+P5*(State.STY-State.STX);
                end;
                State.WIDTH1 := State.WIDTH;
                State.WIDTH := AbsReal(State.STY-State.STX);
            end;
            
            //
            //  NEXT.
            //
            Stage := 3;
            Continue;
        end;
    end;
end;


procedure MCSTEP(var STX : Double;
     var FX : Double;
     var DX : Double;
     var STY : Double;
     var FY : Double;
     var DY : Double;
     var STP : Double;
     const FP : Double;
     const DP : Double;
     var BRACKT : Boolean;
     const STMIN : Double;
     const STMAX : Double;
     var INFO : AlglibInteger);
var
    BOUND : Boolean;
    GAMMA : Double;
    P : Double;
    Q : Double;
    R : Double;
    S : Double;
    SGND : Double;
    STPC : Double;
    STPF : Double;
    STPQ : Double;
    THETA : Double;
begin
    INFO := 0;
    
    //
    //     CHECK THE INPUT PARAMETERS FOR ERRORS.
    //
    if BRACKT and (AP_FP_Less_Eq(STP,Min(STX, STY)) or AP_FP_Greater_Eq(STP,Max(STX, STY))) or AP_FP_Greater_Eq(DX*(STP-STX),0) or AP_FP_Less(STMAX,STMIN) then
    begin
        Exit;
    end;
    
    //
    //     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
    //
    SGND := DP*(DX/AbsReal(DX));
    
    //
    //     FIRST CASE. A HIGHER FUNCTION VALUE.
    //     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
    //     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
    //     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
    //
    if AP_FP_Greater(FP,FX) then
    begin
        INFO := 1;
        BOUND := True;
        THETA := 3*(FX-FP)/(STP-STX)+DX+DP;
        S := Max(AbsReal(THETA), Max(AbsReal(DX), AbsReal(DP)));
        GAMMA := S*Sqrt(AP_Sqr(THETA/S)-DX/S*(DP/S));
        if AP_FP_Less(STP,STX) then
        begin
            GAMMA := -GAMMA;
        end;
        P := GAMMA-DX+THETA;
        Q := GAMMA-DX+GAMMA+DP;
        R := P/Q;
        STPC := STX+R*(STP-STX);
        STPQ := STX+DX/((FX-FP)/(STP-STX)+DX)/2*(STP-STX);
        if AP_FP_Less(AbsReal(STPC-STX),AbsReal(STPQ-STX)) then
        begin
            STPF := STPC;
        end
        else
        begin
            STPF := STPC+(STPQ-STPC)/2;
        end;
        BRACKT := True;
    end
    else
    begin
        if AP_FP_Less(SGND,0) then
        begin
            
            //
            //     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
            //     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
            //     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
            //     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
            //
            INFO := 2;
            BOUND := False;
            THETA := 3*(FX-FP)/(STP-STX)+DX+DP;
            S := Max(AbsReal(THETA), Max(AbsReal(DX), AbsReal(DP)));
            GAMMA := S*SQRT(AP_Sqr(THETA/S)-DX/S*(DP/S));
            if AP_FP_Greater(STP,STX) then
            begin
                GAMMA := -GAMMA;
            end;
            P := GAMMA-DP+THETA;
            Q := GAMMA-DP+GAMMA+DX;
            R := P/Q;
            STPC := STP+R*(STX-STP);
            STPQ := STP+DP/(DP-DX)*(STX-STP);
            if AP_FP_Greater(AbsReal(STPC-STP),AbsReal(STPQ-STP)) then
            begin
                STPF := STPC;
            end
            else
            begin
                STPF := STPQ;
            end;
            BRACKT := True;
        end
        else
        begin
            if AP_FP_Less(AbsReal(DP),AbsReal(DX)) then
            begin
                
                //
                //     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
                //     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
                //     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
                //     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
                //     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
                //     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
                //     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
                //     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
                //
                INFO := 3;
                BOUND := True;
                THETA := 3*(FX-FP)/(STP-STX)+DX+DP;
                S := Max(AbsReal(THETA), Max(AbsReal(DX), AbsReal(DP)));
                
                //
                //        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
                //        TO INFINITY IN THE DIRECTION OF THE STEP.
                //
                GAMMA := S*SQRT(Max(0, AP_Sqr(THETA/S)-DX/S*(DP/S)));
                if AP_FP_Greater(STP,STX) then
                begin
                    GAMMA := -GAMMA;
                end;
                P := GAMMA-DP+THETA;
                Q := GAMMA+(DX-DP)+GAMMA;
                R := P/Q;
                if AP_FP_Less(R,0) and AP_FP_Neq(GAMMA,0) then
                begin
                    STPC := STP+R*(STX-STP);
                end
                else
                begin
                    if AP_FP_Greater(STP,STX) then
                    begin
                        STPC := STMAX;
                    end
                    else
                    begin
                        STPC := STMIN;
                    end;
                end;
                STPQ := STP+DP/(DP-DX)*(STX-STP);
                if BRACKT then
                begin
                    if AP_FP_Less(AbsReal(STP-STPC),AbsReal(STP-STPQ)) then
                    begin
                        STPF := STPC;
                    end
                    else
                    begin
                        STPF := STPQ;
                    end;
                end
                else
                begin
                    if AP_FP_Greater(AbsReal(STP-STPC),AbsReal(STP-STPQ)) then
                    begin
                        STPF := STPC;
                    end
                    else
                    begin
                        STPF := STPQ;
                    end;
                end;
            end
            else
            begin
                
                //
                //     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
                //     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
                //     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
                //     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
                //
                INFO := 4;
                BOUND := False;
                if BRACKT then
                begin
                    THETA := 3*(FP-FY)/(STY-STP)+DY+DP;
                    S := Max(ABSReal(THETA), Max(ABSReal(DY), ABSReal(DP)));
                    GAMMA := S*SQRT(AP_Sqr(THETA/S)-DY/S*(DP/S));
                    if AP_FP_Greater(STP,STY) then
                    begin
                        GAMMA := -GAMMA;
                    end;
                    P := GAMMA-DP+THETA;
                    Q := GAMMA-DP+GAMMA+DY;
                    R := P/Q;
                    STPC := STP+R*(STY-STP);
                    STPF := STPC;
                end
                else
                begin
                    if AP_FP_Greater(STP,STX) then
                    begin
                        STPF := STMAX;
                    end
                    else
                    begin
                        STPF := STMIN;
                    end;
                end;
            end;
        end;
    end;
    
    //
    //     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
    //     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
    //
    if AP_FP_Greater(FP,FX) then
    begin
        STY := STP;
        FY := FP;
        DY := DP;
    end
    else
    begin
        if AP_FP_Less(SGND,0.0) then
        begin
            STY := STX;
            FY := FX;
            DY := DX;
        end;
        STX := STP;
        FX := FP;
        DX := DP;
    end;
    
    //
    //     COMPUTE THE NEW STEP AND SAFEGUARD IT.
    //
    STPF := Min(STMAX, STPF);
    STPF := Max(STMIN, STPF);
    STP := STPF;
    if BRACKT and BOUND then
    begin
        if AP_FP_Greater(STY,STX) then
        begin
            STP := Min(STX+0.66*(STY-STX), STP);
        end
        else
        begin
            STP := Max(STX+0.66*(STY-STX), STP);
        end;
    end;
end;


end.