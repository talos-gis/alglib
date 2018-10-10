unit testodesolverunit;
interface
uses Math, Sysutils, Ap, odesolver;

function TestODESolver(Silent : Boolean):Boolean;
function testodesolverunit_test_silent():Boolean;
function testodesolverunit_test():Boolean;

implementation

procedure Unset2D(var X : TReal2DArray);forward;
procedure Unset1D(var X : TReal1DArray);forward;
procedure UnsetRep(var Rep : ODESolverReport);forward;


(*************************************************************************
Test
*************************************************************************)
function TestODESolver(Silent : Boolean):Boolean;
var
    PassCount : AlglibInteger;
    CurErrors : Boolean;
    RKCKErrors : Boolean;
    WasErrors : Boolean;
    XTbl : TReal1DArray;
    YTbl : TReal2DArray;
    Rep : ODESolverReport;
    XG : TReal1DArray;
    Y : TReal1DArray;
    H : Double;
    Eps : Double;
    Solver : AlglibInteger;
    Pass : AlglibInteger;
    MyNFEV : AlglibInteger;
    V : Double;
    N : AlglibInteger;
    M : AlglibInteger;
    M2 : AlglibInteger;
    I : AlglibInteger;
    Err : Double;
    State : ODESolverState;
begin
    RKCKErrors := False;
    WasErrors := False;
    PassCount := 10;
    
    //
    // simple test: just A*sin(x)+B*cos(x)
    //
    Assert(PassCount>=2);
    Pass:=0;
    while Pass<=PassCount-1 do
    begin
        Solver:=0;
        while Solver<=0 do
        begin
            
            //
            // prepare
            //
            H := 1.0E-2;
            Eps := 1.0E-5;
            if Pass mod 2=0 then
            begin
                Eps := -Eps;
            end;
            SetLength(Y, 2);
            I:=0;
            while I<=1 do
            begin
                Y[I] := 2*RandomReal-1;
                Inc(I);
            end;
            M := 2+RandomInteger(10);
            SetLength(XG, M);
            XG[0] := (M-1)*RandomReal;
            I:=1;
            while I<=M-1 do
            begin
                XG[I] := XG[I-1]+RandomReal;
                Inc(I);
            end;
            V := 2*Pi/(XG[M-1]-XG[0]);
            APVMul(@XG[0], 0, M-1, V);
            if AP_FP_Greater(RandomReal,0.5) then
            begin
                APVMul(@XG[0], 0, M-1, -1);
            end;
            MyNFEV := 0;
            
            //
            // choose solver
            //
            if Solver=0 then
            begin
                ODESolverRKCK(Y, 2, XG, M, Eps, H, State);
            end;
            
            //
            // solve
            //
            while ODESolverIteration(State) do
            begin
                State.DY[0] := State.Y[1];
                State.DY[1] := -State.Y[0];
                MyNFEV := MyNFEV+1;
            end;
            ODESolverResults(State, M2, XTbl, YTbl, Rep);
            
            //
            // check results
            //
            CurErrors := False;
            if Rep.TerminationType<=0 then
            begin
                CurErrors := True;
            end
            else
            begin
                CurErrors := CurErrors or (M2<>M);
                Err := 0;
                I:=0;
                while I<=M-1 do
                begin
                    Err := Max(Err, AbsReal(YTbl[I,0]-(+Y[0]*Cos(XTbl[I]-XTbl[0])+Y[1]*Sin(XTbl[I]-XTbl[0]))));
                    Err := Max(Err, AbsReal(YTbl[I,1]-(-Y[0]*Sin(XTbl[I]-XTbl[0])+Y[1]*Cos(XTbl[I]-XTbl[0]))));
                    Inc(I);
                end;
                CurErrors := CurErrors or AP_FP_Greater(Err,10*AbsReal(Eps));
                CurErrors := CurErrors or (MyNFEV<>Rep.NFEV);
            end;
            if Solver=0 then
            begin
                RKCKErrors := RKCKErrors or CurErrors;
            end;
            Inc(Solver);
        end;
        Inc(Pass);
    end;
    
    //
    // another test:
    //
    //     y(0)   = 0
    //     dy/dx  = f(x,y)
    //     f(x,y) = 0,   x<1
    //              x-1, x>=1
    //
    // with BOTH absolute and fractional tolerances.
    // Starting from zero will be real challenge for
    // fractional tolerance.
    //
    Assert(PassCount>=2);
    Pass:=0;
    while Pass<=PassCount-1 do
    begin
        H := 1.0E-4;
        Eps := 1.0E-4;
        if Pass mod 2=0 then
        begin
            Eps := -Eps;
        end;
        SetLength(Y, 1);
        Y[0] := 0;
        M := 21;
        SetLength(XG, M);
        I:=0;
        while I<=M-1 do
        begin
            XG[I] := AP_Double(2*I)/(M-1);
            Inc(I);
        end;
        MyNFEV := 0;
        ODESolverRKCK(Y, 1, XG, M, Eps, H, State);
        while ODESolverIteration(State) do
        begin
            State.DY[0] := Max(State.X-1, 0);
            MyNFEV := MyNFEV+1;
        end;
        ODESolverResults(State, M2, XTbl, YTbl, Rep);
        if Rep.TerminationType<=0 then
        begin
            RKCKErrors := True;
        end
        else
        begin
            RKCKErrors := RKCKErrors or (M2<>M);
            Err := 0;
            I:=0;
            while I<=M-1 do
            begin
                Err := Max(Err, AbsReal(YTbl[I,0]-AP_Sqr(Max(XG[I]-1, 0))/2));
                Inc(I);
            end;
            RKCKErrors := RKCKErrors or AP_FP_Greater(Err,AbsReal(Eps));
            RKCKErrors := RKCKErrors or (MyNFEV<>Rep.NFEV);
        end;
        Inc(Pass);
    end;
    
    //
    // end
    //
    WasErrors := RKCKErrors;
    if  not Silent then
    begin
        Write(Format('TESTING ODE SOLVER'#13#10'',[]));
        Write(Format('* RK CASH-KARP:                           ',[]));
        if RKCKErrors then
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
procedure UnsetRep(var Rep : ODESolverReport);
begin
    Rep.NFEV := 0;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testodesolverunit_test_silent():Boolean;
begin
    Result := TestODESolver(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testodesolverunit_test():Boolean;
begin
    Result := TestODESolver(False);
end;


end.