unit testgkq;
interface
uses Math, Sysutils, Ap, tsort, hblas, reflections, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, hsschur, evd, gammafunc, gq, gkq;

function TestGKQunit(Silent : Boolean):Boolean;
function testgkq_test_silent():Boolean;
function testgkq_test():Boolean;

implementation

function MapKind(K : AlglibInteger):Double;forward;


(*************************************************************************
Test
*************************************************************************)
function TestGKQunit(Silent : Boolean):Boolean;
var
    PKind : AlglibInteger;
    ErrTol : Double;
    Eps : Double;
    NonStrictErrTol : Double;
    N : AlglibInteger;
    I : AlglibInteger;
    K : AlglibInteger;
    Info : AlglibInteger;
    Err : Double;
    AKind : AlglibInteger;
    BKind : AlglibInteger;
    AlphaC : Double;
    BetaC : Double;
    X1 : TReal1DArray;
    WG1 : TReal1DArray;
    WK1 : TReal1DArray;
    X2 : TReal1DArray;
    WG2 : TReal1DArray;
    WK2 : TReal1DArray;
    Info1 : AlglibInteger;
    Info2 : AlglibInteger;
    SuccessAtLeastOnce : Boolean;
    InTblErrors : Boolean;
    VsTblErrors : Boolean;
    GenErrors : Boolean;
    WasErrors : Boolean;
begin
    InTblErrors := False;
    VsTblErrors := False;
    GenErrors := False;
    WasErrors := False;
    ErrTol := 10000*MachineEpsilon;
    NonStrictErrTol := 1000*ErrTol;
    
    //
    // test recurrence-based Legendre nodes against the precalculated table
    //
    PKind:=0;
    while PKind<=5 do
    begin
        N := 0;
        if PKind=0 then
        begin
            N := 15;
        end;
        if PKind=1 then
        begin
            N := 21;
        end;
        if PKind=2 then
        begin
            N := 31;
        end;
        if PKind=3 then
        begin
            N := 41;
        end;
        if PKind=4 then
        begin
            N := 51;
        end;
        if PKind=5 then
        begin
            N := 61;
        end;
        GKQLegendreCalc(N, Info, X1, WK1, WG1);
        GKQLegendreTbl(N, X2, WK2, WG2, Eps);
        if Info<=0 then
        begin
            GenErrors := True;
            Break;
        end;
        I:=0;
        while I<=N-1 do
        begin
            VsTblErrors := VsTblErrors or AP_FP_Greater(AbsReal(X1[I]-X2[I]),ErrTol);
            VsTblErrors := VsTblErrors or AP_FP_Greater(AbsReal(WK1[I]-WK2[I]),ErrTol);
            VsTblErrors := VsTblErrors or AP_FP_Greater(AbsReal(WG1[I]-WG2[I]),ErrTol);
            Inc(I);
        end;
        Inc(PKind);
    end;
    
    //
    // Test recurrence-baced Gauss-Kronrod nodes against Gauss-only nodes
    // calculated with subroutines from GQ unit.
    //
    K:=1;
    while K<=30 do
    begin
        N := 2*K+1;
        
        //
        // Gauss-Legendre
        //
        Err := 0;
        GKQGenerateGaussLegendre(N, Info1, X1, WK1, WG1);
        GQGenerateGaussLegendre(K, Info2, X2, WG2);
        if (Info1>0) and (Info2>0) then
        begin
            I:=0;
            while I<=K-1 do
            begin
                Err := Max(Err, AbsReal(X1[2*I+1]-X2[I]));
                Err := Max(Err, AbsReal(WG1[2*I+1]-WG2[I]));
                Inc(I);
            end;
        end
        else
        begin
            GenErrors := True;
        end;
        GenErrors := GenErrors or AP_FP_Greater(Err,ErrTol);
        Inc(K);
    end;
    K:=1;
    while K<=15 do
    begin
        N := 2*K+1;
        
        //
        // Gauss-Jacobi
        //
        SuccessAtLeastOnce := False;
        Err := 0;
        AKind:=0;
        while AKind<=9 do
        begin
            BKind:=0;
            while BKind<=9 do
            begin
                AlphaC := MapKind(AKind);
                BetaC := MapKind(BKind);
                GKQGenerateGaussJacobi(N, AlphaC, BetaC, Info1, X1, WK1, WG1);
                GQGenerateGaussJacobi(K, AlphaC, BetaC, Info2, X2, WG2);
                if (Info1>0) and (Info2>0) then
                begin
                    SuccessAtLeastOnce := True;
                    I:=0;
                    while I<=K-1 do
                    begin
                        Err := Max(Err, AbsReal(X1[2*I+1]-X2[I]));
                        Err := Max(Err, AbsReal(WG1[2*I+1]-WG2[I]));
                        Inc(I);
                    end;
                end
                else
                begin
                    GenErrors := GenErrors or (Info1<>-5);
                end;
                Inc(BKind);
            end;
            Inc(AKind);
        end;
        GenErrors := GenErrors or AP_FP_Greater(Err,ErrTol) or  not SuccessAtLeastOnce;
        Inc(K);
    end;
    
    //
    // end
    //
    WasErrors := InTblErrors or VsTblErrors or GenErrors;
    if  not Silent then
    begin
        Write(Format('TESTING GAUSS-KRONROD QUADRATURES'#13#10'',[]));
        Write(Format('FINAL RESULT:                             ',[]));
        if WasErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* PRE-CALCULATED TABLE:                   ',[]));
        if InTblErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* CALCULATED AGAINST THE TABLE:           ',[]));
        if VsTblErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* GENERAL PROPERTIES:                     ',[]));
        if GenErrors then
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
Maps:
    0   =>  -0.9
    1   =>  -0.5
    2   =>  -0.1
    3   =>   0.0
    4   =>  +0.1
    5   =>  +0.5
    6   =>  +0.9
    7   =>  +1.0
    8   =>  +1.5
    9   =>  +2.0
*************************************************************************)
function MapKind(K : AlglibInteger):Double;
begin
    Result := 0;
    if K=0 then
    begin
        Result := -0.9;
    end;
    if K=1 then
    begin
        Result := -0.5;
    end;
    if K=2 then
    begin
        Result := -0.1;
    end;
    if K=3 then
    begin
        Result := 0.0;
    end;
    if K=4 then
    begin
        Result := +0.1;
    end;
    if K=5 then
    begin
        Result := +0.5;
    end;
    if K=6 then
    begin
        Result := +0.9;
    end;
    if K=7 then
    begin
        Result := +1.0;
    end;
    if K=8 then
    begin
        Result := +1.5;
    end;
    if K=9 then
    begin
        Result := +2.0;
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testgkq_test_silent():Boolean;
begin
    Result := TestGKQunit(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testgkq_test():Boolean;
begin
    Result := TestGKQunit(False);
end;


end.