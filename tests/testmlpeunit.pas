unit testmlpeunit;
interface
uses Math, Sysutils, Ap, mlpbase, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, lbfgs, hblas, sblas, ortfac, blas, rotations, bdsvd, svd, xblas, densesolver, mlptrain, tsort, descriptivestatistics, bdss, mlpe;

function TestMLPE(Silent : Boolean):Boolean;
function testmlpeunit_test_silent():Boolean;
function testmlpeunit_test():Boolean;

implementation

procedure CreateEnsemble(var Ensemble : MLPEnsemble;
     NKind : AlglibInteger;
     A1 : Double;
     A2 : Double;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     EC : AlglibInteger);forward;
procedure UnsetEnsemble(var Ensemble : MLPEnsemble);forward;
procedure TestInformational(NKind : AlglibInteger;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     EC : AlglibInteger;
     PassCount : AlglibInteger;
     var Err : Boolean);forward;
procedure TestProcessing(NKind : AlglibInteger;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     EC : AlglibInteger;
     PassCount : AlglibInteger;
     var Err : Boolean);forward;


function TestMLPE(Silent : Boolean):Boolean;
var
    WasErrors : Boolean;
    PassCount : AlglibInteger;
    MaxN : AlglibInteger;
    MaxHid : AlglibInteger;
    NF : AlglibInteger;
    NHid : AlglibInteger;
    NL : AlglibInteger;
    NHid1 : AlglibInteger;
    NHid2 : AlglibInteger;
    EC : AlglibInteger;
    NKind : AlglibInteger;
    AlgType : AlglibInteger;
    TaskType : AlglibInteger;
    Pass : AlglibInteger;
    Ensemble : MLPEnsemble;
    Rep : MLPReport;
    OOBRep : MLPCVReport;
    XY : TReal2DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    NIn : AlglibInteger;
    NOut : AlglibInteger;
    NPoints : AlglibInteger;
    E : Double;
    Info : AlglibInteger;
    NLess : AlglibInteger;
    NAll : AlglibInteger;
    NClasses : AlglibInteger;
    AllSame : Boolean;
    InfErrors : Boolean;
    ProcErrors : Boolean;
    TrnErrors : Boolean;
begin
    WasErrors := False;
    InfErrors := False;
    ProcErrors := False;
    TrnErrors := False;
    PassCount := 10;
    MaxN := 4;
    MaxHid := 4;
    
    //
    // General MLP ensembles tests
    //
    NF:=1;
    while NF<=MaxN do
    begin
        NL:=1;
        while NL<=MaxN do
        begin
            NHid1:=0;
            while NHid1<=MaxHid do
            begin
                NHid2:=0;
                while NHid2<=0 do
                begin
                    NKind:=0;
                    while NKind<=3 do
                    begin
                        EC:=1;
                        while EC<=3 do
                        begin
                            
                            //
                            //  Skip meaningless parameters combinations
                            //
                            if (NKind=1) and (NL<2) then
                            begin
                                Inc(EC);
                                Continue;
                            end;
                            if (NHid1=0) and (NHid2<>0) then
                            begin
                                Inc(EC);
                                Continue;
                            end;
                            
                            //
                            // Tests
                            //
                            TestInformational(NKind, NF, NHid1, NHid2, NL, EC, PassCount, InfErrors);
                            TestProcessing(NKind, NF, NHid1, NHid2, NL, EC, PassCount, ProcErrors);
                            Inc(EC);
                        end;
                        Inc(NKind);
                    end;
                    Inc(NHid2);
                end;
                Inc(NHid1);
            end;
            Inc(NL);
        end;
        Inc(NF);
    end;
    
    //
    // network training must reduce error
    // test on random regression task
    //
    NIn := 3;
    NOut := 2;
    NHid := 5;
    NPoints := 100;
    NLess := 0;
    NAll := 0;
    Pass:=1;
    while Pass<=10 do
    begin
        AlgType:=0;
        while AlgType<=1 do
        begin
            TaskType:=0;
            while TaskType<=1 do
            begin
                if TaskType=0 then
                begin
                    SetLength(XY, NPoints-1+1, NIn+NOut-1+1);
                    I:=0;
                    while I<=NPoints-1 do
                    begin
                        J:=0;
                        while J<=NIn+NOut-1 do
                        begin
                            XY[I,J] := 2*RandomReal-1;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    MLPECreate1(NIn, NHid, NOut, 1+RandomInteger(3), Ensemble);
                end
                else
                begin
                    SetLength(XY, NPoints-1+1, NIn+1);
                    NClasses := 2+RandomInteger(2);
                    I:=0;
                    while I<=NPoints-1 do
                    begin
                        J:=0;
                        while J<=NIn-1 do
                        begin
                            XY[I,J] := 2*RandomReal-1;
                            Inc(J);
                        end;
                        XY[I,NIn] := RandomInteger(NClasses);
                        Inc(I);
                    end;
                    MLPECreateC1(NIn, NHid, NClasses, 1+RandomInteger(3), Ensemble);
                end;
                E := MLPERMSError(Ensemble, XY, NPoints);
                if AlgType=0 then
                begin
                    MLPEBaggingLM(Ensemble, XY, NPoints, 0.001, 1, Info, Rep, OOBRep);
                end
                else
                begin
                    MLPEBaggingLBFGS(Ensemble, XY, NPoints, 0.001, 1, 0.01, 0, Info, Rep, OOBRep);
                end;
                if Info<0 then
                begin
                    TrnErrors := True;
                end
                else
                begin
                    if AP_FP_Less(MLPERMSError(Ensemble, XY, NPoints),E) then
                    begin
                        NLess := NLess+1;
                    end;
                end;
                NAll := NAll+1;
                Inc(TaskType);
            end;
            Inc(AlgType);
        end;
        Inc(Pass);
    end;
    TrnErrors := TrnErrors or AP_FP_Greater(NAll-NLess,0.3*NAll);
    
    //
    // Final report
    //
    WasErrors := InfErrors or ProcErrors or TrnErrors;
    if  not Silent then
    begin
        Write(Format('MLP ENSEMBLE TEST'#13#10'',[]));
        Write(Format('INFORMATIONAL FUNCTIONS:                 ',[]));
        if  not InfErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('BASIC PROCESSING:                        ',[]));
        if  not ProcErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('TRAINING:                                ',[]));
        if  not TrnErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        if WasErrors then
        begin
            Write(Format('TEST SUMMARY: FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('TEST SUMMARY: PASSED'#13#10'',[]));
        end;
        Write(Format(''#13#10''#13#10'',[]));
    end;
    Result :=  not WasErrors;
end;


(*************************************************************************
Network creation
*************************************************************************)
procedure CreateEnsemble(var Ensemble : MLPEnsemble;
     NKind : AlglibInteger;
     A1 : Double;
     A2 : Double;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     EC : AlglibInteger);
begin
    Assert((NIn>0) and (NHid1>=0) and (NHid2>=0) and (NOut>0), 'CreateNetwork error');
    Assert((NHid1<>0) or (NHid2=0), 'CreateNetwork error');
    Assert((NKind<>1) or (NOut>=2), 'CreateNetwork error');
    if NHid1=0 then
    begin
        
        //
        // No hidden layers
        //
        if NKind=0 then
        begin
            MLPECreate0(NIn, NOut, EC, Ensemble);
        end
        else
        begin
            if NKind=1 then
            begin
                MLPECreateC0(NIn, NOut, EC, Ensemble);
            end
            else
            begin
                if NKind=2 then
                begin
                    MLPECreateB0(NIn, NOut, A1, A2, EC, Ensemble);
                end
                else
                begin
                    if NKind=3 then
                    begin
                        MLPECreateR0(NIn, NOut, A1, A2, EC, Ensemble);
                    end;
                end;
            end;
        end;
        Exit;
    end;
    if NHid2=0 then
    begin
        
        //
        // One hidden layer
        //
        if NKind=0 then
        begin
            MLPECreate1(NIn, NHid1, NOut, EC, Ensemble);
        end
        else
        begin
            if NKind=1 then
            begin
                MLPECreateC1(NIn, NHid1, NOut, EC, Ensemble);
            end
            else
            begin
                if NKind=2 then
                begin
                    MLPECreateB1(NIn, NHid1, NOut, A1, A2, EC, Ensemble);
                end
                else
                begin
                    if NKind=3 then
                    begin
                        MLPECreateR1(NIn, NHid1, NOut, A1, A2, EC, Ensemble);
                    end;
                end;
            end;
        end;
        Exit;
    end;
    
    //
    // Two hidden layers
    //
    if NKind=0 then
    begin
        MLPECreate2(NIn, NHid1, NHid2, NOut, EC, Ensemble);
    end
    else
    begin
        if NKind=1 then
        begin
            MLPECreateC2(NIn, NHid1, NHid2, NOut, EC, Ensemble);
        end
        else
        begin
            if NKind=2 then
            begin
                MLPECreateB2(NIn, NHid1, NHid2, NOut, A1, A2, EC, Ensemble);
            end
            else
            begin
                if NKind=3 then
                begin
                    MLPECreateR2(NIn, NHid1, NHid2, NOut, A1, A2, EC, Ensemble);
                end;
            end;
        end;
    end;
end;


(*************************************************************************
Unsets network (initialize it to smallest network possible
*************************************************************************)
procedure UnsetEnsemble(var Ensemble : MLPEnsemble);
begin
    MLPECreate0(1, 1, 1, Ensemble);
end;


(*************************************************************************
Iformational functions test
*************************************************************************)
procedure TestInformational(NKind : AlglibInteger;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     EC : AlglibInteger;
     PassCount : AlglibInteger;
     var Err : Boolean);
var
    Ensemble : MLPEnsemble;
    N1 : AlglibInteger;
    N2 : AlglibInteger;
begin
    CreateEnsemble(Ensemble, NKind, -1.0, 1.0, NIn, NHid1, NHid2, NOut, EC);
    MLPEProperties(Ensemble, N1, N2);
    Err := Err or (N1<>NIn) or (N2<>NOut);
end;


(*************************************************************************
Processing functions test
*************************************************************************)
procedure TestProcessing(NKind : AlglibInteger;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     EC : AlglibInteger;
     PassCount : AlglibInteger;
     var Err : Boolean);
var
    Ensemble : MLPEnsemble;
    Ensemble2 : MLPEnsemble;
    ZeroNet : Boolean;
    A1 : Double;
    A2 : Double;
    Pass : AlglibInteger;
    I : AlglibInteger;
    AllSame : Boolean;
    RLen : AlglibInteger;
    X1 : TReal1DArray;
    X2 : TReal1DArray;
    Y1 : TReal1DArray;
    Y2 : TReal1DArray;
    RA : TReal1DArray;
    RA2 : TReal1DArray;
    V : Double;
begin
    
    //
    // Prepare network
    //
    A1 := 0;
    A2 := 0;
    if NKind=2 then
    begin
        A1 := 1000*RandomReal-500;
        A2 := 2*RandomReal-1;
    end;
    if NKind=3 then
    begin
        A1 := 1000*RandomReal-500;
        A2 := A1+(2*RandomInteger(2)-1)*(0.1+0.9*RandomReal);
    end;
    
    //
    // Initialize arrays
    //
    SetLength(X1, NIn-1+1);
    SetLength(X2, NIn-1+1);
    SetLength(Y1, NOut-1+1);
    SetLength(Y2, NOut-1+1);
    
    //
    // Main cycle
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        CreateEnsemble(Ensemble, NKind, A1, A2, NIn, NHid1, NHid2, NOut, EC);
        
        //
        // Same inputs leads to same outputs
        //
        I:=0;
        while I<=NIn-1 do
        begin
            X1[I] := 2*RandomReal-1;
            X2[I] := X1[I];
            Inc(I);
        end;
        I:=0;
        while I<=NOut-1 do
        begin
            Y1[I] := 2*RandomReal-1;
            Y2[I] := 2*RandomReal-1;
            Inc(I);
        end;
        MLPEProcess(Ensemble, X1, Y1);
        MLPEProcess(Ensemble, X2, Y2);
        AllSame := True;
        I:=0;
        while I<=NOut-1 do
        begin
            AllSame := AllSame and AP_FP_Eq(Y1[I],Y2[I]);
            Inc(I);
        end;
        Err := Err or  not AllSame;
        
        //
        // Same inputs on original network leads to same outputs
        // on copy created using MLPCopy
        //
        UnsetEnsemble(Ensemble2);
        MLPECopy(Ensemble, Ensemble2);
        I:=0;
        while I<=NIn-1 do
        begin
            X1[I] := 2*RandomReal-1;
            X2[I] := X1[I];
            Inc(I);
        end;
        I:=0;
        while I<=NOut-1 do
        begin
            Y1[I] := 2*RandomReal-1;
            Y2[I] := 2*RandomReal-1;
            Inc(I);
        end;
        MLPEProcess(Ensemble, X1, Y1);
        MLPEProcess(Ensemble2, X2, Y2);
        AllSame := True;
        I:=0;
        while I<=NOut-1 do
        begin
            AllSame := AllSame and AP_FP_Eq(Y1[I],Y2[I]);
            Inc(I);
        end;
        Err := Err or  not AllSame;
        
        //
        // Same inputs on original network leads to same outputs
        // on copy created using MLPSerialize
        //
        UnsetEnsemble(Ensemble2);
        MLPESerialize(Ensemble, RA, RLen);
        SetLength(RA2, RLen-1+1);
        I:=0;
        while I<=RLen-1 do
        begin
            RA2[I] := RA[I];
            Inc(I);
        end;
        MLPEUnserialize(RA2, Ensemble2);
        I:=0;
        while I<=NIn-1 do
        begin
            X1[I] := 2*RandomReal-1;
            X2[I] := X1[I];
            Inc(I);
        end;
        I:=0;
        while I<=NOut-1 do
        begin
            Y1[I] := 2*RandomReal-1;
            Y2[I] := 2*RandomReal-1;
            Inc(I);
        end;
        MLPEProcess(Ensemble, X1, Y1);
        MLPEProcess(Ensemble2, X2, Y2);
        AllSame := True;
        I:=0;
        while I<=NOut-1 do
        begin
            AllSame := AllSame and AP_FP_Eq(Y1[I],Y2[I]);
            Inc(I);
        end;
        Err := Err or  not AllSame;
        
        //
        // Different inputs leads to different outputs (non-zero network)
        //
        I:=0;
        while I<=NIn-1 do
        begin
            X1[I] := 2*RandomReal-1;
            X2[I] := 2*RandomReal-1;
            Inc(I);
        end;
        I:=0;
        while I<=NOut-1 do
        begin
            Y1[I] := 2*RandomReal-1;
            Y2[I] := Y1[I];
            Inc(I);
        end;
        MLPEProcess(Ensemble, X1, Y1);
        MLPEProcess(Ensemble, X2, Y2);
        AllSame := True;
        I:=0;
        while I<=NOut-1 do
        begin
            AllSame := AllSame and AP_FP_Eq(Y1[I],Y2[I]);
            Inc(I);
        end;
        Err := Err or AllSame;
        
        //
        // Randomization changes outputs (when inputs are unchanged, non-zero network)
        //
        I:=0;
        while I<=NIn-1 do
        begin
            X1[I] := 2*RandomReal-1;
            X2[I] := 2*RandomReal-1;
            Inc(I);
        end;
        I:=0;
        while I<=NOut-1 do
        begin
            Y1[I] := 2*RandomReal-1;
            Y2[I] := Y1[I];
            Inc(I);
        end;
        MLPECopy(Ensemble, Ensemble2);
        MLPERandomize(Ensemble2);
        MLPEProcess(Ensemble, X1, Y1);
        MLPEProcess(Ensemble2, X1, Y2);
        AllSame := True;
        I:=0;
        while I<=NOut-1 do
        begin
            AllSame := AllSame and AP_FP_Eq(Y1[I],Y2[I]);
            Inc(I);
        end;
        Err := Err or AllSame;
        
        //
        // Normalization properties
        //
        if NKind=1 then
        begin
            
            //
            // Classifier network outputs are normalized
            //
            I:=0;
            while I<=NIn-1 do
            begin
                X1[I] := 2*RandomReal-1;
                Inc(I);
            end;
            MLPEProcess(Ensemble, X1, Y1);
            V := 0;
            I:=0;
            while I<=NOut-1 do
            begin
                V := V+Y1[I];
                Err := Err or AP_FP_Less(Y1[I],0);
                Inc(I);
            end;
            Err := Err or AP_FP_Greater(AbsReal(V-1),1000*MachineEpsilon);
        end;
        if NKind=2 then
        begin
            
            //
            // B-type network outputs are bounded from above/below
            //
            I:=0;
            while I<=NIn-1 do
            begin
                X1[I] := 2*RandomReal-1;
                Inc(I);
            end;
            MLPEProcess(Ensemble, X1, Y1);
            I:=0;
            while I<=NOut-1 do
            begin
                if AP_FP_Greater_Eq(A2,0) then
                begin
                    Err := Err or AP_FP_Less(Y1[I],A1);
                end
                else
                begin
                    Err := Err or AP_FP_Greater(Y1[I],A1);
                end;
                Inc(I);
            end;
        end;
        if NKind=3 then
        begin
            
            //
            // R-type network outputs are within [A1,A2] (or [A2,A1])
            //
            I:=0;
            while I<=NIn-1 do
            begin
                X1[I] := 2*RandomReal-1;
                Inc(I);
            end;
            MLPEProcess(Ensemble, X1, Y1);
            I:=0;
            while I<=NOut-1 do
            begin
                Err := Err or AP_FP_Less(Y1[I],Min(A1, A2)) or AP_FP_Greater(Y1[I],Max(A1, A2));
                Inc(I);
            end;
        end;
        Inc(Pass);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testmlpeunit_test_silent():Boolean;
begin
    Result := TestMLPE(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testmlpeunit_test():Boolean;
begin
    Result := TestMLPE(False);
end;


end.