unit testmlpunit;
interface
uses Math, Sysutils, Ap, mlpbase, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, lbfgs, hblas, sblas, ortfac, blas, rotations, bdsvd, svd, xblas, densesolver, mlptrain;

function TestMLP(Silent : Boolean):Boolean;
function testmlpunit_test_silent():Boolean;
function testmlpunit_test():Boolean;

implementation

procedure CreateNetwork(var Network : MultiLayerPerceptron;
     NKind : AlglibInteger;
     A1 : Double;
     A2 : Double;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger);forward;
procedure UnsetNetwork(var Network : MultiLayerPerceptron);forward;
procedure TestInformational(NKind : AlglibInteger;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     PassCount : AlglibInteger;
     var Err : Boolean);forward;
procedure TestProcessing(NKind : AlglibInteger;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     PassCount : AlglibInteger;
     var Err : Boolean);forward;
procedure TestGradient(NKind : AlglibInteger;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     PassCount : AlglibInteger;
     var Err : Boolean);forward;
procedure TestHessian(NKind : AlglibInteger;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     PassCount : AlglibInteger;
     var Err : Boolean);forward;


function TestMLP(Silent : Boolean):Boolean;
var
    WasErrors : Boolean;
    PassCount : AlglibInteger;
    MaxN : AlglibInteger;
    MaxHid : AlglibInteger;
    Info : AlglibInteger;
    NF : AlglibInteger;
    NHid : AlglibInteger;
    NL : AlglibInteger;
    NHid1 : AlglibInteger;
    NHid2 : AlglibInteger;
    NKind : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    Network : MultiLayerPerceptron;
    Network2 : MultiLayerPerceptron;
    Rep : MLPReport;
    CVRep : MLPCVReport;
    NCount : AlglibInteger;
    XY : TReal2DArray;
    ValXY : TReal2DArray;
    SSize : AlglibInteger;
    ValSize : AlglibInteger;
    AllSame : Boolean;
    InfErrors : Boolean;
    ProcErrors : Boolean;
    GradErrors : Boolean;
    HessErrors : Boolean;
    TrnErrors : Boolean;
begin
    WasErrors := False;
    InfErrors := False;
    ProcErrors := False;
    GradErrors := False;
    HessErrors := False;
    TrnErrors := False;
    PassCount := 10;
    MaxN := 4;
    MaxHid := 4;
    
    //
    // General multilayer network tests
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
                        
                        //
                        //  Skip meaningless parameters combinations
                        //
                        if (NKind=1) and (NL<2) then
                        begin
                            Inc(NKind);
                            Continue;
                        end;
                        if (NHid1=0) and (NHid2<>0) then
                        begin
                            Inc(NKind);
                            Continue;
                        end;
                        
                        //
                        // Tests
                        //
                        TestInformational(NKind, NF, NHid1, NHid2, NL, PassCount, InfErrors);
                        TestProcessing(NKind, NF, NHid1, NHid2, NL, PassCount, ProcErrors);
                        TestGradient(NKind, NF, NHid1, NHid2, NL, PassCount, GradErrors);
                        TestHessian(NKind, NF, NHid1, NHid2, NL, PassCount, HessErrors);
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
    // Test network training on simple XOR problem
    //
    SetLength(XY, 3+1, 2+1);
    XY[0,0] := -1;
    XY[0,1] := -1;
    XY[0,2] := -1;
    XY[1,0] := +1;
    XY[1,1] := -1;
    XY[1,2] := +1;
    XY[2,0] := -1;
    XY[2,1] := +1;
    XY[2,2] := +1;
    XY[3,0] := +1;
    XY[3,1] := +1;
    XY[3,2] := -1;
    MLPCreate1(2, 2, 1, Network);
    MLPTrainLM(Network, XY, 4, 0.001, 10, Info, Rep);
    TrnErrors := TrnErrors or AP_FP_Greater(MLPRMSError(Network, XY, 4),0.1);
    
    //
    // Test CV on random noisy problem
    //
    NCount := 100;
    SetLength(XY, NCount-1+1, 1+1);
    I:=0;
    while I<=NCount-1 do
    begin
        XY[I,0] := 2*RandomReal-1;
        XY[I,1] := RandomInteger(4);
        Inc(I);
    end;
    MLPCreateC0(1, 4, Network);
    MLPKFoldCVLM(Network, XY, NCount, 0.001, 5, 10, Info, Rep, CVRep);
    
    //
    // Final report
    //
    WasErrors := InfErrors or ProcErrors or GradErrors or HessErrors or TrnErrors;
    if  not Silent then
    begin
        Write(Format('MLP TEST'#13#10'',[]));
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
        Write(Format('GRADIENT CALCULATION:                    ',[]));
        if  not GradErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('HESSIAN CALCULATION:                     ',[]));
        if  not HessErrors then
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
procedure CreateNetwork(var Network : MultiLayerPerceptron;
     NKind : AlglibInteger;
     A1 : Double;
     A2 : Double;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger);
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
            MLPCreate0(NIn, NOut, Network);
        end
        else
        begin
            if NKind=1 then
            begin
                MLPCreateC0(NIn, NOut, Network);
            end
            else
            begin
                if NKind=2 then
                begin
                    MLPCreateB0(NIn, NOut, A1, A2, Network);
                end
                else
                begin
                    if NKind=3 then
                    begin
                        MLPCreateR0(NIn, NOut, A1, A2, Network);
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
            MLPCreate1(NIn, NHid1, NOut, Network);
        end
        else
        begin
            if NKind=1 then
            begin
                MLPCreateC1(NIn, NHid1, NOut, Network);
            end
            else
            begin
                if NKind=2 then
                begin
                    MLPCreateB1(NIn, NHid1, NOut, A1, A2, Network);
                end
                else
                begin
                    if NKind=3 then
                    begin
                        MLPCreateR1(NIn, NHid1, NOut, A1, A2, Network);
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
        MLPCreate2(NIn, NHid1, NHid2, NOut, Network);
    end
    else
    begin
        if NKind=1 then
        begin
            MLPCreateC2(NIn, NHid1, NHid2, NOut, Network);
        end
        else
        begin
            if NKind=2 then
            begin
                MLPCreateB2(NIn, NHid1, NHid2, NOut, A1, A2, Network);
            end
            else
            begin
                if NKind=3 then
                begin
                    MLPCreateR2(NIn, NHid1, NHid2, NOut, A1, A2, Network);
                end;
            end;
        end;
    end;
end;


(*************************************************************************
Unsets network (initialize it to smallest network possible
*************************************************************************)
procedure UnsetNetwork(var Network : MultiLayerPerceptron);
begin
    MLPCreate0(1, 1, Network);
end;


(*************************************************************************
Iformational functions test
*************************************************************************)
procedure TestInformational(NKind : AlglibInteger;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     PassCount : AlglibInteger;
     var Err : Boolean);
var
    Network : MultiLayerPerceptron;
    N1 : AlglibInteger;
    N2 : AlglibInteger;
    WCount : AlglibInteger;
begin
    CreateNetwork(Network, NKind, 0.0, 0.0, NIn, NHid1, NHid2, NOut);
    MLPProperties(Network, N1, N2, WCount);
    Err := Err or (N1<>NIn) or (N2<>NOut) or (WCount<=0);
end;


(*************************************************************************
Processing functions test
*************************************************************************)
procedure TestProcessing(NKind : AlglibInteger;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     PassCount : AlglibInteger;
     var Err : Boolean);
var
    Network : MultiLayerPerceptron;
    Network2 : MultiLayerPerceptron;
    N1 : AlglibInteger;
    N2 : AlglibInteger;
    WCount : AlglibInteger;
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
    Assert(PassCount>=2, 'PassCount<2!');
    
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
    CreateNetwork(Network, NKind, A1, A2, NIn, NHid1, NHid2, NOut);
    MLPProperties(Network, N1, N2, WCount);
    
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
        
        //
        // Last run is made on zero network
        //
        MLPRandomizeFull(Network);
        ZeroNet := False;
        if Pass=PassCount then
        begin
            APVMul(@Network.Weights[0], 0, WCount-1, 0);
            ZeroNet := True;
        end;
        
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
        MLPProcess(Network, X1, Y1);
        MLPProcess(Network, X2, Y2);
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
        UnsetNetwork(Network2);
        MLPCopy(Network, Network2);
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
        MLPProcess(Network, X1, Y1);
        MLPProcess(Network2, X2, Y2);
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
        UnsetNetwork(Network2);
        MLPSerialize(Network, RA, RLen);
        SetLength(RA2, RLen-1+1);
        I:=0;
        while I<=RLen-1 do
        begin
            RA2[I] := RA[I];
            Inc(I);
        end;
        MLPUnserialize(RA2, Network2);
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
        MLPProcess(Network, X1, Y1);
        MLPProcess(Network2, X2, Y2);
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
        if  not ZeroNet then
        begin
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
            MLPProcess(Network, X1, Y1);
            MLPProcess(Network, X2, Y2);
            AllSame := True;
            I:=0;
            while I<=NOut-1 do
            begin
                AllSame := AllSame and AP_FP_Eq(Y1[I],Y2[I]);
                Inc(I);
            end;
            Err := Err or AllSame;
        end;
        
        //
        // Randomization changes outputs (when inputs are unchanged, non-zero network)
        //
        if  not ZeroNet then
        begin
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
            MLPCopy(Network, Network2);
            MLPRandomize(Network2);
            MLPProcess(Network, X1, Y1);
            MLPProcess(Network2, X1, Y2);
            AllSame := True;
            I:=0;
            while I<=NOut-1 do
            begin
                AllSame := AllSame and AP_FP_Eq(Y1[I],Y2[I]);
                Inc(I);
            end;
            Err := Err or AllSame;
        end;
        
        //
        // Full randomization changes outputs (when inputs are unchanged, non-zero network)
        //
        if  not ZeroNet then
        begin
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
            MLPCopy(Network, Network2);
            MLPRandomizeFull(Network2);
            MLPProcess(Network, X1, Y1);
            MLPProcess(Network2, X1, Y2);
            AllSame := True;
            I:=0;
            while I<=NOut-1 do
            begin
                AllSame := AllSame and AP_FP_Eq(Y1[I],Y2[I]);
                Inc(I);
            end;
            Err := Err or AllSame;
        end;
        
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
            MLPProcess(Network, X1, Y1);
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
            MLPProcess(Network, X1, Y1);
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
            MLPProcess(Network, X1, Y1);
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
Gradient functions test
*************************************************************************)
procedure TestGradient(NKind : AlglibInteger;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     PassCount : AlglibInteger;
     var Err : Boolean);
var
    Network : MultiLayerPerceptron;
    Network2 : MultiLayerPerceptron;
    N1 : AlglibInteger;
    N2 : AlglibInteger;
    WCount : AlglibInteger;
    ZeroNet : Boolean;
    H : Double;
    ETol : Double;
    A1 : Double;
    A2 : Double;
    Pass : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    AllSame : Boolean;
    ILen : AlglibInteger;
    RLen : AlglibInteger;
    SSize : AlglibInteger;
    XY : TReal2DArray;
    Grad1 : TReal1DArray;
    Grad2 : TReal1DArray;
    X : TReal1DArray;
    Y : TReal1DArray;
    X1 : TReal1DArray;
    X2 : TReal1DArray;
    Y1 : TReal1DArray;
    Y2 : TReal1DArray;
    IA : TInteger1DArray;
    RA : TReal1DArray;
    V : Double;
    E : Double;
    E1 : Double;
    E2 : Double;
    V1 : Double;
    V2 : Double;
    V3 : Double;
    V4 : Double;
    WPrev : Double;
begin
    Assert(PassCount>=2, 'PassCount<2!');
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
    CreateNetwork(Network, NKind, A1, A2, NIn, NHid1, NHid2, NOut);
    MLPProperties(Network, N1, N2, WCount);
    H := 0.0001;
    ETol := 0.01;
    
    //
    // Initialize
    //
    SetLength(X, NIn-1+1);
    SetLength(X1, NIn-1+1);
    SetLength(X2, NIn-1+1);
    SetLength(Y, NOut-1+1);
    SetLength(Y1, NOut-1+1);
    SetLength(Y2, NOut-1+1);
    SetLength(Grad1, WCount-1+1);
    SetLength(Grad2, WCount-1+1);
    
    //
    // Process
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        MLPRandomizeFull(Network);
        
        //
        // Test error/gradient calculation (least squares)
        //
        SetLength(XY, 0+1, NIn+NOut-1+1);
        I:=0;
        while I<=NIn-1 do
        begin
            X[I] := 4*RandomReal-2;
            Inc(I);
        end;
        APVMove(@XY[0][0], 0, NIn-1, @X[0], 0, NIn-1);
        if MLPIsSoftmax(Network) then
        begin
            I:=0;
            while I<=NOut-1 do
            begin
                Y[I] := 0;
                Inc(I);
            end;
            XY[0,NIn] := RandomInteger(NOut);
            Y[Round(XY[0,NIn])] := 1;
        end
        else
        begin
            I:=0;
            while I<=NOut-1 do
            begin
                Y[I] := 4*RandomReal-2;
                Inc(I);
            end;
            APVMove(@XY[0][0], NIn, NIn+NOut-1, @Y[0], 0, NOut-1);
        end;
        MLPGrad(Network, X, Y, E, Grad2);
        MLPProcess(Network, X, Y2);
        APVSub(@Y2[0], 0, NOut-1, @Y[0], 0, NOut-1);
        V := APVDotProduct(@Y2[0], 0, NOut-1, @Y2[0], 0, NOut-1);
        V := V/2;
        Err := Err or AP_FP_Greater(AbsReal((V-E)/V),ETol);
        Err := Err or AP_FP_Greater(AbsReal((MLPError(Network, XY, 1)-V)/V),ETol);
        I:=0;
        while I<=WCount-1 do
        begin
            WPrev := Network.Weights[I];
            Network.Weights[I] := WPrev-2*H;
            MLPProcess(Network, X, Y1);
            APVSub(@Y1[0], 0, NOut-1, @Y[0], 0, NOut-1);
            V1 := APVDotProduct(@Y1[0], 0, NOut-1, @Y1[0], 0, NOut-1);
            V1 := V1/2;
            Network.Weights[I] := WPrev-H;
            MLPProcess(Network, X, Y1);
            APVSub(@Y1[0], 0, NOut-1, @Y[0], 0, NOut-1);
            V2 := APVDotProduct(@Y1[0], 0, NOut-1, @Y1[0], 0, NOut-1);
            V2 := V2/2;
            Network.Weights[I] := WPrev+H;
            MLPProcess(Network, X, Y1);
            APVSub(@Y1[0], 0, NOut-1, @Y[0], 0, NOut-1);
            V3 := APVDotProduct(@Y1[0], 0, NOut-1, @Y1[0], 0, NOut-1);
            V3 := V3/2;
            Network.Weights[I] := WPrev+2*H;
            MLPProcess(Network, X, Y1);
            APVSub(@Y1[0], 0, NOut-1, @Y[0], 0, NOut-1);
            V4 := APVDotProduct(@Y1[0], 0, NOut-1, @Y1[0], 0, NOut-1);
            V4 := V4/2;
            Network.Weights[I] := WPrev;
            Grad1[I] := (V1-8*V2+8*V3-V4)/(12*H);
            if AP_FP_Greater(AbsReal(Grad1[I]),1.0E-3) then
            begin
                Err := Err or AP_FP_Greater(AbsReal((Grad2[I]-Grad1[I])/Grad1[I]),ETol);
            end
            else
            begin
                Err := Err or AP_FP_Greater(AbsReal(Grad2[I]-Grad1[I]),ETol);
            end;
            Inc(I);
        end;
        
        //
        // Test error/gradient calculation (natural).
        // Testing on non-random structure networks
        // (because NKind is representative only in that case).
        //
        SetLength(XY, 0+1, NIn+NOut-1+1);
        I:=0;
        while I<=NIn-1 do
        begin
            X[I] := 4*RandomReal-2;
            Inc(I);
        end;
        APVMove(@XY[0][0], 0, NIn-1, @X[0], 0, NIn-1);
        if MLPIsSoftmax(Network) then
        begin
            I:=0;
            while I<=NOut-1 do
            begin
                Y[I] := 0;
                Inc(I);
            end;
            XY[0,NIn] := RandomInteger(NOut);
            Y[Round(XY[0,NIn])] := 1;
        end
        else
        begin
            I:=0;
            while I<=NOut-1 do
            begin
                Y[I] := 4*RandomReal-2;
                Inc(I);
            end;
            APVMove(@XY[0][0], NIn, NIn+NOut-1, @Y[0], 0, NOut-1);
        end;
        MLPGradN(Network, X, Y, E, Grad2);
        MLPProcess(Network, X, Y2);
        V := 0;
        if NKind<>1 then
        begin
            I:=0;
            while I<=NOut-1 do
            begin
                V := V+0.5*AP_Sqr(Y2[I]-Y[I]);
                Inc(I);
            end;
        end
        else
        begin
            I:=0;
            while I<=NOut-1 do
            begin
                if AP_FP_Neq(Y[I],0) then
                begin
                    if AP_FP_Eq(Y2[I],0) then
                    begin
                        V := V+Y[I]*Ln(MaxRealNumber);
                    end
                    else
                    begin
                        V := V+Y[I]*Ln(Y[I]/Y2[I]);
                    end;
                end;
                Inc(I);
            end;
        end;
        Err := Err or AP_FP_Greater(AbsReal((V-E)/V),ETol);
        Err := Err or AP_FP_Greater(AbsReal((MLPErrorN(Network, XY, 1)-V)/V),ETol);
        I:=0;
        while I<=WCount-1 do
        begin
            WPrev := Network.Weights[I];
            Network.Weights[I] := WPrev+H;
            MLPProcess(Network, X, Y2);
            Network.Weights[I] := WPrev-H;
            MLPProcess(Network, X, Y1);
            Network.Weights[I] := WPrev;
            V := 0;
            if NKind<>1 then
            begin
                J:=0;
                while J<=NOut-1 do
                begin
                    V := V+0.5*(AP_Sqr(Y2[J]-Y[J])-AP_Sqr(Y1[J]-Y[J]))/(2*H);
                    Inc(J);
                end;
            end
            else
            begin
                J:=0;
                while J<=NOut-1 do
                begin
                    if AP_FP_Neq(Y[J],0) then
                    begin
                        if AP_FP_Eq(Y2[J],0) then
                        begin
                            V := V+Y[J]*Ln(MaxRealNumber);
                        end
                        else
                        begin
                            V := V+Y[J]*Ln(Y[J]/Y2[J]);
                        end;
                        if AP_FP_Eq(Y1[J],0) then
                        begin
                            V := V-Y[J]*Ln(MaxRealNumber);
                        end
                        else
                        begin
                            V := V-Y[J]*Ln(Y[J]/Y1[J]);
                        end;
                    end;
                    Inc(J);
                end;
                V := V/(2*H);
            end;
            Grad1[I] := V;
            if AP_FP_Greater(AbsReal(Grad1[I]),1.0E-3) then
            begin
                Err := Err or AP_FP_Greater(AbsReal((Grad2[I]-Grad1[I])/Grad1[I]),ETol);
            end
            else
            begin
                Err := Err or AP_FP_Greater(AbsReal(Grad2[I]-Grad1[I]),ETol);
            end;
            Inc(I);
        end;
        
        //
        // Test gradient calculation: batch (least squares)
        //
        SSize := 1+RandomInteger(10);
        SetLength(XY, SSize-1+1, NIn+NOut-1+1);
        I:=0;
        while I<=WCount-1 do
        begin
            Grad1[I] := 0;
            Inc(I);
        end;
        E1 := 0;
        I:=0;
        while I<=SSize-1 do
        begin
            J:=0;
            while J<=NIn-1 do
            begin
                X1[J] := 4*RandomReal-2;
                Inc(J);
            end;
            APVMove(@XY[I][0], 0, NIn-1, @X1[0], 0, NIn-1);
            if MLPIsSoftmax(Network) then
            begin
                J:=0;
                while J<=NOut-1 do
                begin
                    Y1[J] := 0;
                    Inc(J);
                end;
                XY[I,NIn] := RandomInteger(NOut);
                Y1[Round(XY[I,NIn])] := 1;
            end
            else
            begin
                J:=0;
                while J<=NOut-1 do
                begin
                    Y1[J] := 4*RandomReal-2;
                    Inc(J);
                end;
                APVMove(@XY[I][0], NIn, NIn+NOut-1, @Y1[0], 0, NOut-1);
            end;
            MLPGrad(Network, X1, Y1, V, Grad2);
            E1 := E1+V;
            APVAdd(@Grad1[0], 0, WCount-1, @Grad2[0], 0, WCount-1);
            Inc(I);
        end;
        MLPGradBatch(Network, XY, SSize, E2, Grad2);
        Err := Err or AP_FP_Greater(AbsReal(E1-E2)/E1,0.01);
        I:=0;
        while I<=WCount-1 do
        begin
            if AP_FP_Neq(Grad1[I],0) then
            begin
                Err := Err or AP_FP_Greater(AbsReal((Grad2[I]-Grad1[I])/Grad1[I]),ETol);
            end
            else
            begin
                Err := Err or AP_FP_Neq(Grad2[I],Grad1[I]);
            end;
            Inc(I);
        end;
        
        //
        // Test gradient calculation: batch (natural error func)
        //
        SSize := 1+RandomInteger(10);
        SetLength(XY, SSize-1+1, NIn+NOut-1+1);
        I:=0;
        while I<=WCount-1 do
        begin
            Grad1[I] := 0;
            Inc(I);
        end;
        E1 := 0;
        I:=0;
        while I<=SSize-1 do
        begin
            J:=0;
            while J<=NIn-1 do
            begin
                X1[J] := 4*RandomReal-2;
                Inc(J);
            end;
            APVMove(@XY[I][0], 0, NIn-1, @X1[0], 0, NIn-1);
            if MLPIsSoftmax(Network) then
            begin
                J:=0;
                while J<=NOut-1 do
                begin
                    Y1[J] := 0;
                    Inc(J);
                end;
                XY[I,NIn] := RandomInteger(NOut);
                Y1[Round(XY[I,NIn])] := 1;
            end
            else
            begin
                J:=0;
                while J<=NOut-1 do
                begin
                    Y1[J] := 4*RandomReal-2;
                    Inc(J);
                end;
                APVMove(@XY[I][0], NIn, NIn+NOut-1, @Y1[0], 0, NOut-1);
            end;
            MLPGradN(Network, X1, Y1, V, Grad2);
            E1 := E1+V;
            APVAdd(@Grad1[0], 0, WCount-1, @Grad2[0], 0, WCount-1);
            Inc(I);
        end;
        MLPGradNBatch(Network, XY, SSize, E2, Grad2);
        Err := Err or AP_FP_Greater(AbsReal(E1-E2)/E1,ETol);
        I:=0;
        while I<=WCount-1 do
        begin
            if AP_FP_Neq(Grad1[I],0) then
            begin
                Err := Err or AP_FP_Greater(AbsReal((Grad2[I]-Grad1[I])/Grad1[I]),ETol);
            end
            else
            begin
                Err := Err or AP_FP_Neq(Grad2[I],Grad1[I]);
            end;
            Inc(I);
        end;
        Inc(Pass);
    end;
end;


(*************************************************************************
Hessian functions test
*************************************************************************)
procedure TestHessian(NKind : AlglibInteger;
     NIn : AlglibInteger;
     NHid1 : AlglibInteger;
     NHid2 : AlglibInteger;
     NOut : AlglibInteger;
     PassCount : AlglibInteger;
     var Err : Boolean);
var
    Network : MultiLayerPerceptron;
    Network2 : MultiLayerPerceptron;
    HKind : AlglibInteger;
    N1 : AlglibInteger;
    N2 : AlglibInteger;
    WCount : AlglibInteger;
    ZeroNet : Boolean;
    H : Double;
    ETol : Double;
    Pass : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    AllSame : Boolean;
    ILen : AlglibInteger;
    RLen : AlglibInteger;
    SSize : AlglibInteger;
    A1 : Double;
    A2 : Double;
    XY : TReal2DArray;
    H1 : TReal2DArray;
    H2 : TReal2DArray;
    Grad1 : TReal1DArray;
    Grad2 : TReal1DArray;
    Grad3 : TReal1DArray;
    X : TReal1DArray;
    Y : TReal1DArray;
    X1 : TReal1DArray;
    X2 : TReal1DArray;
    Y1 : TReal1DArray;
    Y2 : TReal1DArray;
    IA : TInteger1DArray;
    RA : TReal1DArray;
    V : Double;
    E : Double;
    E1 : Double;
    E2 : Double;
    V1 : Double;
    V2 : Double;
    V3 : Double;
    V4 : Double;
    WPrev : Double;
begin
    Assert(PassCount>=2, 'PassCount<2!');
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
    CreateNetwork(Network, NKind, A1, A2, NIn, NHid1, NHid2, NOut);
    MLPProperties(Network, N1, N2, WCount);
    H := 0.0001;
    ETol := 0.05;
    
    //
    // Initialize
    //
    SetLength(X, NIn-1+1);
    SetLength(X1, NIn-1+1);
    SetLength(X2, NIn-1+1);
    SetLength(Y, NOut-1+1);
    SetLength(Y1, NOut-1+1);
    SetLength(Y2, NOut-1+1);
    SetLength(Grad1, WCount-1+1);
    SetLength(Grad2, WCount-1+1);
    SetLength(Grad3, WCount-1+1);
    SetLength(H1, WCount-1+1, WCount-1+1);
    SetLength(H2, WCount-1+1, WCount-1+1);
    
    //
    // Process
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        MLPRandomizeFull(Network);
        
        //
        // Test hessian calculation .
        // E1 contains total error (calculated using MLPGrad/MLPGradN)
        // Grad1 contains total gradient (calculated using MLPGrad/MLPGradN)
        // H1 contains Hessian calculated using differences of gradients
        //
        // E2, Grad2 and H2 contains corresponing values calculated using MLPHessianBatch/MLPHessianNBatch
        //
        HKind:=0;
        while HKind<=1 do
        begin
            SSize := 1+RandomInteger(10);
            SetLength(XY, SSize-1+1, NIn+NOut-1+1);
            I:=0;
            while I<=WCount-1 do
            begin
                Grad1[I] := 0;
                Inc(I);
            end;
            I:=0;
            while I<=WCount-1 do
            begin
                J:=0;
                while J<=WCount-1 do
                begin
                    H1[I,J] := 0;
                    Inc(J);
                end;
                Inc(I);
            end;
            E1 := 0;
            I:=0;
            while I<=SSize-1 do
            begin
                
                //
                // X, Y
                //
                J:=0;
                while J<=NIn-1 do
                begin
                    X1[J] := 4*RandomReal-2;
                    Inc(J);
                end;
                APVMove(@XY[I][0], 0, NIn-1, @X1[0], 0, NIn-1);
                if MLPIsSoftmax(Network) then
                begin
                    J:=0;
                    while J<=NOut-1 do
                    begin
                        Y1[J] := 0;
                        Inc(J);
                    end;
                    XY[I,NIn] := RandomInteger(NOut);
                    Y1[Round(XY[I,NIn])] := 1;
                end
                else
                begin
                    J:=0;
                    while J<=NOut-1 do
                    begin
                        Y1[J] := 4*RandomReal-2;
                        Inc(J);
                    end;
                    APVMove(@XY[I][0], NIn, NIn+NOut-1, @Y1[0], 0, NOut-1);
                end;
                
                //
                // E1, Grad1
                //
                if HKind=0 then
                begin
                    MLPGrad(Network, X1, Y1, V, Grad2);
                end
                else
                begin
                    MLPGradN(Network, X1, Y1, V, Grad2);
                end;
                E1 := E1+V;
                APVAdd(@Grad1[0], 0, WCount-1, @Grad2[0], 0, WCount-1);
                
                //
                // H1
                //
                J:=0;
                while J<=WCount-1 do
                begin
                    WPrev := Network.Weights[J];
                    Network.Weights[J] := WPrev-2*H;
                    if HKind=0 then
                    begin
                        MLPGrad(Network, X1, Y1, V, Grad2);
                    end
                    else
                    begin
                        MLPGradN(Network, X1, Y1, V, Grad2);
                    end;
                    Network.Weights[J] := WPrev-H;
                    if HKind=0 then
                    begin
                        MLPGrad(Network, X1, Y1, V, Grad3);
                    end
                    else
                    begin
                        MLPGradN(Network, X1, Y1, V, Grad3);
                    end;
                    APVSub(@Grad2[0], 0, WCount-1, @Grad3[0], 0, WCount-1, 8);
                    Network.Weights[J] := WPrev+H;
                    if HKind=0 then
                    begin
                        MLPGrad(Network, X1, Y1, V, Grad3);
                    end
                    else
                    begin
                        MLPGradN(Network, X1, Y1, V, Grad3);
                    end;
                    APVAdd(@Grad2[0], 0, WCount-1, @Grad3[0], 0, WCount-1, 8);
                    Network.Weights[J] := WPrev+2*H;
                    if HKind=0 then
                    begin
                        MLPGrad(Network, X1, Y1, V, Grad3);
                    end
                    else
                    begin
                        MLPGradN(Network, X1, Y1, V, Grad3);
                    end;
                    APVSub(@Grad2[0], 0, WCount-1, @Grad3[0], 0, WCount-1);
                    V := 1/(12*H);
                    APVAdd(@H1[J][0], 0, WCount-1, @Grad2[0], 0, WCount-1, V);
                    Network.Weights[J] := WPrev;
                    Inc(J);
                end;
                Inc(I);
            end;
            if HKind=0 then
            begin
                MLPHessianBatch(Network, XY, SSize, E2, Grad2, H2);
            end
            else
            begin
                MLPHessianNBatch(Network, XY, SSize, E2, Grad2, H2);
            end;
            Err := Err or AP_FP_Greater(AbsReal(E1-E2)/E1,ETol);
            I:=0;
            while I<=WCount-1 do
            begin
                if AP_FP_Greater(AbsReal(Grad1[I]),1.0E-2) then
                begin
                    Err := Err or AP_FP_Greater(AbsReal((Grad2[I]-Grad1[I])/Grad1[I]),ETol);
                end
                else
                begin
                    Err := Err or AP_FP_Greater(AbsReal(Grad2[I]-Grad1[I]),ETol);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=WCount-1 do
            begin
                J:=0;
                while J<=WCount-1 do
                begin
                    if AP_FP_Greater(AbsReal(H1[I,J]),5.0E-2) then
                    begin
                        Err := Err or AP_FP_Greater(AbsReal((H1[I,J]-H2[I,J])/H1[I,J]),ETol);
                    end
                    else
                    begin
                        Err := Err or AP_FP_Greater(AbsReal(H2[I,J]-H1[I,J]),ETol);
                    end;
                    Inc(J);
                end;
                Inc(I);
            end;
            Inc(HKind);
        end;
        Inc(Pass);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testmlpunit_test_silent():Boolean;
begin
    Result := TestMLP(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testmlpunit_test():Boolean;
begin
    Result := TestMLP(False);
end;


end.