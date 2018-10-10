unit testregressunit;
interface
uses Math, Sysutils, Ap, descriptivestatistics, gammafunc, normaldistr, igammaf, reflections, bidiagonal, qr, lq, blas, rotations, bdsvd, svd, linreg;

function TestLinRegression(Silent : Boolean):Boolean;
function testregressunit_test_silent():Boolean;
function testregressunit_test():Boolean;

implementation

procedure GenerateRandomTask(XL : Double;
     XR : Double;
     RandomX : Boolean;
     YMin : Double;
     YMax : Double;
     SMin : Double;
     SMax : Double;
     N : AlglibInteger;
     var XY : TReal2DArray;
     var S : TReal1DArray);forward;
procedure GenerateTask(A : Double;
     B : Double;
     XL : Double;
     XR : Double;
     RandomX : Boolean;
     SMin : Double;
     SMax : Double;
     N : AlglibInteger;
     var XY : TReal2DArray;
     var S : TReal1DArray);forward;
procedure FillTaskWithY(A : Double;
     B : Double;
     N : AlglibInteger;
     var XY : TReal2DArray;
     var S : TReal1DArray);forward;
function GenerateNormal(Mean : Double; Sigma : Double):Double;forward;
procedure CalculateMV(const X : TReal1DArray;
     N : AlglibInteger;
     var Mean : Double;
     var MeanS : Double;
     var StdDev : Double;
     var StdDevS : Double);forward;
procedure UnsetLR(var LR : LinearModel);forward;


function TestLinRegression(Silent : Boolean):Boolean;
var
    SigmaThreshold : Double;
    MaxN : AlglibInteger;
    MaxM : AlglibInteger;
    PassCount : AlglibInteger;
    EstPassCount : AlglibInteger;
    Threshold : Double;
    N : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    TmpI : AlglibInteger;
    Pass : AlglibInteger;
    EPass : AlglibInteger;
    M : AlglibInteger;
    TaskType : AlglibInteger;
    ModelType : AlglibInteger;
    M1 : AlglibInteger;
    M2 : AlglibInteger;
    N1 : AlglibInteger;
    N2 : AlglibInteger;
    Info : AlglibInteger;
    Info2 : AlglibInteger;
    XY : TReal2DArray;
    XY2 : TReal2DArray;
    S : TReal1DArray;
    S2 : TReal1DArray;
    W2 : TReal1DArray;
    X : TReal1DArray;
    TA : TReal1DArray;
    TB : TReal1DArray;
    TC : TReal1DArray;
    XY0 : TReal1DArray;
    TmpWeights : TReal1DArray;
    W : LinearModel;
    WT : LinearModel;
    WT2 : LinearModel;
    X1 : TReal1DArray;
    X2 : TReal1DArray;
    RA : TReal1DArray;
    RA2 : TReal1DArray;
    Y1 : Double;
    Y2 : Double;
    RLen : AlglibInteger;
    AllSame : Boolean;
    EA : Double;
    EB : Double;
    VarATested : Double;
    VarBTested : Double;
    A : Double;
    B : Double;
    VarA : Double;
    VarB : Double;
    A2 : Double;
    B2 : Double;
    CovAB : Double;
    CorrAB : Double;
    P : Double;
    QCnt : AlglibInteger;
    QTbl : TReal1DArray;
    QVals : TReal1DArray;
    QSigma : TReal1DArray;
    AR : LRReport;
    AR2 : LRReport;
    F : Double;
    FP : Double;
    FM : Double;
    V : Double;
    VV : Double;
    CVRMSError : Double;
    CVAvgError : Double;
    CVAvgRelError : Double;
    RMSError : Double;
    AvgError : Double;
    AvgRelError : Double;
    NonDefect : Boolean;
    SinShift : Double;
    TaskLevel : Double;
    NoiseLevel : Double;
    HStep : Double;
    Sigma : Double;
    Mean : Double;
    MeanS : Double;
    StdDev : Double;
    StdDevS : Double;
    SLCErrors : Boolean;
    SLErrors : Boolean;
    GRCovErrors : Boolean;
    GROptErrors : Boolean;
    GREstErrors : Boolean;
    GROtherErrors : Boolean;
    GRConvErrors : Boolean;
    WasErrors : Boolean;
    i_ : AlglibInteger;
begin
    
    //
    // Primary settings
    //
    MaxN := 40;
    MaxM := 5;
    PassCount := 3;
    EstPassCount := 1000;
    SigmaThreshold := 7;
    Threshold := 1000000*MachineEpsilon;
    SLErrors := False;
    SLCErrors := False;
    GRCovErrors := False;
    GROptErrors := False;
    GREstErrors := False;
    GROtherErrors := False;
    GRConvErrors := False;
    WasErrors := False;
    
    //
    // Quantiles table setup
    //
    QCnt := 5;
    SetLength(QTbl, QCnt-1+1);
    SetLength(QVals, QCnt-1+1);
    SetLength(QSigma, QCnt-1+1);
    QTbl[0] := 0.5;
    QTbl[1] := 0.25;
    QTbl[2] := 0.10;
    QTbl[3] := 0.05;
    QTbl[4] := 0.025;
    I:=0;
    while I<=QCnt-1 do
    begin
        QSigma[I] := Sqrt(QTbl[I]*(1-QTbl[I])/EstPassCount);
        Inc(I);
    end;
    
    //
    // Other setup
    //
    SetLength(TA, EstPassCount-1+1);
    SetLength(TB, EstPassCount-1+1);
    
    //
    // Test straight line regression
    //
    N:=2;
    while N<=MaxN do
    begin
        
        //
        // Fail/pass test
        //
        GenerateRandomTask(-1, 1, False, -1, 1, 1, 2, N, XY, S);
        LRLineS(XY, S, N, Info, A, B, VarA, VarB, CovAB, CorrAB, P);
        SLCErrors := SLCErrors or (Info<>1);
        GenerateRandomTask(+1, 1, False, -1, 1, 1, 2, N, XY, S);
        LRLineS(XY, S, N, Info, A, B, VarA, VarB, CovAB, CorrAB, P);
        SLCErrors := SLCErrors or (Info<>-3);
        GenerateRandomTask(-1, 1, False, -1, 1, -1, -1, N, XY, S);
        LRLineS(XY, S, N, Info, A, B, VarA, VarB, CovAB, CorrAB, P);
        SLCErrors := SLCErrors or (Info<>-2);
        GenerateRandomTask(-1, 1, False, -1, 1, 2, 1, 2, XY, S);
        LRLineS(XY, S, 1, Info, A, B, VarA, VarB, CovAB, CorrAB, P);
        SLCErrors := SLCErrors or (Info<>-1);
        
        //
        // Multipass tests
        //
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // Test S variant against non-S variant
            //
            EA := 2*RandomReal-1;
            EB := 2*RandomReal-1;
            GenerateTask(EA, EB, -5*RandomReal, +5*RandomReal, AP_FP_Greater(RandomReal,0.5), 1, 1, N, XY, S);
            LRLineS(XY, S, N, Info, A, B, VarA, VarB, CovAB, CorrAB, P);
            LRLine(XY, N, Info2, A2, B2);
            if (Info<>1) or (Info2<>1) then
            begin
                SLCErrors := True;
            end
            else
            begin
                SLErrors := SLErrors or AP_FP_Greater(AbsReal(A-A2),Threshold) or AP_FP_Greater(AbsReal(B-B2),Threshold);
            end;
            
            //
            // Test for A/B
            //
            // Generate task with exact, non-perturbed y[i],
            // then make non-zero s[i]
            //
            EA := 2*RandomReal-1;
            EB := 2*RandomReal-1;
            GenerateTask(EA, EB, -5*RandomReal, +5*RandomReal, N>4, 0.0, 0.0, N, XY, S);
            I:=0;
            while I<=N-1 do
            begin
                S[I] := 1+RandomReal;
                Inc(I);
            end;
            LRLineS(XY, S, N, Info, A, B, VarA, VarB, CovAB, CorrAB, P);
            if Info<>1 then
            begin
                SLCErrors := True;
            end
            else
            begin
                SLErrors := SLErrors or AP_FP_Greater(AbsReal(A-EA),0.001) or AP_FP_Greater(AbsReal(B-EB),0.001);
            end;
            
            //
            // Test for VarA, VarB, P (P is being tested only for N>2)
            //
            I:=0;
            while I<=QCnt-1 do
            begin
                QVals[I] := 0;
                Inc(I);
            end;
            EA := 2*RandomReal-1;
            EB := 2*RandomReal-1;
            GenerateTask(EA, EB, -5*RandomReal, +5*RandomReal, N>4, 1.0, 2.0, N, XY, S);
            LRLineS(XY, S, N, Info, A, B, VarA, VarB, CovAB, CorrAB, P);
            if Info<>1 then
            begin
                SLCErrors := True;
                Inc(Pass);
                Continue;
            end;
            VarATested := VarA;
            VarBTested := VarB;
            EPass:=0;
            while EPass<=EstPassCount-1 do
            begin
                
                //
                // Generate
                //
                FillTaskWithY(EA, EB, N, XY, S);
                LRLineS(XY, S, N, Info, A, B, VarA, VarB, CovAB, CorrAB, P);
                if Info<>1 then
                begin
                    SLCErrors := True;
                    Inc(EPass);
                    Continue;
                end;
                
                //
                // A, B, P
                // (P is being tested for uniformity, additional p-tests are below)
                //
                TA[EPass] := A;
                TB[EPass] := B;
                I:=0;
                while I<=QCnt-1 do
                begin
                    if AP_FP_Less_Eq(P,QTbl[I]) then
                    begin
                        QVals[I] := QVals[I]+AP_Double(1)/EstPassCount;
                    end;
                    Inc(I);
                end;
                Inc(EPass);
            end;
            CalculateMV(TA, EstPassCount, Mean, MeanS, StdDev, StdDevS);
            SLErrors := SLErrors or AP_FP_Greater_Eq(AbsReal(Mean-EA)/MeanS,SigmaThreshold);
            SLErrors := SLErrors or AP_FP_Greater_Eq(AbsReal(StdDev-Sqrt(VarATested))/StdDevS,SigmaThreshold);
            CalculateMV(TB, EstPassCount, Mean, MeanS, StdDev, StdDevS);
            SLErrors := SLErrors or AP_FP_Greater_Eq(AbsReal(Mean-EB)/MeanS,SigmaThreshold);
            SLErrors := SLErrors or AP_FP_Greater_Eq(AbsReal(StdDev-Sqrt(VarBTested))/StdDevS,SigmaThreshold);
            if N>2 then
            begin
                I:=0;
                while I<=QCnt-1 do
                begin
                    if AP_FP_Greater(AbsReal(QTbl[I]-QVals[I])/QSigma[I],SigmaThreshold) then
                    begin
                        SLErrors := True;
                    end;
                    Inc(I);
                end;
            end;
            
            //
            // Additional tests for P: correlation with fit quality
            //
            if N>2 then
            begin
                GenerateTask(EA, EB, -5*RandomReal, +5*RandomReal, False, 0.0, 0.0, N, XY, S);
                I:=0;
                while I<=N-1 do
                begin
                    S[I] := 1+RandomReal;
                    Inc(I);
                end;
                LRLineS(XY, S, N, Info, A, B, VarA, VarB, CovAB, CorrAB, P);
                if Info<>1 then
                begin
                    SLCErrors := True;
                    Inc(Pass);
                    Continue;
                end;
                SLErrors := SLErrors or AP_FP_Less(P,0.999);
                GenerateTask(0, 0, -5*RandomReal, +5*RandomReal, False, 1.0, 1.0, N, XY, S);
                I:=0;
                while I<=N-1 do
                begin
                    if I mod 2=0 then
                    begin
                        XY[I,1] := +5.0;
                    end
                    else
                    begin
                        XY[I,1] := -5.0;
                    end;
                    Inc(I);
                end;
                if N mod 2<>0 then
                begin
                    XY[N-1,1] := 0;
                end;
                LRLineS(XY, S, N, Info, A, B, VarA, VarB, CovAB, CorrAB, P);
                if Info<>1 then
                begin
                    SLCErrors := True;
                    Inc(Pass);
                    Continue;
                end;
                SLErrors := SLErrors or AP_FP_Greater(P,0.001);
            end;
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // General regression tests:
    //
    
    //
    // Simple linear tests (small sample, optimum point, covariance)
    //
    N:=3;
    while N<=MaxN do
    begin
        SetLength(S, N-1+1);
        
        //
        // Linear tests:
        // a. random points, sigmas
        // b. no sigmas
        //
        SetLength(XY, N-1+1, 1+1);
        I:=0;
        while I<=N-1 do
        begin
            XY[I,0] := 2*RandomReal-1;
            XY[I,1] := 2*RandomReal-1;
            S[I] := 1+RandomReal;
            Inc(I);
        end;
        LRBuildS(XY, S, N, 1, Info, WT, AR);
        if Info<>1 then
        begin
            GRConvErrors := True;
            Inc(N);
            Continue;
        end;
        LRUnpack(WT, TmpWeights, TmpI);
        LRLineS(XY, S, N, Info2, A, B, VarA, VarB, CovAB, CorrAB, P);
        GROptErrors := GROptErrors or AP_FP_Greater(AbsReal(A-TmpWeights[1]),Threshold);
        GROptErrors := GROptErrors or AP_FP_Greater(AbsReal(B-TmpWeights[0]),Threshold);
        GRCovErrors := GRCovErrors or AP_FP_Greater(AbsReal(VarA-AR.C[1,1]),Threshold);
        GRCovErrors := GRCovErrors or AP_FP_Greater(AbsReal(VarB-AR.C[0,0]),Threshold);
        GRCovErrors := GRCovErrors or AP_FP_Greater(AbsReal(CovAB-AR.C[1,0]),Threshold);
        GRCovErrors := GRCovErrors or AP_FP_Greater(AbsReal(CovAB-AR.C[0,1]),Threshold);
        LRBuild(XY, N, 1, Info, WT, AR);
        if Info<>1 then
        begin
            GRConvErrors := True;
            Inc(N);
            Continue;
        end;
        LRUnpack(WT, TmpWeights, TmpI);
        LRLine(XY, N, Info2, A, B);
        GROptErrors := GROptErrors or AP_FP_Greater(AbsReal(A-TmpWeights[1]),Threshold);
        GROptErrors := GROptErrors or AP_FP_Greater(AbsReal(B-TmpWeights[0]),Threshold);
        Inc(N);
    end;
    
    //
    // S covariance versus S-less covariance.
    // Slightly skewed task, large sample size.
    // Will S-less subroutine estimate covariance matrix good enough?
    //
    N := 1000+RandomInteger(3000);
    Sigma := 0.1+RandomReal*1.9;
    SetLength(XY, N-1+1, 1+1);
    SetLength(S, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        XY[I,0] := 1.5*RandomReal-0.5;
        XY[I,1] := 1.2*XY[I,0]-0.3+GenerateNormal(0, Sigma);
        S[I] := Sigma;
        Inc(I);
    end;
    LRBuild(XY, N, 1, Info, WT, AR);
    LRLineS(XY, S, N, Info2, A, B, VarA, VarB, CovAB, CorrAB, P);
    if (Info<>1) or (Info2<>1) then
    begin
        GRConvErrors := True;
    end
    else
    begin
        GRCovErrors := GRCovErrors or AP_FP_Greater(AbsReal(Ln(AR.C[0,0]/VarB)),Ln(1.2));
        GRCovErrors := GRCovErrors or AP_FP_Greater(AbsReal(Ln(AR.C[1,1]/VarA)),Ln(1.2));
        GRCovErrors := GRCovErrors or AP_FP_Greater(AbsReal(Ln(AR.C[0,1]/CovAB)),Ln(1.2));
        GRCovErrors := GRCovErrors or AP_FP_Greater(AbsReal(Ln(AR.C[1,0]/CovAB)),Ln(1.2));
    end;
    
    //
    // General tests:
    // * basis functions - up to cubic
    // * task types:
    // * data set is noisy sine half-period with random shift
    // * tests:
    //   unpacking/packing
    //   optimality
    //   error estimates
    // * tasks:
    //   0 = noised sine
    //   1 = degenerate task with 1-of-n encoded categorical variables
    //   2 = random task with large variation (for 1-type models)
    //   3 = random task with small variation (for 1-type models)
    //
    //   Additional tasks TODO
    //   specially designed task with defective vectors which leads to
    //   the failure of the fast CV formula.
    //
    //
    ModelType:=0;
    while ModelType<=1 do
    begin
        TaskType:=0;
        while TaskType<=3 do
        begin
            if TaskType=0 then
            begin
                M1 := 1;
                M2 := 3;
            end;
            if TaskType=1 then
            begin
                M1 := 9;
                M2 := 9;
            end;
            if (TaskType=2) or (TaskType=3) then
            begin
                M1 := 9;
                M2 := 9;
            end;
            M:=M1;
            while M<=M2 do
            begin
                if TaskType=0 then
                begin
                    N1 := M+3;
                    N2 := M+20;
                end;
                if TaskType=1 then
                begin
                    N1 := 70+RandomInteger(70);
                    N2 := N1;
                end;
                if (TaskType=2) or (TaskType=3) then
                begin
                    N1 := 100;
                    N2 := N1;
                end;
                N:=N1;
                while N<=N2 do
                begin
                    SetLength(XY, N-1+1, M+1);
                    SetLength(XY0, N-1+1);
                    SetLength(S, N-1+1);
                    HStep := 0.001;
                    NoiseLevel := 0.2;
                    
                    //
                    // Prepare task
                    //
                    if TaskType=0 then
                    begin
                        I:=0;
                        while I<=N-1 do
                        begin
                            XY[I,0] := 2*RandomReal-1;
                            Inc(I);
                        end;
                        I:=0;
                        while I<=N-1 do
                        begin
                            J:=1;
                            while J<=M-1 do
                            begin
                                XY[I,J] := XY[I,0]*XY[I,J-1];
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        SinShift := RandomReal*Pi;
                        I:=0;
                        while I<=N-1 do
                        begin
                            XY0[I] := Sin(SinShift+Pi*0.5*(XY[I,0]+1));
                            XY[I,M] := XY0[I]+NoiseLevel*GenerateNormal(0, 1);
                            Inc(I);
                        end;
                    end;
                    if TaskType=1 then
                    begin
                        Assert(M=9);
                        SetLength(TA, 8+1);
                        TA[0] := 1;
                        TA[1] := 2;
                        TA[2] := 3;
                        TA[3] := 0.25;
                        TA[4] := 0.5;
                        TA[5] := 0.75;
                        TA[6] := 0.06;
                        TA[7] := 0.12;
                        TA[8] := 0.18;
                        I:=0;
                        while I<=N-1 do
                        begin
                            J:=0;
                            while J<=M-1 do
                            begin
                                XY[I,J] := 0;
                                Inc(J);
                            end;
                            XY[I,0+I mod 3] := 1;
                            XY[I,3+I div 3 mod 3] := 1;
                            XY[I,6+I div 9 mod 3] := 1;
                            V := APVDotProduct(@XY[I][0], 0, 8, @TA[0], 0, 8);
                            XY0[I] := V;
                            XY[I,M] := V+NoiseLevel*GenerateNormal(0, 1);
                            Inc(I);
                        end;
                    end;
                    if (TaskType=2) or (TaskType=3) then
                    begin
                        Assert(M=9);
                        SetLength(TA, 8+1);
                        TA[0] := 1;
                        TA[1] := -2;
                        TA[2] := 3;
                        TA[3] := 0.25;
                        TA[4] := -0.5;
                        TA[5] := 0.75;
                        TA[6] := -0.06;
                        TA[7] := 0.12;
                        TA[8] := -0.18;
                        I:=0;
                        while I<=N-1 do
                        begin
                            J:=0;
                            while J<=M-1 do
                            begin
                                if TaskType=2 then
                                begin
                                    XY[I,J] := 1+GenerateNormal(0, 3);
                                end
                                else
                                begin
                                    XY[I,J] := 1+GenerateNormal(0, 0.05);
                                end;
                                Inc(J);
                            end;
                            V := APVDotProduct(@XY[I][0], 0, 8, @TA[0], 0, 8);
                            XY0[I] := V;
                            XY[I,M] := V+NoiseLevel*GenerateNormal(0, 1);
                            Inc(I);
                        end;
                    end;
                    I:=0;
                    while I<=N-1 do
                    begin
                        S[I] := 1+RandomReal;
                        Inc(I);
                    end;
                    
                    //
                    // Solve (using S-variant, non-S-variant is not tested)
                    //
                    if ModelType=0 then
                    begin
                        LRBuildS(XY, S, N, M, Info, WT, AR);
                    end
                    else
                    begin
                        LRBuildZS(XY, S, N, M, Info, WT, AR);
                    end;
                    if Info<>1 then
                    begin
                        GRConvErrors := True;
                        Inc(N);
                        Continue;
                    end;
                    LRUnpack(WT, TmpWeights, TmpI);
                    
                    //
                    // LRProcess test
                    //
                    SetLength(X, M-1+1);
                    V := TmpWeights[M];
                    I:=0;
                    while I<=M-1 do
                    begin
                        X[I] := 2*RandomReal-1;
                        V := V+TmpWeights[I]*X[I];
                        Inc(I);
                    end;
                    GROtherErrors := GROtherErrors or AP_FP_Greater(AbsReal(V-LRProcess(WT, X))/Max(AbsReal(V), 1),Threshold);
                    
                    //
                    // LRPack test
                    //
                    LRPack(TmpWeights, M, WT2);
                    SetLength(X, M-1+1);
                    I:=0;
                    while I<=M-1 do
                    begin
                        X[I] := 2*RandomReal-1;
                        Inc(I);
                    end;
                    V := LRProcess(WT, X);
                    GROtherErrors := GROtherErrors or AP_FP_Greater(AbsReal(V-LRProcess(WT2, X))/AbsReal(V),Threshold);
                    
                    //
                    // Optimality test
                    //
                    K:=0;
                    while K<=M do
                    begin
                        if (ModelType=1) and (K=M) then
                        begin
                            
                            //
                            // 0-type models (with non-zero constant term)
                            // are tested for optimality of all coefficients.
                            //
                            // 1-type models (with zero constant term)
                            // are tested for optimality of non-constant terms only.
                            //
                            Inc(K);
                            Continue;
                        end;
                        F := 0;
                        FP := 0;
                        FM := 0;
                        I:=0;
                        while I<=N-1 do
                        begin
                            V := TmpWeights[M];
                            J:=0;
                            while J<=M-1 do
                            begin
                                V := V+XY[I,J]*TmpWeights[J];
                                Inc(J);
                            end;
                            F := F+AP_Sqr((V-XY[I,M])/S[I]);
                            if K<M then
                            begin
                                VV := XY[I,K];
                            end
                            else
                            begin
                                VV := 1;
                            end;
                            FP := FP+AP_Sqr((V+VV*HStep-XY[I,M])/S[I]);
                            FM := FM+AP_Sqr((V-VV*HStep-XY[I,M])/S[I]);
                            Inc(I);
                        end;
                        GROptErrors := GROptErrors or AP_FP_Greater(F,FP) or AP_FP_Greater(F,FM);
                        Inc(K);
                    end;
                    
                    //
                    // Covariance matrix test:
                    // generate random vector, project coefficients on it,
                    // compare variance of projection with estimate provided
                    // by cov.matrix
                    //
                    SetLength(TA, EstPassCount-1+1);
                    SetLength(TB, M+1);
                    SetLength(TC, M+1);
                    SetLength(XY2, N-1+1, M+1);
                    I:=0;
                    while I<=M do
                    begin
                        TB[I] := GenerateNormal(0, 1);
                        Inc(I);
                    end;
                    EPass:=0;
                    while EPass<=EstPassCount-1 do
                    begin
                        I:=0;
                        while I<=N-1 do
                        begin
                            APVMove(@XY2[I][0], 0, M-1, @XY[I][0], 0, M-1);
                            XY2[I,M] := XY0[I]+S[I]*GenerateNormal(0, 1);
                            Inc(I);
                        end;
                        if ModelType=0 then
                        begin
                            LRBuildS(XY2, S, N, M, Info, WT, AR2);
                        end
                        else
                        begin
                            LRBuildZS(XY2, S, N, M, Info, WT, AR2);
                        end;
                        if Info<>1 then
                        begin
                            TA[EPass] := 0;
                            GRConvErrors := True;
                            Exit;
                        end;
                        LRUnpack(WT, W2, TmpI);
                        V := APVDotProduct(@TB[0], 0, M, @W2[0], 0, M);
                        TA[EPass] := V;
                        Inc(EPass);
                    end;
                    CalculateMV(TA, EstPassCount, Mean, MeanS, StdDev, StdDevS);
                    I:=0;
                    while I<=M do
                    begin
                        V := 0.0;
                        for i_ := 0 to M do
                        begin
                            V := V + TB[i_]*AR.C[i_,I];
                        end;
                        TC[I] := V;
                        Inc(I);
                    end;
                    V := APVDotProduct(@TC[0], 0, M, @TB[0], 0, M);
                    GRCovErrors := GRCovErrors or AP_FP_Greater_Eq(AbsReal((Sqrt(V)-StdDev)/StdDevS),SigmaThreshold);
                    
                    //
                    // Test for the fast CV error:
                    // calculate CV error by definition (leaving out N
                    // points and recalculating solution).
                    //
                    // Test for the training set error
                    //
                    CVRMSError := 0;
                    CVAvgError := 0;
                    CVAvgRelError := 0;
                    RMSError := 0;
                    AvgError := 0;
                    AvgRelError := 0;
                    SetLength(XY2, N-2+1, M+1);
                    SetLength(S2, N-2+1);
                    I:=0;
                    while I<=N-2 do
                    begin
                        APVMove(@XY2[I][0], 0, M, @XY[I+1][0], 0, M);
                        S2[I] := S[I+1];
                        Inc(I);
                    end;
                    I:=0;
                    while I<=N-1 do
                    begin
                        
                        //
                        // Trn
                        //
                        V := APVDotProduct(@XY[I][0], 0, M-1, @TmpWeights[0], 0, M-1);
                        V := V+TmpWeights[M];
                        RMSError := RMSError+AP_Sqr(V-XY[I,M]);
                        AvgError := AvgError+AbsReal(V-XY[I,M]);
                        AvgRelError := AvgRelError+AbsReal((V-XY[I,M])/XY[I,M]);
                        
                        //
                        // CV: non-defect vectors only
                        //
                        NonDefect := True;
                        K:=0;
                        while K<=AR.NCVDefects-1 do
                        begin
                            if AR.CVDefects[K]=I then
                            begin
                                NonDefect := False;
                            end;
                            Inc(K);
                        end;
                        if NonDefect then
                        begin
                            if ModelType=0 then
                            begin
                                LRBuildS(XY2, S2, N-1, M, Info2, WT, AR2);
                            end
                            else
                            begin
                                LRBuildZS(XY2, S2, N-1, M, Info2, WT, AR2);
                            end;
                            if Info2<>1 then
                            begin
                                GRConvErrors := True;
                                Inc(I);
                                Continue;
                            end;
                            LRUnpack(WT, W2, TmpI);
                            V := APVDotProduct(@XY[I][0], 0, M-1, @W2[0], 0, M-1);
                            V := V+W2[M];
                            CVRMSError := CVRMSError+AP_Sqr(V-XY[I,M]);
                            CVAvgError := CVAvgError+AbsReal(V-XY[I,M]);
                            CVAvgRelError := CVAvgRelError+AbsReal((V-XY[I,M])/XY[I,M]);
                        end;
                        
                        //
                        // Next set
                        //
                        if I<>N-1 then
                        begin
                            APVMove(@XY2[I][0], 0, M, @XY[I][0], 0, M);
                            S2[I] := S[I];
                        end;
                        Inc(I);
                    end;
                    CVRMSError := Sqrt(CVRMSError/(N-AR.NCVDefects));
                    CVAvgError := CVAvgError/(N-AR.NCVDefects);
                    CVAvgRelError := CVAvgRelError/(N-AR.NCVDefects);
                    RMSError := Sqrt(RMSError/N);
                    AvgError := AvgError/N;
                    AvgRelError := AvgRelError/N;
                    GREstErrors := GREstErrors or AP_FP_Greater(AbsReal(Ln(AR.CVRMSError/CVRMSError)),Ln(1+1.0E-5));
                    GREstErrors := GREstErrors or AP_FP_Greater(AbsReal(Ln(AR.CVAvgError/CVAvgError)),Ln(1+1.0E-5));
                    GREstErrors := GREstErrors or AP_FP_Greater(AbsReal(Ln(AR.CVAvgRelError/CVAvgRelError)),Ln(1+1.0E-5));
                    GREstErrors := GREstErrors or AP_FP_Greater(AbsReal(Ln(AR.RMSError/RMSError)),Ln(1+1.0E-5));
                    GREstErrors := GREstErrors or AP_FP_Greater(AbsReal(Ln(AR.AvgError/AvgError)),Ln(1+1.0E-5));
                    GREstErrors := GREstErrors or AP_FP_Greater(AbsReal(Ln(AR.AvgRelError/AvgRelError)),Ln(1+1.0E-5));
                    Inc(N);
                end;
                Inc(M);
            end;
            Inc(TaskType);
        end;
        Inc(ModelType);
    end;
    
    //
    // Additional subroutines
    //
    Pass:=1;
    while Pass<=50 do
    begin
        N := 2;
        repeat
            NoiseLevel := RandomReal+0.1;
            TaskLevel := 2*RandomReal-1;
        until AP_FP_Greater(AbsReal(NoiseLevel-TaskLevel),0.05);
        SetLength(XY, 3*N-1+1, 1+1);
        I:=0;
        while I<=N-1 do
        begin
            XY[3*I+0,0] := I;
            XY[3*I+1,0] := I;
            XY[3*I+2,0] := I;
            XY[3*I+0,1] := TaskLevel-NoiseLevel;
            XY[3*I+1,1] := TaskLevel;
            XY[3*I+2,1] := TaskLevel+NoiseLevel;
            Inc(I);
        end;
        LRBuild(XY, 3*N, 1, Info, WT, AR);
        if Info=1 then
        begin
            LRUnpack(WT, TmpWeights, TmpI);
            V := LRRMSError(WT, XY, 3*N);
            GROtherErrors := GROtherErrors or AP_FP_Greater(AbsReal(V-NoiseLevel*Sqrt(AP_Double(2)/3)),Threshold);
            V := LRAvgError(WT, XY, 3*N);
            GROtherErrors := GROtherErrors or AP_FP_Greater(AbsReal(V-NoiseLevel*(AP_Double(2)/3)),Threshold);
            V := LRAvgRelError(WT, XY, 3*N);
            VV := (AbsReal(NoiseLevel/(TaskLevel-NoiseLevel))+AbsReal(NoiseLevel/(TaskLevel+NoiseLevel)))/3;
            GROtherErrors := GROtherErrors or AP_FP_Greater(AbsReal(V-VV),Threshold*VV);
        end
        else
        begin
            GROtherErrors := True;
        end;
        I:=0;
        while I<=N-1 do
        begin
            XY[3*I+0,0] := I;
            XY[3*I+1,0] := I;
            XY[3*I+2,0] := I;
            XY[3*I+0,1] := -NoiseLevel;
            XY[3*I+1,1] := 0;
            XY[3*I+2,1] := +NoiseLevel;
            Inc(I);
        end;
        LRBuild(XY, 3*N, 1, Info, WT, AR);
        if Info=1 then
        begin
            LRUnpack(WT, TmpWeights, TmpI);
            V := LRAvgRelError(WT, XY, 3*N);
            GROtherErrors := GROtherErrors or AP_FP_Greater(AbsReal(V-1),Threshold);
        end
        else
        begin
            GROtherErrors := True;
        end;
        Inc(Pass);
    end;
    Pass:=1;
    while Pass<=10 do
    begin
        M := 1+RandomInteger(5);
        N := 10+RandomInteger(10);
        SetLength(XY, N-1+1, M+1);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M do
            begin
                XY[I,J] := 2*RandomReal-1;
                Inc(J);
            end;
            Inc(I);
        end;
        LRBuild(XY, N, M, Info, W, AR);
        if Info<0 then
        begin
            GROtherErrors := True;
            Break;
        end;
        SetLength(X1, M-1+1);
        SetLength(X2, M-1+1);
        
        //
        // Same inputs on original leads to same outputs
        // on copy created using LRCopy
        //
        UnsetLR(WT);
        LRCopy(W, WT);
        I:=0;
        while I<=M-1 do
        begin
            X1[I] := 2*RandomReal-1;
            X2[I] := X1[I];
            Inc(I);
        end;
        Y1 := LRProcess(W, X1);
        Y2 := LRProcess(WT, X2);
        AllSame := AP_FP_Eq(Y1,Y2);
        GROtherErrors := GROtherErrors or  not AllSame;
        
        //
        // Same inputs on original leads to same outputs
        // on copy created using LRSerialize
        //
        UnsetLR(WT);
        SetLength(RA, 0+1);
        RA[0] := 0;
        RLen := 0;
        LRSerialize(W, RA, RLen);
        SetLength(RA2, RLen-1+1);
        I:=0;
        while I<=RLen-1 do
        begin
            RA2[I] := RA[I];
            Inc(I);
        end;
        LRUnserialize(RA2, WT);
        I:=0;
        while I<=M-1 do
        begin
            X1[I] := 2*RandomReal-1;
            X2[I] := X1[I];
            Inc(I);
        end;
        Y1 := LRProcess(W, X1);
        Y2 := LRProcess(WT, X2);
        AllSame := AP_FP_Eq(Y1,Y2);
        GROtherErrors := GROtherErrors or  not AllSame;
        Inc(Pass);
    end;
    
    //
    // TODO: Degenerate tests (when design matrix and right part are zero)
    //
    
    //
    // Final report
    //
    WasErrors := SLErrors or SLCErrors or GROptErrors or GRCovErrors or GREstErrors or GROtherErrors or GRConvErrors;
    if  not Silent then
    begin
        Write(Format('REGRESSION TEST'#13#10'',[]));
        Write(Format('STRAIGHT LINE REGRESSION:                ',[]));
        if  not SLErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('STRAIGHT LINE REGRESSION CONVERGENCE:    ',[]));
        if  not SLCErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('GENERAL LINEAR REGRESSION:               ',[]));
        if  not (GROptErrors or GRCovErrors or GREstErrors or GROtherErrors or GRConvErrors) then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('* OPTIMALITY:                            ',[]));
        if  not GROptErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('* COV. MATRIX:                           ',[]));
        if  not GRCovErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('* ERROR ESTIMATES:                       ',[]));
        if  not GREstErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('* CONVERGENCE:                           ',[]));
        if  not GRConvErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('* OTHER SUBROUTINES:                     ',[]));
        if  not GROtherErrors then
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
Task generation. Meaningless task, just random numbers.
*************************************************************************)
procedure GenerateRandomTask(XL : Double;
     XR : Double;
     RandomX : Boolean;
     YMin : Double;
     YMax : Double;
     SMin : Double;
     SMax : Double;
     N : AlglibInteger;
     var XY : TReal2DArray;
     var S : TReal1DArray);
var
    I : AlglibInteger;
begin
    SetLength(XY, N-1+1, 1+1);
    SetLength(S, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        if RandomX then
        begin
            XY[I,0] := XL+(XR-XL)*RandomReal;
        end
        else
        begin
            XY[I,0] := XL+(XR-XL)*I/(N-1);
        end;
        XY[I,1] := YMin+(YMax-YMin)*RandomReal;
        S[I] := SMin+(SMax-SMin)*RandomReal;
        Inc(I);
    end;
end;


(*************************************************************************
Task generation.
*************************************************************************)
procedure GenerateTask(A : Double;
     B : Double;
     XL : Double;
     XR : Double;
     RandomX : Boolean;
     SMin : Double;
     SMax : Double;
     N : AlglibInteger;
     var XY : TReal2DArray;
     var S : TReal1DArray);
var
    I : AlglibInteger;
begin
    SetLength(XY, N-1+1, 1+1);
    SetLength(S, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        if RandomX then
        begin
            XY[I,0] := XL+(XR-XL)*RandomReal;
        end
        else
        begin
            XY[I,0] := XL+(XR-XL)*I/(N-1);
        end;
        S[I] := SMin+(SMax-SMin)*RandomReal;
        XY[I,1] := A+B*XY[I,0]+GenerateNormal(0, S[I]);
        Inc(I);
    end;
end;


(*************************************************************************
Task generation.
y[i] are filled based on A, B, X[I], S[I]
*************************************************************************)
procedure FillTaskWithY(A : Double;
     B : Double;
     N : AlglibInteger;
     var XY : TReal2DArray;
     var S : TReal1DArray);
var
    I : AlglibInteger;
begin
    I:=0;
    while I<=N-1 do
    begin
        XY[I,1] := A+B*XY[I,0]+GenerateNormal(0, S[I]);
        Inc(I);
    end;
end;


(*************************************************************************
Normal random numbers
*************************************************************************)
function GenerateNormal(Mean : Double; Sigma : Double):Double;
var
    U : Double;
    V : Double;
    S : Double;
    Sum : Double;
begin
    Result := Mean;
    while True do
    begin
        u := (2*RandomInteger(2)-1)*RandomReal;
        v := (2*RandomInteger(2)-1)*RandomReal;
        sum := u*u+v*v;
        if AP_FP_Less(sum,1) and AP_FP_Greater(sum,0) then
        begin
            sum := sqrt(-2*ln(sum)/sum);
            Result := Sigma*u*sum+Mean;
            Exit;
        end;
    end;
end;


(*************************************************************************
Moments estimates and their errors
*************************************************************************)
procedure CalculateMV(const X : TReal1DArray;
     N : AlglibInteger;
     var Mean : Double;
     var MeanS : Double;
     var StdDev : Double;
     var StdDevS : Double);
var
    I : AlglibInteger;
    V : Double;
    V1 : Double;
    V2 : Double;
    Variance : Double;
begin
    Mean := 0;
    MeanS := 1;
    StdDev := 0;
    StdDevS := 1;
    Variance := 0;
    if N<=1 then
    begin
        Exit;
    end;
    
    //
    // Mean
    //
    I:=0;
    while I<=N-1 do
    begin
        Mean := Mean+X[I];
        Inc(I);
    end;
    Mean := Mean/N;
    
    //
    // Variance (using corrected two-pass algorithm)
    //
    if N<>1 then
    begin
        V1 := 0;
        I:=0;
        while I<=N-1 do
        begin
            V1 := V1+AP_Sqr(X[I]-Mean);
            Inc(I);
        end;
        V2 := 0;
        I:=0;
        while I<=N-1 do
        begin
            V2 := V2+(X[I]-Mean);
            Inc(I);
        end;
        V2 := AP_Sqr(V2)/N;
        Variance := (V1-V2)/(N-1);
        if AP_FP_Less(Variance,0) then
        begin
            Variance := 0;
        end;
        StdDev := Sqrt(Variance);
    end;
    
    //
    // Errors
    //
    MeanS := StdDev/Sqrt(N);
    StdDevS := StdDev*Sqrt(2)/Sqrt(N-1);
end;


(*************************************************************************
Unsets LR
*************************************************************************)
procedure UnsetLR(var LR : LinearModel);
var
    XY : TReal2DArray;
    Info : AlglibInteger;
    Rep : LRReport;
    I : AlglibInteger;
begin
    SetLength(XY, 5+1, 1+1);
    I:=0;
    while I<=5 do
    begin
        XY[I,0] := 0;
        XY[I,1] := 0;
        Inc(I);
    end;
    LRBuild(XY, 6, 1, Info, LR, Rep);
    Assert(Info>0);
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testregressunit_test_silent():Boolean;
begin
    Result := TestLinRegression(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testregressunit_test():Boolean;
begin
    Result := TestLinRegression(False);
end;


end.