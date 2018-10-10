unit testratinterpolation;
interface
uses Math, Sysutils, Ap, tsort, ratinterpolation, blas, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, hblas, sblas, ortfac, rotations, bdsvd, svd, xblas, densesolver, linmin, minlbfgs, minlm, lsfit, ratint;

function TestRI(Silent : Boolean):Boolean;
function testratinterpolation_test_silent():Boolean;
function testratinterpolation_test():Boolean;

implementation

procedure PolDiff2(const X : TReal1DArray;
     F : TReal1DArray;
     N : AlglibInteger;
     T : Double;
     var P : Double;
     var DP : Double;
     var D2P : Double);forward;
procedure BRCUnset(var B : BarycentricInterpolant);forward;
function Is1DSolution(N : AlglibInteger;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     C : Double):Boolean;forward;


function TestRI(Silent : Boolean):Boolean;
var
    WasErrors : Boolean;
    BCErrors : Boolean;
    NPErrors : Boolean;
    FitErrors : Boolean;
    Threshold : Double;
    LipschitzTol : Double;
    MaxN : AlglibInteger;
    PassCount : AlglibInteger;
    B1 : BarycentricInterpolant;
    B2 : BarycentricInterpolant;
    X : TReal1DArray;
    X2 : TReal1DArray;
    Y : TReal1DArray;
    Y2 : TReal1DArray;
    W : TReal1DArray;
    W2 : TReal1DArray;
    XC : TReal1DArray;
    YC : TReal1DArray;
    DC : TInteger1DArray;
    H : Double;
    S1 : Double;
    S2 : Double;
    BSame : Boolean;
    N : AlglibInteger;
    M : AlglibInteger;
    N2 : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    D : AlglibInteger;
    Pass : AlglibInteger;
    Err : Double;
    MaxErr : Double;
    T : Double;
    A : Double;
    B : Double;
    S : Double;
    V : Double;
    V0 : Double;
    V1 : Double;
    V2 : Double;
    V3 : Double;
    D0 : Double;
    D1 : Double;
    D2 : Double;
    Info : AlglibInteger;
    Info2 : AlglibInteger;
    XMin : Double;
    XMax : Double;
    RefRMS : Double;
    RefAvg : Double;
    RefAvgRel : Double;
    RefMax : Double;
    RA : TReal1DArray;
    RA2 : TReal1DArray;
    RALen : AlglibInteger;
    Rep : BarycentricFitReport;
    Rep2 : BarycentricFitReport;
    B3 : BarycentricInterpolant;
    B4 : BarycentricInterpolant;
begin
    NPErrors := False;
    BCErrors := False;
    FitErrors := False;
    WasErrors := False;
    
    //
    // PassCount        number of repeated passes
    // Threshold        error tolerance
    // LipschitzTol     Lipschitz constant increase allowed
    //                  when calculating constant on a twice denser grid
    //
    PassCount := 5;
    MaxN := 15;
    Threshold := 1000000*MachineEpsilon;
    LipschitzTol := 1.3;
    
    //
    // Basic barycentric functions
    //
    N:=1;
    while N<=10 do
    begin
        
        //
        // randomized tests
        //
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // generate weights from polynomial interpolation
            //
            V0 := 1+0.4*RandomReal-0.2;
            V1 := 2*RandomReal-1;
            V2 := 2*RandomReal-1;
            V3 := 2*RandomReal-1;
            SetLength(X, N);
            SetLength(Y, N);
            SetLength(W, N);
            I:=0;
            while I<=N-1 do
            begin
                if N=1 then
                begin
                    X[I] := 0;
                end
                else
                begin
                    X[I] := V0*Cos(I*Pi/(N-1));
                end;
                Y[I] := Sin(V1*X[I])+Cos(V2*X[I])+Exp(V3*X[I]);
                Inc(I);
            end;
            J:=0;
            while J<=N-1 do
            begin
                W[J] := 1;
                K:=0;
                while K<=N-1 do
                begin
                    if K<>J then
                    begin
                        W[J] := W[J]/(X[J]-X[K]);
                    end;
                    Inc(K);
                end;
                Inc(J);
            end;
            BarycentricBuildXYW(X, Y, W, N, B1);
            
            //
            // unpack, then pack again and compare
            //
            BRCUnset(B2);
            BarycentricUnpack(B1, N2, X2, Y2, W2);
            BCErrors := BCErrors or (N2<>N);
            BarycentricBuildXYW(X2, Y2, W2, N2, B2);
            T := 2*RandomReal-1;
            BCErrors := BCErrors or AP_FP_Greater(AbsReal(BarycentricCalc(B1, T)-BarycentricCalc(B2, T)),Threshold);
            
            //
            // serialize, unserialize, compare
            //
            BRCUnset(B2);
            BarycentricSerialize(B1, RA, RALen);
            SetLength(RA2, RALen);
            APVMove(@RA2[0], 0, RALen-1, @RA[0], 0, RALen-1);
            BarycentricUnserialize(RA2, B2);
            T := 2*RandomReal-1;
            BCErrors := BCErrors or AP_FP_Greater(AbsReal(BarycentricCalc(B1, T)-BarycentricCalc(B2, T)),Threshold);
            
            //
            // copy, compare
            //
            BRCUnset(B2);
            BarycentricCopy(B1, B2);
            T := 2*RandomReal-1;
            BCErrors := BCErrors or AP_FP_Greater(AbsReal(BarycentricCalc(B1, T)-BarycentricCalc(B2, T)),Threshold);
            
            //
            // test interpolation properties
            //
            I:=0;
            while I<=N-1 do
            begin
                
                //
                // test interpolation at nodes
                //
                BCErrors := BCErrors or AP_FP_Greater(AbsReal(BarycentricCalc(B1, X[I])-Y[I]),Threshold*AbsReal(Y[I]));
                
                //
                // compare with polynomial interpolation
                //
                T := 2*RandomReal-1;
                PolDiff2(X, Y, N, T, V0, V1, V2);
                BCErrors := BCErrors or AP_FP_Greater(AbsReal(BarycentricCalc(B1, T)-V0),Threshold*Max(AbsReal(V0), 1));
                
                //
                // test continuity between nodes
                // calculate Lipschitz constant on two grids -
                // dense and even more dense. If Lipschitz constant
                // on a denser grid is significantly increased,
                // continuity test is failed
                //
                T := 3.0;
                K := 100;
                S1 := 0;
                J:=0;
                while J<=K-1 do
                begin
                    V1 := X[I]+(T-X[I])*J/K;
                    V2 := X[I]+(T-X[I])*(J+1)/K;
                    S1 := Max(S1, AbsReal(BarycentricCalc(B1, V2)-BarycentricCalc(B1, V1))/AbsReal(V2-V1));
                    Inc(J);
                end;
                K := 2*K;
                S2 := 0;
                J:=0;
                while J<=K-1 do
                begin
                    V1 := X[I]+(T-X[I])*J/K;
                    V2 := X[I]+(T-X[I])*(J+1)/K;
                    S2 := Max(S2, AbsReal(BarycentricCalc(B1, V2)-BarycentricCalc(B1, V1))/AbsReal(V2-V1));
                    Inc(J);
                end;
                BCErrors := BCErrors or AP_FP_Greater(S2,LipschitzTol*S1) and AP_FP_Greater(S1,Threshold*K);
                Inc(I);
            end;
            
            //
            // test differentiation properties
            //
            I:=0;
            while I<=N-1 do
            begin
                T := 2*RandomReal-1;
                PolDiff2(X, Y, N, T, V0, V1, V2);
                D0 := 0;
                D1 := 0;
                D2 := 0;
                BarycentricDiff1(B1, T, D0, D1);
                BCErrors := BCErrors or AP_FP_Greater(AbsReal(V0-D0),Threshold*Max(AbsReal(V0), 1));
                BCErrors := BCErrors or AP_FP_Greater(AbsReal(V1-D1),Threshold*Max(AbsReal(V1), 1));
                D0 := 0;
                D1 := 0;
                D2 := 0;
                BarycentricDiff2(B1, T, D0, D1, D2);
                BCErrors := BCErrors or AP_FP_Greater(AbsReal(V0-D0),Threshold*Max(AbsReal(V0), 1));
                BCErrors := BCErrors or AP_FP_Greater(AbsReal(V1-D1),Threshold*Max(AbsReal(V1), 1));
                BCErrors := BCErrors or AP_FP_Greater(AbsReal(V2-D2),Sqrt(Threshold)*Max(AbsReal(V2), 1));
                Inc(I);
            end;
            
            //
            // test linear translation
            //
            T := 2*RandomReal-1;
            A := 2*RandomReal-1;
            B := 2*RandomReal-1;
            BRCUnset(B2);
            BarycentricCopy(B1, B2);
            BarycentricLinTransX(B2, A, B);
            BCErrors := BCErrors or AP_FP_Greater(AbsReal(BarycentricCalc(B1, A*T+B)-BarycentricCalc(B2, T)),Threshold);
            A := 0;
            B := 2*RandomReal-1;
            BRCUnset(B2);
            BarycentricCopy(B1, B2);
            BarycentricLinTransX(B2, A, B);
            BCErrors := BCErrors or AP_FP_Greater(AbsReal(BarycentricCalc(B1, A*T+B)-BarycentricCalc(B2, T)),Threshold);
            A := 2*RandomReal-1;
            B := 2*RandomReal-1;
            BRCUnset(B2);
            BarycentricCopy(B1, B2);
            BarycentricLinTransY(B2, A, B);
            BCErrors := BCErrors or AP_FP_Greater(AbsReal(A*BarycentricCalc(B1, T)+B-BarycentricCalc(B2, T)),Threshold);
            Inc(Pass);
        end;
        Inc(N);
    end;
    Pass:=0;
    while Pass<=3 do
    begin
        
        //
        // Crash-test: small numbers, large numbers
        //
        SetLength(X, 4);
        SetLength(Y, 4);
        SetLength(W, 4);
        H := 1;
        if Pass mod 2=0 then
        begin
            H := 100*MinRealNumber;
        end;
        if Pass mod 2=1 then
        begin
            H := 0.01*MaxRealNumber;
        end;
        X[0] := 0*H;
        X[1] := 1*H;
        X[2] := 2*H;
        X[3] := 3*H;
        Y[0] := 0*H;
        Y[1] := 1*H;
        Y[2] := 2*H;
        Y[3] := 3*H;
        W[0] := -1/(X[1]-X[0]);
        W[1] := +1*(1/(X[1]-X[0])+1/(X[2]-X[1]));
        W[2] := -1*(1/(X[2]-X[1])+1/(X[3]-X[2]));
        W[3] := +1/(X[3]-X[2]);
        if Pass div 2=0 then
        begin
            V0 := 0;
        end;
        if Pass div 2=1 then
        begin
            V0 := 0.6*H;
        end;
        BarycentricBuildXYW(X, Y, W, 4, B1);
        T := BarycentricCalc(B1, V0);
        D0 := 0;
        D1 := 0;
        D2 := 0;
        BarycentricDiff1(B1, V0, D0, D1);
        BCErrors := BCErrors or AP_FP_Greater(AbsReal(T-V0),Threshold*V0);
        BCErrors := BCErrors or AP_FP_Greater(AbsReal(D0-V0),Threshold*V0);
        BCErrors := BCErrors or AP_FP_Greater(AbsReal(D1-1),1000*Threshold);
        Inc(Pass);
    end;
    
    //
    // crash test: large abscissas, small argument
    //
    // test for errors in D0 is not very strict
    // because renormalization used in Diff1()
    // destroys part of precision.
    //
    SetLength(X, 4);
    SetLength(Y, 4);
    SetLength(W, 4);
    H := 0.01*MaxRealNumber;
    X[0] := 0*H;
    X[1] := 1*H;
    X[2] := 2*H;
    X[3] := 3*H;
    Y[0] := 0*H;
    Y[1] := 1*H;
    Y[2] := 2*H;
    Y[3] := 3*H;
    W[0] := -1/(X[1]-X[0]);
    W[1] := +1*(1/(X[1]-X[0])+1/(X[2]-X[1]));
    W[2] := -1*(1/(X[2]-X[1])+1/(X[3]-X[2]));
    W[3] := +1/(X[3]-X[2]);
    V0 := 100*MinRealNumber;
    BarycentricBuildXYW(X, Y, W, 4, B1);
    T := BarycentricCalc(B1, V0);
    D0 := 0;
    D1 := 0;
    D2 := 0;
    BarycentricDiff1(B1, V0, D0, D1);
    BCErrors := BCErrors or AP_FP_Greater(AbsReal(T),V0*(1+Threshold));
    BCErrors := BCErrors or AP_FP_Greater(AbsReal(D0),V0*(1+Threshold));
    BCErrors := BCErrors or AP_FP_Greater(AbsReal(D1-1),1000*Threshold);
    
    //
    // crash test: test safe barycentric formula
    //
    SetLength(X, 4);
    SetLength(Y, 4);
    SetLength(W, 4);
    H := 2*MinRealNumber;
    X[0] := 0*H;
    X[1] := 1*H;
    X[2] := 2*H;
    X[3] := 3*H;
    Y[0] := 0*H;
    Y[1] := 1*H;
    Y[2] := 2*H;
    Y[3] := 3*H;
    W[0] := -1/(X[1]-X[0]);
    W[1] := +1*(1/(X[1]-X[0])+1/(X[2]-X[1]));
    W[2] := -1*(1/(X[2]-X[1])+1/(X[3]-X[2]));
    W[3] := +1/(X[3]-X[2]);
    V0 := MinRealNumber;
    BarycentricBuildXYW(X, Y, W, 4, B1);
    T := BarycentricCalc(B1, V0);
    BCErrors := BCErrors or AP_FP_Greater(AbsReal(T-V0)/V0,Threshold);
    
    //
    // Testing "No Poles" interpolation
    //
    MaxErr := 0;
    Pass:=1;
    while Pass<=PassCount-1 do
    begin
        SetLength(X, 1);
        SetLength(Y, 1);
        X[0] := 2*RandomReal-1;
        Y[0] := 2*RandomReal-1;
        BarycentricBuildFloaterHormann(X, Y, 1, 1, B1);
        MaxErr := Max(MaxErr, AbsReal(BarycentricCalc(B1, 2*RandomReal-1)-Y[0]));
        Inc(Pass);
    end;
    N:=2;
    while N<=10 do
    begin
        
        //
        // compare interpolant built by subroutine
        // with interpolant built by hands
        //
        SetLength(X, N);
        SetLength(Y, N);
        SetLength(W, N);
        SetLength(W2, N);
        
        //
        // D=1, non-equidistant nodes
        //
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // Initialize X, Y, W
            //
            A := -1-1*RandomReal;
            B := +1+1*RandomReal;
            I:=0;
            while I<=N-1 do
            begin
                X[I] := ArcTan((B-A)*I/(N-1)+A);
                Inc(I);
            end;
            I:=0;
            while I<=N-1 do
            begin
                Y[I] := 2*RandomReal-1;
                Inc(I);
            end;
            W[0] := -1/(X[1]-X[0]);
            S := 1;
            I:=1;
            while I<=N-2 do
            begin
                W[I] := S*(1/(X[I]-X[I-1])+1/(X[I+1]-X[I]));
                S := -S;
                Inc(I);
            end;
            W[N-1] := S/(X[N-1]-X[N-2]);
            I:=0;
            while I<=N-1 do
            begin
                K := RandomInteger(N);
                if K<>I then
                begin
                    T := X[I];
                    X[I] := X[K];
                    X[K] := T;
                    T := Y[I];
                    Y[I] := Y[K];
                    Y[K] := T;
                    T := W[I];
                    W[I] := W[K];
                    W[K] := T;
                end;
                Inc(I);
            end;
            
            //
            // Build and test
            //
            BarycentricBuildFloaterHormann(X, Y, N, 1, B1);
            BarycentricBuildXYW(X, Y, W, N, B2);
            I:=1;
            while I<=2*N do
            begin
                T := A+(B-A)*RandomReal;
                MaxErr := Max(MaxErr, AbsReal(BarycentricCalc(B1, T)-BarycentricCalc(B2, T)));
                Inc(I);
            end;
            Inc(Pass);
        end;
        
        //
        // D = 0, 1, 2. Equidistant nodes.
        //
        D:=0;
        while D<=2 do
        begin
            Pass:=1;
            while Pass<=PassCount do
            begin
                
                //
                // Skip incorrect (N,D) pairs
                //
                if N<2*D then
                begin
                    Inc(Pass);
                    Continue;
                end;
                
                //
                // Initialize X, Y, W
                //
                A := -1-1*RandomReal;
                B := +1+1*RandomReal;
                I:=0;
                while I<=N-1 do
                begin
                    X[I] := (B-A)*I/(N-1)+A;
                    Inc(I);
                end;
                I:=0;
                while I<=N-1 do
                begin
                    Y[I] := 2*RandomReal-1;
                    Inc(I);
                end;
                S := 1;
                if D=0 then
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        W[I] := S;
                        S := -S;
                        Inc(I);
                    end;
                end;
                if D=1 then
                begin
                    W[0] := -S;
                    I:=1;
                    while I<=N-2 do
                    begin
                        W[I] := 2*S;
                        S := -S;
                        Inc(I);
                    end;
                    W[N-1] := S;
                end;
                if D=2 then
                begin
                    W[0] := S;
                    W[1] := -3*S;
                    I:=2;
                    while I<=N-3 do
                    begin
                        W[I] := 4*S;
                        S := -S;
                        Inc(I);
                    end;
                    W[N-2] := 3*S;
                    W[N-1] := -S;
                end;
                
                //
                // Mix
                //
                I:=0;
                while I<=N-1 do
                begin
                    K := RandomInteger(N);
                    if K<>I then
                    begin
                        T := X[I];
                        X[I] := X[K];
                        X[K] := T;
                        T := Y[I];
                        Y[I] := Y[K];
                        Y[K] := T;
                        T := W[I];
                        W[I] := W[K];
                        W[K] := T;
                    end;
                    Inc(I);
                end;
                
                //
                // Build and test
                //
                BarycentricBuildFloaterHormann(X, Y, N, D, B1);
                BarycentricBuildXYW(X, Y, W, N, B2);
                I:=1;
                while I<=2*N do
                begin
                    T := A+(B-A)*RandomReal;
                    MaxErr := Max(MaxErr, AbsReal(BarycentricCalc(B1, T)-BarycentricCalc(B2, T)));
                    Inc(I);
                end;
                Inc(Pass);
            end;
            Inc(D);
        end;
        Inc(N);
    end;
    if AP_FP_Greater(MaxErr,Threshold) then
    begin
        NPErrors := True;
    end;
    
    //
    // Test rational fitting:
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        N:=2;
        while N<=MaxN do
        begin
            
            //
            // N=M+K fitting (i.e. interpolation)
            //
            K:=0;
            while K<=N-1 do
            begin
                SetLength(X, N-K);
                SetLength(Y, N-K);
                SetLength(W, N-K);
                if K>0 then
                begin
                    SetLength(XC, K);
                    SetLength(YC, K);
                    SetLength(DC, K);
                end;
                I:=0;
                while I<=N-K-1 do
                begin
                    X[I] := AP_Double(I)/(N-1);
                    Y[I] := 2*RandomReal-1;
                    W[I] := 1+RandomReal;
                    Inc(I);
                end;
                I:=0;
                while I<=K-1 do
                begin
                    XC[I] := AP_Double((N-K+I))/(N-1);
                    YC[I] := 2*RandomReal-1;
                    DC[I] := 0;
                    Inc(I);
                end;
                BarycentricFitFloaterHormannWC(X, Y, W, N-K, XC, YC, DC, K, N, Info, B1, Rep);
                if Info<=0 then
                begin
                    FitErrors := True;
                end
                else
                begin
                    I:=0;
                    while I<=N-K-1 do
                    begin
                        FitErrors := FitErrors or AP_FP_Greater(AbsReal(BarycentricCalc(B1, X[I])-Y[I]),Threshold);
                        Inc(I);
                    end;
                    I:=0;
                    while I<=K-1 do
                    begin
                        FitErrors := FitErrors or AP_FP_Greater(AbsReal(BarycentricCalc(B1, XC[I])-YC[I]),Threshold);
                        Inc(I);
                    end;
                end;
                Inc(K);
            end;
            
            //
            // Testing constraints on derivatives:
            // * several M's are tried
            // * several K's are tried - 1, 2.
            // * constraints at the ends of the interval
            //
            M:=3;
            while M<=5 do
            begin
                K:=1;
                while K<=2 do
                begin
                    SetLength(X, N);
                    SetLength(Y, N);
                    SetLength(W, N);
                    SetLength(XC, 2);
                    SetLength(YC, 2);
                    SetLength(DC, 2);
                    I:=0;
                    while I<=N-1 do
                    begin
                        X[I] := 2*RandomReal-1;
                        Y[I] := 2*RandomReal-1;
                        W[I] := 1+RandomReal;
                        Inc(I);
                    end;
                    XC[0] := -1;
                    YC[0] := 2*RandomReal-1;
                    DC[0] := 0;
                    XC[1] := +1;
                    YC[1] := 2*RandomReal-1;
                    DC[1] := 0;
                    BarycentricFitFloaterHormannWC(X, Y, W, N, XC, YC, DC, K, M, Info, B1, Rep);
                    if Info<=0 then
                    begin
                        FitErrors := True;
                    end
                    else
                    begin
                        I:=0;
                        while I<=K-1 do
                        begin
                            BarycentricDiff1(B1, XC[I], V0, V1);
                            FitErrors := FitErrors or AP_FP_Greater(AbsReal(V0-YC[I]),Threshold);
                            Inc(I);
                        end;
                    end;
                    Inc(K);
                end;
                Inc(M);
            end;
            Inc(N);
        end;
        Inc(Pass);
    end;
    M:=2;
    while M<=8 do
    begin
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // General fitting
            //
            // interpolating function through M nodes should have
            // greater RMS error than fitting it through the same M nodes
            //
            N := 100;
            SetLength(X2, N);
            SetLength(Y2, N);
            SetLength(W2, N);
            XMin := MaxRealNumber;
            XMax := -MaxRealNumber;
            I:=0;
            while I<=N-1 do
            begin
                X2[I] := 2*Pi*RandomReal;
                Y2[I] := Sin(X2[I]);
                W2[I] := 1;
                XMin := Min(XMin, X2[I]);
                XMax := Max(XMax, X2[I]);
                Inc(I);
            end;
            SetLength(X, M);
            SetLength(Y, M);
            I:=0;
            while I<=M-1 do
            begin
                X[I] := XMin+(XMax-XMin)*I/(M-1);
                Y[I] := Sin(X[I]);
                Inc(I);
            end;
            BarycentricBuildFloaterHormann(X, Y, M, 3, B1);
            BarycentricFitFloaterHormannWC(X2, Y2, W2, N, XC, YC, DC, 0, M, Info, B2, Rep);
            if Info<=0 then
            begin
                FitErrors := True;
            end
            else
            begin
                
                //
                // calculate B1 (interpolant) RMS error, compare with B2 error
                //
                V1 := 0;
                V2 := 0;
                I:=0;
                while I<=N-1 do
                begin
                    V1 := V1+AP_Sqr(BarycentricCalc(B1, X2[I])-Y2[I]);
                    V2 := V2+AP_Sqr(BarycentricCalc(B2, X2[I])-Y2[I]);
                    Inc(I);
                end;
                V1 := Sqrt(V1/N);
                V2 := Sqrt(V2/N);
                FitErrors := FitErrors or AP_FP_Greater(V2,V1);
                FitErrors := FitErrors or AP_FP_Greater(AbsReal(V2-Rep.RMSError),Threshold);
            end;
            
            //
            // compare weighted and non-weighted
            //
            N := 20;
            SetLength(X, N);
            SetLength(Y, N);
            SetLength(W, N);
            I:=0;
            while I<=N-1 do
            begin
                X[I] := 2*RandomReal-1;
                Y[I] := 2*RandomReal-1;
                W[I] := 1;
                Inc(I);
            end;
            BarycentricFitFloaterHormannWC(X, Y, W, N, XC, YC, DC, 0, M, Info, B1, Rep);
            BarycentricFitFloaterHormann(X, Y, N, M, Info2, B2, Rep2);
            if (Info<=0) or (Info2<=0) then
            begin
                FitErrors := True;
            end
            else
            begin
                
                //
                // calculate B1 (interpolant), compare with B2
                // compare RMS errors
                //
                T := 2*RandomReal-1;
                V1 := BarycentricCalc(B1, T);
                V2 := BarycentricCalc(B2, T);
                FitErrors := FitErrors or AP_FP_Neq(V2,V1);
                FitErrors := FitErrors or AP_FP_Neq(Rep.RMSError,Rep2.RMSError);
                FitErrors := FitErrors or AP_FP_Neq(Rep.AvgError,Rep2.AvgError);
                FitErrors := FitErrors or AP_FP_Neq(Rep.AvgRelError,Rep2.AvgRelError);
                FitErrors := FitErrors or AP_FP_Neq(Rep.MaxError,Rep2.MaxError);
            end;
            Inc(Pass);
        end;
        Inc(M);
    end;
    Pass:=1;
    while Pass<=PassCount do
    begin
        Assert(PassCount>=2, 'PassCount should be 2 or greater!');
        
        //
        // solve simple task (all X[] are the same, Y[] are specially
        // calculated to ensure simple form of all types of errors)
        // and check correctness of the errors calculated by subroutines
        //
        // First pass is done with zero Y[], other passes - with random Y[].
        // It should test both ability to correctly calculate errors and
        // ability to not fail while working with zeros :)
        //
        N := 4;
        if Pass=1 then
        begin
            V1 := 0;
            V2 := 0;
            V := 0;
        end
        else
        begin
            V1 := RandomReal;
            V2 := RandomReal;
            V := 1+RandomReal;
        end;
        SetLength(X, 4);
        SetLength(Y, 4);
        SetLength(W, 4);
        X[0] := 0;
        Y[0] := V-V2;
        W[0] := 1;
        X[1] := 0;
        Y[1] := V-V1;
        W[1] := 1;
        X[2] := 0;
        Y[2] := V+V1;
        W[2] := 1;
        X[3] := 0;
        Y[3] := V+V2;
        W[3] := 1;
        RefRms := Sqrt((AP_Sqr(V1)+AP_Sqr(V2))/2);
        RefAvg := (AbsReal(V1)+AbsReal(V2))/2;
        if Pass=1 then
        begin
            RefAvgRel := 0;
        end
        else
        begin
            RefAvgRel := 0.25*(AbsReal(V2)/AbsReal(V-V2)+AbsReal(V1)/AbsReal(V-V1)+AbsReal(V1)/AbsReal(V+V1)+AbsReal(V2)/AbsReal(V+V2));
        end;
        RefMax := Max(V1, V2);
        
        //
        // Test errors correctness
        //
        BarycentricFitFloaterHormann(X, Y, 4, 2, Info, B1, Rep);
        if Info<=0 then
        begin
            FitErrors := True;
        end
        else
        begin
            S := BarycentricCalc(B1, 0);
            FitErrors := FitErrors or AP_FP_Greater(AbsReal(S-V),Threshold);
            FitErrors := FitErrors or AP_FP_Greater(AbsReal(Rep.RMSError-RefRMS),Threshold);
            FitErrors := FitErrors or AP_FP_Greater(AbsReal(Rep.AvgError-RefAvg),Threshold);
            FitErrors := FitErrors or AP_FP_Greater(AbsReal(Rep.AvgRelError-RefAvgRel),Threshold);
            FitErrors := FitErrors or AP_FP_Greater(AbsReal(Rep.MaxError-RefMax),Threshold);
        end;
        Inc(Pass);
    end;
    
    //
    // report
    //
    WasErrors := BCErrors or NPErrors or FitErrors;
    if  not Silent then
    begin
        Write(Format('TESTING RATIONAL INTERPOLATION'#13#10'',[]));
        Write(Format('BASIC BARYCENTRIC FUNCTIONS:             ',[]));
        if BCErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('FLOATER-HORMANN:                         ',[]));
        if NPErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('RATIONAL FITTING:                        ',[]));
        if FitErrors then
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


procedure PolDiff2(const X : TReal1DArray;
     F : TReal1DArray;
     N : AlglibInteger;
     T : Double;
     var P : Double;
     var DP : Double;
     var D2P : Double);
var
    M : AlglibInteger;
    I : AlglibInteger;
    DF : TReal1DArray;
    D2F : TReal1DArray;
begin
    F := DynamicArrayCopy(F);
    N := N-1;
    SetLength(DF, N+1);
    SetLength(D2F, N+1);
    I:=0;
    while I<=N do
    begin
        D2F[I] := 0;
        DF[I] := 0;
        Inc(I);
    end;
    M:=1;
    while M<=N do
    begin
        I:=0;
        while I<=N-M do
        begin
            D2F[I] := ((T-X[I+M])*D2F[I]+(X[I]-T)*D2F[I+1]+2*DF[I]-2*DF[I+1])/(X[I]-X[I+M]);
            DF[I] := ((T-X[I+M])*DF[I]+F[I]+(X[I]-T)*DF[I+1]-F[I+1])/(X[I]-X[I+M]);
            F[I] := ((T-X[I+M])*F[I]+(X[I]-T)*F[I+1])/(X[I]-X[I+M]);
            Inc(I);
        end;
        Inc(M);
    end;
    P := F[0];
    DP := DF[0];
    D2P := D2F[0];
end;


procedure BRCUnset(var B : BarycentricInterpolant);
var
    X : TReal1DArray;
    Y : TReal1DArray;
    W : TReal1DArray;
begin
    SetLength(X, 1);
    SetLength(Y, 1);
    SetLength(W, 1);
    X[0] := 0;
    Y[0] := 0;
    W[0] := 1;
    BarycentricBuildXYW(X, Y, W, 1, B);
end;


(*************************************************************************
Tests whether constant C is solution of 1D LLS problem
*************************************************************************)
function Is1DSolution(N : AlglibInteger;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     C : Double):Boolean;
var
    I : AlglibInteger;
    V : Double;
    S1 : Double;
    S2 : Double;
    S3 : Double;
    Delta : Double;
begin
    Delta := 0.001;
    
    //
    // Test result
    //
    S1 := 0;
    I:=0;
    while I<=N-1 do
    begin
        S1 := S1+AP_Sqr(W[I]*(C-Y[I]));
        Inc(I);
    end;
    S2 := 0;
    S3 := 0;
    I:=0;
    while I<=N-1 do
    begin
        S2 := S2+AP_Sqr(W[I]*(C+Delta-Y[I]));
        S3 := S3+AP_Sqr(W[I]*(C-Delta-Y[I]));
        Inc(I);
    end;
    Result := AP_FP_Greater_Eq(S2,S1) and AP_FP_Greater_Eq(S3,S1);
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testratinterpolation_test_silent():Boolean;
begin
    Result := TestRI(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testratinterpolation_test():Boolean;
begin
    Result := TestRI(False);
end;


end.