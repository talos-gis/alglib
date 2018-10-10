unit testpolintunit;
interface
uses Math, Sysutils, Ap, tsort, ratinterpolation, blas, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, hblas, sblas, ortfac, rotations, bdsvd, svd, xblas, densesolver, linmin, minlbfgs, minlm, lsfit, ratint, taskgen, polint;

function TestPolInt(Silent : Boolean):Boolean;
function testpolintunit_test_silent():Boolean;
function testpolintunit_test():Boolean;

implementation

function InternalPolInt(const x : TReal1DArray;
     F : TReal1DArray;
     n : AlglibInteger;
     t : Double):Double;forward;
procedure BRCUnset(var B : BarycentricInterpolant);forward;


(*************************************************************************
Unit test
*************************************************************************)
function TestPolInt(Silent : Boolean):Boolean;
var
    WasErrors : Boolean;
    IntErrors : Boolean;
    FitErrors : Boolean;
    Threshold : Double;
    X : TReal1DArray;
    Y : TReal1DArray;
    W : TReal1DArray;
    X2 : TReal1DArray;
    Y2 : TReal1DArray;
    W2 : TReal1DArray;
    XFull : TReal1DArray;
    YFull : TReal1DArray;
    A : Double;
    B : Double;
    T : Double;
    I : AlglibInteger;
    K : AlglibInteger;
    XC : TReal1DArray;
    YC : TReal1DArray;
    DC : TInteger1DArray;
    Info : AlglibInteger;
    Info2 : AlglibInteger;
    V : Double;
    V0 : Double;
    V1 : Double;
    V2 : Double;
    S : Double;
    Xmin : Double;
    XMax : Double;
    RefRMS : Double;
    RefAvg : Double;
    RefAvgRel : Double;
    RefMax : Double;
    P : BarycentricInterpolant;
    P1 : BarycentricInterpolant;
    P2 : BarycentricInterpolant;
    Rep : PolynomialFitReport;
    Rep2 : PolynomialFitReport;
    N : AlglibInteger;
    M : AlglibInteger;
    MaxN : AlglibInteger;
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
begin
    WasErrors := False;
    IntErrors := False;
    FitErrors := False;
    MaxN := 5;
    PassCount := 20;
    Threshold := 1.0E8*MachineEpsilon;
    
    //
    // Test equidistant interpolation
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        N:=1;
        while N<=MaxN do
        begin
            
            //
            // prepare task:
            // * equidistant points
            // * random Y
            // * T in [A,B] or near (within 10% of its width)
            //
            repeat
                A := 2*RandomReal-1;
                B := 2*RandomReal-1;
            until AP_FP_Greater(AbsReal(A-B),0.2);
            T := A+(1.2*RandomReal-0.1)*(B-A);
            TaskGenInt1DEquidist(A, B, N, X, Y);
            
            //
            // test "fast" equidistant interpolation (no barycentric model)
            //
            IntErrors := IntErrors or AP_FP_Greater(AbsReal(PolynomialCalcEqDist(A, B, Y, N, T)-InternalPolInt(X, Y, N, T)),Threshold);
            
            //
            // test "slow" equidistant interpolation (create barycentric model)
            //
            BRCUnset(P);
            PolynomialBuild(X, Y, N, P);
            IntErrors := IntErrors or AP_FP_Greater(AbsReal(BarycentricCalc(P, T)-InternalPolInt(X, Y, N, T)),Threshold);
            
            //
            // test "fast" interpolation (create "fast" barycentric model)
            //
            BRCUnset(P);
            PolynomialBuildEqdist(A, B, Y, N, P);
            IntErrors := IntErrors or AP_FP_Greater(AbsReal(BarycentricCalc(P, T)-InternalPolInt(X, Y, N, T)),Threshold);
            Inc(N);
        end;
        Inc(Pass);
    end;
    
    //
    // Test Chebyshev-1 interpolation
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        N:=1;
        while N<=MaxN do
        begin
            
            //
            // prepare task:
            // * equidistant points
            // * random Y
            // * T in [A,B] or near (within 10% of its width)
            //
            repeat
                A := 2*RandomReal-1;
                B := 2*RandomReal-1;
            until AP_FP_Greater(AbsReal(A-B),0.2);
            T := A+(1.2*RandomReal-0.1)*(B-A);
            TaskGenInt1DCheb1(A, B, N, X, Y);
            
            //
            // test "fast" interpolation (no barycentric model)
            //
            IntErrors := IntErrors or AP_FP_Greater(AbsReal(PolynomialCalcCheb1(A, B, Y, N, T)-InternalPolInt(X, Y, N, T)),Threshold);
            
            //
            // test "slow" interpolation (create barycentric model)
            //
            BRCUnset(P);
            PolynomialBuild(X, Y, N, P);
            IntErrors := IntErrors or AP_FP_Greater(AbsReal(BarycentricCalc(P, T)-InternalPolInt(X, Y, N, T)),Threshold);
            
            //
            // test "fast" interpolation (create "fast" barycentric model)
            //
            BRCUnset(P);
            PolynomialBuildCheb1(A, B, Y, N, P);
            IntErrors := IntErrors or AP_FP_Greater(AbsReal(BarycentricCalc(P, T)-InternalPolInt(X, Y, N, T)),Threshold);
            Inc(N);
        end;
        Inc(Pass);
    end;
    
    //
    // Test Chebyshev-2 interpolation
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        N:=1;
        while N<=MaxN do
        begin
            
            //
            // prepare task:
            // * equidistant points
            // * random Y
            // * T in [A,B] or near (within 10% of its width)
            //
            repeat
                A := 2*RandomReal-1;
                B := 2*RandomReal-1;
            until AP_FP_Greater(AbsReal(A-B),0.2);
            T := A+(1.2*RandomReal-0.1)*(B-A);
            TaskGenInt1DCheb2(A, B, N, X, Y);
            
            //
            // test "fast" interpolation (no barycentric model)
            //
            IntErrors := IntErrors or AP_FP_Greater(AbsReal(PolynomialCalcCheb2(A, B, Y, N, T)-InternalPolInt(X, Y, N, T)),Threshold);
            
            //
            // test "slow" interpolation (create barycentric model)
            //
            BRCUnset(P);
            PolynomialBuild(X, Y, N, P);
            IntErrors := IntErrors or AP_FP_Greater(AbsReal(BarycentricCalc(P, T)-InternalPolInt(X, Y, N, T)),Threshold);
            
            //
            // test "fast" interpolation (create "fast" barycentric model)
            //
            BRCUnset(P);
            PolynomialBuildCheb2(A, B, Y, N, P);
            IntErrors := IntErrors or AP_FP_Greater(AbsReal(BarycentricCalc(P, T)-InternalPolInt(X, Y, N, T)),Threshold);
            Inc(N);
        end;
        Inc(Pass);
    end;
    
    //
    // crash-test: ability to solve tasks which will overflow/underflow
    // weights with straightforward implementation
    //
    N:=1;
    while N<=20 do
    begin
        A := -0.1*MaxRealNumber;
        B := +0.1*MaxRealNumber;
        TaskGenInt1DEquidist(A, B, N, X, Y);
        PolynomialBuild(X, Y, N, P);
        I:=0;
        while I<=N-1 do
        begin
            IntErrors := IntErrors or AP_FP_Eq(P.W[I],0);
            Inc(I);
        end;
        Inc(N);
    end;
    
    //
    // Test rational fitting:
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        N:=1;
        while N<=MaxN do
        begin
            
            //
            // N=M+K fitting (i.e. interpolation)
            //
            K:=0;
            while K<=N-1 do
            begin
                TaskGenInt1D(-1, 1, N, XFull, YFull);
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
                    X[I] := XFull[I];
                    Y[I] := YFull[I];
                    W[I] := 1+RandomReal;
                    Inc(I);
                end;
                I:=0;
                while I<=K-1 do
                begin
                    XC[I] := XFull[N-K+I];
                    YC[I] := YFull[N-K+I];
                    DC[I] := 0;
                    Inc(I);
                end;
                PolynomialFitWC(X, Y, W, N-K, XC, YC, DC, K, N, Info, P1, Rep);
                if Info<=0 then
                begin
                    FitErrors := True;
                end
                else
                begin
                    I:=0;
                    while I<=N-K-1 do
                    begin
                        FitErrors := FitErrors or AP_FP_Greater(AbsReal(BarycentricCalc(P1, X[I])-Y[I]),Threshold);
                        Inc(I);
                    end;
                    I:=0;
                    while I<=K-1 do
                    begin
                        FitErrors := FitErrors or AP_FP_Greater(AbsReal(BarycentricCalc(P1, XC[I])-YC[I]),Threshold);
                        Inc(I);
                    end;
                end;
                Inc(K);
            end;
            
            //
            // Testing constraints on derivatives.
            // Special tasks which will always have solution:
            // 1. P(0)=YC[0]
            // 2. P(0)=YC[0], P'(0)=YC[1]
            //
            if N>1 then
            begin
                M:=3;
                while M<=5 do
                begin
                    K:=1;
                    while K<=2 do
                    begin
                        TaskGenInt1D(-1, 1, N, X, Y);
                        SetLength(W, N);
                        SetLength(XC, 2);
                        SetLength(YC, 2);
                        SetLength(DC, 2);
                        I:=0;
                        while I<=N-1 do
                        begin
                            W[I] := 1+RandomReal;
                            Inc(I);
                        end;
                        XC[0] := 0;
                        YC[0] := 2*RandomReal-1;
                        DC[0] := 0;
                        XC[1] := 0;
                        YC[1] := 2*RandomReal-1;
                        DC[1] := 1;
                        PolynomialFitWC(X, Y, W, N, XC, YC, DC, K, M, Info, P1, Rep);
                        if Info<=0 then
                        begin
                            FitErrors := True;
                        end
                        else
                        begin
                            BarycentricDiff1(P1, 0.0, V0, V1);
                            FitErrors := FitErrors or AP_FP_Greater(AbsReal(V0-YC[0]),Threshold);
                            if K=2 then
                            begin
                                FitErrors := FitErrors or AP_FP_Greater(AbsReal(V1-YC[1]),Threshold);
                            end;
                        end;
                        Inc(K);
                    end;
                    Inc(M);
                end;
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
            XMin := 0;
            XMax := 2*Pi;
            I:=0;
            while I<=N-1 do
            begin
                X2[I] := 2*Pi*RandomReal;
                Y2[I] := Sin(X2[I]);
                W2[I] := 1;
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
            PolynomialBuild(X, Y, M, P1);
            PolynomialFitWC(X2, Y2, W2, N, XC, YC, DC, 0, M, Info, P2, Rep);
            if Info<=0 then
            begin
                FitErrors := True;
            end
            else
            begin
                
                //
                // calculate P1 (interpolant) RMS error, compare with P2 error
                //
                V1 := 0;
                V2 := 0;
                I:=0;
                while I<=N-1 do
                begin
                    V1 := V1+AP_Sqr(BarycentricCalc(P1, X2[I])-Y2[I]);
                    V2 := V2+AP_Sqr(BarycentricCalc(P2, X2[I])-Y2[I]);
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
            PolynomialFitWC(X, Y, W, N, XC, YC, DC, 0, M, Info, P1, Rep);
            PolynomialFit(X, Y, N, M, Info2, P2, Rep2);
            if (Info<=0) or (Info2<=0) then
            begin
                FitErrors := True;
            end
            else
            begin
                
                //
                // calculate P1 (interpolant), compare with P2 error
                // compare RMS errors
                //
                T := 2*RandomReal-1;
                V1 := BarycentricCalc(P1, T);
                V2 := BarycentricCalc(P2, T);
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
        PolynomialFit(X, Y, 4, 1, Info, P, Rep);
        if Info<=0 then
        begin
            FitErrors := True;
        end
        else
        begin
            S := BarycentricCalc(P, 0);
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
    WasErrors := IntErrors or FitErrors;
    if  not Silent then
    begin
        Write(Format('TESTING POLYNOMIAL INTERPOLATION AND FITTING'#13#10'',[]));
        
        //
        // Normal tests
        //
        Write(Format('INTERPOLATION TEST:                      ',[]));
        if IntErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('FITTING TEST:                            ',[]));
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


function InternalPolInt(const x : TReal1DArray;
     F : TReal1DArray;
     n : AlglibInteger;
     t : Double):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    F := DynamicArrayCopy(F);
    N := N-1;
    j:=0;
    while j<=n-1 do
    begin
        i:=j+1;
        while i<=n do
        begin
            F[i] := ((t-x[j])*F[i]-(t-x[i])*F[j])/(x[i]-x[j]);
            Inc(i);
        end;
        Inc(j);
    end;
    Result := F[n];
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
Silent unit test
*************************************************************************)
function testpolintunit_test_silent():Boolean;
begin
    Result := TestPolInt(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testpolintunit_test():Boolean;
begin
    Result := TestPolInt(False);
end;


end.