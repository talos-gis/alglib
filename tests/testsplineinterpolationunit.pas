unit testsplineinterpolationunit;
interface
uses Math, Sysutils, Ap, spline3, blas, trinverse, cholesky, spdsolve, lbfgs, minlm, reflections, bidiagonal, qr, lq, rotations, bdsvd, svd, lu, trlinsolve, rcond, leastsquares, lsfit, spline1d;

function TestSplineInterpolation(Silent : Boolean):Boolean;
function testsplineinterpolationunit_test_silent():Boolean;
function testsplineinterpolationunit_test():Boolean;

implementation

procedure LConst(A : Double;
     B : Double;
     const C : Spline1DInterpolant;
     LStep : Double;
     var L0 : Double;
     var L1 : Double;
     var L2 : Double);forward;
function TestUnpack(const C : Spline1DInterpolant;
     const X : TReal1DArray):Boolean;forward;
procedure UnsetSpline1D(var C : Spline1DInterpolant);forward;
procedure Unset1D(var X : TReal1DArray);forward;
function Is1DSolution(N : AlglibInteger;
     const Y : TReal1DArray;
     const W : TReal1DArray;
     C : Double):Boolean;forward;


function TestSplineInterpolation(Silent : Boolean):Boolean;
var
    WasErrors : Boolean;
    CSErrors : Boolean;
    HSErrors : Boolean;
    ASErrors : Boolean;
    LSErrors : Boolean;
    DSErrors : Boolean;
    UPErrors : Boolean;
    CPErrors : Boolean;
    LTErrors : Boolean;
    IErrors : Boolean;
    FitErrors : Boolean;
    NonStrictThreshold : Double;
    Threshold : Double;
    PassCount : AlglibInteger;
    LStep : Double;
    H : Double;
    MaxN : AlglibInteger;
    N : AlglibInteger;
    M : AlglibInteger;
    I : AlglibInteger;
    K : AlglibInteger;
    Pass : AlglibInteger;
    BLType : AlglibInteger;
    BRType : AlglibInteger;
    SType : AlglibInteger;
    X : TReal1DArray;
    Y : TReal1DArray;
    W : TReal1DArray;
    W2 : TReal1DArray;
    Y2 : TReal1DArray;
    D : TReal1DArray;
    XC : TReal1DArray;
    YC : TReal1DArray;
    DC : TInteger1DArray;
    C : Spline1DInterpolant;
    C2 : Spline1DInterpolant;
    Info : AlglibInteger;
    Info1 : AlglibInteger;
    Info2 : AlglibInteger;
    A : Double;
    B : Double;
    BL : Double;
    BR : Double;
    T : Double;
    SA : Double;
    SB : Double;
    V : Double;
    V1 : Double;
    V2 : Double;
    L10 : Double;
    L11 : Double;
    L12 : Double;
    L20 : Double;
    L21 : Double;
    L22 : Double;
    P0 : Double;
    P1 : Double;
    P2 : Double;
    S : Double;
    DS : Double;
    D2S : Double;
    S2 : Double;
    DS2 : Double;
    D2S2 : Double;
    VL : Double;
    VM : Double;
    VR : Double;
    Err : Double;
    Rep : Spline1DFitReport;
    Rep2 : Spline1DFitReport;
    RefRMS : Double;
    RefAvg : Double;
    RefAvgRel : Double;
    RefMax : Double;
    RA : TReal1DArray;
    RA2 : TReal1DArray;
    RALen : AlglibInteger;
begin
    WasErrors := False;
    PassCount := 20;
    LStep := 0.005;
    H := 0.00001;
    MaxN := 10;
    Threshold := 10000*MachineEpsilon;
    NonStrictThreshold := 0.00001;
    LSErrors := False;
    CSErrors := False;
    HSErrors := False;
    ASErrors := False;
    DSErrors := False;
    CPErrors := False;
    UPErrors := False;
    LTErrors := False;
    IErrors := False;
    FitErrors := False;
    
    //
    // General test: linear, cubic, Hermite, Akima
    //
    N:=2;
    while N<=MaxN do
    begin
        SetLength(X, N-1+1);
        SetLength(Y, N-1+1);
        SetLength(D, N-1+1);
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // Prepare task
            //
            A := -1-RandomReal;
            B := +1+RandomReal;
            BL := 2*RandomReal-1;
            BR := 2*RandomReal-1;
            I:=0;
            while I<=N-1 do
            begin
                X[I] := 0.5*(B+A)+0.5*(B-A)*Cos(PI*(2*i+1)/(2*n));
                if I=0 then
                begin
                    X[I] := A;
                end;
                if I=N-1 then
                begin
                    X[I] := B;
                end;
                Y[I] := Cos(1.3*Pi*X[I]+0.4);
                D[I] := -1.3*Pi*Sin(1.3*Pi*X[I]+0.4);
                Inc(I);
            end;
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
                    T := D[I];
                    D[I] := D[K];
                    D[K] := T;
                end;
                Inc(I);
            end;
            
            //
            // Build linear spline
            // Test for general interpolation scheme properties:
            // * values at nodes
            // * continuous function
            // Test for specific properties is implemented below.
            //
            Spline1DBuildLinear(X, Y, N, C);
            Err := 0;
            I:=0;
            while I<=N-1 do
            begin
                Err := Max(Err, AbsReal(Y[I]-Spline1DCalc(C, X[I])));
                Inc(I);
            end;
            LSErrors := LSErrors or AP_FP_Greater(Err,Threshold);
            LConst(A, B, C, LStep, L10, L11, L12);
            LConst(A, B, C, LStep/3, L20, L21, L22);
            LSErrors := LSErrors or AP_FP_Greater(L20/L10,1.2);
            
            //
            // Build cubic spline.
            // Test for interpolation scheme properties:
            // * values at nodes
            // * boundary conditions
            // * continuous function
            // * continuous first derivative
            // * continuous second derivative
            //
            BLType:=0;
            while BLType<=2 do
            begin
                BRType:=0;
                while BRType<=2 do
                begin
                    Spline1DBuildCubic(X, Y, N, BLType, BL, BRType, BR, C);
                    Err := 0;
                    I:=0;
                    while I<=N-1 do
                    begin
                        Err := Max(Err, AbsReal(Y[I]-Spline1DCalc(C, X[I])));
                        Inc(I);
                    end;
                    CSErrors := CSErrors or AP_FP_Greater(Err,Threshold);
                    Err := 0;
                    if BLType=0 then
                    begin
                        Spline1DDiff(C, A-H, S, DS, D2S);
                        Spline1DDiff(C, A+H, S2, DS2, D2S2);
                        T := (D2S2-D2S)/(2*H);
                        Err := Max(Err, AbsReal(T));
                    end;
                    if BLType=1 then
                    begin
                        T := (Spline1DCalc(C, A+H)-Spline1DCalc(C, A-H))/(2*H);
                        Err := Max(Err, AbsReal(BL-T));
                    end;
                    if BLType=2 then
                    begin
                        T := (Spline1DCalc(C, A+H)-2*Spline1DCalc(C, A)+Spline1DCalc(C, A-H))/AP_Sqr(H);
                        Err := Max(Err, AbsReal(BL-T));
                    end;
                    if BRType=0 then
                    begin
                        Spline1DDiff(C, B-H, S, DS, D2S);
                        Spline1DDiff(C, B+H, S2, DS2, D2S2);
                        T := (D2S2-D2S)/(2*H);
                        Err := Max(Err, AbsReal(T));
                    end;
                    if BRType=1 then
                    begin
                        T := (Spline1DCalc(C, B+H)-Spline1DCalc(C, B-H))/(2*H);
                        Err := Max(Err, AbsReal(BR-T));
                    end;
                    if BRType=2 then
                    begin
                        T := (Spline1DCalc(C, B+H)-2*Spline1DCalc(C, B)+Spline1DCalc(C, B-H))/AP_Sqr(H);
                        Err := Max(Err, AbsReal(BR-T));
                    end;
                    CSErrors := CSErrors or AP_FP_Greater(Err,1.0E-3);
                    LConst(A, B, C, LStep, L10, L11, L12);
                    LConst(A, B, C, LStep/3, L20, L21, L22);
                    CSErrors := CSErrors or AP_FP_Greater(L20/L10,1.2) and AP_FP_Greater(L10,1.0E-6);
                    CSErrors := CSErrors or AP_FP_Greater(L21/L11,1.2) and AP_FP_Greater(L11,1.0E-6);
                    CSErrors := CSErrors or AP_FP_Greater(L22/L12,1.2) and AP_FP_Greater(L12,1.0E-6);
                    Inc(BRType);
                end;
                Inc(BLType);
            end;
            
            //
            // Build Hermite spline.
            // Test for interpolation scheme properties:
            // * values and derivatives at nodes
            // * continuous function
            // * continuous first derivative
            //
            Spline1DBuildHermite(X, Y, D, N, C);
            Err := 0;
            I:=0;
            while I<=N-1 do
            begin
                Err := Max(Err, AbsReal(Y[I]-Spline1DCalc(C, X[I])));
                Inc(I);
            end;
            HSErrors := HSErrors or AP_FP_Greater(Err,Threshold);
            Err := 0;
            I:=0;
            while I<=N-1 do
            begin
                T := (Spline1DCalc(C, X[I]+H)-Spline1DCalc(C, X[I]-H))/(2*H);
                Err := Max(Err, AbsReal(D[I]-T));
                Inc(I);
            end;
            HSErrors := HSErrors or AP_FP_Greater(Err,1.0E-3);
            LConst(A, B, C, LStep, L10, L11, L12);
            LConst(A, B, C, LStep/3, L20, L21, L22);
            HSErrors := HSErrors or AP_FP_Greater(L20/L10,1.2);
            HSErrors := HSErrors or AP_FP_Greater(L21/L11,1.2);
            
            //
            // Build Akima spline
            // Test for general interpolation scheme properties:
            // * values at nodes
            // * continuous function
            // * continuous first derivative
            // Test for specific properties is implemented below.
            //
            if N>=5 then
            begin
                Spline1DBuildAkima(X, Y, N, C);
                Err := 0;
                I:=0;
                while I<=N-1 do
                begin
                    Err := Max(Err, AbsReal(Y[I]-Spline1DCalc(C, X[I])));
                    Inc(I);
                end;
                ASErrors := ASErrors or AP_FP_Greater(Err,Threshold);
                LConst(A, B, C, LStep, L10, L11, L12);
                LConst(A, B, C, LStep/3, L20, L21, L22);
                HSErrors := HSErrors or AP_FP_Greater(L20/L10,1.2);
                HSErrors := HSErrors or AP_FP_Greater(L21/L11,1.2);
            end;
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // Special linear spline test:
    // test for linearity between x[i] and x[i+1]
    //
    N:=2;
    while N<=MaxN do
    begin
        SetLength(X, N-1+1);
        SetLength(Y, N-1+1);
        
        //
        // Prepare task
        //
        A := -1;
        B := +1;
        I:=0;
        while I<=N-1 do
        begin
            X[I] := A+(B-A)*I/(N-1);
            Y[I] := 2*RandomReal-1;
            Inc(I);
        end;
        Spline1DBuildLinear(X, Y, N, C);
        
        //
        // Test
        //
        Err := 0;
        K:=0;
        while K<=N-2 do
        begin
            A := X[K];
            B := X[K+1];
            Pass:=1;
            while Pass<=PassCount do
            begin
                T := A+(B-A)*RandomReal;
                V := Y[K]+(T-A)/(B-A)*(Y[K+1]-Y[K]);
                Err := Max(Err, AbsReal(Spline1DCalc(C, T)-V));
                Inc(Pass);
            end;
            Inc(K);
        end;
        LSErrors := LSErrors or AP_FP_Greater(Err,Threshold);
        Inc(N);
    end;
    
    //
    // Special Akima test: test outlier sensitivity
    // Spline value at (x[i], x[i+1]) should depend from
    // f[i-2], f[i-1], f[i], f[i+1], f[i+2], f[i+3] only.
    //
    N:=5;
    while N<=MaxN do
    begin
        SetLength(X, N-1+1);
        SetLength(Y, N-1+1);
        SetLength(Y2, N-1+1);
        
        //
        // Prepare unperturbed Akima spline
        //
        A := -1;
        B := +1;
        I:=0;
        while I<=N-1 do
        begin
            X[I] := A+(B-A)*I/(N-1);
            Y[I] := Cos(1.3*Pi*X[I]+0.4);
            Inc(I);
        end;
        Spline1DBuildAkima(X, Y, N, C);
        
        //
        // Process perturbed tasks
        //
        Err := 0;
        K:=0;
        while K<=N-1 do
        begin
            APVMove(@Y2[0], 0, N-1, @Y[0], 0, N-1);
            Y2[K] := 5;
            Spline1DBuildAkima(X, Y2, N, C2);
            
            //
            // Test left part independence
            //
            if K-3>=1 then
            begin
                A := -1;
                B := X[K-3];
                Pass:=1;
                while Pass<=PassCount do
                begin
                    T := A+(B-A)*RandomReal;
                    Err := Max(Err, AbsReal(Spline1DCalc(C, T)-Spline1DCalc(C2, T)));
                    Inc(Pass);
                end;
            end;
            
            //
            // Test right part independence
            //
            if K+3<=N-2 then
            begin
                A := X[K+3];
                B := +1;
                Pass:=1;
                while Pass<=PassCount do
                begin
                    T := A+(B-A)*RandomReal;
                    Err := Max(Err, AbsReal(Spline1DCalc(C, T)-Spline1DCalc(C2, T)));
                    Inc(Pass);
                end;
            end;
            Inc(K);
        end;
        ASErrors := ASErrors or AP_FP_Greater(Err,Threshold);
        Inc(N);
    end;
    
    //
    // Differentiation, copy/serialize/unpack test
    //
    N:=2;
    while N<=MaxN do
    begin
        SetLength(X, N-1+1);
        SetLength(Y, N-1+1);
        
        //
        // Prepare cubic spline
        //
        A := -1-RandomReal;
        B := +1+RandomReal;
        I:=0;
        while I<=N-1 do
        begin
            X[I] := A+(B-A)*I/(N-1);
            Y[I] := Cos(1.3*Pi*X[I]+0.4);
            Inc(I);
        end;
        Spline1DBuildCubic(X, Y, N, 2, 0.0, 2, 0.0, C);
        
        //
        // Test diff
        //
        Err := 0;
        Pass:=1;
        while Pass<=PassCount do
        begin
            T := A+(B-A)*RandomReal;
            Spline1DDiff(C, T, S, DS, D2S);
            VL := Spline1DCalc(C, T-H);
            VM := Spline1DCalc(C, T);
            VR := Spline1DCalc(C, T+H);
            Err := Max(Err, AbsReal(S-VM));
            Err := Max(Err, AbsReal(DS-(VR-VL)/(2*H)));
            Err := Max(Err, AbsReal(D2S-(VR-2*VM+VL)/AP_Sqr(H)));
            Inc(Pass);
        end;
        DSErrors := DSErrors or AP_FP_Greater(Err,0.001);
        
        //
        // Test copy
        //
        UnsetSpline1D(C2);
        Spline1DCopy(C, C2);
        Err := 0;
        Pass:=1;
        while Pass<=PassCount do
        begin
            T := A+(B-A)*RandomReal;
            Err := Max(Err, AbsReal(Spline1DCalc(C, T)-Spline1DCalc(C2, T)));
            Inc(Pass);
        end;
        CPErrors := CPErrors or AP_FP_Greater(Err,Threshold);
        
        //
        // Test serialize/deserialize
        //
        UnsetSpline1D(C2);
        Unset1D(RA2);
        Spline1DSerialize(C, RA, RALen);
        SetLength(RA2, RALen);
        APVMove(@RA2[0], 0, RALen-1, @RA[0], 0, RALen-1);
        Spline1DUnserialize(RA2, C2);
        Err := 0;
        Pass:=1;
        while Pass<=PassCount do
        begin
            T := A+(B-A)*RandomReal;
            Err := Max(Err, AbsReal(Spline1DCalc(C, T)-Spline1DCalc(C2, T)));
            Inc(Pass);
        end;
        CPErrors := CPErrors or AP_FP_Greater(Err,Threshold);
        
        //
        // Test unpack
        //
        UPErrors := UPErrors or  not TestUnpack(C, X);
        
        //
        // Test lin.trans.
        //
        Err := 0;
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // LinTransX, general A
            //
            SA := 4*RandomReal-2;
            SB := 2*RandomReal-1;
            T := A+(B-A)*RandomReal;
            Spline1DCopy(C, C2);
            Spline1DLinTransX(C2, SA, SB);
            Err := Max(Err, AbsReal(Spline1DCalc(C, T)-Spline1DCalc(C2, (T-SB)/SA)));
            
            //
            // LinTransX, special case: A=0
            //
            SB := 2*RandomReal-1;
            T := A+(B-A)*RandomReal;
            Spline1DCopy(C, C2);
            Spline1DLinTransX(C2, 0, SB);
            Err := Max(Err, AbsReal(Spline1DCalc(C, SB)-Spline1DCalc(C2, T)));
            
            //
            // LinTransY
            //
            SA := 2*RandomReal-1;
            SB := 2*RandomReal-1;
            T := A+(B-A)*RandomReal;
            Spline1DCopy(C, C2);
            Spline1DLinTransY(C2, SA, SB);
            Err := Max(Err, AbsReal(SA*Spline1DCalc(C, T)+SB-Spline1DCalc(C2, T)));
            Inc(Pass);
        end;
        LTErrors := LTErrors or AP_FP_Greater(Err,Threshold);
        Inc(N);
    end;
    
    //
    // Testing integration. Two tests are performed:
    // * approximate test (well behaved smooth function, many points,
    //   integration inside [a,b])
    // * exact test (integration of parabola, outside of [a,b]
    //
    Err := 0;
    N:=20;
    while N<=35 do
    begin
        SetLength(X, N-1+1);
        SetLength(Y, N-1+1);
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // Prepare cubic spline
            //
            A := -1-0.2*RandomReal;
            B := +1+0.2*RandomReal;
            I:=0;
            while I<=N-1 do
            begin
                X[I] := A+(B-A)*I/(N-1);
                Y[I] := Sin(Pi*X[I]+0.4)+Exp(X[I]);
                Inc(I);
            end;
            BL := Pi*Cos(Pi*A+0.4)+Exp(A);
            BR := Pi*Cos(Pi*B+0.4)+Exp(B);
            Spline1DBuildCubic(X, Y, N, 1, BL, 1, BR, C);
            
            //
            // Test
            //
            T := A+(B-A)*RandomReal;
            V := -Cos(Pi*A+0.4)/Pi+Exp(A);
            V := -Cos(Pi*T+0.4)/Pi+Exp(T)-V;
            V := V-Spline1DIntegrate(C, T);
            Err := Max(Err, AbsReal(V));
            Inc(Pass);
        end;
        Inc(N);
    end;
    IErrors := IErrors or AP_FP_Greater(Err,0.001);
    P0 := 2*RandomReal-1;
    P1 := 2*RandomReal-1;
    P2 := 2*RandomReal-1;
    A := -RandomReal-0.5;
    B := +RandomReal+0.5;
    N := 2;
    SetLength(X, N);
    SetLength(Y, N);
    SetLength(D, N);
    X[0] := A;
    Y[0] := P0+P1*A+P2*AP_Sqr(A);
    D[0] := P1+2*P2*A;
    X[1] := B;
    Y[1] := P0+P1*B+P2*AP_Sqr(B);
    D[1] := P1+2*P2*B;
    Spline1DBuildHermite(X, Y, D, N, C);
    BL := Min(A, B)-AbsReal(B-A);
    BR := Min(A, B)+AbsReal(B-A);
    Err := 0;
    Pass:=1;
    while Pass<=100 do
    begin
        T := BL+(BR-BL)*RandomReal;
        V := P0*T+P1*AP_Sqr(T)/2+P2*AP_Sqr(T)*T/3-(P0*A+P1*AP_Sqr(A)/2+P2*AP_Sqr(A)*A/3);
        V := V-Spline1DIntegrate(C, T);
        Err := Max(Err, AbsReal(V));
        Inc(Pass);
    end;
    IErrors := IErrors or AP_FP_Greater(Err,Threshold);
    
    //
    // Test fitting.
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        
        //
        // Cubic splines
        // Ability to handle boundary constraints (1-4 constraints on F, dF/dx).
        //
        M:=4;
        while M<=8 do
        begin
            K:=1;
            while K<=4 do
            begin
                if K>=M then
                begin
                    Inc(K);
                    Continue;
                end;
                N := 100;
                SetLength(X, N);
                SetLength(Y, N);
                SetLength(W, N);
                SetLength(XC, 4);
                SetLength(YC, 4);
                SetLength(DC, 4);
                SA := 1+RandomReal;
                SB := 2*RandomReal-1;
                I:=0;
                while I<=N-1 do
                begin
                    X[I] := SA*RandomReal+SB;
                    Y[I] := 2*RandomReal-1;
                    W[I] := 1+RandomReal;
                    Inc(I);
                end;
                XC[0] := SB;
                YC[0] := 2*RandomReal-1;
                DC[0] := 0;
                XC[1] := SB;
                YC[1] := 2*RandomReal-1;
                DC[1] := 1;
                XC[2] := SA+SB;
                YC[2] := 2*RandomReal-1;
                DC[2] := 0;
                XC[3] := SA+SB;
                YC[3] := 2*RandomReal-1;
                DC[3] := 1;
                Spline1DFitCubicWC(X, Y, W, N, XC, YC, DC, K, M, Info, C, Rep);
                if Info<=0 then
                begin
                    FitErrors := True;
                end
                else
                begin
                    
                    //
                    // Check that constraints are satisfied
                    //
                    I:=0;
                    while I<=K-1 do
                    begin
                        Spline1DDiff(C, XC[I], S, DS, D2S);
                        if DC[I]=0 then
                        begin
                            FitErrors := FitErrors or AP_FP_Greater(AbsReal(S-YC[I]),Threshold);
                        end;
                        if DC[I]=1 then
                        begin
                            FitErrors := FitErrors or AP_FP_Greater(AbsReal(DS-YC[I]),Threshold);
                        end;
                        if DC[I]=2 then
                        begin
                            FitErrors := FitErrors or AP_FP_Greater(AbsReal(D2S-YC[I]),Threshold);
                        end;
                        Inc(I);
                    end;
                end;
                Inc(K);
            end;
            Inc(M);
        end;
        
        //
        // Cubic splines
        // Ability to handle one internal constraint
        //
        M:=4;
        while M<=8 do
        begin
            N := 100;
            SetLength(X, N);
            SetLength(Y, N);
            SetLength(W, N);
            SetLength(XC, 1);
            SetLength(YC, 1);
            SetLength(DC, 1);
            SA := 1+RandomReal;
            SB := 2*RandomReal-1;
            I:=0;
            while I<=N-1 do
            begin
                X[I] := SA*RandomReal+SB;
                Y[I] := 2*RandomReal-1;
                W[I] := 1+RandomReal;
                Inc(I);
            end;
            XC[0] := SA*RandomReal+SB;
            YC[0] := 2*RandomReal-1;
            DC[0] := RandomInteger(2);
            Spline1DFitCubicWC(X, Y, W, N, XC, YC, DC, 1, M, Info, C, Rep);
            if Info<=0 then
            begin
                FitErrors := True;
            end
            else
            begin
                
                //
                // Check that constraints are satisfied
                //
                Spline1DDiff(C, XC[0], S, DS, D2S);
                if DC[0]=0 then
                begin
                    FitErrors := FitErrors or AP_FP_Greater(AbsReal(S-YC[0]),Threshold);
                end;
                if DC[0]=1 then
                begin
                    FitErrors := FitErrors or AP_FP_Greater(AbsReal(DS-YC[0]),Threshold);
                end;
                if DC[0]=2 then
                begin
                    FitErrors := FitErrors or AP_FP_Greater(AbsReal(D2S-YC[0]),Threshold);
                end;
            end;
            Inc(M);
        end;
        
        //
        // Hermite splines
        // Ability to handle boundary constraints (1-4 constraints on F, dF/dx).
        //
        M:=4;
        while M<=8 do
        begin
            K:=1;
            while K<=4 do
            begin
                if K>=M then
                begin
                    Inc(K);
                    Continue;
                end;
                if M mod 2<>0 then
                begin
                    Inc(K);
                    Continue;
                end;
                N := 100;
                SetLength(X, N);
                SetLength(Y, N);
                SetLength(W, N);
                SetLength(XC, 4);
                SetLength(YC, 4);
                SetLength(DC, 4);
                SA := 1+RandomReal;
                SB := 2*RandomReal-1;
                I:=0;
                while I<=N-1 do
                begin
                    X[I] := SA*RandomReal+SB;
                    Y[I] := 2*RandomReal-1;
                    W[I] := 1+RandomReal;
                    Inc(I);
                end;
                XC[0] := SB;
                YC[0] := 2*RandomReal-1;
                DC[0] := 0;
                XC[1] := SB;
                YC[1] := 2*RandomReal-1;
                DC[1] := 1;
                XC[2] := SA+SB;
                YC[2] := 2*RandomReal-1;
                DC[2] := 0;
                XC[3] := SA+SB;
                YC[3] := 2*RandomReal-1;
                DC[3] := 1;
                Spline1DFitHermiteWC(X, Y, W, N, XC, YC, DC, K, M, Info, C, Rep);
                if Info<=0 then
                begin
                    FitErrors := True;
                end
                else
                begin
                    
                    //
                    // Check that constraints are satisfied
                    //
                    I:=0;
                    while I<=K-1 do
                    begin
                        Spline1DDiff(C, XC[I], S, DS, D2S);
                        if DC[I]=0 then
                        begin
                            FitErrors := FitErrors or AP_FP_Greater(AbsReal(S-YC[I]),Threshold);
                        end;
                        if DC[I]=1 then
                        begin
                            FitErrors := FitErrors or AP_FP_Greater(AbsReal(DS-YC[I]),Threshold);
                        end;
                        if DC[I]=2 then
                        begin
                            FitErrors := FitErrors or AP_FP_Greater(AbsReal(D2S-YC[I]),Threshold);
                        end;
                        Inc(I);
                    end;
                end;
                Inc(K);
            end;
            Inc(M);
        end;
        
        //
        // Hermite splines
        // Ability to handle one internal constraint
        //
        M:=4;
        while M<=8 do
        begin
            if M mod 2<>0 then
            begin
                Inc(M);
                Continue;
            end;
            N := 100;
            SetLength(X, N);
            SetLength(Y, N);
            SetLength(W, N);
            SetLength(XC, 1);
            SetLength(YC, 1);
            SetLength(DC, 1);
            SA := 1+RandomReal;
            SB := 2*RandomReal-1;
            I:=0;
            while I<=N-1 do
            begin
                X[I] := SA*RandomReal+SB;
                Y[I] := 2*RandomReal-1;
                W[I] := 1+RandomReal;
                Inc(I);
            end;
            XC[0] := SA*RandomReal+SB;
            YC[0] := 2*RandomReal-1;
            DC[0] := RandomInteger(2);
            Spline1DFitHermiteWC(X, Y, W, N, XC, YC, DC, 1, M, Info, C, Rep);
            if Info<=0 then
            begin
                FitErrors := True;
            end
            else
            begin
                
                //
                // Check that constraints are satisfied
                //
                Spline1DDiff(C, XC[0], S, DS, D2S);
                if DC[0]=0 then
                begin
                    FitErrors := FitErrors or AP_FP_Greater(AbsReal(S-YC[0]),Threshold);
                end;
                if DC[0]=1 then
                begin
                    FitErrors := FitErrors or AP_FP_Greater(AbsReal(DS-YC[0]),Threshold);
                end;
                if DC[0]=2 then
                begin
                    FitErrors := FitErrors or AP_FP_Greater(AbsReal(D2S-YC[0]),Threshold);
                end;
            end;
            Inc(M);
        end;
        Inc(Pass);
    end;
    M:=4;
    while M<=8 do
    begin
        SType:=0;
        while SType<=1 do
        begin
            Pass:=1;
            while Pass<=PassCount do
            begin
                if (SType=1) and (M mod 2<>0) then
                begin
                    Inc(Pass);
                    Continue;
                end;
                
                //
                // cubic/Hermite spline fitting:
                // * generate "template spline" C2
                // * generate 2*N points from C2, such that result of
                //   ideal fit should be equal to C2
                // * fit, store in C
                // * compare C and C2
                //
                SA := 1+RandomReal;
                SB := 2*RandomReal-1;
                if SType=0 then
                begin
                    SetLength(X, M-2);
                    SetLength(Y, M-2);
                    I:=0;
                    while I<=M-2-1 do
                    begin
                        X[I] := SA*I/(M-2-1)+SB;
                        Y[I] := 2*RandomReal-1;
                        Inc(I);
                    end;
                    Spline1DBuildCubic(X, Y, M-2, 1, 2*RandomReal-1, 1, 2*RandomReal-1, C2);
                end;
                if SType=1 then
                begin
                    SetLength(X, M div 2);
                    SetLength(Y, M div 2);
                    SetLength(D, M div 2);
                    I:=0;
                    while I<=M div 2-1 do
                    begin
                        X[I] := SA*I/(M div 2-1)+SB;
                        Y[I] := 2*RandomReal-1;
                        D[I] := 2*RandomReal-1;
                        Inc(I);
                    end;
                    Spline1DBuildHermite(X, Y, D, M div 2, C2);
                end;
                N := 50;
                SetLength(X, 2*N);
                SetLength(Y, 2*N);
                SetLength(W, 2*N);
                I:=0;
                while I<=N-1 do
                begin
                    
                    //
                    // "if i=0" and "if i=1" are needed to
                    // synchronize interval size for C2 and
                    // spline being fitted (i.e. C).
                    //
                    T := RandomReal;
                    X[I] := SA*RandomReal+SB;
                    if I=0 then
                    begin
                        X[I] := SB;
                    end;
                    if I=1 then
                    begin
                        X[I] := SA+SB;
                    end;
                    V := Spline1DCalc(C2, X[I]);
                    Y[I] := V+T;
                    W[I] := 1+RandomReal;
                    X[N+I] := X[I];
                    Y[N+I] := V-T;
                    W[N+I] := W[I];
                    Inc(I);
                end;
                if SType=0 then
                begin
                    Spline1DFitCubicWC(X, Y, W, 2*N, XC, YC, DC, 0, M, Info, C, Rep);
                end;
                if SType=1 then
                begin
                    Spline1DFitHermiteWC(X, Y, W, 2*N, XC, YC, DC, 0, M, Info, C, Rep);
                end;
                if Info<=0 then
                begin
                    FitErrors := True;
                end
                else
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        V := SA*RandomReal+SB;
                        FitErrors := FitErrors or AP_FP_Greater(AbsReal(Spline1DCalc(C, V)-Spline1DCalc(C2, V)),Threshold);
                        Inc(I);
                    end;
                end;
                Inc(Pass);
            end;
            Inc(SType);
        end;
        Inc(M);
    end;
    M:=4;
    while M<=8 do
    begin
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // prepare points/weights
            //
            SA := 1+RandomReal;
            SB := 2*RandomReal-1;
            N := 10+RandomInteger(10);
            SetLength(X, N);
            SetLength(Y, N);
            SetLength(W, N);
            I:=0;
            while I<=N-1 do
            begin
                X[I] := SA*RandomReal+SB;
                Y[I] := 2*RandomReal-1;
                W[I] := 1;
                Inc(I);
            end;
            
            //
            // Fit cubic with unity weights, without weights, then compare
            //
            if M>=4 then
            begin
                Spline1DFitCubicWC(X, Y, W, N, XC, YC, DC, 0, M, Info1, C, Rep);
                Spline1DFitCubic(X, Y, N, M, Info2, C2, Rep2);
                if (Info1<=0) or (Info2<=0) then
                begin
                    FitErrors := True;
                end
                else
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        V := SA*RandomReal+SB;
                        FitErrors := FitErrors or AP_FP_Neq(Spline1DCalc(C, V),Spline1DCalc(C2, V));
                        FitErrors := FitErrors or AP_FP_Neq(Rep.TaskRCond,Rep2.TaskRCond);
                        FitErrors := FitErrors or AP_FP_Neq(Rep.RMSError,Rep2.RMSError);
                        FitErrors := FitErrors or AP_FP_Neq(Rep.AvgError,Rep2.AvgError);
                        FitErrors := FitErrors or AP_FP_Neq(Rep.AvgRelError,Rep2.AvgRelError);
                        FitErrors := FitErrors or AP_FP_Neq(Rep.MaxError,Rep2.MaxError);
                        Inc(I);
                    end;
                end;
            end;
            
            //
            // Fit Hermite with unity weights, without weights, then compare
            //
            if (M>=4) and (M mod 2=0) then
            begin
                Spline1DFitHermiteWC(X, Y, W, N, XC, YC, DC, 0, M, Info1, C, Rep);
                Spline1DFitHermite(X, Y, N, M, Info2, C2, Rep2);
                if (Info1<=0) or (Info2<=0) then
                begin
                    FitErrors := True;
                end
                else
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        V := SA*RandomReal+SB;
                        FitErrors := FitErrors or AP_FP_Neq(Spline1DCalc(C, V),Spline1DCalc(C2, V));
                        FitErrors := FitErrors or AP_FP_Neq(Rep.TaskRCond,Rep2.TaskRCond);
                        FitErrors := FitErrors or AP_FP_Neq(Rep.RMSError,Rep2.RMSError);
                        FitErrors := FitErrors or AP_FP_Neq(Rep.AvgError,Rep2.AvgError);
                        FitErrors := FitErrors or AP_FP_Neq(Rep.AvgRelError,Rep2.AvgRelError);
                        FitErrors := FitErrors or AP_FP_Neq(Rep.MaxError,Rep2.MaxError);
                        Inc(I);
                    end;
                end;
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
        // Test cubic fitting
        //
        Spline1DFitCubic(X, Y, 4, 4, Info, C, Rep);
        if Info<=0 then
        begin
            FitErrors := True;
        end
        else
        begin
            S := Spline1DCalc(C, 0);
            FitErrors := FitErrors or AP_FP_Greater(AbsReal(S-V),Threshold);
            FitErrors := FitErrors or AP_FP_Greater(AbsReal(Rep.RMSError-RefRMS),Threshold);
            FitErrors := FitErrors or AP_FP_Greater(AbsReal(Rep.AvgError-RefAvg),Threshold);
            FitErrors := FitErrors or AP_FP_Greater(AbsReal(Rep.AvgRelError-RefAvgRel),Threshold);
            FitErrors := FitErrors or AP_FP_Greater(AbsReal(Rep.MaxError-RefMax),Threshold);
        end;
        
        //
        // Test cubic fitting
        //
        Spline1DFitHermite(X, Y, 4, 4, Info, C, Rep);
        if Info<=0 then
        begin
            FitErrors := True;
        end
        else
        begin
            S := Spline1DCalc(C, 0);
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
    WasErrors := LSErrors or CSErrors or HSErrors or ASErrors or DSErrors or CPErrors or UPErrors or LTErrors or IErrors or FitErrors;
    if  not Silent then
    begin
        Write(Format('TESTING SPLINE INTERPOLATION'#13#10'',[]));
        
        //
        // Normal tests
        //
        Write(Format('LINEAR SPLINE TEST:                      ',[]));
        if LSErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('CUBIC SPLINE TEST:                       ',[]));
        if CSErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('HERMITE SPLINE TEST:                     ',[]));
        if HSErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('AKIMA SPLINE TEST:                       ',[]));
        if ASErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('DIFFERENTIATION TEST:                    ',[]));
        if DSErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('COPY/SERIALIZATION TEST:                 ',[]));
        if CPErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('UNPACK TEST:                             ',[]));
        if UPErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('LIN.TRANS. TEST:                         ',[]));
        if LTErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('INTEGRATION TEST:                        ',[]));
        if IErrors then
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


(*************************************************************************
Lipschitz constants for spline inself, first and second derivatives.
*************************************************************************)
procedure LConst(A : Double;
     B : Double;
     const C : Spline1DInterpolant;
     LStep : Double;
     var L0 : Double;
     var L1 : Double;
     var L2 : Double);
var
    T : Double;
    VL : Double;
    VM : Double;
    VR : Double;
    PrevF : Double;
    PrevD : Double;
    PrevD2 : Double;
    F : Double;
    D : Double;
    D2 : Double;
begin
    L0 := 0;
    L1 := 0;
    L2 := 0;
    T := A-0.1;
    VL := Spline1DCalc(C, T-2*LStep);
    VM := Spline1DCalc(C, T-LStep);
    VR := Spline1DCalc(C, T);
    F := VM;
    D := (VR-VL)/(2*LStep);
    D2 := (VR-2*VM+VL)/AP_Sqr(LStep);
    while AP_FP_Less_Eq(T,B+0.1) do
    begin
        PrevF := F;
        PrevD := D;
        PrevD2 := D2;
        VL := VM;
        VM := VR;
        VR := Spline1DCalc(C, T+LStep);
        F := VM;
        D := (VR-VL)/(2*LStep);
        D2 := (VR-2*VM+VL)/AP_Sqr(LStep);
        L0 := Max(L0, AbsReal((F-PrevF)/LStep));
        L1 := Max(L1, AbsReal((D-PrevD)/LStep));
        L2 := Max(L2, AbsReal((D2-PrevD2)/LStep));
        T := T+LStep;
    end;
end;


(*************************************************************************
Unpack testing
*************************************************************************)
function TestUnpack(const C : Spline1DInterpolant;
     const X : TReal1DArray):Boolean;
var
    I : AlglibInteger;
    N : AlglibInteger;
    Err : Double;
    T : Double;
    V1 : Double;
    V2 : Double;
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    Tbl : TReal2DArray;
begin
    PassCount := 20;
    Err := 0;
    Spline1DUnpack(C, N, Tbl);
    I:=0;
    while I<=N-2 do
    begin
        Pass:=1;
        while Pass<=PassCount do
        begin
            T := RandomReal*(Tbl[I,1]-Tbl[I,0]);
            V1 := Tbl[I,2]+T*Tbl[I,3]+AP_Sqr(T)*Tbl[I,4]+T*AP_Sqr(T)*Tbl[I,5];
            V2 := Spline1DCalc(C, Tbl[I,0]+T);
            Err := Max(Err, AbsReal(V1-V2));
            Inc(Pass);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        Err := Max(Err, AbsReal(X[I]-Tbl[I,0]));
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        Err := Max(Err, AbsReal(X[I+1]-Tbl[I,1]));
        Inc(I);
    end;
    Result := AP_FP_Less(Err,100*MachineEpsilon);
end;


(*************************************************************************
Unset spline, i.e. initialize it with random garbage
*************************************************************************)
procedure UnsetSpline1D(var C : Spline1DInterpolant);
var
    X : TReal1DArray;
    Y : TReal1DArray;
    D : TReal1DArray;
begin
    SetLength(X, 2);
    SetLength(Y, 2);
    SetLength(D, 2);
    X[0] := -1;
    Y[0] := RandomReal;
    D[0] := RandomReal;
    X[1] := 1;
    Y[1] := RandomReal;
    D[1] := RandomReal;
    Spline1DBuildHermite(X, Y, D, 2, C);
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
function testsplineinterpolationunit_test_silent():Boolean;
begin
    Result := TestSplineInterpolation(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testsplineinterpolationunit_test():Boolean;
begin
    Result := TestSplineInterpolation(False);
end;


end.