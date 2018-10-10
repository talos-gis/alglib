unit testpsplineunit;
interface
uses Math, Sysutils, Ap, spline3, blas, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond, matinv, hblas, sblas, ortfac, rotations, bdsvd, svd, xblas, densesolver, linmin, minlbfgs, minlm, lsfit, apserv, spline1d, tsort, hsschur, evd, gammafunc, gq, gkq, autogk, pspline;

function TestPSplineInterpolation(Silent : Boolean):Boolean;
function testpsplineunit_test_silent():Boolean;
function testpsplineunit_test():Boolean;

implementation

procedure UnsetP2(var P : PSpline2Interpolant);forward;
procedure UnsetP3(var P : PSpline3Interpolant);forward;
procedure Unset1D(var X : TReal1DArray);forward;


function TestPSplineInterpolation(Silent : Boolean):Boolean;
var
    WasErrors : Boolean;
    P2Errors : Boolean;
    P3Errors : Boolean;
    NonStrictThreshold : Double;
    Threshold : Double;
    PassCount : AlglibInteger;
    LStep : Double;
    H : Double;
    MaxN : AlglibInteger;
    Periodicity : AlglibInteger;
    SKind : AlglibInteger;
    PKind : AlglibInteger;
    Periodic : Boolean;
    A : Double;
    B : Double;
    N : AlglibInteger;
    TmpN : AlglibInteger;
    I : AlglibInteger;
    K : AlglibInteger;
    VX : Double;
    VY : Double;
    VZ : Double;
    VX2 : Double;
    VY2 : Double;
    VZ2 : Double;
    VDX : Double;
    VDY : Double;
    VDZ : Double;
    VDX2 : Double;
    VDY2 : Double;
    VDZ2 : Double;
    VD2X : Double;
    VD2Y : Double;
    VD2Z : Double;
    VD2X2 : Double;
    VD2Y2 : Double;
    VD2Z2 : Double;
    V0 : Double;
    V1 : Double;
    X : TReal1DArray;
    Y : TReal1DArray;
    Z : TReal1DArray;
    T : TReal1DArray;
    T2 : TReal1DArray;
    T3 : TReal1DArray;
    XY : TReal2DArray;
    XYZ : TReal2DArray;
    P2 : PSpline2Interpolant;
    P3 : PSpline3Interpolant;
    S : Spline1DInterpolant;
    i_ : AlglibInteger;
begin
    WasErrors := False;
    PassCount := 20;
    LStep := 0.005;
    H := 0.00001;
    MaxN := 10;
    Threshold := 10000*MachineEpsilon;
    NonStrictThreshold := 0.00001;
    P2Errors := False;
    P3Errors := False;
    
    //
    // Test basic properties of 2- and 3-dimensional splines:
    // * PSpline2ParameterValues() properties
    // * values at nodes
    // * for periodic splines - periodicity properties
    //
    // Variables used:
    // * N              points count
    // * SKind          spline
    // * PKind          parameterization
    // * Periodicity    whether we have periodic spline or not
    //
    N:=2;
    while N<=MaxN do
    begin
        SKind:=0;
        while SKind<=2 do
        begin
            PKind:=0;
            while PKind<=2 do
            begin
                Periodicity:=0;
                while Periodicity<=1 do
                begin
                    Periodic := Periodicity=1;
                    
                    //
                    // skip unsupported combinations of parameters
                    //
                    if Periodic and (N<3) then
                    begin
                        Inc(Periodicity);
                        Continue;
                    end;
                    if Periodic and (SKind=0) then
                    begin
                        Inc(Periodicity);
                        Continue;
                    end;
                    if (N<5) and (SKind=0) then
                    begin
                        Inc(Periodicity);
                        Continue;
                    end;
                    
                    //
                    // init
                    //
                    SetLength(XY, N, 2);
                    SetLength(XYZ, N, 3);
                    TaskGenInt1DEquidist(-1, +1, N, T2, X);
                    for i_ := 0 to N-1 do
                    begin
                        XY[i_,0] := X[i_];
                    end;
                    for i_ := 0 to N-1 do
                    begin
                        XYZ[i_,0] := X[i_];
                    end;
                    TaskGenInt1DEquidist(-1, +1, N, T2, Y);
                    for i_ := 0 to N-1 do
                    begin
                        XY[i_,1] := Y[i_];
                    end;
                    for i_ := 0 to N-1 do
                    begin
                        XYZ[i_,1] := Y[i_];
                    end;
                    TaskGenInt1DEquidist(-1, +1, N, T2, Z);
                    for i_ := 0 to N-1 do
                    begin
                        XYZ[i_,2] := Z[i_];
                    end;
                    UnsetP2(P2);
                    UnsetP3(P3);
                    if Periodic then
                    begin
                        PSpline2BuildPeriodic(XY, N, SKind, PKind, P2);
                        PSpline3BuildPeriodic(XYZ, N, SKind, PKind, P3);
                    end
                    else
                    begin
                        PSpline2Build(XY, N, SKind, PKind, P2);
                        PSpline3Build(XYZ, N, SKind, PKind, P3);
                    end;
                    
                    //
                    // PSpline2ParameterValues() properties
                    //
                    PSpline2ParameterValues(P2, TmpN, T2);
                    if TmpN<>N then
                    begin
                        P2Errors := True;
                        Inc(Periodicity);
                        Continue;
                    end;
                    PSpline3ParameterValues(P3, TmpN, T3);
                    if TmpN<>N then
                    begin
                        P3Errors := True;
                        Inc(Periodicity);
                        Continue;
                    end;
                    P2Errors := P2Errors or AP_FP_Neq(T2[0],0);
                    P3Errors := P3Errors or AP_FP_Neq(T3[0],0);
                    I:=1;
                    while I<=N-1 do
                    begin
                        P2Errors := P2Errors or AP_FP_Less_Eq(T2[I],T2[I-1]);
                        P3Errors := P3Errors or AP_FP_Less_Eq(T3[I],T3[I-1]);
                        Inc(I);
                    end;
                    if Periodic then
                    begin
                        P2Errors := P2Errors or AP_FP_Greater_Eq(T2[N-1],1);
                        P3Errors := P3Errors or AP_FP_Greater_Eq(T3[N-1],1);
                    end
                    else
                    begin
                        P2Errors := P2Errors or AP_FP_Neq(T2[N-1],1);
                        P3Errors := P3Errors or AP_FP_Neq(T3[N-1],1);
                    end;
                    
                    //
                    // Now we have parameter values stored at T,
                    // and want to test whether the actully correspond to
                    // points
                    //
                    I:=0;
                    while I<=N-1 do
                    begin
                        
                        //
                        // 2-dimensional test
                        //
                        PSpline2Calc(P2, T2[I], VX, VY);
                        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VX-X[I]),Threshold);
                        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VY-Y[I]),Threshold);
                        
                        //
                        // 3-dimensional test
                        //
                        PSpline3Calc(P3, T3[I], VX, VY, VZ);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VX-X[I]),Threshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VY-Y[I]),Threshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VZ-Z[I]),Threshold);
                        Inc(I);
                    end;
                    
                    //
                    // Test periodicity (if needed)
                    //
                    if Periodic then
                    begin
                        
                        //
                        // periodicity at nodes
                        //
                        I:=0;
                        while I<=N-1 do
                        begin
                            
                            //
                            // 2-dimensional test
                            //
                            PSpline2Calc(P2, T2[I]+RandomInteger(10)-5, VX, VY);
                            P2Errors := P2Errors or AP_FP_Greater(AbsReal(VX-X[I]),Threshold);
                            P2Errors := P2Errors or AP_FP_Greater(AbsReal(VY-Y[I]),Threshold);
                            PSpline2Diff(P2, T2[I]+RandomInteger(10)-5, VX, VDX, VY, VDY);
                            P2Errors := P2Errors or AP_FP_Greater(AbsReal(VX-X[I]),Threshold);
                            P2Errors := P2Errors or AP_FP_Greater(AbsReal(VY-Y[I]),Threshold);
                            PSpline2Diff2(P2, T2[I]+RandomInteger(10)-5, VX, VDX, VD2X, VY, VDY, VD2Y);
                            P2Errors := P2Errors or AP_FP_Greater(AbsReal(VX-X[I]),Threshold);
                            P2Errors := P2Errors or AP_FP_Greater(AbsReal(VY-Y[I]),Threshold);
                            
                            //
                            // 3-dimensional test
                            //
                            PSpline3Calc(P3, T3[I]+RandomInteger(10)-5, VX, VY, VZ);
                            P3Errors := P3Errors or AP_FP_Greater(AbsReal(VX-X[I]),Threshold);
                            P3Errors := P3Errors or AP_FP_Greater(AbsReal(VY-Y[I]),Threshold);
                            P3Errors := P3Errors or AP_FP_Greater(AbsReal(VZ-Z[I]),Threshold);
                            PSpline3Diff(P3, T3[I]+RandomInteger(10)-5, VX, VDX, VY, VDY, VZ, VDZ);
                            P3Errors := P3Errors or AP_FP_Greater(AbsReal(VX-X[I]),Threshold);
                            P3Errors := P3Errors or AP_FP_Greater(AbsReal(VY-Y[I]),Threshold);
                            P3Errors := P3Errors or AP_FP_Greater(AbsReal(VZ-Z[I]),Threshold);
                            PSpline3Diff2(P3, T3[I]+RandomInteger(10)-5, VX, VDX, VD2X, VY, VDY, VD2Y, VZ, VDZ, VD2Z);
                            P3Errors := P3Errors or AP_FP_Greater(AbsReal(VX-X[I]),Threshold);
                            P3Errors := P3Errors or AP_FP_Greater(AbsReal(VY-Y[I]),Threshold);
                            P3Errors := P3Errors or AP_FP_Greater(AbsReal(VZ-Z[I]),Threshold);
                            Inc(I);
                        end;
                        
                        //
                        // periodicity between nodes
                        //
                        V0 := RandomReal;
                        PSpline2Calc(P2, V0, VX, VY);
                        PSpline2Calc(P2, V0+RandomInteger(10)-5, VX2, VY2);
                        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VX-VX2),Threshold);
                        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VY-VY2),Threshold);
                        PSpline3Calc(P3, V0, VX, VY, VZ);
                        PSpline3Calc(P3, V0+RandomInteger(10)-5, VX2, VY2, VZ2);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VX-VX2),Threshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VY-VY2),Threshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VZ-VZ2),Threshold);
                        
                        //
                        // near-boundary test for continuity of function values and derivatives:
                        // 2-dimensional curve
                        //
                        Assert((SKind=1) or (SKind=2), 'TEST: unexpected spline type!');
                        V0 := 100*MachineEpsilon;
                        V1 := 1-V0;
                        PSpline2Calc(P2, V0, VX, VY);
                        PSpline2Calc(P2, V1, VX2, VY2);
                        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VX-VX2),Threshold);
                        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VY-VY2),Threshold);
                        PSpline2Diff(P2, V0, VX, VDX, VY, VDY);
                        PSpline2Diff(P2, V1, VX2, VDX2, VY2, VDY2);
                        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VX-VX2),Threshold);
                        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VY-VY2),Threshold);
                        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VDX-VDX2),NonStrictThreshold);
                        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VDY-VDY2),NonStrictThreshold);
                        PSpline2Diff2(P2, V0, VX, VDX, VD2X, VY, VDY, VD2Y);
                        PSpline2Diff2(P2, V1, VX2, VDX2, VD2X2, VY2, VDY2, VD2Y2);
                        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VX-VX2),Threshold);
                        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VY-VY2),Threshold);
                        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VDX-VDX2),NonStrictThreshold);
                        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VDY-VDY2),NonStrictThreshold);
                        if SKind=2 then
                        begin
                            
                            //
                            // second derivative test only for cubic splines
                            //
                            P2Errors := P2Errors or AP_FP_Greater(AbsReal(VD2X-VD2X2),NonStrictThreshold);
                            P2Errors := P2Errors or AP_FP_Greater(AbsReal(VD2Y-VD2Y2),NonStrictThreshold);
                        end;
                        
                        //
                        // near-boundary test for continuity of function values and derivatives:
                        // 3-dimensional curve
                        //
                        Assert((SKind=1) or (SKind=2), 'TEST: unexpected spline type!');
                        V0 := 100*MachineEpsilon;
                        V1 := 1-V0;
                        PSpline3Calc(P3, V0, VX, VY, VZ);
                        PSpline3Calc(P3, V1, VX2, VY2, VZ2);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VX-VX2),Threshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VY-VY2),Threshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VZ-VZ2),Threshold);
                        PSpline3Diff(P3, V0, VX, VDX, VY, VDY, VZ, VDZ);
                        PSpline3Diff(P3, V1, VX2, VDX2, VY2, VDY2, VZ2, VDZ2);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VX-VX2),Threshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VY-VY2),Threshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VZ-VZ2),Threshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VDX-VDX2),NonStrictThreshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VDY-VDY2),NonStrictThreshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VDZ-VDZ2),NonStrictThreshold);
                        PSpline3Diff2(P3, V0, VX, VDX, VD2X, VY, VDY, VD2Y, VZ, VDZ, VD2Z);
                        PSpline3Diff2(P3, V1, VX2, VDX2, VD2X2, VY2, VDY2, VD2Y2, VZ2, VDZ2, VD2Z2);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VX-VX2),Threshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VY-VY2),Threshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VZ-VZ2),Threshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VDX-VDX2),NonStrictThreshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VDY-VDY2),NonStrictThreshold);
                        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VDZ-VDZ2),NonStrictThreshold);
                        if SKind=2 then
                        begin
                            
                            //
                            // second derivative test only for cubic splines
                            //
                            P3Errors := P3Errors or AP_FP_Greater(AbsReal(VD2X-VD2X2),NonStrictThreshold);
                            P3Errors := P3Errors or AP_FP_Greater(AbsReal(VD2Y-VD2Y2),NonStrictThreshold);
                            P3Errors := P3Errors or AP_FP_Greater(AbsReal(VD2Z-VD2Z2),NonStrictThreshold);
                        end;
                    end;
                    Inc(Periodicity);
                end;
                Inc(PKind);
            end;
            Inc(SKind);
        end;
        Inc(N);
    end;
    
    //
    // Test differentiation, tangents, calculation between nodes.
    //
    // Because differentiation is done in parameterization/spline/periodicity
    // oblivious manner, we don't have to test all possible combinations
    // of spline types and parameterizations.
    //
    // Actually we test special combination with properties which allow us
    // to easily solve this problem:
    // * 2 (3) variables
    // * first variable is sampled from equidistant grid on [0,1]
    // * other variables are random
    // * uniform parameterization is used
    // * periodicity - none
    // * spline type - any (we use cubic splines)
    // Same problem allows us to test calculation BETWEEN nodes.
    //
    N:=2;
    while N<=MaxN do
    begin
        
        //
        // init
        //
        SetLength(XY, N, 2);
        SetLength(XYZ, N, 3);
        TaskGenInt1DEquidist(0, +1, N, T, X);
        for i_ := 0 to N-1 do
        begin
            XY[i_,0] := X[i_];
        end;
        for i_ := 0 to N-1 do
        begin
            XYZ[i_,0] := X[i_];
        end;
        TaskGenInt1DEquidist(0, +1, N, T, Y);
        for i_ := 0 to N-1 do
        begin
            XY[i_,1] := Y[i_];
        end;
        for i_ := 0 to N-1 do
        begin
            XYZ[i_,1] := Y[i_];
        end;
        TaskGenInt1DEquidist(0, +1, N, T, Z);
        for i_ := 0 to N-1 do
        begin
            XYZ[i_,2] := Z[i_];
        end;
        UnsetP2(P2);
        UnsetP3(P3);
        PSpline2Build(XY, N, 2, 0, P2);
        PSpline3Build(XYZ, N, 2, 0, P3);
        
        //
        // Test 2D/3D spline:
        // * build non-parametric cubic spline from T and X/Y
        // * calculate its value and derivatives at V0
        // * compare with Spline2Calc/Spline2Diff/Spline2Diff2
        // Because of task properties both variants should
        // return same answer.
        //
        V0 := RandomReal;
        Spline1DBuildCubic(T, X, N, 0, 0.0, 0, 0.0, S);
        Spline1DDiff(S, V0, VX2, VDX2, VD2X2);
        Spline1DBuildCubic(T, Y, N, 0, 0.0, 0, 0.0, S);
        Spline1DDiff(S, V0, VY2, VDY2, VD2Y2);
        Spline1DBuildCubic(T, Z, N, 0, 0.0, 0, 0.0, S);
        Spline1DDiff(S, V0, VZ2, VDZ2, VD2Z2);
        
        //
        // 2D test
        //
        PSpline2Calc(P2, V0, VX, VY);
        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VX-VX2),Threshold);
        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VY-VY2),Threshold);
        PSpline2Diff(P2, V0, VX, VDX, VY, VDY);
        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VX-VX2),Threshold);
        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VY-VY2),Threshold);
        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VDX-VDX2),Threshold);
        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VDY-VDY2),Threshold);
        PSpline2Diff2(P2, V0, VX, VDX, VD2X, VY, VDY, VD2Y);
        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VX-VX2),Threshold);
        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VY-VY2),Threshold);
        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VDX-VDX2),Threshold);
        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VDY-VDY2),Threshold);
        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VD2X-VD2X2),Threshold);
        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VD2Y-VD2Y2),Threshold);
        
        //
        // 3D test
        //
        PSpline3Calc(P3, V0, VX, VY, VZ);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VX-VX2),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VY-VY2),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VZ-VZ2),Threshold);
        PSpline3Diff(P3, V0, VX, VDX, VY, VDY, VZ, VDZ);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VX-VX2),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VY-VY2),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VZ-VZ2),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VDX-VDX2),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VDY-VDY2),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VDZ-VDZ2),Threshold);
        PSpline3Diff2(P3, V0, VX, VDX, VD2X, VY, VDY, VD2Y, VZ, VDZ, VD2Z);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VX-VX2),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VY-VY2),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VZ-VZ2),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VDX-VDX2),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VDY-VDY2),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VDZ-VDZ2),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VD2X-VD2X2),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VD2Y-VD2Y2),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VD2Z-VD2Z2),Threshold);
        
        //
        // Test tangents for 2D/3D
        //
        PSpline2Tangent(P2, V0, VX, VY);
        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VX-VDX2/SafePythag2(VDX2, VDY2)),Threshold);
        P2Errors := P2Errors or AP_FP_Greater(AbsReal(VY-VDY2/SafePythag2(VDX2, VDY2)),Threshold);
        PSpline3Tangent(P3, V0, VX, VY, VZ);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VX-VDX2/SafePythag3(VDX2, VDY2, VDZ2)),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VY-VDY2/SafePythag3(VDX2, VDY2, VDZ2)),Threshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(VZ-VDZ2/SafePythag3(VDX2, VDY2, VDZ2)),Threshold);
        Inc(N);
    end;
    
    //
    // Arc length test.
    //
    // Simple problem with easy solution (points on a straight line with
    // uniform parameterization).
    //
    N:=2;
    while N<=MaxN do
    begin
        SetLength(XY, N, 2);
        SetLength(XYZ, N, 3);
        I:=0;
        while I<=N-1 do
        begin
            XY[I,0] := I;
            XY[I,1] := I;
            XYZ[I,0] := I;
            XYZ[I,1] := I;
            XYZ[I,2] := I;
            Inc(I);
        end;
        PSpline2Build(XY, N, 1, 0, P2);
        PSpline3Build(XYZ, N, 1, 0, P3);
        A := RandomReal;
        B := RandomReal;
        P2Errors := P2Errors or AP_FP_Greater(AbsReal(PSpline2ArcLength(P2, A, B)-(B-A)*Sqrt(2)*(N-1)),NonStrictThreshold);
        P3Errors := P3Errors or AP_FP_Greater(AbsReal(PSpline3ArcLength(P3, A, B)-(B-A)*Sqrt(3)*(N-1)),NonStrictThreshold);
        Inc(N);
    end;
    
    //
    // report
    //
    WasErrors := P2Errors or P3Errors;
    if  not Silent then
    begin
        Write(Format('TESTING SPLINE INTERPOLATION'#13#10'',[]));
        
        //
        // Normal tests
        //
        Write(Format('2D TEST:                                 ',[]));
        if P2Errors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('3D TEST:                                 ',[]));
        if P3Errors then
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
Unset spline, i.e. initialize it with random garbage
*************************************************************************)
procedure UnsetP2(var P : PSpline2Interpolant);
var
    XY : TReal2DArray;
begin
    SetLength(XY, 2, 2);
    XY[0,0] := -1;
    XY[0,1] := -1;
    XY[1,0] := +1;
    XY[1,1] := +1;
    PSpline2Build(XY, 2, 1, 0, P);
end;


(*************************************************************************
Unset spline, i.e. initialize it with random garbage
*************************************************************************)
procedure UnsetP3(var P : PSpline3Interpolant);
var
    XY : TReal2DArray;
begin
    SetLength(XY, 2, 3);
    XY[0,0] := -1;
    XY[0,1] := -1;
    XY[0,2] := -1;
    XY[1,0] := +1;
    XY[1,1] := +1;
    XY[1,2] := +1;
    PSpline3Build(XY, 2, 1, 0, P);
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
Silent unit test
*************************************************************************)
function testpsplineunit_test_silent():Boolean;
begin
    Result := TestPSplineInterpolation(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testpsplineunit_test():Boolean;
begin
    Result := TestPSplineInterpolation(False);
end;


end.