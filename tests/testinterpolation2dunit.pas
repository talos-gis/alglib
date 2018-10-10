unit testinterpolation2dunit;
interface
uses Math, Sysutils, Ap, spline3, blas, reflections, creflections, hqrnd, matgen, trinverse, ablasf, ablas, trfac, bidiagonal, qr, lq, rotations, bdsvd, svd, trlinsolve, safesolve, rcond, tsort, xblas, densesolver, lbfgs, minlm, leastsquares, lsfit, spline1d, spline2d;

function Test2DInterpolation(Silent : Boolean):Boolean;
function testinterpolation2dunit_test_silent():Boolean;
function testinterpolation2dunit_test():Boolean;

implementation

procedure LConst(const C : Spline2DInterpolant;
     const LX : TReal1DArray;
     const LY : TReal1DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     LStep : Double;
     var LC : Double;
     var LCX : Double;
     var LCY : Double;
     var LCXY : Double);forward;
procedure TwoDNumDer(const C : Spline2DInterpolant;
     X : Double;
     Y : Double;
     H : Double;
     var F : Double;
     var FX : Double;
     var FY : Double;
     var FXY : Double);forward;
function TestUnpack(const C : Spline2DInterpolant;
     const LX : TReal1DArray;
     const LY : TReal1DArray):Boolean;forward;
function TestLinTrans(const C : Spline2DInterpolant;
     AX : Double;
     BX : Double;
     AY : Double;
     BY : Double):Boolean;forward;
procedure UnsetSpline2D(var C : Spline2DInterpolant);forward;


function Test2DInterpolation(Silent : Boolean):Boolean;
var
    WasErrors : Boolean;
    BLErrors : Boolean;
    BCErrors : Boolean;
    DSErrors : Boolean;
    CPErrors : Boolean;
    UPErrors : Boolean;
    LTErrors : Boolean;
    SYErrors : Boolean;
    RLErrors : Boolean;
    RCErrors : Boolean;
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    JobType : AlglibInteger;
    LStep : Double;
    H : Double;
    X : TReal1DArray;
    Y : TReal1DArray;
    C : Spline2DInterpolant;
    C2 : Spline2DInterpolant;
    LX : TReal1DArray;
    LY : TReal1DArray;
    F : TReal2DArray;
    FR : TReal2DArray;
    FT : TReal2DArray;
    AX : Double;
    AY : Double;
    BX : Double;
    BY : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    N : AlglibInteger;
    M : AlglibInteger;
    N2 : AlglibInteger;
    M2 : AlglibInteger;
    Err : Double;
    T : Double;
    T1 : Double;
    T2 : Double;
    L1 : Double;
    L1X : Double;
    L1Y : Double;
    L1XY : Double;
    L2 : Double;
    L2X : Double;
    L2Y : Double;
    L2XY : Double;
    FM : Double;
    F1 : Double;
    F2 : Double;
    F3 : Double;
    F4 : Double;
    V1 : Double;
    V1X : Double;
    V1Y : Double;
    V1XY : Double;
    V2 : Double;
    V2X : Double;
    V2Y : Double;
    V2XY : Double;
    MF : Double;
    RA : TReal1DArray;
    RA2 : TReal1DArray;
    RALen : AlglibInteger;
begin
    WasErrors := False;
    PassCount := 10;
    H := 0.00001;
    LStep := 0.001;
    BLErrors := False;
    BCErrors := False;
    DSErrors := False;
    CPErrors := False;
    UPErrors := False;
    LTErrors := False;
    SYErrors := False;
    RLErrors := False;
    RCErrors := False;
    
    //
    // Test: bilinear, bicubic
    //
    N:=2;
    while N<=7 do
    begin
        M:=2;
        while M<=7 do
        begin
            SetLength(X, N-1+1);
            SetLength(Y, M-1+1);
            SetLength(LX, 2*N-2+1);
            SetLength(LY, 2*M-2+1);
            SetLength(F, M-1+1, N-1+1);
            SetLength(FT, N-1+1, M-1+1);
            Pass:=1;
            while Pass<=PassCount do
            begin
                
                //
                // Prepare task:
                // * X and Y stores grid
                // * F stores function values
                // * LX and LY stores twice dense grid (for Lipschitz testing)
                //
                AX := -1-RandomReal;
                BX := +1+RandomReal;
                AY := -1-RandomReal;
                BY := +1+RandomReal;
                J:=0;
                while J<=N-1 do
                begin
                    X[J] := 0.5*(BX+AX)-0.5*(BX-AX)*Cos(PI*(2*J+1)/(2*n));
                    if J=0 then
                    begin
                        X[J] := AX;
                    end;
                    if J=N-1 then
                    begin
                        X[J] := BX;
                    end;
                    LX[2*J] := X[J];
                    if J>0 then
                    begin
                        LX[2*J-1] := 0.5*(X[J]+X[J-1]);
                    end;
                    Inc(J);
                end;
                J:=0;
                while J<=N-1 do
                begin
                    K := RandomInteger(N);
                    if K<>J then
                    begin
                        T := X[J];
                        X[J] := X[K];
                        X[K] := T;
                    end;
                    Inc(J);
                end;
                I:=0;
                while I<=M-1 do
                begin
                    Y[I] := 0.5*(BY+AY)-0.5*(BY-AY)*Cos(PI*(2*i+1)/(2*M));
                    if I=0 then
                    begin
                        Y[I] := AY;
                    end;
                    if I=M-1 then
                    begin
                        Y[I] := BY;
                    end;
                    LY[2*I] := Y[I];
                    if I>0 then
                    begin
                        LY[2*I-1] := 0.5*(Y[I]+Y[I-1]);
                    end;
                    Inc(I);
                end;
                I:=0;
                while I<=M-1 do
                begin
                    K := RandomInteger(M);
                    if K<>I then
                    begin
                        T := Y[I];
                        Y[I] := Y[K];
                        Y[K] := T;
                    end;
                    Inc(I);
                end;
                I:=0;
                while I<=M-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        F[I,J] := Exp(0.6*X[J])-Exp(-0.3*Y[I]+0.08*X[J])+2*Cos(Pi*(X[J]+1.2*Y[I]))+0.1*Cos(20*X[J]+15*Y[I]);
                        Inc(J);
                    end;
                    Inc(I);
                end;
                
                //
                // Test bilinear interpolation:
                // * interpolation at the nodes
                // * linearity
                // * continuity
                // * differentiation in the inner points
                //
                Spline2DBuildBilinear(X, Y, F, M, N, C);
                Err := 0;
                I:=0;
                while I<=M-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        Err := Max(Err, AbsReal(F[I,J]-Spline2DCalc(C, X[J], Y[I])));
                        Inc(J);
                    end;
                    Inc(I);
                end;
                BLErrors := BLErrors or AP_FP_Greater(Err,10000*MachineEpsilon);
                Err := 0;
                I:=0;
                while I<=M-2 do
                begin
                    J:=0;
                    while J<=N-2 do
                    begin
                        
                        //
                        // Test for linearity between grid points
                        // (test point - geometric center of the cell)
                        //
                        FM := Spline2DCalc(C, LX[2*J+1], LY[2*I+1]);
                        F1 := Spline2DCalc(C, LX[2*J], LY[2*I]);
                        F2 := Spline2DCalc(C, LX[2*J+2], LY[2*I]);
                        F3 := Spline2DCalc(C, LX[2*J+2], LY[2*I+2]);
                        F4 := Spline2DCalc(C, LX[2*J], LY[2*I+2]);
                        Err := Max(Err, AbsReal(0.25*(F1+F2+F3+F4)-FM));
                        Inc(J);
                    end;
                    Inc(I);
                end;
                BLErrors := BLErrors or AP_FP_Greater(Err,10000*MachineEpsilon);
                LConst(C, LX, LY, M, N, LStep, L1, L1X, L1Y, L1XY);
                LConst(C, LX, LY, M, N, LStep/3, L2, L2X, L2Y, L2XY);
                BLErrors := BLErrors or AP_FP_Greater(L2/L1,1.2);
                Err := 0;
                I:=0;
                while I<=M-2 do
                begin
                    J:=0;
                    while J<=N-2 do
                    begin
                        Spline2DDiff(C, LX[2*J+1], LY[2*I+1], V1, V1X, V1Y, V1XY);
                        TwoDNumDer(C, LX[2*J+1], LY[2*I+1], H, V2, V2X, V2Y, V2XY);
                        Err := Max(Err, AbsReal(V1-V2));
                        Err := Max(Err, AbsReal(V1X-V2X));
                        Err := Max(Err, AbsReal(V1Y-V2Y));
                        Err := Max(Err, AbsReal(V1XY-V2XY));
                        Inc(J);
                    end;
                    Inc(I);
                end;
                DSErrors := DSErrors or AP_FP_Greater(Err,1.0E-3);
                UPErrors := UPErrors or  not TestUnpack(C, LX, LY);
                LTErrors := LTErrors or  not TestLinTrans(C, AX, BX, AY, BY);
                
                //
                // Test bicubic interpolation.
                // * interpolation at the nodes
                // * smoothness
                // * differentiation
                //
                Spline2DBuildBicubic(X, Y, F, M, N, C);
                Err := 0;
                I:=0;
                while I<=M-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        Err := Max(Err, AbsReal(F[I,J]-Spline2DCalc(C, X[J], Y[I])));
                        Inc(J);
                    end;
                    Inc(I);
                end;
                BCErrors := BCErrors or AP_FP_Greater(Err,10000*MachineEpsilon);
                LConst(C, LX, LY, M, N, LStep, L1, L1X, L1Y, L1XY);
                LConst(C, LX, LY, M, N, LStep/3, L2, L2X, L2Y, L2XY);
                BCErrors := BCErrors or AP_FP_Greater(L2/L1,1.2);
                BCErrors := BCErrors or AP_FP_Greater(L2X/L1X,1.2);
                BCErrors := BCErrors or AP_FP_Greater(L2Y/L1Y,1.2);
                if AP_FP_Greater(L2XY,0.01) and AP_FP_Greater(L1XY,0.01) then
                begin
                    
                    //
                    // Cross-derivative continuity is tested only when
                    // bigger than 0.01. When the task size is too
                    // small, the d2F/dXdY is nearly zero and Lipschitz
                    // constant ratio is ill-conditioned.
                    //
                    BCErrors := BCErrors or AP_FP_Greater(L2XY/L1XY,1.2);
                end;
                Err := 0;
                I:=0;
                while I<=2*M-2 do
                begin
                    J:=0;
                    while J<=2*N-2 do
                    begin
                        Spline2DDiff(C, LX[J], LY[I], V1, V1X, V1Y, V1XY);
                        TwoDNumDer(C, LX[J], LY[I], H, V2, V2X, V2Y, V2XY);
                        Err := Max(Err, AbsReal(V1-V2));
                        Err := Max(Err, AbsReal(V1X-V2X));
                        Err := Max(Err, AbsReal(V1Y-V2Y));
                        Err := Max(Err, AbsReal(V1XY-V2XY));
                        Inc(J);
                    end;
                    Inc(I);
                end;
                DSErrors := DSErrors or AP_FP_Greater(Err,1.0E-3);
                UPErrors := UPErrors or  not TestUnpack(C, LX, LY);
                LTErrors := LTErrors or  not TestLinTrans(C, AX, BX, AY, BY);
                
                //
                // Copy/Serialise test
                //
                if AP_FP_Greater(RandomReal,0.5) then
                begin
                    Spline2DBuildBicubic(X, Y, F, M, N, C);
                end
                else
                begin
                    Spline2DBuildBilinear(X, Y, F, M, N, C);
                end;
                UnsetSpline2D(C2);
                Spline2DCopy(C, C2);
                Err := 0;
                I:=1;
                while I<=5 do
                begin
                    T1 := AX+(BX-AX)*RandomReal;
                    T2 := AY+(BY-AY)*RandomReal;
                    Err := Max(Err, AbsReal(Spline2DCalc(C, T1, T2)-Spline2DCalc(C2, T1, T2)));
                    Inc(I);
                end;
                CPErrors := CPErrors or AP_FP_Greater(Err,10000*MachineEpsilon);
                UnsetSpline2D(C2);
                Spline2DSerialize(C, RA, RALen);
                SetLength(RA2, RALen);
                APVMove(@RA2[0], 0, RALen-1, @RA[0], 0, RALen-1);
                Spline2DUnserialize(RA2, C2);
                Err := 0;
                I:=1;
                while I<=5 do
                begin
                    T1 := AX+(BX-AX)*RandomReal;
                    T2 := AY+(BY-AY)*RandomReal;
                    Err := Max(Err, AbsReal(Spline2DCalc(C, T1, T2)-Spline2DCalc(C2, T1, T2)));
                    Inc(I);
                end;
                CPErrors := CPErrors or AP_FP_Greater(Err,10000*MachineEpsilon);
                
                //
                // Special symmetry test
                //
                Err := 0;
                JobType:=0;
                while JobType<=1 do
                begin
                    
                    //
                    // Prepare
                    //
                    I:=0;
                    while I<=M-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            FT[J,I] := F[I,J];
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    if JobType=0 then
                    begin
                        Spline2DBuildBilinear(X, Y, F, M, N, C);
                        Spline2DBuildBilinear(Y, X, FT, N, M, C2);
                    end
                    else
                    begin
                        Spline2DBuildBicubic(X, Y, F, M, N, C);
                        Spline2DBuildBicubic(Y, X, FT, N, M, C2);
                    end;
                    
                    //
                    // Test
                    //
                    I:=1;
                    while I<=10 do
                    begin
                        T1 := AX+(BX-AX)*RandomReal;
                        T2 := AY+(BY-AY)*RandomReal;
                        Err := Max(Err, AbsReal(Spline2DCalc(C, T1, T2)-Spline2DCalc(C2, T2, T1)));
                        Inc(I);
                    end;
                    Inc(JobType);
                end;
                SYErrors := SYErrors or AP_FP_Greater(Err,10000*MachineEpsilon);
                Inc(Pass);
            end;
            Inc(M);
        end;
        Inc(N);
    end;
    
    //
    // Test resample
    //
    M:=2;
    while M<=6 do
    begin
        N:=2;
        while N<=6 do
        begin
            SetLength(F, M-1+1, N-1+1);
            SetLength(X, N-1+1);
            SetLength(Y, M-1+1);
            J:=0;
            while J<=N-1 do
            begin
                X[J] := AP_Double(J)/(N-1);
                Inc(J);
            end;
            I:=0;
            while I<=M-1 do
            begin
                Y[I] := AP_Double(I)/(M-1);
                Inc(I);
            end;
            I:=0;
            while I<=M-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    F[I,J] := Exp(0.6*X[J])-Exp(-0.3*Y[I]+0.08*X[J])+2*Cos(Pi*(X[J]+1.2*Y[I]))+0.1*Cos(20*X[J]+15*Y[I]);
                    Inc(J);
                end;
                Inc(I);
            end;
            M2:=2;
            while M2<=6 do
            begin
                N2:=2;
                while N2<=6 do
                begin
                    Pass:=1;
                    while Pass<=PassCount do
                    begin
                        JobType:=0;
                        while JobType<=1 do
                        begin
                            if JobType=0 then
                            begin
                                Spline2DResampleBilinear(F, M, N, FR, M2, N2);
                                Spline2DBuildBilinear(X, Y, F, M, N, C);
                            end;
                            if JobType=1 then
                            begin
                                Spline2DResampleBicubic(F, M, N, FR, M2, N2);
                                Spline2DBuildBicubic(X, Y, F, M, N, C);
                            end;
                            Err := 0;
                            MF := 0;
                            I:=0;
                            while I<=M2-1 do
                            begin
                                J:=0;
                                while J<=N2-1 do
                                begin
                                    V1 := Spline2DCalc(C, AP_Double(J)/(N2-1), AP_Double(I)/(M2-1));
                                    V2 := FR[I,J];
                                    Err := Max(Err, AbsReal(V1-V2));
                                    MF := Max(MF, AbsReal(V1));
                                    Inc(J);
                                end;
                                Inc(I);
                            end;
                            if JobType=0 then
                            begin
                                RLErrors := RLErrors or AP_FP_Greater(Err/MF,10000*MachineEpsilon);
                            end;
                            if JobType=1 then
                            begin
                                RCErrors := RCErrors or AP_FP_Greater(Err/MF,10000*MachineEpsilon);
                            end;
                            Inc(JobType);
                        end;
                        Inc(Pass);
                    end;
                    Inc(N2);
                end;
                Inc(M2);
            end;
            Inc(N);
        end;
        Inc(M);
    end;
    
    //
    // report
    //
    WasErrors := BLErrors or BCErrors or DSErrors or CPErrors or UPErrors or LTErrors or SYErrors or RLErrors or RCErrors;
    if  not Silent then
    begin
        Write(Format('TESTING 2D INTERPOLATION'#13#10'',[]));
        
        //
        // Normal tests
        //
        Write(Format('BILINEAR TEST:                           ',[]));
        if BLErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('BICUBIC TEST:                            ',[]));
        if BCErrors then
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
        Write(Format('COPY/SERIALIZE TEST:                     ',[]));
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
        Write(Format('SPECIAL SYMMETRY TEST:                   ',[]));
        if SYErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('BILINEAR RESAMPLING TEST:                ',[]));
        if RLErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('BICUBIC RESAMPLING TEST:                 ',[]));
        if RCErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        
        //
        // Summary
        //
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
procedure LConst(const C : Spline2DInterpolant;
     const LX : TReal1DArray;
     const LY : TReal1DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     LStep : Double;
     var LC : Double;
     var LCX : Double;
     var LCY : Double;
     var LCXY : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
    F1 : Double;
    F2 : Double;
    F3 : Double;
    F4 : Double;
    FX1 : Double;
    FX2 : Double;
    FX3 : Double;
    FX4 : Double;
    FY1 : Double;
    FY2 : Double;
    FY3 : Double;
    FY4 : Double;
    FXY1 : Double;
    FXY2 : Double;
    FXY3 : Double;
    FXY4 : Double;
    S2LStep : Double;
begin
    LC := 0;
    LCX := 0;
    LCY := 0;
    LCXY := 0;
    S2LStep := Sqrt(2)*LStep;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            
            //
            // Calculate
            //
            TwoDNumDer(C, LX[J]-LStep/2, LY[I]-LStep/2, LStep/4, F1, FX1, FY1, FXY1);
            TwoDNumDer(C, LX[J]+LStep/2, LY[I]-LStep/2, LStep/4, F2, FX2, FY2, FXY2);
            TwoDNumDer(C, LX[J]+LStep/2, LY[I]+LStep/2, LStep/4, F3, FX3, FY3, FXY3);
            TwoDNumDer(C, LX[J]-LStep/2, LY[I]+LStep/2, LStep/4, F4, FX4, FY4, FXY4);
            
            //
            // Lipschitz constant for the function itself
            //
            LC := Max(LC, AbsReal((F1-F2)/LStep));
            LC := Max(LC, AbsReal((F2-F3)/LStep));
            LC := Max(LC, AbsReal((F3-F4)/LStep));
            LC := Max(LC, AbsReal((F4-F1)/LStep));
            LC := Max(LC, AbsReal((F1-F3)/S2LStep));
            LC := Max(LC, AbsReal((F2-F4)/S2LStep));
            
            //
            // Lipschitz constant for the first derivative
            //
            LCX := Max(LCX, AbsReal((FX1-FX2)/LStep));
            LCX := Max(LCX, AbsReal((FX2-FX3)/LStep));
            LCX := Max(LCX, AbsReal((FX3-FX4)/LStep));
            LCX := Max(LCX, AbsReal((FX4-FX1)/LStep));
            LCX := Max(LCX, AbsReal((FX1-FX3)/S2LStep));
            LCX := Max(LCX, AbsReal((FX2-FX4)/S2LStep));
            
            //
            // Lipschitz constant for the first derivative
            //
            LCY := Max(LCY, AbsReal((FY1-FY2)/LStep));
            LCY := Max(LCY, AbsReal((FY2-FY3)/LStep));
            LCY := Max(LCY, AbsReal((FY3-FY4)/LStep));
            LCY := Max(LCY, AbsReal((FY4-FY1)/LStep));
            LCY := Max(LCY, AbsReal((FY1-FY3)/S2LStep));
            LCY := Max(LCY, AbsReal((FY2-FY4)/S2LStep));
            
            //
            // Lipschitz constant for the cross-derivative
            //
            LCXY := Max(LCXY, AbsReal((FXY1-FXY2)/LStep));
            LCXY := Max(LCXY, AbsReal((FXY2-FXY3)/LStep));
            LCXY := Max(LCXY, AbsReal((FXY3-FXY4)/LStep));
            LCXY := Max(LCXY, AbsReal((FXY4-FXY1)/LStep));
            LCXY := Max(LCXY, AbsReal((FXY1-FXY3)/S2LStep));
            LCXY := Max(LCXY, AbsReal((FXY2-FXY4)/S2LStep));
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Numerical differentiation.
*************************************************************************)
procedure TwoDNumDer(const C : Spline2DInterpolant;
     X : Double;
     Y : Double;
     H : Double;
     var F : Double;
     var FX : Double;
     var FY : Double;
     var FXY : Double);
begin
    F := Spline2DCalc(C, X, Y);
    FX := (Spline2DCalc(C, X+H, Y)-Spline2DCalc(C, X-H, Y))/(2*H);
    FY := (Spline2DCalc(C, X, Y+H)-Spline2DCalc(C, X, Y-H))/(2*H);
    FXY := (Spline2DCalc(C, X+H, Y+H)-Spline2DCalc(C, X-H, Y+H)-Spline2DCalc(C, X+H, Y-H)+Spline2DCalc(C, X-H, Y-H))/AP_Sqr(2*H);
end;


(*************************************************************************
Unpack test
*************************************************************************)
function TestUnpack(const C : Spline2DInterpolant;
     const LX : TReal1DArray;
     const LY : TReal1DArray):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    N : AlglibInteger;
    M : AlglibInteger;
    CI : AlglibInteger;
    CJ : AlglibInteger;
    P : AlglibInteger;
    Err : Double;
    TX : Double;
    TY : Double;
    V1 : Double;
    V2 : Double;
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    Tbl : TReal2DArray;
begin
    PassCount := 20;
    Err := 0;
    Spline2DUnpack(C, M, N, Tbl);
    I:=0;
    while I<=M-2 do
    begin
        J:=0;
        while J<=N-2 do
        begin
            Pass:=1;
            while Pass<=PassCount do
            begin
                P := (N-1)*I+J;
                TX := (0.001+0.999*RandomReal)*(Tbl[P,1]-Tbl[P,0]);
                TY := (0.001+0.999*RandomReal)*(Tbl[P,3]-Tbl[P,2]);
                
                //
                // Interpolation properties
                //
                V1 := 0;
                CI:=0;
                while CI<=3 do
                begin
                    CJ:=0;
                    while CJ<=3 do
                    begin
                        V1 := V1+Tbl[P,4+CI*4+CJ]*Power(TX, CI)*Power(TY, CJ);
                        Inc(CJ);
                    end;
                    Inc(CI);
                end;
                V2 := Spline2DCalc(C, Tbl[P,0]+TX, Tbl[P,2]+TY);
                Err := Max(Err, AbsReal(V1-V2));
                
                //
                // Grid correctness
                //
                Err := Max(Err, AbsReal(LX[2*J]-Tbl[P,0]));
                Err := Max(Err, AbsReal(LX[2*(J+1)]-Tbl[P,1]));
                Err := Max(Err, AbsReal(LY[2*I]-Tbl[P,2]));
                Err := Max(Err, AbsReal(LY[2*(I+1)]-Tbl[P,3]));
                Inc(Pass);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    Result := AP_FP_Less(Err,10000*MachineEpsilon);
end;


(*************************************************************************
LinTrans test
*************************************************************************)
function TestLinTrans(const C : Spline2DInterpolant;
     AX : Double;
     BX : Double;
     AY : Double;
     BY : Double):Boolean;
var
    Err : Double;
    A1 : Double;
    A2 : Double;
    B1 : Double;
    B2 : Double;
    TX : Double;
    TY : Double;
    VX : Double;
    VY : Double;
    V1 : Double;
    V2 : Double;
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    XJob : AlglibInteger;
    YJob : AlglibInteger;
    C2 : Spline2DInterpolant;
begin
    PassCount := 5;
    Err := 0;
    XJob:=0;
    while XJob<=1 do
    begin
        YJob:=0;
        while YJob<=1 do
        begin
            Pass:=1;
            while Pass<=PassCount do
            begin
                
                //
                // Prepare
                //
                repeat
                    A1 := 2*RandomReal-1;
                until AP_FP_Neq(A1,0);
                A1 := A1*XJob;
                B1 := 2*RandomReal-1;
                repeat
                    A2 := 2*RandomReal-1;
                until AP_FP_Neq(A2,0);
                A2 := A2*YJob;
                B2 := 2*RandomReal-1;
                
                //
                // Test XY
                //
                Spline2DCopy(C, C2);
                Spline2DLinTransXY(C2, A1, B1, A2, B2);
                TX := AX+RandomReal*(BX-AX);
                TY := AY+RandomReal*(BY-AY);
                if XJob=0 then
                begin
                    TX := B1;
                    VX := AX+RandomReal*(BX-AX);
                end
                else
                begin
                    VX := (TX-B1)/A1;
                end;
                if YJob=0 then
                begin
                    TY := B2;
                    VY := AY+RandomReal*(BY-AY);
                end
                else
                begin
                    VY := (TY-B2)/A2;
                end;
                V1 := Spline2DCalc(C, TX, TY);
                V2 := Spline2DCalc(C2, VX, VY);
                Err := Max(Err, AbsReal(V1-V2));
                
                //
                // Test F
                //
                Spline2DCopy(C, C2);
                Spline2DLinTransF(C2, A1, B1);
                TX := AX+RandomReal*(BX-AX);
                TY := AY+RandomReal*(BY-AY);
                V1 := Spline2DCalc(C, TX, TY);
                V2 := Spline2DCalc(C2, TX, TY);
                Err := Max(Err, AbsReal(A1*V1+B1-V2));
                Inc(Pass);
            end;
            Inc(YJob);
        end;
        Inc(XJob);
    end;
    Result := AP_FP_Less(Err,10000*MachineEpsilon);
end;


(*************************************************************************
Unset spline, i.e. initialize it with random garbage
*************************************************************************)
procedure UnsetSpline2D(var C : Spline2DInterpolant);
var
    X : TReal1DArray;
    Y : TReal1DArray;
    F : TReal2DArray;
begin
    SetLength(X, 2);
    SetLength(Y, 2);
    SetLength(F, 2, 2);
    X[0] := -1;
    X[1] := +1;
    Y[0] := -1;
    Y[1] := +1;
    F[0,0] := 0;
    F[0,1] := 0;
    F[1,0] := 0;
    F[1,1] := 0;
    Spline2DBuildBilinear(X, Y, F, 2, 2, C);
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testinterpolation2dunit_test_silent():Boolean;
begin
    Result := Test2DInterpolation(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testinterpolation2dunit_test():Boolean;
begin
    Result := Test2DInterpolation(False);
end;


end.