unit testidwunit;
interface
uses Math, Sysutils, Ap, tsort, nearestneighbor, reflections, hblas, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, bdsvd, svd, hqrnd, matgen, trfac, trlinsolve, safesolve, rcond, xblas, densesolver, idwint;

function TestIDW(Silent : Boolean):Boolean;
function testidwunit_test_silent():Boolean;
function testidwunit_test():Boolean;

implementation

procedure Unset2D(var A : TComplex2DArray);forward;
procedure Unset1D(var A : TReal1DArray);forward;
procedure TestXY(const XY : TReal2DArray;
     const N : AlglibInteger;
     const NX : AlglibInteger;
     const D : AlglibInteger;
     const NQ : AlglibInteger;
     const NW : AlglibInteger;
     var IDWErrors : Boolean);forward;
procedure TestRXY(const XY : TReal2DArray;
     const N : AlglibInteger;
     const NX : AlglibInteger;
     const R : Double;
     var IDWErrors : Boolean);forward;
procedure TestDegree(const N : AlglibInteger;
     const NX : AlglibInteger;
     const D : AlglibInteger;
     const DTask : AlglibInteger;
     var IDWErrors : Boolean);forward;
procedure TestNoisy(var IDWErrors : Boolean);forward;


(*************************************************************************
Testing IDW interpolation
*************************************************************************)
function TestIDW(Silent : Boolean):Boolean;
var
    XY : TReal2DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    VX : Double;
    VY : Double;
    VZ : Double;
    D : AlglibInteger;
    DTask : AlglibInteger;
    NX : AlglibInteger;
    N : AlglibInteger;
    NQ : AlglibInteger;
    NW : AlglibInteger;
    SmallN : AlglibInteger;
    LargeN : AlglibInteger;
    WasErrors : Boolean;
    IDWErrors : Boolean;
begin
    IDWErrors := False;
    SmallN := 256;
    LargeN := 1024;
    NQ := 10;
    NW := 18;
    
    //
    // Simple test:
    // * F = x^3 + sin(pi*y)*z^2 - (x+y)^2
    // * space is either R1=[-1,+1] (other dimensions are
    //   fixed at 0), R1^2 or R1^3.
    //* D = -1, 0, 1, 2
    //
    NX:=1;
    while NX<=2 do
    begin
        SetLength(XY, LargeN, NX+1);
        I:=0;
        while I<=LargeN-1 do
        begin
            J:=0;
            while J<=NX-1 do
            begin
                XY[I,J] := 2*RandomReal-1;
                Inc(J);
            end;
            if NX>=1 then
            begin
                VX := XY[I,0];
            end
            else
            begin
                VX := 0;
            end;
            if NX>=2 then
            begin
                VY := XY[I,1];
            end
            else
            begin
                VY := 0;
            end;
            if NX>=3 then
            begin
                VZ := XY[I,2];
            end
            else
            begin
                VZ := 0;
            end;
            XY[I,NX] := VX*VX*VX+Sin(Pi*VY)*AP_Sqr(VZ)-AP_Sqr(VX+VY);
            Inc(I);
        end;
        D:=-1;
        while D<=2 do
        begin
            TestXY(XY, LargeN, NX, D, NQ, NW, IDWErrors);
            Inc(D);
        end;
        Inc(NX);
    end;
    
    //
    // Another simple test:
    // * five points in 2D - (0,0), (0,1), (1,0), (-1,0) (0,-1)
    // * F is random
    // * D = -1, 0, 1, 2
    //
    NX := 2;
    SetLength(XY, 5, NX+1);
    XY[0,0] := 0;
    XY[0,1] := 0;
    XY[0,2] := 2*RandomReal-1;
    XY[1,0] := 1;
    XY[1,1] := 0;
    XY[1,2] := 2*RandomReal-1;
    XY[2,0] := 0;
    XY[2,1] := 1;
    XY[2,2] := 2*RandomReal-1;
    XY[3,0] := -1;
    XY[3,1] := 0;
    XY[3,2] := 2*RandomReal-1;
    XY[4,0] := 0;
    XY[4,1] := -1;
    XY[4,2] := 2*RandomReal-1;
    D:=-1;
    while D<=2 do
    begin
        TestXY(XY, 5, NX, D, NQ, NW, IDWErrors);
        Inc(D);
    end;
    
    //
    // Degree test.
    //
    // F is either:
    // * constant (DTask=0)
    // * linear (DTask=1)
    // * quadratic (DTask=2)
    //
    // Nodal functions are either
    // * constant (D=0)
    // * linear (D=1)
    // * quadratic (D=2)
    //
    // When DTask<=D, we can interpolate without errors.
    // When DTask>D, we MUST have errors.
    //
    NX:=1;
    while NX<=3 do
    begin
        D:=0;
        while D<=2 do
        begin
            DTask:=0;
            while DTask<=2 do
            begin
                TestDegree(SmallN, NX, D, DTask, IDWErrors);
                Inc(DTask);
            end;
            Inc(D);
        end;
        Inc(NX);
    end;
    
    //
    // Noisy test
    //
    TestNoisy(IDWErrors);
    
    //
    // report
    //
    WasErrors := IDWErrors;
    if  not Silent then
    begin
        Write(Format('TESTING INVERSE DISTANCE WEIGHTING'#13#10'',[]));
        Write(Format('* IDW:                                   ',[]));
        if  not IDWErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
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
Unsets 2D array.
*************************************************************************)
procedure Unset2D(var A : TComplex2DArray);
begin
    SetLength(A, 0+1, 0+1);
    A[0,0] := C_Complex(2*RandomReal-1);
end;


(*************************************************************************
Unsets 1D array.
*************************************************************************)
procedure Unset1D(var A : TReal1DArray);
begin
    SetLength(A, 0+1);
    A[0] := 2*RandomReal-1;
end;


(*************************************************************************
Testing IDW:
* generate model using N/NX/D/NQ/NW
* test basic properties
*************************************************************************)
procedure TestXY(const XY : TReal2DArray;
     const N : AlglibInteger;
     const NX : AlglibInteger;
     const D : AlglibInteger;
     const NQ : AlglibInteger;
     const NW : AlglibInteger;
     var IDWErrors : Boolean);
var
    Threshold : Double;
    LipschitzStep : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    V : Double;
    V1 : Double;
    V2 : Double;
    T : Double;
    L1 : Double;
    L2 : Double;
    Z1 : IDWInterpolant;
    X : TReal1DArray;
begin
    Threshold := 1000*MachineEpsilon;
    LipschitzStep := 0.001;
    SetLength(X, NX);
    
    //
    // build
    //
    IDWBuildModifiedShepard(XY, N, NX, D, NQ, NW, Z1);
    
    //
    // first, test interpolation properties at nodes
    //
    I:=0;
    while I<=N-1 do
    begin
        APVMove(@X[0], 0, NX-1, @XY[I][0], 0, NX-1);
        IDWErrors := IDWErrors or AP_FP_Neq(IDWCalc(Z1, X),XY[I,NX]);
        Inc(I);
    end;
    
    //
    // test Lipschitz continuity
    //
    I1 := RandomInteger(N);
    repeat
        I2 := RandomInteger(N);
    until I2<>I1;
    L1 := 0;
    T := 0;
    while AP_FP_Less(T,1) do
    begin
        V := 1-T;
        APVMove(@X[0], 0, NX-1, @XY[I1][0], 0, NX-1, V);
        V := T;
        APVAdd(@X[0], 0, NX-1, @XY[I2][0], 0, NX-1, V);
        V1 := IDWCalc(Z1, X);
        V := 1-(T+LipschitzStep);
        APVMove(@X[0], 0, NX-1, @XY[I1][0], 0, NX-1, V);
        V := T+LipschitzStep;
        APVAdd(@X[0], 0, NX-1, @XY[I2][0], 0, NX-1, V);
        V2 := IDWCalc(Z1, X);
        L1 := Max(L1, AbsReal(V2-V1)/LipschitzStep);
        T := T+LipschitzStep;
    end;
    L2 := 0;
    T := 0;
    while AP_FP_Less(T,1) do
    begin
        V := 1-T;
        APVMove(@X[0], 0, NX-1, @XY[I1][0], 0, NX-1, V);
        V := T;
        APVAdd(@X[0], 0, NX-1, @XY[I2][0], 0, NX-1, V);
        V1 := IDWCalc(Z1, X);
        V := 1-(T+LipschitzStep/3);
        APVMove(@X[0], 0, NX-1, @XY[I1][0], 0, NX-1, V);
        V := T+LipschitzStep/3;
        APVAdd(@X[0], 0, NX-1, @XY[I2][0], 0, NX-1, V);
        V2 := IDWCalc(Z1, X);
        L2 := Max(L2, AbsReal(V2-V1)/(LipschitzStep/3));
        T := T+LipschitzStep/3;
    end;
    IDWErrors := IDWErrors or AP_FP_Greater(L2,2.0*L1);
end;


(*************************************************************************
Testing IDW:
* generate model using R-based model
* test basic properties
*************************************************************************)
procedure TestRXY(const XY : TReal2DArray;
     const N : AlglibInteger;
     const NX : AlglibInteger;
     const R : Double;
     var IDWErrors : Boolean);
var
    Threshold : Double;
    LipschitzStep : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    V : Double;
    V1 : Double;
    V2 : Double;
    T : Double;
    L1 : Double;
    L2 : Double;
    Z1 : IDWInterpolant;
    X : TReal1DArray;
begin
    Threshold := 1000*MachineEpsilon;
    LipschitzStep := 0.001;
    SetLength(X, NX);
    
    //
    // build
    //
    IDWBuildModifiedShepardR(XY, N, NX, R, Z1);
    
    //
    // first, test interpolation properties at nodes
    //
    I:=0;
    while I<=N-1 do
    begin
        APVMove(@X[0], 0, NX-1, @XY[I][0], 0, NX-1);
        IDWErrors := IDWErrors or AP_FP_Neq(IDWCalc(Z1, X),XY[I,NX]);
        Inc(I);
    end;
    
    //
    // test Lipschitz continuity
    //
    I1 := RandomInteger(N);
    repeat
        I2 := RandomInteger(N);
    until I2<>I1;
    L1 := 0;
    T := 0;
    while AP_FP_Less(T,1) do
    begin
        V := 1-T;
        APVMove(@X[0], 0, NX-1, @XY[I1][0], 0, NX-1, V);
        V := T;
        APVAdd(@X[0], 0, NX-1, @XY[I2][0], 0, NX-1, V);
        V1 := IDWCalc(Z1, X);
        V := 1-(T+LipschitzStep);
        APVMove(@X[0], 0, NX-1, @XY[I1][0], 0, NX-1, V);
        V := T+LipschitzStep;
        APVAdd(@X[0], 0, NX-1, @XY[I2][0], 0, NX-1, V);
        V2 := IDWCalc(Z1, X);
        L1 := Max(L1, AbsReal(V2-V1)/LipschitzStep);
        T := T+LipschitzStep;
    end;
    L2 := 0;
    T := 0;
    while AP_FP_Less(T,1) do
    begin
        V := 1-T;
        APVMove(@X[0], 0, NX-1, @XY[I1][0], 0, NX-1, V);
        V := T;
        APVAdd(@X[0], 0, NX-1, @XY[I2][0], 0, NX-1, V);
        V1 := IDWCalc(Z1, X);
        V := 1-(T+LipschitzStep/3);
        APVMove(@X[0], 0, NX-1, @XY[I1][0], 0, NX-1, V);
        V := T+LipschitzStep/3;
        APVAdd(@X[0], 0, NX-1, @XY[I2][0], 0, NX-1, V);
        V2 := IDWCalc(Z1, X);
        L2 := Max(L2, AbsReal(V2-V1)/(LipschitzStep/3));
        T := T+LipschitzStep/3;
    end;
    IDWErrors := IDWErrors or AP_FP_Greater(L2,2.0*L1);
end;


(*************************************************************************
Testing degree properties

F is either:
* constant (DTask=0)
* linear (DTask=1)
* quadratic (DTask=2)

Nodal functions are either
* constant (D=0)
* linear (D=1)
* quadratic (D=2)

When DTask<=D, we can interpolate without errors.
When DTask>D, we MUST have errors.
*************************************************************************)
procedure TestDegree(const N : AlglibInteger;
     const NX : AlglibInteger;
     const D : AlglibInteger;
     const DTask : AlglibInteger;
     var IDWErrors : Boolean);
var
    Threshold : Double;
    NQ : AlglibInteger;
    NW : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    C0 : Double;
    C1 : TReal1DArray;
    C2 : TReal2DArray;
    X : TReal1DArray;
    XY : TReal2DArray;
    Z1 : IDWInterpolant;
    V1 : Double;
    V2 : Double;
begin
    Threshold := 1.0E6*MachineEpsilon;
    NQ := 2*(NX*NX+NX+1);
    NW := 10;
    Assert(NQ<=N, 'TestDegree: internal error');
    
    //
    // prepare model
    //
    C0 := 2*RandomReal-1;
    SetLength(C1, NX);
    I:=0;
    while I<=NX-1 do
    begin
        C1[I] := 2*RandomReal-1;
        Inc(I);
    end;
    SetLength(C2, NX, NX);
    I:=0;
    while I<=NX-1 do
    begin
        J:=I+1;
        while J<=NX-1 do
        begin
            C2[I,J] := 2*RandomReal-1;
            C2[J,I] := C2[I,J];
            Inc(J);
        end;
        repeat
            C2[I,I] := 2*RandomReal-1;
        until AP_FP_Greater(AbsReal(C2[I,I]),0.3);
        Inc(I);
    end;
    
    //
    // prepare points
    //
    SetLength(XY, N, NX+1);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=NX-1 do
        begin
            XY[I,J] := 4*RandomReal-2;
            Inc(J);
        end;
        XY[I,NX] := C0;
        if DTask>=1 then
        begin
            V := APVDotProduct(@C1[0], 0, NX-1, @XY[I][0], 0, NX-1);
            XY[I,NX] := XY[I,NX]+V;
        end;
        if DTask=2 then
        begin
            J:=0;
            while J<=NX-1 do
            begin
                V := APVDotProduct(@C2[J][0], 0, NX-1, @XY[I][0], 0, NX-1);
                XY[I,NX] := XY[I,NX]+XY[I,J]*V;
                Inc(J);
            end;
        end;
        Inc(I);
    end;
    
    //
    // build interpolant, calculate value at random point
    //
    IDWBuildModifiedShepard(XY, N, NX, D, NQ, NW, Z1);
    SetLength(X, NX);
    I:=0;
    while I<=NX-1 do
    begin
        X[I] := 4*RandomReal-2;
        Inc(I);
    end;
    V1 := IDWCalc(Z1, X);
    
    //
    // calculate model value at the same point
    //
    V2 := C0;
    if DTask>=1 then
    begin
        V := APVDotProduct(@C1[0], 0, NX-1, @X[0], 0, NX-1);
        V2 := V2+V;
    end;
    if DTask=2 then
    begin
        J:=0;
        while J<=NX-1 do
        begin
            V := APVDotProduct(@C2[J][0], 0, NX-1, @X[0], 0, NX-1);
            V2 := V2+X[J]*V;
            Inc(J);
        end;
    end;
    
    //
    // Compare
    //
    if DTask<=D then
    begin
        IDWErrors := IDWErrors or AP_FP_Greater(AbsReal(V2-V1),Threshold);
    end
    else
    begin
        IDWErrors := IDWErrors or AP_FP_Less(AbsReal(V2-V1),Threshold);
    end;
end;


(*************************************************************************
Noisy test:
 * F = x^2 + y^2 + z^2 + noise on [-1,+1]^3
 * space is either R1=[-1,+1] (other dimensions are
   fixed at 0), R1^2 or R1^3.
 * D = 1, 2
 * 4096 points is used for function generation,
   4096 points - for testing
 * RMS error of "noisy" model on test set must be
   lower than RMS error of interpolation model.
*************************************************************************)
procedure TestNoisy(var IDWErrors : Boolean);
var
    NoiseLevel : Double;
    NQ : AlglibInteger;
    NW : AlglibInteger;
    D : AlglibInteger;
    NX : AlglibInteger;
    NTrn : AlglibInteger;
    NTst : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    T : Double;
    V1 : Double;
    V2 : Double;
    VE : Double;
    XY : TReal2DArray;
    X : TReal1DArray;
    Z1 : IDWInterpolant;
    Z2 : IDWInterpolant;
    Rms1 : Double;
    Rms2 : Double;
begin
    NQ := 20;
    NW := 40;
    NoiseLevel := 0.2;
    NTrn := 256;
    NTst := 1024;
    D:=1;
    while D<=2 do
    begin
        NX:=1;
        while NX<=2 do
        begin
            
            //
            // prepare dataset
            //
            SetLength(XY, NTrn, NX+1);
            I:=0;
            while I<=NTrn-1 do
            begin
                V := NoiseLevel*(2*RandomReal-1);
                J:=0;
                while J<=NX-1 do
                begin
                    T := 2*RandomReal-1;
                    V := V+AP_Sqr(T);
                    XY[I,J] := T;
                    Inc(J);
                end;
                XY[I,NX] := V;
                Inc(I);
            end;
            
            //
            // build interpolants
            //
            IDWBuildModifiedShepard(XY, NTrn, NX, D, NQ, NW, Z1);
            IDWBuildNoisy(XY, NTrn, NX, D, NQ, NW, Z2);
            
            //
            // calculate RMS errors
            //
            SetLength(X, NX);
            Rms1 := 0;
            Rms2 := 0;
            I:=0;
            while I<=NTst-1 do
            begin
                VE := 0;
                J:=0;
                while J<=NX-1 do
                begin
                    T := 2*RandomReal-1;
                    X[J] := T;
                    VE := VE+AP_Sqr(T);
                    Inc(J);
                end;
                V1 := IDWCalc(Z1, X);
                V2 := IDWCalc(Z2, X);
                Rms1 := Rms1+AP_Sqr(V1-VE);
                Rms2 := Rms2+AP_Sqr(V2-VE);
                Inc(I);
            end;
            IDWErrors := IDWErrors or AP_FP_Greater(Rms2,Rms1);
            Inc(NX);
        end;
        Inc(D);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testidwunit_test_silent():Boolean;
begin
    Result := TestIDW(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testidwunit_test():Boolean;
begin
    Result := TestIDW(False);
end;


end.