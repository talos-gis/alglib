unit testevdunit;
interface
uses Math, Sysutils, Ap, hblas, reflections, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, hsschur, evd;

function TestEVD(Silent : Boolean):Boolean;
function testevdunit_test_silent():Boolean;
function testevdunit_test():Boolean;

implementation

procedure RMatrixFillSparseA(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Sparcity : Double);forward;
procedure CMatrixFillSparseA(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Sparcity : Double);forward;
procedure RMatrixSymmetricSplit(const A : TReal2DArray;
     N : AlglibInteger;
     var AL : TReal2DArray;
     var AU : TReal2DArray);forward;
procedure CMatrixHermitianSplit(const A : TComplex2DArray;
     N : AlglibInteger;
     var AL : TComplex2DArray;
     var AU : TComplex2DArray);forward;
procedure Unset2D(var A : TReal2DArray);forward;
procedure CUnset2D(var A : TComplex2DArray);forward;
procedure Unset1D(var A : TReal1DArray);forward;
procedure CUnset1D(var A : TComplex1DArray);forward;
function TDTestProduct(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     const Z : TReal2DArray;
     const Lambda : TReal1DArray):Double;forward;
function TestProduct(const A : TReal2DArray;
     N : AlglibInteger;
     const Z : TReal2DArray;
     const Lambda : TReal1DArray):Double;forward;
function TestOrt(const Z : TReal2DArray; N : AlglibInteger):Double;forward;
function TestCProduct(const A : TComplex2DArray;
     N : AlglibInteger;
     const Z : TComplex2DArray;
     const Lambda : TReal1DArray):Double;forward;
function TestCOrt(const Z : TComplex2DArray; N : AlglibInteger):Double;forward;
procedure TestSEVDProblem(const A : TReal2DArray;
     const AL : TReal2DArray;
     const AU : TReal2DArray;
     N : AlglibInteger;
     Threshold : Double;
     var SErrors : Boolean;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger);forward;
procedure TestHEVDProblem(const A : TComplex2DArray;
     const AL : TComplex2DArray;
     const AU : TComplex2DArray;
     N : AlglibInteger;
     Threshold : Double;
     var HErrors : Boolean;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger);forward;
procedure TestSEVDBIProblem(const AFull : TReal2DArray;
     const AL : TReal2DArray;
     const AU : TReal2DArray;
     N : AlglibInteger;
     DistVals : Boolean;
     Threshold : Double;
     var SErrors : Boolean;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger);forward;
procedure TestHEVDBIProblem(const AFull : TComplex2DArray;
     const AL : TComplex2DArray;
     const AU : TComplex2DArray;
     N : AlglibInteger;
     DistVals : Boolean;
     Threshold : Double;
     var HErrors : Boolean;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger);forward;
procedure TestTDEVDProblem(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     Threshold : Double;
     var TDErrors : Boolean;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger);forward;
procedure TestTDEVDBIProblem(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     DistVals : Boolean;
     Threshold : Double;
     var SErrors : Boolean;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger);forward;
procedure TestNSEVDProblem(const A : TReal2DArray;
     N : AlglibInteger;
     Threshold : Double;
     var NSErrors : Boolean;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger);forward;
procedure TestEVDSet(const N : AlglibInteger;
     const Threshold : Double;
     const BIThreshold : Double;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger;
     var NSErrors : Boolean;
     var SErrors : Boolean;
     var HErrors : Boolean;
     var TDErrors : Boolean;
     var SBIErrors : Boolean;
     var HBIErrors : Boolean;
     var TDBIErrors : Boolean);forward;


(*************************************************************************
Testing symmetric EVD subroutine
*************************************************************************)
function TestEVD(Silent : Boolean):Boolean;
var
    RA : TReal2DArray;
    N : AlglibInteger;
    J : AlglibInteger;
    FailC : AlglibInteger;
    Runs : AlglibInteger;
    FailR : Double;
    FailThreshold : Double;
    Threshold : Double;
    BIThreshold : Double;
    WasErrors : Boolean;
    NSErrors : Boolean;
    SErrors : Boolean;
    HErrors : Boolean;
    TDErrors : Boolean;
    SBIErrors : Boolean;
    HBIErrors : Boolean;
    TDBIErrors : Boolean;
    WFailed : Boolean;
begin
    FailThreshold := 0.005;
    Threshold := 100000*MachineEpsilon;
    BIThreshold := 1.0E-6;
    NSErrors := False;
    SErrors := False;
    HErrors := False;
    TDErrors := False;
    SBIErrors := False;
    HBIErrors := False;
    TDBIErrors := False;
    FailC := 0;
    Runs := 0;
    
    //
    // Test problems
    //
    N:=1;
    while N<=ABLASBlockSize(RA) do
    begin
        TestEVDSet(N, Threshold, BIThreshold, FailC, Runs, NSErrors, SErrors, HErrors, TDErrors, SBIErrors, HBIErrors, TDBIErrors);
        Inc(N);
    end;
    J:=2;
    while J<=3 do
    begin
        N:=J*ABLASBlockSize(RA)-1;
        while N<=J*ABLASBlockSize(RA)+1 do
        begin
            TestEVDSet(N, Threshold, BIThreshold, FailC, Runs, NSErrors, SErrors, HErrors, TDErrors, SBIErrors, HBIErrors, TDBIErrors);
            Inc(N);
        end;
        Inc(J);
    end;
    
    //
    // report
    //
    WFailed := AP_FP_Greater(AP_Double(FailC)/Runs,FailThreshold);
    WasErrors := NSErrors or SErrors or HErrors or TDErrors or SBIErrors or HBIErrors or TDBIErrors or WFailed;
    if  not Silent then
    begin
        Write(Format('TESTING EVD UNIT'#13#10'',[]));
        Write(Format('NS ERRORS:                               ',[]));
        if  not NSErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('S ERRORS:                                ',[]));
        if  not SErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('H ERRORS:                                ',[]));
        if  not HErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('TD ERRORS:                               ',[]));
        if  not TDErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('SBI ERRORS:                              ',[]));
        if  not SBIErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('HBI ERRORS:                              ',[]));
        if  not HBIErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('TDBI ERRORS:                             ',[]));
        if  not TDBIErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('FAILURE THRESHOLD:                       ',[]));
        if  not WFailed then
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
Sparse fill
*************************************************************************)
procedure RMatrixFillSparseA(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Sparcity : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if AP_FP_Greater_Eq(RandomReal,Sparcity) then
            begin
                A[I,J] := 2*RandomReal-1;
            end
            else
            begin
                A[I,J] := 0;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Sparse fill
*************************************************************************)
procedure CMatrixFillSparseA(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Sparcity : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if AP_FP_Greater_Eq(RandomReal,Sparcity) then
            begin
                A[I,J].X := 2*RandomReal-1;
                A[I,J].Y := 2*RandomReal-1;
            end
            else
            begin
                A[I,J] := C_Complex(0);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Copies A to AL (lower half) and AU (upper half), filling unused parts by
random garbage.
*************************************************************************)
procedure RMatrixSymmetricSplit(const A : TReal2DArray;
     N : AlglibInteger;
     var AL : TReal2DArray;
     var AU : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    I:=0;
    while I<=N-1 do
    begin
        J:=I+1;
        while J<=N-1 do
        begin
            AL[I,J] := 2*RandomReal-1;
            AL[J,I] := A[I,J];
            AU[I,J] := A[I,J];
            AU[J,I] := 2*RandomReal-1;
            Inc(J);
        end;
        AL[I,I] := A[I,I];
        AU[I,I] := A[I,I];
        Inc(I);
    end;
end;


(*************************************************************************
Copies A to AL (lower half) and AU (upper half), filling unused parts by
random garbage.
*************************************************************************)
procedure CMatrixHermitianSplit(const A : TComplex2DArray;
     N : AlglibInteger;
     var AL : TComplex2DArray;
     var AU : TComplex2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    I:=0;
    while I<=N-1 do
    begin
        J:=I+1;
        while J<=N-1 do
        begin
            AL[I,J] := C_Complex(2*RandomReal-1);
            AL[J,I] := Conj(A[I,J]);
            AU[I,J] := A[I,J];
            AU[J,I] := C_Complex(2*RandomReal-1);
            Inc(J);
        end;
        AL[I,I] := A[I,I];
        AU[I,I] := A[I,I];
        Inc(I);
    end;
end;


(*************************************************************************
Unsets 2D array.
*************************************************************************)
procedure Unset2D(var A : TReal2DArray);
begin
    SetLength(A, 0+1, 0+1);
    A[0,0] := 2*RandomReal-1;
end;


(*************************************************************************
Unsets 2D array.
*************************************************************************)
procedure CUnset2D(var A : TComplex2DArray);
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
Unsets 1D array.
*************************************************************************)
procedure CUnset1D(var A : TComplex1DArray);
begin
    SetLength(A, 0+1);
    A[0] := C_Complex(2*RandomReal-1);
end;


(*************************************************************************
Tests Z*Lambda*Z' against tridiag(D,E).
Returns relative error.
*************************************************************************)
function TDTestProduct(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     const Z : TReal2DArray;
     const Lambda : TReal1DArray):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    MX : Double;
begin
    Result := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            
            //
            // Calculate V = A[i,j], A = Z*Lambda*Z'
            //
            V := 0;
            K:=0;
            while K<=N-1 do
            begin
                V := V+Z[I,K]*Lambda[K]*Z[J,K];
                Inc(K);
            end;
            
            //
            // Compare
            //
            if AbsInt(I-J)=0 then
            begin
                Result := Max(Result, AbsReal(V-D[I]));
            end;
            if AbsInt(I-J)=1 then
            begin
                Result := Max(Result, AbsReal(V-E[Min(I, J)]));
            end;
            if AbsInt(I-J)>1 then
            begin
                Result := Max(Result, AbsReal(V));
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    MX := 0;
    I:=0;
    while I<=N-1 do
    begin
        MX := Max(MX, AbsReal(D[I]));
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        MX := Max(MX, AbsReal(E[I]));
        Inc(I);
    end;
    if AP_FP_Eq(MX,0) then
    begin
        MX := 1;
    end;
    Result := Result/MX;
end;


(*************************************************************************
Tests Z*Lambda*Z' against A
Returns relative error.
*************************************************************************)
function TestProduct(const A : TReal2DArray;
     N : AlglibInteger;
     const Z : TReal2DArray;
     const Lambda : TReal1DArray):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    MX : Double;
begin
    Result := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            
            //
            // Calculate V = A[i,j], A = Z*Lambda*Z'
            //
            V := 0;
            K:=0;
            while K<=N-1 do
            begin
                V := V+Z[I,K]*Lambda[K]*Z[J,K];
                Inc(K);
            end;
            
            //
            // Compare
            //
            Result := Max(Result, AbsReal(V-A[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    MX := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            MX := Max(MX, AbsReal(A[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(MX,0) then
    begin
        MX := 1;
    end;
    Result := Result/MX;
end;


(*************************************************************************
Tests Z*Z' against diag(1...1)
Returns absolute error.
*************************************************************************)
function TestOrt(const Z : TReal2DArray; N : AlglibInteger):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
begin
    Result := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + Z[i_,I]*Z[i_,J];
            end;
            if I=J then
            begin
                V := V-1;
            end;
            Result := Max(Result, AbsReal(V));
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Tests Z*Lambda*Z' against A
Returns relative error.
*************************************************************************)
function TestCProduct(const A : TComplex2DArray;
     N : AlglibInteger;
     const Z : TComplex2DArray;
     const Lambda : TReal1DArray):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Complex;
    MX : Double;
begin
    Result := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            
            //
            // Calculate V = A[i,j], A = Z*Lambda*Z'
            //
            V := C_Complex(0);
            K:=0;
            while K<=N-1 do
            begin
                V := C_Add(V,C_Mul(C_MulR(Z[I,K],Lambda[K]),Conj(Z[J,K])));
                Inc(K);
            end;
            
            //
            // Compare
            //
            Result := Max(Result, AbsComplex(C_Sub(V,A[I,J])));
            Inc(J);
        end;
        Inc(I);
    end;
    MX := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            MX := Max(MX, AbsComplex(A[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(MX,0) then
    begin
        MX := 1;
    end;
    Result := Result/MX;
end;


(*************************************************************************
Tests Z*Z' against diag(1...1)
Returns absolute error.
*************************************************************************)
function TestCOrt(const Z : TComplex2DArray; N : AlglibInteger):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Complex;
    i_ : AlglibInteger;
begin
    Result := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                V := C_Add(V,C_Mul(Z[i_,I],Conj(Z[i_,J])));
            end;
            if I=J then
            begin
                V := C_SubR(V,1);
            end;
            Result := Max(Result, AbsComplex(V));
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Tests SEVD problem
*************************************************************************)
procedure TestSEVDProblem(const A : TReal2DArray;
     const AL : TReal2DArray;
     const AU : TReal2DArray;
     N : AlglibInteger;
     Threshold : Double;
     var SErrors : Boolean;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger);
var
    Lambda : TReal1DArray;
    LambdaRef : TReal1DArray;
    Z : TReal2DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
begin
    
    //
    // Test simple EVD: values and full vectors, lower A
    //
    Unset1D(LambdaRef);
    Unset2D(Z);
    Runs := Runs+1;
    if  not SMatrixEVD(AL, N, 1, False, LambdaRef, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    SErrors := SErrors or AP_FP_Greater(TestProduct(A, N, Z, LambdaRef),Threshold);
    SErrors := SErrors or AP_FP_Greater(TestOrt(Z, N),Threshold);
    I:=0;
    while I<=N-2 do
    begin
        if AP_FP_Less(LambdaRef[I+1],LambdaRef[I]) then
        begin
            SErrors := True;
            Exit;
        end;
        Inc(I);
    end;
    
    //
    // Test simple EVD: values and full vectors, upper A
    //
    Unset1D(Lambda);
    Unset2D(Z);
    Runs := Runs+1;
    if  not SMatrixEVD(AU, N, 1, True, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    SErrors := SErrors or AP_FP_Greater(TestProduct(A, N, Z, Lambda),Threshold);
    SErrors := SErrors or AP_FP_Greater(TestOrt(Z, N),Threshold);
    I:=0;
    while I<=N-2 do
    begin
        if AP_FP_Less(Lambda[I+1],Lambda[I]) then
        begin
            SErrors := True;
            Exit;
        end;
        Inc(I);
    end;
    
    //
    // Test simple EVD: values only, lower A
    //
    Unset1D(Lambda);
    Unset2D(Z);
    Runs := Runs+1;
    if  not SMatrixEVD(AL, N, 0, False, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[I]-LambdaRef[I]),Threshold);
        Inc(I);
    end;
    
    //
    // Test simple EVD: values only, upper A
    //
    Unset1D(Lambda);
    Unset2D(Z);
    Runs := Runs+1;
    if  not SMatrixEVD(AU, N, 0, True, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[I]-LambdaRef[I]),Threshold);
        Inc(I);
    end;
end;


(*************************************************************************
Tests SEVD problem
*************************************************************************)
procedure TestHEVDProblem(const A : TComplex2DArray;
     const AL : TComplex2DArray;
     const AU : TComplex2DArray;
     N : AlglibInteger;
     Threshold : Double;
     var HErrors : Boolean;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger);
var
    Lambda : TReal1DArray;
    LambdaRef : TReal1DArray;
    Z : TComplex2DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Complex;
begin
    
    //
    // Test simple EVD: values and full vectors, lower A
    //
    Unset1D(LambdaRef);
    CUnset2D(Z);
    Runs := Runs+1;
    if  not HMatrixEVD(AL, N, 1, False, LambdaRef, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    HErrors := HErrors or AP_FP_Greater(TestCProduct(A, N, Z, LambdaRef),Threshold);
    HErrors := HErrors or AP_FP_Greater(TestCOrt(Z, N),Threshold);
    I:=0;
    while I<=N-2 do
    begin
        if AP_FP_Less(LambdaRef[I+1],LambdaRef[I]) then
        begin
            HErrors := True;
            Exit;
        end;
        Inc(I);
    end;
    
    //
    // Test simple EVD: values and full vectors, upper A
    //
    Unset1D(Lambda);
    CUnset2D(Z);
    Runs := Runs+1;
    if  not HMatrixEVD(AU, N, 1, True, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    HErrors := HErrors or AP_FP_Greater(TestCProduct(A, N, Z, Lambda),Threshold);
    HErrors := HErrors or AP_FP_Greater(TestCOrt(Z, N),Threshold);
    I:=0;
    while I<=N-2 do
    begin
        if AP_FP_Less(Lambda[I+1],Lambda[I]) then
        begin
            HErrors := True;
            Exit;
        end;
        Inc(I);
    end;
    
    //
    // Test simple EVD: values only, lower A
    //
    Unset1D(Lambda);
    CUnset2D(Z);
    Runs := Runs+1;
    if  not HMatrixEVD(AL, N, 0, False, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        HErrors := HErrors or AP_FP_Greater(AbsReal(Lambda[I]-LambdaRef[I]),Threshold);
        Inc(I);
    end;
    
    //
    // Test simple EVD: values only, upper A
    //
    Unset1D(Lambda);
    CUnset2D(Z);
    Runs := Runs+1;
    if  not HMatrixEVD(AU, N, 0, True, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        HErrors := HErrors or AP_FP_Greater(AbsReal(Lambda[I]-LambdaRef[I]),Threshold);
        Inc(I);
    end;
end;


(*************************************************************************
Tests EVD problem

DistVals    -   is True, when eigenvalues are distinct. Is False, when we
                are solving sparse task with  lots  of  zero  eigenvalues.
                In such cases some tests related to the  eigenvectors  are
                not performed.
*************************************************************************)
procedure TestSEVDBIProblem(const AFull : TReal2DArray;
     const AL : TReal2DArray;
     const AU : TReal2DArray;
     N : AlglibInteger;
     DistVals : Boolean;
     Threshold : Double;
     var SErrors : Boolean;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger);
var
    Lambda : TReal1DArray;
    LambdaRef : TReal1DArray;
    Z : TReal2DArray;
    ZRef : TReal2DArray;
    A1 : TReal2DArray;
    A2 : TReal2DArray;
    AR : TReal2DArray;
    WSucc : Boolean;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    M : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    V : Double;
    A : Double;
    B : Double;
    i_ : AlglibInteger;
begin
    SetLength(LambdaRef, N-1+1);
    SetLength(ZRef, N-1+1, N-1+1);
    SetLength(A1, N-1+1, N-1+1);
    SetLength(A2, N-1+1, N-1+1);
    
    //
    // Reference EVD
    //
    Runs := Runs+1;
    if  not SMatrixEVD(AFull, N, 1, True, LambdaRef, ZRef) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    
    //
    // Select random interval boundaries.
    // If there are non-distinct eigenvalues at the boundaries,
    // we move indexes further until values splits. It is done to
    // avoid situations where we can't get definite answer.
    //
    I1 := RandomInteger(N);
    I2 := I1+RandomInteger(N-I1);
    while I1>0 do
    begin
        if AP_FP_Greater(AbsReal(LambdaRef[I1-1]-LambdaRef[I1]),10*Threshold) then
        begin
            Break;
        end;
        I1 := I1-1;
    end;
    while I2<N-1 do
    begin
        if AP_FP_Greater(AbsReal(LambdaRef[I2+1]-LambdaRef[I2]),10*Threshold) then
        begin
            Break;
        end;
        I2 := I2+1;
    end;
    
    //
    // Select A, B
    //
    if I1>0 then
    begin
        A := 0.5*(LambdaRef[I1]+LambdaRef[I1-1]);
    end
    else
    begin
        A := LambdaRef[0]-1;
    end;
    if I2<N-1 then
    begin
        B := 0.5*(LambdaRef[I2]+LambdaRef[I2+1]);
    end
    else
    begin
        B := LambdaRef[N-1]+1;
    end;
    
    //
    // Test interval, no vectors, lower A
    //
    Unset1D(Lambda);
    Unset2D(Z);
    Runs := Runs+1;
    if  not SMatrixEVDR(AL, N, 0, False, A, B, M, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    if M<>I2-I1+1 then
    begin
        FailC := FailC+1;
        Exit;
    end;
    K:=0;
    while K<=M-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    
    //
    // Test interval, no vectors, upper A
    //
    Unset1D(Lambda);
    Unset2D(Z);
    Runs := Runs+1;
    if  not SMatrixEVDR(AU, N, 0, True, A, B, M, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    if M<>I2-I1+1 then
    begin
        FailC := FailC+1;
        Exit;
    end;
    K:=0;
    while K<=M-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    
    //
    // Test indexes, no vectors, lower A
    //
    Unset1D(Lambda);
    Unset2D(Z);
    Runs := Runs+1;
    if  not SMatrixEVDI(AL, N, 0, False, I1, I2, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    M := I2-I1+1;
    K:=0;
    while K<=M-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    
    //
    // Test indexes, no vectors, upper A
    //
    Unset1D(Lambda);
    Unset2D(Z);
    Runs := Runs+1;
    if  not SMatrixEVDI(AU, N, 0, True, I1, I2, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    M := I2-I1+1;
    K:=0;
    while K<=M-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    
    //
    // Test interval, vectors, lower A
    //
    Unset1D(Lambda);
    Unset2D(Z);
    Runs := Runs+1;
    if  not SMatrixEVDR(AL, N, 1, False, A, B, M, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    if M<>I2-I1+1 then
    begin
        FailC := FailC+1;
        Exit;
    end;
    K:=0;
    while K<=M-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    if DistVals then
    begin
        
        //
        // Distinct eigenvalues, test vectors
        //
        J:=0;
        while J<=M-1 do
        begin
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + Z[i_,J]*ZRef[i_,I1+J];
            end;
            if AP_FP_Less(V,0) then
            begin
                for i_ := 0 to N-1 do
                begin
                    Z[i_,J] := -1*Z[i_,J];
                end;
            end;
            Inc(J);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                SErrors := SErrors or AP_FP_Greater(AbsReal(Z[I,J]-ZRef[I,I1+J]),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    
    //
    // Test interval, vectors, upper A
    //
    Unset1D(Lambda);
    Unset2D(Z);
    Runs := Runs+1;
    if  not SMatrixEVDR(AU, N, 1, True, A, B, M, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    if M<>I2-I1+1 then
    begin
        FailC := FailC+1;
        Exit;
    end;
    K:=0;
    while K<=M-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    if DistVals then
    begin
        
        //
        // Distinct eigenvalues, test vectors
        //
        J:=0;
        while J<=M-1 do
        begin
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + Z[i_,J]*ZRef[i_,I1+J];
            end;
            if AP_FP_Less(V,0) then
            begin
                for i_ := 0 to N-1 do
                begin
                    Z[i_,J] := -1*Z[i_,J];
                end;
            end;
            Inc(J);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                SErrors := SErrors or AP_FP_Greater(AbsReal(Z[I,J]-ZRef[I,I1+J]),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    
    //
    // Test indexes, vectors, lower A
    //
    Unset1D(Lambda);
    Unset2D(Z);
    Runs := Runs+1;
    if  not SMatrixEVDI(AL, N, 1, False, I1, I2, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    M := I2-I1+1;
    K:=0;
    while K<=M-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    if DistVals then
    begin
        
        //
        // Distinct eigenvalues, test vectors
        //
        J:=0;
        while J<=M-1 do
        begin
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + Z[i_,J]*ZRef[i_,I1+J];
            end;
            if AP_FP_Less(V,0) then
            begin
                for i_ := 0 to N-1 do
                begin
                    Z[i_,J] := -1*Z[i_,J];
                end;
            end;
            Inc(J);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                SErrors := SErrors or AP_FP_Greater(AbsReal(Z[I,J]-ZRef[I,I1+J]),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    
    //
    // Test indexes, vectors, upper A
    //
    Unset1D(Lambda);
    Unset2D(Z);
    Runs := Runs+1;
    if  not SMatrixEVDI(AU, N, 1, True, I1, I2, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    M := I2-I1+1;
    K:=0;
    while K<=M-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    if DistVals then
    begin
        
        //
        // Distinct eigenvalues, test vectors
        //
        J:=0;
        while J<=M-1 do
        begin
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + Z[i_,J]*ZRef[i_,I1+J];
            end;
            if AP_FP_Less(V,0) then
            begin
                for i_ := 0 to N-1 do
                begin
                    Z[i_,J] := -1*Z[i_,J];
                end;
            end;
            Inc(J);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                SErrors := SErrors or AP_FP_Greater(AbsReal(Z[I,J]-ZRef[I,I1+J]),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Tests EVD problem

DistVals    -   is True, when eigenvalues are distinct. Is False, when we
                are solving sparse task with  lots  of  zero  eigenvalues.
                In such cases some tests related to the  eigenvectors  are
                not performed.
*************************************************************************)
procedure TestHEVDBIProblem(const AFull : TComplex2DArray;
     const AL : TComplex2DArray;
     const AU : TComplex2DArray;
     N : AlglibInteger;
     DistVals : Boolean;
     Threshold : Double;
     var HErrors : Boolean;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger);
var
    Lambda : TReal1DArray;
    LambdaRef : TReal1DArray;
    Z : TComplex2DArray;
    ZRef : TComplex2DArray;
    A1 : TComplex2DArray;
    A2 : TComplex2DArray;
    AR : TComplex2DArray;
    WSucc : Boolean;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    M : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    V : Complex;
    A : Double;
    B : Double;
    i_ : AlglibInteger;
begin
    SetLength(LambdaRef, N-1+1);
    SetLength(ZRef, N-1+1, N-1+1);
    SetLength(A1, N-1+1, N-1+1);
    SetLength(A2, N-1+1, N-1+1);
    
    //
    // Reference EVD
    //
    Runs := Runs+1;
    if  not HMatrixEVD(AFull, N, 1, True, LambdaRef, ZRef) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    
    //
    // Select random interval boundaries.
    // If there are non-distinct eigenvalues at the boundaries,
    // we move indexes further until values splits. It is done to
    // avoid situations where we can't get definite answer.
    //
    I1 := RandomInteger(N);
    I2 := I1+RandomInteger(N-I1);
    while I1>0 do
    begin
        if AP_FP_Greater(AbsReal(LambdaRef[I1-1]-LambdaRef[I1]),10*Threshold) then
        begin
            Break;
        end;
        I1 := I1-1;
    end;
    while I2<N-1 do
    begin
        if AP_FP_Greater(AbsReal(LambdaRef[I2+1]-LambdaRef[I2]),10*Threshold) then
        begin
            Break;
        end;
        I2 := I2+1;
    end;
    
    //
    // Select A, B
    //
    if I1>0 then
    begin
        A := 0.5*(LambdaRef[I1]+LambdaRef[I1-1]);
    end
    else
    begin
        A := LambdaRef[0]-1;
    end;
    if I2<N-1 then
    begin
        B := 0.5*(LambdaRef[I2]+LambdaRef[I2+1]);
    end
    else
    begin
        B := LambdaRef[N-1]+1;
    end;
    
    //
    // Test interval, no vectors, lower A
    //
    Unset1D(Lambda);
    CUnset2D(Z);
    Runs := Runs+1;
    if  not HMatrixEVDR(AL, N, 0, False, A, B, M, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    if M<>I2-I1+1 then
    begin
        FailC := FailC+1;
        Exit;
    end;
    K:=0;
    while K<=M-1 do
    begin
        HErrors := HErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    
    //
    // Test interval, no vectors, upper A
    //
    Unset1D(Lambda);
    CUnset2D(Z);
    Runs := Runs+1;
    if  not HMatrixEVDR(AU, N, 0, True, A, B, M, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    if M<>I2-I1+1 then
    begin
        FailC := FailC+1;
        Exit;
    end;
    K:=0;
    while K<=M-1 do
    begin
        HErrors := HErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    
    //
    // Test indexes, no vectors, lower A
    //
    Unset1D(Lambda);
    CUnset2D(Z);
    Runs := Runs+1;
    if  not HMatrixEVDI(AL, N, 0, False, I1, I2, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    M := I2-I1+1;
    K:=0;
    while K<=M-1 do
    begin
        HErrors := HErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    
    //
    // Test indexes, no vectors, upper A
    //
    Unset1D(Lambda);
    CUnset2D(Z);
    Runs := Runs+1;
    if  not HMatrixEVDI(AU, N, 0, True, I1, I2, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    M := I2-I1+1;
    K:=0;
    while K<=M-1 do
    begin
        HErrors := HErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    
    //
    // Test interval, vectors, lower A
    //
    Unset1D(Lambda);
    CUnset2D(Z);
    Runs := Runs+1;
    if  not HMatrixEVDR(AL, N, 1, False, A, B, M, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    if M<>I2-I1+1 then
    begin
        FailC := FailC+1;
        Exit;
    end;
    K:=0;
    while K<=M-1 do
    begin
        HErrors := HErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    if DistVals then
    begin
        
        //
        // Distinct eigenvalues, test vectors
        //
        J:=0;
        while J<=M-1 do
        begin
            V := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                V := C_Add(V,C_Mul(Z[i_,J],Conj(ZRef[i_,I1+J])));
            end;
            V := Conj(C_DivR(V,AbsComplex(V)));
            for i_ := 0 to N-1 do
            begin
                Z[i_,J] := C_Mul(V, Z[i_,J]);
            end;
            Inc(J);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                HErrors := HErrors or AP_FP_Greater(AbsComplex(C_Sub(Z[I,J],ZRef[I,I1+J])),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    
    //
    // Test interval, vectors, upper A
    //
    Unset1D(Lambda);
    CUnset2D(Z);
    Runs := Runs+1;
    if  not HMatrixEVDR(AU, N, 1, True, A, B, M, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    if M<>I2-I1+1 then
    begin
        FailC := FailC+1;
        Exit;
    end;
    K:=0;
    while K<=M-1 do
    begin
        HErrors := HErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    if DistVals then
    begin
        
        //
        // Distinct eigenvalues, test vectors
        //
        J:=0;
        while J<=M-1 do
        begin
            V := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                V := C_Add(V,C_Mul(Z[i_,J],Conj(ZRef[i_,I1+J])));
            end;
            V := Conj(C_DivR(V,AbsComplex(V)));
            for i_ := 0 to N-1 do
            begin
                Z[i_,J] := C_Mul(V, Z[i_,J]);
            end;
            Inc(J);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                HErrors := HErrors or AP_FP_Greater(AbsComplex(C_Sub(Z[I,J],ZRef[I,I1+J])),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    
    //
    // Test indexes, vectors, lower A
    //
    Unset1D(Lambda);
    CUnset2D(Z);
    Runs := Runs+1;
    if  not HMatrixEVDI(AL, N, 1, False, I1, I2, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    M := I2-I1+1;
    K:=0;
    while K<=M-1 do
    begin
        HErrors := HErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    if DistVals then
    begin
        
        //
        // Distinct eigenvalues, test vectors
        //
        J:=0;
        while J<=M-1 do
        begin
            V := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                V := C_Add(V,C_Mul(Z[i_,J],Conj(ZRef[i_,I1+J])));
            end;
            V := Conj(C_DivR(V,AbsComplex(V)));
            for i_ := 0 to N-1 do
            begin
                Z[i_,J] := C_Mul(V, Z[i_,J]);
            end;
            Inc(J);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                HErrors := HErrors or AP_FP_Greater(AbsComplex(C_Sub(Z[I,J],ZRef[I,I1+J])),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    
    //
    // Test indexes, vectors, upper A
    //
    Unset1D(Lambda);
    CUnset2D(Z);
    Runs := Runs+1;
    if  not HMatrixEVDI(AU, N, 1, True, I1, I2, Lambda, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    M := I2-I1+1;
    K:=0;
    while K<=M-1 do
    begin
        HErrors := HErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    if DistVals then
    begin
        
        //
        // Distinct eigenvalues, test vectors
        //
        J:=0;
        while J<=M-1 do
        begin
            V := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                V := C_Add(V,C_Mul(Z[i_,J],Conj(ZRef[i_,I1+J])));
            end;
            V := Conj(C_DivR(V,AbsComplex(V)));
            for i_ := 0 to N-1 do
            begin
                Z[i_,J] := C_Mul(V, Z[i_,J]);
            end;
            Inc(J);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                HErrors := HErrors or AP_FP_Greater(AbsComplex(C_Sub(Z[I,J],ZRef[I,I1+J])),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Tests EVD problem
*************************************************************************)
procedure TestTDEVDProblem(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     Threshold : Double;
     var TDErrors : Boolean;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger);
var
    Lambda : TReal1DArray;
    EE : TReal1DArray;
    Lambda2 : TReal1DArray;
    Z : TReal2DArray;
    ZRef : TReal2DArray;
    A1 : TReal2DArray;
    A2 : TReal2DArray;
    WSucc : Boolean;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
begin
    SetLength(Lambda, N-1+1);
    SetLength(Lambda2, N-1+1);
    SetLength(ZRef, N-1+1, N-1+1);
    SetLength(A1, N-1+1, N-1+1);
    SetLength(A2, N-1+1, N-1+1);
    if N>1 then
    begin
        SetLength(EE, N-2+1);
    end;
    
    //
    // Test simple EVD: values and full vectors
    //
    I:=0;
    while I<=N-1 do
    begin
        Lambda[I] := D[I];
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        EE[I] := E[I];
        Inc(I);
    end;
    Runs := Runs+1;
    WSucc := SMatrixTDEVD(Lambda, EE, N, 2, Z);
    if  not WSucc then
    begin
        FailC := FailC+1;
        Exit;
    end;
    TDErrors := TDErrors or AP_FP_Greater(TDTestProduct(D, E, N, Z, Lambda),Threshold);
    TDErrors := TDErrors or AP_FP_Greater(TestOrt(Z, N),Threshold);
    I:=0;
    while I<=N-2 do
    begin
        if AP_FP_Less(Lambda[I+1],Lambda[I]) then
        begin
            TDErrors := True;
            Exit;
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            ZRef[I,J] := Z[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Test values only variant
    //
    I:=0;
    while I<=N-1 do
    begin
        Lambda2[I] := D[I];
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        EE[I] := E[I];
        Inc(I);
    end;
    Runs := Runs+1;
    WSucc := SMatrixTDEVD(Lambda2, EE, N, 0, Z);
    if  not WSucc then
    begin
        FailC := FailC+1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        TDErrors := TDErrors or AP_FP_Greater(AbsReal(Lambda2[I]-Lambda[I]),Threshold);
        Inc(I);
    end;
    
    //
    // Test multiplication variant
    //
    I:=0;
    while I<=N-1 do
    begin
        Lambda2[I] := D[I];
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        EE[I] := E[I];
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            A1[I,J] := 2*RandomReal-1;
            A2[I,J] := A1[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
    Runs := Runs+1;
    WSucc := SMatrixTDEVD(Lambda2, EE, N, 1, A1);
    if  not WSucc then
    begin
        FailC := FailC+1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        TDErrors := TDErrors or AP_FP_Greater(AbsReal(Lambda2[I]-Lambda[I]),Threshold);
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + A2[I,i_]*ZRef[i_,J];
            end;
            
            //
            // next line is a bit complicated because
            // depending on algorithm used we can get either
            // z or -z as eigenvector. so we compare result
            // with both A*ZRef and -A*ZRef
            //
            TDErrors := TDErrors or AP_FP_Greater(AbsReal(V-A1[I,J]),Threshold) and AP_FP_Greater(AbsReal(V+A1[I,J]),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Test first row variant
    //
    I:=0;
    while I<=N-1 do
    begin
        Lambda2[I] := D[I];
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        EE[I] := E[I];
        Inc(I);
    end;
    Runs := Runs+1;
    WSucc := SMatrixTDEVD(Lambda2, EE, N, 3, Z);
    if  not WSucc then
    begin
        FailC := FailC+1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        TDErrors := TDErrors or AP_FP_Greater(AbsReal(Lambda2[I]-Lambda[I]),Threshold);
        
        //
        // next line is a bit complicated because
        // depending on algorithm used we can get either
        // z or -z as eigenvector. so we compare result
        // with both z and -z
        //
        TDErrors := TDErrors or AP_FP_Greater(AbsReal(Z[0,I]-ZRef[0,I]),Threshold) and AP_FP_Greater(AbsReal(Z[0,I]+ZRef[0,I]),Threshold);
        Inc(I);
    end;
end;


(*************************************************************************
Tests EVD problem

DistVals    -   is True, when eigenvalues are distinct. Is False, when we
                are solving sparse task with  lots  of  zero  eigenvalues.
                In such cases some tests related to the  eigenvectors  are
                not performed.
*************************************************************************)
procedure TestTDEVDBIProblem(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     DistVals : Boolean;
     Threshold : Double;
     var SErrors : Boolean;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger);
var
    Lambda : TReal1DArray;
    LambdaRef : TReal1DArray;
    Z : TReal2DArray;
    ZRef : TReal2DArray;
    A1 : TReal2DArray;
    A2 : TReal2DArray;
    AR : TReal2DArray;
    WSucc : Boolean;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    M : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    V : Double;
    A : Double;
    B : Double;
    i_ : AlglibInteger;
begin
    SetLength(LambdaRef, N-1+1);
    SetLength(ZRef, N-1+1, N-1+1);
    SetLength(A1, N-1+1, N-1+1);
    SetLength(A2, N-1+1, N-1+1);
    
    //
    // Reference EVD
    //
    SetLength(LambdaRef, N);
    APVMove(@LambdaRef[0], 0, N-1, @D[0], 0, N-1);
    Runs := Runs+1;
    if  not SMatrixTDEVD(LambdaRef, E, N, 2, ZRef) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    
    //
    // Select random interval boundaries.
    // If there are non-distinct eigenvalues at the boundaries,
    // we move indexes further until values splits. It is done to
    // avoid situations where we can't get definite answer.
    //
    I1 := RandomInteger(N);
    I2 := I1+RandomInteger(N-I1);
    while I1>0 do
    begin
        if AP_FP_Greater(AbsReal(LambdaRef[I1-1]-LambdaRef[I1]),10*Threshold) then
        begin
            Break;
        end;
        I1 := I1-1;
    end;
    while I2<N-1 do
    begin
        if AP_FP_Greater(AbsReal(LambdaRef[I2+1]-LambdaRef[I2]),10*Threshold) then
        begin
            Break;
        end;
        I2 := I2+1;
    end;
    
    //
    // Test different combinations
    //
    
    //
    // Select A, B
    //
    if I1>0 then
    begin
        A := 0.5*(LambdaRef[I1]+LambdaRef[I1-1]);
    end
    else
    begin
        A := LambdaRef[0]-1;
    end;
    if I2<N-1 then
    begin
        B := 0.5*(LambdaRef[I2]+LambdaRef[I2+1]);
    end
    else
    begin
        B := LambdaRef[N-1]+1;
    end;
    
    //
    // Test interval, no vectors
    //
    SetLength(Lambda, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        Lambda[I] := D[I];
        Inc(I);
    end;
    Runs := Runs+1;
    if  not SMatrixTDEVDR(Lambda, E, N, 0, A, B, M, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    if M<>I2-I1+1 then
    begin
        FailC := FailC+1;
        Exit;
    end;
    K:=0;
    while K<=M-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    
    //
    // Test indexes, no vectors
    //
    SetLength(Lambda, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        Lambda[I] := D[I];
        Inc(I);
    end;
    Runs := Runs+1;
    if  not SMatrixTDEVDI(Lambda, E, N, 0, I1, I2, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    M := I2-I1+1;
    K:=0;
    while K<=M-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    
    //
    // Test interval, transform vectors
    //
    SetLength(Lambda, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        Lambda[I] := D[I];
        Inc(I);
    end;
    SetLength(A1, N-1+1, N-1+1);
    SetLength(A2, N-1+1, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            A1[I,J] := 2*RandomReal-1;
            A2[I,J] := A1[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
    Runs := Runs+1;
    if  not SMatrixTDEVDR(Lambda, E, N, 1, A, B, M, A1) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    if M<>I2-I1+1 then
    begin
        FailC := FailC+1;
        Exit;
    end;
    K:=0;
    while K<=M-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    if DistVals then
    begin
        SetLength(AR, N-1+1, M-1+1);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                V := 0.0;
                for i_ := 0 to N-1 do
                begin
                    V := V + A2[I,i_]*ZRef[i_,I1+J];
                end;
                AR[I,J] := V;
                Inc(J);
            end;
            Inc(I);
        end;
        J:=0;
        while J<=M-1 do
        begin
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + A1[i_,J]*AR[i_,J];
            end;
            if AP_FP_Less(V,0) then
            begin
                for i_ := 0 to N-1 do
                begin
                    AR[i_,J] := -1*AR[i_,J];
                end;
            end;
            Inc(J);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                SErrors := SErrors or AP_FP_Greater(AbsReal(A1[I,J]-AR[I,J]),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    
    //
    // Test indexes, transform vectors
    //
    SetLength(Lambda, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        Lambda[I] := D[I];
        Inc(I);
    end;
    SetLength(A1, N-1+1, N-1+1);
    SetLength(A2, N-1+1, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            A1[I,J] := 2*RandomReal-1;
            A2[I,J] := A1[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
    Runs := Runs+1;
    if  not SMatrixTDEVDI(Lambda, E, N, 1, I1, I2, A1) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    M := I2-I1+1;
    K:=0;
    while K<=M-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    if DistVals then
    begin
        SetLength(AR, N-1+1, M-1+1);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                V := 0.0;
                for i_ := 0 to N-1 do
                begin
                    V := V + A2[I,i_]*ZRef[i_,I1+J];
                end;
                AR[I,J] := V;
                Inc(J);
            end;
            Inc(I);
        end;
        J:=0;
        while J<=M-1 do
        begin
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + A1[i_,J]*AR[i_,J];
            end;
            if AP_FP_Less(V,0) then
            begin
                for i_ := 0 to N-1 do
                begin
                    AR[i_,J] := -1*AR[i_,J];
                end;
            end;
            Inc(J);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                SErrors := SErrors or AP_FP_Greater(AbsReal(A1[I,J]-AR[I,J]),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    
    //
    // Test interval, do not transform vectors
    //
    SetLength(Lambda, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        Lambda[I] := D[I];
        Inc(I);
    end;
    SetLength(Z, 0+1, 0+1);
    Runs := Runs+1;
    if  not SMatrixTDEVDR(Lambda, E, N, 2, A, B, M, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    if M<>I2-I1+1 then
    begin
        FailC := FailC+1;
        Exit;
    end;
    K:=0;
    while K<=M-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    if DistVals then
    begin
        J:=0;
        while J<=M-1 do
        begin
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + Z[i_,J]*ZRef[i_,I1+J];
            end;
            if AP_FP_Less(V,0) then
            begin
                for i_ := 0 to N-1 do
                begin
                    Z[i_,J] := -1*Z[i_,J];
                end;
            end;
            Inc(J);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                SErrors := SErrors or AP_FP_Greater(AbsReal(Z[I,J]-ZRef[I,I1+J]),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    
    //
    // Test indexes, do not transform vectors
    //
    SetLength(Lambda, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        Lambda[I] := D[I];
        Inc(I);
    end;
    SetLength(Z, 0+1, 0+1);
    Runs := Runs+1;
    if  not SMatrixTDEVDI(Lambda, E, N, 2, I1, I2, Z) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    M := I2-I1+1;
    K:=0;
    while K<=M-1 do
    begin
        SErrors := SErrors or AP_FP_Greater(AbsReal(Lambda[K]-LambdaRef[I1+K]),Threshold);
        Inc(K);
    end;
    if DistVals then
    begin
        J:=0;
        while J<=M-1 do
        begin
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + Z[i_,J]*ZRef[i_,I1+J];
            end;
            if AP_FP_Less(V,0) then
            begin
                for i_ := 0 to N-1 do
                begin
                    Z[i_,J] := -1*Z[i_,J];
                end;
            end;
            Inc(J);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                SErrors := SErrors or AP_FP_Greater(AbsReal(Z[I,J]-ZRef[I,I1+J]),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Non-symmetric problem
*************************************************************************)
procedure TestNSEVDProblem(const A : TReal2DArray;
     N : AlglibInteger;
     Threshold : Double;
     var NSErrors : Boolean;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger);
var
    MX : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    VJob : AlglibInteger;
    NeedL : Boolean;
    NeedR : Boolean;
    WR0 : TReal1DArray;
    WI0 : TReal1DArray;
    WR1 : TReal1DArray;
    WI1 : TReal1DArray;
    WR0S : TReal1DArray;
    WI0S : TReal1DArray;
    WR1S : TReal1DArray;
    WI1S : TReal1DArray;
    VL : TReal2DArray;
    VR : TReal2DArray;
    Vec1R : TReal1DArray;
    Vec1I : TReal1DArray;
    Vec2R : TReal1DArray;
    Vec2I : TReal1DArray;
    Vec3R : TReal1DArray;
    Vec3I : TReal1DArray;
    CurWR : Double;
    CurWI : Double;
    VT : Double;
    Tmp : Double;
    i_ : AlglibInteger;
begin
    SetLength(Vec1R, N-1+1);
    SetLength(Vec2R, N-1+1);
    SetLength(Vec3R, N-1+1);
    SetLength(Vec1I, N-1+1);
    SetLength(Vec2I, N-1+1);
    SetLength(Vec3I, N-1+1);
    SetLength(WR0S, N-1+1);
    SetLength(WR1S, N-1+1);
    SetLength(WI0S, N-1+1);
    SetLength(WI1S, N-1+1);
    MX := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if AP_FP_Greater(AbsReal(A[I,J]),MX) then
            begin
                MX := AbsReal(A[I,J]);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(MX,0) then
    begin
        MX := 1;
    end;
    
    //
    // Load values-only
    //
    Runs := Runs+1;
    if  not RMatrixEVD(A, N, 0, WR0, WI0, VL, VR) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    
    //
    // Test different jobs
    //
    VJob:=1;
    while VJob<=3 do
    begin
        NeedR := (VJob=1) or (VJob=3);
        NeedL := (VJob=2) or (VJob=3);
        Runs := Runs+1;
        if  not RMatrixEVD(A, N, VJob, WR1, WI1, VL, VR) then
        begin
            FailC := FailC+1;
            Exit;
        end;
        
        //
        // Test values:
        // 1. sort by real part
        // 2. test
        //
        APVMove(@WR0S[0], 0, N-1, @WR0[0], 0, N-1);
        APVMove(@WI0S[0], 0, N-1, @WI0[0], 0, N-1);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-2-I do
            begin
                if AP_FP_Greater(WR0S[J],WR0S[J+1]) then
                begin
                    Tmp := WR0S[J];
                    WR0S[J] := WR0S[J+1];
                    WR0S[J+1] := Tmp;
                    Tmp := WI0S[J];
                    WI0S[J] := WI0S[J+1];
                    WI0S[J+1] := Tmp;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        APVMove(@WR1S[0], 0, N-1, @WR1[0], 0, N-1);
        APVMove(@WI1S[0], 0, N-1, @WI1[0], 0, N-1);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-2-I do
            begin
                if AP_FP_Greater(WR1S[J],WR1S[J+1]) then
                begin
                    Tmp := WR1S[J];
                    WR1S[J] := WR1S[J+1];
                    WR1S[J+1] := Tmp;
                    Tmp := WI1S[J];
                    WI1S[J] := WI1S[J+1];
                    WI1S[J+1] := Tmp;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            NSErrors := NSErrors or AP_FP_Greater(AbsReal(WR0S[I]-WR1S[I]),Threshold);
            NSErrors := NSErrors or AP_FP_Greater(AbsReal(WI0S[I]-WI1S[I]),Threshold);
            Inc(I);
        end;
        
        //
        // Test right vectors
        //
        if NeedR then
        begin
            K := 0;
            while K<=N-1 do
            begin
                if AP_FP_Eq(WI1[K],0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        Vec1R[i_] := VR[i_,K];
                    end;
                    I:=0;
                    while I<=N-1 do
                    begin
                        Vec1I[I] := 0;
                        Inc(I);
                    end;
                    CurWR := WR1[K];
                    CurWI := 0;
                end;
                if AP_FP_Greater(WI1[K],0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        Vec1R[i_] := VR[i_,K];
                    end;
                    for i_ := 0 to N-1 do
                    begin
                        Vec1I[i_] := VR[i_,K+1];
                    end;
                    CurWR := WR1[K];
                    CurWI := WI1[K];
                end;
                if AP_FP_Less(WI1[K],0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        Vec1R[i_] := VR[i_,K-1];
                    end;
                    for i_ := 0 to N-1 do
                    begin
                        Vec1I[i_] := -VR[i_,K];
                    end;
                    CurWR := WR1[K];
                    CurWI := WI1[K];
                end;
                I:=0;
                while I<=N-1 do
                begin
                    VT := APVDotProduct(@A[I][0], 0, N-1, @Vec1R[0], 0, N-1);
                    Vec2R[I] := VT;
                    VT := APVDotProduct(@A[I][0], 0, N-1, @Vec1I[0], 0, N-1);
                    Vec2I[I] := VT;
                    Inc(I);
                end;
                APVMove(@Vec3R[0], 0, N-1, @Vec1R[0], 0, N-1, CurWR);
                APVSub(@Vec3R[0], 0, N-1, @Vec1I[0], 0, N-1, CurWI);
                APVMove(@Vec3I[0], 0, N-1, @Vec1R[0], 0, N-1, CurWI);
                APVAdd(@Vec3I[0], 0, N-1, @Vec1I[0], 0, N-1, CurWR);
                I:=0;
                while I<=N-1 do
                begin
                    NSErrors := NSErrors or AP_FP_Greater(AbsReal(Vec2R[I]-Vec3R[I]),Threshold);
                    NSErrors := NSErrors or AP_FP_Greater(AbsReal(Vec2I[I]-Vec3I[I]),Threshold);
                    Inc(I);
                end;
                K := K+1;
            end;
        end;
        
        //
        // Test left vectors
        //
        if NeedL then
        begin
            K := 0;
            while K<=N-1 do
            begin
                if AP_FP_Eq(WI1[K],0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        Vec1R[i_] := VL[i_,K];
                    end;
                    I:=0;
                    while I<=N-1 do
                    begin
                        Vec1I[I] := 0;
                        Inc(I);
                    end;
                    CurWR := WR1[K];
                    CurWI := 0;
                end;
                if AP_FP_Greater(WI1[K],0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        Vec1R[i_] := VL[i_,K];
                    end;
                    for i_ := 0 to N-1 do
                    begin
                        Vec1I[i_] := VL[i_,K+1];
                    end;
                    CurWR := WR1[K];
                    CurWI := WI1[K];
                end;
                if AP_FP_Less(WI1[K],0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        Vec1R[i_] := VL[i_,K-1];
                    end;
                    for i_ := 0 to N-1 do
                    begin
                        Vec1I[i_] := -VL[i_,K];
                    end;
                    CurWR := WR1[K];
                    CurWI := WI1[K];
                end;
                J:=0;
                while J<=N-1 do
                begin
                    VT := 0.0;
                    for i_ := 0 to N-1 do
                    begin
                        VT := VT + Vec1R[i_]*A[i_,J];
                    end;
                    Vec2R[J] := VT;
                    VT := 0.0;
                    for i_ := 0 to N-1 do
                    begin
                        VT := VT + Vec1I[i_]*A[i_,J];
                    end;
                    Vec2I[J] := -VT;
                    Inc(J);
                end;
                APVMove(@Vec3R[0], 0, N-1, @Vec1R[0], 0, N-1, CurWR);
                APVAdd(@Vec3R[0], 0, N-1, @Vec1I[0], 0, N-1, CurWI);
                APVMove(@Vec3I[0], 0, N-1, @Vec1R[0], 0, N-1, CurWI);
                APVSub(@Vec3I[0], 0, N-1, @Vec1I[0], 0, N-1, CurWR);
                I:=0;
                while I<=N-1 do
                begin
                    NSErrors := NSErrors or AP_FP_Greater(AbsReal(Vec2R[I]-Vec3R[I]),Threshold);
                    NSErrors := NSErrors or AP_FP_Greater(AbsReal(Vec2I[I]-Vec3I[I]),Threshold);
                    Inc(I);
                end;
                K := K+1;
            end;
        end;
        Inc(VJob);
    end;
end;


(*************************************************************************
Testing EVD subroutines for one N

NOTES:
* BIThreshold is a threshold for bisection-and-inverse-iteration subroutines.
  special threshold is needed because these subroutines may have much more
  larger error than QR-based algorithms.
*************************************************************************)
procedure TestEVDSet(const N : AlglibInteger;
     const Threshold : Double;
     const BIThreshold : Double;
     var FailC : AlglibInteger;
     var Runs : AlglibInteger;
     var NSErrors : Boolean;
     var SErrors : Boolean;
     var HErrors : Boolean;
     var TDErrors : Boolean;
     var SBIErrors : Boolean;
     var HBIErrors : Boolean;
     var TDBIErrors : Boolean);
var
    RA : TReal2DArray;
    RAL : TReal2DArray;
    RAU : TReal2DArray;
    CA : TComplex2DArray;
    CAL : TComplex2DArray;
    CAU : TComplex2DArray;
    D : TReal1DArray;
    E : TReal1DArray;
    Pass : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    MKind : AlglibInteger;
begin
    
    //
    // Test symmetric problems
    //
    
    //
    // Test symmetric problem: zero, random, sparse matrices.
    //
    SetLength(RA, N, N);
    SetLength(RAL, N, N);
    SetLength(RAU, N, N);
    SetLength(CA, N, N);
    SetLength(CAL, N, N);
    SetLength(CAU, N, N);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            RA[I,J] := 0;
            CA[I,J] := C_Complex(0);
            Inc(J);
        end;
        Inc(I);
    end;
    RMatrixSymmetricSplit(RA, N, RAL, RAU);
    CMatrixHermitianSplit(CA, N, CAL, CAU);
    TestSEVDProblem(RA, RAL, RAU, N, Threshold, SErrors, FailC, Runs);
    TestHEVDProblem(CA, CAL, CAU, N, Threshold, HErrors, FailC, Runs);
    TestSEVDBIProblem(RA, RAL, RAU, N, False, BIThreshold, SBIErrors, FailC, Runs);
    TestHEVDBIProblem(CA, CAL, CAU, N, False, BIThreshold, HBIErrors, FailC, Runs);
    I:=0;
    while I<=N-1 do
    begin
        J:=I+1;
        while J<=N-1 do
        begin
            RA[I,J] := 2*RandomReal-1;
            CA[I,J].X := 2*RandomReal-1;
            CA[I,J].Y := 2*RandomReal-1;
            RA[J,I] := RA[I,J];
            CA[J,I] := Conj(CA[I,J]);
            Inc(J);
        end;
        RA[I,I] := 2*RandomReal-1;
        CA[I,I] := C_Complex(2*RandomReal-1);
        Inc(I);
    end;
    RMatrixSymmetricSplit(RA, N, RAL, RAU);
    CMatrixHermitianSplit(CA, N, CAL, CAU);
    TestSEVDProblem(RA, RAL, RAU, N, Threshold, SErrors, FailC, Runs);
    TestHEVDProblem(CA, CAL, CAU, N, Threshold, HErrors, FailC, Runs);
    TestSEVDBIProblem(RA, RAL, RAU, N, True, BIThreshold, SBIErrors, FailC, Runs);
    TestHEVDBIProblem(CA, CAL, CAU, N, True, BIThreshold, HBIErrors, FailC, Runs);
    RMatrixFillSparseA(RA, N, N, 0.995);
    CMatrixFillSparseA(CA, N, N, 0.995);
    I:=0;
    while I<=N-1 do
    begin
        J:=I+1;
        while J<=N-1 do
        begin
            RA[J,I] := RA[I,J];
            CA[J,I] := Conj(CA[I,J]);
            Inc(J);
        end;
        CA[I,I].Y := 0;
        Inc(I);
    end;
    RMatrixSymmetricSplit(RA, N, RAL, RAU);
    CMatrixHermitianSplit(CA, N, CAL, CAU);
    TestSEVDProblem(RA, RAL, RAU, N, Threshold, SErrors, FailC, Runs);
    TestHEVDProblem(CA, CAL, CAU, N, Threshold, HErrors, FailC, Runs);
    TestSEVDBIProblem(RA, RAL, RAU, N, False, BIThreshold, SBIErrors, FailC, Runs);
    TestHEVDBIProblem(CA, CAL, CAU, N, False, BIThreshold, HBIErrors, FailC, Runs);
    
    //
    // testing tridiagonal problems
    //
    MKind:=0;
    while MKind<=4 do
    begin
        SetLength(D, N);
        if N>1 then
        begin
            SetLength(E, N-1);
        end;
        if MKind=0 then
        begin
            
            //
            // Zero matrix
            //
            I:=0;
            while I<=N-1 do
            begin
                D[I] := 0;
                Inc(I);
            end;
            I:=0;
            while I<=N-2 do
            begin
                E[I] := 0;
                Inc(I);
            end;
        end;
        if MKind=1 then
        begin
            
            //
            // Diagonal matrix
            //
            I:=0;
            while I<=N-1 do
            begin
                D[I] := 2*RandomReal-1;
                Inc(I);
            end;
            I:=0;
            while I<=N-2 do
            begin
                E[I] := 0;
                Inc(I);
            end;
        end;
        if MKind=2 then
        begin
            
            //
            // Off-diagonal matrix
            //
            I:=0;
            while I<=N-1 do
            begin
                D[I] := 0;
                Inc(I);
            end;
            I:=0;
            while I<=N-2 do
            begin
                E[I] := 2*RandomReal-1;
                Inc(I);
            end;
        end;
        if MKind=3 then
        begin
            
            //
            // Dense matrix with blocks
            //
            I:=0;
            while I<=N-1 do
            begin
                D[I] := 2*RandomReal-1;
                Inc(I);
            end;
            I:=0;
            while I<=N-2 do
            begin
                E[I] := 2*RandomReal-1;
                Inc(I);
            end;
            J := 1;
            I := 2;
            while J<=N-2 do
            begin
                E[J] := 0;
                J := J+I;
                I := I+1;
            end;
        end;
        if MKind=4 then
        begin
            
            //
            // dense matrix
            //
            I:=0;
            while I<=N-1 do
            begin
                D[I] := 2*RandomReal-1;
                Inc(I);
            end;
            I:=0;
            while I<=N-2 do
            begin
                E[I] := 2*RandomReal-1;
                Inc(I);
            end;
        end;
        TestTDEVDProblem(D, E, N, Threshold, TDErrors, FailC, Runs);
        TestTDEVDBIProblem(D, E, N, (MKind=1) or (MKind=2) or (MKind=4), BIThreshold, TDBIErrors, FailC, Runs);
        Inc(MKind);
    end;
    
    //
    // Test non-symmetric problems
    //
    
    //
    // Test non-symmetric problems: zero, random, sparse matrices.
    //
    SetLength(RA, N, N);
    SetLength(CA, N, N);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            RA[I,J] := 0;
            CA[I,J] := C_Complex(0);
            Inc(J);
        end;
        Inc(I);
    end;
    TestNSEVDProblem(RA, N, Threshold, NSErrors, FailC, Runs);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            RA[I,J] := 2*RandomReal-1;
            CA[I,J].X := 2*RandomReal-1;
            CA[I,J].Y := 2*RandomReal-1;
            Inc(J);
        end;
        Inc(I);
    end;
    TestNSEVDProblem(RA, N, Threshold, NSErrors, FailC, Runs);
    RMatrixFillSparseA(RA, N, N, 0.995);
    CMatrixFillSparseA(CA, N, N, 0.995);
    TestNSEVDProblem(RA, N, Threshold, NSErrors, FailC, Runs);
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testevdunit_test_silent():Boolean;
begin
    Result := TestEVD(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testevdunit_test():Boolean;
begin
    Result := TestEVD(False);
end;


end.