unit testhevdunit;
interface
uses Math, Sysutils, Ap, blas, rotations, tdevd, cblas, creflections, hblas, htridiagonal, hevd;

function TestHEVD(Silent : Boolean):Boolean;
function testhevdunit_test_silent():Boolean;
function testhevdunit_test():Boolean;

implementation

procedure Unset2D(var A : TComplex2DArray);forward;
procedure Unset1D(var A : TReal1DArray);forward;
function TestProduct(const A : TComplex2DArray;
     N : AlglibInteger;
     const Z : TComplex2DArray;
     const Lambda : TReal1DArray):Double;forward;
function TestOrt(const Z : TComplex2DArray; N : AlglibInteger):Double;forward;
procedure TestEVDProblem(const A : TComplex2DArray;
     const AL : TComplex2DArray;
     const AU : TComplex2DArray;
     N : AlglibInteger;
     var MatErr : Double;
     var ValErr : Double;
     var OrtErr : Double;
     var WNSorted : Boolean;
     var FailC : AlglibInteger);forward;


(*************************************************************************
Testing symmetric EVD subroutine
*************************************************************************)
function TestHEVD(Silent : Boolean):Boolean;
var
    A : TComplex2DArray;
    AL : TComplex2DArray;
    AU : TComplex2DArray;
    Z : TComplex2DArray;
    Pass : AlglibInteger;
    N : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    MKind : AlglibInteger;
    PassCount : AlglibInteger;
    MaxN : AlglibInteger;
    MatErr : Double;
    ValErr : Double;
    OrtErr : Double;
    WNSorted : Boolean;
    FailC : AlglibInteger;
    Runs : AlglibInteger;
    FailR : Double;
    FailThreshold : Double;
    Threshold : Double;
    WasErrors : Boolean;
    WFailed : Boolean;
begin
    FailThreshold := 0.005;
    Threshold := 1000*MachineEpsilon;
    MatErr := 0;
    ValErr := 0;
    OrtErr := 0;
    WNSorted := False;
    WFailed := False;
    FailC := 0;
    Runs := 0;
    MaxN := 20;
    PassCount := 10;
    
    //
    // Main cycle
    //
    N:=1;
    while N<=MaxN do
    begin
        
        //
        // Prepare
        //
        SetLength(A, N-1+1, N-1+1);
        SetLength(AL, N-1+1, N-1+1);
        SetLength(AU, N-1+1, N-1+1);
        I:=0;
        while I<=N-1 do
        begin
            J:=I+1;
            while J<=N-1 do
            begin
                
                //
                // A
                //
                A[I,J].X := 2*RandomReal-1;
                A[I,J].Y := 2*RandomReal-1;
                A[J,I] := Conj(A[I,J]);
                
                //
                // A lower
                //
                AL[I,J].X := 2*RandomReal-1;
                AL[I,J].Y := 2*RandomReal-1;
                AL[J,I] := A[J,I];
                
                //
                // A upper
                //
                AU[I,J] := A[I,J];
                AU[J,I].X := 2*RandomReal-1;
                AU[J,I].Y := 2*RandomReal-1;
                Inc(J);
            end;
            A[I,I] := C_Complex(2*RandomReal-1);
            AL[I,I] := A[I,I];
            AU[I,I] := A[I,I];
            Inc(I);
        end;
        
        //
        // Test
        //
        TestEVDProblem(A, AL, AU, N, MatErr, ValErr, OrtErr, WNSorted, FailC);
        Runs := Runs+1;
        Inc(N);
    end;
    
    //
    // report
    //
    FailR := AP_Double(FailC)/Runs;
    WFailed := AP_FP_Greater(FailR,FailThreshold);
    WasErrors := AP_FP_Greater(MatErr,Threshold) or AP_FP_Greater(ValErr,Threshold) or AP_FP_Greater(OrtErr,Threshold) or WNSorted or WFailed;
    if  not Silent then
    begin
        Write(Format('TESTING HERMITIAN EVD'#13#10'',[]));
        Write(Format('EVD matrix error:                        %5.4e'#13#10'',[
            MatErr]));
        Write(Format('EVD values error (different variants):   %5.4e'#13#10'',[
            ValErr]));
        Write(Format('EVD orthogonality error:                 %5.4e'#13#10'',[
            OrtErr]));
        Write(Format('Eigen values order:                      ',[]));
        if  not WNSorted then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('Always converged:                        ',[]));
        if  not WFailed then
        begin
            Write(Format('YES'#13#10'',[]));
        end
        else
        begin
            Write(Format('NO'#13#10'',[]));
            Write(Format('Fail ratio:                              %5.3f'#13#10'',[
                FailR]));
        end;
        Write(Format('Threshold:                               %5.4e'#13#10'',[
            Threshold]));
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
Tests Z*Lambda*Z' against tridiag(D,E).
Returns relative error.
*************************************************************************)
function TestProduct(const A : TComplex2DArray;
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
function TestOrt(const Z : TComplex2DArray; N : AlglibInteger):Double;
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
Tests EVD problem
*************************************************************************)
procedure TestEVDProblem(const A : TComplex2DArray;
     const AL : TComplex2DArray;
     const AU : TComplex2DArray;
     N : AlglibInteger;
     var MatErr : Double;
     var ValErr : Double;
     var OrtErr : Double;
     var WNSorted : Boolean;
     var FailC : AlglibInteger);
var
    Lambda : TReal1DArray;
    LambdaRef : TReal1DArray;
    Z : TComplex2DArray;
    WSucc : Boolean;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
begin
    
    //
    // Test simple EVD: values and full vectors, lower A
    //
    Unset1D(LambdaRef);
    Unset2D(Z);
    WSucc := HMatrixEVD(AL, N, 1, False, LambdaRef, Z);
    if  not WSucc then
    begin
        FailC := FailC+1;
        Exit;
    end;
    MatErr := Max(MatErr, TestProduct(A, N, Z, LambdaRef));
    OrtErr := Max(OrtErr, TestOrt(Z, N));
    I:=0;
    while I<=N-2 do
    begin
        if AP_FP_Less(LambdaRef[I+1],LambdaRef[I]) then
        begin
            WNSorted := True;
        end;
        Inc(I);
    end;
    
    //
    // Test simple EVD: values and full vectors, upper A
    //
    Unset1D(Lambda);
    Unset2D(Z);
    WSucc := HMatrixEVD(AU, N, 1, True, Lambda, Z);
    if  not WSucc then
    begin
        FailC := FailC+1;
        Exit;
    end;
    MatErr := Max(MatErr, TestProduct(A, N, Z, Lambda));
    OrtErr := Max(OrtErr, TestOrt(Z, N));
    I:=0;
    while I<=N-2 do
    begin
        if AP_FP_Less(Lambda[I+1],Lambda[I]) then
        begin
            WNSorted := True;
        end;
        Inc(I);
    end;
    
    //
    // Test simple EVD: values only, lower A
    //
    Unset1D(Lambda);
    Unset2D(Z);
    WSucc := HMatrixEVD(AL, N, 0, False, Lambda, Z);
    if  not WSucc then
    begin
        FailC := FailC+1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        ValErr := Max(ValErr, AbsReal(Lambda[I]-LambdaRef[I]));
        Inc(I);
    end;
    
    //
    // Test simple EVD: values only, upper A
    //
    Unset1D(Lambda);
    Unset2D(Z);
    WSucc := HMatrixEVD(AU, N, 0, True, Lambda, Z);
    if  not WSucc then
    begin
        FailC := FailC+1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        ValErr := Max(ValErr, AbsReal(Lambda[I]-LambdaRef[I]));
        Inc(I);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testhevdunit_test_silent():Boolean;
begin
    Result := TestHEVD(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testhevdunit_test():Boolean;
begin
    Result := TestHEVD(False);
end;


end.