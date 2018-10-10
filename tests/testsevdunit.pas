unit testsevdunit;
interface
uses Math, Sysutils, Ap, blas, rotations, tdevd, sblas, reflections, tridiagonal, sevd;

function TestSEVD(Silent : Boolean):Boolean;
function testsevdunit_test_silent():Boolean;
function testsevdunit_test():Boolean;

implementation

procedure Unset2D(var A : TReal2DArray);forward;
procedure Unset1D(var A : TReal1DArray);forward;
function TestProduct(const A : TReal2DArray;
     N : AlglibInteger;
     const Z : TReal2DArray;
     const Lambda : TReal1DArray):Double;forward;
function TestOrt(const Z : TReal2DArray; N : AlglibInteger):Double;forward;
procedure TestEVDProblem(const A : TReal2DArray;
     const AL : TReal2DArray;
     const AU : TReal2DArray;
     N : AlglibInteger;
     var MatErr : Double;
     var ValErr : Double;
     var OrtErr : Double;
     var WNSorted : Boolean;
     var FailC : AlglibInteger);forward;


(*************************************************************************
Testing symmetric EVD subroutine
*************************************************************************)
function TestSEVD(Silent : Boolean):Boolean;
var
    A : TReal2DArray;
    AL : TReal2DArray;
    AU : TReal2DArray;
    Z : TReal2DArray;
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
                A[I,J] := 2*RandomReal-1;
                A[J,I] := A[I,J];
                
                //
                // A lower
                //
                AL[I,J] := 2*RandomReal-1;
                AL[J,I] := A[I,J];
                
                //
                // A upper
                //
                AU[I,J] := A[I,J];
                AU[J,I] := 2*RandomReal-1;
                Inc(J);
            end;
            A[I,I] := 2*RandomReal-1;
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
        Write(Format('TESTING SYMMETRIC EVD'#13#10'',[]));
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
procedure Unset2D(var A : TReal2DArray);
begin
    SetLength(A, 0+1, 0+1);
    A[0,0] := 2*RandomReal-1;
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
Tests EVD problem
*************************************************************************)
procedure TestEVDProblem(const A : TReal2DArray;
     const AL : TReal2DArray;
     const AU : TReal2DArray;
     N : AlglibInteger;
     var MatErr : Double;
     var ValErr : Double;
     var OrtErr : Double;
     var WNSorted : Boolean;
     var FailC : AlglibInteger);
var
    Lambda : TReal1DArray;
    LambdaRef : TReal1DArray;
    Z : TReal2DArray;
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
    WSucc := SMatrixEVD(AL, N, 1, False, LambdaRef, Z);
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
    WSucc := SMatrixEVD(AU, N, 1, True, Lambda, Z);
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
    WSucc := SMatrixEVD(AL, N, 0, False, Lambda, Z);
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
    WSucc := SMatrixEVD(AU, N, 0, True, Lambda, Z);
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
function testsevdunit_test_silent():Boolean;
begin
    Result := TestSEVD(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testsevdunit_test():Boolean;
begin
    Result := TestSEVD(False);
end;


end.