unit testtdevdunit;
interface
uses Math, Sysutils, Ap, blas, rotations, tdevd;

function TestTDEVD(Silent : Boolean):Boolean;
function testtdevdunit_test_silent():Boolean;
function testtdevdunit_test():Boolean;

implementation

procedure FillDE(var D : TReal1DArray;
     var E : TReal1DArray;
     N : AlglibInteger;
     FillType : AlglibInteger);forward;
function TestProduct(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     const Z : TReal2DArray;
     const Lambda : TReal1DArray):Double;forward;
function TestOrt(const Z : TReal2DArray; N : AlglibInteger):Double;forward;
procedure TestEVDProblem(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     var MatErr : Double;
     var ValErr : Double;
     var OrtErr : Double;
     var WNSorted : Boolean;
     var FailC : AlglibInteger);forward;


(*************************************************************************
Testing bidiagonal SVD decomposition subroutine
*************************************************************************)
function TestTDEVD(Silent : Boolean):Boolean;
var
    D : TReal1DArray;
    E : TReal1DArray;
    Pass : AlglibInteger;
    N : AlglibInteger;
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
        SetLength(D, N-1+1);
        if N>1 then
        begin
            SetLength(E, N-2+1);
        end;
        
        //
        // Different tasks
        //
        MKind:=0;
        while MKind<=4 do
        begin
            FillDE(D, E, N, MKind);
            TestEVDProblem(D, E, N, MatErr, ValErr, OrtErr, WNSorted, FailC);
            Runs := Runs+1;
            Inc(MKind);
        end;
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
        Write(Format('TESTING TRIDIAGONAL EVD'#13#10'',[]));
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
Fills D and E
*************************************************************************)
procedure FillDE(var D : TReal1DArray;
     var E : TReal1DArray;
     N : AlglibInteger;
     FillType : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    if FillType=0 then
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
        Exit;
    end;
    if FillType=1 then
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
        Exit;
    end;
    if FillType=2 then
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
        Exit;
    end;
    if FillType=3 then
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
        Exit;
    end;
    
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


(*************************************************************************
Tests Z*Lambda*Z' against tridiag(D,E).
Returns relative error.
*************************************************************************)
function TestProduct(const D : TReal1DArray;
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
procedure TestEVDProblem(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     var MatErr : Double;
     var ValErr : Double;
     var OrtErr : Double;
     var WNSorted : Boolean;
     var FailC : AlglibInteger);
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
    WSucc := SMatrixTDEVD(Lambda, EE, N, 2, Z);
    if  not WSucc then
    begin
        FailC := FailC+1;
        Exit;
    end;
    MatErr := Max(MatErr, TestProduct(D, E, N, Z, Lambda));
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
    WSucc := SMatrixTDEVD(Lambda2, EE, N, 0, Z);
    if  not WSucc then
    begin
        FailC := FailC+1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        ValErr := Max(ValErr, AbsReal(Lambda2[I]-Lambda[I]));
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
    WSucc := SMatrixTDEVD(Lambda2, EE, N, 1, A1);
    if  not WSucc then
    begin
        FailC := FailC+1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        ValErr := Max(ValErr, AbsReal(Lambda2[I]-Lambda[I]));
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
            MatErr := Max(MatErr, AbsReal(V-A1[I,J]));
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
    WSucc := SMatrixTDEVD(Lambda2, EE, N, 3, Z);
    if  not WSucc then
    begin
        FailC := FailC+1;
        Exit;
    end;
    I:=0;
    while I<=N-1 do
    begin
        ValErr := Max(ValErr, AbsReal(Lambda2[I]-Lambda[I]));
        MatErr := Max(MatErr, AbsReal(Z[0,I]-ZRef[0,I]));
        Inc(I);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testtdevdunit_test_silent():Boolean;
begin
    Result := TestTDEVD(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testtdevdunit_test():Boolean;
begin
    Result := TestTDEVD(False);
end;


end.