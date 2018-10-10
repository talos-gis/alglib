unit testschurunit;
interface
uses Math, Sysutils, Ap, hblas, reflections, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, hsschur, schur;

function TestSchur(Silent : Boolean):Boolean;
function testschurunit_test_silent():Boolean;
function testschurunit_test():Boolean;

implementation

procedure FillSparseA(var A : TReal2DArray;
     N : AlglibInteger;
     Sparcity : Double);forward;
procedure TestSchurProblem(const A : TReal2DArray;
     N : AlglibInteger;
     var MatErr : Double;
     var OrtErr : Double;
     var ErrStruct : Boolean;
     var WFailed : Boolean);forward;


(*************************************************************************
Testing Schur decomposition subroutine
*************************************************************************)
function TestSchur(Silent : Boolean):Boolean;
var
    A : TReal2DArray;
    N : AlglibInteger;
    MaxN : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    WasErrors : Boolean;
    ErrStruct : Boolean;
    WFailed : Boolean;
    MatErr : Double;
    OrtErr : Double;
    Threshold : Double;
begin
    MatErr := 0;
    OrtErr := 0;
    ErrStruct := False;
    WFailed := False;
    WasErrors := False;
    MaxN := 70;
    PassCount := 1;
    Threshold := 5*100*MachineEpsilon;
    SetLength(A, MaxN-1+1, MaxN-1+1);
    
    //
    // zero matrix, several cases
    //
    I:=0;
    while I<=MaxN-1 do
    begin
        J:=0;
        while J<=MaxN-1 do
        begin
            A[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    N:=1;
    while N<=MaxN do
    begin
        if (N>30) and (N mod 2=0) then
        begin
            Inc(N);
            Continue;
        end;
        TestSchurProblem(A, N, MatErr, OrtErr, ErrStruct, WFailed);
        Inc(N);
    end;
    
    //
    // Dense matrix
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        N:=1;
        while N<=MaxN do
        begin
            if (N>30) and (N mod 2=0) then
            begin
                Inc(N);
                Continue;
            end;
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := 2*RandomReal-1;
                    Inc(J);
                end;
                Inc(I);
            end;
            TestSchurProblem(A, N, MatErr, OrtErr, ErrStruct, WFailed);
            Inc(N);
        end;
        Inc(Pass);
    end;
    
    //
    // Sparse matrices, very sparse matrices, incredible sparse matrices
    //
    Pass:=1;
    while Pass<=1 do
    begin
        N:=1;
        while N<=MaxN do
        begin
            if (N>30) and (N mod 3<>0) then
            begin
                Inc(N);
                Continue;
            end;
            FillSparseA(A, N, 0.8);
            TestSchurProblem(A, N, MatErr, OrtErr, ErrStruct, WFailed);
            FillSparseA(A, N, 0.9);
            TestSchurProblem(A, N, MatErr, OrtErr, ErrStruct, WFailed);
            FillSparseA(A, N, 0.95);
            TestSchurProblem(A, N, MatErr, OrtErr, ErrStruct, WFailed);
            FillSparseA(A, N, 0.997);
            TestSchurProblem(A, N, MatErr, OrtErr, ErrStruct, WFailed);
            Inc(N);
        end;
        Inc(Pass);
    end;
    
    //
    // report
    //
    WasErrors := AP_FP_Greater(MatErr,Threshold) or AP_FP_Greater(OrtErr,Threshold) or ErrStruct or WFailed;
    if  not Silent then
    begin
        Write(Format('TESTING SCHUR DECOMPOSITION'#13#10'',[]));
        Write(Format('Schur decomposition error:               %5.4e'#13#10'',[
            MatErr]));
        Write(Format('Schur orthogonality error:               %5.4e'#13#10'',[
            OrtErr]));
        Write(Format('T matrix structure:                      ',[]));
        if  not ErrStruct then
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
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
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


procedure FillSparseA(var A : TReal2DArray;
     N : AlglibInteger;
     Sparcity : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    I:=0;
    while I<=N-1 do
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


procedure TestSchurProblem(const A : TReal2DArray;
     N : AlglibInteger;
     var MatErr : Double;
     var OrtErr : Double;
     var ErrStruct : Boolean;
     var WFailed : Boolean);
var
    S : TReal2DArray;
    T : TReal2DArray;
    SR : TReal1DArray;
    ASTC : TReal1DArray;
    SASTC : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    LocErr : Double;
    i_ : AlglibInteger;
begin
    SetLength(SR, N-1+1);
    SetLength(ASTC, N-1+1);
    SetLength(SASTC, N-1+1);
    
    //
    // Schur decomposition, convergence test
    //
    SetLength(T, N-1+1, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            T[I,J] := A[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
    if  not RMatrixSchur(T, N, S) then
    begin
        WFailed := True;
        Exit;
    end;
    
    //
    // decomposition error
    //
    LocErr := 0;
    J:=0;
    while J<=N-1 do
    begin
        APVMove(@SR[0], 0, N-1, @S[J][0], 0, N-1);
        K:=0;
        while K<=N-1 do
        begin
            V := APVDotProduct(@T[K][0], 0, N-1, @SR[0], 0, N-1);
            ASTC[K] := V;
            Inc(K);
        end;
        K:=0;
        while K<=N-1 do
        begin
            V := APVDotProduct(@S[K][0], 0, N-1, @ASTC[0], 0, N-1);
            SASTC[K] := V;
            Inc(K);
        end;
        K:=0;
        while K<=N-1 do
        begin
            LocErr := Max(LocErr, AbsReal(SASTC[K]-A[K,J]));
            Inc(K);
        end;
        Inc(J);
    end;
    MatErr := Max(MatErr, LocErr);
    
    //
    // orthogonality error
    //
    LocErr := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + S[i_,I]*S[i_,J];
            end;
            if I<>J then
            begin
                LocErr := Max(LocErr, AbsReal(V));
            end
            else
            begin
                LocErr := Max(LocErr, AbsReal(V-1));
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    OrtErr := Max(OrtErr, LocErr);
    
    //
    // T matrix structure
    //
    J:=0;
    while J<=N-1 do
    begin
        I:=J+2;
        while I<=N-1 do
        begin
            if AP_FP_Neq(T[I,J],0) then
            begin
                ErrStruct := True;
            end;
            Inc(I);
        end;
        Inc(J);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testschurunit_test_silent():Boolean;
begin
    Result := TestSchur(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testschurunit_test():Boolean;
begin
    Result := TestSchur(False);
end;


end.