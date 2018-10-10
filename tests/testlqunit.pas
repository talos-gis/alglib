unit testlqunit;
interface
uses Math, Sysutils, Ap, reflections, lq;

function TestLQ(Silent : Boolean):Boolean;
function testlqunit_test_silent():Boolean;
function testlqunit_test():Boolean;

implementation

var
    Threshold : Double;
    StructErrors : Boolean;
    DecompErrors : Boolean;
    OtherErrors : Boolean;

procedure FillSparseA(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Sparcity : Double);forward;
procedure MakeACopy(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);forward;
procedure TestProblem(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger);forward;


(*************************************************************************
Main unittest subroutine
*************************************************************************)
function TestLQ(Silent : Boolean):Boolean;
var
    ShortMN : AlglibInteger;
    MaxMN : AlglibInteger;
    GPassCount : AlglibInteger;
    A : TReal2DArray;
    M : AlglibInteger;
    N : AlglibInteger;
    GPass : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    WasErrors : Boolean;
begin
    DecompErrors := False;
    OtherErrors := False;
    StructErrors := False;
    WasErrors := False;
    ShortMN := 5;
    MaxMN := 15;
    GPassCount := 5;
    Threshold := 5*100*MachineEpsilon;
    SetLength(A, MaxMN-1+1, MaxMN-1+1);
    
    //
    // Different problems
    //
    GPass:=1;
    while GPass<=GPassCount do
    begin
        
        //
        // zero matrix, several cases
        //
        I:=0;
        while I<=MaxMN-1 do
        begin
            J:=0;
            while J<=MaxMN-1 do
            begin
                A[I,J] := 0;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=1;
        while I<=MaxMN do
        begin
            J:=1;
            while J<=MaxMN do
            begin
                TestProblem(A, I, J);
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Long dense matrix
        //
        I:=0;
        while I<=MaxMN-1 do
        begin
            J:=0;
            while J<=ShortMN-1 do
            begin
                A[I,J] := 2*RandomReal-1;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=ShortMN+1;
        while I<=MaxMN do
        begin
            TestProblem(A, I, ShortMN);
            Inc(I);
        end;
        I:=0;
        while I<=ShortMN-1 do
        begin
            J:=0;
            while J<=MaxMN-1 do
            begin
                A[I,J] := 2*RandomReal-1;
                Inc(J);
            end;
            Inc(I);
        end;
        J:=ShortMN+1;
        while J<=MaxMN do
        begin
            TestProblem(A, ShortMN, J);
            Inc(J);
        end;
        
        //
        // Dense matrices
        //
        M:=1;
        while M<=MaxMN do
        begin
            N:=1;
            while N<=MaxMN do
            begin
                I:=0;
                while I<=M-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        A[I,J] := 2*RandomReal-1;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                TestProblem(A, M, N);
                Inc(N);
            end;
            Inc(M);
        end;
        
        //
        // Sparse matrices, very sparse matrices, incredible sparse matrices
        //
        M:=1;
        while M<=MaxMN do
        begin
            N:=1;
            while N<=MaxMN do
            begin
                FillSparseA(A, M, N, 0.8);
                TestProblem(A, M, N);
                FillSparseA(A, M, N, 0.9);
                TestProblem(A, M, N);
                FillSparseA(A, M, N, 0.95);
                TestProblem(A, M, N);
                Inc(N);
            end;
            Inc(M);
        end;
        Inc(GPass);
    end;
    
    //
    // report
    //
    WasErrors := StructErrors or DecompErrors or OtherErrors;
    if  not Silent then
    begin
        Write(Format('TESTING RMatrixLQ'#13#10'',[]));
        Write(Format('STRUCTURAL ERRORS:                       ',[]));
        if  not StructErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('DECOMPOSITION ERRORS:                    ',[]));
        if  not DecompErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('OTHER ERRORS:                            ',[]));
        if  not OtherErrors then
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
procedure FillSparseA(var A : TReal2DArray;
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
Copy
*************************************************************************)
procedure MakeACopy(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    SetLength(B, M-1+1, N-1+1);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            B[I,J] := A[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Problem testing
*************************************************************************)
procedure TestProblem(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    MX : Double;
    B : TReal2DArray;
    TauB : TReal1DArray;
    Q : TReal2DArray;
    L : TReal2DArray;
    Q2 : TReal2DArray;
    V : Double;
    i_ : AlglibInteger;
begin
    
    //
    // MX - estimate of the matrix norm
    //
    MX := 0;
    I:=0;
    while I<=M-1 do
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
    // Test decompose-and-unpack error
    //
    MakeACopy(A, M, N, B);
    RMatrixLQ(B, M, N, TauB);
    RMatrixLQUnpackQ(B, M, N, TauB, N, Q);
    RMatrixLQUnpackL(B, M, N, L);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + L[I,i_]*Q[i_,J];
            end;
            DecompErrors := DecompErrors or AP_FP_Greater_Eq(AbsReal(V-A[I,J]),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        J:=Min(I, N-1)+1;
        while J<=N-1 do
        begin
            StructErrors := StructErrors or AP_FP_Neq(L[I,J],0);
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := APVDotProduct(@Q[I][0], 0, N-1, @Q[J][0], 0, N-1);
            if I=J then
            begin
                StructErrors := StructErrors or AP_FP_Greater_Eq(AbsReal(V-1),Threshold);
            end
            else
            begin
                StructErrors := StructErrors or AP_FP_Greater_Eq(AbsReal(V),Threshold);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Test for other errors
    //
    K:=1;
    while K<=N-1 do
    begin
        RMatrixLQUnpackQ(B, M, N, TauB, K, Q2);
        I:=0;
        while I<=K-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                OtherErrors := OtherErrors or AP_FP_Neq(Q2[I,J],Q[I,J]);
                Inc(J);
            end;
            Inc(I);
        end;
        Inc(K);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testlqunit_test_silent():Boolean;
begin
    Result := TestLQ(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testlqunit_test():Boolean;
begin
    Result := TestLQ(False);
end;


end.