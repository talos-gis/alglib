unit testcqrunit;
interface
uses Math, Sysutils, Ap, creflections, cqr;

function TestCQR(Silent : Boolean):Boolean;
function testcqrunit_test_silent():Boolean;
function testcqrunit_test():Boolean;

implementation

var
    Threshold : Double;
    StructErrors : Boolean;
    DecompErrors : Boolean;
    OtherErrors : Boolean;

procedure FillSparseA(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Sparcity : Double);forward;
procedure MakeACopy(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TComplex2DArray);forward;
procedure TestProblem(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger);forward;


(*************************************************************************
Main unittest subroutine
*************************************************************************)
function TestCQR(Silent : Boolean):Boolean;
var
    ShortMN : AlglibInteger;
    MaxMN : AlglibInteger;
    GPassCount : AlglibInteger;
    A : TComplex2DArray;
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
                A[I,J] := C_Complex(0);
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
                A[I,J].X := 2*RandomReal-1;
                A[I,J].Y := 2*RandomReal-1;
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
                A[I,J].X := 2*RandomReal-1;
                A[I,J].Y := 2*RandomReal-1;
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
                        A[I,J].X := 2*RandomReal-1;
                        A[I,J].Y := 2*RandomReal-1;
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
        Write(Format('TESTING CMatrixQR'#13#10'',[]));
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
procedure FillSparseA(var A : TComplex2DArray;
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
Copy
*************************************************************************)
procedure MakeACopy(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TComplex2DArray);
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
procedure TestProblem(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    MX : Double;
    B : TComplex2DArray;
    TauB : TComplex1DArray;
    Q : TComplex2DArray;
    R : TComplex2DArray;
    Q2 : TComplex2DArray;
    V : Complex;
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
            if AP_FP_Greater(AbsComplex(A[I,J]),MX) then
            begin
                MX := AbsComplex(A[I,J]);
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
    CMatrixQR(B, M, N, TauB);
    CMatrixQRUnpackQ(B, M, N, TauB, M, Q);
    CMatrixQRUnpackR(B, M, N, R);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := C_Complex(0.0);
            for i_ := 0 to M-1 do
            begin
                V := C_Add(V,C_Mul(Q[I,i_],R[i_,J]));
            end;
            DecompErrors := DecompErrors or AP_FP_Greater_Eq(AbsComplex(C_Sub(V,A[I,J])),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=Min(I, N-1)-1 do
        begin
            StructErrors := StructErrors or C_NotEqualR(R[I,J],0);
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=M-1 do
        begin
            V := C_Complex(0.0);
            for i_ := 0 to M-1 do
            begin
                V := C_Add(V,C_Mul(Q[I,i_],Conj(Q[J,i_])));
            end;
            if I=J then
            begin
                StructErrors := StructErrors or AP_FP_Greater_Eq(AbsComplex(C_SubR(V,1)),Threshold);
            end
            else
            begin
                StructErrors := StructErrors or AP_FP_Greater_Eq(AbsComplex(V),Threshold);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Test for other errors
    //
    K:=1;
    while K<=M-1 do
    begin
        CMatrixQRUnpackQ(B, M, N, TauB, K, Q2);
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=K-1 do
            begin
                OtherErrors := OtherErrors or AP_FP_Greater(AbsComplex(C_Sub(Q2[I,J],Q[I,J])),10*MachineEpsilon);
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
function testcqrunit_test_silent():Boolean;
begin
    Result := TestCQR(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testcqrunit_test():Boolean;
begin
    Result := TestCQR(False);
end;


end.