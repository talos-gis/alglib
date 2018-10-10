unit testcluunit;
interface
uses Math, Sysutils, Ap, clu;

function TestCLU(Silent : Boolean):Boolean;
function testcluunit_test_silent():Boolean;
function testcluunit_test():Boolean;

implementation

procedure TestCLUProblem(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var DiffPU : Double;
     var LUErr : Double);forward;


function TestCLU(Silent : Boolean):Boolean;
var
    A : TComplex2DArray;
    M : AlglibInteger;
    N : AlglibInteger;
    MaxMN : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    Pass : AlglibInteger;
    WasErrors : Boolean;
    DiffErr : Double;
    LUErr : Double;
    Threshold : Double;
begin
    DiffErr := 0;
    LUErr := 0;
    WasErrors := False;
    MaxMN := 50;
    Threshold := MachineEpsilon*MaxMN;
    SetLength(A, MaxMN+1, MaxMN+1);
    
    //
    // zero matrix, several cases
    //
    I:=1;
    while I<=MaxMN do
    begin
        J:=1;
        while J<=MaxMN do
        begin
            A[I,J] := C_Complex(0);
            Inc(J);
        end;
        Inc(I);
    end;
    I:=1;
    while I<=5 do
    begin
        J:=1;
        while J<=5 do
        begin
            TestCLUProblem(A, I, J, DiffErr, LUErr);
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Long non-zero matrix
    //
    I:=1;
    while I<=MaxMN do
    begin
        J:=1;
        while J<=5 do
        begin
            A[I,J].X := 2*RandomReal-1;
            A[I,J].Y := 2*RandomReal-1;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=1;
    while I<=MaxMN do
    begin
        J:=1;
        while J<=5 do
        begin
            TestCLUProblem(A, I, J, DiffErr, LUErr);
            Inc(J);
        end;
        Inc(I);
    end;
    I:=1;
    while I<=5 do
    begin
        J:=1;
        while J<=MaxMN do
        begin
            A[I,J].X := 2*RandomReal-1;
            A[I,J].Y := 2*RandomReal-1;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=1;
    while I<=5 do
    begin
        J:=1;
        while J<=MaxMN do
        begin
            TestCLUProblem(A, I, J, DiffErr, LUErr);
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Sparse matrices
    //
    M:=1;
    while M<=10 do
    begin
        N:=1;
        while N<=10 do
        begin
            Pass:=1;
            while Pass<=5 do
            begin
                I:=1;
                while I<=M do
                begin
                    J:=1;
                    while J<=N do
                    begin
                        if AP_FP_Greater(RandomReal,0.8) then
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
                TestCLUProblem(A, M, N, DiffErr, LUErr);
                Inc(Pass);
            end;
            Inc(N);
        end;
        Inc(M);
    end;
    
    //
    // Dense matrices
    //
    M:=1;
    while M<=10 do
    begin
        N:=1;
        while N<=10 do
        begin
            I:=1;
            while I<=M do
            begin
                J:=1;
                while J<=N do
                begin
                    A[I,J].X := 2*RandomReal-1;
                    A[I,J].Y := 2*RandomReal-1;
                    Inc(J);
                end;
                Inc(I);
            end;
            TestCLUProblem(A, M, N, DiffErr, LUErr);
            Inc(N);
        end;
        Inc(M);
    end;
    
    //
    // report
    //
    WasErrors := AP_FP_Greater(DiffErr,Threshold) or AP_FP_Greater(LUErr,Threshold);
    if  not Silent then
    begin
        Write(Format('TESTING COMPLEX LU DECOMPOSITION'#13#10'',[]));
        Write(Format('Difference (normal/packed/0-based LU):   %5.4e'#13#10'',[
            DiffErr]));
        Write(Format('LU decomposition error:                  %5.4e'#13#10'',[
            LUErr]));
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


procedure TestCLUProblem(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var DiffPU : Double;
     var LUErr : Double);
var
    T1 : TComplex2DArray;
    T2 : TComplex2DArray;
    T3 : TComplex2DArray;
    IT1 : TInteger1DArray;
    IT2 : TInteger1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Complex;
    A0 : TComplex2DArray;
    P0 : TInteger1DArray;
    MX : Double;
    i_ : AlglibInteger;
begin
    MX := 0;
    I:=1;
    while I<=M do
    begin
        J:=1;
        while J<=N do
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
    // Compare LU and unpacked LU
    //
    SetLength(T1, M+1, N+1);
    I:=1;
    while I<=M do
    begin
        for i_ := 1 to N do
        begin
            T1[I,i_] := A[I,i_];
        end;
        Inc(I);
    end;
    ComplexLUDecomposition(T1, M, N, IT1);
    ComplexLUDecompositionUnpacked(A, M, N, T2, T3, IT2);
    I:=1;
    while I<=M do
    begin
        J:=1;
        while J<=Min(M, N) do
        begin
            if I>J then
            begin
                DiffPU := Max(DiffPU, AbsComplex(C_Sub(T1[I,J],T2[I,J]))/MX);
            end;
            if I=J then
            begin
                DiffPU := Max(DiffPU, AbsComplex(C_RSub(1,T2[I,J]))/MX);
            end;
            if I<J then
            begin
                DiffPU := Max(DiffPU, AbsComplex(C_RSub(0,T2[I,J]))/MX);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=1;
    while I<=Min(M, N) do
    begin
        J:=1;
        while J<=N do
        begin
            if I>J then
            begin
                DiffPU := Max(DiffPU, AbsComplex(C_RSub(0,T3[I,J]))/MX);
            end;
            if I<=J then
            begin
                DiffPU := Max(DiffPU, AbsComplex(C_Sub(T1[I,J],T3[I,J]))/MX);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=1;
    while I<=Min(M, N) do
    begin
        DiffPU := Max(DiffPU, AbsComplex(C_Complex(IT1[I]-IT2[I])));
        Inc(I);
    end;
    
    //
    // Test unpacked LU
    //
    ComplexLUDecompositionUnpacked(A, M, N, T1, T2, IT1);
    SetLength(T3, M+1, N+1);
    K := Min(M, N);
    I:=1;
    while I<=M do
    begin
        J:=1;
        while J<=N do
        begin
            V := C_Complex(0.0);
            for i_ := 1 to K do
            begin
                V := C_Add(V,C_Mul(T1[I,i_],T2[i_,J]));
            end;
            T3[I,J] := V;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=Min(M, N);
    while I>=1 do
    begin
        if I<>IT1[I] then
        begin
            J:=1;
            while J<=N do
            begin
                V := T3[I,J];
                T3[I,J] := T3[IT1[I],J];
                T3[IT1[I],J] := V;
                Inc(J);
            end;
        end;
        Dec(I);
    end;
    I:=1;
    while I<=M do
    begin
        J:=1;
        while J<=N do
        begin
            LUErr := Max(LUErr, AbsComplex(C_Sub(A[I,J],T3[I,J]))/MX);
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Test 0-based LU
    //
    SetLength(T1, M+1, N+1);
    I:=1;
    while I<=M do
    begin
        for i_ := 1 to N do
        begin
            T1[I,i_] := A[I,i_];
        end;
        Inc(I);
    end;
    ComplexLUDecomposition(T1, M, N, IT1);
    SetLength(A0, M-1+1, N-1+1);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            A0[I,J] := A[I+1,J+1];
            Inc(J);
        end;
        Inc(I);
    end;
    CMatrixLU(A0, M, N, P0);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            DiffPU := Max(DiffPU, AbsComplex(C_Sub(A0[I,J],T1[I+1,J+1])));
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=Min(M-1, N-1) do
    begin
        DiffPU := Max(DiffPU, AbsComplex(C_Complex(P0[I]+1-IT1[I+1])));
        Inc(I);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testcluunit_test_silent():Boolean;
begin
    Result := TestCLU(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testcluunit_test():Boolean;
begin
    Result := TestCLU(False);
end;


end.