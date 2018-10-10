unit testcdetunit;
interface
uses Math, Sysutils, Ap, clu, cdet;

function TestCDet(Silent : Boolean):Boolean;
function testcdetunit_test_silent():Boolean;
function testcdetunit_test():Boolean;

implementation

var
    DetErrors : Boolean;

procedure FillSparseA(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Sparcity : Double);forward;
procedure TestProblem(const A : TComplex2DArray; N : AlglibInteger);forward;
procedure MakeACopy(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TComplex2DArray);forward;
function DetTriangle(A : TComplex2DArray;
     const N : AlglibInteger):Complex;forward;


(*************************************************************************
Main unittest subroutine
*************************************************************************)
function TestCDet(Silent : Boolean):Boolean;
var
    MaxN : AlglibInteger;
    GPassCount : AlglibInteger;
    Threshold : Double;
    A : TComplex2DArray;
    N : AlglibInteger;
    GPass : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    WasErrors : Boolean;
begin
    DetErrors := False;
    WasErrors := False;
    MaxN := 8;
    GPassCount := 5;
    Threshold := 5*100*MachineEpsilon;
    SetLength(A, MaxN-1+1, MaxN-1+1);
    
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
        while I<=MaxN-1 do
        begin
            J:=0;
            while J<=MaxN-1 do
            begin
                A[I,J] := C_Complex(0);
                Inc(J);
            end;
            Inc(I);
        end;
        I:=1;
        while I<=MaxN do
        begin
            TestProblem(A, I);
            Inc(I);
        end;
        
        //
        // Dense matrices
        //
        N:=1;
        while N<=MaxN do
        begin
            I:=0;
            while I<=N-1 do
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
            TestProblem(A, N);
            Inc(N);
        end;
        Inc(GPass);
    end;
    
    //
    // report
    //
    WasErrors := DetErrors;
    if  not Silent then
    begin
        Write(Format('TESTING DET'#13#10'',[]));
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
Sparse matrix
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
Problem testing
*************************************************************************)
procedure TestProblem(const A : TComplex2DArray; N : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
    B : TComplex2DArray;
    C : TComplex2DArray;
    Pivots : TInteger1DArray;
    V1 : Complex;
    V2 : Complex;
    VE : Complex;
begin
    SetLength(B, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        J:=1;
        while J<=N do
        begin
            B[I,J] := A[I-1,J-1];
            Inc(J);
        end;
        Inc(I);
    end;
    VE := DetTriangle(B, N);
    V1 := CMatrixDet(A, N);
    MakeACopy(A, N, N, C);
    CMatrixLU(C, N, N, Pivots);
    V2 := CMatrixLUDet(C, Pivots, N);
    if C_NotEqualR(VE,0) then
    begin
        DetErrors := DetErrors or AP_FP_Greater(AbsComplex(C_DivR(C_Sub(V1,VE),Max(AbsComplex(VE), 1))),1.0E-9);
        DetErrors := DetErrors or AP_FP_Greater(AbsComplex(C_DivR(C_Sub(V1,VE),Max(AbsComplex(VE), 1))),1.0E-9);
    end
    else
    begin
        DetErrors := DetErrors or C_NotEqual(V1,VE);
        DetErrors := DetErrors or C_NotEqual(V2,VE);
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
Basic det subroutine
*************************************************************************)
function DetTriangle(A : TComplex2DArray; const N : AlglibInteger):Complex;
var
    i : AlglibInteger;
    j : AlglibInteger;
    k : AlglibInteger;
    l : AlglibInteger;
    f : AlglibInteger;
    z : AlglibInteger;
    t : Complex;
    m1 : Complex;
begin
    A := DynamicArrayCopy(A);
    Result := C_Complex(1);
    k := 1;
    repeat
        m1 := C_Complex(0);
        i := k;
        while i<=n do
        begin
            t := a[i,k];
            if AP_FP_Greater(AbsComplex(t),AbsComplex(m1)) then
            begin
                m1 := t;
                j := i;
            end;
            i := i+1;
        end;
        if AP_FP_Eq(AbsComplex(m1),0) then
        begin
            Result := C_Complex(0);
            k := n+1;
        end
        else
        begin
            if j<>k then
            begin
                Result := C_Opposite(Result);
                l := k;
                while l<=n do
                begin
                    t := a[j,l];
                    a[j,l] := a[k,l];
                    a[k,l] := t;
                    l := l+1;
                end;
            end;
            f := k+1;
            while f<=n do
            begin
                t := C_Div(a[f,k],m1);
                z := k+1;
                while z<=n do
                begin
                    a[f,z] := C_Sub(a[f,z],C_Mul(t,a[k,z]));
                    z := z+1;
                end;
                f := f+1;
            end;
            Result := C_Mul(Result,a[k,k]);
        end;
        k := k+1;
    until  not (k<=n);
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testcdetunit_test_silent():Boolean;
begin
    Result := TestCDet(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testcdetunit_test():Boolean;
begin
    Result := TestCDet(False);
end;


end.