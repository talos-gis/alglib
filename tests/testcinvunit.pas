unit testcinvunit;
interface
uses Math, Sysutils, Ap, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, ctrinverse, cinverse;

function TestCInv(Silent : Boolean):Boolean;
function testcinvunit_test_silent():Boolean;
function testcinvunit_test():Boolean;

implementation

var
    InvErrors : Boolean;
    Threshold : Double;

procedure TestProblem(const A : TComplex2DArray; N : AlglibInteger);forward;
procedure MakeACopy(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TComplex2DArray);forward;
procedure MakeACopyOldMem(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TComplex2DArray);forward;
function MatrixDiff(const A : TComplex2DArray;
     const B : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger):Double;forward;
procedure InternalMatrixMatrixMultiply(const A : TComplex2DArray;
     AI1 : AlglibInteger;
     AI2 : AlglibInteger;
     AJ1 : AlglibInteger;
     AJ2 : AlglibInteger;
     TransA : Boolean;
     const B : TComplex2DArray;
     BI1 : AlglibInteger;
     BI2 : AlglibInteger;
     BJ1 : AlglibInteger;
     BJ2 : AlglibInteger;
     TransB : Boolean;
     var C : TComplex2DArray;
     CI1 : AlglibInteger;
     CI2 : AlglibInteger;
     CJ1 : AlglibInteger;
     CJ2 : AlglibInteger);forward;


(*************************************************************************
Main unittest subroutine
*************************************************************************)
function TestCInv(Silent : Boolean):Boolean;
var
    MaxN : AlglibInteger;
    GPassCount : AlglibInteger;
    A : TComplex2DArray;
    N : AlglibInteger;
    GPass : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    WasErrors : Boolean;
begin
    InvErrors := False;
    Threshold := 5*1000*MachineEpsilon;
    WasErrors := False;
    MaxN := 10;
    GPassCount := 30;
    
    //
    // Different problems
    //
    N:=1;
    while N<=MaxN do
    begin
        SetLength(A, N-1+1, N-1+1);
        GPass:=1;
        while GPass<=GPassCount do
        begin
            
            //
            // diagonal matrix, several cases
            //
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := C_Complex(0);
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=N-1 do
            begin
                A[I,I].X := 2*RandomReal-1;
                A[I,I].Y := 2*RandomReal-1;
                Inc(I);
            end;
            TestProblem(A, N);
            
            //
            // shifted diagonal matrix, several cases
            //
            K := RandomInteger(N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := C_Complex(0);
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=N-1 do
            begin
                A[I,(I+K) mod N].X := 2*RandomReal-1;
                A[I,(I+K) mod N].Y := 2*RandomReal-1;
                Inc(I);
            end;
            TestProblem(A, N);
            
            //
            // Dense matrices
            //
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
            Inc(GPass);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    WasErrors := InvErrors;
    if  not Silent then
    begin
        Write(Format('TESTING COMPLEX INVERSE'#13#10'',[]));
        Write(Format('INVERSE ERRORS                           ',[]));
        if InvErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
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
Problem testing
*************************************************************************)
procedure TestProblem(const A : TComplex2DArray; N : AlglibInteger);
var
    B : TComplex2DArray;
    BLU : TComplex2DArray;
    T1 : TComplex2DArray;
    P : TInteger1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Complex;
begin
    
    //
    // Decomposition
    //
    MakeACopy(A, N, N, B);
    CMatrixInverse(B, N);
    MakeACopy(A, N, N, BLU);
    CMatrixLU(BLU, N, N, P);
    CMatrixLUInverse(BLU, P, N);
    
    //
    // Test
    //
    SetLength(T1, N-1+1, N-1+1);
    InternalMatrixMatrixMultiply(A, 0, N-1, 0, N-1, False, B, 0, N-1, 0, N-1, False, T1, 0, N-1, 0, N-1);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := T1[I,J];
            if I=J then
            begin
                V := C_SubR(V,1);
            end;
            InvErrors := InvErrors or AP_FP_Greater(AbsComplex(V),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
    InternalMatrixMatrixMultiply(A, 0, N-1, 0, N-1, False, BLU, 0, N-1, 0, N-1, False, T1, 0, N-1, 0, N-1);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := T1[I,J];
            if I=J then
            begin
                V := C_SubR(V,1);
            end;
            InvErrors := InvErrors or AP_FP_Greater(AbsComplex(V),Threshold);
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
Copy
*************************************************************************)
procedure MakeACopyOldMem(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TComplex2DArray);
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
            B[I,J] := A[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Diff
*************************************************************************)
function MatrixDiff(const A : TComplex2DArray;
     const B : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Result := 0;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            Result := Max(Result, AbsComplex(C_Sub(B[I,J],A[I,J])));
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Matrix multiplication
*************************************************************************)
procedure InternalMatrixMatrixMultiply(const A : TComplex2DArray;
     AI1 : AlglibInteger;
     AI2 : AlglibInteger;
     AJ1 : AlglibInteger;
     AJ2 : AlglibInteger;
     TransA : Boolean;
     const B : TComplex2DArray;
     BI1 : AlglibInteger;
     BI2 : AlglibInteger;
     BJ1 : AlglibInteger;
     BJ2 : AlglibInteger;
     TransB : Boolean;
     var C : TComplex2DArray;
     CI1 : AlglibInteger;
     CI2 : AlglibInteger;
     CJ1 : AlglibInteger;
     CJ2 : AlglibInteger);
var
    ARows : AlglibInteger;
    ACols : AlglibInteger;
    BRows : AlglibInteger;
    BCols : AlglibInteger;
    CRows : AlglibInteger;
    CCols : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    L : AlglibInteger;
    R : AlglibInteger;
    V : Complex;
    WORK : TComplex1DArray;
    Beta : Complex;
    Alpha : Complex;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // Pre-setup
    //
    K := Max(AI2-AI1+1, AJ2-AJ1+1);
    K := Max(K, BI2-BI1+1);
    K := Max(K, BJ2-BJ1+1);
    SetLength(WORK, K+1);
    Beta := C_Complex(0);
    Alpha := C_Complex(1);
    
    //
    // Setup
    //
    if  not TransA then
    begin
        ARows := AI2-AI1+1;
        ACols := AJ2-AJ1+1;
    end
    else
    begin
        ARows := AJ2-AJ1+1;
        ACols := AI2-AI1+1;
    end;
    if  not TransB then
    begin
        BRows := BI2-BI1+1;
        BCols := BJ2-BJ1+1;
    end
    else
    begin
        BRows := BJ2-BJ1+1;
        BCols := BI2-BI1+1;
    end;
    Assert(ACols=BRows, 'MatrixMatrixMultiply: incorrect matrix sizes!');
    if (ARows<=0) or (ACols<=0) or (BRows<=0) or (BCols<=0) then
    begin
        Exit;
    end;
    CRows := ARows;
    CCols := BCols;
    
    //
    // Test WORK
    //
    I := Max(ARows, ACols);
    I := Max(BRows, I);
    I := Max(I, BCols);
    Work[1] := C_Complex(0);
    Work[I] := C_Complex(0);
    
    //
    // Prepare C
    //
    if C_EqualR(Beta,0) then
    begin
        I:=CI1;
        while I<=CI2 do
        begin
            J:=CJ1;
            while J<=CJ2 do
            begin
                C[I,J] := C_Complex(0);
                Inc(J);
            end;
            Inc(I);
        end;
    end
    else
    begin
        I:=CI1;
        while I<=CI2 do
        begin
            for i_ := CJ1 to CJ2 do
            begin
                C[I,i_] := C_Mul(Beta, C[I,i_]);
            end;
            Inc(I);
        end;
    end;
    
    //
    // A*B
    //
    if  not TransA and  not TransB then
    begin
        L:=AI1;
        while L<=AI2 do
        begin
            R:=BI1;
            while R<=BI2 do
            begin
                V := C_Mul(Alpha,A[L,AJ1+R-BI1]);
                K := CI1+L-AI1;
                i1_ := (BJ1) - (CJ1);
                for i_ := CJ1 to CJ2 do
                begin
                    C[K,i_] := C_Add(C[K,i_], C_Mul(V, B[R,i_+i1_]));
                end;
                Inc(R);
            end;
            Inc(L);
        end;
        Exit;
    end;
    
    //
    // A*B'
    //
    if  not TransA and TransB then
    begin
        if ARows*ACols<BRows*BCols then
        begin
            R:=BI1;
            while R<=BI2 do
            begin
                L:=AI1;
                while L<=AI2 do
                begin
                    i1_ := (BJ1)-(AJ1);
                    V := C_Complex(0.0);
                    for i_ := AJ1 to AJ2 do
                    begin
                        V := C_Add(V,C_Mul(A[L,i_],B[R,i_+i1_]));
                    end;
                    C[CI1+L-AI1,CJ1+R-BI1] := C_Add(C[CI1+L-AI1,CJ1+R-BI1],C_Mul(Alpha,V));
                    Inc(L);
                end;
                Inc(R);
            end;
            Exit;
        end
        else
        begin
            L:=AI1;
            while L<=AI2 do
            begin
                R:=BI1;
                while R<=BI2 do
                begin
                    i1_ := (BJ1)-(AJ1);
                    V := C_Complex(0.0);
                    for i_ := AJ1 to AJ2 do
                    begin
                        V := C_Add(V,C_Mul(A[L,i_],B[R,i_+i1_]));
                    end;
                    C[CI1+L-AI1,CJ1+R-BI1] := C_Add(C[CI1+L-AI1,CJ1+R-BI1],C_Mul(Alpha,V));
                    Inc(R);
                end;
                Inc(L);
            end;
            Exit;
        end;
    end;
    
    //
    // A'*B
    //
    if TransA and  not TransB then
    begin
        L:=AJ1;
        while L<=AJ2 do
        begin
            R:=BI1;
            while R<=BI2 do
            begin
                V := C_Mul(Alpha,A[AI1+R-BI1,L]);
                K := CI1+L-AJ1;
                i1_ := (BJ1) - (CJ1);
                for i_ := CJ1 to CJ2 do
                begin
                    C[K,i_] := C_Add(C[K,i_], C_Mul(V, B[R,i_+i1_]));
                end;
                Inc(R);
            end;
            Inc(L);
        end;
        Exit;
    end;
    
    //
    // A'*B'
    //
    if TransA and TransB then
    begin
        if ARows*ACols<BRows*BCols then
        begin
            R:=BI1;
            while R<=BI2 do
            begin
                I:=1;
                while I<=CRows do
                begin
                    WORK[I] := C_Complex(0.0);
                    Inc(I);
                end;
                L:=AI1;
                while L<=AI2 do
                begin
                    V := C_Mul(Alpha,B[R,BJ1+L-AI1]);
                    K := CJ1+R-BI1;
                    i1_ := (AJ1) - (1);
                    for i_ := 1 to CRows do
                    begin
                        WORK[i_] := C_Add(WORK[i_], C_Mul(V, A[L,i_+i1_]));
                    end;
                    Inc(L);
                end;
                i1_ := (1) - (CI1);
                for i_ := CI1 to CI2 do
                begin
                    C[i_,K] := C_Add(C[i_,K], WORK[i_+i1_]);
                end;
                Inc(R);
            end;
            Exit;
        end
        else
        begin
            L:=AJ1;
            while L<=AJ2 do
            begin
                K := AI2-AI1+1;
                i1_ := (AI1) - (1);
                for i_ := 1 to K do
                begin
                    WORK[i_] := A[i_+i1_,L];
                end;
                R:=BI1;
                while R<=BI2 do
                begin
                    i1_ := (BJ1)-(1);
                    V := C_Complex(0.0);
                    for i_ := 1 to K do
                    begin
                        V := C_Add(V,C_Mul(WORK[i_],B[R,i_+i1_]));
                    end;
                    C[CI1+L-AJ1,CJ1+R-BI1] := C_Add(C[CI1+L-AJ1,CJ1+R-BI1],C_Mul(Alpha,V));
                    Inc(R);
                end;
                Inc(L);
            end;
            Exit;
        end;
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testcinvunit_test_silent():Boolean;
begin
    Result := TestCInv(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testcinvunit_test():Boolean;
begin
    Result := TestCInv(False);
end;


end.