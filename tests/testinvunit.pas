unit testinvunit;
interface
uses Math, Sysutils, Ap, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trinverse, inv;

function TestInv(Silent : Boolean):Boolean;
function testinvunit_test_silent():Boolean;
function testinvunit_test():Boolean;

implementation

var
    InvErrors : Boolean;
    Threshold : Double;

procedure TestProblem(const A : TReal2DArray; N : AlglibInteger);forward;
procedure MakeACopy(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);forward;
procedure MakeACopyOldMem(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);forward;
function MatrixDiff(const A : TReal2DArray;
     const B : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger):Double;forward;
procedure InternalMatrixMatrixMultiply(const A : TReal2DArray;
     AI1 : AlglibInteger;
     AI2 : AlglibInteger;
     AJ1 : AlglibInteger;
     AJ2 : AlglibInteger;
     TransA : Boolean;
     const B : TReal2DArray;
     BI1 : AlglibInteger;
     BI2 : AlglibInteger;
     BJ1 : AlglibInteger;
     BJ2 : AlglibInteger;
     TransB : Boolean;
     var C : TReal2DArray;
     CI1 : AlglibInteger;
     CI2 : AlglibInteger;
     CJ1 : AlglibInteger;
     CJ2 : AlglibInteger);forward;


(*************************************************************************
Main unittest subroutine
*************************************************************************)
function TestInv(Silent : Boolean):Boolean;
var
    MaxN : AlglibInteger;
    GPassCount : AlglibInteger;
    A : TReal2DArray;
    N : AlglibInteger;
    GPass : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    WasErrors : Boolean;
begin
    InvErrors := False;
    Threshold := 1.0E7*MachineEpsilon;
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
                    A[I,J] := 0;
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=N-1 do
            begin
                A[I,I] := 2*RandomReal-1;
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
                    A[I,J] := 0;
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=N-1 do
            begin
                A[I,(I+K) mod N] := 2*RandomReal-1;
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
                    A[I,J] := 2*RandomReal-1;
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
        Write(Format('TESTING INVERSE'#13#10'',[]));
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
procedure TestProblem(const A : TReal2DArray; N : AlglibInteger);
var
    B : TReal2DArray;
    BLU : TReal2DArray;
    T1 : TReal2DArray;
    P : TInteger1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
begin
    
    //
    // Decomposition
    //
    MakeACopy(A, N, N, B);
    RMatrixInverse(B, N);
    MakeACopy(A, N, N, BLU);
    RMatrixLU(BLU, N, N, P);
    RMatrixLUInverse(BLU, P, N);
    
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
                V := V-1;
            end;
            InvErrors := InvErrors or AP_FP_Greater(AbsReal(V),Threshold);
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
                V := V-1;
            end;
            InvErrors := InvErrors or AP_FP_Greater(AbsReal(V),Threshold);
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
Copy
*************************************************************************)
procedure MakeACopyOldMem(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);
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
function MatrixDiff(const A : TReal2DArray;
     const B : TReal2DArray;
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
            Result := Max(Result, AbsReal(B[I,J]-A[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Matrix multiplication
*************************************************************************)
procedure InternalMatrixMatrixMultiply(const A : TReal2DArray;
     AI1 : AlglibInteger;
     AI2 : AlglibInteger;
     AJ1 : AlglibInteger;
     AJ2 : AlglibInteger;
     TransA : Boolean;
     const B : TReal2DArray;
     BI1 : AlglibInteger;
     BI2 : AlglibInteger;
     BJ1 : AlglibInteger;
     BJ2 : AlglibInteger;
     TransB : Boolean;
     var C : TReal2DArray;
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
    V : Double;
    WORK : TReal1DArray;
    Beta : Double;
    Alpha : Double;
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
    Beta := 0;
    Alpha := 1;
    
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
    Work[1] := 0;
    Work[I] := 0;
    
    //
    // Prepare C
    //
    if AP_FP_Eq(Beta,0) then
    begin
        I:=CI1;
        while I<=CI2 do
        begin
            J:=CJ1;
            while J<=CJ2 do
            begin
                C[I,J] := 0;
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
            APVMul(@C[I][0], CJ1, CJ2, Beta);
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
                V := Alpha*A[L,AJ1+R-BI1];
                K := CI1+L-AI1;
                APVAdd(@C[K][0], CJ1, CJ2, @B[R][0], BJ1, BJ2, V);
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
                    V := APVDotProduct(@A[L][0], AJ1, AJ2, @B[R][0], BJ1, BJ2);
                    C[CI1+L-AI1,CJ1+R-BI1] := C[CI1+L-AI1,CJ1+R-BI1]+Alpha*V;
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
                    V := APVDotProduct(@A[L][0], AJ1, AJ2, @B[R][0], BJ1, BJ2);
                    C[CI1+L-AI1,CJ1+R-BI1] := C[CI1+L-AI1,CJ1+R-BI1]+Alpha*V;
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
                V := Alpha*A[AI1+R-BI1,L];
                K := CI1+L-AJ1;
                APVAdd(@C[K][0], CJ1, CJ2, @B[R][0], BJ1, BJ2, V);
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
                    WORK[I] := 0.0;
                    Inc(I);
                end;
                L:=AI1;
                while L<=AI2 do
                begin
                    V := Alpha*B[R,BJ1+L-AI1];
                    K := CJ1+R-BI1;
                    APVAdd(@WORK[0], 1, CRows, @A[L][0], AJ1, AJ2, V);
                    Inc(L);
                end;
                i1_ := (1) - (CI1);
                for i_ := CI1 to CI2 do
                begin
                    C[i_,K] := C[i_,K] + WORK[i_+i1_];
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
                    V := APVDotProduct(@WORK[0], 1, K, @B[R][0], BJ1, BJ2);
                    C[CI1+L-AJ1,CJ1+R-BI1] := C[CI1+L-AJ1,CJ1+R-BI1]+Alpha*V;
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
function testinvunit_test_silent():Boolean;
begin
    Result := TestInv(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testinvunit_test():Boolean;
begin
    Result := TestInv(False);
end;


end.