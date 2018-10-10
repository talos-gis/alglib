unit testbdunit;
interface
uses Math, Sysutils, Ap, reflections, bidiagonal;

function TestBD(Silent : Boolean):Boolean;
function testbdunit_test_silent():Boolean;
function testbdunit_test():Boolean;

implementation

var
    DecompErrors : Boolean;
    PropErrors : Boolean;
    PartErrors : Boolean;
    MulErrors : Boolean;
    Threshold : Double;

procedure FillSparseA(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Sparcity : Double);forward;
procedure TestProblem(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger);forward;
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
function TestBD(Silent : Boolean):Boolean;
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
    PropErrors := False;
    MulErrors := False;
    PartErrors := False;
    Threshold := 5*100*MachineEpsilon;
    WasErrors := False;
    ShortMN := 5;
    MaxMN := 10;
    GPassCount := 5;
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
    WasErrors := DecompErrors or PropErrors or PartErrors or MulErrors;
    if  not Silent then
    begin
        Write(Format('TESTING 2BIDIAGONAL'#13#10'',[]));
        Write(Format('DECOMPOSITION ERRORS                     ',[]));
        if DecompErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('MATRIX PROPERTIES                        ',[]));
        if PropErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('PARTIAL UNPACKING                        ',[]));
        if PartErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('MULTIPLICATION TEST                      ',[]));
        if MulErrors then
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
Sparse matrix
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
    T : TReal2DArray;
    PT : TReal2DArray;
    Q : TReal2DArray;
    R : TReal2DArray;
    BD : TReal2DArray;
    X : TReal2DArray;
    R1 : TReal2DArray;
    R2 : TReal2DArray;
    TauP : TReal1DArray;
    TauQ : TReal1DArray;
    D : TReal1DArray;
    E : TReal1DArray;
    Up : Boolean;
    V : Double;
    MTSize : AlglibInteger;
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
    // Bidiagonal decomposition error
    //
    MakeACopy(A, M, N, T);
    RMatrixBD(T, M, N, TauQ, TauP);
    RMatrixBDUnpackQ(T, M, N, TauQ, M, Q);
    RMatrixBDUnpackPT(T, M, N, TauP, N, PT);
    RMatrixBDUnpackDiagonals(T, M, N, Up, D, E);
    SetLength(BD, M-1+1, N-1+1);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            BD[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=Min(M, N)-1 do
    begin
        BD[I,I] := D[I];
        Inc(I);
    end;
    if Up then
    begin
        I:=0;
        while I<=Min(M, N)-2 do
        begin
            BD[I,I+1] := E[I];
            Inc(I);
        end;
    end
    else
    begin
        I:=0;
        while I<=Min(M, N)-2 do
        begin
            BD[I+1,I] := E[I];
            Inc(I);
        end;
    end;
    SetLength(R, M-1+1, N-1+1);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := 0.0;
            for i_ := 0 to M-1 do
            begin
                V := V + Q[I,i_]*BD[i_,J];
            end;
            R[I,J] := V;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + R[I,i_]*PT[i_,J];
            end;
            DecompErrors := DecompErrors or AP_FP_Greater(AbsReal(V-A[I,J]),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Orthogonality test for Q/PT
    //
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=M-1 do
        begin
            V := 0.0;
            for i_ := 0 to M-1 do
            begin
                V := V + Q[i_,I]*Q[i_,J];
            end;
            if I=J then
            begin
                PropErrors := PropErrors or AP_FP_Greater(AbsReal(V-1),Threshold);
            end
            else
            begin
                PropErrors := PropErrors or AP_FP_Greater(AbsReal(V),Threshold);
            end;
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
            V := APVDotProduct(@PT[I][0], 0, N-1, @PT[J][0], 0, N-1);
            if I=J then
            begin
                PropErrors := PropErrors or AP_FP_Greater(AbsReal(V-1),Threshold);
            end
            else
            begin
                PropErrors := PropErrors or AP_FP_Greater(AbsReal(V),Threshold);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Partial unpacking test
    //
    K:=1;
    while K<=M-1 do
    begin
        RMatrixBDUnpackQ(T, M, N, TauQ, K, R);
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=K-1 do
            begin
                PartErrors := PartErrors or AP_FP_Greater(AbsReal(R[I,J]-Q[I,J]),10*MachineEpsilon);
                Inc(J);
            end;
            Inc(I);
        end;
        Inc(K);
    end;
    K:=1;
    while K<=N-1 do
    begin
        RMatrixBDUnpackPT(T, M, N, TauP, K, R);
        I:=0;
        while I<=K-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                PartErrors := PartErrors or AP_FP_Neq(R[I,J]-PT[I,J],0);
                Inc(J);
            end;
            Inc(I);
        end;
        Inc(K);
    end;
    
    //
    // Multiplication test
    //
    SetLength(X, Max(M, N)-1+1, Max(M, N)-1+1);
    SetLength(R, Max(M, N)-1+1, Max(M, N)-1+1);
    SetLength(R1, Max(M, N)-1+1, Max(M, N)-1+1);
    SetLength(R2, Max(M, N)-1+1, Max(M, N)-1+1);
    I:=0;
    while I<=Max(M, N)-1 do
    begin
        J:=0;
        while J<=Max(M, N)-1 do
        begin
            X[I,J] := 2*RandomReal-1;
            Inc(J);
        end;
        Inc(I);
    end;
    MTSize := 1+RandomInteger(Max(M, N));
    MakeACopyOldMem(X, MTSize, M, R);
    InternalMatrixMatrixMultiply(R, 0, MTSize-1, 0, M-1, False, Q, 0, M-1, 0, M-1, False, R1, 0, MTSize-1, 0, M-1);
    MakeACopyOldMem(X, MTSize, M, R2);
    RMatrixBDMultiplyByQ(T, M, N, TauQ, R2, MTSize, M, True, False);
    MulErrors := MulErrors or AP_FP_Greater(MatrixDiff(R1, R2, MTSize, M),Threshold);
    MakeACopyOldMem(X, MTSize, M, R);
    InternalMatrixMatrixMultiply(R, 0, MTSize-1, 0, M-1, False, Q, 0, M-1, 0, M-1, True, R1, 0, MTSize-1, 0, M-1);
    MakeACopyOldMem(X, MTSize, M, R2);
    RMatrixBDMultiplyByQ(T, M, N, TauQ, R2, MTSize, M, True, True);
    MulErrors := MulErrors or AP_FP_Greater(MatrixDiff(R1, R2, MTSize, M),Threshold);
    MakeACopyOldMem(X, M, MTSize, R);
    InternalMatrixMatrixMultiply(Q, 0, M-1, 0, M-1, False, R, 0, M-1, 0, MTSize-1, False, R1, 0, M-1, 0, MTSize-1);
    MakeACopyOldMem(X, M, MTSize, R2);
    RMatrixBDMultiplyByQ(T, M, N, TauQ, R2, M, MTSize, False, False);
    MulErrors := MulErrors or AP_FP_Greater(MatrixDiff(R1, R2, M, MTSize),Threshold);
    MakeACopyOldMem(X, M, MTSize, R);
    InternalMatrixMatrixMultiply(Q, 0, M-1, 0, M-1, True, R, 0, M-1, 0, MTSize-1, False, R1, 0, M-1, 0, MTSize-1);
    MakeACopyOldMem(X, M, MTSize, R2);
    RMatrixBDMultiplyByQ(T, M, N, TauQ, R2, M, MTSize, False, True);
    MulErrors := MulErrors or AP_FP_Greater(MatrixDiff(R1, R2, M, MTSize),Threshold);
    MakeACopyOldMem(X, MTSize, N, R);
    InternalMatrixMatrixMultiply(R, 0, MTSize-1, 0, N-1, False, PT, 0, N-1, 0, N-1, True, R1, 0, MTSize-1, 0, N-1);
    MakeACopyOldMem(X, MTSize, N, R2);
    RMatrixBDMultiplyByP(T, M, N, TauP, R2, MTSize, N, True, False);
    MulErrors := MulErrors or AP_FP_Greater(MatrixDiff(R1, R2, MTSize, N),Threshold);
    MakeACopyOldMem(X, MTSize, N, R);
    InternalMatrixMatrixMultiply(R, 0, MTSize-1, 0, N-1, False, PT, 0, N-1, 0, N-1, False, R1, 0, MTSize-1, 0, N-1);
    MakeACopyOldMem(X, MTSize, N, R2);
    RMatrixBDMultiplyByP(T, M, N, TauP, R2, MTSize, N, True, True);
    MulErrors := MulErrors or AP_FP_Greater(MatrixDiff(R1, R2, MTSize, N),Threshold);
    MakeACopyOldMem(X, N, MTSize, R);
    InternalMatrixMatrixMultiply(PT, 0, N-1, 0, N-1, True, R, 0, N-1, 0, MTSize-1, False, R1, 0, N-1, 0, MTSize-1);
    MakeACopyOldMem(X, N, MTSize, R2);
    RMatrixBDMultiplyByP(T, M, N, TauP, R2, N, MTSize, False, False);
    MulErrors := MulErrors or AP_FP_Greater(MatrixDiff(R1, R2, N, MTSize),Threshold);
    MakeACopyOldMem(X, N, MTSize, R);
    InternalMatrixMatrixMultiply(PT, 0, N-1, 0, N-1, False, R, 0, N-1, 0, MTSize-1, False, R1, 0, N-1, 0, MTSize-1);
    MakeACopyOldMem(X, N, MTSize, R2);
    RMatrixBDMultiplyByP(T, M, N, TauP, R2, N, MTSize, False, True);
    MulErrors := MulErrors or AP_FP_Greater(MatrixDiff(R1, R2, N, MTSize),Threshold);
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
function testbdunit_test_silent():Boolean;
begin
    Result := TestBD(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testbdunit_test():Boolean;
begin
    Result := TestBD(False);
end;


end.