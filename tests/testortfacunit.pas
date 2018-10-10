unit testortfacunit;
interface
uses Math, Sysutils, Ap, hblas, reflections, creflections, sblas, ablasf, ablas, ortfac;

function TestOrtFac(Silent : Boolean):Boolean;
function testortfacunit_test_silent():Boolean;
function testortfacunit_test():Boolean;

implementation

function RMatrixDiff(const A : TReal2DArray;
     const B : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger):Double;forward;
procedure RMatrixMakeACopy(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);forward;
procedure CMatrixMakeACopy(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TComplex2DArray);forward;
procedure RMatrixFillSparseA(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Sparcity : Double);forward;
procedure CMatrixFillSparseA(var A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Sparcity : Double);forward;
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
procedure TestRQRProblem(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Threshold : Double;
     var QRErrors : Boolean);forward;
procedure TestCQRProblem(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Threshold : Double;
     var QRErrors : Boolean);forward;
procedure TestRLQProblem(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Threshold : Double;
     var LQErrors : Boolean);forward;
procedure TestCLQProblem(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Threshold : Double;
     var LQErrors : Boolean);forward;
procedure TestRBDProblem(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Threshold : Double;
     var BDErrors : Boolean);forward;
procedure TestRHessProblem(const A : TReal2DArray;
     N : AlglibInteger;
     Threshold : Double;
     var HessErrors : Boolean);forward;
procedure TestRTDProblem(const A : TReal2DArray;
     N : AlglibInteger;
     Threshold : Double;
     var TDErrors : Boolean);forward;
procedure TestCTDProblem(const A : TComplex2DArray;
     N : AlglibInteger;
     Threshold : Double;
     var TDErrors : Boolean);forward;


(*************************************************************************
Main unittest subroutine
*************************************************************************)
function TestOrtFac(Silent : Boolean):Boolean;
var
    MaxMN : AlglibInteger;
    Threshold : Double;
    PassCount : AlglibInteger;
    MX : AlglibInteger;
    RA : TReal2DArray;
    CA : TComplex2DArray;
    M : AlglibInteger;
    N : AlglibInteger;
    Pass : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    RQRErrors : Boolean;
    RLQErrors : Boolean;
    CQRErrors : Boolean;
    CLQErrors : Boolean;
    RBDErrors : Boolean;
    RHessErrors : Boolean;
    RTDErrors : Boolean;
    CTDErrors : Boolean;
    WasErrors : Boolean;
begin
    WasErrors := False;
    RQRErrors := False;
    RLQErrors := False;
    CQRErrors := False;
    CLQErrors := False;
    RBDErrors := False;
    RHessErrors := False;
    RTDErrors := False;
    CTDErrors := False;
    MaxMN := 3*ABLASBlockSize(RA)+1;
    PassCount := 1;
    Threshold := 5*1000*MachineEpsilon;
    
    //
    // Different problems
    //
    MX:=1;
    while MX<=MaxMN do
    begin
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // Rectangular factorizations: QR, LQ, bidiagonal
            // Matrix types: zero, dense, sparse
            //
            N := 1+RandomInteger(MX);
            M := 1+RandomInteger(MX);
            if AP_FP_Greater(RandomReal,0.5) then
            begin
                N := MX;
            end
            else
            begin
                M := MX;
            end;
            SetLength(RA, M, N);
            SetLength(CA, M, N);
            I:=0;
            while I<=M-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    RA[I,J] := 0;
                    CA[I,J] := C_Complex(0);
                    Inc(J);
                end;
                Inc(I);
            end;
            TestRQRProblem(RA, M, N, Threshold, RQRErrors);
            TestRLQProblem(RA, M, N, Threshold, RLQErrors);
            TestCQRProblem(CA, M, N, Threshold, CQRErrors);
            TestCLQProblem(CA, M, N, Threshold, CLQErrors);
            TestRBDProblem(RA, M, N, Threshold, RBDErrors);
            I:=0;
            while I<=M-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    RA[I,J] := 2*RandomReal-1;
                    CA[I,J].X := 2*RandomReal-1;
                    CA[I,J].Y := 2*RandomReal-1;
                    Inc(J);
                end;
                Inc(I);
            end;
            TestRQRProblem(RA, M, N, Threshold, RQRErrors);
            TestRLQProblem(RA, M, N, Threshold, RLQErrors);
            TestCQRProblem(CA, M, N, Threshold, CQRErrors);
            TestCLQProblem(CA, M, N, Threshold, CLQErrors);
            TestRBDProblem(RA, M, N, Threshold, RBDErrors);
            RMatrixFillSparseA(RA, M, N, 0.95);
            CMatrixFillSparseA(CA, M, N, 0.95);
            TestRQRProblem(RA, M, N, Threshold, RQRErrors);
            TestRLQProblem(RA, M, N, Threshold, RLQErrors);
            TestCQRProblem(CA, M, N, Threshold, CQRErrors);
            TestCLQProblem(CA, M, N, Threshold, CLQErrors);
            TestRBDProblem(RA, M, N, Threshold, RBDErrors);
            
            //
            // Square factorizations: Hessenberg, tridiagonal
            // Matrix types: zero, dense, sparse
            //
            SetLength(RA, MX, MX);
            SetLength(CA, MX, MX);
            I:=0;
            while I<=MX-1 do
            begin
                J:=0;
                while J<=MX-1 do
                begin
                    RA[I,J] := 0;
                    CA[I,J] := C_Complex(0);
                    Inc(J);
                end;
                Inc(I);
            end;
            TestRHessProblem(RA, MX, Threshold, RHessErrors);
            I:=0;
            while I<=MX-1 do
            begin
                J:=0;
                while J<=MX-1 do
                begin
                    RA[I,J] := 2*RandomReal-1;
                    CA[I,J].X := 2*RandomReal-1;
                    CA[I,J].Y := 2*RandomReal-1;
                    Inc(J);
                end;
                Inc(I);
            end;
            TestRHessProblem(RA, MX, Threshold, RHessErrors);
            RMatrixFillSparseA(RA, MX, MX, 0.95);
            CMatrixFillSparseA(CA, MX, MX, 0.95);
            TestRHessProblem(RA, MX, Threshold, RHessErrors);
            
            //
            // Symetric factorizations: tridiagonal
            // Matrix types: zero, dense, sparse
            //
            SetLength(RA, MX, MX);
            SetLength(CA, MX, MX);
            I:=0;
            while I<=MX-1 do
            begin
                J:=0;
                while J<=MX-1 do
                begin
                    RA[I,J] := 0;
                    CA[I,J] := C_Complex(0);
                    Inc(J);
                end;
                Inc(I);
            end;
            TestRTDProblem(RA, MX, Threshold, RTDErrors);
            TestCTDProblem(CA, MX, Threshold, CTDErrors);
            I:=0;
            while I<=MX-1 do
            begin
                J:=I;
                while J<=MX-1 do
                begin
                    RA[I,J] := 2*RandomReal-1;
                    CA[I,J].X := 2*RandomReal-1;
                    CA[I,J].Y := 2*RandomReal-1;
                    RA[J,I] := RA[I,J];
                    CA[J,I] := Conj(CA[I,J]);
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=MX-1 do
            begin
                CA[I,I] := C_Complex(2*RandomReal-1);
                Inc(I);
            end;
            TestRTDProblem(RA, MX, Threshold, RTDErrors);
            TestCTDProblem(CA, MX, Threshold, CTDErrors);
            RMatrixFillSparseA(RA, MX, MX, 0.95);
            CMatrixFillSparseA(CA, MX, MX, 0.95);
            I:=0;
            while I<=MX-1 do
            begin
                J:=I;
                while J<=MX-1 do
                begin
                    RA[J,I] := RA[I,J];
                    CA[J,I] := Conj(CA[I,J]);
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=MX-1 do
            begin
                CA[I,I] := C_Complex(2*RandomReal-1);
                Inc(I);
            end;
            TestRTDProblem(RA, MX, Threshold, RTDErrors);
            TestCTDProblem(CA, MX, Threshold, CTDErrors);
            Inc(Pass);
        end;
        Inc(MX);
    end;
    
    //
    // report
    //
    WasErrors := RQRErrors or RLQErrors or CQRErrors or CLQErrors or RBDErrors or RHessErrors or RTDErrors or CTDErrors;
    if  not Silent then
    begin
        Write(Format('TESTING ORTFAC UNIT'#13#10'',[]));
        Write(Format('RQR ERRORS:                              ',[]));
        if  not RQRErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('RLQ ERRORS:                              ',[]));
        if  not RLQErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('CQR ERRORS:                              ',[]));
        if  not CQRErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('CLQ ERRORS:                              ',[]));
        if  not CLQErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('RBD ERRORS:                              ',[]));
        if  not RBDErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('RHESS ERRORS:                            ',[]));
        if  not RHessErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('RTD ERRORS:                              ',[]));
        if  not RTDErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('CTD ERRORS:                              ',[]));
        if  not CTDErrors then
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
Diff
*************************************************************************)
function RMatrixDiff(const A : TReal2DArray;
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
Copy
*************************************************************************)
procedure RMatrixMakeACopy(const A : TReal2DArray;
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
procedure CMatrixMakeACopy(const A : TComplex2DArray;
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
Sparse fill
*************************************************************************)
procedure RMatrixFillSparseA(var A : TReal2DArray;
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
Sparse fill
*************************************************************************)
procedure CMatrixFillSparseA(var A : TComplex2DArray;
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
Problem testing
*************************************************************************)
procedure TestRQRProblem(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Threshold : Double;
     var QRErrors : Boolean);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    B : TReal2DArray;
    TauB : TReal1DArray;
    Q : TReal2DArray;
    R : TReal2DArray;
    Q2 : TReal2DArray;
    V : Double;
    i_ : AlglibInteger;
begin
    
    //
    // Test decompose-and-unpack error
    //
    RMatrixMakeACopy(A, M, N, B);
    RMatrixQR(B, M, N, TauB);
    RMatrixQRUnpackQ(B, M, N, TauB, M, Q);
    RMatrixQRUnpackR(B, M, N, R);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := 0.0;
            for i_ := 0 to M-1 do
            begin
                V := V + Q[I,i_]*R[i_,J];
            end;
            QRErrors := QRErrors or AP_FP_Greater(AbsReal(V-A[I,J]),Threshold);
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
            QRErrors := QRErrors or AP_FP_Neq(R[I,J],0);
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
            V := APVDotProduct(@Q[I][0], 0, M-1, @Q[J][0], 0, M-1);
            if I=J then
            begin
                V := V-1;
            end;
            QRErrors := QRErrors or AP_FP_Greater_Eq(AbsReal(V),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Test for other errors
    //
    K := 1+RandomInteger(M);
    RMatrixQRUnpackQ(B, M, N, TauB, K, Q2);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=K-1 do
        begin
            QRErrors := QRErrors or AP_FP_Greater(AbsReal(Q2[I,J]-Q[I,J]),10*MachineEpsilon);
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Problem testing
*************************************************************************)
procedure TestCQRProblem(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Threshold : Double;
     var QRErrors : Boolean);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    B : TComplex2DArray;
    TauB : TComplex1DArray;
    Q : TComplex2DArray;
    R : TComplex2DArray;
    Q2 : TComplex2DArray;
    V : Complex;
    i_ : AlglibInteger;
begin
    
    //
    // Test decompose-and-unpack error
    //
    CMatrixMakeACopy(A, M, N, B);
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
            QRErrors := QRErrors or AP_FP_Greater(AbsComplex(C_Sub(V,A[I,J])),Threshold);
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
            QRErrors := QRErrors or C_NotEqualR(R[I,J],0);
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
                V := C_SubR(V,1);
            end;
            QRErrors := QRErrors or AP_FP_Greater_Eq(AbsComplex(V),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Test for other errors
    //
    K := 1+RandomInteger(M);
    CMatrixQRUnpackQ(B, M, N, TauB, K, Q2);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=K-1 do
        begin
            QRErrors := QRErrors or AP_FP_Greater(AbsComplex(C_Sub(Q2[I,J],Q[I,J])),10*MachineEpsilon);
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Problem testing
*************************************************************************)
procedure TestRLQProblem(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Threshold : Double;
     var LQErrors : Boolean);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    B : TReal2DArray;
    TauB : TReal1DArray;
    Q : TReal2DArray;
    L : TReal2DArray;
    Q2 : TReal2DArray;
    V : Double;
    i_ : AlglibInteger;
begin
    
    //
    // Test decompose-and-unpack error
    //
    RMatrixMakeACopy(A, M, N, B);
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
            LQErrors := LQErrors or AP_FP_Greater_Eq(AbsReal(V-A[I,J]),Threshold);
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
            LQErrors := LQErrors or AP_FP_Neq(L[I,J],0);
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
                V := V-1;
            end;
            LQErrors := LQErrors or AP_FP_Greater_Eq(AbsReal(V),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Test for other errors
    //
    K := 1+RandomInteger(N);
    RMatrixLQUnpackQ(B, M, N, TauB, K, Q2);
    I:=0;
    while I<=K-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            LQErrors := LQErrors or AP_FP_Greater(AbsReal(Q2[I,J]-Q[I,J]),10*MachineEpsilon);
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Problem testing
*************************************************************************)
procedure TestCLQProblem(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Threshold : Double;
     var LQErrors : Boolean);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    B : TComplex2DArray;
    TauB : TComplex1DArray;
    Q : TComplex2DArray;
    L : TComplex2DArray;
    Q2 : TComplex2DArray;
    V : Complex;
    i_ : AlglibInteger;
begin
    
    //
    // Test decompose-and-unpack error
    //
    CMatrixMakeACopy(A, M, N, B);
    CMatrixLQ(B, M, N, TauB);
    CMatrixLQUnpackQ(B, M, N, TauB, N, Q);
    CMatrixLQUnpackL(B, M, N, L);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                V := C_Add(V,C_Mul(L[I,i_],Q[i_,J]));
            end;
            LQErrors := LQErrors or AP_FP_Greater_Eq(AbsComplex(C_Sub(V,A[I,J])),Threshold);
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
            LQErrors := LQErrors or C_NotEqualR(L[I,J],0);
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
            V := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                V := C_Add(V,C_Mul(Q[I,i_],Conj(Q[J,i_])));
            end;
            if I=J then
            begin
                V := C_SubR(V,1);
            end;
            LQErrors := LQErrors or AP_FP_Greater_Eq(AbsComplex(V),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Test for other errors
    //
    K := 1+RandomInteger(N);
    CMatrixLQUnpackQ(B, M, N, TauB, K, Q2);
    I:=0;
    while I<=K-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            LQErrors := LQErrors or AP_FP_Greater(AbsComplex(C_Sub(Q2[I,J],Q[I,J])),10*MachineEpsilon);
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Problem testing
*************************************************************************)
procedure TestRBDProblem(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Threshold : Double;
     var BDErrors : Boolean);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
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
    // Bidiagonal decomposition error
    //
    RMatrixMakeACopy(A, M, N, T);
    RMatrixBD(T, M, N, TauQ, TauP);
    RMatrixBDUnpackQ(T, M, N, TauQ, M, Q);
    RMatrixBDUnpackPT(T, M, N, TauP, N, PT);
    RMatrixBDUnpackDiagonals(T, M, N, Up, D, E);
    SetLength(BD, M, N);
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
    SetLength(R, M, N);
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
            BDErrors := BDErrors or AP_FP_Greater(AbsReal(V-A[I,J]),Threshold);
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
                BDErrors := BDErrors or AP_FP_Greater(AbsReal(V-1),Threshold);
            end
            else
            begin
                BDErrors := BDErrors or AP_FP_Greater(AbsReal(V),Threshold);
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
                BDErrors := BDErrors or AP_FP_Greater(AbsReal(V-1),Threshold);
            end
            else
            begin
                BDErrors := BDErrors or AP_FP_Greater(AbsReal(V),Threshold);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Partial unpacking test
    //
    K := 1+RandomInteger(M);
    RMatrixBDUnpackQ(T, M, N, TauQ, K, R);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=K-1 do
        begin
            BDErrors := BDErrors or AP_FP_Greater(AbsReal(R[I,J]-Q[I,J]),10*MachineEpsilon);
            Inc(J);
        end;
        Inc(I);
    end;
    K := 1+RandomInteger(N);
    RMatrixBDUnpackPT(T, M, N, TauP, K, R);
    I:=0;
    while I<=K-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            BDErrors := BDErrors or AP_FP_Neq(R[I,J]-PT[I,J],0);
            Inc(J);
        end;
        Inc(I);
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
    RMatrixMakeACopy(X, MTSize, M, R);
    InternalMatrixMatrixMultiply(R, 0, MTSize-1, 0, M-1, False, Q, 0, M-1, 0, M-1, False, R1, 0, MTSize-1, 0, M-1);
    RMatrixMakeACopy(X, MTSize, M, R2);
    RMatrixBDMultiplyByQ(T, M, N, TauQ, R2, MTSize, M, True, False);
    BDErrors := BDErrors or AP_FP_Greater(RMatrixDiff(R1, R2, MTSize, M),Threshold);
    RMatrixMakeACopy(X, MTSize, M, R);
    InternalMatrixMatrixMultiply(R, 0, MTSize-1, 0, M-1, False, Q, 0, M-1, 0, M-1, True, R1, 0, MTSize-1, 0, M-1);
    RMatrixMakeACopy(X, MTSize, M, R2);
    RMatrixBDMultiplyByQ(T, M, N, TauQ, R2, MTSize, M, True, True);
    BDErrors := BDErrors or AP_FP_Greater(RMatrixDiff(R1, R2, MTSize, M),Threshold);
    RMatrixMakeACopy(X, M, MTSize, R);
    InternalMatrixMatrixMultiply(Q, 0, M-1, 0, M-1, False, R, 0, M-1, 0, MTSize-1, False, R1, 0, M-1, 0, MTSize-1);
    RMatrixMakeACopy(X, M, MTSize, R2);
    RMatrixBDMultiplyByQ(T, M, N, TauQ, R2, M, MTSize, False, False);
    BDErrors := BDErrors or AP_FP_Greater(RMatrixDiff(R1, R2, M, MTSize),Threshold);
    RMatrixMakeACopy(X, M, MTSize, R);
    InternalMatrixMatrixMultiply(Q, 0, M-1, 0, M-1, True, R, 0, M-1, 0, MTSize-1, False, R1, 0, M-1, 0, MTSize-1);
    RMatrixMakeACopy(X, M, MTSize, R2);
    RMatrixBDMultiplyByQ(T, M, N, TauQ, R2, M, MTSize, False, True);
    BDErrors := BDErrors or AP_FP_Greater(RMatrixDiff(R1, R2, M, MTSize),Threshold);
    RMatrixMakeACopy(X, MTSize, N, R);
    InternalMatrixMatrixMultiply(R, 0, MTSize-1, 0, N-1, False, PT, 0, N-1, 0, N-1, True, R1, 0, MTSize-1, 0, N-1);
    RMatrixMakeACopy(X, MTSize, N, R2);
    RMatrixBDMultiplyByP(T, M, N, TauP, R2, MTSize, N, True, False);
    BDErrors := BDErrors or AP_FP_Greater(RMatrixDiff(R1, R2, MTSize, N),Threshold);
    RMatrixMakeACopy(X, MTSize, N, R);
    InternalMatrixMatrixMultiply(R, 0, MTSize-1, 0, N-1, False, PT, 0, N-1, 0, N-1, False, R1, 0, MTSize-1, 0, N-1);
    RMatrixMakeACopy(X, MTSize, N, R2);
    RMatrixBDMultiplyByP(T, M, N, TauP, R2, MTSize, N, True, True);
    BDErrors := BDErrors or AP_FP_Greater(RMatrixDiff(R1, R2, MTSize, N),Threshold);
    RMatrixMakeACopy(X, N, MTSize, R);
    InternalMatrixMatrixMultiply(PT, 0, N-1, 0, N-1, True, R, 0, N-1, 0, MTSize-1, False, R1, 0, N-1, 0, MTSize-1);
    RMatrixMakeACopy(X, N, MTSize, R2);
    RMatrixBDMultiplyByP(T, M, N, TauP, R2, N, MTSize, False, False);
    BDErrors := BDErrors or AP_FP_Greater(RMatrixDiff(R1, R2, N, MTSize),Threshold);
    RMatrixMakeACopy(X, N, MTSize, R);
    InternalMatrixMatrixMultiply(PT, 0, N-1, 0, N-1, False, R, 0, N-1, 0, MTSize-1, False, R1, 0, N-1, 0, MTSize-1);
    RMatrixMakeACopy(X, N, MTSize, R2);
    RMatrixBDMultiplyByP(T, M, N, TauP, R2, N, MTSize, False, True);
    BDErrors := BDErrors or AP_FP_Greater(RMatrixDiff(R1, R2, N, MTSize),Threshold);
end;


(*************************************************************************
Problem testing
*************************************************************************)
procedure TestRHessProblem(const A : TReal2DArray;
     N : AlglibInteger;
     Threshold : Double;
     var HessErrors : Boolean);
var
    B : TReal2DArray;
    H : TReal2DArray;
    Q : TReal2DArray;
    T1 : TReal2DArray;
    T2 : TReal2DArray;
    Tau : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
begin
    RMatrixMakeACopy(A, N, N, B);
    
    //
    // Decomposition
    //
    RMatrixHessenberg(B, N, Tau);
    RMatrixHessenbergUnpackQ(B, N, Tau, Q);
    RMatrixHessenbergUnpackH(B, N, H);
    
    //
    // Matrix properties
    //
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + Q[i_,I]*Q[i_,J];
            end;
            if I=J then
            begin
                V := V-1;
            end;
            HessErrors := HessErrors or AP_FP_Greater(AbsReal(V),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=I-2 do
        begin
            HessErrors := HessErrors or AP_FP_Neq(H[I,J],0);
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Decomposition error
    //
    SetLength(T1, N, N);
    SetLength(T2, N, N);
    InternalMatrixMatrixMultiply(Q, 0, N-1, 0, N-1, False, H, 0, N-1, 0, N-1, False, T1, 0, N-1, 0, N-1);
    InternalMatrixMatrixMultiply(T1, 0, N-1, 0, N-1, False, Q, 0, N-1, 0, N-1, True, T2, 0, N-1, 0, N-1);
    HessErrors := HessErrors or AP_FP_Greater(RMatrixDiff(T2, A, N, N),Threshold);
end;


(*************************************************************************
Tridiagonal tester
*************************************************************************)
procedure TestRTDProblem(const A : TReal2DArray;
     N : AlglibInteger;
     Threshold : Double;
     var TDErrors : Boolean);
var
    I : AlglibInteger;
    J : AlglibInteger;
    UA : TReal2DArray;
    LA : TReal2DArray;
    T : TReal2DArray;
    Q : TReal2DArray;
    T2 : TReal2DArray;
    T3 : TReal2DArray;
    Tau : TReal1DArray;
    D : TReal1DArray;
    E : TReal1DArray;
    V : Double;
    i_ : AlglibInteger;
begin
    SetLength(UA, N-1+1, N-1+1);
    SetLength(LA, N-1+1, N-1+1);
    SetLength(T, N-1+1, N-1+1);
    SetLength(Q, N-1+1, N-1+1);
    SetLength(T2, N-1+1, N-1+1);
    SetLength(T3, N-1+1, N-1+1);
    
    //
    // fill
    //
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            UA[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        J:=I;
        while J<=N-1 do
        begin
            UA[I,J] := A[I,J];
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
            LA[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=I do
        begin
            LA[I,J] := A[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Test 2tridiagonal: upper
    //
    SMatrixTD(UA, N, True, Tau, D, E);
    SMatrixTDUnpackQ(UA, N, True, Tau, Q);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            T[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        T[I,I] := D[I];
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        T[I,I+1] := E[I];
        T[I+1,I] := E[I];
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
                V := V + Q[i_,I]*A[i_,J];
            end;
            T2[I,J] := V;
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
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + T2[I,i_]*Q[i_,J];
            end;
            T3[I,J] := V;
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
            TDErrors := TDErrors or AP_FP_Greater(AbsReal(T3[I,J]-T[I,J]),Threshold);
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
                V := V-1;
            end;
            TDErrors := TDErrors or AP_FP_Greater(AbsReal(V),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Test 2tridiagonal: lower
    //
    SMatrixTD(LA, N, False, Tau, D, E);
    SMatrixTDUnpackQ(LA, N, False, Tau, Q);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            T[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        T[I,I] := D[I];
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        T[I,I+1] := E[I];
        T[I+1,I] := E[I];
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
                V := V + Q[i_,I]*A[i_,J];
            end;
            T2[I,J] := V;
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
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + T2[I,i_]*Q[i_,J];
            end;
            T3[I,J] := V;
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
            TDErrors := TDErrors or AP_FP_Greater(AbsReal(T3[I,J]-T[I,J]),Threshold);
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
                V := V-1;
            end;
            TDErrors := TDErrors or AP_FP_Greater(AbsReal(V),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Hermitian problem tester
*************************************************************************)
procedure TestCTDProblem(const A : TComplex2DArray;
     N : AlglibInteger;
     Threshold : Double;
     var TDErrors : Boolean);
var
    I : AlglibInteger;
    J : AlglibInteger;
    UA : TComplex2DArray;
    LA : TComplex2DArray;
    T : TComplex2DArray;
    Q : TComplex2DArray;
    T2 : TComplex2DArray;
    T3 : TComplex2DArray;
    Tau : TComplex1DArray;
    D : TReal1DArray;
    E : TReal1DArray;
    V : Complex;
    i_ : AlglibInteger;
begin
    SetLength(UA, N-1+1, N-1+1);
    SetLength(LA, N-1+1, N-1+1);
    SetLength(T, N-1+1, N-1+1);
    SetLength(Q, N-1+1, N-1+1);
    SetLength(T2, N-1+1, N-1+1);
    SetLength(T3, N-1+1, N-1+1);
    
    //
    // fill
    //
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            UA[I,J] := C_Complex(0);
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        J:=I;
        while J<=N-1 do
        begin
            UA[I,J] := A[I,J];
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
            LA[I,J] := C_Complex(0);
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=I do
        begin
            LA[I,J] := A[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Test 2tridiagonal: upper
    //
    HMatrixTD(UA, N, True, Tau, D, E);
    HMatrixTDUnpackQ(UA, N, True, Tau, Q);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            T[I,J] := C_Complex(0);
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        T[I,I] := C_Complex(D[I]);
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        T[I,I+1] := C_Complex(E[I]);
        T[I+1,I] := C_Complex(E[I]);
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                V := C_Add(V,C_Mul(Conj(Q[i_,I]),A[i_,J]));
            end;
            T2[I,J] := V;
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
            V := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                V := C_Add(V,C_Mul(T2[I,i_],Q[i_,J]));
            end;
            T3[I,J] := V;
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
            TDErrors := TDErrors or AP_FP_Greater(AbsComplex(C_Sub(T3[I,J],T[I,J])),Threshold);
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
            V := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                V := C_Add(V,C_Mul(Q[I,i_],Conj(Q[J,i_])));
            end;
            if I=J then
            begin
                V := C_SubR(V,1);
            end;
            TDErrors := TDErrors or AP_FP_Greater(AbsComplex(V),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Test 2tridiagonal: lower
    //
    HMatrixTD(LA, N, False, Tau, D, E);
    HMatrixTDUnpackQ(LA, N, False, Tau, Q);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            T[I,J] := C_Complex(0);
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        T[I,I] := C_Complex(D[I]);
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        T[I,I+1] := C_Complex(E[I]);
        T[I+1,I] := C_Complex(E[I]);
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                V := C_Add(V,C_Mul(Conj(Q[i_,I]),A[i_,J]));
            end;
            T2[I,J] := V;
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
            V := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                V := C_Add(V,C_Mul(T2[I,i_],Q[i_,J]));
            end;
            T3[I,J] := V;
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
            TDErrors := TDErrors or AP_FP_Greater(AbsComplex(C_Sub(T3[I,J],T[I,J])),Threshold);
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
            V := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                V := C_Add(V,C_Mul(Q[I,i_],Conj(Q[J,i_])));
            end;
            if I=J then
            begin
                V := C_SubR(V,1);
            end;
            TDErrors := TDErrors or AP_FP_Greater(AbsComplex(V),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testortfacunit_test_silent():Boolean;
begin
    Result := TestOrtFac(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testortfacunit_test():Boolean;
begin
    Result := TestOrtFac(False);
end;


end.