unit testsvdunit;
interface
uses Math, Sysutils, Ap, hblas, reflections, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, bdsvd, svd;

function TestSVD(Silent : Boolean):Boolean;
function testsvdunit_test_silent():Boolean;
function testsvdunit_test():Boolean;

implementation

var
    FailCount : AlglibInteger;
    SuccCount : AlglibInteger;

procedure FillSparseA(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Sparcity : Double);forward;
procedure GetSVDError(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const U : TReal2DArray;
     const W : TReal1DArray;
     const VT : TReal2DArray;
     var MatErr : Double;
     var OrtErr : Double;
     var WSorted : Boolean);forward;
procedure TestSVDProblem(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var MatErr : Double;
     var OrtErr : Double;
     var OtherErr : Double;
     var WSorted : Boolean;
     var WFailed : Boolean);forward;


(*************************************************************************
Testing SVD decomposition subroutine
*************************************************************************)
function TestSVD(Silent : Boolean):Boolean;
var
    A : TReal2DArray;
    M : AlglibInteger;
    N : AlglibInteger;
    MaxMN : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    GPass : AlglibInteger;
    Pass : AlglibInteger;
    WasErrors : Boolean;
    WSorted : Boolean;
    WFailed : Boolean;
    MatErr : Double;
    OrtErr : Double;
    OtherErr : Double;
    Threshold : Double;
    FailThreshold : Double;
    FailR : Double;
begin
    FailCount := 0;
    SuccCount := 0;
    MatErr := 0;
    OrtErr := 0;
    OtherErr := 0;
    WSorted := True;
    WFailed := False;
    WasErrors := False;
    MaxMN := 30;
    Threshold := 5*100*MachineEpsilon;
    FailThreshold := 5.0E-3;
    SetLength(A, MaxMN-1+1, MaxMN-1+1);
    
    //
    // TODO: div by zero fail, convergence fail
    //
    GPass:=1;
    while GPass<=1 do
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
        while I<=Min(5, MaxMN) do
        begin
            J:=1;
            while J<=Min(5, MaxMN) do
            begin
                TestSVDProblem(A, I, J, MatErr, OrtErr, OtherErr, WSorted, WFailed);
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
            while J<=Min(5, MaxMN)-1 do
            begin
                A[I,J] := 2*RandomReal-1;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=1;
        while I<=MaxMN do
        begin
            J:=1;
            while J<=Min(5, MaxMN) do
            begin
                TestSVDProblem(A, I, J, MatErr, OrtErr, OtherErr, WSorted, WFailed);
                Inc(J);
            end;
            Inc(I);
        end;
        I:=0;
        while I<=Min(5, MaxMN)-1 do
        begin
            J:=0;
            while J<=MaxMN-1 do
            begin
                A[I,J] := 2*RandomReal-1;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=1;
        while I<=Min(5, MaxMN) do
        begin
            J:=1;
            while J<=MaxMN do
            begin
                TestSVDProblem(A, I, J, MatErr, OrtErr, OtherErr, WSorted, WFailed);
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Dense matrices
        //
        M:=1;
        while M<=Min(10, MaxMN) do
        begin
            N:=1;
            while N<=Min(10, MaxMN) do
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
                TestSVDProblem(A, M, N, MatErr, OrtErr, OtherErr, WSorted, WFailed);
                Inc(N);
            end;
            Inc(M);
        end;
        
        //
        // Sparse matrices, very sparse matrices, incredible sparse matrices
        //
        M:=1;
        while M<=10 do
        begin
            N:=1;
            while N<=10 do
            begin
                Pass:=1;
                while Pass<=2 do
                begin
                    FillSparseA(A, M, N, 0.8);
                    TestSVDProblem(A, M, N, MatErr, OrtErr, OtherErr, WSorted, WFailed);
                    FillSparseA(A, M, N, 0.9);
                    TestSVDProblem(A, M, N, MatErr, OrtErr, OtherErr, WSorted, WFailed);
                    FillSparseA(A, M, N, 0.95);
                    TestSVDProblem(A, M, N, MatErr, OrtErr, OtherErr, WSorted, WFailed);
                    Inc(Pass);
                end;
                Inc(N);
            end;
            Inc(M);
        end;
        Inc(GPass);
    end;
    
    //
    // report
    //
    FailR := AP_Double(FailCount)/(SuccCount+FailCount);
    WasErrors := AP_FP_Greater(MatErr,Threshold) or AP_FP_Greater(OrtErr,Threshold) or AP_FP_Greater(OtherErr,Threshold) or  not WSorted or AP_FP_Greater(FailR,FailThreshold);
    if  not Silent then
    begin
        Write(Format('TESTING SVD DECOMPOSITION'#13#10'',[]));
        Write(Format('SVD decomposition error:                 %5.4e'#13#10'',[
            MatErr]));
        Write(Format('SVD orthogonality error:                 %5.4e'#13#10'',[
            OrtErr]));
        Write(Format('SVD with different parameters error:     %5.4e'#13#10'',[
            OtherErr]));
        Write(Format('Singular values order:                   ',[]));
        if WSorted then
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


procedure GetSVDError(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     const U : TReal2DArray;
     const W : TReal1DArray;
     const VT : TReal2DArray;
     var MatErr : Double;
     var OrtErr : Double;
     var WSorted : Boolean);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    MinMN : AlglibInteger;
    LocErr : Double;
    SM : Double;
    i_ : AlglibInteger;
begin
    MinMN := Min(M, N);
    
    //
    // decomposition error
    //
    LocErr := 0;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            SM := 0;
            K:=0;
            while K<=MinMN-1 do
            begin
                SM := SM+W[K]*U[I,K]*VT[K,J];
                Inc(K);
            end;
            LocErr := Max(LocErr, AbsReal(A[I,J]-SM));
            Inc(J);
        end;
        Inc(I);
    end;
    MatErr := Max(MatErr, LocErr);
    
    //
    // orthogonality error
    //
    LocErr := 0;
    I:=0;
    while I<=MinMN-1 do
    begin
        J:=I;
        while J<=MinMN-1 do
        begin
            SM := 0.0;
            for i_ := 0 to M-1 do
            begin
                SM := SM + U[i_,I]*U[i_,J];
            end;
            if I<>J then
            begin
                LocErr := Max(LocErr, AbsReal(SM));
            end
            else
            begin
                LocErr := Max(LocErr, AbsReal(SM-1));
            end;
            SM := APVDotProduct(@VT[I][0], 0, N-1, @VT[J][0], 0, N-1);
            if I<>J then
            begin
                LocErr := Max(LocErr, AbsReal(SM));
            end
            else
            begin
                LocErr := Max(LocErr, AbsReal(SM-1));
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    OrtErr := Max(OrtErr, LocErr);
    
    //
    // values order error
    //
    I:=1;
    while I<=MinMN-1 do
    begin
        if AP_FP_Greater(W[I],W[I-1]) then
        begin
            WSorted := False;
        end;
        Inc(I);
    end;
end;


procedure TestSVDProblem(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var MatErr : Double;
     var OrtErr : Double;
     var OtherErr : Double;
     var WSorted : Boolean;
     var WFailed : Boolean);
var
    U : TReal2DArray;
    VT : TReal2DArray;
    U2 : TReal2DArray;
    VT2 : TReal2DArray;
    W : TReal1DArray;
    W2 : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    UJob : AlglibInteger;
    VTJob : AlglibInteger;
    MemJob : AlglibInteger;
    UCheck : AlglibInteger;
    VTCheck : AlglibInteger;
    V : Double;
    MX : Double;
begin
    
    //
    // Main SVD test
    //
    if  not RMatrixSVD(A, M, N, 2, 2, 2, W, U, VT) then
    begin
        FailCount := FailCount+1;
        WFailed := True;
        Exit;
    end;
    GetSVDError(A, M, N, U, W, VT, MatErr, OrtErr, WSorted);
    
    //
    // Additional SVD tests
    //
    UJob:=0;
    while UJob<=2 do
    begin
        VTJob:=0;
        while VTJob<=2 do
        begin
            MemJob:=0;
            while MemJob<=2 do
            begin
                if  not RMatrixSVD(A, M, N, UJob, VTJob, MemJob, W2, U2, VT2) then
                begin
                    FailCount := FailCount+1;
                    WFailed := True;
                    Exit;
                end;
                UCheck := 0;
                if UJob=1 then
                begin
                    UCheck := Min(M, N);
                end;
                if UJob=2 then
                begin
                    UCheck := M;
                end;
                VTCheck := 0;
                if VTJob=1 then
                begin
                    VTCheck := Min(M, N);
                end;
                if VTJob=2 then
                begin
                    VTCheck := N;
                end;
                I:=0;
                while I<=M-1 do
                begin
                    J:=0;
                    while J<=UCheck-1 do
                    begin
                        OtherErr := Max(OtherErr, AbsReal(U[I,J]-U2[I,J]));
                        Inc(J);
                    end;
                    Inc(I);
                end;
                I:=0;
                while I<=VTCheck-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        OtherErr := Max(OtherErr, AbsReal(VT[I,J]-VT2[I,J]));
                        Inc(J);
                    end;
                    Inc(I);
                end;
                I:=0;
                while I<=Min(M, N)-1 do
                begin
                    OtherErr := Max(OtherErr, AbsReal(W[I]-W2[I]));
                    Inc(I);
                end;
                Inc(MemJob);
            end;
            Inc(VTJob);
        end;
        Inc(UJob);
    end;
    
    //
    // update counter
    //
    SuccCount := SuccCount+1;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testsvdunit_test_silent():Boolean;
begin
    Result := TestSVD(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testsvdunit_test():Boolean;
begin
    Result := TestSVD(False);
end;


end.