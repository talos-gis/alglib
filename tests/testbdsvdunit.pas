unit testbdsvdunit;
interface
uses Math, Sysutils, Ap, rotations, bdsvd;

function TestBDSVD(Silent : Boolean):Boolean;
function testbdsvdunit_test_silent():Boolean;
function testbdsvdunit_test():Boolean;

implementation

var
    FailCount : AlglibInteger;
    SuccCount : AlglibInteger;

procedure FillIdentity(var A : TReal2DArray; N : AlglibInteger);forward;
procedure FillSparseDE(var D : TReal1DArray;
     var E : TReal1DArray;
     N : AlglibInteger;
     Sparcity : Double);forward;
procedure GetBDSVDError(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const U : TReal2DArray;
     const C : TReal2DArray;
     const W : TReal1DArray;
     const VT : TReal2DArray;
     var MatErr : Double;
     var OrtErr : Double;
     var WSorted : Boolean);forward;
procedure TestBDSVDProblem(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     var MatErr : Double;
     var OrtErr : Double;
     var WSorted : Boolean;
     var WFailed : Boolean);forward;


(*************************************************************************
Testing bidiagonal SVD decomposition subroutine
*************************************************************************)
function TestBDSVD(Silent : Boolean):Boolean;
var
    D : TReal1DArray;
    E : TReal1DArray;
    MEmpty : TReal2DArray;
    N : AlglibInteger;
    MaxN : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    GPass : AlglibInteger;
    Pass : AlglibInteger;
    WasErrors : Boolean;
    WSorted : Boolean;
    WFailed : Boolean;
    FailCase : Boolean;
    MatErr : Double;
    OrtErr : Double;
    Threshold : Double;
    FailThreshold : Double;
    FailR : Double;
begin
    FailCount := 0;
    SuccCount := 0;
    MatErr := 0;
    OrtErr := 0;
    WSorted := True;
    WFailed := False;
    WasErrors := False;
    MaxN := 15;
    Threshold := 5*100*MachineEpsilon;
    FailThreshold := 1.0E-2;
    SetLength(D, MaxN-1+1);
    SetLength(E, MaxN-2+1);
    
    //
    // special case: fail matrix
    //
    N := 5;
    d[0] := -8.27448347422711894000e-01;
    d[1] := -8.16705832087160854600e-01;
    d[2] := -2.53974358904729382800e-17;
    d[3] := -1.24626684881972815700e+00;
    d[4] := -4.64744131545637651000e-01;
    e[0] := -3.25785088656270038800e-01;
    e[1] := -1.03732413708914436580e-01;
    e[2] := -9.57365642262031357700e-02;
    e[3] := -2.71564153973817390400e-01;
    FailCase := RMatrixBDSVD(D, E, N, True, False, MEmpty, 0, MEmpty, 0, MEmpty, 0);
    
    //
    // special case: zero divide matrix
    // unfixed LAPACK routine should fail on this problem
    //
    N := 7;
    d[0] := -6.96462904751731892700e-01;
    d[1] := 0.00000000000000000000e+00;
    d[2] := -5.73827770385971991400e-01;
    d[3] := -6.62562624399371191700e-01;
    d[4] := 5.82737148001782223600e-01;
    d[5] := 3.84825263580925003300e-01;
    d[6] := 9.84087420830525472200e-01;
    e[0] := -7.30307931760612871800e-02;
    e[1] := -2.30079042939542843800e-01;
    e[2] := -6.87824621739351216300e-01;
    e[3] := -1.77306437707837570600e-02;
    e[4] := 1.78285126526551632000e-15;
    e[5] := -4.89434737751289969400e-02;
    RMatrixBDSVD(D, E, N, True, False, MEmpty, 0, MEmpty, 0, MEmpty, 0);
    
    //
    // zero matrix, several cases
    //
    I:=0;
    while I<=MaxN-1 do
    begin
        D[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=MaxN-2 do
    begin
        E[I] := 0;
        Inc(I);
    end;
    N:=1;
    while N<=MaxN do
    begin
        TestBDSVDProblem(D, E, N, MatErr, OrtErr, WSorted, WFailed);
        Inc(N);
    end;
    
    //
    // Dense matrix
    //
    N:=1;
    while N<=MaxN do
    begin
        Pass:=1;
        while Pass<=10 do
        begin
            I:=0;
            while I<=MaxN-1 do
            begin
                D[I] := 2*RandomReal-1;
                Inc(I);
            end;
            I:=0;
            while I<=MaxN-2 do
            begin
                E[I] := 2*RandomReal-1;
                Inc(I);
            end;
            TestBDSVDProblem(D, E, N, MatErr, OrtErr, WSorted, WFailed);
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // Sparse matrices, very sparse matrices, incredible sparse matrices
    //
    N:=1;
    while N<=MaxN do
    begin
        Pass:=1;
        while Pass<=10 do
        begin
            FillSparseDE(D, E, N, 0.5);
            TestBDSVDProblem(D, E, N, MatErr, OrtErr, WSorted, WFailed);
            FillSparseDE(D, E, N, 0.8);
            TestBDSVDProblem(D, E, N, MatErr, OrtErr, WSorted, WFailed);
            FillSparseDE(D, E, N, 0.9);
            TestBDSVDProblem(D, E, N, MatErr, OrtErr, WSorted, WFailed);
            FillSparseDE(D, E, N, 0.95);
            TestBDSVDProblem(D, E, N, MatErr, OrtErr, WSorted, WFailed);
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    FailR := AP_Double(FailCount)/(SuccCount+FailCount);
    WasErrors := AP_FP_Greater(MatErr,Threshold) or AP_FP_Greater(OrtErr,Threshold) or  not WSorted or AP_FP_Greater(FailR,FailThreshold);
    if  not Silent then
    begin
        Write(Format('TESTING BIDIAGONAL SVD DECOMPOSITION'#13#10'',[]));
        Write(Format('SVD decomposition error:                 %5.4e'#13#10'',[
            MatErr]));
        Write(Format('SVD orthogonality error:                 %5.4e'#13#10'',[
            OrtErr]));
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
        Write(Format('Fail matrix test:                        ',[]));
        if  not FailCase then
        begin
            Write(Format('AS EXPECTED'#13#10'',[]));
        end
        else
        begin
            Write(Format('CONVERGED (UNEXPECTED)'#13#10'',[]));
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


procedure FillIdentity(var A : TReal2DArray; N : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    SetLength(A, N-1+1, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if I=J then
            begin
                A[I,J] := 1;
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


procedure FillSparseDE(var D : TReal1DArray;
     var E : TReal1DArray;
     N : AlglibInteger;
     Sparcity : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    SetLength(D, N-1+1);
    SetLength(E, Max(0, N-2)+1);
    I:=0;
    while I<=N-1 do
    begin
        if AP_FP_Greater_Eq(RandomReal,Sparcity) then
        begin
            D[I] := 2*RandomReal-1;
        end
        else
        begin
            D[I] := 0;
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        if AP_FP_Greater_Eq(RandomReal,Sparcity) then
        begin
            E[I] := 2*RandomReal-1;
        end
        else
        begin
            E[I] := 0;
        end;
        Inc(I);
    end;
end;


procedure GetBDSVDError(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     const U : TReal2DArray;
     const C : TReal2DArray;
     const W : TReal1DArray;
     const VT : TReal2DArray;
     var MatErr : Double;
     var OrtErr : Double;
     var WSorted : Boolean);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    LocErr : Double;
    SM : Double;
    i_ : AlglibInteger;
begin
    
    //
    // decomposition error
    //
    LocErr := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            SM := 0;
            K:=0;
            while K<=N-1 do
            begin
                SM := SM+W[K]*U[I,K]*VT[K,J];
                Inc(K);
            end;
            if IsUpper then
            begin
                if I=J then
                begin
                    LocErr := Max(LocErr, AbsReal(D[I]-SM));
                end
                else
                begin
                    if I=J-1 then
                    begin
                        LocErr := Max(LocErr, AbsReal(E[I]-SM));
                    end
                    else
                    begin
                        LocErr := Max(LocErr, AbsReal(SM));
                    end;
                end;
            end
            else
            begin
                if I=J then
                begin
                    LocErr := Max(LocErr, AbsReal(D[I]-SM));
                end
                else
                begin
                    if I-1=J then
                    begin
                        LocErr := Max(LocErr, AbsReal(E[J]-SM));
                    end
                    else
                    begin
                        LocErr := Max(LocErr, AbsReal(SM));
                    end;
                end;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    MatErr := Max(MatErr, LocErr);
    
    //
    // check for C = U'
    // we consider it as decomposition error
    //
    LocErr := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            LocErr := Max(LocErr, AbsReal(U[I,J]-C[J,I]));
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
    while I<=N-1 do
    begin
        J:=I;
        while J<=N-1 do
        begin
            SM := 0.0;
            for i_ := 0 to N-1 do
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
    while I<=N-1 do
    begin
        if AP_FP_Greater(W[I],W[I-1]) then
        begin
            WSorted := False;
        end;
        Inc(I);
    end;
end;


procedure TestBDSVDProblem(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     var MatErr : Double;
     var OrtErr : Double;
     var WSorted : Boolean;
     var WFailed : Boolean);
var
    U : TReal2DArray;
    VT : TReal2DArray;
    C : TReal2DArray;
    W : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    MX : Double;
begin
    MX := 0;
    I:=0;
    while I<=N-1 do
    begin
        if AP_FP_Greater(AbsReal(D[I]),MX) then
        begin
            MX := AbsReal(D[I]);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        if AP_FP_Greater(AbsReal(E[I]),MX) then
        begin
            MX := AbsReal(E[I]);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(MX,0) then
    begin
        MX := 1;
    end;
    
    //
    // Upper BDSVD tests
    //
    SetLength(W, N-1+1);
    FillIdentity(U, N);
    FillIdentity(VT, N);
    FillIdentity(C, N);
    I:=0;
    while I<=N-1 do
    begin
        W[I] := D[I];
        Inc(I);
    end;
    if  not RMatrixBDSVD(W, E, N, True, False, U, N, C, N, VT, N) then
    begin
        FailCount := FailCount+1;
        WFailed := True;
        Exit;
    end;
    GetBDSVDError(D, E, N, True, U, C, W, VT, MatErr, OrtErr, WSorted);
    FillIdentity(U, N);
    FillIdentity(VT, N);
    FillIdentity(C, N);
    I:=0;
    while I<=N-1 do
    begin
        W[I] := D[I];
        Inc(I);
    end;
    if  not RMatrixBDSVD(W, E, N, True, True, U, N, C, N, VT, N) then
    begin
        FailCount := FailCount+1;
        WFailed := True;
        Exit;
    end;
    GetBDSVDError(D, E, N, True, U, C, W, VT, MatErr, OrtErr, WSorted);
    
    //
    // Lower BDSVD tests
    //
    SetLength(W, N-1+1);
    FillIdentity(U, N);
    FillIdentity(VT, N);
    FillIdentity(C, N);
    I:=0;
    while I<=N-1 do
    begin
        W[I] := D[I];
        Inc(I);
    end;
    if  not RMatrixBDSVD(W, E, N, False, False, U, N, C, N, VT, N) then
    begin
        FailCount := FailCount+1;
        WFailed := True;
        Exit;
    end;
    GetBDSVDError(D, E, N, False, U, C, W, VT, MatErr, OrtErr, WSorted);
    FillIdentity(U, N);
    FillIdentity(VT, N);
    FillIdentity(C, N);
    I:=0;
    while I<=N-1 do
    begin
        W[I] := D[I];
        Inc(I);
    end;
    if  not RMatrixBDSVD(W, E, N, False, True, U, N, C, N, VT, N) then
    begin
        FailCount := FailCount+1;
        WFailed := True;
        Exit;
    end;
    GetBDSVDError(D, E, N, False, U, C, W, VT, MatErr, OrtErr, WSorted);
    
    //
    // update counter
    //
    SuccCount := SuccCount+1;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testbdsvdunit_test_silent():Boolean;
begin
    Result := TestBDSVD(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testbdsvdunit_test():Boolean;
begin
    Result := TestBDSVD(False);
end;


end.