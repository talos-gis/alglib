unit testspdgevdunit;
interface
uses Math, Sysutils, Ap, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, sblas, blas, trlinsolve, safesolve, rcond, matinv, hblas, ortfac, rotations, hsschur, evd, spdgevd;

function TestSPDGEVD(Silent : Boolean):Boolean;
function testspdgevdunit_test_silent():Boolean;
function testspdgevdunit_test():Boolean;

implementation

(*************************************************************************
Testing bidiagonal SVD decomposition subroutine
*************************************************************************)
function TestSPDGEVD(Silent : Boolean):Boolean;
var
    Pass : AlglibInteger;
    N : AlglibInteger;
    PassCount : AlglibInteger;
    MaxN : AlglibInteger;
    ATask : AlglibInteger;
    BTask : AlglibInteger;
    D : TReal1DArray;
    T1 : TReal1DArray;
    A : TReal2DArray;
    B : TReal2DArray;
    AFull : TReal2DArray;
    BFull : TReal2DArray;
    L : TReal2DArray;
    Z : TReal2DArray;
    IsUpperA : Boolean;
    IsUpperB : Boolean;
    I : AlglibInteger;
    J : AlglibInteger;
    MinIJ : AlglibInteger;
    V : Double;
    V1 : Double;
    V2 : Double;
    CW : Boolean;
    Err : Double;
    ValErr : Double;
    Threshold : Double;
    WasErrors : Boolean;
    WFailed : Boolean;
    WNSorted : Boolean;
    i_ : AlglibInteger;
begin
    Threshold := 10000*MachineEpsilon;
    ValErr := 0;
    WFailed := False;
    WNSorted := False;
    MaxN := 20;
    PassCount := 5;
    
    //
    // Main cycle
    //
    N:=1;
    while N<=MaxN do
    begin
        Pass:=1;
        while Pass<=PassCount do
        begin
            ATask:=0;
            while ATask<=1 do
            begin
                BTask:=0;
                while BTask<=1 do
                begin
                    IsUpperA := ATask=0;
                    IsUpperB := BTask=0;
                    
                    //
                    // Initialize A, B, AFull, BFull
                    //
                    SetLength(T1, N-1+1);
                    SetLength(A, N-1+1, N-1+1);
                    SetLength(B, N-1+1, N-1+1);
                    SetLength(AFull, N-1+1, N-1+1);
                    SetLength(BFull, N-1+1, N-1+1);
                    SetLength(L, N-1+1, N-1+1);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A[I,J] := 2*RandomReal-1;
                            A[J,I] := A[I,J];
                            AFull[I,J] := A[I,J];
                            AFull[J,I] := A[I,J];
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=I+1;
                        while J<=N-1 do
                        begin
                            L[I,J] := RandomReal;
                            L[J,I] := L[I,J];
                            Inc(J);
                        end;
                        L[I,I] := 1.5+RandomReal;
                        Inc(I);
                    end;
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            MinIJ := Min(I, J);
                            V := 0.0;
                            for i_ := 0 to MinIJ do
                            begin
                                V := V + L[I,i_]*L[i_,J];
                            end;
                            B[I,J] := V;
                            B[J,I] := V;
                            BFull[I,J] := V;
                            BFull[J,I] := V;
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
                            if IsUpperA then
                            begin
                                if J<I then
                                begin
                                    A[I,J] := 2*RandomReal-1;
                                end;
                            end
                            else
                            begin
                                if I<J then
                                begin
                                    A[I,J] := 2*RandomReal-1;
                                end;
                            end;
                            if IsUpperB then
                            begin
                                if J<I then
                                begin
                                    B[I,J] := 2*RandomReal-1;
                                end;
                            end
                            else
                            begin
                                if I<J then
                                begin
                                    B[I,J] := 2*RandomReal-1;
                                end;
                            end;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    
                    //
                    // Problem 1
                    //
                    if  not SMatrixGEVD(A, N, IsUpperA, B, IsUpperB, 1, 1, D, Z) then
                    begin
                        WFailed := True;
                        Inc(BTask);
                        Continue;
                    end;
                    Err := 0;
                    J:=0;
                    while J<=N-1 do
                    begin
                        I:=0;
                        while I<=N-1 do
                        begin
                            V1 := 0.0;
                            for i_ := 0 to N-1 do
                            begin
                                V1 := V1 + AFull[I,i_]*Z[i_,J];
                            end;
                            V2 := 0.0;
                            for i_ := 0 to N-1 do
                            begin
                                V2 := V2 + BFull[I,i_]*Z[i_,J];
                            end;
                            Err := Max(Err, AbsReal(V1-D[J]*V2));
                            Inc(I);
                        end;
                        Inc(J);
                    end;
                    ValErr := Max(Err, ValErr);
                    
                    //
                    // Problem 2
                    //
                    if  not SMatrixGEVD(A, N, IsUpperA, B, IsUpperB, 1, 2, D, Z) then
                    begin
                        WFailed := True;
                        Inc(BTask);
                        Continue;
                    end;
                    Err := 0;
                    J:=0;
                    while J<=N-1 do
                    begin
                        I:=0;
                        while I<=N-1 do
                        begin
                            V1 := 0.0;
                            for i_ := 0 to N-1 do
                            begin
                                V1 := V1 + BFull[I,i_]*Z[i_,J];
                            end;
                            T1[I] := V1;
                            Inc(I);
                        end;
                        I:=0;
                        while I<=N-1 do
                        begin
                            V2 := APVDotProduct(@AFull[I][0], 0, N-1, @T1[0], 0, N-1);
                            Err := Max(Err, AbsReal(V2-D[J]*Z[I,J]));
                            Inc(I);
                        end;
                        Inc(J);
                    end;
                    ValErr := Max(Err, ValErr);
                    
                    //
                    // Test problem 3
                    //
                    if  not SMatrixGEVD(A, N, IsUpperA, B, IsUpperB, 1, 3, D, Z) then
                    begin
                        WFailed := True;
                        Inc(BTask);
                        Continue;
                    end;
                    Err := 0;
                    J:=0;
                    while J<=N-1 do
                    begin
                        I:=0;
                        while I<=N-1 do
                        begin
                            V1 := 0.0;
                            for i_ := 0 to N-1 do
                            begin
                                V1 := V1 + AFull[I,i_]*Z[i_,J];
                            end;
                            T1[I] := V1;
                            Inc(I);
                        end;
                        I:=0;
                        while I<=N-1 do
                        begin
                            V2 := APVDotProduct(@BFull[I][0], 0, N-1, @T1[0], 0, N-1);
                            Err := Max(Err, AbsReal(V2-D[J]*Z[I,J]));
                            Inc(I);
                        end;
                        Inc(J);
                    end;
                    ValErr := Max(Err, ValErr);
                    Inc(BTask);
                end;
                Inc(ATask);
            end;
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    WasErrors := AP_FP_Greater(ValErr,Threshold) or WFailed or WNSorted;
    if  not Silent then
    begin
        Write(Format('TESTING SYMMETRIC GEVD'#13#10'',[]));
        Write(Format('Av-lambdav error (generalized):          %5.4e'#13#10'',[
            ValErr]));
        Write(Format('Eigen values order:                      ',[]));
        if  not WNSorted then
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


(*************************************************************************
Silent unit test
*************************************************************************)
function testspdgevdunit_test_silent():Boolean;
begin
    Result := TestSPDGEVD(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testspdgevdunit_test():Boolean;
begin
    Result := TestSPDGEVD(False);
end;


end.