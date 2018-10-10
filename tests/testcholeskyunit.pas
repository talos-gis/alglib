unit testcholeskyunit;
interface
uses Math, Sysutils, Ap, cholesky;

function TestCholesky(Silent : Boolean):Boolean;
function testcholeskyunit_test_silent():Boolean;
function testcholeskyunit_test():Boolean;

implementation

function TestCholesky(Silent : Boolean):Boolean;
var
    L : TReal2DArray;
    A : TReal2DArray;
    N : AlglibInteger;
    Pass : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    MinIJ : AlglibInteger;
    UpperIn : Boolean;
    CR : Boolean;
    V : Double;
    Err : Double;
    WF : Boolean;
    WasErrors : Boolean;
    PassCount : AlglibInteger;
    MaxN : AlglibInteger;
    HTask : AlglibInteger;
    Threshold : Double;
    i_ : AlglibInteger;
begin
    Err := 0;
    WF := False;
    PassCount := 10;
    MaxN := 20;
    Threshold := 1000*MachineEpsilon;
    WasErrors := False;
    
    //
    // Test
    //
    N:=1;
    while N<=MaxN do
    begin
        SetLength(L, N-1+1, N-1+1);
        SetLength(A, N-1+1, N-1+1);
        HTask:=0;
        while HTask<=1 do
        begin
            Pass:=1;
            while Pass<=PassCount do
            begin
                UpperIn := HTask=0;
                
                //
                // Prepare task:
                // * A contains upper (or lower) half of SPD matrix
                // * L contains its Cholesky factor (upper or lower)
                //
                I:=0;
                while I<=N-1 do
                begin
                    J:=I+1;
                    while J<=N-1 do
                    begin
                        L[I,J] := 2*RandomReal-1;
                        L[J,I] := L[I,J];
                        Inc(J);
                    end;
                    L[I,I] := 0.8+RandomReal;
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
                        A[I,J] := V;
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
                        if UpperIn then
                        begin
                            if J<I then
                            begin
                                A[I,J] := 0;
                                L[I,J] := 0;
                            end;
                        end
                        else
                        begin
                            if I<J then
                            begin
                                A[I,J] := 0;
                                L[I,J] := 0;
                            end;
                        end;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                
                //
                // decomposition
                //
                if  not SPDMatrixCholesky(A, N, UpperIn) then
                begin
                    WF := True;
                    Inc(Pass);
                    Continue;
                end;
                
                //
                // Test
                //
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        Err := Max(Err, AbsReal(A[I,J]-L[I,J]));
                        Inc(J);
                    end;
                    Inc(I);
                end;
                Inc(Pass);
            end;
            Inc(HTask);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    WasErrors := AP_FP_Greater(Err,Threshold) or WF;
    if  not Silent then
    begin
        Write(Format('TESTING CHOLESKY DECOMPOSITION'#13#10'',[]));
        Write(Format('ERROR:                                   %5.4e'#13#10'',[
            Err]));
        Write(Format('ALWAYS SUCCEDED:                         ',[]));
        if WF then
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
Silent unit test
*************************************************************************)
function testcholeskyunit_test_silent():Boolean;
begin
    Result := TestCholesky(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testcholeskyunit_test():Boolean;
begin
    Result := TestCholesky(False);
end;


end.