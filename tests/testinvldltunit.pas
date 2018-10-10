unit testinvldltunit;
interface
uses Math, Sysutils, Ap, sblas, ldlt, sinverse;

function TestInvLDLT(Silent : Boolean):Boolean;
function testinvldltunit_test_silent():Boolean;
function testinvldltunit_test():Boolean;

implementation

procedure GenerateMatrix(var A : TReal2DArray;
     N : AlglibInteger;
     Task : AlglibInteger);forward;
procedure RestoreMatrix(var A : TReal2DArray;
     N : AlglibInteger;
     UpperIn : Boolean);forward;


function TestInvLDLT(Silent : Boolean):Boolean;
var
    A : TReal2DArray;
    A2 : TReal2DArray;
    A3 : TReal2DArray;
    P : TInteger1DArray;
    N : AlglibInteger;
    Pass : AlglibInteger;
    MTask : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    MinIJ : AlglibInteger;
    UpperIn : Boolean;
    CR : Boolean;
    V : Double;
    Err : Double;
    WasErrors : Boolean;
    PassCount : AlglibInteger;
    MaxN : AlglibInteger;
    HTask : AlglibInteger;
    Threshold : Double;
    i_ : AlglibInteger;
begin
    Err := 0;
    PassCount := 10;
    MaxN := 20;
    Threshold := 10000000*MachineEpsilon;
    WasErrors := False;
    
    //
    // Test
    //
    N:=1;
    while N<=MaxN do
    begin
        SetLength(A, N-1+1, N-1+1);
        SetLength(A2, N-1+1, N-1+1);
        SetLength(A3, N-1+1, N-1+1);
        MTask:=2;
        while MTask<=2 do
        begin
            HTask:=0;
            while HTask<=1 do
            begin
                Pass:=1;
                while Pass<=PassCount do
                begin
                    UpperIn := HTask=0;
                    
                    //
                    // Prepare task:
                    // * A contains symmetric matrix
                    // * A2, A3 contains its upper (or lower) half
                    //
                    GenerateMatrix(A, N, MTask);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A2[I,J] := A[I,J];
                            A3[I,J] := A[I,J];
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
                                    A2[I,J] := 0;
                                    A3[I,J] := 0;
                                end;
                            end
                            else
                            begin
                                if I<J then
                                begin
                                    A2[I,J] := 0;
                                    A3[I,J] := 0;
                                end;
                            end;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    
                    //
                    // Test 1: inv(A2)
                    //
                    SMatrixInverse(A2, N, UpperIn);
                    RestoreMatrix(A2, N, UpperIn);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            V := 0.0;
                            for i_ := 0 to N-1 do
                            begin
                                V := V + A[I,i_]*A2[i_,J];
                            end;
                            if I=J then
                            begin
                                V := V-1;
                            end;
                            Err := Max(Err, AbsReal(V));
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    
                    //
                    // Test 2: inv(LDLt(A3))
                    //
                    SMatrixLDLT(A3, N, UpperIn, P);
                    SMatrixLDLTInverse(A3, P, N, UpperIn);
                    RestoreMatrix(A3, N, UpperIn);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            V := 0.0;
                            for i_ := 0 to N-1 do
                            begin
                                V := V + A[I,i_]*A3[i_,J];
                            end;
                            if I=J then
                            begin
                                V := V-1;
                            end;
                            Err := Max(Err, AbsReal(V));
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    Inc(Pass);
                end;
                Inc(HTask);
            end;
            Inc(MTask);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    WasErrors := AP_FP_Greater(Err,Threshold);
    if  not Silent then
    begin
        Write(Format('TESTING LDLT INVERSE'#13#10'',[]));
        Write(Format('ERROR:                                   %5.4e'#13#10'',[
            Err]));
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


procedure GenerateMatrix(var A : TReal2DArray;
     N : AlglibInteger;
     Task : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    if Task=0 then
    begin
        
        //
        // Zero matrix
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
    end;
    if Task=1 then
    begin
        
        //
        // Sparse matrix
        //
        I:=0;
        while I<=N-1 do
        begin
            J:=I+1;
            while J<=N-1 do
            begin
                if AP_FP_Greater(RandomReal,0.95) then
                begin
                    A[I,J] := 2*RandomReal-1;
                end
                else
                begin
                    A[I,J] := 0;
                end;
                A[J,I] := A[I,J];
                Inc(J);
            end;
            if AP_FP_Greater(RandomReal,0.95) then
            begin
                A[I,I] := (2*RandomInteger(2)-1)*(0.8+RandomReal);
            end
            else
            begin
                A[I,I] := 0;
            end;
            Inc(I);
        end;
    end;
    if Task=2 then
    begin
        
        //
        // Dense matrix
        //
        I:=0;
        while I<=N-1 do
        begin
            J:=I+1;
            while J<=N-1 do
            begin
                A[I,J] := 2*RandomReal-1;
                A[J,I] := A[I,J];
                Inc(J);
            end;
            A[I,I] := (2*RandomInteger(2)-1)*(0.7+RandomReal);
            Inc(I);
        end;
    end;
end;


procedure RestoreMatrix(var A : TReal2DArray;
     N : AlglibInteger;
     UpperIn : Boolean);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
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
                    A[I,J] := A[J,I];
                end;
            end
            else
            begin
                if I<J then
                begin
                    A[I,J] := A[J,I];
                end;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testinvldltunit_test_silent():Boolean;
begin
    Result := TestInvLDLT(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testinvldltunit_test():Boolean;
begin
    Result := TestInvLDLT(False);
end;


end.