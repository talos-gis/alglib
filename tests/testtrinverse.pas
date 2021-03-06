unit testtrinverse;
interface
uses Math, Sysutils, Ap, trinverse;

function TestTRInv(Silent : Boolean):Boolean;
function testtrinverse_test_silent():Boolean;
function testtrinverse_test():Boolean;

implementation

(*************************************************************************
Main unittest subroutine
*************************************************************************)
function TestTRInv(Silent : Boolean):Boolean;
var
    ShortMN : AlglibInteger;
    MaxN : AlglibInteger;
    PassCount : AlglibInteger;
    Threshold : Double;
    A : TReal2DArray;
    B : TReal2DArray;
    N : AlglibInteger;
    Pass : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    Task : AlglibInteger;
    IsUpper : Boolean;
    Isunit : Boolean;
    V : Double;
    InvFailed : Boolean;
    InvErrors : Boolean;
    StructErrors : Boolean;
    WasErrors : Boolean;
    i_ : AlglibInteger;
begin
    InvFailed := False;
    InvErrors := False;
    StructErrors := False;
    WasErrors := False;
    MaxN := 15;
    PassCount := 5;
    Threshold := 5*100*MachineEpsilon;
    
    //
    // Test
    //
    N:=1;
    while N<=MaxN do
    begin
        SetLength(A, N-1+1, N-1+1);
        SetLength(B, N-1+1, N-1+1);
        Task:=0;
        while Task<=3 do
        begin
            Pass:=1;
            while Pass<=PassCount do
            begin
                
                //
                // Determine task
                //
                IsUpper := Task mod 2=0;
                Isunit := Task div 2 mod 2=0;
                
                //
                // Generate matrix
                //
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        if I=J then
                        begin
                            A[I,J] := 1.5+RandomReal;
                        end
                        else
                        begin
                            A[I,J] := 2*RandomReal-1;
                        end;
                        B[I,J] := A[I,J];
                        Inc(J);
                    end;
                    Inc(I);
                end;
                
                //
                // Inverse
                //
                if  not RMatrixTRInverse(B, N, IsUpper, Isunit) then
                begin
                    InvFailed := True;
                    Inc(Pass);
                    Continue;
                end;
                
                //
                // Structural test
                //
                if Isunit then
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        StructErrors := StructErrors or AP_FP_Neq(A[I,I],B[I,I]);
                        Inc(I);
                    end;
                end;
                if IsUpper then
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=I-1 do
                        begin
                            StructErrors := StructErrors or AP_FP_Neq(A[I,J],B[I,J]);
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                end
                else
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=I+1;
                        while J<=N-1 do
                        begin
                            StructErrors := StructErrors or AP_FP_Neq(A[I,J],B[I,J]);
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                end;
                
                //
                // Inverse test
                //
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        if (J<I) and IsUpper or (J>I) and  not IsUpper then
                        begin
                            A[I,J] := 0;
                            B[I,J] := 0;
                        end;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                if Isunit then
                begin
                    I:=0;
                    while I<=N-1 do
                    begin
                        A[I,I] := 1;
                        B[I,I] := 1;
                        Inc(I);
                    end;
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
                            V := V + A[I,i_]*B[i_,J];
                        end;
                        if J<>I then
                        begin
                            InvErrors := InvErrors or AP_FP_Greater(AbsReal(V),Threshold);
                        end
                        else
                        begin
                            InvErrors := InvErrors or AP_FP_Greater(AbsReal(V-1),Threshold);
                        end;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                Inc(Pass);
            end;
            Inc(Task);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    WasErrors := InvErrors or StructErrors or InvFailed;
    if  not Silent then
    begin
        Write(Format('TESTING TRIANGULAR INVERSE'#13#10'',[]));
        if InvFailed then
        begin
            Write(Format('SOME INVERSIONS FAILED'#13#10'',[]));
        end;
        Write(Format('INVERSION TEST:                          ',[]));
        if InvErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('STRUCTURE TEST:                          ',[]));
        if StructErrors then
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
function testtrinverse_test_silent():Boolean;
begin
    Result := TestTRInv(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testtrinverse_test():Boolean;
begin
    Result := TestTRInv(False);
end;


end.