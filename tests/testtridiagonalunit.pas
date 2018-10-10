unit testtridiagonalunit;
interface
uses Math, Sysutils, Ap, sblas, reflections, tridiagonal;

function TestTridiagonal(Silent : Boolean):Boolean;
function testtridiagonalunit_test_silent():Boolean;
function testtridiagonalunit_test():Boolean;

implementation

procedure TestSTDProblem(const A : TReal2DArray;
     N : AlglibInteger;
     var MatErr : Double;
     var OrtErr : Double);forward;


function TestTridiagonal(Silent : Boolean):Boolean;
var
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    MaxN : AlglibInteger;
    MatErr : Double;
    OrtErr : Double;
    N : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    A : TReal2DArray;
    Threshold : Double;
    WasErrors : Boolean;
begin
    MatErr := 0;
    OrtErr := 0;
    WasErrors := False;
    Threshold := 1000*MachineEpsilon;
    MaxN := 15;
    PassCount := 20;
    N:=1;
    while N<=MaxN do
    begin
        SetLength(A, N-1+1, N-1+1);
        
        //
        // Test zero matrix
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
        TestSTDProblem(A, N, MatErr, OrtErr);
        
        //
        // Test other matrix types
        //
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // Test dense matrix
            //
            I:=0;
            while I<=N-1 do
            begin
                A[I,I] := 2*RandomReal-1;
                J:=I+1;
                while J<=N-1 do
                begin
                    A[I,J] := 2*RandomReal-1;
                    A[J,I] := A[I,J];
                    Inc(J);
                end;
                Inc(I);
            end;
            TestSTDProblem(A, N, MatErr, OrtErr);
            
            //
            // Diagonal matrix
            //
            I:=0;
            while I<=N-1 do
            begin
                A[I,I] := 2*RandomReal-1;
                J:=I+1;
                while J<=N-1 do
                begin
                    A[I,J] := 0;
                    A[J,I] := 0;
                    Inc(J);
                end;
                Inc(I);
            end;
            TestSTDProblem(A, N, MatErr, OrtErr);
            
            //
            // sparse matrix
            //
            I:=0;
            while I<=N-1 do
            begin
                A[I,I] := 0;
                J:=I+1;
                while J<=N-1 do
                begin
                    A[I,J] := 0;
                    A[J,I] := 0;
                    Inc(J);
                end;
                Inc(I);
            end;
            K:=1;
            while K<=2 do
            begin
                I := RandomInteger(N);
                J := RandomInteger(N);
                if I=J then
                begin
                    A[I,J] := 2*RandomReal-1;
                end
                else
                begin
                    A[I,J] := 2*RandomReal-1;
                    A[J,I] := A[I,J];
                end;
                Inc(K);
            end;
            TestSTDProblem(A, N, MatErr, OrtErr);
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    WasErrors := AP_FP_Greater(MatErr,Threshold) or AP_FP_Greater(OrtErr,Threshold);
    if  not Silent then
    begin
        Write(Format('TESTING SYMMETRIC TO TRIDIAGONAL'#13#10'',[]));
        Write(Format('Matrix error:                            %5.4e'#13#10'',[
            MatErr]));
        Write(Format('Q orthogonality error:                   %5.4e'#13#10'',[
            OrtErr]));
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


procedure TestSTDProblem(const A : TReal2DArray;
     N : AlglibInteger;
     var MatErr : Double;
     var OrtErr : Double);
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
            MatErr := Max(MatErr, AbsComplex(C_Complex(T3[I,J]-T[I,J])));
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
                OrtErr := Max(OrtErr, AbsComplex(C_Complex(V-1)));
            end
            else
            begin
                OrtErr := Max(OrtErr, AbsComplex(C_Complex(V)));
            end;
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
            MatErr := Max(MatErr, AbsComplex(C_Complex(T3[I,J]-T[I,J])));
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
                OrtErr := Max(OrtErr, AbsComplex(C_Complex(V-1)));
            end
            else
            begin
                OrtErr := Max(OrtErr, AbsComplex(C_Complex(V)));
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testtridiagonalunit_test_silent():Boolean;
begin
    Result := TestTridiagonal(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testtridiagonalunit_test():Boolean;
begin
    Result := TestTridiagonal(False);
end;


end.