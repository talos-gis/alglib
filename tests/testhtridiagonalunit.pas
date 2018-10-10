unit testhtridiagonalunit;
interface
uses Math, Sysutils, Ap, cblas, creflections, hblas, htridiagonal;

function TestHTridiagonal(Silent : Boolean):Boolean;
function testhtridiagonalunit_test_silent():Boolean;
function testhtridiagonalunit_test():Boolean;

implementation

procedure TestHTDProblem(const A : TComplex2DArray;
     N : AlglibInteger;
     var MatErr : Double;
     var OrtErr : Double);forward;


function TestHTridiagonal(Silent : Boolean):Boolean;
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
    A : TComplex2DArray;
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
                A[I,J] := C_Complex(0);
                Inc(J);
            end;
            Inc(I);
        end;
        TestHTDProblem(A, N, MatErr, OrtErr);
        
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
                A[I,I] := C_Complex(2*RandomReal-1);
                J:=I+1;
                while J<=N-1 do
                begin
                    A[I,J].X := 2*RandomReal-1;
                    A[I,J].Y := 2*RandomReal-1;
                    A[J,I] := Conj(A[I,J]);
                    Inc(J);
                end;
                Inc(I);
            end;
            TestHTDProblem(A, N, MatErr, OrtErr);
            
            //
            // Diagonal matrix
            //
            I:=0;
            while I<=N-1 do
            begin
                A[I,I] := C_Complex(2*RandomReal-1);
                J:=I+1;
                while J<=N-1 do
                begin
                    A[I,J] := C_Complex(0);
                    A[J,I] := C_Complex(0);
                    Inc(J);
                end;
                Inc(I);
            end;
            TestHTDProblem(A, N, MatErr, OrtErr);
            
            //
            // sparse matrix
            //
            I:=0;
            while I<=N-1 do
            begin
                A[I,I] := C_Complex(0);
                J:=I+1;
                while J<=N-1 do
                begin
                    A[I,J] := C_Complex(0);
                    A[J,I] := C_Complex(0);
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
                    A[I,J] := C_Complex(2*RandomReal-1);
                end
                else
                begin
                    A[I,J].X := 2*RandomReal-1;
                    A[I,J].Y := 2*RandomReal-1;
                    A[J,I] := Conj(A[I,J]);
                end;
                Inc(K);
            end;
            TestHTDProblem(A, N, MatErr, OrtErr);
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
        Write(Format('TESTING HERMITIAN TO TRIDIAGONAL'#13#10'',[]));
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


procedure TestHTDProblem(const A : TComplex2DArray;
     N : AlglibInteger;
     var MatErr : Double;
     var OrtErr : Double);
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
            MatErr := Max(MatErr, AbsComplex(C_Sub(T3[I,J],T[I,J])));
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
                OrtErr := Max(OrtErr, AbsComplex(C_SubR(V,1)));
            end
            else
            begin
                OrtErr := Max(OrtErr, AbsComplex(V));
            end;
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
            MatErr := Max(MatErr, AbsComplex(C_Sub(T3[I,J],T[I,J])));
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
                OrtErr := Max(OrtErr, AbsComplex(C_SubR(V,1)));
            end
            else
            begin
                OrtErr := Max(OrtErr, AbsComplex(V));
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testhtridiagonalunit_test_silent():Boolean;
begin
    Result := TestHTridiagonal(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testhtridiagonalunit_test():Boolean;
begin
    Result := TestHTridiagonal(False);
end;


end.