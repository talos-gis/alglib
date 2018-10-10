unit testhblasunit;
interface
uses Math, Sysutils, Ap, hblas;

function TestHBLAS(Silent : Boolean):Boolean;
function testhblasunit_test_silent():Boolean;
function testhblasunit_test():Boolean;

implementation

function TestHBLAS(Silent : Boolean):Boolean;
var
    A : TComplex2DArray;
    UA : TComplex2DArray;
    LA : TComplex2DArray;
    X : TComplex1DArray;
    Y1 : TComplex1DArray;
    Y2 : TComplex1DArray;
    Y3 : TComplex1DArray;
    N : AlglibInteger;
    MaxN : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    GPass : AlglibInteger;
    WasErrors : Boolean;
    MVErr : Double;
    Threshold : Double;
    Alpha : Complex;
    V : Complex;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    MVErr := 0;
    WasErrors := False;
    MaxN := 10;
    Threshold := 1000*MachineEpsilon;
    
    //
    // Test MV
    //
    N:=2;
    while N<=MaxN do
    begin
        SetLength(A, N+1, N+1);
        SetLength(UA, N+1, N+1);
        SetLength(LA, N+1, N+1);
        SetLength(X, N+1);
        SetLength(Y1, N+1);
        SetLength(Y2, N+1);
        SetLength(Y3, N+1);
        
        //
        // fill A, UA, LA
        //
        I:=1;
        while I<=N do
        begin
            A[I,I].X := 2*RandomReal-1;
            A[I,I].Y := 0;
            J:=I+1;
            while J<=N do
            begin
                A[I,J].X := 2*RandomReal-1;
                A[I,J].Y := 2*RandomReal-1;
                A[J,I] := Conj(A[I,J]);
                Inc(J);
            end;
            Inc(I);
        end;
        I:=1;
        while I<=N do
        begin
            J:=1;
            while J<=N do
            begin
                UA[I,J] := C_Complex(0);
                Inc(J);
            end;
            Inc(I);
        end;
        I:=1;
        while I<=N do
        begin
            J:=I;
            while J<=N do
            begin
                UA[I,J] := A[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
        I:=1;
        while I<=N do
        begin
            J:=1;
            while J<=N do
            begin
                LA[I,J] := C_Complex(0);
                Inc(J);
            end;
            Inc(I);
        end;
        I:=1;
        while I<=N do
        begin
            J:=1;
            while J<=I do
            begin
                LA[I,J] := A[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // test on different I1, I2
        //
        I1:=1;
        while I1<=N do
        begin
            I2:=I1;
            while I2<=N do
            begin
                
                //
                // Fill X, choose Alpha
                //
                I:=1;
                while I<=I2-I1+1 do
                begin
                    X[I].X := 2*RandomReal-1;
                    X[I].Y := 2*RandomReal-1;
                    Inc(I);
                end;
                Alpha.X := 2*RandomReal-1;
                Alpha.Y := 2*RandomReal-1;
                
                //
                // calculate A*x, UA*x, LA*x
                //
                I:=I1;
                while I<=I2 do
                begin
                    i1_ := (1)-(I1);
                    V := C_Complex(0.0);
                    for i_ := I1 to I2 do
                    begin
                        V := C_Add(V,C_Mul(A[I,i_],X[i_+i1_]));
                    end;
                    Y1[I-I1+1] := C_Mul(Alpha,V);
                    Inc(I);
                end;
                HermitianMatrixVectorMultiply(UA, True, I1, I2, X, Alpha, Y2);
                HermitianMatrixVectorMultiply(LA, False, I1, I2, X, Alpha, Y3);
                
                //
                // Calculate error
                //
                for i_ := 1 to I2-I1+1 do
                begin
                    Y2[i_] := C_Sub(Y2[i_], Y1[i_]);
                end;
                V := C_Complex(0.0);
                for i_ := 1 to I2-I1+1 do
                begin
                    V := C_Add(V,C_Mul(Y2[i_],Conj(Y2[i_])));
                end;
                MVErr := Max(MVErr, Sqrt(AbsComplex(V)));
                for i_ := 1 to I2-I1+1 do
                begin
                    Y3[i_] := C_Sub(Y3[i_], Y1[i_]);
                end;
                V := C_Complex(0.0);
                for i_ := 1 to I2-I1+1 do
                begin
                    V := C_Add(V,C_Mul(Y3[i_],Conj(Y3[i_])));
                end;
                MVErr := Max(MVErr, Sqrt(AbsComplex(V)));
                Inc(I2);
            end;
            Inc(I1);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    WasErrors := AP_FP_Greater(MVErr,Threshold);
    if  not Silent then
    begin
        Write(Format('TESTING HERMITIAN BLAS'#13#10'',[]));
        Write(Format('MV error:                                %5.4e'#13#10'',[
            MVErr]));
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
function testhblasunit_test_silent():Boolean;
begin
    Result := TestHBLAS(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testhblasunit_test():Boolean;
begin
    Result := TestHBLAS(False);
end;


end.