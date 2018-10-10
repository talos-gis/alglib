unit testsblasunit;
interface
uses Math, Sysutils, Ap, sblas;

function TestSBLAS(Silent : Boolean):Boolean;
function testsblasunit_test_silent():Boolean;
function testsblasunit_test():Boolean;

implementation

function TestSBLAS(Silent : Boolean):Boolean;
var
    A : TReal2DArray;
    UA : TReal2DArray;
    LA : TReal2DArray;
    X : TReal1DArray;
    Y1 : TReal1DArray;
    Y2 : TReal1DArray;
    Y3 : TReal1DArray;
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
    Alpha : Double;
    V : Double;
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
            A[I,I] := 2*RandomReal-1;
            J:=I+1;
            while J<=N do
            begin
                A[I,J] := 2*RandomReal-1;
                A[J,I] := A[I,J];
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
                UA[I,J] := 0;
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
                LA[I,J] := 0;
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
                    X[I] := 2*RandomReal-1;
                    Inc(I);
                end;
                Alpha := 2*RandomReal-1;
                
                //
                // calculate A*x, UA*x, LA*x
                //
                I:=I1;
                while I<=I2 do
                begin
                    V := APVDotProduct(@A[I][0], I1, I2, @X[0], 1, I2-I1+1);
                    Y1[I-I1+1] := Alpha*V;
                    Inc(I);
                end;
                SymmetricMatrixVectorMultiply(UA, True, I1, I2, X, Alpha, Y2);
                SymmetricMatrixVectorMultiply(LA, False, I1, I2, X, Alpha, Y3);
                
                //
                // Calculate error
                //
                APVSub(@Y2[0], 1, I2-I1+1, @Y1[0], 1, I2-I1+1);
                V := APVDotProduct(@Y2[0], 1, I2-I1+1, @Y2[0], 1, I2-I1+1);
                MVErr := Max(MVErr, Sqrt(V));
                APVSub(@Y3[0], 1, I2-I1+1, @Y1[0], 1, I2-I1+1);
                V := APVDotProduct(@Y3[0], 1, I2-I1+1, @Y3[0], 1, I2-I1+1);
                MVErr := Max(MVErr, Sqrt(V));
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
        Write(Format('TESTING SYMMETRIC BLAS'#13#10'',[]));
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
function testsblasunit_test_silent():Boolean;
begin
    Result := TestSBLAS(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testsblasunit_test():Boolean;
begin
    Result := TestSBLAS(False);
end;


end.