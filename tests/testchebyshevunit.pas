unit testchebyshevunit;
interface
uses Math, Sysutils, Ap, chebyshev;

function TestChebyshev(Silent : Boolean):Boolean;
function testchebyshevunit_test_silent():Boolean;
function testchebyshevunit_test():Boolean;

implementation

function TestChebyshev(Silent : Boolean):Boolean;
var
    Err : Double;
    SumErr : Double;
    CErr : Double;
    FErr : Double;
    Threshold : Double;
    X : Double;
    V : Double;
    T : Double;
    Pass : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    N : AlglibInteger;
    MaxN : AlglibInteger;
    C : TReal1DArray;
    P1 : TReal1DArray;
    P2 : TReal1DArray;
    A : TReal2DArray;
    WasErrors : Boolean;
begin
    Err := 0;
    SumErr := 0;
    CErr := 0;
    FErr := 0;
    Threshold := 1.0E-9;
    WasErrors := False;
    
    //
    // Testing Chebyshev polynomials of the first kind
    //
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, 0, 0.00)-1));
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, 0, 0.33)-1));
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, 0, -0.42)-1));
    X := 0.2;
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, 1, X)-0.2));
    X := 0.4;
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, 1, X)-0.4));
    X := 0.6;
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, 1, X)-0.6));
    X := 0.8;
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, 1, X)-0.8));
    X := 1.0;
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, 1, X)-1.0));
    X := 0.2;
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, 2, X)+0.92));
    X := 0.4;
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, 2, X)+0.68));
    X := 0.6;
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, 2, X)+0.28));
    X := 0.8;
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, 2, X)-0.28));
    X := 1.0;
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, 2, X)-1.00));
    N := 10;
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, N, 0.2)-0.4284556288));
    N := 11;
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, N, 0.2)+0.7996160205));
    N := 12;
    Err := Max(Err, AbsReal(ChebyshevCalculate(1, N, 0.2)+0.7483020370));
    
    //
    // Testing Chebyshev polynomials of the second kind
    //
    N := 0;
    Err := Max(Err, AbsReal(ChebyshevCalculate(2, N, 0.2)-1.0000000000));
    N := 1;
    Err := Max(Err, AbsReal(ChebyshevCalculate(2, N, 0.2)-0.4000000000));
    N := 2;
    Err := Max(Err, AbsReal(ChebyshevCalculate(2, N, 0.2)+0.8400000000));
    N := 3;
    Err := Max(Err, AbsReal(ChebyshevCalculate(2, N, 0.2)+0.7360000000));
    N := 4;
    Err := Max(Err, AbsReal(ChebyshevCalculate(2, N, 0.2)-0.5456000000));
    N := 10;
    Err := Max(Err, AbsReal(ChebyshevCalculate(2, N, 0.2)-0.6128946176));
    N := 11;
    Err := Max(Err, AbsReal(ChebyshevCalculate(2, N, 0.2)+0.6770370970));
    N := 12;
    Err := Max(Err, AbsReal(ChebyshevCalculate(2, N, 0.2)+0.8837094564));
    
    //
    // Testing Clenshaw summation
    //
    MaxN := 20;
    SetLength(C, MaxN+1);
    K:=1;
    while K<=2 do
    begin
        Pass:=1;
        while Pass<=10 do
        begin
            X := 2*RandomReal-1;
            V := 0;
            N:=0;
            while N<=MaxN do
            begin
                C[N] := 2*RandomReal-1;
                V := V+ChebyshevCalculate(K, N, X)*C[N];
                SumErr := Max(SumErr, AbsReal(V-ChebyshevSum(C, K, N, X)));
                Inc(N);
            end;
            Inc(Pass);
        end;
        Inc(K);
    end;
    
    //
    // Testing coefficients
    //
    ChebyshevCoefficients(0, C);
    CErr := Max(CErr, AbsReal(C[0]-1));
    ChebyshevCoefficients(1, C);
    CErr := Max(CErr, AbsReal(C[0]-0));
    CErr := Max(CErr, AbsReal(C[1]-1));
    ChebyshevCoefficients(2, C);
    CErr := Max(CErr, AbsReal(C[0]+1));
    CErr := Max(CErr, AbsReal(C[1]-0));
    CErr := Max(CErr, AbsReal(C[2]-2));
    ChebyshevCoefficients(3, C);
    CErr := Max(CErr, AbsReal(C[0]-0));
    CErr := Max(CErr, AbsReal(C[1]+3));
    CErr := Max(CErr, AbsReal(C[2]-0));
    CErr := Max(CErr, AbsReal(C[3]-4));
    ChebyshevCoefficients(4, C);
    CErr := Max(CErr, AbsReal(C[0]-1));
    CErr := Max(CErr, AbsReal(C[1]-0));
    CErr := Max(CErr, AbsReal(C[2]+8));
    CErr := Max(CErr, AbsReal(C[3]-0));
    CErr := Max(CErr, AbsReal(C[4]-8));
    ChebyshevCoefficients(9, C);
    CErr := Max(CErr, AbsReal(C[0]-0));
    CErr := Max(CErr, AbsReal(C[1]-9));
    CErr := Max(CErr, AbsReal(C[2]-0));
    CErr := Max(CErr, AbsReal(C[3]+120));
    CErr := Max(CErr, AbsReal(C[4]-0));
    CErr := Max(CErr, AbsReal(C[5]-432));
    CErr := Max(CErr, AbsReal(C[6]-0));
    CErr := Max(CErr, AbsReal(C[7]+576));
    CErr := Max(CErr, AbsReal(C[8]-0));
    CErr := Max(CErr, AbsReal(C[9]-256));
    
    //
    // Testing FromChebyshev
    //
    MaxN := 10;
    SetLength(A, MaxN+1, MaxN+1);
    I:=0;
    while I<=MaxN do
    begin
        J:=0;
        while J<=MaxN do
        begin
            A[I,J] := 0;
            Inc(J);
        end;
        ChebyshevCoefficients(I, C);
        APVMove(@A[I][0], 0, I, @C[0], 0, I);
        Inc(I);
    end;
    SetLength(C, MaxN+1);
    SetLength(P1, MaxN+1);
    N:=0;
    while N<=MaxN do
    begin
        Pass:=1;
        while Pass<=10 do
        begin
            I:=0;
            while I<=N do
            begin
                P1[I] := 0;
                Inc(I);
            end;
            I:=0;
            while I<=N do
            begin
                C[I] := 2*RandomReal-1;
                V := C[I];
                APVAdd(@P1[0], 0, I, @A[I][0], 0, I, V);
                Inc(I);
            end;
            FromChebyshev(C, N, P2);
            I:=0;
            while I<=N do
            begin
                FErr := Max(FErr, AbsReal(P1[I]-P2[I]));
                Inc(I);
            end;
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // Reporting
    //
    WasErrors := AP_FP_Greater(Err,Threshold) or AP_FP_Greater(SumErr,Threshold) or AP_FP_Greater(CErr,Threshold) or AP_FP_Greater(FErr,Threshold);
    if  not Silent then
    begin
        Write(Format('TESTING CALCULATION OF THE CHEBYSHEV POLYNOMIALS'#13#10'',[]));
        Write(Format('Max error against table                   %5.3e'#13#10'',[
            Err]));
        Write(Format('Summation error                           %5.3e'#13#10'',[
            SumErr]));
        Write(Format('Coefficients error                        %5.3e'#13#10'',[
            CErr]));
        Write(Format('FrobChebyshev error                       %5.3e'#13#10'',[
            FErr]));
        Write(Format('Threshold                                 %5.3e'#13#10'',[
            Threshold]));
        if  not WasErrors then
        begin
            Write(Format('TEST PASSED'#13#10'',[]));
        end
        else
        begin
            Write(Format('TEST FAILED'#13#10'',[]));
        end;
    end;
    Result :=  not WasErrors;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testchebyshevunit_test_silent():Boolean;
begin
    Result := TestChebyshev(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testchebyshevunit_test():Boolean;
begin
    Result := TestChebyshev(False);
end;


end.