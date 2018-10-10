unit testhermiteunit;
interface
uses Math, Sysutils, Ap, hermite;

function TestHermiteCalculate(Silent : Boolean):Boolean;
function testhermiteunit_test_silent():Boolean;
function testhermiteunit_test():Boolean;

implementation

function TestHermiteCalculate(Silent : Boolean):Boolean;
var
    Err : Double;
    SumErr : Double;
    CErr : Double;
    Threshold : Double;
    N : AlglibInteger;
    MaxN : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    Pass : AlglibInteger;
    C : TReal1DArray;
    X : Double;
    V : Double;
    WasErrors : Boolean;
begin
    Err := 0;
    SumErr := 0;
    CErr := 0;
    Threshold := 1.0E-9;
    WasErrors := False;
    
    //
    // Testing Hermite polynomials
    //
    N := 0;
    Err := Max(Err, AbsReal(HermiteCalculate(N, 1)-1));
    N := 1;
    Err := Max(Err, AbsReal(HermiteCalculate(N, 1)-2));
    N := 2;
    Err := Max(Err, AbsReal(HermiteCalculate(N, 1)-2));
    N := 3;
    Err := Max(Err, AbsReal(HermiteCalculate(N, 1)+4));
    N := 4;
    Err := Max(Err, AbsReal(HermiteCalculate(N, 1)+20));
    N := 5;
    Err := Max(Err, AbsReal(HermiteCalculate(N, 1)+8));
    N := 6;
    Err := Max(Err, AbsReal(HermiteCalculate(N, 1)-184));
    N := 7;
    Err := Max(Err, AbsReal(HermiteCalculate(N, 1)-464));
    N := 11;
    Err := Max(Err, AbsReal(HermiteCalculate(N, 1)-230848));
    N := 12;
    Err := Max(Err, AbsReal(HermiteCalculate(N, 1)-280768));
    
    //
    // Testing Clenshaw summation
    //
    MaxN := 10;
    SetLength(C, MaxN+1);
    Pass:=1;
    while Pass<=10 do
    begin
        X := 2*RandomReal-1;
        V := 0;
        N:=0;
        while N<=MaxN do
        begin
            C[N] := 2*RandomReal-1;
            V := V+HermiteCalculate(N, X)*C[N];
            SumErr := Max(SumErr, AbsReal(V-HermiteSum(C, N, X)));
            Inc(N);
        end;
        Inc(Pass);
    end;
    
    //
    // Testing coefficients
    //
    HermiteCoefficients(0, C);
    CErr := Max(CErr, AbsReal(C[0]-1));
    HermiteCoefficients(1, C);
    CErr := Max(CErr, AbsReal(C[0]-0));
    CErr := Max(CErr, AbsReal(C[1]-2));
    HermiteCoefficients(2, C);
    CErr := Max(CErr, AbsReal(C[0]+2));
    CErr := Max(CErr, AbsReal(C[1]-0));
    CErr := Max(CErr, AbsReal(C[2]-4));
    HermiteCoefficients(3, C);
    CErr := Max(CErr, AbsReal(C[0]-0));
    CErr := Max(CErr, AbsReal(C[1]+12));
    CErr := Max(CErr, AbsReal(C[2]-0));
    CErr := Max(CErr, AbsReal(C[3]-8));
    HermiteCoefficients(4, C);
    CErr := Max(CErr, AbsReal(C[0]-12));
    CErr := Max(CErr, AbsReal(C[1]-0));
    CErr := Max(CErr, AbsReal(C[2]+48));
    CErr := Max(CErr, AbsReal(C[3]-0));
    CErr := Max(CErr, AbsReal(C[4]-16));
    HermiteCoefficients(5, C);
    CErr := Max(CErr, AbsReal(C[0]-0));
    CErr := Max(CErr, AbsReal(C[1]-120));
    CErr := Max(CErr, AbsReal(C[2]-0));
    CErr := Max(CErr, AbsReal(C[3]+160));
    CErr := Max(CErr, AbsReal(C[4]-0));
    CErr := Max(CErr, AbsReal(C[5]-32));
    HermiteCoefficients(6, C);
    CErr := Max(CErr, AbsReal(C[0]+120));
    CErr := Max(CErr, AbsReal(C[1]-0));
    CErr := Max(CErr, AbsReal(C[2]-720));
    CErr := Max(CErr, AbsReal(C[3]-0));
    CErr := Max(CErr, AbsReal(C[4]+480));
    CErr := Max(CErr, AbsReal(C[5]-0));
    CErr := Max(CErr, AbsReal(C[6]-64));
    
    //
    // Reporting
    //
    WasErrors := AP_FP_Greater(Err,Threshold) or AP_FP_Greater(SumErr,Threshold) or AP_FP_Greater(CErr,Threshold);
    if  not Silent then
    begin
        Write(Format('TESTING CALCULATION OF THE HERMITE POLYNOMIALS'#13#10'',[]));
        Write(Format('Max error                                 %5.3e'#13#10'',[
            Err]));
        Write(Format('Summation error                           %5.3e'#13#10'',[
            SumErr]));
        Write(Format('Coefficients error                        %5.3e'#13#10'',[
            CErr]));
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
function testhermiteunit_test_silent():Boolean;
begin
    Result := TestHermiteCalculate(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testhermiteunit_test():Boolean;
begin
    Result := TestHermiteCalculate(False);
end;


end.