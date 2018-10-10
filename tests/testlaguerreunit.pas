unit testlaguerreunit;
interface
uses Math, Sysutils, Ap, laguerre;

function TestLaguerreCalculate(Silent : Boolean):Boolean;
function testlaguerreunit_test_silent():Boolean;
function testlaguerreunit_test():Boolean;

implementation

function TestLaguerreCalculate(Silent : Boolean):Boolean;
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
    // Testing Laguerre polynomials
    //
    N := 0;
    Err := Max(Err, AbsReal(LaguerreCalculate(N, 0.5)-1.0000000000));
    N := 1;
    Err := Max(Err, AbsReal(LaguerreCalculate(N, 0.5)-0.5000000000));
    N := 2;
    Err := Max(Err, AbsReal(LaguerreCalculate(N, 0.5)-0.1250000000));
    N := 3;
    Err := Max(Err, AbsReal(LaguerreCalculate(N, 0.5)+0.1458333333));
    N := 4;
    Err := Max(Err, AbsReal(LaguerreCalculate(N, 0.5)+0.3307291667));
    N := 5;
    Err := Max(Err, AbsReal(LaguerreCalculate(N, 0.5)+0.4455729167));
    N := 6;
    Err := Max(Err, AbsReal(LaguerreCalculate(N, 0.5)+0.5041449653));
    N := 7;
    Err := Max(Err, AbsReal(LaguerreCalculate(N, 0.5)+0.5183392237));
    N := 8;
    Err := Max(Err, AbsReal(LaguerreCalculate(N, 0.5)+0.4983629984));
    N := 9;
    Err := Max(Err, AbsReal(LaguerreCalculate(N, 0.5)+0.4529195204));
    N := 10;
    Err := Max(Err, AbsReal(LaguerreCalculate(N, 0.5)+0.3893744141));
    N := 11;
    Err := Max(Err, AbsReal(LaguerreCalculate(N, 0.5)+0.3139072988));
    N := 12;
    Err := Max(Err, AbsReal(LaguerreCalculate(N, 0.5)+0.2316496389));
    
    //
    // Testing Clenshaw summation
    //
    MaxN := 20;
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
            V := V+LaguerreCalculate(N, X)*C[N];
            SumErr := Max(SumErr, AbsReal(V-LaguerreSum(C, N, X)));
            Inc(N);
        end;
        Inc(Pass);
    end;
    
    //
    // Testing coefficients
    //
    LaguerreCoefficients(0, C);
    CErr := Max(CErr, AbsReal(C[0]-1));
    LaguerreCoefficients(1, C);
    CErr := Max(CErr, AbsReal(C[0]-1));
    CErr := Max(CErr, AbsReal(C[1]+1));
    LaguerreCoefficients(2, C);
    CErr := Max(CErr, AbsReal(C[0]-AP_Double(2)/2));
    CErr := Max(CErr, AbsReal(C[1]+AP_Double(4)/2));
    CErr := Max(CErr, AbsReal(C[2]-AP_Double(1)/2));
    LaguerreCoefficients(3, C);
    CErr := Max(CErr, AbsReal(C[0]-AP_Double(6)/6));
    CErr := Max(CErr, AbsReal(C[1]+AP_Double(18)/6));
    CErr := Max(CErr, AbsReal(C[2]-AP_Double(9)/6));
    CErr := Max(CErr, AbsReal(C[3]+AP_Double(1)/6));
    LaguerreCoefficients(4, C);
    CErr := Max(CErr, AbsReal(C[0]-AP_Double(24)/24));
    CErr := Max(CErr, AbsReal(C[1]+AP_Double(96)/24));
    CErr := Max(CErr, AbsReal(C[2]-AP_Double(72)/24));
    CErr := Max(CErr, AbsReal(C[3]+AP_Double(16)/24));
    CErr := Max(CErr, AbsReal(C[4]-AP_Double(1)/24));
    LaguerreCoefficients(5, C);
    CErr := Max(CErr, AbsReal(C[0]-AP_Double(120)/120));
    CErr := Max(CErr, AbsReal(C[1]+AP_Double(600)/120));
    CErr := Max(CErr, AbsReal(C[2]-AP_Double(600)/120));
    CErr := Max(CErr, AbsReal(C[3]+AP_Double(200)/120));
    CErr := Max(CErr, AbsReal(C[4]-AP_Double(25)/120));
    CErr := Max(CErr, AbsReal(C[5]+AP_Double(1)/120));
    LaguerreCoefficients(6, C);
    CErr := Max(CErr, AbsReal(C[0]-AP_Double(720)/720));
    CErr := Max(CErr, AbsReal(C[1]+AP_Double(4320)/720));
    CErr := Max(CErr, AbsReal(C[2]-AP_Double(5400)/720));
    CErr := Max(CErr, AbsReal(C[3]+AP_Double(2400)/720));
    CErr := Max(CErr, AbsReal(C[4]-AP_Double(450)/720));
    CErr := Max(CErr, AbsReal(C[5]+AP_Double(36)/720));
    CErr := Max(CErr, AbsReal(C[6]-AP_Double(1)/720));
    
    //
    // Reporting
    //
    WasErrors := AP_FP_Greater(Err,Threshold) or AP_FP_Greater(SumErr,Threshold) or AP_FP_Greater(CErr,Threshold);
    if  not Silent then
    begin
        Write(Format('TESTING CALCULATION OF THE LAGUERRE POLYNOMIALS'#13#10'',[]));
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
function testlaguerreunit_test_silent():Boolean;
begin
    Result := TestLaguerreCalculate(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testlaguerreunit_test():Boolean;
begin
    Result := TestLaguerreCalculate(False);
end;


end.