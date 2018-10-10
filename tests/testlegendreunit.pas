unit testlegendreunit;
interface
uses Math, Sysutils, Ap, legendre;

function TestLegendreCalculate(Silent : Boolean):Boolean;
function testlegendreunit_test_silent():Boolean;
function testlegendreunit_test():Boolean;

implementation

function TestLegendreCalculate(Silent : Boolean):Boolean;
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
    T : Double;
    WasErrors : Boolean;
begin
    Err := 0;
    SumErr := 0;
    CErr := 0;
    Threshold := 1.0E-9;
    WasErrors := False;
    
    //
    // Testing Legendre polynomials values
    //
    N:=0;
    while N<=10 do
    begin
        LegendreCoefficients(N, C);
        Pass:=1;
        while Pass<=10 do
        begin
            X := 2*RandomReal-1;
            V := LegendreCalculate(N, X);
            T := 1;
            I:=0;
            while I<=N do
            begin
                V := V-C[I]*T;
                T := T*X;
                Inc(I);
            end;
            Err := Max(Err, AbsReal(V));
            Inc(Pass);
        end;
        Inc(N);
    end;
    
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
            V := V+LegendreCalculate(N, X)*C[N];
            SumErr := Max(SumErr, AbsReal(V-LegendreSum(C, N, X)));
            Inc(N);
        end;
        Inc(Pass);
    end;
    
    //
    // Testing coefficients
    //
    LegendreCoefficients(0, C);
    CErr := Max(CErr, AbsReal(C[0]-1));
    LegendreCoefficients(1, C);
    CErr := Max(CErr, AbsReal(C[0]-0));
    CErr := Max(CErr, AbsReal(C[1]-1));
    LegendreCoefficients(2, C);
    CErr := Max(CErr, AbsReal(C[0]+AP_Double(1)/2));
    CErr := Max(CErr, AbsReal(C[1]-0));
    CErr := Max(CErr, AbsReal(C[2]-AP_Double(3)/2));
    LegendreCoefficients(3, C);
    CErr := Max(CErr, AbsReal(C[0]-0));
    CErr := Max(CErr, AbsReal(C[1]+AP_Double(3)/2));
    CErr := Max(CErr, AbsReal(C[2]-0));
    CErr := Max(CErr, AbsReal(C[3]-AP_Double(5)/2));
    LegendreCoefficients(4, C);
    CErr := Max(CErr, AbsReal(C[0]-AP_Double(3)/8));
    CErr := Max(CErr, AbsReal(C[1]-0));
    CErr := Max(CErr, AbsReal(C[2]+AP_Double(30)/8));
    CErr := Max(CErr, AbsReal(C[3]-0));
    CErr := Max(CErr, AbsReal(C[4]-AP_Double(35)/8));
    LegendreCoefficients(9, C);
    CErr := Max(CErr, AbsReal(C[0]-0));
    CErr := Max(CErr, AbsReal(C[1]-AP_Double(315)/128));
    CErr := Max(CErr, AbsReal(C[2]-0));
    CErr := Max(CErr, AbsReal(C[3]+AP_Double(4620)/128));
    CErr := Max(CErr, AbsReal(C[4]-0));
    CErr := Max(CErr, AbsReal(C[5]-AP_Double(18018)/128));
    CErr := Max(CErr, AbsReal(C[6]-0));
    CErr := Max(CErr, AbsReal(C[7]+AP_Double(25740)/128));
    CErr := Max(CErr, AbsReal(C[8]-0));
    CErr := Max(CErr, AbsReal(C[9]-AP_Double(12155)/128));
    
    //
    // Reporting
    //
    WasErrors := AP_FP_Greater(Err,Threshold) or AP_FP_Greater(SumErr,Threshold) or AP_FP_Greater(CErr,Threshold);
    if  not Silent then
    begin
        Write(Format('TESTING CALCULATION OF THE LEGENDRE POLYNOMIALS'#13#10'',[]));
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
function testlegendreunit_test_silent():Boolean;
begin
    Result := TestLegendreCalculate(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testlegendreunit_test():Boolean;
begin
    Result := TestLegendreCalculate(False);
end;


end.