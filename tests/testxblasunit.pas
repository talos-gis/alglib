unit testxblasunit;
interface
uses Math, Sysutils, Ap, tsort, xblas;

function TestXBlas(Silent : Boolean):Boolean;
function testxblasunit_test_silent():Boolean;
function testxblasunit_test():Boolean;

implementation

const
    XChunk = 1048576;
    XChunkCount = 4;

function TestXBlas(Silent : Boolean):Boolean;
var
    ApproxErrors : Boolean;
    ExactnessErrors : Boolean;
    WasErrors : Boolean;
    ApproxThreshold : Double;
    MaxN : AlglibInteger;
    PassCount : AlglibInteger;
    N : AlglibInteger;
    I : AlglibInteger;
    Pass : AlglibInteger;
    V1 : Double;
    V2 : Double;
    V2Err : Double;
    X : TReal1DArray;
    Y : TReal1DArray;
    Temp : TReal1DArray;
    B : Double;
    S : Double;
begin
    ApproxErrors := False;
    ExactnessErrors := False;
    WasErrors := False;
    ApproxThreshold := 1000*MachineEpsilon;
    MaxN := 1000;
    PassCount := 10;
    
    //
    // tests:
    // 1. ability to calculate dot product
    // 2. higher precision
    //
    N:=1;
    while N<=MaxN do
    begin
        Pass:=1;
        while Pass<=PassCount do
        begin
            SetLength(X, N);
            SetLength(Y, N);
            SetLength(Temp, N);
            
            //
            //  ability to approximately calculate dot product
            //
            I:=0;
            while I<=N-1 do
            begin
                if AP_FP_Greater(RandomReal,0.2) then
                begin
                    X[I] := 2*RandomReal-1;
                end
                else
                begin
                    X[I] := 0;
                end;
                if AP_FP_Greater(RandomReal,0.2) then
                begin
                    Y[I] := 2*RandomReal-1;
                end
                else
                begin
                    Y[I] := 0;
                end;
                Inc(I);
            end;
            V1 := APVDotProduct(@X[0], 0, N-1, @Y[0], 0, N-1);
            XDot(X, Y, N, Temp, V2, V2Err);
            ApproxErrors := ApproxErrors or AP_FP_Greater(AbsReal(V1-V2),ApproxThreshold);
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // test of precision
    //
    N := 10000;
    SetLength(X, N);
    SetLength(Y, N);
    SetLength(Temp, N);
    Pass:=0;
    while Pass<=PassCount-1 do
    begin
        
        //
        // First test: X + X + ... + X - X - X - ... - X = 1*X
        //
        S := Exp(Max(Pass, 50));
        if (Pass=PassCount-1) and (Pass>1) then
        begin
            S := MaxRealNumber;
        end;
        Y[0] := (2*RandomReal-1)*S*Sqrt(2*RandomReal);
        I:=1;
        while I<=N-1 do
        begin
            Y[I] := Y[0];
            Inc(I);
        end;
        I:=0;
        while I<=N div 2-1 do
        begin
            X[I] := 1;
            Inc(I);
        end;
        I:=N div 2;
        while I<=N-2 do
        begin
            X[I] := -1;
            Inc(I);
        end;
        X[N-1] := 0;
        XDot(X, Y, N, Temp, V2, V2Err);
        ExactnessErrors := ExactnessErrors or AP_FP_Less(V2Err,0);
        ExactnessErrors := ExactnessErrors or AP_FP_Greater(V2Err,4*MachineEpsilon*AbsReal(Y[0]));
        ExactnessErrors := ExactnessErrors or AP_FP_Greater(AbsReal(V2-Y[0]),V2Err);
        
        //
        // First test: X + X + ... + X = N*X
        //
        S := Exp(Max(Pass, 50));
        if (Pass=PassCount-1) and (Pass>1) then
        begin
            S := MaxRealNumber;
        end;
        Y[0] := (2*RandomReal-1)*S*Sqrt(2*RandomReal);
        I:=1;
        while I<=N-1 do
        begin
            Y[I] := Y[0];
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            X[I] := 1;
            Inc(I);
        end;
        XDot(X, Y, N, Temp, V2, V2Err);
        ExactnessErrors := ExactnessErrors or AP_FP_Less(V2Err,0);
        ExactnessErrors := ExactnessErrors or AP_FP_Greater(V2Err,4*MachineEpsilon*AbsReal(Y[0])*N);
        ExactnessErrors := ExactnessErrors or AP_FP_Greater(AbsReal(V2-N*Y[0]),V2Err);
        Inc(Pass);
    end;
    
    //
    // report
    //
    WasErrors := ApproxErrors or ExactnessErrors;
    if  not Silent then
    begin
        Write(Format('TESTING XBLAS'#13#10'',[]));
        Write(Format('APPROX.TESTS:                            ',[]));
        if ApproxErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('EXACT TESTS:                             ',[]));
        if ExactnessErrors then
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
    
    //
    // end
    //
    Result :=  not WasErrors;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testxblasunit_test_silent():Boolean;
begin
    Result := TestXBlas(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testxblasunit_test():Boolean;
begin
    Result := TestXBlas(False);
end;


end.