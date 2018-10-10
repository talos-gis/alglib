unit testxblasunit;
interface
uses Math, Sysutils, Ap, xblas;

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
    RV1 : Double;
    RV2 : Double;
    RV2Err : Double;
    CV1 : Complex;
    CV2 : Complex;
    CV2Err : Double;
    CV2ErrX : Double;
    CV2ErrY : Double;
    RX : TReal1DArray;
    RY : TReal1DArray;
    CX : TComplex1DArray;
    CY : TComplex1DArray;
    Temp : TReal1DArray;
    B : Double;
    S : Double;
    i_ : AlglibInteger;
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
            
            //
            //  ability to approximately calculate real dot product
            //
            SetLength(RX, N);
            SetLength(RY, N);
            SetLength(Temp, N);
            I:=0;
            while I<=N-1 do
            begin
                if AP_FP_Greater(RandomReal,0.2) then
                begin
                    RX[I] := 2*RandomReal-1;
                end
                else
                begin
                    RX[I] := 0;
                end;
                if AP_FP_Greater(RandomReal,0.2) then
                begin
                    RY[I] := 2*RandomReal-1;
                end
                else
                begin
                    RY[I] := 0;
                end;
                Inc(I);
            end;
            RV1 := APVDotProduct(@RX[0], 0, N-1, @RY[0], 0, N-1);
            XDot(RX, RY, N, Temp, RV2, RV2Err);
            ApproxErrors := ApproxErrors or AP_FP_Greater(AbsReal(RV1-RV2),ApproxThreshold);
            
            //
            //  ability to approximately calculate complex dot product
            //
            SetLength(CX, N);
            SetLength(CY, N);
            SetLength(Temp, 2*N);
            I:=0;
            while I<=N-1 do
            begin
                if AP_FP_Greater(RandomReal,0.2) then
                begin
                    CX[I].X := 2*RandomReal-1;
                    CX[I].Y := 2*RandomReal-1;
                end
                else
                begin
                    CX[I] := C_Complex(0);
                end;
                if AP_FP_Greater(RandomReal,0.2) then
                begin
                    CY[I].X := 2*RandomReal-1;
                    CY[I].Y := 2*RandomReal-1;
                end
                else
                begin
                    CY[I] := C_Complex(0);
                end;
                Inc(I);
            end;
            CV1 := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                CV1 := C_Add(CV1,C_Mul(CX[i_],CY[i_]));
            end;
            XCDot(CX, CY, N, Temp, CV2, CV2Err);
            ApproxErrors := ApproxErrors or AP_FP_Greater(AbsComplex(C_Sub(CV1,CV2)),ApproxThreshold);
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // test of precision: real
    //
    N := 50000;
    SetLength(RX, N);
    SetLength(RY, N);
    SetLength(Temp, N);
    Pass:=0;
    while Pass<=PassCount-1 do
    begin
        Assert(N mod 2=0);
        
        //
        // First test: X + X + ... + X - X - X - ... - X = 1*X
        //
        S := Exp(Max(Pass, 50));
        if (Pass=PassCount-1) and (Pass>1) then
        begin
            S := MaxRealNumber;
        end;
        RY[0] := (2*RandomReal-1)*S*Sqrt(2*RandomReal);
        I:=1;
        while I<=N-1 do
        begin
            RY[I] := RY[0];
            Inc(I);
        end;
        I:=0;
        while I<=N div 2-1 do
        begin
            RX[I] := 1;
            Inc(I);
        end;
        I:=N div 2;
        while I<=N-2 do
        begin
            RX[I] := -1;
            Inc(I);
        end;
        RX[N-1] := 0;
        XDot(RX, RY, N, Temp, RV2, RV2Err);
        ExactnessErrors := ExactnessErrors or AP_FP_Less(RV2Err,0);
        ExactnessErrors := ExactnessErrors or AP_FP_Greater(RV2Err,4*MachineEpsilon*AbsReal(RY[0]));
        ExactnessErrors := ExactnessErrors or AP_FP_Greater(AbsReal(RV2-RY[0]),RV2Err);
        
        //
        // First test: X + X + ... + X = N*X
        //
        S := Exp(Max(Pass, 50));
        if (Pass=PassCount-1) and (Pass>1) then
        begin
            S := MaxRealNumber;
        end;
        RY[0] := (2*RandomReal-1)*S*Sqrt(2*RandomReal);
        I:=1;
        while I<=N-1 do
        begin
            RY[I] := RY[0];
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            RX[I] := 1;
            Inc(I);
        end;
        XDot(RX, RY, N, Temp, RV2, RV2Err);
        ExactnessErrors := ExactnessErrors or AP_FP_Less(RV2Err,0);
        ExactnessErrors := ExactnessErrors or AP_FP_Greater(RV2Err,4*MachineEpsilon*AbsReal(RY[0])*N);
        ExactnessErrors := ExactnessErrors or AP_FP_Greater(AbsReal(RV2-N*RY[0]),RV2Err);
        Inc(Pass);
    end;
    
    //
    // test of precision: complex
    //
    N := 50000;
    SetLength(CX, N);
    SetLength(CY, N);
    SetLength(Temp, 2*N);
    Pass:=0;
    while Pass<=PassCount-1 do
    begin
        Assert(N mod 2=0);
        
        //
        // First test: X + X + ... + X - X - X - ... - X = 1*X
        //
        S := Exp(Max(Pass, 50));
        if (Pass=PassCount-1) and (Pass>1) then
        begin
            S := MaxRealNumber;
        end;
        CY[0].X := (2*RandomReal-1)*S*Sqrt(2*RandomReal);
        CY[0].Y := (2*RandomReal-1)*S*Sqrt(2*RandomReal);
        I:=1;
        while I<=N-1 do
        begin
            CY[I] := CY[0];
            Inc(I);
        end;
        I:=0;
        while I<=N div 2-1 do
        begin
            CX[I] := C_Complex(1);
            Inc(I);
        end;
        I:=N div 2;
        while I<=N-2 do
        begin
            CX[I] := C_Complex(-1);
            Inc(I);
        end;
        CX[N-1] := C_Complex(0);
        XCDot(CX, CY, N, Temp, CV2, CV2Err);
        ExactnessErrors := ExactnessErrors or AP_FP_Less(CV2Err,0);
        ExactnessErrors := ExactnessErrors or AP_FP_Greater(CV2Err,4*MachineEpsilon*AbsComplex(CY[0]));
        ExactnessErrors := ExactnessErrors or AP_FP_Greater(AbsComplex(C_Sub(CV2,CY[0])),CV2Err);
        
        //
        // First test: X + X + ... + X = N*X
        //
        S := Exp(Max(Pass, 50));
        if (Pass=PassCount-1) and (Pass>1) then
        begin
            S := MaxRealNumber;
        end;
        CY[0] := C_Complex((2*RandomReal-1)*S*Sqrt(2*RandomReal));
        I:=1;
        while I<=N-1 do
        begin
            CY[I] := CY[0];
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            CX[I] := C_Complex(1);
            Inc(I);
        end;
        XCDot(CX, CY, N, Temp, CV2, CV2Err);
        ExactnessErrors := ExactnessErrors or AP_FP_Less(CV2Err,0);
        ExactnessErrors := ExactnessErrors or AP_FP_Greater(CV2Err,4*MachineEpsilon*AbsComplex(CY[0])*N);
        ExactnessErrors := ExactnessErrors or AP_FP_Greater(AbsComplex(C_Sub(CV2,C_MulR(CY[0],N))),CV2Err);
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