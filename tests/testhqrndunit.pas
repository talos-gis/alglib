unit testhqrndunit;
interface
uses Math, Sysutils, Ap, hqrnd;

procedure CalculateMV(const X : TReal1DArray;
     N : AlglibInteger;
     var Mean : Double;
     var MeanS : Double;
     var StdDev : Double;
     var StdDevS : Double);
function TestHQRND(Silent : Boolean):Boolean;
function testhqrndunit_test_silent():Boolean;
function testhqrndunit_test():Boolean;

implementation

procedure UnsetState(var State : HQRNDState);forward;


procedure CalculateMV(const X : TReal1DArray;
     N : AlglibInteger;
     var Mean : Double;
     var MeanS : Double;
     var StdDev : Double;
     var StdDevS : Double);
var
    I : AlglibInteger;
    V : Double;
    V1 : Double;
    V2 : Double;
    Variance : Double;
begin
    Mean := 0;
    MeanS := 1;
    StdDev := 0;
    StdDevS := 1;
    Variance := 0;
    if N<=1 then
    begin
        Exit;
    end;
    
    //
    // Mean
    //
    I:=0;
    while I<=N-1 do
    begin
        Mean := Mean+X[I];
        Inc(I);
    end;
    Mean := Mean/N;
    
    //
    // Variance (using corrected two-pass algorithm)
    //
    if N<>1 then
    begin
        V1 := 0;
        I:=0;
        while I<=N-1 do
        begin
            V1 := V1+AP_Sqr(X[I]-Mean);
            Inc(I);
        end;
        V2 := 0;
        I:=0;
        while I<=N-1 do
        begin
            V2 := V2+(X[I]-Mean);
            Inc(I);
        end;
        V2 := AP_Sqr(V2)/N;
        Variance := (V1-V2)/(N-1);
        if AP_FP_Less(Variance,0) then
        begin
            Variance := 0;
        end;
        StdDev := Sqrt(Variance);
    end;
    
    //
    // Errors
    //
    MeanS := StdDev/Sqrt(N);
    StdDevS := StdDev*Sqrt(2)/Sqrt(N-1);
end;


function TestHQRND(Silent : Boolean):Boolean;
var
    WasErrors : Boolean;
    SampleSize : AlglibInteger;
    SigmaThreshold : Double;
    PassCount : AlglibInteger;
    N : AlglibInteger;
    I : AlglibInteger;
    Pass : AlglibInteger;
    S1 : AlglibInteger;
    S2 : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    R1 : Double;
    R2 : Double;
    X : TReal1DArray;
    Mean : Double;
    MeanS : Double;
    StdDev : Double;
    StdDevS : Double;
    Lambda : Double;
    SeedErrors : Boolean;
    URErrors : Boolean;
    URSigmaErr : Double;
    UIErrors : Boolean;
    UISigmaErr : Double;
    NormErrors : Boolean;
    NormSigmaErr : Double;
    ExpErrors : Boolean;
    ExpSigmaErr : Double;
    State : HQRNDState;
begin
    WasErrors := False;
    SigmaThreshold := 7;
    SampleSize := 100000;
    PassCount := 50;
    SetLength(X, SampleSize-1+1);
    
    //
    // Test seed errors
    //
    SeedErrors := False;
    Pass:=1;
    while Pass<=PassCount do
    begin
        S1 := 1+RandomInteger(32000);
        S2 := 1+RandomInteger(32000);
        UnsetState(State);
        HQRNDSeed(S1, S2, State);
        I1 := HQRNDUniformI(100, State);
        UnsetState(State);
        HQRNDSeed(S1, S2, State);
        I2 := HQRNDUniformI(100, State);
        SeedErrors := SeedErrors or (I1<>I2);
        UnsetState(State);
        HQRNDSeed(S1, S2, State);
        R1 := HQRNDUniformR(State);
        UnsetState(State);
        HQRNDSeed(S1, S2, State);
        R2 := HQRNDUniformR(State);
        SeedErrors := SeedErrors or AP_FP_Neq(R1,R2);
        Inc(Pass);
    end;
    
    //
    // Test HQRNDRandomize() and real uniform generator
    //
    UnsetState(State);
    HQRNDRandomize(State);
    URErrors := False;
    URSigmaErr := 0;
    I:=0;
    while I<=SampleSize-1 do
    begin
        X[I] := HQRNDUniformR(State);
        Inc(I);
    end;
    I:=0;
    while I<=SampleSize-1 do
    begin
        URErrors := URErrors or AP_FP_Less_Eq(X[I],0) or AP_FP_Greater_Eq(X[I],1);
        Inc(I);
    end;
    CalculateMV(X, SampleSize, Mean, MeanS, StdDev, StdDevS);
    if AP_FP_Neq(MeanS,0) then
    begin
        URSigmaErr := Max(URSigmaErr, AbsReal((Mean-0.5)/MeanS));
    end
    else
    begin
        URErrors := True;
    end;
    if AP_FP_Neq(StdDevS,0) then
    begin
        URSigmaErr := Max(URSigmaErr, AbsReal((StdDev-Sqrt(AP_Double(1)/12))/StdDevS));
    end
    else
    begin
        URErrors := True;
    end;
    URErrors := URErrors or AP_FP_Greater(URSigmaErr,SigmaThreshold);
    
    //
    // Test HQRNDRandomize() and integer uniform
    //
    UnsetState(State);
    HQRNDRandomize(State);
    UIErrors := False;
    UISigmaErr := 0;
    N:=2;
    while N<=10 do
    begin
        I:=0;
        while I<=SampleSize-1 do
        begin
            X[I] := HQRNDUniformI(N, State);
            Inc(I);
        end;
        I:=0;
        while I<=SampleSize-1 do
        begin
            UIErrors := UIErrors or AP_FP_Less(X[I],0) or AP_FP_Greater_Eq(X[I],N);
            Inc(I);
        end;
        CalculateMV(X, SampleSize, Mean, MeanS, StdDev, StdDevS);
        if AP_FP_Neq(MeanS,0) then
        begin
            UISigmaErr := Max(UISigmaErr, AbsReal((Mean-0.5*(N-1))/MeanS));
        end
        else
        begin
            UIErrors := True;
        end;
        if AP_FP_Neq(StdDevS,0) then
        begin
            UISigmaErr := Max(UISigmaErr, AbsReal((StdDev-Sqrt((AP_Sqr(N)-1)/12))/StdDevS));
        end
        else
        begin
            UIErrors := True;
        end;
        Inc(N);
    end;
    UIErrors := UIErrors or AP_FP_Greater(UISigmaErr,SigmaThreshold);
    
    //
    // Special 'close-to-limit' test on uniformity of integers
    // (straightforward implementation like 'RND mod N' will return
    //  non-uniform numbers for N=2/3*LIMIT)
    //
    UnsetState(State);
    HQRNDRandomize(State);
    UIErrors := False;
    UISigmaErr := 0;
    N := Round(2.0/3.0*2147483563.0);
    I:=0;
    while I<=SampleSize-1 do
    begin
        X[I] := HQRNDUniformI(N, State);
        Inc(I);
    end;
    I:=0;
    while I<=SampleSize-1 do
    begin
        UIErrors := UIErrors or AP_FP_Less(X[I],0) or AP_FP_Greater_Eq(X[I],N);
        Inc(I);
    end;
    CalculateMV(X, SampleSize, Mean, MeanS, StdDev, StdDevS);
    if AP_FP_Neq(MeanS,0) then
    begin
        UISigmaErr := Max(UISigmaErr, AbsReal((Mean-0.5*(N-1))/MeanS));
    end
    else
    begin
        UIErrors := True;
    end;
    if AP_FP_Neq(StdDevS,0) then
    begin
        UISigmaErr := Max(UISigmaErr, AbsReal((StdDev-Sqrt((AP_Sqr(N)-1)/12))/StdDevS));
    end
    else
    begin
        UIErrors := True;
    end;
    UIErrors := UIErrors or AP_FP_Greater(UISigmaErr,SigmaThreshold);
    
    //
    // Test normal
    //
    UnsetState(State);
    HQRNDRandomize(State);
    NormErrors := False;
    NormSigmaErr := 0;
    I := 0;
    while I<SampleSize do
    begin
        HQRNDNormal2(State, R1, R2);
        X[I] := R1;
        if I+1<SampleSize then
        begin
            X[I+1] := R2;
        end;
        I := I+2;
    end;
    CalculateMV(X, SampleSize, Mean, MeanS, StdDev, StdDevS);
    if AP_FP_Neq(MeanS,0) then
    begin
        NormSigmaErr := Max(NormSigmaErr, AbsReal((Mean-0)/MeanS));
    end
    else
    begin
        NormErrors := True;
    end;
    if AP_FP_Neq(StdDevS,0) then
    begin
        NormSigmaErr := Max(NormSigmaErr, AbsReal((StdDev-1)/StdDevS));
    end
    else
    begin
        NormErrors := True;
    end;
    NormErrors := NormErrors or AP_FP_Greater(NormSigmaErr,SigmaThreshold);
    
    //
    // Test exponential
    //
    UnsetState(State);
    HQRNDRandomize(State);
    ExpErrors := False;
    ExpSigmaErr := 0;
    Lambda := 2+5*RandomReal;
    I:=0;
    while I<=SampleSize-1 do
    begin
        X[I] := HQRNDExponential(Lambda, State);
        Inc(I);
    end;
    I:=0;
    while I<=SampleSize-1 do
    begin
        UIErrors := UIErrors or AP_FP_Less(X[I],0);
        Inc(I);
    end;
    CalculateMV(X, SampleSize, Mean, MeanS, StdDev, StdDevS);
    if AP_FP_Neq(MeanS,0) then
    begin
        ExpSigmaErr := Max(ExpSigmaErr, AbsReal((Mean-1.0/Lambda)/MeanS));
    end
    else
    begin
        ExpErrors := True;
    end;
    if AP_FP_Neq(StdDevS,0) then
    begin
        ExpSigmaErr := Max(ExpSigmaErr, AbsReal((StdDev-1.0/Lambda)/StdDevS));
    end
    else
    begin
        ExpErrors := True;
    end;
    ExpErrors := ExpErrors or AP_FP_Greater(ExpSigmaErr,SigmaThreshold);
    
    //
    // Final report
    //
    WasErrors := SeedErrors or URErrors or UIErrors or NormErrors or ExpErrors;
    if  not Silent then
    begin
        Write(Format('RNG TEST'#13#10'',[]));
        Write(Format('SEED TEST:                               ',[]));
        if  not SeedErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('UNIFORM CONTINUOUS:                      ',[]));
        if  not URErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('UNIFORM INTEGER:                         ',[]));
        if  not UIErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('NORMAL:                                  ',[]));
        if  not NormErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('EXPONENTIAL:                             ',[]));
        if  not ExpErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        if WasErrors then
        begin
            Write(Format('TEST SUMMARY: FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('TEST SUMMARY: PASSED'#13#10'',[]));
        end;
        Write(Format(''#13#10''#13#10'',[]));
    end;
    Result :=  not WasErrors;
end;


(*************************************************************************
Unsets HQRNDState structure
*************************************************************************)
procedure UnsetState(var State : HQRNDState);
begin
    State.S1 := 0;
    State.S2 := 0;
    State.V := 0;
    State.MagicV := 0;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testhqrndunit_test_silent():Boolean;
begin
    Result := TestHQRND(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testhqrndunit_test():Boolean;
begin
    Result := TestHQRND(False);
end;


end.