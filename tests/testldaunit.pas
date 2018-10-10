unit testldaunit;
interface
uses Math, Sysutils, Ap, blas, rotations, tdevd, sblas, reflections, tridiagonal, sevd, cholesky, spdinverse, lda;

function TestLDA(Silent : Boolean):Boolean;
function testldaunit_test_silent():Boolean;
function testldaunit_test():Boolean;

implementation

procedure GenSimpleSet(NFeatures : AlglibInteger;
     NClasses : AlglibInteger;
     NSamples : AlglibInteger;
     Axis : AlglibInteger;
     var XY : TReal2DArray);forward;
procedure GenDeg1Set(NFeatures : AlglibInteger;
     NClasses : AlglibInteger;
     NSamples : AlglibInteger;
     Axis : AlglibInteger;
     var XY : TReal2DArray);forward;
function GenerateNormal(Mean : Double; Sigma : Double):Double;forward;
function TestWN(const XY : TReal2DArray;
     const WN : TReal2DArray;
     NS : AlglibInteger;
     NF : AlglibInteger;
     NC : AlglibInteger;
     NDeg : AlglibInteger):Boolean;forward;
function CalcJ(NF : AlglibInteger;
     const ST : TReal2DArray;
     const SW : TReal2DArray;
     const W : TReal1DArray;
     var P : Double;
     var Q : Double):Double;forward;
procedure FisherS(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NFeatures : AlglibInteger;
     NClasses : AlglibInteger;
     var ST : TReal2DArray;
     var SW : TReal2DArray);forward;


function TestLDA(Silent : Boolean):Boolean;
var
    MaxNF : AlglibInteger;
    MaxNS : AlglibInteger;
    MaxNC : AlglibInteger;
    PassCount : AlglibInteger;
    LDANErrors : Boolean;
    LDA1Errors : Boolean;
    WasErrors : Boolean;
    NF : AlglibInteger;
    NC : AlglibInteger;
    NS : AlglibInteger;
    I : AlglibInteger;
    Info : AlglibInteger;
    Pass : AlglibInteger;
    Axis : AlglibInteger;
    XY : TReal2DArray;
    WN : TReal2DArray;
    W1 : TReal1DArray;
begin
    
    //
    // Primary settings
    //
    MaxNF := 10;
    MaxNS := 1000;
    MaxNC := 5;
    PassCount := 1;
    WasErrors := False;
    LDANErrors := False;
    LDA1Errors := False;
    
    //
    // General tests
    //
    NF:=1;
    while NF<=MaxNF do
    begin
        NC:=2;
        while NC<=MaxNC do
        begin
            Pass:=1;
            while Pass<=PassCount do
            begin
                
                //
                // Simple test for LDA-N/LDA-1
                //
                Axis := RandomInteger(NF);
                NS := MaxNS div 2+RandomInteger(MaxNS div 2);
                GenSimpleSet(NF, NC, NS, Axis, XY);
                FisherLDAN(XY, NS, NF, NC, Info, WN);
                if Info<>1 then
                begin
                    LDANErrors := True;
                    Inc(Pass);
                    Continue;
                end;
                LDANErrors := LDANErrors or  not TestWN(XY, WN, NS, NF, NC, 0);
                LDANErrors := LDANErrors or AP_FP_Less_Eq(AbsReal(WN[Axis,0]),0.75);
                FisherLDA(XY, NS, NF, NC, Info, W1);
                I:=0;
                while I<=NF-1 do
                begin
                    LDA1Errors := LDA1Errors or AP_FP_Neq(W1[I],WN[I,0]);
                    Inc(I);
                end;
                
                //
                // Degenerate test for LDA-N
                //
                if NF>=3 then
                begin
                    NS := MaxNS div 2+RandomInteger(MaxNS div 2);
                    
                    //
                    // there are two duplicate features,
                    // axis is oriented along non-duplicate feature
                    //
                    Axis := RandomInteger(NF-2);
                    GenDeg1Set(NF, NC, NS, Axis, XY);
                    FisherLDAN(XY, NS, NF, NC, Info, WN);
                    if Info<>2 then
                    begin
                        LDANErrors := True;
                        Inc(Pass);
                        Continue;
                    end;
                    LDANErrors := LDANErrors or AP_FP_Less_Eq(WN[Axis,0],0.75);
                    FisherLDA(XY, NS, NF, NC, Info, W1);
                    I:=0;
                    while I<=NF-1 do
                    begin
                        LDA1Errors := LDA1Errors or AP_FP_Neq(W1[I],WN[I,0]);
                        Inc(I);
                    end;
                end;
                Inc(Pass);
            end;
            Inc(NC);
        end;
        Inc(NF);
    end;
    
    //
    // Final report
    //
    WasErrors := LDANErrors or LDA1Errors;
    if  not Silent then
    begin
        Write(Format('LDA TEST'#13#10'',[]));
        Write(Format('FISHER LDA-N:                            ',[]));
        if  not LDANErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('FISHER LDA-1:                            ',[]));
        if  not LDA1Errors then
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
Generates 'simple' set - a sequence of unit 'balls' at (0,0), (1,0), (2,0)
and so on.
*************************************************************************)
procedure GenSimpleSet(NFeatures : AlglibInteger;
     NClasses : AlglibInteger;
     NSamples : AlglibInteger;
     Axis : AlglibInteger;
     var XY : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    C : AlglibInteger;
    V : Double;
begin
    Assert((Axis>=0) and (Axis<NFeatures), 'GenSimpleSet: wrong Axis!');
    SetLength(XY, NSamples-1+1, NFeatures+1);
    I:=0;
    while I<=NSamples-1 do
    begin
        J:=0;
        while J<=NFeatures-1 do
        begin
            XY[I,J] := GenerateNormal(0.0, 1.0);
            Inc(J);
        end;
        C := I mod NClasses;
        XY[I,Axis] := XY[I,Axis]+C;
        XY[I,NFeatures] := C;
        Inc(I);
    end;
end;


(*************************************************************************
Generates 'degenerate' set #1.
NFeatures>=3.
*************************************************************************)
procedure GenDeg1Set(NFeatures : AlglibInteger;
     NClasses : AlglibInteger;
     NSamples : AlglibInteger;
     Axis : AlglibInteger;
     var XY : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    C : AlglibInteger;
    V : Double;
begin
    Assert((Axis>=0) and (Axis<NFeatures), 'GenDeg1Set: wrong Axis!');
    Assert(NFeatures>=3, 'GenDeg1Set: wrong NFeatures!');
    SetLength(XY, NSamples-1+1, NFeatures+1);
    if Axis>=NFeatures-2 then
    begin
        Axis := NFeatures-3;
    end;
    I:=0;
    while I<=NSamples-1 do
    begin
        J:=0;
        while J<=NFeatures-2 do
        begin
            XY[I,J] := GenerateNormal(0.0, 1.0);
            Inc(J);
        end;
        XY[I,NFeatures-1] := XY[I,NFeatures-2];
        C := I mod NClasses;
        XY[I,Axis] := XY[I,Axis]+C;
        XY[I,NFeatures] := C;
        Inc(I);
    end;
end;


(*************************************************************************
Normal random number
*************************************************************************)
function GenerateNormal(Mean : Double; Sigma : Double):Double;
var
    U : Double;
    V : Double;
    S : Double;
    Sum : Double;
begin
    Result := Mean;
    while True do
    begin
        u := (2*RandomInteger(2)-1)*RandomReal;
        v := (2*RandomInteger(2)-1)*RandomReal;
        sum := u*u+v*v;
        if AP_FP_Less(sum,1) and AP_FP_Greater(sum,0) then
        begin
            sum := sqrt(-2*ln(sum)/sum);
            Result := Sigma*u*sum+Mean;
            Exit;
        end;
    end;
end;


(*************************************************************************
Tests WN for correctness
*************************************************************************)
function TestWN(const XY : TReal2DArray;
     const WN : TReal2DArray;
     NS : AlglibInteger;
     NF : AlglibInteger;
     NC : AlglibInteger;
     NDeg : AlglibInteger):Boolean;
var
    ST : TReal2DArray;
    SW : TReal2DArray;
    A : TReal2DArray;
    Z : TReal2DArray;
    TX : TReal1DArray;
    JP : TReal1DArray;
    JQ : TReal1DArray;
    WORK : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    WPrev : Double;
    TOL : Double;
    P : Double;
    Q : Double;
    i_ : AlglibInteger;
begin
    TOL := 10000;
    Result := True;
    FisherS(XY, NS, NF, NC, ST, SW);
    
    //
    // Test for decreasing of J
    //
    SetLength(TX, NF-1+1);
    SetLength(JP, NF-1+1);
    SetLength(JQ, NF-1+1);
    J:=0;
    while J<=NF-1 do
    begin
        for i_ := 0 to NF-1 do
        begin
            TX[i_] := WN[i_,J];
        end;
        V := CalcJ(NF, ST, SW, TX, P, Q);
        JP[J] := P;
        JQ[J] := Q;
        Inc(J);
    end;
    I:=1;
    while I<=NF-1-NDeg do
    begin
        Result := Result and AP_FP_Greater_Eq(JP[I-1]/JQ[I-1],(1-TOL*MachineEpsilon)*JP[I]/JQ[I]);
        Inc(I);
    end;
    I:=NF-1-NDeg+1;
    while I<=NF-1 do
    begin
        Result := Result and AP_FP_Less_Eq(JP[I],TOL*MachineEpsilon*JP[0]);
        Inc(I);
    end;
    
    //
    // Test for J optimality
    //
    for i_ := 0 to NF-1 do
    begin
        TX[i_] := WN[i_,0];
    end;
    V := CalcJ(NF, ST, SW, TX, P, Q);
    I:=0;
    while I<=NF-1 do
    begin
        WPrev := TX[I];
        TX[I] := WPrev+0.01;
        Result := Result and AP_FP_Greater_Eq(V,(1-TOL*MachineEpsilon)*CalcJ(NF, ST, SW, TX, P, Q));
        TX[I] := WPrev-0.01;
        Result := Result and AP_FP_Greater_Eq(V,(1-TOL*MachineEpsilon)*CalcJ(NF, ST, SW, TX, P, Q));
        TX[I] := WPrev;
        Inc(I);
    end;
    
    //
    // Test for linear independence of W
    //
    SetLength(WORK, NF+1);
    SetLength(A, NF-1+1, NF-1+1);
    MatrixMatrixMultiply(WN, 0, NF-1, 0, NF-1, False, WN, 0, NF-1, 0, NF-1, True, 1.0, A, 0, NF-1, 0, NF-1, 0.0, WORK);
    if SMatrixEVD(A, NF, 1, True, TX, Z) then
    begin
        Result := Result and AP_FP_Greater(TX[0],TX[NF-1]*1000*MachineEpsilon);
    end;
    
    //
    // Test for other properties
    //
    J:=0;
    while J<=NF-1 do
    begin
        V := 0.0;
        for i_ := 0 to NF-1 do
        begin
            V := V + WN[i_,J]*WN[i_,J];
        end;
        V := Sqrt(V);
        Result := Result and AP_FP_Less_Eq(AbsReal(V-1),1000*MachineEpsilon);
        V := 0;
        I:=0;
        while I<=NF-1 do
        begin
            V := V+WN[I,J];
            Inc(I);
        end;
        Result := Result and AP_FP_Greater_Eq(V,0);
        Inc(J);
    end;
end;


(*************************************************************************
Calculates J
*************************************************************************)
function CalcJ(NF : AlglibInteger;
     const ST : TReal2DArray;
     const SW : TReal2DArray;
     const W : TReal1DArray;
     var P : Double;
     var Q : Double):Double;
var
    TX : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
begin
    SetLength(TX, NF-1+1);
    I:=0;
    while I<=NF-1 do
    begin
        V := APVDotProduct(@ST[I][0], 0, NF-1, @W[0], 0, NF-1);
        TX[I] := V;
        Inc(I);
    end;
    V := APVDotProduct(@W[0], 0, NF-1, @TX[0], 0, NF-1);
    P := V;
    I:=0;
    while I<=NF-1 do
    begin
        V := APVDotProduct(@SW[I][0], 0, NF-1, @W[0], 0, NF-1);
        TX[I] := V;
        Inc(I);
    end;
    V := APVDotProduct(@W[0], 0, NF-1, @TX[0], 0, NF-1);
    Q := V;
    Result := P/Q;
end;


(*************************************************************************
Calculates ST/SW
*************************************************************************)
procedure FisherS(const XY : TReal2DArray;
     NPoints : AlglibInteger;
     NFeatures : AlglibInteger;
     NClasses : AlglibInteger;
     var ST : TReal2DArray;
     var SW : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    C : TInteger1DArray;
    Mu : TReal1DArray;
    MuC : TReal2DArray;
    NC : TInteger1DArray;
    TF : TReal1DArray;
    WORK : TReal1DArray;
begin
    
    //
    // Prepare temporaries
    //
    SetLength(TF, NFeatures-1+1);
    SetLength(WORK, NFeatures+1);
    
    //
    // Convert class labels from reals to integers (just for convenience)
    //
    SetLength(C, NPoints-1+1);
    I:=0;
    while I<=NPoints-1 do
    begin
        C[I] := Round(XY[I,NFeatures]);
        Inc(I);
    end;
    
    //
    // Calculate class sizes and means
    //
    SetLength(Mu, NFeatures-1+1);
    SetLength(MuC, NClasses-1+1, NFeatures-1+1);
    SetLength(NC, NClasses-1+1);
    J:=0;
    while J<=NFeatures-1 do
    begin
        Mu[J] := 0;
        Inc(J);
    end;
    I:=0;
    while I<=NClasses-1 do
    begin
        NC[I] := 0;
        J:=0;
        while J<=NFeatures-1 do
        begin
            MuC[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=NPoints-1 do
    begin
        APVAdd(@Mu[0], 0, NFeatures-1, @XY[I][0], 0, NFeatures-1);
        APVAdd(@MuC[C[I]][0], 0, NFeatures-1, @XY[I][0], 0, NFeatures-1);
        NC[C[I]] := NC[C[I]]+1;
        Inc(I);
    end;
    I:=0;
    while I<=NClasses-1 do
    begin
        V := AP_Double(1)/NC[I];
        APVMul(@MuC[I][0], 0, NFeatures-1, V);
        Inc(I);
    end;
    V := AP_Double(1)/NPoints;
    APVMul(@Mu[0], 0, NFeatures-1, V);
    
    //
    // Create ST matrix
    //
    SetLength(ST, NFeatures-1+1, NFeatures-1+1);
    I:=0;
    while I<=NFeatures-1 do
    begin
        J:=0;
        while J<=NFeatures-1 do
        begin
            ST[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    K:=0;
    while K<=NPoints-1 do
    begin
        APVMove(@TF[0], 0, NFeatures-1, @XY[K][0], 0, NFeatures-1);
        APVSub(@TF[0], 0, NFeatures-1, @Mu[0], 0, NFeatures-1);
        I:=0;
        while I<=NFeatures-1 do
        begin
            V := TF[I];
            APVAdd(@ST[I][0], 0, NFeatures-1, @TF[0], 0, NFeatures-1, V);
            Inc(I);
        end;
        Inc(K);
    end;
    
    //
    // Create SW matrix
    //
    SetLength(SW, NFeatures-1+1, NFeatures-1+1);
    I:=0;
    while I<=NFeatures-1 do
    begin
        J:=0;
        while J<=NFeatures-1 do
        begin
            SW[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    K:=0;
    while K<=NPoints-1 do
    begin
        APVMove(@TF[0], 0, NFeatures-1, @XY[K][0], 0, NFeatures-1);
        APVSub(@TF[0], 0, NFeatures-1, @MuC[C[K]][0], 0, NFeatures-1);
        I:=0;
        while I<=NFeatures-1 do
        begin
            V := TF[I];
            APVAdd(@SW[I][0], 0, NFeatures-1, @TF[0], 0, NFeatures-1, V);
            Inc(I);
        end;
        Inc(K);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testldaunit_test_silent():Boolean;
begin
    Result := TestLDA(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testldaunit_test():Boolean;
begin
    Result := TestLDA(False);
end;


end.