unit testforestunit;
interface
uses Math, Sysutils, Ap, tsort, descriptivestatistics, bdss, dforest;

function TestForest(Silent : Boolean):Boolean;
function testforestunit_test_silent():Boolean;
function testforestunit_test():Boolean;

implementation

procedure TestProcessing(var Err : Boolean);forward;
procedure BasicTest1(NVars : AlglibInteger;
     NClasses : AlglibInteger;
     PassCount : AlglibInteger;
     var Err : Boolean);forward;
procedure BasicTest2(var Err : Boolean);forward;
procedure BasicTest3(var Err : Boolean);forward;
procedure BasicTest4(var Err : Boolean);forward;
procedure BasicTest5(var Err : Boolean);forward;
function RNormal():Double;forward;
function RSphere(var XY : TReal2DArray;
     N : AlglibInteger;
     I : AlglibInteger):Double;forward;
procedure UnsetDF(var DF : DecisionForest);forward;


function TestForest(Silent : Boolean):Boolean;
var
    NCMax : AlglibInteger;
    NVMax : AlglibInteger;
    PassCount : AlglibInteger;
    NVars : AlglibInteger;
    NClasses : AlglibInteger;
    WasErrors : Boolean;
    BasicErrors : Boolean;
    ProcErrors : Boolean;
    I : AlglibInteger;
    J : AlglibInteger;
begin
    
    //
    // Primary settings
    //
    NVMax := 4;
    NCMax := 3;
    PassCount := 10;
    BasicErrors := False;
    ProcErrors := False;
    WasErrors := False;
    
    //
    // Tests
    //
    TestProcessing(ProcErrors);
    NVars:=1;
    while NVars<=NVMax do
    begin
        NClasses:=1;
        while NClasses<=NCMax do
        begin
            BasicTest1(NVars, NClasses, PassCount, BasicErrors);
            Inc(NClasses);
        end;
        Inc(NVars);
    end;
    BasicTest2(BasicErrors);
    BasicTest3(BasicErrors);
    BasicTest4(BasicErrors);
    BasicTest5(BasicErrors);
    
    //
    // Final report
    //
    WasErrors := BasicErrors or ProcErrors;
    if  not Silent then
    begin
        Write(Format('RANDOM FOREST TEST'#13#10'',[]));
        Write(Format('TOTAL RESULTS:                           ',[]));
        if  not WasErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('* PROCESSING FUNCTIONS:                  ',[]));
        if  not ProcErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('* BASIC TESTS:                           ',[]));
        if  not BasicErrors then
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
Processing functions test
*************************************************************************)
procedure TestProcessing(var Err : Boolean);
var
    NVars : AlglibInteger;
    NClasses : AlglibInteger;
    NSample : AlglibInteger;
    NTrees : AlglibInteger;
    NFeatures : AlglibInteger;
    Flags : AlglibInteger;
    DF1 : DecisionForest;
    DF2 : DecisionForest;
    NPoints : AlglibInteger;
    XY : TReal2DArray;
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    AllSame : Boolean;
    RLen : AlglibInteger;
    Info : AlglibInteger;
    Rep : DFReport;
    X1 : TReal1DArray;
    X2 : TReal1DArray;
    Y1 : TReal1DArray;
    Y2 : TReal1DArray;
    RA : TReal1DArray;
    RA2 : TReal1DArray;
    V : Double;
begin
    PassCount := 100;
    
    //
    // Main cycle
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        
        //
        // initialize parameters
        //
        NVars := 1+RandomInteger(5);
        NClasses := 1+RandomInteger(3);
        NTrees := 1+RandomInteger(4);
        NFeatures := 1+RandomInteger(NVars);
        Flags := 0;
        if AP_FP_Greater(RandomReal,0.5) then
        begin
            Flags := Flags+2;
        end;
        
        //
        // Initialize arrays and data
        //
        NPoints := 10+RandomInteger(50);
        NSample := Max(10, RandomInteger(NPoints));
        SetLength(X1, NVars-1+1);
        SetLength(X2, NVars-1+1);
        SetLength(Y1, NClasses-1+1);
        SetLength(Y2, NClasses-1+1);
        SetLength(XY, NPoints-1+1, NVars+1);
        I:=0;
        while I<=NPoints-1 do
        begin
            J:=0;
            while J<=NVars-1 do
            begin
                if J mod 2=0 then
                begin
                    XY[I,J] := 2*RandomReal-1;
                end
                else
                begin
                    XY[I,J] := RandomInteger(2);
                end;
                Inc(J);
            end;
            if NClasses=1 then
            begin
                XY[I,NVars] := 2*RandomReal-1;
            end
            else
            begin
                XY[I,NVars] := RandomInteger(NClasses);
            end;
            Inc(I);
        end;
        
        //
        // create forest
        //
        DFBuildInternal(XY, NPoints, NVars, NClasses, NTrees, NSample, NFeatures, Flags, Info, DF1, Rep);
        if Info<=0 then
        begin
            Err := True;
            Exit;
        end;
        
        //
        // Same inputs leads to same outputs
        //
        I:=0;
        while I<=NVars-1 do
        begin
            X1[I] := 2*RandomReal-1;
            X2[I] := X1[I];
            Inc(I);
        end;
        I:=0;
        while I<=NClasses-1 do
        begin
            Y1[I] := 2*RandomReal-1;
            Y2[I] := 2*RandomReal-1;
            Inc(I);
        end;
        DFProcess(DF1, X1, Y1);
        DFProcess(DF1, X2, Y2);
        AllSame := True;
        I:=0;
        while I<=NClasses-1 do
        begin
            AllSame := AllSame and AP_FP_Eq(Y1[I],Y2[I]);
            Inc(I);
        end;
        Err := Err or  not AllSame;
        
        //
        // Same inputs on original forest leads to same outputs
        // on copy created using DFCopy
        //
        UnsetDF(DF2);
        DFCopy(DF1, DF2);
        I:=0;
        while I<=NVars-1 do
        begin
            X1[I] := 2*RandomReal-1;
            X2[I] := X1[I];
            Inc(I);
        end;
        I:=0;
        while I<=NClasses-1 do
        begin
            Y1[I] := 2*RandomReal-1;
            Y2[I] := 2*RandomReal-1;
            Inc(I);
        end;
        DFProcess(DF1, X1, Y1);
        DFProcess(DF2, X2, Y2);
        AllSame := True;
        I:=0;
        while I<=NClasses-1 do
        begin
            AllSame := AllSame and AP_FP_Eq(Y1[I],Y2[I]);
            Inc(I);
        end;
        Err := Err or  not AllSame;
        
        //
        // Same inputs on original forest leads to same outputs
        // on copy created using DFSerialize
        //
        UnsetDF(DF2);
        SetLength(RA, 0+1);
        RA[0] := 0;
        RLen := 0;
        DFSerialize(DF1, RA, RLen);
        SetLength(RA2, RLen-1+1);
        I:=0;
        while I<=RLen-1 do
        begin
            RA2[I] := RA[I];
            Inc(I);
        end;
        DFUnserialize(RA2, DF2);
        I:=0;
        while I<=NVars-1 do
        begin
            X1[I] := 2*RandomReal-1;
            X2[I] := X1[I];
            Inc(I);
        end;
        I:=0;
        while I<=NClasses-1 do
        begin
            Y1[I] := 2*RandomReal-1;
            Y2[I] := 2*RandomReal-1;
            Inc(I);
        end;
        DFProcess(DF1, X1, Y1);
        DFProcess(DF2, X2, Y2);
        AllSame := True;
        I:=0;
        while I<=NClasses-1 do
        begin
            AllSame := AllSame and AP_FP_Eq(Y1[I],Y2[I]);
            Inc(I);
        end;
        Err := Err or  not AllSame;
        
        //
        // Normalization properties
        //
        if NClasses>1 then
        begin
            I:=0;
            while I<=NVars-1 do
            begin
                X1[I] := 2*RandomReal-1;
                Inc(I);
            end;
            DFProcess(DF1, X1, Y1);
            V := 0;
            I:=0;
            while I<=NClasses-1 do
            begin
                V := V+Y1[I];
                Err := Err or AP_FP_Less(Y1[I],0);
                Inc(I);
            end;
            Err := Err or AP_FP_Greater(AbsReal(V-1),1000*MachineEpsilon);
        end;
        Inc(Pass);
    end;
end;


(*************************************************************************
Basic test:  one-tree forest built using full sample must remember all the
training cases
*************************************************************************)
procedure BasicTest1(NVars : AlglibInteger;
     NClasses : AlglibInteger;
     PassCount : AlglibInteger;
     var Err : Boolean);
var
    Pass : AlglibInteger;
    XY : TReal2DArray;
    NPoints : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    S : Double;
    Info : AlglibInteger;
    DF : DecisionForest;
    X : TReal1DArray;
    Y : TReal1DArray;
    Rep : DFReport;
    HasSame : Boolean;
begin
    if NClasses=1 then
    begin
        
        //
        // only classification tasks
        //
        Exit;
    end;
    Pass:=1;
    while Pass<=PassCount do
    begin
        
        //
        // select number of points
        //
        if (Pass<=3) and (PassCount>3) then
        begin
            NPoints := Pass;
        end
        else
        begin
            NPoints := 100+RandomInteger(100);
        end;
        
        //
        // Prepare task
        //
        SetLength(XY, NPoints-1+1, NVars+1);
        SetLength(X, NVars-1+1);
        SetLength(Y, NClasses-1+1);
        I:=0;
        while I<=NPoints-1 do
        begin
            J:=0;
            while J<=NVars-1 do
            begin
                XY[I,J] := 2*RandomReal-1;
                Inc(J);
            end;
            XY[I,NVars] := RandomInteger(NClasses);
            Inc(I);
        end;
        
        //
        // Test
        //
        DFBuildInternal(XY, NPoints, NVars, NClasses, 1, NPoints, 1, 1, Info, DF, Rep);
        if Info<=0 then
        begin
            Err := True;
            Exit;
        end;
        I:=0;
        while I<=NPoints-1 do
        begin
            APVMove(@X[0], 0, NVars-1, @XY[I][0], 0, NVars-1);
            DFProcess(DF, X, Y);
            S := 0;
            J:=0;
            while J<=NClasses-1 do
            begin
                if AP_FP_Less(Y[J],0) then
                begin
                    Err := True;
                    Exit;
                end;
                S := S+Y[J];
                Inc(J);
            end;
            if AP_FP_Greater(AbsReal(S-1),1000*MachineEpsilon) then
            begin
                Err := True;
                Exit;
            end;
            if AP_FP_Greater(AbsReal(Y[Round(XY[I,NVars])]-1),1000*MachineEpsilon) then
            begin
                
                //
                // not an error if there exists such K,J that XY[K,J]=XY[I,J]
                // (may be we just can't distinguish two tied values).
                //
                // definitely error otherwise.
                //
                HasSame := False;
                K:=0;
                while K<=NPoints-1 do
                begin
                    if K<>I then
                    begin
                        J:=0;
                        while J<=NVars-1 do
                        begin
                            if AP_FP_Eq(XY[K,J],XY[I,J]) then
                            begin
                                HasSame := True;
                            end;
                            Inc(J);
                        end;
                    end;
                    Inc(K);
                end;
                if  not HasSame then
                begin
                    Err := True;
                    Exit;
                end;
            end;
            Inc(I);
        end;
        Inc(Pass);
    end;
end;


(*************************************************************************
Basic test:  tests generalization ability on a simple noisy classification
task:
* 0<x<1 - P(class=0)=1
* 1<x<2 - P(class=0)=2-x
* 2<x<3 - P(class=0)=0
*************************************************************************)
procedure BasicTest2(var Err : Boolean);
var
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    XY : TReal2DArray;
    NPoints : AlglibInteger;
    NTrees : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    S : Double;
    Info : AlglibInteger;
    DF : DecisionForest;
    X : TReal1DArray;
    Y : TReal1DArray;
    Rep : DFReport;
    HasSame : Boolean;
begin
    PassCount := 1;
    Pass:=1;
    while Pass<=PassCount do
    begin
        
        //
        // select npoints and ntrees
        //
        NPoints := 3000;
        NTrees := 50;
        
        //
        // Prepare task
        //
        SetLength(XY, NPoints-1+1, 1+1);
        SetLength(X, 0+1);
        SetLength(Y, 1+1);
        I:=0;
        while I<=NPoints-1 do
        begin
            XY[I,0] := 3*RandomReal;
            if AP_FP_Less_Eq(XY[I,0],1) then
            begin
                XY[I,1] := 0;
            end
            else
            begin
                if AP_FP_Less_Eq(XY[I,0],2) then
                begin
                    if AP_FP_Less(RandomReal,XY[I,0]-1) then
                    begin
                        XY[I,1] := 1;
                    end
                    else
                    begin
                        XY[I,1] := 0;
                    end;
                end
                else
                begin
                    XY[I,1] := 1;
                end;
            end;
            Inc(I);
        end;
        
        //
        // Test
        //
        DFBuildInternal(XY, NPoints, 1, 2, NTrees, Round(0.05*NPoints), 1, 0, Info, DF, Rep);
        if Info<=0 then
        begin
            Err := True;
            Exit;
        end;
        X[0] := 0.0;
        while AP_FP_Less_Eq(X[0],3.0) do
        begin
            DFProcess(DF, X, Y);
            
            //
            // Test for basic properties
            //
            S := 0;
            J:=0;
            while J<=1 do
            begin
                if AP_FP_Less(Y[J],0) then
                begin
                    Err := True;
                    Exit;
                end;
                S := S+Y[J];
                Inc(J);
            end;
            if AP_FP_Greater(AbsReal(S-1),1000*MachineEpsilon) then
            begin
                Err := True;
                Exit;
            end;
            
            //
            // test for good correlation with results
            //
            if AP_FP_Less(X[0],1) then
            begin
                Err := Err or AP_FP_Less(Y[0],0.8);
            end;
            if AP_FP_Greater_Eq(X[0],1) and AP_FP_Less_Eq(X[0],2) then
            begin
                Err := Err or AP_FP_Greater(AbsReal(Y[1]-(X[0]-1)),0.5);
            end;
            if AP_FP_Greater(X[0],2) then
            begin
                Err := Err or AP_FP_Less(Y[1],0.8);
            end;
            X[0] := X[0]+0.01;
        end;
        Inc(Pass);
    end;
end;


(*************************************************************************
Basic test:  tests  generalization ability on a simple classification task
(no noise):
* |x|<1, |y|<1
* x^2+y^2<=0.25 - P(class=0)=1
* x^2+y^2>0.25  - P(class=0)=0
*************************************************************************)
procedure BasicTest3(var Err : Boolean);
var
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    XY : TReal2DArray;
    NPoints : AlglibInteger;
    NTrees : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    S : Double;
    Info : AlglibInteger;
    DF : DecisionForest;
    X : TReal1DArray;
    Y : TReal1DArray;
    Rep : DFReport;
    TestGridSize : AlglibInteger;
    R : Double;
begin
    PassCount := 1;
    TestGridSize := 50;
    Pass:=1;
    while Pass<=PassCount do
    begin
        
        //
        // select npoints and ntrees
        //
        NPoints := 2000;
        NTrees := 100;
        
        //
        // Prepare task
        //
        SetLength(XY, NPoints-1+1, 2+1);
        SetLength(X, 1+1);
        SetLength(Y, 1+1);
        I:=0;
        while I<=NPoints-1 do
        begin
            XY[I,0] := 2*RandomReal-1;
            XY[I,1] := 2*RandomReal-1;
            if AP_FP_Less_Eq(AP_Sqr(XY[I,0])+AP_Sqr(XY[I,1]),0.25) then
            begin
                XY[I,2] := 0;
            end
            else
            begin
                XY[I,2] := 1;
            end;
            Inc(I);
        end;
        
        //
        // Test
        //
        DFBuildInternal(XY, NPoints, 2, 2, NTrees, Round(0.1*NPoints), 1, 0, Info, DF, Rep);
        if Info<=0 then
        begin
            Err := True;
            Exit;
        end;
        I:=-TestGridSize div 2;
        while I<=TestGridSize div 2 do
        begin
            J:=-TestGridSize div 2;
            while J<=TestGridSize div 2 do
            begin
                X[0] := AP_Double(I)/(TestGridSize div 2);
                X[1] := AP_Double(J)/(TestGridSize div 2);
                DFProcess(DF, X, Y);
                
                //
                // Test for basic properties
                //
                S := 0;
                K:=0;
                while K<=1 do
                begin
                    if AP_FP_Less(Y[K],0) then
                    begin
                        Err := True;
                        Exit;
                    end;
                    S := S+Y[K];
                    Inc(K);
                end;
                if AP_FP_Greater(AbsReal(S-1),1000*MachineEpsilon) then
                begin
                    Err := True;
                    Exit;
                end;
                
                //
                // test for good correlation with results
                //
                R := Sqrt(AP_Sqr(X[0])+AP_Sqr(X[1]));
                if AP_FP_Less(R,0.5*0.5) then
                begin
                    Err := Err or AP_FP_Less(Y[0],0.6);
                end;
                if AP_FP_Greater(R,0.5*1.5) then
                begin
                    Err := Err or AP_FP_Less(Y[1],0.6);
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        Inc(Pass);
    end;
end;


(*************************************************************************
Basic test: simple regression task without noise:
* |x|<1, |y|<1
* F(x,y) = x^2+y
*************************************************************************)
procedure BasicTest4(var Err : Boolean);
var
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    XY : TReal2DArray;
    NPoints : AlglibInteger;
    NTrees : AlglibInteger;
    NS : AlglibInteger;
    StrongC : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    S : Double;
    Info : AlglibInteger;
    DF : DecisionForest;
    DF2 : DecisionForest;
    X : TReal1DArray;
    Y : TReal1DArray;
    Rep : DFReport;
    Rep2 : DFReport;
    TestGridSize : AlglibInteger;
    MaxErr : Double;
    MaxErr2 : Double;
    AvgErr : Double;
    AvgErr2 : Double;
    Cnt : AlglibInteger;
    EY : Double;
begin
    PassCount := 1;
    TestGridSize := 50;
    Pass:=1;
    while Pass<=PassCount do
    begin
        
        //
        // select npoints and ntrees
        //
        NPoints := 5000;
        NTrees := 100;
        NS := Round(0.1*NPoints);
        StrongC := 1;
        
        //
        // Prepare task
        //
        SetLength(XY, NPoints-1+1, 2+1);
        SetLength(X, 1+1);
        SetLength(Y, 0+1);
        I:=0;
        while I<=NPoints-1 do
        begin
            XY[I,0] := 2*RandomReal-1;
            XY[I,1] := 2*RandomReal-1;
            XY[I,2] := AP_Sqr(XY[I,0])+XY[I,1];
            Inc(I);
        end;
        
        //
        // Test
        //
        DFBuildInternal(XY, NPoints, 2, 1, NTrees, NS, 1, 0, Info, DF, Rep);
        if Info<=0 then
        begin
            Err := True;
            Exit;
        end;
        DFBuildInternal(XY, NPoints, 2, 1, NTrees, NS, 1, StrongC, Info, DF2, Rep2);
        if Info<=0 then
        begin
            Err := True;
            Exit;
        end;
        MaxErr := 0;
        MaxErr2 := 0;
        AvgErr := 0;
        AvgErr2 := 0;
        Cnt := 0;
        I:=Round(-0.7*TestGridSize/2);
        while I<=Round(0.7*TestGridSize/2) do
        begin
            J:=Round(-0.7*TestGridSize/2);
            while J<=Round(0.7*TestGridSize/2) do
            begin
                X[0] := AP_Double(I)/(TestGridSize div 2);
                X[1] := AP_Double(J)/(TestGridSize div 2);
                EY := AP_Sqr(X[0])+X[1];
                DFProcess(DF, X, Y);
                MaxErr := Max(MaxErr, AbsReal(Y[0]-EY));
                AvgErr := AvgErr+AbsReal(Y[0]-EY);
                DFProcess(DF2, X, Y);
                MaxErr2 := Max(MaxErr2, AbsReal(Y[0]-EY));
                AvgErr2 := AvgErr2+AbsReal(Y[0]-EY);
                Cnt := Cnt+1;
                Inc(J);
            end;
            Inc(I);
        end;
        AvgErr := AvgErr/Cnt;
        AvgErr2 := AvgErr2/Cnt;
        Err := Err or AP_FP_Greater(MaxErr,0.2);
        Err := Err or AP_FP_Greater(MaxErr2,0.2);
        Err := Err or AP_FP_Greater(AvgErr,0.1);
        Err := Err or AP_FP_Greater(AvgErr2,0.1);
        Inc(Pass);
    end;
end;


(*************************************************************************
Basic test: extended variable selection leads to better results.

Next task CAN be solved without EVS but it is very unlikely. With EVS
it can be easily and exactly solved.

Task matrix:
    1 0 0 0 ... 0   0
    0 1 0 0 ... 0   1
    0 0 1 0 ... 0   2
    0 0 0 1 ... 0   3
    0 0 0 0 ... 1   N-1
*************************************************************************)
procedure BasicTest5(var Err : Boolean);
var
    XY : TReal2DArray;
    NVars : AlglibInteger;
    NPoints : AlglibInteger;
    NFeatures : AlglibInteger;
    NSample : AlglibInteger;
    NTrees : AlglibInteger;
    EVS : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    EFlag : Boolean;
    Info : AlglibInteger;
    DF : DecisionForest;
    X : TReal1DArray;
    Y : TReal1DArray;
    Rep : DFReport;
begin
    
    //
    // select npoints and ntrees
    //
    NPoints := 50;
    NVars := NPoints;
    NTrees := 1;
    NSample := NPoints;
    EVS := 2;
    NFeatures := 1;
    
    //
    // Prepare task
    //
    SetLength(XY, NPoints-1+1, NVars+1);
    SetLength(X, NVars-1+1);
    SetLength(Y, 0+1);
    I:=0;
    while I<=NPoints-1 do
    begin
        J:=0;
        while J<=NVars-1 do
        begin
            XY[I,J] := 0;
            Inc(J);
        end;
        XY[I,I] := 1;
        XY[I,NVars] := I;
        Inc(I);
    end;
    
    //
    // Without EVS
    //
    DFBuildInternal(XY, NPoints, NVars, 1, NTrees, NSample, NFeatures, 0, Info, DF, Rep);
    if Info<=0 then
    begin
        Err := True;
        Exit;
    end;
    EFlag := False;
    I:=0;
    while I<=NPoints-1 do
    begin
        APVMove(@X[0], 0, NVars-1, @XY[I][0], 0, NVars-1);
        DFProcess(DF, X, Y);
        if AP_FP_Greater(AbsReal(Y[0]-XY[I,NVars]),1000*MachineEpsilon) then
        begin
            EFlag := True;
        end;
        Inc(I);
    end;
    if  not EFlag then
    begin
        Err := True;
        Exit;
    end;
    
    //
    // With EVS
    //
    DFBuildInternal(XY, NPoints, NVars, 1, NTrees, NSample, NFeatures, EVS, Info, DF, Rep);
    if Info<=0 then
    begin
        Err := True;
        Exit;
    end;
    EFlag := False;
    I:=0;
    while I<=NPoints-1 do
    begin
        APVMove(@X[0], 0, NVars-1, @XY[I][0], 0, NVars-1);
        DFProcess(DF, X, Y);
        if AP_FP_Greater(AbsReal(Y[0]-XY[I,NVars]),1000*MachineEpsilon) then
        begin
            EFlag := True;
        end;
        Inc(I);
    end;
    if EFlag then
    begin
        Err := True;
        Exit;
    end;
end;


(*************************************************************************
Random normal number
*************************************************************************)
function RNormal():Double;
var
    U : Double;
    V : Double;
    S : Double;
    X1 : Double;
    X2 : Double;
begin
    while True do
    begin
        U := 2*RandomReal-1;
        V := 2*RandomReal-1;
        S := AP_Sqr(u)+AP_Sqr(v);
        if AP_FP_Greater(S,0) and AP_FP_Less(S,1) then
        begin
            S := Sqrt(-2*Ln(S)/S);
            X1 := U*S;
            X2 := V*S;
            Break;
        end;
    end;
    Result := X1;
end;


(*************************************************************************
Random point from sphere
*************************************************************************)
function RSphere(var XY : TReal2DArray;
     N : AlglibInteger;
     I : AlglibInteger):Double;
var
    J : AlglibInteger;
    V : Double;
begin
    J:=0;
    while J<=N-1 do
    begin
        XY[I,J] := RNormal;
        Inc(J);
    end;
    V := APVDotProduct(@XY[I][0], 0, N-1, @XY[I][0], 0, N-1);
    V := RandomReal/Sqrt(V);
    APVMul(@XY[I][0], 0, N-1, V);
end;


(*************************************************************************
Unsets DF
*************************************************************************)
procedure UnsetDF(var DF : DecisionForest);
var
    XY : TReal2DArray;
    Info : AlglibInteger;
    Rep : DFReport;
begin
    SetLength(XY, 0+1, 1+1);
    XY[0,0] := 0;
    XY[0,1] := 0;
    DFBuildInternal(XY, 1, 1, 1, 1, 1, 1, 0, Info, DF, Rep);
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testforestunit_test_silent():Boolean;
begin
    Result := TestForest(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testforestunit_test():Boolean;
begin
    Result := TestForest(False);
end;


end.