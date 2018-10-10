unit testbdssunit;
interface
uses Math, Sysutils, Ap, tsort, descriptivestatistics, bdss;

function TestBDSS(Silent : Boolean):Boolean;
function testbdssunit_test_silent():Boolean;
function testbdssunit_test():Boolean;

implementation

procedure Unset2D(var A : TComplex2DArray);forward;
procedure Unset1D(var A : TReal1DArray);forward;
procedure Unset1DI(var A : TInteger1DArray);forward;
procedure TestSortResults(const ASorted : TReal1DArray;
     const P1 : TInteger1DArray;
     const P2 : TInteger1DArray;
     const AOriginal : TReal1DArray;
     N : AlglibInteger;
     var WasErrors : Boolean);forward;


(*************************************************************************
Testing BDSS operations
*************************************************************************)
function TestBDSS(Silent : Boolean):Boolean;
var
    N : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    MaxN : AlglibInteger;
    MaxNQ : AlglibInteger;
    A : TReal1DArray;
    A0 : TReal1DArray;
    AT : TReal1DArray;
    P : TReal2DArray;
    Thresholds : TReal1DArray;
    NI : AlglibInteger;
    C : TInteger1DArray;
    P1 : TInteger1DArray;
    P2 : TInteger1DArray;
    Ties : TInteger1DArray;
    PT1 : TInteger1DArray;
    PT2 : TInteger1DArray;
    TieCount : AlglibInteger;
    C1 : AlglibInteger;
    C0 : AlglibInteger;
    NC : AlglibInteger;
    Tmp : TReal1DArray;
    PAL : Double;
    PBL : Double;
    PAR : Double;
    PBR : Double;
    CVE : Double;
    CVR : Double;
    Info : AlglibInteger;
    Threshold : Double;
    TieBuf : TInteger1DArray;
    CntBuf : TInteger1DArray;
    RMS : Double;
    CVRMS : Double;
    WasErrors : Boolean;
    TiesErrors : Boolean;
    Split2Errors : Boolean;
    OptimalSplitKErrors : Boolean;
    SplitKErrors : Boolean;
begin
    WasErrors := False;
    TiesErrors := False;
    Split2Errors := False;
    SplitKErrors := False;
    OptimalSplitKErrors := False;
    MaxN := 100;
    MaxNQ := 49;
    PassCount := 10;
    
    //
    // Test ties
    //
    N:=1;
    while N<=MaxN do
    begin
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // untied data, test DSTie
            //
            Unset1DI(P1);
            Unset1DI(P2);
            Unset1DI(PT1);
            Unset1DI(PT2);
            SetLength(A, N-1+1);
            SetLength(A0, N-1+1);
            SetLength(AT, N-1+1);
            SetLength(Tmp, N-1+1);
            A[0] := 2*RandomReal-1;
            Tmp[0] := RandomReal;
            I:=1;
            while I<=N-1 do
            begin
                
                //
                // A is randomly permuted
                //
                A[I] := A[I-1]+0.1*RandomReal+0.1;
                Tmp[I] := RandomReal;
                Inc(I);
            end;
            TagSortFastR(Tmp, A, N);
            I:=0;
            while I<=N-1 do
            begin
                A0[I] := A[I];
                AT[I] := A[I];
                Inc(I);
            end;
            DSTie(A0, N, Ties, TieCount, P1, P2);
            TagSort(AT, N, PT1, PT2);
            I:=0;
            while I<=N-1 do
            begin
                TiesErrors := TiesErrors or (P1[I]<>PT1[I]);
                TiesErrors := TiesErrors or (P2[I]<>PT2[I]);
                Inc(I);
            end;
            TiesErrors := TiesErrors or (TieCount<>N);
            if TieCount=N then
            begin
                I:=0;
                while I<=N do
                begin
                    TiesErrors := TiesErrors or (Ties[I]<>I);
                    Inc(I);
                end;
            end;
            
            //
            // tied data, test DSTie
            //
            Unset1DI(P1);
            Unset1DI(P2);
            Unset1DI(PT1);
            Unset1DI(PT2);
            SetLength(A, N-1+1);
            SetLength(A0, N-1+1);
            SetLength(AT, N-1+1);
            C1 := 0;
            C0 := 0;
            I:=0;
            while I<=N-1 do
            begin
                A[I] := RandomInteger(2);
                if AP_FP_Eq(A[I],0) then
                begin
                    C0 := C0+1;
                end
                else
                begin
                    C1 := C1+1;
                end;
                A0[I] := A[I];
                AT[I] := A[I];
                Inc(I);
            end;
            DSTie(A0, N, Ties, TieCount, P1, P2);
            TagSort(AT, N, PT1, PT2);
            I:=0;
            while I<=N-1 do
            begin
                TiesErrors := TiesErrors or (P1[I]<>PT1[I]);
                TiesErrors := TiesErrors or (P2[I]<>PT2[I]);
                Inc(I);
            end;
            if (C0=0) or (C1=0) then
            begin
                TiesErrors := TiesErrors or (TieCount<>1);
                if TieCount=1 then
                begin
                    TiesErrors := TiesErrors or (Ties[0]<>0);
                    TiesErrors := TiesErrors or (Ties[1]<>N);
                end;
            end
            else
            begin
                TiesErrors := TiesErrors or (TieCount<>2);
                if TieCount=2 then
                begin
                    TiesErrors := TiesErrors or (Ties[0]<>0);
                    TiesErrors := TiesErrors or (Ties[1]<>C0);
                    TiesErrors := TiesErrors or (Ties[2]<>N);
                end;
            end;
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // split-2
    //
    
    //
    // General tests for different N's
    //
    N:=1;
    while N<=MaxN do
    begin
        SetLength(A, N-1+1);
        SetLength(C, N-1+1);
        
        //
        // one-tie test
        //
        if N mod 2=0 then
        begin
            I:=0;
            while I<=N-1 do
            begin
                A[I] := N;
                C[I] := I mod 2;
                Inc(I);
            end;
            DSOptimalSplit2(A, C, N, Info, Threshold, PAL, PBL, PAR, PBR, CVE);
            if Info<>-3 then
            begin
                Split2Errors := True;
                Inc(N);
                Continue;
            end;
        end;
        
        //
        // two-tie test
        //
        
        //
        // test #1
        //
        if N>1 then
        begin
            I:=0;
            while I<=N-1 do
            begin
                A[I] := I div ((N+1) div 2);
                C[I] := I div ((N+1) div 2);
                Inc(I);
            end;
            DSOptimalSplit2(A, C, N, Info, Threshold, PAL, PBL, PAR, PBR, CVE);
            if Info<>1 then
            begin
                Split2Errors := True;
                Inc(N);
                Continue;
            end;
            Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(Threshold-0.5),100*MachineEpsilon);
            Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(PAL-1),100*MachineEpsilon);
            Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(PBL-0),100*MachineEpsilon);
            Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(PAR-0),100*MachineEpsilon);
            Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(PBR-1),100*MachineEpsilon);
        end;
        Inc(N);
    end;
    
    //
    // Special "CREDIT"-test (transparency coefficient)
    //
    N := 110;
    SetLength(A, N-1+1);
    SetLength(C, N-1+1);
    A[0] := 0.000;
    C[0] := 0;
    A[1] := 0.000;
    C[1] := 0;
    A[2] := 0.000;
    C[2] := 0;
    A[3] := 0.000;
    C[3] := 0;
    A[4] := 0.000;
    C[4] := 0;
    A[5] := 0.000;
    C[5] := 0;
    A[6] := 0.000;
    C[6] := 0;
    A[7] := 0.000;
    C[7] := 1;
    A[8] := 0.000;
    C[8] := 0;
    A[9] := 0.000;
    C[9] := 1;
    A[10] := 0.000;
    C[10] := 0;
    A[11] := 0.000;
    C[11] := 0;
    A[12] := 0.000;
    C[12] := 0;
    A[13] := 0.000;
    C[13] := 0;
    A[14] := 0.000;
    C[14] := 0;
    A[15] := 0.000;
    C[15] := 0;
    A[16] := 0.000;
    C[16] := 0;
    A[17] := 0.000;
    C[17] := 0;
    A[18] := 0.000;
    C[18] := 0;
    A[19] := 0.000;
    C[19] := 0;
    A[20] := 0.000;
    C[20] := 0;
    A[21] := 0.000;
    C[21] := 0;
    A[22] := 0.000;
    C[22] := 1;
    A[23] := 0.000;
    C[23] := 0;
    A[24] := 0.000;
    C[24] := 0;
    A[25] := 0.000;
    C[25] := 0;
    A[26] := 0.000;
    C[26] := 0;
    A[27] := 0.000;
    C[27] := 1;
    A[28] := 0.000;
    C[28] := 0;
    A[29] := 0.000;
    C[29] := 1;
    A[30] := 0.000;
    C[30] := 0;
    A[31] := 0.000;
    C[31] := 1;
    A[32] := 0.000;
    C[32] := 0;
    A[33] := 0.000;
    C[33] := 1;
    A[34] := 0.000;
    C[34] := 0;
    A[35] := 0.030;
    C[35] := 0;
    A[36] := 0.030;
    C[36] := 0;
    A[37] := 0.050;
    C[37] := 0;
    A[38] := 0.070;
    C[38] := 1;
    A[39] := 0.110;
    C[39] := 0;
    A[40] := 0.110;
    C[40] := 1;
    A[41] := 0.120;
    C[41] := 0;
    A[42] := 0.130;
    C[42] := 0;
    A[43] := 0.140;
    C[43] := 0;
    A[44] := 0.140;
    C[44] := 0;
    A[45] := 0.140;
    C[45] := 0;
    A[46] := 0.150;
    C[46] := 0;
    A[47] := 0.150;
    C[47] := 0;
    A[48] := 0.170;
    C[48] := 0;
    A[49] := 0.190;
    C[49] := 1;
    A[50] := 0.200;
    C[50] := 0;
    A[51] := 0.200;
    C[51] := 0;
    A[52] := 0.250;
    C[52] := 0;
    A[53] := 0.250;
    C[53] := 0;
    A[54] := 0.260;
    C[54] := 0;
    A[55] := 0.270;
    C[55] := 0;
    A[56] := 0.280;
    C[56] := 0;
    A[57] := 0.310;
    C[57] := 0;
    A[58] := 0.310;
    C[58] := 0;
    A[59] := 0.330;
    C[59] := 0;
    A[60] := 0.330;
    C[60] := 0;
    A[61] := 0.340;
    C[61] := 0;
    A[62] := 0.340;
    C[62] := 0;
    A[63] := 0.370;
    C[63] := 0;
    A[64] := 0.380;
    C[64] := 1;
    A[65] := 0.380;
    C[65] := 0;
    A[66] := 0.410;
    C[66] := 0;
    A[67] := 0.460;
    C[67] := 0;
    A[68] := 0.520;
    C[68] := 0;
    A[69] := 0.530;
    C[69] := 0;
    A[70] := 0.540;
    C[70] := 0;
    A[71] := 0.560;
    C[71] := 0;
    A[72] := 0.560;
    C[72] := 0;
    A[73] := 0.570;
    C[73] := 0;
    A[74] := 0.600;
    C[74] := 0;
    A[75] := 0.600;
    C[75] := 0;
    A[76] := 0.620;
    C[76] := 0;
    A[77] := 0.650;
    C[77] := 0;
    A[78] := 0.660;
    C[78] := 0;
    A[79] := 0.680;
    C[79] := 0;
    A[80] := 0.700;
    C[80] := 0;
    A[81] := 0.750;
    C[81] := 0;
    A[82] := 0.770;
    C[82] := 0;
    A[83] := 0.770;
    C[83] := 0;
    A[84] := 0.770;
    C[84] := 0;
    A[85] := 0.790;
    C[85] := 0;
    A[86] := 0.810;
    C[86] := 0;
    A[87] := 0.840;
    C[87] := 0;
    A[88] := 0.860;
    C[88] := 0;
    A[89] := 0.870;
    C[89] := 0;
    A[90] := 0.890;
    C[90] := 0;
    A[91] := 0.900;
    C[91] := 1;
    A[92] := 0.900;
    C[92] := 0;
    A[93] := 0.910;
    C[93] := 0;
    A[94] := 0.940;
    C[94] := 0;
    A[95] := 0.950;
    C[95] := 0;
    A[96] := 0.952;
    C[96] := 0;
    A[97] := 0.970;
    C[97] := 0;
    A[98] := 0.970;
    C[98] := 0;
    A[99] := 0.980;
    C[99] := 0;
    A[100] := 1.000;
    C[100] := 0;
    A[101] := 1.000;
    C[101] := 0;
    A[102] := 1.000;
    C[102] := 0;
    A[103] := 1.000;
    C[103] := 0;
    A[104] := 1.000;
    C[104] := 0;
    A[105] := 1.020;
    C[105] := 0;
    A[106] := 1.090;
    C[106] := 0;
    A[107] := 1.130;
    C[107] := 0;
    A[108] := 1.840;
    C[108] := 0;
    A[109] := 2.470;
    C[109] := 0;
    DSOptimalSplit2(A, C, N, Info, Threshold, PAL, PBL, PAR, PBR, CVE);
    if Info<>1 then
    begin
        Split2Errors := True;
    end
    else
    begin
        Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(Threshold-0.195),100*MachineEpsilon);
        Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(PAL-0.80),0.02);
        Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(PBL-0.20),0.02);
        Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(PAR-0.97),0.02);
        Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(PBR-0.03),0.02);
    end;
    
    //
    // split-2 fast
    //
    
    //
    // General tests for different N's
    //
    N:=1;
    while N<=MaxN do
    begin
        SetLength(A, N-1+1);
        SetLength(C, N-1+1);
        SetLength(TieBuf, N+1);
        SetLength(CntBuf, 3+1);
        
        //
        // one-tie test
        //
        if N mod 2=0 then
        begin
            I:=0;
            while I<=N-1 do
            begin
                A[I] := N;
                C[I] := I mod 2;
                Inc(I);
            end;
            DSOptimalSplit2Fast(A, C, TieBuf, CntBuf, N, 2, 0.00, Info, Threshold, RMS, CVRMS);
            if Info<>-3 then
            begin
                Split2Errors := True;
                Inc(N);
                Continue;
            end;
        end;
        
        //
        // two-tie test
        //
        
        //
        // test #1
        //
        if N>1 then
        begin
            I:=0;
            while I<=N-1 do
            begin
                A[I] := I div ((N+1) div 2);
                C[I] := I div ((N+1) div 2);
                Inc(I);
            end;
            DSOptimalSplit2Fast(A, C, TieBuf, CntBuf, N, 2, 0.00, Info, Threshold, RMS, CVRMS);
            if Info<>1 then
            begin
                Split2Errors := True;
                Inc(N);
                Continue;
            end;
            Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(Threshold-0.5),100*MachineEpsilon);
            Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(RMS-0),100*MachineEpsilon);
            if N=2 then
            begin
                Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(CVRMS-0.5),100*MachineEpsilon);
            end
            else
            begin
                if N=3 then
                begin
                    Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(CVRMS-Sqrt((2*0+2*0+2*0.25)/6)),100*MachineEpsilon);
                end
                else
                begin
                    Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(CVRMS),100*MachineEpsilon);
                end;
            end;
        end;
        Inc(N);
    end;
    
    //
    // special tests
    //
    N := 10;
    SetLength(A, N-1+1);
    SetLength(C, N-1+1);
    SetLength(TieBuf, N+1);
    SetLength(CntBuf, 2*3-1+1);
    I:=0;
    while I<=N-1 do
    begin
        A[I] := I;
        if I<=N-3 then
        begin
            C[I] := 0;
        end
        else
        begin
            C[I] := I-(N-3);
        end;
        Inc(I);
    end;
    DSOptimalSplit2Fast(A, C, TieBuf, CntBuf, N, 3, 0.00, Info, Threshold, RMS, CVRMS);
    if Info<>1 then
    begin
        Split2Errors := True;
    end
    else
    begin
        Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(Threshold-(N-2.5)),100*MachineEpsilon);
        Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(RMS-Sqrt((0.25+0.25+0.25+0.25)/(3*N))),100*MachineEpsilon);
        Split2Errors := Split2Errors or AP_FP_Greater(AbsReal(CVRMS-Sqrt(AP_Double((1+1+1+1))/(3*N))),100*MachineEpsilon);
    end;
    
    //
    // Optimal split-K
    //
    
    //
    // General tests for different N's
    //
    N:=1;
    while N<=MaxNQ do
    begin
        SetLength(A, N-1+1);
        SetLength(C, N-1+1);
        
        //
        // one-tie test
        //
        if N mod 2=0 then
        begin
            I:=0;
            while I<=N-1 do
            begin
                A[I] := Pass;
                C[I] := I mod 2;
                Inc(I);
            end;
            DSOptimalSplitK(A, C, N, 2, 2+RandomInteger(5), Info, Thresholds, NI, CVE);
            if Info<>-3 then
            begin
                OptimalSplitKErrors := True;
                Inc(N);
                Continue;
            end;
        end;
        
        //
        // two-tie test
        //
        
        //
        // test #1
        //
        if N>1 then
        begin
            C0 := 0;
            C1 := 0;
            I:=0;
            while I<=N-1 do
            begin
                A[I] := I div ((N+1) div 2);
                C[I] := I div ((N+1) div 2);
                if C[I]=0 then
                begin
                    C0 := C0+1;
                end;
                if C[I]=1 then
                begin
                    C1 := C1+1;
                end;
                Inc(I);
            end;
            DSOptimalSplitK(A, C, N, 2, 2+RandomInteger(5), Info, Thresholds, NI, CVE);
            if Info<>1 then
            begin
                OptimalSplitKErrors := True;
                Inc(N);
                Continue;
            end;
            OptimalSplitKErrors := OptimalSplitKErrors or (NI<>2);
            OptimalSplitKErrors := OptimalSplitKErrors or AP_FP_Greater(AbsReal(Thresholds[0]-0.5),100*MachineEpsilon);
            OptimalSplitKErrors := OptimalSplitKErrors or AP_FP_Greater(AbsReal(CVE-(-C0*Ln(AP_Double(C0)/(C0+1))-C1*Ln(AP_Double(C1)/(C1+1)))),100*MachineEpsilon);
        end;
        
        //
        // test #2
        //
        if N>2 then
        begin
            C0 := 1+RandomInteger(N-1);
            C1 := N-C0;
            I:=0;
            while I<=N-1 do
            begin
                if I<C0 then
                begin
                    A[I] := 0;
                    C[I] := 0;
                end
                else
                begin
                    A[I] := 1;
                    C[I] := 1;
                end;
                Inc(I);
            end;
            DSOptimalSplitK(A, C, N, 2, 2+RandomInteger(5), Info, Thresholds, NI, CVE);
            if Info<>1 then
            begin
                OptimalSplitKErrors := True;
                Inc(N);
                Continue;
            end;
            OptimalSplitKErrors := OptimalSplitKErrors or (NI<>2);
            OptimalSplitKErrors := OptimalSplitKErrors or AP_FP_Greater(AbsReal(Thresholds[0]-0.5),100*MachineEpsilon);
            OptimalSplitKErrors := OptimalSplitKErrors or AP_FP_Greater(AbsReal(CVE-(-C0*Ln(AP_Double(C0)/(C0+1))-C1*Ln(AP_Double(C1)/(C1+1)))),100*MachineEpsilon);
        end;
        
        //
        // multi-tie test
        //
        if N>=16 then
        begin
            
            //
            // Multi-tie test.
            //
            // First NC-1 ties have C0 entries, remaining NC-th tie
            // have C1 entries.
            //
            NC := Round(Sqrt(N));
            C0 := N div NC;
            C1 := N-C0*(NC-1);
            I:=0;
            while I<=NC-2 do
            begin
                J:=C0*I;
                while J<=C0*(I+1)-1 do
                begin
                    A[J] := J;
                    C[J] := I;
                    Inc(J);
                end;
                Inc(I);
            end;
            J:=C0*(NC-1);
            while J<=N-1 do
            begin
                A[J] := J;
                C[J] := NC-1;
                Inc(J);
            end;
            DSOptimalSplitK(A, C, N, NC, NC+RandomInteger(NC), Info, Thresholds, NI, CVE);
            if Info<>1 then
            begin
                OptimalSplitKErrors := True;
                Inc(N);
                Continue;
            end;
            OptimalSplitKErrors := OptimalSplitKErrors or (NI<>NC);
            if NI=NC then
            begin
                I:=0;
                while I<=NC-2 do
                begin
                    OptimalSplitKErrors := OptimalSplitKErrors or AP_FP_Greater(AbsReal(Thresholds[I]-(C0*(I+1)-1+0.5)),100*MachineEpsilon);
                    Inc(I);
                end;
                CVR := -((NC-1)*C0*Ln(AP_Double(C0)/(C0+NC-1))+C1*Ln(AP_Double(C1)/(C1+NC-1)));
                OptimalSplitKErrors := OptimalSplitKErrors or AP_FP_Greater(AbsReal(CVE-CVR),100*MachineEpsilon);
            end;
        end;
        Inc(N);
    end;
    
    //
    // Non-optimal split-K
    //
    
    //
    // General tests for different N's
    //
    N:=1;
    while N<=MaxNQ do
    begin
        SetLength(A, N-1+1);
        SetLength(C, N-1+1);
        
        //
        // one-tie test
        //
        if N mod 2=0 then
        begin
            I:=0;
            while I<=N-1 do
            begin
                A[I] := Pass;
                C[I] := I mod 2;
                Inc(I);
            end;
            DSSplitK(A, C, N, 2, 2+RandomInteger(5), Info, Thresholds, NI, CVE);
            if Info<>-3 then
            begin
                SplitKErrors := True;
                Inc(N);
                Continue;
            end;
        end;
        
        //
        // two-tie test
        //
        
        //
        // test #1
        //
        if N>1 then
        begin
            C0 := 0;
            C1 := 0;
            I:=0;
            while I<=N-1 do
            begin
                A[I] := I div ((N+1) div 2);
                C[I] := I div ((N+1) div 2);
                if C[I]=0 then
                begin
                    C0 := C0+1;
                end;
                if C[I]=1 then
                begin
                    C1 := C1+1;
                end;
                Inc(I);
            end;
            DSSplitK(A, C, N, 2, 2+RandomInteger(5), Info, Thresholds, NI, CVE);
            if Info<>1 then
            begin
                SplitKErrors := True;
                Inc(N);
                Continue;
            end;
            SplitKErrors := SplitKErrors or (NI<>2);
            if NI=2 then
            begin
                SplitKErrors := SplitKErrors or AP_FP_Greater(AbsReal(Thresholds[0]-0.5),100*MachineEpsilon);
                SplitKErrors := SplitKErrors or AP_FP_Greater(AbsReal(CVE-(-C0*Ln(AP_Double(C0)/(C0+1))-C1*Ln(AP_Double(C1)/(C1+1)))),100*MachineEpsilon);
            end;
        end;
        
        //
        // test #2
        //
        if N>2 then
        begin
            C0 := 1+RandomInteger(N-1);
            C1 := N-C0;
            I:=0;
            while I<=N-1 do
            begin
                if I<C0 then
                begin
                    A[I] := 0;
                    C[I] := 0;
                end
                else
                begin
                    A[I] := 1;
                    C[I] := 1;
                end;
                Inc(I);
            end;
            DSSplitK(A, C, N, 2, 2+RandomInteger(5), Info, Thresholds, NI, CVE);
            if Info<>1 then
            begin
                SplitKErrors := True;
                Inc(N);
                Continue;
            end;
            SplitKErrors := SplitKErrors or (NI<>2);
            if NI=2 then
            begin
                SplitKErrors := SplitKErrors or AP_FP_Greater(AbsReal(Thresholds[0]-0.5),100*MachineEpsilon);
                SplitKErrors := SplitKErrors or AP_FP_Greater(AbsReal(CVE-(-C0*Ln(AP_Double(C0)/(C0+1))-C1*Ln(AP_Double(C1)/(C1+1)))),100*MachineEpsilon);
            end;
        end;
        
        //
        // multi-tie test
        //
        C0:=4;
        while C0<=N do
        begin
            if (N mod C0=0) and (N div C0<=C0) and (N div C0>1) then
            begin
                NC := N div C0;
                I:=0;
                while I<=NC-1 do
                begin
                    J:=C0*I;
                    while J<=C0*(I+1)-1 do
                    begin
                        A[J] := J;
                        C[J] := I;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                DSSplitK(A, C, N, NC, NC+RandomInteger(NC), Info, Thresholds, NI, CVE);
                if Info<>1 then
                begin
                    SplitKErrors := True;
                    Inc(C0);
                    Continue;
                end;
                SplitKErrors := SplitKErrors or (NI<>NC);
                if NI=NC then
                begin
                    I:=0;
                    while I<=NC-2 do
                    begin
                        SplitKErrors := SplitKErrors or AP_FP_Greater(AbsReal(Thresholds[I]-(C0*(I+1)-1+0.5)),100*MachineEpsilon);
                        Inc(I);
                    end;
                    CVR := -NC*C0*Ln(AP_Double(C0)/(C0+NC-1));
                    SplitKErrors := SplitKErrors or AP_FP_Greater(AbsReal(CVE-CVR),100*MachineEpsilon);
                end;
            end;
            Inc(C0);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    WasErrors := TiesErrors or Split2Errors or OptimalSplitKErrors or SplitKErrors;
    if  not Silent then
    begin
        Write(Format('TESTING BASIC DATASET SUBROUTINES'#13#10'',[]));
        Write(Format('TIES:                               ',[]));
        if  not TiesErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('SPLIT-2:                            ',[]));
        if  not Split2Errors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('OPTIMAL SPLIT-K:                    ',[]));
        if  not OptimalSplitKErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('SPLIT-K:                            ',[]));
        if  not SplitKErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
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
    Result :=  not WasErrors;
end;


(*************************************************************************
Unsets 2D array.
*************************************************************************)
procedure Unset2D(var A : TComplex2DArray);
begin
    SetLength(A, 0+1, 0+1);
    A[0,0] := C_Complex(2*RandomReal-1);
end;


(*************************************************************************
Unsets 1D array.
*************************************************************************)
procedure Unset1D(var A : TReal1DArray);
begin
    SetLength(A, 0+1);
    A[0] := 2*RandomReal-1;
end;


(*************************************************************************
Unsets 1D array.
*************************************************************************)
procedure Unset1DI(var A : TInteger1DArray);
begin
    SetLength(A, 0+1);
    A[0] := RandomInteger(3)-1;
end;


procedure TestSortResults(const ASorted : TReal1DArray;
     const P1 : TInteger1DArray;
     const P2 : TInteger1DArray;
     const AOriginal : TReal1DArray;
     N : AlglibInteger;
     var WasErrors : Boolean);
var
    I : AlglibInteger;
    A2 : TReal1DArray;
    T : Double;
    F : TInteger1DArray;
begin
    SetLength(A2, N-1+1);
    SetLength(F, N-1+1);
    
    //
    // is set ordered?
    //
    I:=0;
    while I<=N-2 do
    begin
        WasErrors := WasErrors or AP_FP_Greater(ASorted[I],ASorted[I+1]);
        Inc(I);
    end;
    
    //
    // P1 correctness
    //
    I:=0;
    while I<=N-1 do
    begin
        WasErrors := WasErrors or AP_FP_Neq(ASorted[I],AOriginal[P1[I]]);
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        F[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        F[P1[I]] := F[P1[I]]+1;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        WasErrors := WasErrors or (F[I]<>1);
        Inc(I);
    end;
    
    //
    // P2 correctness
    //
    I:=0;
    while I<=N-1 do
    begin
        A2[I] := AOriginal[I];
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        if P2[I]<>I then
        begin
            T := A2[I];
            A2[I] := A2[P2[I]];
            A2[P2[I]] := T;
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        WasErrors := WasErrors or AP_FP_Neq(ASorted[I],A2[I]);
        Inc(I);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testbdssunit_test_silent():Boolean;
begin
    Result := TestBDSS(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testbdssunit_test():Boolean;
begin
    Result := TestBDSS(False);
end;


end.