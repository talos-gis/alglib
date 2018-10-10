unit testnsevdunit;
interface
uses Math, Sysutils, Ap, blas, reflections, rotations, hsschur, hessenberg, nsevd;

function TestNonSymmetricEVD(Silent : Boolean):Boolean;
function testnsevdunit_test_silent():Boolean;
function testnsevdunit_test():Boolean;

implementation

procedure FillSparseA(var A : TReal2DArray;
     N : AlglibInteger;
     Sparcity : Double);forward;
procedure TestNSEVDProblem(const A : TReal2DArray;
     N : AlglibInteger;
     var VecErr : Double;
     var ValOnlyDiff : Double;
     var WFailed : Boolean);forward;


function TestNonSymmetricEVD(Silent : Boolean):Boolean;
var
    A : TReal2DArray;
    N : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    GPass : AlglibInteger;
    WasErrors : Boolean;
    WFailed : Boolean;
    VecErr : Double;
    ValOnlyDiff : Double;
    Threshold : Double;
begin
    VecErr := 0;
    ValOnlyDiff := 0;
    WFailed := False;
    WasErrors := False;
    Threshold := 1000*MachineEpsilon;
    
    //
    // First set: N = 1..10
    //
    N:=1;
    while N<=10 do
    begin
        SetLength(A, N-1+1, N-1+1);
        
        //
        // zero matrix
        //
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                A[I,J] := 0;
                Inc(J);
            end;
            Inc(I);
        end;
        TestNSEVDProblem(A, N, VecErr, ValOnlyDiff, WFailed);
        
        //
        // Dense and sparse matrices
        //
        GPass:=1;
        while GPass<=1 do
        begin
            
            //
            // Dense matrix
            //
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := 2*RandomReal-1;
                    Inc(J);
                end;
                Inc(I);
            end;
            TestNSEVDProblem(A, N, VecErr, ValOnlyDiff, WFailed);
            
            //
            // Very matrix
            //
            FillSparseA(A, N, 0.98);
            TestNSEVDProblem(A, N, VecErr, ValOnlyDiff, WFailed);
            
            //
            // Incredible sparse matrix
            //
            FillSparseA(A, N, 0.995);
            TestNSEVDProblem(A, N, VecErr, ValOnlyDiff, WFailed);
            Inc(GPass);
        end;
        Inc(N);
    end;
    
    //
    // Second set: N = 70..72
    //
    N:=70;
    while N<=72 do
    begin
        SetLength(A, N-1+1, N-1+1);
        
        //
        // zero matrix
        //
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                A[I,J] := 0;
                Inc(J);
            end;
            Inc(I);
        end;
        TestNSEVDProblem(A, N, VecErr, ValOnlyDiff, WFailed);
        
        //
        // Dense and sparse matrices
        //
        GPass:=1;
        while GPass<=1 do
        begin
            
            //
            // Dense matrix
            //
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := 2*RandomReal-1;
                    Inc(J);
                end;
                Inc(I);
            end;
            TestNSEVDProblem(A, N, VecErr, ValOnlyDiff, WFailed);
            
            //
            // Very matrix
            //
            FillSparseA(A, N, 0.98);
            TestNSEVDProblem(A, N, VecErr, ValOnlyDiff, WFailed);
            
            //
            // Incredible sparse matrix
            //
            FillSparseA(A, N, 0.995);
            TestNSEVDProblem(A, N, VecErr, ValOnlyDiff, WFailed);
            Inc(GPass);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    WasErrors := AP_FP_Greater(ValOnlyDiff,1000*Threshold) or AP_FP_Greater(VecErr,Threshold) or WFailed;
    if  not Silent then
    begin
        Write(Format('TESTING NONSYMMETTRIC EVD'#13#10'',[]));
        Write(Format('Av-lambdav error:                        %5.4e'#13#10'',[
            VecErr]));
        Write(Format('Values only difference:                  %5.4e'#13#10'',[
            ValOnlyDiff]));
        Write(Format('Always converged:                        ',[]));
        if  not WFailed then
        begin
            Write(Format('YES'#13#10'',[]));
        end
        else
        begin
            Write(Format('NO'#13#10'',[]));
        end;
        Write(Format('Threshold:                               %5.4e'#13#10'',[
            Threshold]));
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


procedure FillSparseA(var A : TReal2DArray;
     N : AlglibInteger;
     Sparcity : Double);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if AP_FP_Greater_Eq(RandomReal,Sparcity) then
            begin
                A[I,J] := 2*RandomReal-1;
            end
            else
            begin
                A[I,J] := 0;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


procedure TestNSEVDProblem(const A : TReal2DArray;
     N : AlglibInteger;
     var VecErr : Double;
     var ValOnlyDiff : Double;
     var WFailed : Boolean);
var
    MX : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    VJob : AlglibInteger;
    NeedL : Boolean;
    NeedR : Boolean;
    WR0 : TReal1DArray;
    WI0 : TReal1DArray;
    WR1 : TReal1DArray;
    WI1 : TReal1DArray;
    WR0S : TReal1DArray;
    WI0S : TReal1DArray;
    WR1S : TReal1DArray;
    WI1S : TReal1DArray;
    VL : TReal2DArray;
    VR : TReal2DArray;
    Vec1R : TReal1DArray;
    Vec1I : TReal1DArray;
    Vec2R : TReal1DArray;
    Vec2I : TReal1DArray;
    Vec3R : TReal1DArray;
    Vec3I : TReal1DArray;
    CurWR : Double;
    CurWI : Double;
    VT : Double;
    Tmp : Double;
    i_ : AlglibInteger;
begin
    SetLength(Vec1R, N-1+1);
    SetLength(Vec2R, N-1+1);
    SetLength(Vec3R, N-1+1);
    SetLength(Vec1I, N-1+1);
    SetLength(Vec2I, N-1+1);
    SetLength(Vec3I, N-1+1);
    SetLength(WR0S, N-1+1);
    SetLength(WR1S, N-1+1);
    SetLength(WI0S, N-1+1);
    SetLength(WI1S, N-1+1);
    MX := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if AP_FP_Greater(AbsReal(A[I,J]),MX) then
            begin
                MX := AbsReal(A[I,J]);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(MX,0) then
    begin
        MX := 1;
    end;
    
    //
    // Load values-only
    //
    if  not RMatrixEVD(A, N, 0, WR0, WI0, VL, VR) then
    begin
        WFailed := False;
        Exit;
    end;
    
    //
    // Test different jobs
    //
    VJob:=1;
    while VJob<=3 do
    begin
        NeedR := (VJob=1) or (VJob=3);
        NeedL := (VJob=2) or (VJob=3);
        if  not RMatrixEVD(A, N, VJob, WR1, WI1, VL, VR) then
        begin
            WFailed := False;
            Exit;
        end;
        
        //
        // Test values:
        // 1. sort by real part
        // 2. test
        //
        APVMove(@WR0S[0], 0, N-1, @WR0[0], 0, N-1);
        APVMove(@WI0S[0], 0, N-1, @WI0[0], 0, N-1);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-2-I do
            begin
                if AP_FP_Greater(WR0S[J],WR0S[J+1]) then
                begin
                    Tmp := WR0S[J];
                    WR0S[J] := WR0S[J+1];
                    WR0S[J+1] := Tmp;
                    Tmp := WI0S[J];
                    WI0S[J] := WI0S[J+1];
                    WI0S[J+1] := Tmp;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        APVMove(@WR1S[0], 0, N-1, @WR1[0], 0, N-1);
        APVMove(@WI1S[0], 0, N-1, @WI1[0], 0, N-1);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-2-I do
            begin
                if AP_FP_Greater(WR1S[J],WR1S[J+1]) then
                begin
                    Tmp := WR1S[J];
                    WR1S[J] := WR1S[J+1];
                    WR1S[J+1] := Tmp;
                    Tmp := WI1S[J];
                    WI1S[J] := WI1S[J+1];
                    WI1S[J+1] := Tmp;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            ValOnlyDiff := Max(ValOnlyDiff, AbsReal(WR0S[I]-WR1S[I]));
            ValOnlyDiff := Max(ValOnlyDiff, AbsReal(WI0S[I]-WI1S[I]));
            Inc(I);
        end;
        
        //
        // Test right vectors
        //
        if NeedR then
        begin
            K := 0;
            while K<=N-1 do
            begin
                if AP_FP_Eq(WI1[K],0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        Vec1R[i_] := VR[i_,K];
                    end;
                    I:=0;
                    while I<=N-1 do
                    begin
                        Vec1I[I] := 0;
                        Inc(I);
                    end;
                    CurWR := WR1[K];
                    CurWI := 0;
                end;
                if AP_FP_Greater(WI1[K],0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        Vec1R[i_] := VR[i_,K];
                    end;
                    for i_ := 0 to N-1 do
                    begin
                        Vec1I[i_] := VR[i_,K+1];
                    end;
                    CurWR := WR1[K];
                    CurWI := WI1[K];
                end;
                if AP_FP_Less(WI1[K],0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        Vec1R[i_] := VR[i_,K-1];
                    end;
                    for i_ := 0 to N-1 do
                    begin
                        Vec1I[i_] := -VR[i_,K];
                    end;
                    CurWR := WR1[K];
                    CurWI := WI1[K];
                end;
                I:=0;
                while I<=N-1 do
                begin
                    VT := APVDotProduct(@A[I][0], 0, N-1, @Vec1R[0], 0, N-1);
                    Vec2R[I] := VT;
                    VT := APVDotProduct(@A[I][0], 0, N-1, @Vec1I[0], 0, N-1);
                    Vec2I[I] := VT;
                    Inc(I);
                end;
                APVMove(@Vec3R[0], 0, N-1, @Vec1R[0], 0, N-1, CurWR);
                APVSub(@Vec3R[0], 0, N-1, @Vec1I[0], 0, N-1, CurWI);
                APVMove(@Vec3I[0], 0, N-1, @Vec1R[0], 0, N-1, CurWI);
                APVAdd(@Vec3I[0], 0, N-1, @Vec1I[0], 0, N-1, CurWR);
                I:=0;
                while I<=N-1 do
                begin
                    VecErr := Max(VecErr, AbsReal(Vec2R[I]-Vec3R[I]));
                    VecErr := Max(VecErr, AbsReal(Vec2I[I]-Vec3I[I]));
                    Inc(I);
                end;
                K := K+1;
            end;
        end;
        
        //
        // Test left vectors
        //
        if NeedL then
        begin
            K := 0;
            while K<=N-1 do
            begin
                if AP_FP_Eq(WI1[K],0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        Vec1R[i_] := VL[i_,K];
                    end;
                    I:=0;
                    while I<=N-1 do
                    begin
                        Vec1I[I] := 0;
                        Inc(I);
                    end;
                    CurWR := WR1[K];
                    CurWI := 0;
                end;
                if AP_FP_Greater(WI1[K],0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        Vec1R[i_] := VL[i_,K];
                    end;
                    for i_ := 0 to N-1 do
                    begin
                        Vec1I[i_] := VL[i_,K+1];
                    end;
                    CurWR := WR1[K];
                    CurWI := WI1[K];
                end;
                if AP_FP_Less(WI1[K],0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        Vec1R[i_] := VL[i_,K-1];
                    end;
                    for i_ := 0 to N-1 do
                    begin
                        Vec1I[i_] := -VL[i_,K];
                    end;
                    CurWR := WR1[K];
                    CurWI := WI1[K];
                end;
                J:=0;
                while J<=N-1 do
                begin
                    VT := 0.0;
                    for i_ := 0 to N-1 do
                    begin
                        VT := VT + Vec1R[i_]*A[i_,J];
                    end;
                    Vec2R[J] := VT;
                    VT := 0.0;
                    for i_ := 0 to N-1 do
                    begin
                        VT := VT + Vec1I[i_]*A[i_,J];
                    end;
                    Vec2I[J] := -VT;
                    Inc(J);
                end;
                APVMove(@Vec3R[0], 0, N-1, @Vec1R[0], 0, N-1, CurWR);
                APVAdd(@Vec3R[0], 0, N-1, @Vec1I[0], 0, N-1, CurWI);
                APVMove(@Vec3I[0], 0, N-1, @Vec1R[0], 0, N-1, CurWI);
                APVSub(@Vec3I[0], 0, N-1, @Vec1I[0], 0, N-1, CurWR);
                I:=0;
                while I<=N-1 do
                begin
                    VecErr := Max(VecErr, AbsReal(Vec2R[I]-Vec3R[I]));
                    VecErr := Max(VecErr, AbsReal(Vec2I[I]-Vec3I[I]));
                    Inc(I);
                end;
                K := K+1;
            end;
        end;
        Inc(VJob);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testnsevdunit_test_silent():Boolean;
begin
    Result := TestNonSymmetricEVD(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testnsevdunit_test():Boolean;
begin
    Result := TestNonSymmetricEVD(False);
end;


end.