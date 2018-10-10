unit testnearestneighborunit;
interface
uses Math, Sysutils, Ap, tsort, nearestneighbor;

function TestNearestNeighbor(Silent : Boolean):Boolean;
function testnearestneighborunit_test_silent():Boolean;
function testnearestneighborunit_test():Boolean;

implementation

procedure Unset2D(var A : TComplex2DArray);forward;
procedure Unset1D(var A : TReal1DArray);forward;
function KDTResultsDifferent(const RefXY : TReal2DArray;
     NTotal : AlglibInteger;
     const QX : TReal2DArray;
     const QXY : TReal2DArray;
     const QT : TInteger1DArray;
     N : AlglibInteger;
     NX : AlglibInteger;
     NY : AlglibInteger):Boolean;forward;
function VNorm(const X : TReal1DArray;
     N : AlglibInteger;
     NormType : AlglibInteger):Double;forward;
procedure TestKDTUniform(const XY : TReal2DArray;
     const N : AlglibInteger;
     const NX : AlglibInteger;
     const NY : AlglibInteger;
     const NormType : AlglibInteger;
     var KDTErrors : Boolean);forward;


(*************************************************************************
Testing Nearest Neighbor Search
*************************************************************************)
function TestNearestNeighbor(Silent : Boolean):Boolean;
var
    XY : TReal2DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    NormType : AlglibInteger;
    NX : AlglibInteger;
    NY : AlglibInteger;
    N : AlglibInteger;
    SmallN : AlglibInteger;
    LargeN : AlglibInteger;
    PassCount : AlglibInteger;
    Pass : AlglibInteger;
    WasErrors : Boolean;
    KDTErrors : Boolean;
begin
    KDTErrors := False;
    PassCount := 2;
    SmallN := 256;
    LargeN := 2048;
    NY := 3;
    
    //
    //
    //
    Pass:=1;
    while Pass<=PassCount do
    begin
        NormType:=0;
        while NormType<=2 do
        begin
            NX:=1;
            while NX<=3 do
            begin
                
                //
                // Test in hypercube
                //
                SetLength(XY, LargeN, NX+NY);
                I:=0;
                while I<=LargeN-1 do
                begin
                    J:=0;
                    while J<=NX+NY-1 do
                    begin
                        XY[I,J] := 10*RandomReal-5;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                N:=1;
                while N<=10 do
                begin
                    TestKDTUniform(XY, N, NX, RandomInteger(NY+1), NormType, KDTErrors);
                    Inc(N);
                end;
                TestKDTUniform(XY, LargeN, NX, RandomInteger(NY+1), NormType, KDTErrors);
                
                //
                // Test clustered (2*N points, pairs of equal points)
                //
                SetLength(XY, 2*SmallN, NX+NY);
                I:=0;
                while I<=SmallN-1 do
                begin
                    J:=0;
                    while J<=NX+NY-1 do
                    begin
                        XY[2*I+0,J] := 10*RandomReal-5;
                        XY[2*I+1,J] := XY[2*I+0,J];
                        Inc(J);
                    end;
                    Inc(I);
                end;
                TestKDTUniform(XY, 2*SmallN, NX, RandomInteger(NY+1), NormType, KDTErrors);
                
                //
                // Test degenerate case: all points are same except for one
                //
                SetLength(XY, SmallN, NX+NY);
                V := RandomReal;
                I:=0;
                while I<=SmallN-2 do
                begin
                    J:=0;
                    while J<=NX+NY-1 do
                    begin
                        XY[I,J] := V;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                J:=0;
                while J<=NX+NY-1 do
                begin
                    XY[SmallN-1,J] := 10*RandomReal-5;
                    Inc(J);
                end;
                TestKDTUniform(XY, SmallN, NX, RandomInteger(NY+1), NormType, KDTErrors);
                Inc(NX);
            end;
            Inc(NormType);
        end;
        Inc(Pass);
    end;
    
    //
    // report
    //
    WasErrors := KDTErrors;
    if  not Silent then
    begin
        Write(Format('TESTING NEAREST NEIGHBOR SEARCH'#13#10'',[]));
        Write(Format('* KD TREES:                              ',[]));
        if  not KDTErrors then
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
Compare results from different queries:
* X     just X-values
* XY    X-values and Y-values
* XT    X-values and tag values
*************************************************************************)
function KDTResultsDifferent(const RefXY : TReal2DArray;
     NTotal : AlglibInteger;
     const QX : TReal2DArray;
     const QXY : TReal2DArray;
     const QT : TInteger1DArray;
     N : AlglibInteger;
     NX : AlglibInteger;
     NY : AlglibInteger):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Result := False;
    I:=0;
    while I<=N-1 do
    begin
        if (QT[I]<0) or (QT[I]>=NTotal) then
        begin
            Result := True;
            Exit;
        end;
        J:=0;
        while J<=NX-1 do
        begin
            Result := Result or AP_FP_Neq(QX[I,J],RefXY[QT[I],J]);
            Result := Result or AP_FP_Neq(QXY[I,J],RefXY[QT[I],J]);
            Inc(J);
        end;
        J:=0;
        while J<=NY-1 do
        begin
            Result := Result or AP_FP_Neq(QXY[I,NX+J],RefXY[QT[I],NX+J]);
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Returns norm
*************************************************************************)
function VNorm(const X : TReal1DArray;
     N : AlglibInteger;
     NormType : AlglibInteger):Double;
var
    I : AlglibInteger;
begin
    Result := RandomReal;
    if NormType=0 then
    begin
        Result := 0;
        I:=0;
        while I<=N-1 do
        begin
            Result := Max(Result, AbsReal(X[I]));
            Inc(I);
        end;
        Exit;
    end;
    if NormType=1 then
    begin
        Result := 0;
        I:=0;
        while I<=N-1 do
        begin
            Result := Result+AbsReal(X[I]);
            Inc(I);
        end;
        Exit;
    end;
    if NormType=2 then
    begin
        Result := 0;
        I:=0;
        while I<=N-1 do
        begin
            Result := Result+AP_Sqr(X[I]);
            Inc(I);
        end;
        Result := Sqrt(Result);
        Exit;
    end;
end;


(*************************************************************************
Testing Nearest Neighbor Search on uniformly distributed hypercube

NormType: 0, 1, 2
D: space dimension
N: points count
*************************************************************************)
procedure TestKDTUniform(const XY : TReal2DArray;
     const N : AlglibInteger;
     const NX : AlglibInteger;
     const NY : AlglibInteger;
     const NormType : AlglibInteger;
     var KDTErrors : Boolean);
var
    ErrTol : Double;
    Tags : TInteger1DArray;
    PtX : TReal1DArray;
    TmpX : TReal1DArray;
    TmpB : TBoolean1DArray;
    TreeX : KDTree;
    TreeXY : KDTree;
    TreeXT : KDTree;
    QX : TReal2DArray;
    QXY : TReal2DArray;
    QTags : TInteger1DArray;
    QR : TReal1DArray;
    KX : AlglibInteger;
    KXY : AlglibInteger;
    KT : AlglibInteger;
    KR : AlglibInteger;
    Eps : Double;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    Task : AlglibInteger;
    IsEqual : Boolean;
    R : Double;
    Q : AlglibInteger;
    QCount : AlglibInteger;
begin
    QCount := 10;
    
    //
    // Tol - roundoff error tolerance (for '>=' comparisons)
    //
    ErrTol := 100000*MachineEpsilon;
    
    //
    // fill tags
    //
    SetLength(Tags, N);
    I:=0;
    while I<=N-1 do
    begin
        Tags[I] := I;
        Inc(I);
    end;
    
    //
    // build trees
    //
    KDTreeBuild(XY, N, NX, 0, NormType, TreeX);
    KDTreeBuild(XY, N, NX, NY, NormType, TreeXY);
    KDTreeBuildTagged(XY, Tags, N, NX, 0, NormType, TreeXT);
    
    //
    // allocate arrays
    //
    SetLength(TmpX, NX);
    SetLength(TmpB, N);
    SetLength(QX, N, NX);
    SetLength(QXY, N, NX+NY);
    SetLength(QTags, N);
    SetLength(QR, N);
    SetLength(PtX, NX);
    
    //
    // test general K-NN queries (with self-matches):
    // * compare results from different trees (must be equal) and
    //   check that correct (value,tag) pairs are returned
    // * test results from XT tree - let R be radius of query result.
    //   then all points not in result must be not closer than R.
    //
    Q:=1;
    while Q<=QCount do
    begin
        
        //
        // Select K: 1..N
        //
        if AP_FP_Greater(RandomReal,0.5) then
        begin
            K := 1+RandomInteger(N);
        end
        else
        begin
            K := 1;
        end;
        
        //
        // Select point (either one of the points, or random)
        //
        if AP_FP_Greater(RandomReal,0.5) then
        begin
            I := RandomInteger(N);
            APVMove(@PtX[0], 0, NX-1, @XY[I][0], 0, NX-1);
        end
        else
        begin
            I:=0;
            while I<=NX-1 do
            begin
                PtX[I] := 2*RandomReal-1;
                Inc(I);
            end;
        end;
        
        //
        // Test:
        // * consistency of results from different queries
        // * points in query are IN the R-sphere (or at the boundary),
        //   and points not in query are outside of the R-sphere (or at the boundary)
        // * distances are correct and are ordered
        //
        KX := KDTreeQueryKNN(TreeX, PtX, K, True);
        KXY := KDTreeQueryKNN(TreeXY, PtX, K, True);
        KT := KDTreeQueryKNN(TreeXT, PtX, K, True);
        if (KX<>K) or (KXY<>K) or (KT<>K) then
        begin
            KDTErrors := True;
            Exit;
        end;
        KX := 0;
        KXY := 0;
        KT := 0;
        KDTreeQueryResultsX(TreeX, QX, KX);
        KDTreeQueryResultsXY(TreeXY, QXY, KXY);
        KDTreeQueryResultsTags(TreeXT, QTags, KT);
        KDTreeQueryResultsDistances(TreeXT, QR, KR);
        if (KX<>K) or (KXY<>K) or (KT<>K) or (KR<>K) then
        begin
            KDTErrors := True;
            Exit;
        end;
        KDTErrors := KDTErrors or KDTResultsDifferent(XY, N, QX, QXY, QTags, K, NX, NY);
        I:=0;
        while I<=N-1 do
        begin
            TmpB[I] := True;
            Inc(I);
        end;
        R := 0;
        I:=0;
        while I<=K-1 do
        begin
            TmpB[QTags[I]] := False;
            APVMove(@TmpX[0], 0, NX-1, @PtX[0], 0, NX-1);
            APVSub(@TmpX[0], 0, NX-1, @QX[I][0], 0, NX-1);
            R := Max(R, VNorm(TmpX, NX, NormType));
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            if TmpB[I] then
            begin
                APVMove(@TmpX[0], 0, NX-1, @PtX[0], 0, NX-1);
                APVSub(@TmpX[0], 0, NX-1, @XY[I][0], 0, NX-1);
                KDTErrors := KDTErrors or AP_FP_Less(VNorm(TmpX, NX, NormType),R*(1-ErrTol));
            end;
            Inc(I);
        end;
        I:=0;
        while I<=K-2 do
        begin
            KDTErrors := KDTErrors or AP_FP_Greater(QR[I],QR[I+1]);
            Inc(I);
        end;
        I:=0;
        while I<=K-1 do
        begin
            APVMove(@TmpX[0], 0, NX-1, @PtX[0], 0, NX-1);
            APVSub(@TmpX[0], 0, NX-1, @XY[QTags[I]][0], 0, NX-1);
            KDTErrors := KDTErrors or AP_FP_Greater(AbsReal(VNorm(TmpX, NX, NormType)-QR[I]),ErrTol);
            Inc(I);
        end;
        Inc(Q);
    end;
    
    //
    // test general approximate K-NN queries (with self-matches):
    // * compare results from different trees (must be equal) and
    //   check that correct (value,tag) pairs are returned
    // * test results from XT tree - let R be radius of query result.
    //   then all points not in result must be not closer than R/(1+Eps).
    //
    Q:=1;
    while Q<=QCount do
    begin
        
        //
        // Select K: 1..N
        //
        if AP_FP_Greater(RandomReal,0.5) then
        begin
            K := 1+RandomInteger(N);
        end
        else
        begin
            K := 1;
        end;
        
        //
        // Select Eps
        //
        Eps := 0.5+RandomReal;
        
        //
        // Select point (either one of the points, or random)
        //
        if AP_FP_Greater(RandomReal,0.5) then
        begin
            I := RandomInteger(N);
            APVMove(@PtX[0], 0, NX-1, @XY[I][0], 0, NX-1);
        end
        else
        begin
            I:=0;
            while I<=NX-1 do
            begin
                PtX[I] := 2*RandomReal-1;
                Inc(I);
            end;
        end;
        
        //
        // Test:
        // * consistency of results from different queries
        // * points in query are IN the R-sphere (or at the boundary),
        //   and points not in query are outside of the R-sphere (or at the boundary)
        // * distances are correct and are ordered
        //
        KX := KDTreeQueryAKNN(TreeX, PtX, K, True, Eps);
        KXY := KDTreeQueryAKNN(TreeXY, PtX, K, True, Eps);
        KT := KDTreeQueryAKNN(TreeXT, PtX, K, True, Eps);
        if (KX<>K) or (KXY<>K) or (KT<>K) then
        begin
            KDTErrors := True;
            Exit;
        end;
        KX := 0;
        KXY := 0;
        KT := 0;
        KDTreeQueryResultsX(TreeX, QX, KX);
        KDTreeQueryResultsXY(TreeXY, QXY, KXY);
        KDTreeQueryResultsTags(TreeXT, QTags, KT);
        KDTreeQueryResultsDistances(TreeXT, QR, KR);
        if (KX<>K) or (KXY<>K) or (KT<>K) or (KR<>K) then
        begin
            KDTErrors := True;
            Exit;
        end;
        KDTErrors := KDTErrors or KDTResultsDifferent(XY, N, QX, QXY, QTags, K, NX, NY);
        I:=0;
        while I<=N-1 do
        begin
            TmpB[I] := True;
            Inc(I);
        end;
        R := 0;
        I:=0;
        while I<=K-1 do
        begin
            TmpB[QTags[I]] := False;
            APVMove(@TmpX[0], 0, NX-1, @PtX[0], 0, NX-1);
            APVSub(@TmpX[0], 0, NX-1, @QX[I][0], 0, NX-1);
            R := Max(R, VNorm(TmpX, NX, NormType));
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            if TmpB[I] then
            begin
                APVMove(@TmpX[0], 0, NX-1, @PtX[0], 0, NX-1);
                APVSub(@TmpX[0], 0, NX-1, @XY[I][0], 0, NX-1);
                KDTErrors := KDTErrors or AP_FP_Less(VNorm(TmpX, NX, NormType),R*(1-ErrTol)/(1+Eps));
            end;
            Inc(I);
        end;
        I:=0;
        while I<=K-2 do
        begin
            KDTErrors := KDTErrors or AP_FP_Greater(QR[I],QR[I+1]);
            Inc(I);
        end;
        I:=0;
        while I<=K-1 do
        begin
            APVMove(@TmpX[0], 0, NX-1, @PtX[0], 0, NX-1);
            APVSub(@TmpX[0], 0, NX-1, @XY[QTags[I]][0], 0, NX-1);
            KDTErrors := KDTErrors or AP_FP_Greater(AbsReal(VNorm(TmpX, NX, NormType)-QR[I]),ErrTol);
            Inc(I);
        end;
        Inc(Q);
    end;
    
    //
    // test general R-NN queries  (with self-matches):
    // * compare results from different trees (must be equal) and
    //   check that correct (value,tag) pairs are returned
    // * test results from XT tree - let R be radius of query result.
    //   then all points not in result must be not closer than R.
    //
    Q:=1;
    while Q<=QCount do
    begin
        
        //
        // Select R
        //
        if AP_FP_Greater(RandomReal,0.3) then
        begin
            R := Max(RandomReal, MachineEpsilon);
        end
        else
        begin
            R := MachineEpsilon;
        end;
        
        //
        // Select point (either one of the points, or random)
        //
        if AP_FP_Greater(RandomReal,0.5) then
        begin
            I := RandomInteger(N);
            APVMove(@PtX[0], 0, NX-1, @XY[I][0], 0, NX-1);
        end
        else
        begin
            I:=0;
            while I<=NX-1 do
            begin
                PtX[I] := 2*RandomReal-1;
                Inc(I);
            end;
        end;
        
        //
        // Test:
        // * consistency of results from different queries
        // * points in query are IN the R-sphere (or at the boundary),
        //   and points not in query are outside of the R-sphere (or at the boundary)
        // * distances are correct and are ordered
        //
        KX := KDTreeQueryRNN(TreeX, PtX, R, True);
        KXY := KDTreeQueryRNN(TreeXY, PtX, R, True);
        KT := KDTreeQueryRNN(TreeXT, PtX, R, True);
        if (KXY<>KX) or (KT<>KX) then
        begin
            KDTErrors := True;
            Exit;
        end;
        KX := 0;
        KXY := 0;
        KT := 0;
        KDTreeQueryResultsX(TreeX, QX, KX);
        KDTreeQueryResultsXY(TreeXY, QXY, KXY);
        KDTreeQueryResultsTags(TreeXT, QTags, KT);
        KDTreeQueryResultsDistances(TreeXT, QR, KR);
        if (KXY<>KX) or (KT<>KX) or (KR<>KX) then
        begin
            KDTErrors := True;
            Exit;
        end;
        KDTErrors := KDTErrors or KDTResultsDifferent(XY, N, QX, QXY, QTags, KX, NX, NY);
        I:=0;
        while I<=N-1 do
        begin
            TmpB[I] := True;
            Inc(I);
        end;
        I:=0;
        while I<=KX-1 do
        begin
            TmpB[QTags[I]] := False;
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            APVMove(@TmpX[0], 0, NX-1, @PtX[0], 0, NX-1);
            APVSub(@TmpX[0], 0, NX-1, @XY[I][0], 0, NX-1);
            if TmpB[I] then
            begin
                KDTErrors := KDTErrors or AP_FP_Less(VNorm(TmpX, NX, NormType),R*(1-ErrTol));
            end
            else
            begin
                KDTErrors := KDTErrors or AP_FP_Greater(VNorm(TmpX, NX, NormType),R*(1+ErrTol));
            end;
            Inc(I);
        end;
        I:=0;
        while I<=KX-2 do
        begin
            KDTErrors := KDTErrors or AP_FP_Greater(QR[I],QR[I+1]);
            Inc(I);
        end;
        Inc(Q);
    end;
    
    //
    // Test self-matching:
    // * self-match - nearest neighbor of each point in XY is the point itself
    // * no self-match - nearest neighbor is NOT the point itself
    //
    if N>1 then
    begin
        
        //
        // test for N=1 have non-general form, but it is not really needed
        //
        Task:=0;
        while Task<=1 do
        begin
            I:=0;
            while I<=N-1 do
            begin
                APVMove(@PtX[0], 0, NX-1, @XY[I][0], 0, NX-1);
                KX := KDTreeQueryKNN(TreeX, PtX, 1, Task=0);
                KDTreeQueryResultsX(TreeX, QX, KX);
                if KX<>1 then
                begin
                    KDTErrors := True;
                    Exit;
                end;
                IsEqual := True;
                J:=0;
                while J<=NX-1 do
                begin
                    IsEqual := IsEqual and AP_FP_Eq(QX[0,J],PtX[J]);
                    Inc(J);
                end;
                if Task=0 then
                begin
                    KDTErrors := KDTErrors or  not IsEqual;
                end
                else
                begin
                    KDTErrors := KDTErrors or IsEqual;
                end;
                Inc(I);
            end;
            Inc(Task);
        end;
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testnearestneighborunit_test_silent():Boolean;
begin
    Result := TestNearestNeighbor(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testnearestneighborunit_test():Boolean;
begin
    Result := TestNearestNeighbor(False);
end;


end.