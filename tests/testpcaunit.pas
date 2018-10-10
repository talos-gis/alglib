unit testpcaunit;
interface
uses Math, Sysutils, Ap, hblas, reflections, creflections, sblas, ablasf, ablas, ortfac, blas, rotations, bdsvd, svd, descriptivestatistics, pca;

function TestPCA(Silent : Boolean):Boolean;
function testpcaunit_test_silent():Boolean;
function testpcaunit_test():Boolean;

implementation

procedure CalculateMV(const X : TReal1DArray;
     N : AlglibInteger;
     var Mean : Double;
     var MeanS : Double;
     var StdDev : Double;
     var StdDevS : Double);forward;


function TestPCA(Silent : Boolean):Boolean;
var
    PassCount : AlglibInteger;
    MaxN : AlglibInteger;
    MaxM : AlglibInteger;
    Threshold : Double;
    M : AlglibInteger;
    N : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    Info : AlglibInteger;
    Means : TReal1DArray;
    S : TReal1DArray;
    T2 : TReal1DArray;
    T3 : TReal1DArray;
    V : TReal2DArray;
    X : TReal2DArray;
    T : Double;
    H : Double;
    TMean : Double;
    TMeanS : Double;
    TStdDev : Double;
    TStdDevS : Double;
    TMean2 : Double;
    TMeanS2 : Double;
    TStdDev2 : Double;
    TStdDevS2 : Double;
    PCAConvErrors : Boolean;
    PCAOrtErrors : Boolean;
    PCAVarErrors : Boolean;
    PCAOptErrors : Boolean;
    WasErrors : Boolean;
    i_ : AlglibInteger;
begin
    
    //
    // Primary settings
    //
    MaxM := 10;
    MaxN := 100;
    PassCount := 1;
    Threshold := 1000*MachineEpsilon;
    WasErrors := False;
    PCAConvErrors := False;
    PCAOrtErrors := False;
    PCAVarErrors := False;
    PCAOptErrors := False;
    
    //
    // Test 1: N random points in M-dimensional space
    //
    M:=1;
    while M<=MaxM do
    begin
        N:=1;
        while N<=MaxN do
        begin
            
            //
            // Generate task
            //
            SetLength(X, N-1+1, M-1+1);
            SetLength(Means, M-1+1);
            J:=0;
            while J<=M-1 do
            begin
                Means[J] := 1.5*RandomReal-0.75;
                Inc(J);
            end;
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=M-1 do
                begin
                    X[I,J] := Means[J]+(2*RandomReal-1);
                    Inc(J);
                end;
                Inc(I);
            end;
            
            //
            // Solve
            //
            PCABuildBasis(X, N, M, Info, S, V);
            if Info<>1 then
            begin
                PCAConvErrors := True;
                Inc(N);
                Continue;
            end;
            
            //
            // Orthogonality test
            //
            I:=0;
            while I<=M-1 do
            begin
                J:=0;
                while J<=M-1 do
                begin
                    T := 0.0;
                    for i_ := 0 to M-1 do
                    begin
                        T := T + V[i_,I]*V[i_,J];
                    end;
                    if I=J then
                    begin
                        T := T-1;
                    end;
                    PCAOrtErrors := PCAOrtErrors or AP_FP_Greater(AbsReal(T),Threshold);
                    Inc(J);
                end;
                Inc(I);
            end;
            
            //
            // Variance test
            //
            SetLength(T2, N-1+1);
            K:=0;
            while K<=M-1 do
            begin
                I:=0;
                while I<=N-1 do
                begin
                    T := 0.0;
                    for i_ := 0 to M-1 do
                    begin
                        T := T + X[I,i_]*V[i_,K];
                    end;
                    T2[I] := T;
                    Inc(I);
                end;
                CalculateMV(T2, N, TMean, TMeanS, TStdDev, TStdDevS);
                if N<>1 then
                begin
                    T := AP_Sqr(TStdDev)*N/(N-1);
                end
                else
                begin
                    T := 0;
                end;
                PCAVarErrors := PCAVarErrors or AP_FP_Greater(AbsReal(T-S[K]),Threshold);
                Inc(K);
            end;
            K:=0;
            while K<=M-2 do
            begin
                PCAVarErrors := PCAVarErrors or AP_FP_Less(S[K],S[K+1]);
                Inc(K);
            end;
            
            //
            // Optimality: different perturbations in V[..,0] can't
            // increase variance of projection - can only decrease.
            //
            SetLength(T2, N-1+1);
            SetLength(T3, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                T := 0.0;
                for i_ := 0 to M-1 do
                begin
                    T := T + X[I,i_]*V[i_,0];
                end;
                T2[I] := T;
                Inc(I);
            end;
            CalculateMV(T2, N, TMean, TMeanS, TStdDev, TStdDevS);
            K:=0;
            while K<=2*M-1 do
            begin
                H := 0.001;
                if K mod 2<>0 then
                begin
                    H := -H;
                end;
                APVMove(@T3[0], 0, N-1, @T2[0], 0, N-1);
                for i_ := 0 to N-1 do
                begin
                    T3[i_] := T3[i_] + H*X[i_,K div 2];
                end;
                T := 0;
                J:=0;
                while J<=M-1 do
                begin
                    if J<>K div 2 then
                    begin
                        T := T+AP_Sqr(V[J,0]);
                    end
                    else
                    begin
                        T := T+AP_Sqr(V[J,0]+H);
                    end;
                    Inc(J);
                end;
                T := 1/Sqrt(T);
                APVMul(@T3[0], 0, N-1, T);
                CalculateMV(T3, N, TMean2, TMeanS2, TStdDev2, TStdDevS2);
                PCAOptErrors := PCAOptErrors or AP_FP_Greater(TStdDev2,TStdDev+Threshold);
                Inc(K);
            end;
            Inc(N);
        end;
        Inc(M);
    end;
    
    //
    // Special test for N=0
    //
    M:=1;
    while M<=MaxM do
    begin
        
        //
        // Solve
        //
        PCABuildBasis(X, 0, M, Info, S, V);
        if Info<>1 then
        begin
            PCAConvErrors := True;
            Inc(M);
            Continue;
        end;
        
        //
        // Orthogonality test
        //
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                T := 0.0;
                for i_ := 0 to M-1 do
                begin
                    T := T + V[i_,I]*V[i_,J];
                end;
                if I=J then
                begin
                    T := T-1;
                end;
                PCAOrtErrors := PCAOrtErrors or AP_FP_Greater(AbsReal(T),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
        Inc(M);
    end;
    
    //
    // Final report
    //
    WasErrors := PCAConvErrors or PCAOrtErrors or PCAVarErrors or PCAOptErrors;
    if  not Silent then
    begin
        Write(Format('PCA TEST'#13#10'',[]));
        Write(Format('TOTAL RESULTS:                           ',[]));
        if  not WasErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('* CONVERGENCE                            ',[]));
        if  not PCAConvErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('* ORTOGONALITY                           ',[]));
        if  not PCAOrtErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('* VARIANCE REPORT                        ',[]));
        if  not PCAVarErrors then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('* OPTIMALITY                             ',[]));
        if  not PCAOptErrors then
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
Moments estimates and their errors
*************************************************************************)
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
        Variance := (V1-V2)/N;
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


(*************************************************************************
Silent unit test
*************************************************************************)
function testpcaunit_test_silent():Boolean;
begin
    Result := TestPCA(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testpcaunit_test():Boolean;
begin
    Result := TestPCA(False);
end;


end.