unit testtdevdbiunit;
interface
uses Math, Sysutils, Ap, blas, tdbisinv;

function TestTDEVDBI(Silent : Boolean):Boolean;
function testtdevdbiunit_test_silent():Boolean;
function testtdevdbiunit_test():Boolean;

implementation

function RefEVD(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     var Lambda : TReal1DArray;
     var Z : TReal2DArray):Boolean;forward;
function TDBITridiagonalQLIEigenValuesAndVectors(var d : TReal1DArray;
     e : TReal1DArray;
     n : AlglibInteger;
     var z : TReal2DArray):Boolean;forward;
procedure FillDE(var D : TReal1DArray;
     var E : TReal1DArray;
     N : AlglibInteger;
     FillType : AlglibInteger);forward;
function TestProduct(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     const Z : TReal2DArray;
     const Lambda : TReal1DArray):Double;forward;
function TestOrt(const Z : TReal2DArray; N : AlglibInteger):Double;forward;
procedure TestEVDProblem(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     var ValErr : Double;
     var VecErr : Double;
     var WNSorted : Boolean;
     var FailC : AlglibInteger);forward;


(*************************************************************************
Testing EVD, BI
*************************************************************************)
function TestTDEVDBI(Silent : Boolean):Boolean;
var
    D : TReal1DArray;
    E : TReal1DArray;
    Pass : AlglibInteger;
    N : AlglibInteger;
    I : AlglibInteger;
    MKind : AlglibInteger;
    PassCount : AlglibInteger;
    MaxN : AlglibInteger;
    ValErr : Double;
    VecErr : Double;
    WNSorted : Boolean;
    FailC : AlglibInteger;
    Runs : AlglibInteger;
    FailR : Double;
    FailThreshold : Double;
    Threshold : Double;
    WasErrors : Boolean;
    WFailed : Boolean;
    WSpecialF : Boolean;
    M : AlglibInteger;
    Z : TReal2DArray;
begin
    FailThreshold := 0.005;
    Threshold := 1.0E7*MachineEpsilon;
    ValErr := 0;
    VecErr := 0;
    WNSorted := False;
    WFailed := False;
    WSpecialF := False;
    FailC := 0;
    Runs := 0;
    MaxN := 15;
    PassCount := 5;
    
    //
    // Main cycle
    //
    N:=1;
    while N<=MaxN do
    begin
        
        //
        // Different tasks
        //
        MKind:=1;
        while MKind<=4 do
        begin
            FillDE(D, E, N, MKind);
            TestEVDProblem(D, E, N, ValErr, VecErr, WNSorted, FailC);
            Runs := Runs+1;
            Inc(MKind);
        end;
        
        //
        // Special tests
        //
        FillDE(D, E, N, 0);
        if  not SMatrixTDEVDR(D, E, N, 0, -1.0, 0.0, M, Z) then
        begin
            WSpecialF := True;
            Inc(N);
            Continue;
        end;
        WSpecialF := WSpecialF or (M<>N);
        FillDE(D, E, N, 0);
        if  not SMatrixTDEVDR(D, E, N, 0, 0.0, 1.0, M, Z) then
        begin
            WSpecialF := True;
            Inc(N);
            Continue;
        end;
        WSpecialF := WSpecialF or (M<>0);
        I:=0;
        while I<=N-1 do
        begin
            FillDE(D, E, N, 0);
            if  not SMatrixTDEVDI(D, E, N, 0, I, I, Z) then
            begin
                WSpecialF := True;
                Inc(I);
                Continue;
            end;
            WSpecialF := WSpecialF or AP_FP_Neq(D[0],0);
            Inc(I);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    FailR := AP_Double(FailC)/Runs;
    WFailed := AP_FP_Greater(FailR,FailThreshold);
    WasErrors := AP_FP_Greater(ValErr,Threshold) or AP_FP_Greater(VecErr,Threshold) or WNSorted or WFailed or WSpecialF;
    if  not Silent then
    begin
        Write(Format('TESTING TRIDIAGONAL BISECTION AND INVERSE ITERATION EVD'#13#10'',[]));
        Write(Format('EVD values error (different variants):   %5.4e'#13#10'',[
            ValErr]));
        Write(Format('EVD vectors error:                       %5.4e'#13#10'',[
            VecErr]));
        Write(Format('Eigen values order:                      ',[]));
        if  not WNSorted then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('Special tests:                           ',[]));
        if  not WSpecialF then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('Always successfully converged:           ',[]));
        if  not WFailed then
        begin
            Write(Format('YES'#13#10'',[]));
        end
        else
        begin
            Write(Format('NO'#13#10'',[]));
            Write(Format('Fail ratio:                              %5.3f'#13#10'',[
                FailR]));
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


(*************************************************************************
Reference EVD
*************************************************************************)
function RefEVD(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     var Lambda : TReal1DArray;
     var Z : TReal2DArray):Boolean;
var
    Z1 : TReal2DArray;
    D1 : TReal1DArray;
    E1 : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
begin
    SetLength(Z1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        J:=1;
        while J<=N do
        begin
            Z1[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=1;
    while I<=N do
    begin
        Z1[I,I] := 1;
        Inc(I);
    end;
    SetLength(D1, N+1);
    I:=1;
    while I<=N do
    begin
        D1[I] := D[I-1];
        Inc(I);
    end;
    SetLength(E1, N+1);
    I:=2;
    while I<=N do
    begin
        E1[I] := E[I-2];
        Inc(I);
    end;
    Result := TDBITridiagonalQLIEigenValuesAndVectors(D1, E1, N, Z1);
    if Result then
    begin
        
        //
        // copy
        //
        SetLength(Lambda, N-1+1);
        I:=0;
        while I<=N-1 do
        begin
            Lambda[I] := D1[I+1];
            Inc(I);
        end;
        SetLength(Z, N-1+1, N-1+1);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                Z[I,J] := Z1[I+1,J+1];
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Use Selection Sort to minimize swaps of eigenvectors
        //
        I:=0;
        while I<=N-2 do
        begin
            K := I;
            J:=I+1;
            while J<=N-1 do
            begin
                if AP_FP_Less(Lambda[J],Lambda[K]) then
                begin
                    K := J;
                end;
                Inc(J);
            end;
            if K<>I then
            begin
                V := Lambda[I];
                Lambda[I] := Lambda[K];
                Lambda[K] := V;
                J:=0;
                while J<=N-1 do
                begin
                    V := Z[J,I];
                    Z[J,I] := Z[J,K];
                    Z[J,K] := V;
                    Inc(J);
                end;
            end;
            Inc(I);
        end;
    end;
end;


function TDBITridiagonalQLIEigenValuesAndVectors(var d : TReal1DArray;
     e : TReal1DArray;
     n : AlglibInteger;
     var z : TReal2DArray):Boolean;
var
    m : AlglibInteger;
    l : AlglibInteger;
    iter : AlglibInteger;
    i : AlglibInteger;
    k : AlglibInteger;
    s : Double;
    r : Double;
    p : Double;
    g : Double;
    f : Double;
    dd : Double;
    c : Double;
    b : Double;
begin
    e := DynamicArrayCopy(e);
    Result := True;
    if n=1 then
    begin
        Exit;
    end;
    i:=2;
    while i<=n do
    begin
        e[i-1] := e[i];
        Inc(i);
    end;
    e[n] := 0.0;
    l:=1;
    while l<=n do
    begin
        iter := 0;
        repeat
            m:=l;
            while m<=n-1 do
            begin
                dd := AbsReal(d[m])+AbsReal(d[m+1]);
                if AP_FP_Eq(AbsReal(e[m])+dd,dd) then
                begin
                    Break;
                end;
                Inc(m);
            end;
            if m<>l then
            begin
                if iter=30 then
                begin
                    Result := False;
                    Exit;
                end;
                iter := iter+1;
                g := (d[l+1]-d[l])/(2*e[l]);
                if AP_FP_Less(AbsReal(g),1) then
                begin
                    r := Sqrt(1+AP_Sqr(g));
                end
                else
                begin
                    r := AbsReal(g)*Sqrt(1+AP_Sqr(1/g));
                end;
                if AP_FP_Less(g,0) then
                begin
                    g := d[m]-d[l]+e[l]/(g-r);
                end
                else
                begin
                    g := d[m]-d[l]+e[l]/(g+r);
                end;
                s := 1;
                c := 1;
                p := 0;
                i:=m-1;
                while i>=l do
                begin
                    f := s*e[i];
                    b := c*e[i];
                    if AP_FP_Less(AbsReal(f),AbsReal(g)) then
                    begin
                        r := AbsReal(g)*Sqrt(1+AP_Sqr(f/g));
                    end
                    else
                    begin
                        r := AbsReal(f)*Sqrt(1+AP_Sqr(g/f));
                    end;
                    e[i+1] := r;
                    if AP_FP_Eq(r,0) then
                    begin
                        d[i+1] := d[i+1]-p;
                        e[m] := 0;
                        Break;
                    end;
                    s := f/r;
                    c := g/r;
                    g := d[i+1]-p;
                    r := (d[i]-g)*s+2.0*c*b;
                    p := s*r;
                    d[i+1] := g+p;
                    g := c*r-b;
                    k:=1;
                    while k<=n do
                    begin
                        f := z[k,i+1];
                        z[k,i+1] := s*z[k,i]+c*f;
                        z[k,i] := c*z[k,i]-s*f;
                        Inc(k);
                    end;
                    Dec(i);
                end;
                if AP_FP_Eq(r,0) and (i>=1) then
                begin
                    Continue;
                end;
                d[l] := d[l]-p;
                e[l] := g;
                e[m] := 0.0;
            end;
        until m=l;
        Inc(l);
    end;
end;


(*************************************************************************
Fills D and E
*************************************************************************)
procedure FillDE(var D : TReal1DArray;
     var E : TReal1DArray;
     N : AlglibInteger;
     FillType : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    SetLength(D, N-1+1);
    if N>1 then
    begin
        SetLength(E, N-2+1);
    end;
    if FillType=0 then
    begin
        
        //
        // Zero matrix
        //
        I:=0;
        while I<=N-1 do
        begin
            D[I] := 0;
            Inc(I);
        end;
        I:=0;
        while I<=N-2 do
        begin
            E[I] := 0;
            Inc(I);
        end;
        Exit;
    end;
    if FillType=1 then
    begin
        
        //
        // Diagonal matrix
        //
        I:=0;
        while I<=N-1 do
        begin
            D[I] := 2*RandomReal-1;
            Inc(I);
        end;
        I:=0;
        while I<=N-2 do
        begin
            E[I] := 0;
            Inc(I);
        end;
        Exit;
    end;
    if FillType=2 then
    begin
        
        //
        // Off-diagonal matrix
        //
        I:=0;
        while I<=N-1 do
        begin
            D[I] := 0;
            Inc(I);
        end;
        I:=0;
        while I<=N-2 do
        begin
            E[I] := 2*RandomReal-1;
            Inc(I);
        end;
        Exit;
    end;
    if FillType=3 then
    begin
        
        //
        // Dense matrix with blocks
        //
        I:=0;
        while I<=N-1 do
        begin
            D[I] := 2*RandomReal-1;
            Inc(I);
        end;
        I:=0;
        while I<=N-2 do
        begin
            E[I] := 2*RandomReal-1;
            Inc(I);
        end;
        J := 1;
        I := 2;
        while J<=N-2 do
        begin
            E[J] := 0;
            J := J+I;
            I := I+1;
        end;
        Exit;
    end;
    
    //
    // dense matrix
    //
    I:=0;
    while I<=N-1 do
    begin
        D[I] := 2*RandomReal-1;
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        E[I] := 2*RandomReal-1;
        Inc(I);
    end;
end;


(*************************************************************************
Tests Z*Lambda*Z' against tridiag(D,E).
Returns relative error.
*************************************************************************)
function TestProduct(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     const Z : TReal2DArray;
     const Lambda : TReal1DArray):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
    MX : Double;
begin
    Result := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            
            //
            // Calculate V = A[i,j], A = Z*Lambda*Z'
            //
            V := 0;
            K:=0;
            while K<=N-1 do
            begin
                V := V+Z[I,K]*Lambda[K]*Z[J,K];
                Inc(K);
            end;
            
            //
            // Compare
            //
            if AbsInt(I-J)=0 then
            begin
                Result := Max(Result, AbsReal(V-D[I]));
            end;
            if AbsInt(I-J)=1 then
            begin
                Result := Max(Result, AbsReal(V-E[Min(I, J)]));
            end;
            if AbsInt(I-J)>1 then
            begin
                Result := Max(Result, AbsReal(V));
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    MX := 0;
    I:=0;
    while I<=N-1 do
    begin
        MX := Max(MX, AbsReal(D[I]));
        Inc(I);
    end;
    I:=0;
    while I<=N-2 do
    begin
        MX := Max(MX, AbsReal(E[I]));
        Inc(I);
    end;
    if AP_FP_Eq(MX,0) then
    begin
        MX := 1;
    end;
    Result := Result/MX;
end;


(*************************************************************************
Tests Z*Z' against diag(1...1)
Returns absolute error.
*************************************************************************)
function TestOrt(const Z : TReal2DArray; N : AlglibInteger):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
begin
    Result := 0;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := 0.0;
            for i_ := 0 to N-1 do
            begin
                V := V + Z[i_,I]*Z[i_,J];
            end;
            if I=J then
            begin
                V := V-1;
            end;
            Result := Max(Result, AbsReal(V));
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Tests EVD problem
*************************************************************************)
procedure TestEVDProblem(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     var ValErr : Double;
     var VecErr : Double;
     var WNSorted : Boolean;
     var FailC : AlglibInteger);
var
    Lambda : TReal1DArray;
    LambdaRef : TReal1DArray;
    Z : TReal2DArray;
    ZRef : TReal2DArray;
    A1 : TReal2DArray;
    A2 : TReal2DArray;
    AR : TReal2DArray;
    WSucc : Boolean;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    M : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    V : Double;
    A : Double;
    B : Double;
    i_ : AlglibInteger;
begin
    SetLength(LambdaRef, N-1+1);
    SetLength(ZRef, N-1+1, N-1+1);
    SetLength(A1, N-1+1, N-1+1);
    SetLength(A2, N-1+1, N-1+1);
    
    //
    // Reference EVD
    //
    if  not RefEVD(D, E, N, LambdaRef, ZRef) then
    begin
        FailC := FailC+1;
        Exit;
    end;
    
    //
    // Test different combinations
    //
    I1:=0;
    while I1<=N-1 do
    begin
        I2:=I1;
        while I2<=N-1 do
        begin
            
            //
            // Select A, B
            //
            if I1>0 then
            begin
                A := 0.5*(LambdaRef[I1]+LambdaRef[I1-1]);
            end
            else
            begin
                A := LambdaRef[0]-1;
            end;
            if I2<N-1 then
            begin
                B := 0.5*(LambdaRef[I2]+LambdaRef[I2+1]);
            end
            else
            begin
                B := LambdaRef[N-1]+1;
            end;
            
            //
            // Test interval, no vectors
            //
            SetLength(Lambda, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                Lambda[I] := D[I];
                Inc(I);
            end;
            if  not SMatrixTDEVDR(Lambda, E, N, 0, A, B, M, Z) then
            begin
                FailC := FailC+1;
                Exit;
            end;
            if M<>I2-I1+1 then
            begin
                FailC := FailC+1;
                Exit;
            end;
            K:=0;
            while K<=M-1 do
            begin
                ValErr := Max(ValErr, AbsReal(Lambda[K]-LambdaRef[I1+K]));
                Inc(K);
            end;
            
            //
            // Test indexes, no vectors
            //
            SetLength(Lambda, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                Lambda[I] := D[I];
                Inc(I);
            end;
            if  not SMatrixTDEVDI(Lambda, E, N, 0, I1, I2, Z) then
            begin
                FailC := FailC+1;
                Exit;
            end;
            M := I2-I1+1;
            K:=0;
            while K<=M-1 do
            begin
                ValErr := Max(ValErr, AbsReal(Lambda[K]-LambdaRef[I1+K]));
                Inc(K);
            end;
            
            //
            // Test interval, transform vectors
            //
            SetLength(Lambda, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                Lambda[I] := D[I];
                Inc(I);
            end;
            SetLength(A1, N-1+1, N-1+1);
            SetLength(A2, N-1+1, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A1[I,J] := 2*RandomReal-1;
                    A2[I,J] := A1[I,J];
                    Inc(J);
                end;
                Inc(I);
            end;
            if  not SMatrixTDEVDR(Lambda, E, N, 1, A, B, M, A1) then
            begin
                FailC := FailC+1;
                Exit;
            end;
            if M<>I2-I1+1 then
            begin
                FailC := FailC+1;
                Exit;
            end;
            K:=0;
            while K<=M-1 do
            begin
                ValErr := Max(ValErr, AbsReal(Lambda[K]-LambdaRef[I1+K]));
                Inc(K);
            end;
            SetLength(AR, N-1+1, M-1+1);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=M-1 do
                begin
                    V := 0.0;
                    for i_ := 0 to N-1 do
                    begin
                        V := V + A2[I,i_]*ZRef[i_,I1+J];
                    end;
                    AR[I,J] := V;
                    Inc(J);
                end;
                Inc(I);
            end;
            J:=0;
            while J<=M-1 do
            begin
                V := 0.0;
                for i_ := 0 to N-1 do
                begin
                    V := V + A1[i_,J]*AR[i_,J];
                end;
                if AP_FP_Less(V,0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        AR[i_,J] := -1*AR[i_,J];
                    end;
                end;
                Inc(J);
            end;
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=M-1 do
                begin
                    VecErr := Max(VecErr, AbsReal(A1[I,J]-AR[I,J]));
                    Inc(J);
                end;
                Inc(I);
            end;
            
            //
            // Test indexes, transform vectors
            //
            SetLength(Lambda, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                Lambda[I] := D[I];
                Inc(I);
            end;
            SetLength(A1, N-1+1, N-1+1);
            SetLength(A2, N-1+1, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A1[I,J] := 2*RandomReal-1;
                    A2[I,J] := A1[I,J];
                    Inc(J);
                end;
                Inc(I);
            end;
            if  not SMatrixTDEVDI(Lambda, E, N, 1, I1, I2, A1) then
            begin
                FailC := FailC+1;
                Exit;
            end;
            M := I2-I1+1;
            K:=0;
            while K<=M-1 do
            begin
                ValErr := Max(ValErr, AbsReal(Lambda[K]-LambdaRef[I1+K]));
                Inc(K);
            end;
            SetLength(AR, N-1+1, M-1+1);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=M-1 do
                begin
                    V := 0.0;
                    for i_ := 0 to N-1 do
                    begin
                        V := V + A2[I,i_]*ZRef[i_,I1+J];
                    end;
                    AR[I,J] := V;
                    Inc(J);
                end;
                Inc(I);
            end;
            J:=0;
            while J<=M-1 do
            begin
                V := 0.0;
                for i_ := 0 to N-1 do
                begin
                    V := V + A1[i_,J]*AR[i_,J];
                end;
                if AP_FP_Less(V,0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        AR[i_,J] := -1*AR[i_,J];
                    end;
                end;
                Inc(J);
            end;
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=M-1 do
                begin
                    VecErr := Max(VecErr, AbsReal(A1[I,J]-AR[I,J]));
                    Inc(J);
                end;
                Inc(I);
            end;
            
            //
            // Test interval, do not transform vectors
            //
            SetLength(Lambda, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                Lambda[I] := D[I];
                Inc(I);
            end;
            SetLength(Z, 0+1, 0+1);
            if  not SMatrixTDEVDR(Lambda, E, N, 2, A, B, M, Z) then
            begin
                FailC := FailC+1;
                Exit;
            end;
            if M<>I2-I1+1 then
            begin
                FailC := FailC+1;
                Exit;
            end;
            K:=0;
            while K<=M-1 do
            begin
                ValErr := Max(ValErr, AbsReal(Lambda[K]-LambdaRef[I1+K]));
                Inc(K);
            end;
            J:=0;
            while J<=M-1 do
            begin
                V := 0.0;
                for i_ := 0 to N-1 do
                begin
                    V := V + Z[i_,J]*ZRef[i_,I1+J];
                end;
                if AP_FP_Less(V,0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        Z[i_,J] := -1*Z[i_,J];
                    end;
                end;
                Inc(J);
            end;
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=M-1 do
                begin
                    VecErr := Max(VecErr, AbsReal(Z[I,J]-ZRef[I,I1+J]));
                    Inc(J);
                end;
                Inc(I);
            end;
            
            //
            // Test interval, do not transform vectors
            //
            SetLength(Lambda, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                Lambda[I] := D[I];
                Inc(I);
            end;
            SetLength(Z, 0+1, 0+1);
            if  not SMatrixTDEVDI(Lambda, E, N, 2, I1, I2, Z) then
            begin
                FailC := FailC+1;
                Exit;
            end;
            M := I2-I1+1;
            K:=0;
            while K<=M-1 do
            begin
                ValErr := Max(ValErr, AbsReal(Lambda[K]-LambdaRef[I1+K]));
                Inc(K);
            end;
            J:=0;
            while J<=M-1 do
            begin
                V := 0.0;
                for i_ := 0 to N-1 do
                begin
                    V := V + Z[i_,J]*ZRef[i_,I1+J];
                end;
                if AP_FP_Less(V,0) then
                begin
                    for i_ := 0 to N-1 do
                    begin
                        Z[i_,J] := -1*Z[i_,J];
                    end;
                end;
                Inc(J);
            end;
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=M-1 do
                begin
                    VecErr := Max(VecErr, AbsReal(Z[I,J]-ZRef[I,I1+J]));
                    Inc(J);
                end;
                Inc(I);
            end;
            Inc(I2);
        end;
        Inc(I1);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testtdevdbiunit_test_silent():Boolean;
begin
    Result := TestTDEVDBI(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testtdevdbiunit_test():Boolean;
begin
    Result := TestTDEVDBI(False);
end;


end.