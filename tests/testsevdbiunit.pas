unit testsevdbiunit;
interface
uses Math, Sysutils, Ap, sblas, reflections, tridiagonal, blas, tdbisinv, sbisinv;

function TestSEVDBI(Silent : Boolean):Boolean;
function testsevdbiunit_test_silent():Boolean;
function testsevdbiunit_test():Boolean;

implementation

procedure Unset2D(var A : TReal2DArray);forward;
procedure Unset1D(var A : TReal1DArray);forward;
function RefEVD(A : TReal2DArray;
     N : AlglibInteger;
     var Lambda : TReal1DArray;
     var Z : TReal2DArray):Boolean;forward;
function SBITridiagonalQLIEigenValuesAndVectors(var d : TReal1DArray;
     e : TReal1DArray;
     n : AlglibInteger;
     var z : TReal2DArray):Boolean;forward;
function TestProduct(const D : TReal1DArray;
     const E : TReal1DArray;
     N : AlglibInteger;
     const Z : TReal2DArray;
     const Lambda : TReal1DArray):Double;forward;
function TestOrt(const Z : TReal2DArray; N : AlglibInteger):Double;forward;
procedure TestEVDProblem(const AFull : TReal2DArray;
     const AL : TReal2DArray;
     const AU : TReal2DArray;
     N : AlglibInteger;
     var ValErr : Double;
     var VecErr : Double;
     var WNSorted : Boolean;
     var FailC : AlglibInteger);forward;


(*************************************************************************
Testing symmetric EVD, BI
*************************************************************************)
function TestSEVDBI(Silent : Boolean):Boolean;
var
    A : TReal2DArray;
    AL : TReal2DArray;
    AU : TReal2DArray;
    Pass : AlglibInteger;
    N : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
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
    M : AlglibInteger;
    Z : TReal2DArray;
begin
    FailThreshold := 0.005;
    Threshold := 1000*MachineEpsilon;
    ValErr := 0;
    VecErr := 0;
    WNSorted := False;
    WFailed := False;
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
        // Prepare
        //
        SetLength(A, N-1+1, N-1+1);
        SetLength(AL, N-1+1, N-1+1);
        SetLength(AU, N-1+1, N-1+1);
        I:=0;
        while I<=N-1 do
        begin
            J:=I+1;
            while J<=N-1 do
            begin
                
                //
                // A
                //
                A[I,J] := 2*RandomReal-1;
                A[J,I] := A[I,J];
                
                //
                // A lower
                //
                AL[I,J] := 2*RandomReal-1;
                AL[J,I] := A[I,J];
                
                //
                // A upper
                //
                AU[I,J] := A[I,J];
                AU[J,I] := 2*RandomReal-1;
                Inc(J);
            end;
            A[I,I] := 2*RandomReal-1;
            AL[I,I] := A[I,I];
            AU[I,I] := A[I,I];
            Inc(I);
        end;
        TestEVDProblem(A, AL, AU, N, ValErr, VecErr, WNSorted, FailC);
        Runs := Runs+1;
        Inc(N);
    end;
    
    //
    // report
    //
    FailR := AP_Double(FailC)/Runs;
    WFailed := AP_FP_Greater(FailR,FailThreshold);
    WasErrors := AP_FP_Greater(ValErr,Threshold) or AP_FP_Greater(VecErr,Threshold) or WNSorted or WFailed;
    if  not Silent then
    begin
        Write(Format('TESTING SYMMETRIC BISECTION AND INVERSE ITERATION EVD'#13#10'',[]));
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
Unsets 2D array.
*************************************************************************)
procedure Unset2D(var A : TReal2DArray);
begin
    SetLength(A, 0+1, 0+1);
    A[0,0] := 2*RandomReal-1;
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
Reference EVD
*************************************************************************)
function RefEVD(A : TReal2DArray;
     N : AlglibInteger;
     var Lambda : TReal1DArray;
     var Z : TReal2DArray):Boolean;
var
    Z2 : TReal2DArray;
    Z1 : TReal2DArray;
    D : TReal1DArray;
    E : TReal1DArray;
    Tau : TReal1DArray;
    D1 : TReal1DArray;
    E1 : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    V : Double;
begin
    A := DynamicArrayCopy(A);
    
    //
    // to tridiagonal
    //
    SMatrixTD(A, N, True, Tau, D, E);
    SMatrixTDUnpackQ(A, N, True, Tau, Z2);
    
    //
    // TDEVD
    //
    SetLength(Z1, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        J:=1;
        while J<=N do
        begin
            Z1[I,J] := Z2[I-1,J-1];
            Inc(J);
        end;
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
    Result := SBITridiagonalQLIEigenValuesAndVectors(D1, E1, N, Z1);
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


function SBITridiagonalQLIEigenValuesAndVectors(var d : TReal1DArray;
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
procedure TestEVDProblem(const AFull : TReal2DArray;
     const AL : TReal2DArray;
     const AU : TReal2DArray;
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
    if  not RefEVD(AFull, N, LambdaRef, ZRef) then
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
            // Test interval, no vectors, lower A
            //
            Unset1D(Lambda);
            Unset2D(Z);
            if  not SMatrixEVDR(AL, N, 0, False, A, B, M, Lambda, Z) then
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
            // Test interval, no vectors, upper A
            //
            Unset1D(Lambda);
            Unset2D(Z);
            if  not SMatrixEVDR(AU, N, 0, True, A, B, M, Lambda, Z) then
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
            // Test indexes, no vectors, lower A
            //
            Unset1D(Lambda);
            Unset2D(Z);
            if  not SMatrixEVDI(AL, N, 0, False, I1, I2, Lambda, Z) then
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
            // Test indexes, no vectors, upper A
            //
            Unset1D(Lambda);
            Unset2D(Z);
            if  not SMatrixEVDI(AU, N, 0, True, I1, I2, Lambda, Z) then
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
            // Test interval, do not transform vectors, lower A
            //
            Unset1D(Lambda);
            Unset2D(Z);
            if  not SMatrixEVDR(AL, N, 1, False, A, B, M, Lambda, Z) then
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
            // Test interval, do not transform vectors, upper A
            //
            Unset1D(Lambda);
            Unset2D(Z);
            if  not SMatrixEVDR(AU, N, 1, True, A, B, M, Lambda, Z) then
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
            // Test indexes, do not transform vectors, lower A
            //
            Unset1D(Lambda);
            Unset2D(Z);
            if  not SMatrixEVDI(AL, N, 1, False, I1, I2, Lambda, Z) then
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
            
            //
            // Test indexes, do not transform vectors, upper A
            //
            Unset1D(Lambda);
            Unset2D(Z);
            if  not SMatrixEVDI(AU, N, 1, True, I1, I2, Lambda, Z) then
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
function testsevdbiunit_test_silent():Boolean;
begin
    Result := TestSEVDBI(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testsevdbiunit_test():Boolean;
begin
    Result := TestSEVDBI(False);
end;


end.