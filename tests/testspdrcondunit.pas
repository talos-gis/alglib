unit testspdrcondunit;
interface
uses Math, Sysutils, Ap, trlinsolve, cholesky, estnorm, spdrcond;

function TestRCondCholesky(Silent : Boolean):Boolean;
function testspdrcondunit_test_silent():Boolean;
function testspdrcondunit_test():Boolean;

implementation

procedure MakeACopy(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);forward;
function InvMatTR(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsunitTriangular : Boolean):Boolean;forward;
function InvMatLU(var A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger):Boolean;forward;
procedure MatLU(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);forward;
function InvMat(var A : TReal2DArray; N : AlglibInteger):Boolean;forward;
procedure RefRCond(const A : TReal2DArray;
     N : AlglibInteger;
     var RC1 : Double;
     var RCInf : Double);forward;


function TestRCondCholesky(Silent : Boolean):Boolean;
var
    L : TReal2DArray;
    A : TReal2DArray;
    A2 : TReal2DArray;
    MinIJ : AlglibInteger;
    P : TInteger1DArray;
    N : AlglibInteger;
    MaxN : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    WasErrors : Boolean;
    Err50 : Boolean;
    Err90 : Boolean;
    ErrLess : Boolean;
    ERC1 : Double;
    ERCInf : Double;
    Q50 : TReal1DArray;
    Q90 : TReal1DArray;
    V : Double;
    Threshold50 : Double;
    Threshold90 : Double;
    HTask : AlglibInteger;
    UpperIn : Boolean;
    i_ : AlglibInteger;
begin
    WasErrors := False;
    Err50 := False;
    Err90 := False;
    ErrLess := False;
    MaxN := 10;
    PassCount := 100;
    Threshold50 := 0.5;
    Threshold90 := 0.1;
    SetLength(Q50, 1+1);
    SetLength(Q90, 1+1);
    
    //
    // process
    //
    N:=1;
    while N<=MaxN do
    begin
        SetLength(L, N-1+1, N-1+1);
        SetLength(A, N-1+1, N-1+1);
        SetLength(A2, N-1+1, N-1+1);
        I:=0;
        while I<=1 do
        begin
            Q50[I] := 0;
            Q90[I] := 0;
            Inc(I);
        end;
        HTask:=0;
        while HTask<=1 do
        begin
            Pass:=1;
            while Pass<=PassCount do
            begin
                UpperIn := HTask=0;
                
                //
                // Prepare task:
                // * A contains full matrix A
                // * L contains its Cholesky factor (upper or lower)
                // * A2 contains upper/lower triangle of A
                //
                I:=0;
                while I<=N-1 do
                begin
                    J:=I+1;
                    while J<=N-1 do
                    begin
                        L[I,J] := 2*RandomReal-1;
                        L[J,I] := L[I,J];
                        Inc(J);
                    end;
                    L[I,I] := 1.1+RandomReal;
                    Inc(I);
                end;
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        MinIJ := Min(I, J);
                        V := 0.0;
                        for i_ := 0 to MinIJ do
                        begin
                            V := V + L[I,i_]*L[i_,J];
                        end;
                        A[I,J] := V;
                        A2[I,J] := V;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        if UpperIn then
                        begin
                            if J<I then
                            begin
                                L[I,J] := 0;
                                A2[I,J] := 0;
                            end;
                        end
                        else
                        begin
                            if I<J then
                            begin
                                L[I,J] := 0;
                                A2[I,J] := 0;
                            end;
                        end;
                        Inc(J);
                    end;
                    Inc(I);
                end;
                RefRCond(A, N, ERC1, ERCInf);
                
                //
                // normal
                //
                V := 1/SPDMatrixRCond(A, N, UpperIn);
                if AP_FP_Greater_Eq(V,Threshold50*ERC1) then
                begin
                    Q50[0] := Q50[0]+AP_Double(1)/PassCount;
                end;
                if AP_FP_Greater_Eq(V,Threshold90*ERC1) then
                begin
                    Q90[0] := Q90[0]+AP_Double(1)/PassCount;
                end;
                ErrLess := ErrLess or AP_FP_Greater(V,ERC1*1.001);
                
                //
                // Cholesky
                //
                V := 1/SPDMatrixCholeskyRCond(L, N, UpperIn);
                if AP_FP_Greater_Eq(V,Threshold50*ERC1) then
                begin
                    Q50[1] := Q50[1]+AP_Double(1)/PassCount;
                end;
                if AP_FP_Greater_Eq(V,Threshold90*ERC1) then
                begin
                    Q90[1] := Q90[1]+AP_Double(1)/PassCount;
                end;
                ErrLess := ErrLess or AP_FP_Greater(V,ERC1*1.001);
                Inc(Pass);
            end;
            Inc(HTask);
        end;
        I:=0;
        while I<=1 do
        begin
            Err50 := Err50 or AP_FP_Less(Q50[I],0.50);
            Err90 := Err90 or AP_FP_Less(Q90[I],0.90);
            Inc(I);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    WasErrors := Err50 or Err90 or ErrLess;
    if  not Silent then
    begin
        Write(Format('TESTING RCOND (SYMMETRIC)'#13#10'',[]));
        Write(Format('50%% within [0.5*cond,cond]:              ',[]));
        if Err50 or ErrLess then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('90%% within [0.1*cond,cond]               ',[]));
        if Err90 or ErrLess then
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
    Result :=  not WasErrors;
end;


(*************************************************************************
Copy
*************************************************************************)
procedure MakeACopy(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    SetLength(B, M-1+1, N-1+1);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            B[I,J] := A[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
triangular inverse
*************************************************************************)
function InvMatTR(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsunitTriangular : Boolean):Boolean;
var
    NOunit : Boolean;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Double;
    AJJ : Double;
    T : TReal1DArray;
    i_ : AlglibInteger;
begin
    Result := True;
    SetLength(T, N-1+1);
    
    //
    // Test the input parameters.
    //
    NOunit :=  not IsunitTriangular;
    if IsUpper then
    begin
        
        //
        // Compute inverse of upper triangular matrix.
        //
        J:=0;
        while J<=N-1 do
        begin
            if NOunit then
            begin
                if AP_FP_Eq(A[J,J],0) then
                begin
                    Result := False;
                    Exit;
                end;
                A[J,J] := 1/A[J,J];
                AJJ := -A[J,J];
            end
            else
            begin
                AJJ := -1;
            end;
            
            //
            // Compute elements 1:j-1 of j-th column.
            //
            if J>0 then
            begin
                for i_ := 0 to J-1 do
                begin
                    T[i_] := A[i_,J];
                end;
                I:=0;
                while I<=J-1 do
                begin
                    if I<J-1 then
                    begin
                        V := APVDotProduct(@A[I][0], I+1, J-1, @T[0], I+1, J-1);
                    end
                    else
                    begin
                        V := 0;
                    end;
                    if NOunit then
                    begin
                        A[I,J] := V+A[I,I]*T[I];
                    end
                    else
                    begin
                        A[I,J] := V+T[I];
                    end;
                    Inc(I);
                end;
                for i_ := 0 to J-1 do
                begin
                    A[i_,J] := AJJ*A[i_,J];
                end;
            end;
            Inc(J);
        end;
    end
    else
    begin
        
        //
        // Compute inverse of lower triangular matrix.
        //
        J:=N-1;
        while J>=0 do
        begin
            if NOunit then
            begin
                if AP_FP_Eq(A[J,J],0) then
                begin
                    Result := False;
                    Exit;
                end;
                A[J,J] := 1/A[J,J];
                AJJ := -A[J,J];
            end
            else
            begin
                AJJ := -1;
            end;
            if J<N-1 then
            begin
                
                //
                // Compute elements j+1:n of j-th column.
                //
                for i_ := J+1 to N-1 do
                begin
                    T[i_] := A[i_,J];
                end;
                I:=J+1;
                while I<=N-1 do
                begin
                    if I>J+1 then
                    begin
                        V := APVDotProduct(@A[I][0], J+1, I-1, @T[0], J+1, I-1);
                    end
                    else
                    begin
                        V := 0;
                    end;
                    if NOunit then
                    begin
                        A[I,J] := V+A[I,I]*T[I];
                    end
                    else
                    begin
                        A[I,J] := V+T[I];
                    end;
                    Inc(I);
                end;
                for i_ := J+1 to N-1 do
                begin
                    A[i_,J] := AJJ*A[i_,J];
                end;
            end;
            Dec(J);
        end;
    end;
end;


(*************************************************************************
LU inverse
*************************************************************************)
function InvMatLU(var A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger):Boolean;
var
    WORK : TReal1DArray;
    I : AlglibInteger;
    IWS : AlglibInteger;
    J : AlglibInteger;
    JB : AlglibInteger;
    JJ : AlglibInteger;
    JP : AlglibInteger;
    V : Double;
    i_ : AlglibInteger;
begin
    Result := True;
    
    //
    // Quick return if possible
    //
    if N=0 then
    begin
        Exit;
    end;
    SetLength(WORK, N-1+1);
    
    //
    // Form inv(U)
    //
    if  not InvMatTR(A, N, True, False) then
    begin
        Result := False;
        Exit;
    end;
    
    //
    // Solve the equation inv(A)*L = inv(U) for inv(A).
    //
    J:=N-1;
    while J>=0 do
    begin
        
        //
        // Copy current column of L to WORK and replace with zeros.
        //
        I:=J+1;
        while I<=N-1 do
        begin
            WORK[I] := A[I,J];
            A[I,J] := 0;
            Inc(I);
        end;
        
        //
        // Compute current column of inv(A).
        //
        if J<N-1 then
        begin
            I:=0;
            while I<=N-1 do
            begin
                V := APVDotProduct(@A[I][0], J+1, N-1, @WORK[0], J+1, N-1);
                A[I,J] := A[I,J]-V;
                Inc(I);
            end;
        end;
        Dec(J);
    end;
    
    //
    // Apply column interchanges.
    //
    J:=N-2;
    while J>=0 do
    begin
        JP := Pivots[J];
        if JP<>J then
        begin
            for i_ := 0 to N-1 do
            begin
                WORK[i_] := A[i_,J];
            end;
            for i_ := 0 to N-1 do
            begin
                A[i_,J] := A[i_,JP];
            end;
            for i_ := 0 to N-1 do
            begin
                A[i_,JP] := WORK[i_];
            end;
        end;
        Dec(J);
    end;
end;


(*************************************************************************
LU decomposition
*************************************************************************)
procedure MatLU(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    JP : AlglibInteger;
    T1 : TReal1DArray;
    s : Double;
    i_ : AlglibInteger;
begin
    SetLength(Pivots, Min(M-1, N-1)+1);
    SetLength(T1, Max(M-1, N-1)+1);
    Assert((M>=0) and (N>=0), 'Error in LUDecomposition: incorrect function arguments');
    
    //
    // Quick return if possible
    //
    if (M=0) or (N=0) then
    begin
        Exit;
    end;
    J:=0;
    while J<=Min(M-1, N-1) do
    begin
        
        //
        // Find pivot and test for singularity.
        //
        JP := J;
        I:=J+1;
        while I<=M-1 do
        begin
            if AP_FP_Greater(AbsReal(A[I,J]),AbsReal(A[JP,J])) then
            begin
                JP := I;
            end;
            Inc(I);
        end;
        Pivots[J] := JP;
        if AP_FP_Neq(A[JP,J],0) then
        begin
            
            //
            //Apply the interchange to rows
            //
            if JP<>J then
            begin
                APVMove(@T1[0], 0, N-1, @A[J][0], 0, N-1);
                APVMove(@A[J][0], 0, N-1, @A[JP][0], 0, N-1);
                APVMove(@A[JP][0], 0, N-1, @T1[0], 0, N-1);
            end;
            
            //
            //Compute elements J+1:M of J-th column.
            //
            if J<M then
            begin
                JP := J+1;
                S := 1/A[J,J];
                for i_ := JP to M-1 do
                begin
                    A[i_,J] := S*A[i_,J];
                end;
            end;
        end;
        if J<Min(M, N)-1 then
        begin
            
            //
            //Update trailing submatrix.
            //
            JP := J+1;
            I:=J+1;
            while I<=M-1 do
            begin
                S := A[I,J];
                APVSub(@A[I][0], JP, N-1, @A[J][0], JP, N-1, S);
                Inc(I);
            end;
        end;
        Inc(J);
    end;
end;


(*************************************************************************
Matrix inverse
*************************************************************************)
function InvMat(var A : TReal2DArray; N : AlglibInteger):Boolean;
var
    Pivots : TInteger1DArray;
begin
    MatLU(A, N, N, Pivots);
    Result := InvMatLU(A, Pivots, N);
end;


(*************************************************************************
reference RCond
*************************************************************************)
procedure RefRCond(const A : TReal2DArray;
     N : AlglibInteger;
     var RC1 : Double;
     var RCInf : Double);
var
    InvA : TReal2DArray;
    Nrm1A : Double;
    NrmInfA : Double;
    Nrm1InvA : Double;
    NrmInfInvA : Double;
    V : Double;
    K : AlglibInteger;
    I : AlglibInteger;
begin
    
    //
    // inv A
    //
    MakeACopy(A, N, N, InvA);
    if  not InvMat(InvA, N) then
    begin
        RC1 := 0;
        RCInf := 0;
        Exit;
    end;
    
    //
    // norm A
    //
    Nrm1A := 0;
    NrmInfA := 0;
    K:=0;
    while K<=N-1 do
    begin
        V := 0;
        I:=0;
        while I<=N-1 do
        begin
            V := V+AbsReal(A[I,K]);
            Inc(I);
        end;
        Nrm1A := Max(Nrm1A, V);
        V := 0;
        I:=0;
        while I<=N-1 do
        begin
            V := V+AbsReal(A[K,I]);
            Inc(I);
        end;
        NrmInfA := Max(NrmInfA, V);
        Inc(K);
    end;
    
    //
    // norm inv A
    //
    Nrm1InvA := 0;
    NrmInfInvA := 0;
    K:=0;
    while K<=N-1 do
    begin
        V := 0;
        I:=0;
        while I<=N-1 do
        begin
            V := V+AbsReal(InvA[I,K]);
            Inc(I);
        end;
        Nrm1InvA := Max(Nrm1InvA, V);
        V := 0;
        I:=0;
        while I<=N-1 do
        begin
            V := V+AbsReal(InvA[K,I]);
            Inc(I);
        end;
        NrmInfInvA := Max(NrmInfInvA, V);
        Inc(K);
    end;
    
    //
    // result
    //
    RC1 := Nrm1InvA*Nrm1A;
    RCInf := NrmInfInvA*NrmInfA;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testspdrcondunit_test_silent():Boolean;
begin
    Result := TestRCondCholesky(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testspdrcondunit_test():Boolean;
begin
    Result := TestRCondCholesky(False);
end;


end.