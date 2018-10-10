unit testrcondunit;
interface
uses Math, Sysutils, Ap, lu, trlinsolve, rcond;

function TestRCond(Silent : Boolean):Boolean;
function testrcondunit_test_silent():Boolean;
function testrcondunit_test():Boolean;

implementation

procedure MakeACopy(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);forward;
procedure GenerateRandomMatrixCond(var A0 : TReal2DArray;
     N : AlglibInteger;
     C : Double);forward;
function InvMatTR(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsunitTriangular : Boolean):Boolean;forward;
function InvMatLU(var A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger):Boolean;forward;
function InvMat(var A : TReal2DArray; N : AlglibInteger):Boolean;forward;
procedure RefRCond(const A : TReal2DArray;
     N : AlglibInteger;
     var RC1 : Double;
     var RCInf : Double);forward;


function TestRCond(Silent : Boolean):Boolean;
var
    A : TReal2DArray;
    LUA : TReal2DArray;
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
begin
    WasErrors := False;
    Err50 := False;
    Err90 := False;
    ErrLess := False;
    MaxN := 10;
    PassCount := 100;
    Threshold50 := 0.5;
    Threshold90 := 0.1;
    SetLength(Q50, 3+1);
    SetLength(Q90, 3+1);
    
    //
    // process
    //
    N:=1;
    while N<=MaxN do
    begin
        SetLength(A, N-1+1, N-1+1);
        I:=0;
        while I<=3 do
        begin
            Q50[I] := 0;
            Q90[I] := 0;
            Inc(I);
        end;
        Pass:=1;
        while Pass<=PassCount do
        begin
            GenerateRandomMatrixCond(A, N, Exp(RandomReal*Ln(1000)));
            MakeACopy(A, N, N, LUA);
            RMatrixLU(LUA, N, N, P);
            RefRCond(A, N, ERC1, ERCInf);
            
            //
            // 1-norm, normal
            //
            V := 1/RMatrixRCond1(A, N);
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
            // 1-norm, LU
            //
            V := 1/RMatrixLURCond1(LUA, N);
            if AP_FP_Greater_Eq(V,Threshold50*ERC1) then
            begin
                Q50[1] := Q50[1]+AP_Double(1)/PassCount;
            end;
            if AP_FP_Greater_Eq(V,Threshold90*ERC1) then
            begin
                Q90[1] := Q90[1]+AP_Double(1)/PassCount;
            end;
            ErrLess := ErrLess or AP_FP_Greater(V,ERC1*1.001);
            
            //
            // Inf-norm, normal
            //
            V := 1/RMatrixRCondInf(A, N);
            if AP_FP_Greater_Eq(V,Threshold50*ERCInf) then
            begin
                Q50[2] := Q50[2]+AP_Double(1)/PassCount;
            end;
            if AP_FP_Greater_Eq(V,Threshold90*ERCInf) then
            begin
                Q90[2] := Q90[2]+AP_Double(1)/PassCount;
            end;
            ErrLess := ErrLess or AP_FP_Greater(V,ERCInf*1.001);
            
            //
            // Inf-norm, LU
            //
            V := 1/RMatrixLURCondInf(LUA, N);
            if AP_FP_Greater_Eq(V,Threshold50*ERCInf) then
            begin
                Q50[3] := Q50[3]+AP_Double(1)/PassCount;
            end;
            if AP_FP_Greater_Eq(V,Threshold90*ERCInf) then
            begin
                Q90[3] := Q90[3]+AP_Double(1)/PassCount;
            end;
            ErrLess := ErrLess or AP_FP_Greater(V,ERCInf*1.001);
            Inc(Pass);
        end;
        I:=0;
        while I<=3 do
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
        Write(Format('TESTING RCOND (REAL)'#13#10'',[]));
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
Generate matrix with given condition number C (2-norm)
*************************************************************************)
procedure GenerateRandomMatrixCond(var A0 : TReal2DArray;
     N : AlglibInteger;
     C : Double);
var
    T : Double;
    Lambda : Double;
    S : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    U1 : Double;
    U2 : Double;
    W : TReal1DArray;
    V : TReal1DArray;
    SM : Double;
    L1 : Double;
    L2 : Double;
    A : TReal2DArray;
begin
    if (N<=0) or AP_FP_Less(C,1) then
    begin
        Exit;
    end;
    SetLength(A, N+1, N+1);
    SetLength(A0, N-1+1, N-1+1);
    if N=1 then
    begin
        A0[0,0] := 1;
        Exit;
    end;
    SetLength(W, N+1);
    SetLength(V, N+1);
    
    //
    // N>=2, prepare A
    //
    I:=1;
    while I<=N do
    begin
        J:=1;
        while J<=N do
        begin
            A[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    L1 := 0;
    L2 := Ln(1/C);
    A[1,1] := Exp(L1);
    I:=2;
    while I<=N-1 do
    begin
        A[I,I] := Exp(RandomReal*(L2-L1)+L1);
        Inc(I);
    end;
    A[N,N] := Exp(L2);
    
    //
    // Multiply A by random Q from the right (using Stewart algorithm)
    //
    S:=2;
    while S<=N do
    begin
        
        //
        // Prepare v and Lambda = v'*v
        //
        repeat
            I := 1;
            while i<=S do
            begin
                U1 := 2*RandomReal-1;
                U2 := 2*RandomReal-1;
                SM := U1*U1+U2*U2;
                if AP_FP_Eq(SM,0) or AP_FP_Greater(SM,1) then
                begin
                    Continue;
                end;
                SM := Sqrt(-2*Ln(SM)/SM);
                V[I] := U1*SM;
                if I+1<=S then
                begin
                    V[I+1] := U2*SM;
                end;
                I := I+2;
            end;
            Lambda := APVDotProduct(@V[0], 1, S, @V[0], 1, S);
        until AP_FP_Neq(Lambda,0);
        Lambda := 2/Lambda;
        
        //
        // A * (I - 2 vv'/v'v ) =
        //   = A - (2/v'v) * A * v * v' =
        //   = A - (2/v'v) * w * v'
        //  where w = Av
        //
        // Note that size of A is SxS, not SxN
        // due to A structure!!!
        //
        I:=1;
        while I<=S do
        begin
            T := APVDotProduct(@A[I][0], 1, S, @V[0], 1, S);
            W[I] := T;
            Inc(I);
        end;
        I:=1;
        while I<=S do
        begin
            T := W[I]*Lambda;
            APVSub(@A[I][0], 1, S, @V[0], 1, S, T);
            Inc(I);
        end;
        Inc(S);
    end;
    
    //
    // Multiply A by random Q from the left (using Stewart algorithm)
    //
    S:=2;
    while S<=N do
    begin
        
        //
        // Prepare v and Lambda = v'*v
        //
        repeat
            I := 1;
            while i<=S do
            begin
                U1 := 2*RandomReal-1;
                U2 := 2*RandomReal-1;
                SM := U1*U1+U2*U2;
                if AP_FP_Eq(SM,0) or AP_FP_Greater(SM,1) then
                begin
                    Continue;
                end;
                SM := Sqrt(-2*Ln(SM)/SM);
                V[I] := U1*SM;
                if I+1<=S then
                begin
                    V[I+1] := U2*SM;
                end;
                I := I+2;
            end;
            Lambda := APVDotProduct(@V[0], 1, S, @V[0], 1, S);
        until AP_FP_Neq(Lambda,0);
        Lambda := 2/Lambda;
        
        //
        // (I - 2 vv'/v'v ) * A =
        //   = A - (2/v'v) * v * v' * A =
        //   = A - (2/v'v) * v * w
        //  where w = v'A
        //
        // Note that size of A is SxN, not SxS!!!
        //
        I:=1;
        while I<=N do
        begin
            W[I] := 0;
            Inc(I);
        end;
        I:=1;
        while I<=S do
        begin
            T := V[I];
            APVAdd(@W[0], 1, N, @A[I][0], 1, N, T);
            Inc(I);
        end;
        I:=1;
        while I<=S do
        begin
            T := V[I]*Lambda;
            APVSub(@A[I][0], 1, N, @W[0], 1, N, T);
            Inc(I);
        end;
        Inc(S);
    end;
    
    //
    //
    //
    I:=1;
    while I<=N do
    begin
        J:=1;
        while J<=N do
        begin
            A0[I-1,J-1] := A[I,J];
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
Matrix inverse
*************************************************************************)
function InvMat(var A : TReal2DArray; N : AlglibInteger):Boolean;
var
    Pivots : TInteger1DArray;
begin
    RMatrixLU(A, N, N, Pivots);
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
function testrcondunit_test_silent():Boolean;
begin
    Result := TestRCond(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testrcondunit_test():Boolean;
begin
    Result := TestRCond(False);
end;


end.