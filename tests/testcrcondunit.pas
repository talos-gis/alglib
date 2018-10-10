unit testcrcondunit;
interface
uses Math, Sysutils, Ap, clu, ctrlinsolve, crcond;

function TestCRCond(Silent : Boolean):Boolean;
function testcrcondunit_test_silent():Boolean;
function testcrcondunit_test():Boolean;

implementation

procedure MakeACopy(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TComplex2DArray);forward;
procedure GenerateRandomMatrix(var A0 : TComplex2DArray;
     N : AlglibInteger);forward;
function InvMatTR(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsunitTriangular : Boolean):Boolean;forward;
function InvMatLU(var A : TComplex2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger):Boolean;forward;
function InvMat(var A : TComplex2DArray; N : AlglibInteger):Boolean;forward;
procedure RefRCond(const A : TComplex2DArray;
     N : AlglibInteger;
     var RC1 : Double;
     var RCInf : Double);forward;


function TestCRCond(Silent : Boolean):Boolean;
var
    A : TComplex2DArray;
    LUA : TComplex2DArray;
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
            GenerateRandomMatrix(A, N);
            MakeACopy(A, N, N, LUA);
            CMatrixLU(LUA, N, N, P);
            RefRCond(A, N, ERC1, ERCInf);
            
            //
            // 1-norm, normal
            //
            V := 1/CMatrixRCond1(A, N);
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
            V := 1/CMatrixLURCond1(LUA, N);
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
            V := 1/CMatrixRCondInf(A, N);
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
            V := 1/CMatrixLURCondInf(LUA, N);
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
        Write(Format('TESTING RCOND (COMPLEX)'#13#10'',[]));
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
procedure MakeACopy(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TComplex2DArray);
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
procedure GenerateRandomMatrix(var A0 : TComplex2DArray; N : AlglibInteger);
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
            A0[I,J].X := 2*RandomReal-1;
            A0[I,J].Y := 2*RandomReal-1;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
triangular inverse
*************************************************************************)
function InvMatTR(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsunitTriangular : Boolean):Boolean;
var
    NOunit : Boolean;
    I : AlglibInteger;
    J : AlglibInteger;
    V : Complex;
    AJJ : Complex;
    T : TComplex1DArray;
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
                if C_EqualR(A[J,J],0) then
                begin
                    Result := False;
                    Exit;
                end;
                A[J,J] := C_RDiv(1,A[J,J]);
                AJJ := C_Opposite(A[J,J]);
            end
            else
            begin
                AJJ := C_Complex(-1);
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
                        V := C_Complex(0.0);
                        for i_ := I+1 to J-1 do
                        begin
                            V := C_Add(V,C_Mul(A[I,i_],T[i_]));
                        end;
                    end
                    else
                    begin
                        V := C_Complex(0);
                    end;
                    if NOunit then
                    begin
                        A[I,J] := C_Add(V,C_Mul(A[I,I],T[I]));
                    end
                    else
                    begin
                        A[I,J] := C_Add(V,T[I]);
                    end;
                    Inc(I);
                end;
                for i_ := 0 to J-1 do
                begin
                    A[i_,J] := C_Mul(AJJ, A[i_,J]);
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
                if C_EqualR(A[J,J],0) then
                begin
                    Result := False;
                    Exit;
                end;
                A[J,J] := C_RDiv(1,A[J,J]);
                AJJ := C_Opposite(A[J,J]);
            end
            else
            begin
                AJJ := C_Complex(-1);
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
                        V := C_Complex(0.0);
                        for i_ := J+1 to I-1 do
                        begin
                            V := C_Add(V,C_Mul(A[I,i_],T[i_]));
                        end;
                    end
                    else
                    begin
                        V := C_Complex(0);
                    end;
                    if NOunit then
                    begin
                        A[I,J] := C_Add(V,C_Mul(A[I,I],T[I]));
                    end
                    else
                    begin
                        A[I,J] := C_Add(V,T[I]);
                    end;
                    Inc(I);
                end;
                for i_ := J+1 to N-1 do
                begin
                    A[i_,J] := C_Mul(AJJ, A[i_,J]);
                end;
            end;
            Dec(J);
        end;
    end;
end;


(*************************************************************************
LU inverse
*************************************************************************)
function InvMatLU(var A : TComplex2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger):Boolean;
var
    WORK : TComplex1DArray;
    I : AlglibInteger;
    IWS : AlglibInteger;
    J : AlglibInteger;
    JB : AlglibInteger;
    JJ : AlglibInteger;
    JP : AlglibInteger;
    V : Complex;
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
            A[I,J] := C_Complex(0);
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
                V := C_Complex(0.0);
                for i_ := J+1 to N-1 do
                begin
                    V := C_Add(V,C_Mul(A[I,i_],WORK[i_]));
                end;
                A[I,J] := C_Sub(A[I,J],V);
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
function InvMat(var A : TComplex2DArray; N : AlglibInteger):Boolean;
var
    Pivots : TInteger1DArray;
begin
    CMatrixLU(A, N, N, Pivots);
    Result := InvMatLU(A, Pivots, N);
end;


(*************************************************************************
reference RCond
*************************************************************************)
procedure RefRCond(const A : TComplex2DArray;
     N : AlglibInteger;
     var RC1 : Double;
     var RCInf : Double);
var
    InvA : TComplex2DArray;
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
            V := V+AbsComplex(A[I,K]);
            Inc(I);
        end;
        Nrm1A := Max(Nrm1A, V);
        V := 0;
        I:=0;
        while I<=N-1 do
        begin
            V := V+AbsComplex(A[K,I]);
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
            V := V+AbsComplex(InvA[I,K]);
            Inc(I);
        end;
        Nrm1InvA := Max(Nrm1InvA, V);
        V := 0;
        I:=0;
        while I<=N-1 do
        begin
            V := V+AbsComplex(InvA[K,I]);
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
function testcrcondunit_test_silent():Boolean;
begin
    Result := TestCRCond(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testcrcondunit_test():Boolean;
begin
    Result := TestCRCond(False);
end;


end.