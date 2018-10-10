unit testrcondunit;
interface
uses Math, Sysutils, Ap, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac, trlinsolve, safesolve, rcond;

function TestRCond(Silent : Boolean):Boolean;
function testrcondunit_test_silent():Boolean;
function testrcondunit_test():Boolean;

implementation

const
    Threshold50 = 0.25;
    Threshold90 = 0.10;

procedure RMatrixMakeACopy(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);forward;
procedure RMatrixDropHalf(var A : TReal2DArray;
     N : AlglibInteger;
     DropLower : Boolean);forward;
procedure CMatrixDropHalf(var A : TComplex2DArray;
     N : AlglibInteger;
     DropLower : Boolean);forward;
procedure RMatrixGenZero(var A0 : TReal2DArray; N : AlglibInteger);forward;
function RMatrixInvMatTR(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsunitTriangular : Boolean):Boolean;forward;
function RMatrixInvMatLU(var A : TReal2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger):Boolean;forward;
function RMatrixInvMat(var A : TReal2DArray;
     N : AlglibInteger):Boolean;forward;
procedure RMatrixRefRCond(const A : TReal2DArray;
     N : AlglibInteger;
     var RC1 : Double;
     var RCInf : Double);forward;
procedure CMatrixMakeACopy(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TComplex2DArray);forward;
procedure CMatrixGenZero(var A0 : TComplex2DArray; N : AlglibInteger);forward;
function CMatrixInvMatTR(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsunitTriangular : Boolean):Boolean;forward;
function CMatrixInvMatLU(var A : TComplex2DArray;
     const Pivots : TInteger1DArray;
     N : AlglibInteger):Boolean;forward;
function CMatrixInvMat(var A : TComplex2DArray;
     N : AlglibInteger):Boolean;forward;
procedure CMatrixRefRCond(const A : TComplex2DArray;
     N : AlglibInteger;
     var RC1 : Double;
     var RCInf : Double);forward;
function TestRMatrixTRRCond(MaxN : AlglibInteger;
     PassCount : AlglibInteger):Boolean;forward;
function TestCMatrixTRRCond(MaxN : AlglibInteger;
     PassCount : AlglibInteger):Boolean;forward;
function TestRMatrixRCond(MaxN : AlglibInteger;
     PassCount : AlglibInteger):Boolean;forward;
function TestSPDMatrixRCond(MaxN : AlglibInteger;
     PassCount : AlglibInteger):Boolean;forward;
function TestCMatrixRCond(MaxN : AlglibInteger;
     PassCount : AlglibInteger):Boolean;forward;
function TestHPDMatrixRCond(MaxN : AlglibInteger;
     PassCount : AlglibInteger):Boolean;forward;


function TestRCond(Silent : Boolean):Boolean;
var
    MaxN : AlglibInteger;
    PassCount : AlglibInteger;
    WasErrors : Boolean;
    RTRErr : Boolean;
    CTRErr : Boolean;
    RErr : Boolean;
    CErr : Boolean;
    SPDErr : Boolean;
    HPDErr : Boolean;
begin
    MaxN := 10;
    PassCount := 100;
    
    //
    // report
    //
    RTRErr :=  not TestRMatrixTRRCond(MaxN, PassCount);
    CTRErr :=  not TestCMatrixTRRCond(MaxN, PassCount);
    RErr :=  not TestRMatrixRCond(MaxN, PassCount);
    CErr :=  not TestCMatrixRCond(MaxN, PassCount);
    SPDErr :=  not TestSPDMatrixRCond(MaxN, PassCount);
    HPDErr :=  not TestHPDMatrixRCond(MaxN, PassCount);
    WasErrors := RTRErr or CTRErr or RErr or CErr or SPDErr or HPDErr;
    if  not Silent then
    begin
        Write(Format('TESTING RCOND'#13#10'',[]));
        Write(Format('REAL TRIANGULAR:                         ',[]));
        if  not RTRErr then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('COMPLEX TRIANGULAR:                      ',[]));
        if  not CTRErr then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('REAL:                                    ',[]));
        if  not RErr then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('SPD:                                     ',[]));
        if  not SPDErr then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('HPD:                                     ',[]));
        if  not HPDErr then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('COMPLEX:                                 ',[]));
        if  not CErr then
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
Copy
*************************************************************************)
procedure RMatrixMakeACopy(const A : TReal2DArray;
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
Drops upper or lower half of the matrix - fills it by special pattern
which may be used later to ensure that this part wasn't changed
*************************************************************************)
procedure RMatrixDropHalf(var A : TReal2DArray;
     N : AlglibInteger;
     DropLower : Boolean);
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
            if DropLower and (I>J) or  not DropLower and (I<J) then
            begin
                A[I,J] := 1+2*I+3*J;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Drops upper or lower half of the matrix - fills it by special pattern
which may be used later to ensure that this part wasn't changed
*************************************************************************)
procedure CMatrixDropHalf(var A : TComplex2DArray;
     N : AlglibInteger;
     DropLower : Boolean);
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
            if DropLower and (I>J) or  not DropLower and (I<J) then
            begin
                A[I,J] := C_Complex(1+2*I+3*J);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Generate matrix with given condition number C (2-norm)
*************************************************************************)
procedure RMatrixGenZero(var A0 : TReal2DArray; N : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    SetLength(A0, N, N);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            A0[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
triangular inverse
*************************************************************************)
function RMatrixInvMatTR(var A : TReal2DArray;
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
function RMatrixInvMatLU(var A : TReal2DArray;
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
    if  not RMatrixInvMatTR(A, N, True, False) then
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
function RMatrixInvMat(var A : TReal2DArray; N : AlglibInteger):Boolean;
var
    Pivots : TInteger1DArray;
begin
    RMatrixLU(A, N, N, Pivots);
    Result := RMatrixInvMatLU(A, Pivots, N);
end;


(*************************************************************************
reference RCond
*************************************************************************)
procedure RMatrixRefRCond(const A : TReal2DArray;
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
    RMatrixMakeACopy(A, N, N, InvA);
    if  not RMatrixInvMat(InvA, N) then
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
Copy
*************************************************************************)
procedure CMatrixMakeACopy(const A : TComplex2DArray;
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
procedure CMatrixGenZero(var A0 : TComplex2DArray; N : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    SetLength(A0, N, N);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            A0[I,J] := C_Complex(0);
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
triangular inverse
*************************************************************************)
function CMatrixInvMatTR(var A : TComplex2DArray;
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
function CMatrixInvMatLU(var A : TComplex2DArray;
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
    if  not CMatrixInvMatTR(A, N, True, False) then
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
function CMatrixInvMat(var A : TComplex2DArray; N : AlglibInteger):Boolean;
var
    Pivots : TInteger1DArray;
begin
    CMatrixLU(A, N, N, Pivots);
    Result := CMatrixInvMatLU(A, Pivots, N);
end;


(*************************************************************************
reference RCond
*************************************************************************)
procedure CMatrixRefRCond(const A : TComplex2DArray;
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
    CMatrixMakeACopy(A, N, N, InvA);
    if  not CMatrixInvMat(InvA, N) then
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
Returns True for successful test, False - for failed test
*************************************************************************)
function TestRMatrixTRRCond(MaxN : AlglibInteger;
     PassCount : AlglibInteger):Boolean;
var
    A : TReal2DArray;
    EA : TReal2DArray;
    P : TInteger1DArray;
    N : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
    Pass : AlglibInteger;
    Err50 : Boolean;
    Err90 : Boolean;
    ErrSpec : Boolean;
    ErrLess : Boolean;
    ERC1 : Double;
    ERCInf : Double;
    Q50 : TReal1DArray;
    Q90 : TReal1DArray;
    V : Double;
    IsUpper : Boolean;
    IsUnit : Boolean;
begin
    Err50 := False;
    Err90 := False;
    ErrLess := False;
    ErrSpec := False;
    SetLength(Q50, 2);
    SetLength(Q90, 2);
    N:=1;
    while N<=MaxN do
    begin
        
        //
        // special test for zero matrix
        //
        RMatrixGenZero(A, N);
        ErrSpec := ErrSpec or AP_FP_Neq(RMatrixTRRCond1(A, N, AP_FP_Greater(RandomReal,0.5), False),0);
        ErrSpec := ErrSpec or AP_FP_Neq(RMatrixTRRCondInf(A, N, AP_FP_Greater(RandomReal,0.5), False),0);
        
        //
        // general test
        //
        SetLength(A, N, N);
        I:=0;
        while I<=1 do
        begin
            Q50[I] := 0;
            Q90[I] := 0;
            Inc(I);
        end;
        Pass:=1;
        while Pass<=PassCount do
        begin
            IsUpper := AP_FP_Greater(RandomReal,0.5);
            IsUnit := AP_FP_Greater(RandomReal,0.5);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := RandomReal-0.5;
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=N-1 do
            begin
                A[I,I] := 1+RandomReal;
                Inc(I);
            end;
            RMatrixMakeACopy(A, N, N, EA);
            I:=0;
            while I<=N-1 do
            begin
                if IsUpper then
                begin
                    J1 := 0;
                    J2 := I-1;
                end
                else
                begin
                    J1 := I+1;
                    J2 := N-1;
                end;
                J:=J1;
                while J<=J2 do
                begin
                    EA[I,J] := 0;
                    Inc(J);
                end;
                if IsUnit then
                begin
                    EA[I,I] := 1;
                end;
                Inc(I);
            end;
            RMatrixRefRCond(EA, N, ERC1, ERCInf);
            
            //
            // 1-norm
            //
            V := 1/RMatrixTRRCond1(A, N, IsUpper, IsUnit);
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
            // Inf-norm
            //
            V := 1/RMatrixTRRCondInf(A, N, IsUpper, IsUnit);
            if AP_FP_Greater_Eq(V,Threshold50*ERCInf) then
            begin
                Q50[1] := Q50[1]+AP_Double(1)/PassCount;
            end;
            if AP_FP_Greater_Eq(V,Threshold90*ERCInf) then
            begin
                Q90[1] := Q90[1]+AP_Double(1)/PassCount;
            end;
            ErrLess := ErrLess or AP_FP_Greater(V,ERCInf*1.001);
            Inc(Pass);
        end;
        I:=0;
        while I<=1 do
        begin
            Err50 := Err50 or AP_FP_Less(Q50[I],0.50);
            Err90 := Err90 or AP_FP_Less(Q90[I],0.90);
            Inc(I);
        end;
        
        //
        // degenerate matrix test
        //
        if N>=3 then
        begin
            SetLength(A, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := 0.0;
                    Inc(J);
                end;
                Inc(I);
            end;
            A[0,0] := 1;
            A[N-1,N-1] := 1;
            ErrSpec := ErrSpec or AP_FP_Neq(RMatrixTRRCond1(A, N, AP_FP_Greater(RandomReal,0.5), False),0);
            ErrSpec := ErrSpec or AP_FP_Neq(RMatrixTRRCondInf(A, N, AP_FP_Greater(RandomReal,0.5), False),0);
        end;
        
        //
        // near-degenerate matrix test
        //
        if N>=2 then
        begin
            SetLength(A, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := 0.0;
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=N-1 do
            begin
                A[I,I] := 1;
                Inc(I);
            end;
            I := RandomInteger(N);
            A[I,I] := 0.1*MaxRealNumber;
            ErrSpec := ErrSpec or AP_FP_Neq(RMatrixTRRCond1(A, N, AP_FP_Greater(RandomReal,0.5), False),0);
            ErrSpec := ErrSpec or AP_FP_Neq(RMatrixTRRCondInf(A, N, AP_FP_Greater(RandomReal,0.5), False),0);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    Result :=  not (Err50 or Err90 or ErrLess or ErrSpec);
end;


(*************************************************************************
Returns True for successful test, False - for failed test
*************************************************************************)
function TestCMatrixTRRCond(MaxN : AlglibInteger;
     PassCount : AlglibInteger):Boolean;
var
    A : TComplex2DArray;
    EA : TComplex2DArray;
    P : TInteger1DArray;
    N : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
    Pass : AlglibInteger;
    Err50 : Boolean;
    Err90 : Boolean;
    ErrSpec : Boolean;
    ErrLess : Boolean;
    ERC1 : Double;
    ERCInf : Double;
    Q50 : TReal1DArray;
    Q90 : TReal1DArray;
    V : Double;
    IsUpper : Boolean;
    IsUnit : Boolean;
begin
    Err50 := False;
    Err90 := False;
    ErrLess := False;
    ErrSpec := False;
    SetLength(Q50, 2);
    SetLength(Q90, 2);
    N:=1;
    while N<=MaxN do
    begin
        
        //
        // special test for zero matrix
        //
        CMatrixGenZero(A, N);
        ErrSpec := ErrSpec or AP_FP_Neq(CMatrixTRRCond1(A, N, AP_FP_Greater(RandomReal,0.5), False),0);
        ErrSpec := ErrSpec or AP_FP_Neq(CMatrixTRRCondInf(A, N, AP_FP_Greater(RandomReal,0.5), False),0);
        
        //
        // general test
        //
        SetLength(A, N, N);
        I:=0;
        while I<=1 do
        begin
            Q50[I] := 0;
            Q90[I] := 0;
            Inc(I);
        end;
        Pass:=1;
        while Pass<=PassCount do
        begin
            IsUpper := AP_FP_Greater(RandomReal,0.5);
            IsUnit := AP_FP_Greater(RandomReal,0.5);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J].X := RandomReal-0.5;
                    A[I,J].Y := RandomReal-0.5;
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=N-1 do
            begin
                A[I,I].X := 1+RandomReal;
                A[I,I].Y := 1+RandomReal;
                Inc(I);
            end;
            CMatrixMakeACopy(A, N, N, EA);
            I:=0;
            while I<=N-1 do
            begin
                if IsUpper then
                begin
                    J1 := 0;
                    J2 := I-1;
                end
                else
                begin
                    J1 := I+1;
                    J2 := N-1;
                end;
                J:=J1;
                while J<=J2 do
                begin
                    EA[I,J] := C_Complex(0);
                    Inc(J);
                end;
                if IsUnit then
                begin
                    EA[I,I] := C_Complex(1);
                end;
                Inc(I);
            end;
            CMatrixRefRCond(EA, N, ERC1, ERCInf);
            
            //
            // 1-norm
            //
            V := 1/CMatrixTRRCond1(A, N, IsUpper, IsUnit);
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
            // Inf-norm
            //
            V := 1/CMatrixTRRCondInf(A, N, IsUpper, IsUnit);
            if AP_FP_Greater_Eq(V,Threshold50*ERCInf) then
            begin
                Q50[1] := Q50[1]+AP_Double(1)/PassCount;
            end;
            if AP_FP_Greater_Eq(V,Threshold90*ERCInf) then
            begin
                Q90[1] := Q90[1]+AP_Double(1)/PassCount;
            end;
            ErrLess := ErrLess or AP_FP_Greater(V,ERCInf*1.001);
            Inc(Pass);
        end;
        I:=0;
        while I<=1 do
        begin
            Err50 := Err50 or AP_FP_Less(Q50[I],0.50);
            Err90 := Err90 or AP_FP_Less(Q90[I],0.90);
            Inc(I);
        end;
        
        //
        // degenerate matrix test
        //
        if N>=3 then
        begin
            SetLength(A, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := C_Complex(0.0);
                    Inc(J);
                end;
                Inc(I);
            end;
            A[0,0] := C_Complex(1);
            A[N-1,N-1] := C_Complex(1);
            ErrSpec := ErrSpec or AP_FP_Neq(CMatrixTRRCond1(A, N, AP_FP_Greater(RandomReal,0.5), False),0);
            ErrSpec := ErrSpec or AP_FP_Neq(CMatrixTRRCondInf(A, N, AP_FP_Greater(RandomReal,0.5), False),0);
        end;
        
        //
        // near-degenerate matrix test
        //
        if N>=2 then
        begin
            SetLength(A, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := C_Complex(0.0);
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=N-1 do
            begin
                A[I,I] := C_Complex(1);
                Inc(I);
            end;
            I := RandomInteger(N);
            A[I,I] := C_Complex(0.1*MaxRealNumber);
            ErrSpec := ErrSpec or AP_FP_Neq(CMatrixTRRCond1(A, N, AP_FP_Greater(RandomReal,0.5), False),0);
            ErrSpec := ErrSpec or AP_FP_Neq(CMatrixTRRCondInf(A, N, AP_FP_Greater(RandomReal,0.5), False),0);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    Result :=  not (Err50 or Err90 or ErrLess or ErrSpec);
end;


(*************************************************************************
Returns True for successful test, False - for failed test
*************************************************************************)
function TestRMatrixRCond(MaxN : AlglibInteger;
     PassCount : AlglibInteger):Boolean;
var
    A : TReal2DArray;
    LUA : TReal2DArray;
    P : TInteger1DArray;
    N : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    Pass : AlglibInteger;
    Err50 : Boolean;
    Err90 : Boolean;
    ErrSpec : Boolean;
    ErrLess : Boolean;
    ERC1 : Double;
    ERCInf : Double;
    Q50 : TReal1DArray;
    Q90 : TReal1DArray;
    V : Double;
begin
    Err50 := False;
    Err90 := False;
    ErrLess := False;
    ErrSpec := False;
    SetLength(Q50, 3+1);
    SetLength(Q90, 3+1);
    N:=1;
    while N<=MaxN do
    begin
        
        //
        // special test for zero matrix
        //
        RMatrixGenZero(A, N);
        RMatrixMakeACopy(A, N, N, LUA);
        RMatrixLU(LUA, N, N, P);
        ErrSpec := ErrSpec or AP_FP_Neq(RMatrixRCond1(A, N),0);
        ErrSpec := ErrSpec or AP_FP_Neq(RMatrixRCondInf(A, N),0);
        ErrSpec := ErrSpec or AP_FP_Neq(RMatrixLURCond1(LUA, N),0);
        ErrSpec := ErrSpec or AP_FP_Neq(RMatrixLURCondInf(LUA, N),0);
        
        //
        // general test
        //
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
            RMatrixRndCond(N, Exp(RandomReal*Ln(1000)), A);
            RMatrixMakeACopy(A, N, N, LUA);
            RMatrixLU(LUA, N, N, P);
            RMatrixRefRCond(A, N, ERC1, ERCInf);
            
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
        
        //
        // degenerate matrix test
        //
        if N>=3 then
        begin
            SetLength(A, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := 0.0;
                    Inc(J);
                end;
                Inc(I);
            end;
            A[0,0] := 1;
            A[N-1,N-1] := 1;
            ErrSpec := ErrSpec or AP_FP_Neq(RMatrixRCond1(A, N),0);
            ErrSpec := ErrSpec or AP_FP_Neq(RMatrixRCondInf(A, N),0);
            ErrSpec := ErrSpec or AP_FP_Neq(RMatrixLURCond1(A, N),0);
            ErrSpec := ErrSpec or AP_FP_Neq(RMatrixLURCondInf(A, N),0);
        end;
        
        //
        // near-degenerate matrix test
        //
        if N>=2 then
        begin
            SetLength(A, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := 0.0;
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=N-1 do
            begin
                A[I,I] := 1;
                Inc(I);
            end;
            I := RandomInteger(N);
            A[I,I] := 0.1*MaxRealNumber;
            ErrSpec := ErrSpec or AP_FP_Neq(RMatrixRCond1(A, N),0);
            ErrSpec := ErrSpec or AP_FP_Neq(RMatrixRCondInf(A, N),0);
            ErrSpec := ErrSpec or AP_FP_Neq(RMatrixLURCond1(A, N),0);
            ErrSpec := ErrSpec or AP_FP_Neq(RMatrixLURCondInf(A, N),0);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    Result :=  not (Err50 or Err90 or ErrLess or ErrSpec);
end;


(*************************************************************************
Returns True for successful test, False - for failed test
*************************************************************************)
function TestSPDMatrixRCond(MaxN : AlglibInteger;
     PassCount : AlglibInteger):Boolean;
var
    A : TReal2DArray;
    CHA : TReal2DArray;
    P : TInteger1DArray;
    N : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    Pass : AlglibInteger;
    Err50 : Boolean;
    Err90 : Boolean;
    ErrSpec : Boolean;
    ErrLess : Boolean;
    IsUpper : Boolean;
    ERC1 : Double;
    ERCInf : Double;
    Q50 : TReal1DArray;
    Q90 : TReal1DArray;
    V : Double;
begin
    Err50 := False;
    Err90 := False;
    ErrLess := False;
    ErrSpec := False;
    SetLength(Q50, 2);
    SetLength(Q90, 2);
    N:=1;
    while N<=MaxN do
    begin
        IsUpper := AP_FP_Greater(RandomReal,0.5);
        
        //
        // general test
        //
        SetLength(A, N, N);
        I:=0;
        while I<=1 do
        begin
            Q50[I] := 0;
            Q90[I] := 0;
            Inc(I);
        end;
        Pass:=1;
        while Pass<=PassCount do
        begin
            SPDMatrixRndCond(N, Exp(RandomReal*Ln(1000)), A);
            RMatrixRefRCond(A, N, ERC1, ERCInf);
            RMatrixDropHalf(A, N, IsUpper);
            RMatrixMakeACopy(A, N, N, CHA);
            SPDMatrixCholesky(CHA, N, IsUpper);
            
            //
            // normal
            //
            V := 1/SPDMatrixRCond(A, N, IsUpper);
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
            V := 1/SPDMatrixCholeskyRCond(CHA, N, IsUpper);
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
        I:=0;
        while I<=1 do
        begin
            Err50 := Err50 or AP_FP_Less(Q50[I],0.50);
            Err90 := Err90 or AP_FP_Less(Q90[I],0.90);
            Inc(I);
        end;
        
        //
        // degenerate matrix test
        //
        if N>=3 then
        begin
            SetLength(A, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := 0.0;
                    Inc(J);
                end;
                Inc(I);
            end;
            A[0,0] := 1;
            A[N-1,N-1] := 1;
            ErrSpec := ErrSpec or AP_FP_Neq(SPDMatrixRCond(A, N, IsUpper),-1);
            ErrSpec := ErrSpec or AP_FP_Neq(SPDMatrixCholeskyRCond(A, N, IsUpper),0);
        end;
        
        //
        // near-degenerate matrix test
        //
        if N>=2 then
        begin
            SetLength(A, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := 0.0;
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=N-1 do
            begin
                A[I,I] := 1;
                Inc(I);
            end;
            I := RandomInteger(N);
            A[I,I] := 0.1*MaxRealNumber;
            ErrSpec := ErrSpec or AP_FP_Neq(SPDMatrixRCond(A, N, IsUpper),0);
            ErrSpec := ErrSpec or AP_FP_Neq(SPDMatrixCholeskyRCond(A, N, IsUpper),0);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    Result :=  not (Err50 or Err90 or ErrLess or ErrSpec);
end;


(*************************************************************************
Returns True for successful test, False - for failed test
*************************************************************************)
function TestCMatrixRCond(MaxN : AlglibInteger;
     PassCount : AlglibInteger):Boolean;
var
    A : TComplex2DArray;
    LUA : TComplex2DArray;
    P : TInteger1DArray;
    N : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    Pass : AlglibInteger;
    Err50 : Boolean;
    Err90 : Boolean;
    ErrLess : Boolean;
    ErrSpec : Boolean;
    ERC1 : Double;
    ERCInf : Double;
    Q50 : TReal1DArray;
    Q90 : TReal1DArray;
    V : Double;
begin
    SetLength(Q50, 3+1);
    SetLength(Q90, 3+1);
    Err50 := False;
    Err90 := False;
    ErrLess := False;
    ErrSpec := False;
    
    //
    // process
    //
    N:=1;
    while N<=MaxN do
    begin
        
        //
        // special test for zero matrix
        //
        CMatrixGenZero(A, N);
        CMatrixMakeACopy(A, N, N, LUA);
        CMatrixLU(LUA, N, N, P);
        ErrSpec := ErrSpec or AP_FP_Neq(CMatrixRCond1(A, N),0);
        ErrSpec := ErrSpec or AP_FP_Neq(CMatrixRCondInf(A, N),0);
        ErrSpec := ErrSpec or AP_FP_Neq(CMatrixLURCond1(LUA, N),0);
        ErrSpec := ErrSpec or AP_FP_Neq(CMatrixLURCondInf(LUA, N),0);
        
        //
        // general test
        //
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
            CMatrixRndCond(N, Exp(RandomReal*Ln(1000)), A);
            CMatrixMakeACopy(A, N, N, LUA);
            CMatrixLU(LUA, N, N, P);
            CMatrixRefRCond(A, N, ERC1, ERCInf);
            
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
        
        //
        // degenerate matrix test
        //
        if N>=3 then
        begin
            SetLength(A, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := C_Complex(0.0);
                    Inc(J);
                end;
                Inc(I);
            end;
            A[0,0] := C_Complex(1);
            A[N-1,N-1] := C_Complex(1);
            ErrSpec := ErrSpec or AP_FP_Neq(CMatrixRCond1(A, N),0);
            ErrSpec := ErrSpec or AP_FP_Neq(cMatrixRCondInf(A, N),0);
            ErrSpec := ErrSpec or AP_FP_Neq(CMatrixLURCond1(A, N),0);
            ErrSpec := ErrSpec or AP_FP_Neq(CMatrixLURCondInf(A, N),0);
        end;
        
        //
        // near-degenerate matrix test
        //
        if N>=2 then
        begin
            SetLength(A, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := C_Complex(0.0);
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=N-1 do
            begin
                A[I,I] := C_Complex(1);
                Inc(I);
            end;
            I := RandomInteger(N);
            A[I,I] := C_Complex(0.1*MaxRealNumber);
            ErrSpec := ErrSpec or AP_FP_Neq(CMatrixRCond1(A, N),0);
            ErrSpec := ErrSpec or AP_FP_Neq(cMatrixRCondInf(A, N),0);
            ErrSpec := ErrSpec or AP_FP_Neq(CMatrixLURCond1(A, N),0);
            ErrSpec := ErrSpec or AP_FP_Neq(CMatrixLURCondInf(A, N),0);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    Result :=  not (Err50 or Err90 or ErrLess or ErrSpec);
end;


(*************************************************************************
Returns True for successful test, False - for failed test
*************************************************************************)
function TestHPDMatrixRCond(MaxN : AlglibInteger;
     PassCount : AlglibInteger):Boolean;
var
    A : TComplex2DArray;
    CHA : TComplex2DArray;
    P : TInteger1DArray;
    N : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    Pass : AlglibInteger;
    Err50 : Boolean;
    Err90 : Boolean;
    ErrSpec : Boolean;
    ErrLess : Boolean;
    IsUpper : Boolean;
    ERC1 : Double;
    ERCInf : Double;
    Q50 : TReal1DArray;
    Q90 : TReal1DArray;
    V : Double;
begin
    Err50 := False;
    Err90 := False;
    ErrLess := False;
    ErrSpec := False;
    SetLength(Q50, 2);
    SetLength(Q90, 2);
    N:=1;
    while N<=MaxN do
    begin
        IsUpper := AP_FP_Greater(RandomReal,0.5);
        
        //
        // general test
        //
        SetLength(A, N, N);
        I:=0;
        while I<=1 do
        begin
            Q50[I] := 0;
            Q90[I] := 0;
            Inc(I);
        end;
        Pass:=1;
        while Pass<=PassCount do
        begin
            HPDMatrixRndCond(N, Exp(RandomReal*Ln(1000)), A);
            CMatrixRefRCond(A, N, ERC1, ERCInf);
            CMatrixDropHalf(A, N, IsUpper);
            CMatrixMakeACopy(A, N, N, CHA);
            HPDMatrixCholesky(CHA, N, IsUpper);
            
            //
            // normal
            //
            V := 1/HPDMatrixRCond(A, N, IsUpper);
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
            V := 1/HPDMatrixCholeskyRCond(CHA, N, IsUpper);
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
        I:=0;
        while I<=1 do
        begin
            Err50 := Err50 or AP_FP_Less(Q50[I],0.50);
            Err90 := Err90 or AP_FP_Less(Q90[I],0.90);
            Inc(I);
        end;
        
        //
        // degenerate matrix test
        //
        if N>=3 then
        begin
            SetLength(A, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := C_Complex(0.0);
                    Inc(J);
                end;
                Inc(I);
            end;
            A[0,0] := C_Complex(1);
            A[N-1,N-1] := C_Complex(1);
            ErrSpec := ErrSpec or AP_FP_Neq(HPDMatrixRCond(A, N, IsUpper),-1);
            ErrSpec := ErrSpec or AP_FP_Neq(HPDMatrixCholeskyRCond(A, N, IsUpper),0);
        end;
        
        //
        // near-degenerate matrix test
        //
        if N>=2 then
        begin
            SetLength(A, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := C_Complex(0.0);
                    Inc(J);
                end;
                Inc(I);
            end;
            I:=0;
            while I<=N-1 do
            begin
                A[I,I] := C_Complex(1);
                Inc(I);
            end;
            I := RandomInteger(N);
            A[I,I] := C_Complex(0.1*MaxRealNumber);
            ErrSpec := ErrSpec or AP_FP_Neq(HPDMatrixRCond(A, N, IsUpper),0);
            ErrSpec := ErrSpec or AP_FP_Neq(HPDMatrixCholeskyRCond(A, N, IsUpper),0);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    Result :=  not (Err50 or Err90 or ErrLess or ErrSpec);
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