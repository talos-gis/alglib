unit testtrfacunit;
interface
uses Math, Sysutils, Ap, reflections, creflections, hqrnd, matgen, ablasf, ablas, trfac;

function TestTRFAC(Silent : Boolean):Boolean;
function testtrfacunit_test_silent():Boolean;
function testtrfacunit_test():Boolean;

implementation

procedure TestCLUProblem(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Threshold : Double;
     var Err : Boolean;
     var PropErr : Boolean);forward;
procedure TestRLUProblem(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Threshold : Double;
     var Err : Boolean;
     var PropErr : Boolean);forward;


function TestTRFAC(Silent : Boolean):Boolean;
var
    RA : TReal2DArray;
    RAL : TReal2DArray;
    RAU : TReal2DArray;
    CA : TComplex2DArray;
    CAL : TComplex2DArray;
    CAU : TComplex2DArray;
    M : AlglibInteger;
    N : AlglibInteger;
    MX : AlglibInteger;
    MaxMN : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    MinIJ : AlglibInteger;
    Pass : AlglibInteger;
    VC : Complex;
    VR : Double;
    WasErrors : Boolean;
    SPDErr : Boolean;
    HPDErr : Boolean;
    RErr : Boolean;
    CErr : Boolean;
    PropErr : Boolean;
    Threshold : Double;
    i_ : AlglibInteger;
begin
    RErr := False;
    SPDErr := False;
    CErr := False;
    HPDErr := False;
    PropErr := False;
    WasErrors := False;
    MaxMN := 4*ABLASBlockSize(RA)+1;
    Threshold := 1000*MachineEpsilon*MaxMN;
    
    //
    // test LU
    //
    MX:=1;
    while MX<=MaxMN do
    begin
        
        //
        // Initialize N/M, both are <=MX,
        // at least one of them is exactly equal to MX
        //
        N := 1+RandomInteger(MX);
        M := 1+RandomInteger(MX);
        if AP_FP_Greater(RandomReal,0.5) then
        begin
            N := MX;
        end
        else
        begin
            M := MX;
        end;
        
        //
        // First, test on zero matrix
        //
        SetLength(RA, M, N);
        SetLength(CA, M, N);
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                RA[I,J] := 0;
                CA[I,J] := C_Complex(0);
                Inc(J);
            end;
            Inc(I);
        end;
        TestCLUProblem(CA, M, N, Threshold, CErr, PropErr);
        TestRLUProblem(RA, M, N, Threshold, RErr, PropErr);
        
        //
        // Second, random matrix with moderate condition number
        //
        SetLength(RA, M, N);
        SetLength(CA, M, N);
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                RA[I,J] := 0;
                CA[I,J] := C_Complex(0);
                Inc(J);
            end;
            Inc(I);
        end;
        I:=0;
        while I<=Min(M, N)-1 do
        begin
            RA[I,I] := 1+10*RandomReal;
            CA[I,I] := C_Complex(1+10*RandomReal);
            Inc(I);
        end;
        CMatrixRndOrthogonalFromTheLeft(CA, M, N);
        CMatrixRndOrthogonalFromTheRight(CA, M, N);
        RMatrixRndOrthogonalFromTheLeft(RA, M, N);
        RMatrixRndOrthogonalFromTheRight(RA, M, N);
        TestCLUProblem(CA, M, N, Threshold, CErr, PropErr);
        TestRLUProblem(RA, M, N, Threshold, RErr, PropErr);
        Inc(MX);
    end;
    
    //
    // Test Cholesky
    //
    N:=1;
    while N<=MaxMN do
    begin
        
        //
        // Load CA (HPD matrix with low condition number),
        //      CAL and CAU - its lower and upper triangles
        //
        HPDMatrixRndCond(N, 1+50*RandomReal, CA);
        SetLength(CAL, N, N);
        SetLength(CAU, N, N);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                CAL[I,J] := C_Complex(I);
                CAU[I,J] := C_Complex(J);
                Inc(J);
            end;
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            for i_ := 0 to I do
            begin
                CAL[I,i_] := CA[I,i_];
            end;
            for i_ := I to N-1 do
            begin
                CAU[I,i_] := CA[I,i_];
            end;
            Inc(I);
        end;
        
        //
        // Test HPDMatrixCholesky:
        // 1. it must leave upper (lower) part unchanged
        // 2. max(A-L*L^H) must be small
        //
        if HPDMatrixCholesky(CAL, N, False) then
        begin
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    if J>I then
                    begin
                        HPDErr := HPDErr or C_NotEqualR(CAL[I,J],I);
                    end
                    else
                    begin
                        VC := C_Complex(0.0);
                        for i_ := 0 to J do
                        begin
                            VC := C_Add(VC,C_Mul(CAL[I,i_],Conj(CAL[J,i_])));
                        end;
                        HPDErr := HPDErr or AP_FP_Greater(AbsComplex(C_Sub(CA[I,J],VC)),Threshold);
                    end;
                    Inc(J);
                end;
                Inc(I);
            end;
        end
        else
        begin
            HPDErr := True;
        end;
        if HPDMatrixCholesky(CAU, N, True) then
        begin
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    if J<I then
                    begin
                        HPDErr := HPDErr or C_NotEqualR(CAU[I,J],J);
                    end
                    else
                    begin
                        VC := C_Complex(0.0);
                        for i_ := 0 to I do
                        begin
                            VC := C_Add(VC,C_Mul(Conj(CAU[i_,I]),CAU[i_,J]));
                        end;
                        HPDErr := HPDErr or AP_FP_Greater(AbsComplex(C_Sub(CA[I,J],VC)),Threshold);
                    end;
                    Inc(J);
                end;
                Inc(I);
            end;
        end
        else
        begin
            HPDErr := True;
        end;
        
        //
        // Load RA (SPD matrix with low condition number),
        //      RAL and RAU - its lower and upper triangles
        //
        SPDMatrixRndCond(N, 1+50*RandomReal, RA);
        SetLength(RAL, N, N);
        SetLength(RAU, N, N);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                RAL[I,J] := I;
                RAU[I,J] := J;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            APVMove(@RAL[I][0], 0, I, @RA[I][0], 0, I);
            APVMove(@RAU[I][0], I, N-1, @RA[I][0], I, N-1);
            Inc(I);
        end;
        
        //
        // Test SPDMatrixCholesky:
        // 1. it must leave upper (lower) part unchanged
        // 2. max(A-L*L^H) must be small
        //
        if SPDMatrixCholesky(RAL, N, False) then
        begin
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    if J>I then
                    begin
                        SPDErr := SPDErr or AP_FP_Neq(RAL[I,J],I);
                    end
                    else
                    begin
                        VR := APVDotProduct(@RAL[I][0], 0, J, @RAL[J][0], 0, J);
                        SPDErr := SPDErr or AP_FP_Greater(AbsReal(RA[I,J]-VR),Threshold);
                    end;
                    Inc(J);
                end;
                Inc(I);
            end;
        end
        else
        begin
            SPDErr := True;
        end;
        if SPDMatrixCholesky(RAU, N, True) then
        begin
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    if J<I then
                    begin
                        SPDErr := SPDErr or AP_FP_Neq(RAU[I,J],J);
                    end
                    else
                    begin
                        VR := 0.0;
                        for i_ := 0 to I do
                        begin
                            VR := VR + RAU[i_,I]*RAU[i_,J];
                        end;
                        SPDErr := SPDErr or AP_FP_Greater(AbsReal(RA[I,J]-VR),Threshold);
                    end;
                    Inc(J);
                end;
                Inc(I);
            end;
        end
        else
        begin
            SPDErr := True;
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    WasErrors := RErr or SPDErr or CErr or HPDErr or PropErr;
    if  not Silent then
    begin
        Write(Format('TESTING TRIANGULAR FACTORIZATIONS'#13#10'',[]));
        Write(Format('* REAL:                                  ',[]));
        if RErr then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* SPD:                                   ',[]));
        if SPDErr then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* COMPLEX:                               ',[]));
        if CErr then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* HPD:                                   ',[]));
        if HPDErr then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* OTHER PROPERTIES:                      ',[]));
        if PropErr then
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


procedure TestCLUProblem(const A : TComplex2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Threshold : Double;
     var Err : Boolean;
     var PropErr : Boolean);
var
    CA : TComplex2DArray;
    CL : TComplex2DArray;
    CU : TComplex2DArray;
    CA2 : TComplex2DArray;
    CT : TComplex1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    MinMN : AlglibInteger;
    V : Complex;
    P : TInteger1DArray;
    i_ : AlglibInteger;
begin
    MinMN := Min(M, N);
    
    //
    // PLU test
    //
    SetLength(CA, M, N);
    I:=0;
    while I<=M-1 do
    begin
        for i_ := 0 to N-1 do
        begin
            CA[I,i_] := A[I,i_];
        end;
        Inc(I);
    end;
    CMatrixPLU(CA, M, N, P);
    I:=0;
    while I<=MinMN-1 do
    begin
        if (P[I]<I) or (P[I]>=M) then
        begin
            PropErr := False;
            Exit;
        end;
        Inc(I);
    end;
    SetLength(CL, M, MinMN);
    J:=0;
    while J<=MinMN-1 do
    begin
        I:=0;
        while I<=J-1 do
        begin
            CL[I,J] := C_Complex(0.0);
            Inc(I);
        end;
        CL[J,J] := C_Complex(1.0);
        I:=J+1;
        while I<=M-1 do
        begin
            CL[I,J] := CA[I,J];
            Inc(I);
        end;
        Inc(J);
    end;
    SetLength(CU, MinMN, N);
    I:=0;
    while I<=MinMN-1 do
    begin
        J:=0;
        while J<=I-1 do
        begin
            CU[I,J] := C_Complex(0.0);
            Inc(J);
        end;
        J:=I;
        while J<=N-1 do
        begin
            CU[I,J] := CA[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
    SetLength(CA2, M, N);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := C_Complex(0.0);
            for i_ := 0 to MinMN-1 do
            begin
                V := C_Add(V,C_Mul(CL[I,i_],CU[i_,J]));
            end;
            CA2[I,J] := V;
            Inc(J);
        end;
        Inc(I);
    end;
    SetLength(CT, N);
    I:=MinMN-1;
    while I>=0 do
    begin
        if I<>P[I] then
        begin
            for i_ := 0 to N-1 do
            begin
                CT[i_] := CA2[I,i_];
            end;
            for i_ := 0 to N-1 do
            begin
                CA2[I,i_] := CA2[P[I],i_];
            end;
            for i_ := 0 to N-1 do
            begin
                CA2[P[I],i_] := CT[i_];
            end;
        end;
        Dec(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            Err := Err or AP_FP_Greater(AbsComplex(C_Sub(A[I,J],CA2[I,J])),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // LUP test
    //
    SetLength(CA, M, N);
    I:=0;
    while I<=M-1 do
    begin
        for i_ := 0 to N-1 do
        begin
            CA[I,i_] := A[I,i_];
        end;
        Inc(I);
    end;
    CMatrixLUP(CA, M, N, P);
    I:=0;
    while I<=MinMN-1 do
    begin
        if (P[I]<I) or (P[I]>=N) then
        begin
            PropErr := False;
            Exit;
        end;
        Inc(I);
    end;
    SetLength(CL, M, MinMN);
    J:=0;
    while J<=MinMN-1 do
    begin
        I:=0;
        while I<=J-1 do
        begin
            CL[I,J] := C_Complex(0.0);
            Inc(I);
        end;
        I:=J;
        while I<=M-1 do
        begin
            CL[I,J] := CA[I,J];
            Inc(I);
        end;
        Inc(J);
    end;
    SetLength(CU, MinMN, N);
    I:=0;
    while I<=MinMN-1 do
    begin
        J:=0;
        while J<=I-1 do
        begin
            CU[I,J] := C_Complex(0.0);
            Inc(J);
        end;
        CU[I,I] := C_Complex(1.0);
        J:=I+1;
        while J<=N-1 do
        begin
            CU[I,J] := CA[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
    SetLength(CA2, M, N);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := C_Complex(0.0);
            for i_ := 0 to MinMN-1 do
            begin
                V := C_Add(V,C_Mul(CL[I,i_],CU[i_,J]));
            end;
            CA2[I,J] := V;
            Inc(J);
        end;
        Inc(I);
    end;
    SetLength(CT, M);
    I:=MinMN-1;
    while I>=0 do
    begin
        if I<>P[I] then
        begin
            for i_ := 0 to M-1 do
            begin
                CT[i_] := CA2[i_,I];
            end;
            for i_ := 0 to M-1 do
            begin
                CA2[i_,I] := CA2[i_,P[I]];
            end;
            for i_ := 0 to M-1 do
            begin
                CA2[i_,P[I]] := CT[i_];
            end;
        end;
        Dec(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            Err := Err or AP_FP_Greater(AbsComplex(C_Sub(A[I,J],CA2[I,J])),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
end;


procedure TestRLUProblem(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     Threshold : Double;
     var Err : Boolean;
     var PropErr : Boolean);
var
    CA : TReal2DArray;
    CL : TReal2DArray;
    CU : TReal2DArray;
    CA2 : TReal2DArray;
    CT : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    MinMN : AlglibInteger;
    V : Double;
    P : TInteger1DArray;
    i_ : AlglibInteger;
begin
    MinMN := Min(M, N);
    
    //
    // PLU test
    //
    SetLength(CA, M, N);
    I:=0;
    while I<=M-1 do
    begin
        APVMove(@CA[I][0], 0, N-1, @A[I][0], 0, N-1);
        Inc(I);
    end;
    RMatrixPLU(CA, M, N, P);
    I:=0;
    while I<=MinMN-1 do
    begin
        if (P[I]<I) or (P[I]>=M) then
        begin
            PropErr := False;
            Exit;
        end;
        Inc(I);
    end;
    SetLength(CL, M, MinMN);
    J:=0;
    while J<=MinMN-1 do
    begin
        I:=0;
        while I<=J-1 do
        begin
            CL[I,J] := 0.0;
            Inc(I);
        end;
        CL[J,J] := 1.0;
        I:=J+1;
        while I<=M-1 do
        begin
            CL[I,J] := CA[I,J];
            Inc(I);
        end;
        Inc(J);
    end;
    SetLength(CU, MinMN, N);
    I:=0;
    while I<=MinMN-1 do
    begin
        J:=0;
        while J<=I-1 do
        begin
            CU[I,J] := 0.0;
            Inc(J);
        end;
        J:=I;
        while J<=N-1 do
        begin
            CU[I,J] := CA[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
    SetLength(CA2, M, N);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := 0.0;
            for i_ := 0 to MinMN-1 do
            begin
                V := V + CL[I,i_]*CU[i_,J];
            end;
            CA2[I,J] := V;
            Inc(J);
        end;
        Inc(I);
    end;
    SetLength(CT, N);
    I:=MinMN-1;
    while I>=0 do
    begin
        if I<>P[I] then
        begin
            APVMove(@CT[0], 0, N-1, @CA2[I][0], 0, N-1);
            APVMove(@CA2[I][0], 0, N-1, @CA2[P[I]][0], 0, N-1);
            APVMove(@CA2[P[I]][0], 0, N-1, @CT[0], 0, N-1);
        end;
        Dec(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            Err := Err or AP_FP_Greater(AbsReal(A[I,J]-CA2[I,J]),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // LUP test
    //
    SetLength(CA, M, N);
    I:=0;
    while I<=M-1 do
    begin
        APVMove(@CA[I][0], 0, N-1, @A[I][0], 0, N-1);
        Inc(I);
    end;
    RMatrixLUP(CA, M, N, P);
    I:=0;
    while I<=MinMN-1 do
    begin
        if (P[I]<I) or (P[I]>=N) then
        begin
            PropErr := False;
            Exit;
        end;
        Inc(I);
    end;
    SetLength(CL, M, MinMN);
    J:=0;
    while J<=MinMN-1 do
    begin
        I:=0;
        while I<=J-1 do
        begin
            CL[I,J] := 0.0;
            Inc(I);
        end;
        I:=J;
        while I<=M-1 do
        begin
            CL[I,J] := CA[I,J];
            Inc(I);
        end;
        Inc(J);
    end;
    SetLength(CU, MinMN, N);
    I:=0;
    while I<=MinMN-1 do
    begin
        J:=0;
        while J<=I-1 do
        begin
            CU[I,J] := 0.0;
            Inc(J);
        end;
        CU[I,I] := 1.0;
        J:=I+1;
        while J<=N-1 do
        begin
            CU[I,J] := CA[I,J];
            Inc(J);
        end;
        Inc(I);
    end;
    SetLength(CA2, M, N);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            V := 0.0;
            for i_ := 0 to MinMN-1 do
            begin
                V := V + CL[I,i_]*CU[i_,J];
            end;
            CA2[I,J] := V;
            Inc(J);
        end;
        Inc(I);
    end;
    SetLength(CT, M);
    I:=MinMN-1;
    while I>=0 do
    begin
        if I<>P[I] then
        begin
            for i_ := 0 to M-1 do
            begin
                CT[i_] := CA2[i_,I];
            end;
            for i_ := 0 to M-1 do
            begin
                CA2[i_,I] := CA2[i_,P[I]];
            end;
            for i_ := 0 to M-1 do
            begin
                CA2[i_,P[I]] := CT[i_];
            end;
        end;
        Dec(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            Err := Err or AP_FP_Greater(AbsReal(A[I,J]-CA2[I,J]),Threshold);
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testtrfacunit_test_silent():Boolean;
begin
    Result := TestTRFAC(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testtrfacunit_test():Boolean;
begin
    Result := TestTRFAC(False);
end;


end.