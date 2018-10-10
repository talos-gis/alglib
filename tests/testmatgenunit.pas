unit testmatgenunit;
interface
uses Math, Sysutils, Ap, reflections, creflections, hqrnd, matgen;

function TestMatGen(Silent : Boolean):Boolean;
function IsSPD(A : TReal2DArray; N : AlglibInteger; IsUpper : Boolean):Boolean;
function ObsoleteSVDDecomposition(var a : TReal2DArray;
     m : AlglibInteger;
     n : AlglibInteger;
     var w : TReal1DArray;
     var v : TReal2DArray):Boolean;
function testmatgenunit_test_silent():Boolean;
function testmatgenunit_test():Boolean;

implementation

const
    MaxSVDIterations = 60;

procedure Unset2D(var A : TReal2DArray);forward;
procedure Unset2DC(var A : TComplex2DArray);forward;
function IsHPD(A : TComplex2DArray; N : AlglibInteger):Boolean;forward;
function SVDCond(const a : TReal2DArray; N : AlglibInteger):Double;forward;
function ExtSign(a : Double; b : Double):Double;forward;
function MyMax(a : Double; b : Double):Double;forward;
function Pythag(A : Double; B : Double):Double;forward;


function TestMatGen(Silent : Boolean):Boolean;
var
    A : TReal2DArray;
    B : TReal2DArray;
    U : TReal2DArray;
    V : TReal2DArray;
    CA : TComplex2DArray;
    CB : TComplex2DArray;
    R1 : TReal2DArray;
    R2 : TReal2DArray;
    C1 : TComplex2DArray;
    C2 : TComplex2DArray;
    W : TReal1DArray;
    N : AlglibInteger;
    MaxN : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    WasErrors : Boolean;
    Cond : Double;
    Threshold : Double;
    VT : Double;
    CT : Complex;
    MinW : Double;
    MaxW : Double;
    SErr : Boolean;
    HErr : Boolean;
    SPDErr : Boolean;
    HPDErr : Boolean;
    RErr : Boolean;
    CErr : Boolean;
    i_ : AlglibInteger;
begin
    RErr := False;
    CErr := False;
    SErr := False;
    HErr := False;
    SPDErr := False;
    HPDErr := False;
    WasErrors := False;
    MaxN := 20;
    PassCount := 15;
    Threshold := 1000*MachineEpsilon;
    
    //
    // Testing orthogonal
    //
    N:=1;
    while N<=MaxN do
    begin
        Pass:=1;
        while Pass<=PassCount do
        begin
            SetLength(R1, N-1+1, 2*N-1+1);
            SetLength(R2, 2*N-1+1, N-1+1);
            SetLength(C1, N-1+1, 2*N-1+1);
            SetLength(C2, 2*N-1+1, N-1+1);
            
            //
            // Random orthogonal, real
            //
            Unset2D(A);
            Unset2D(B);
            RMatrixRndOrthogonal(N, A);
            RMatrixRndOrthogonal(N, B);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    
                    //
                    // orthogonality test
                    //
                    VT := APVDotProduct(@A[I][0], 0, N-1, @A[J][0], 0, N-1);
                    if I=J then
                    begin
                        RErr := RErr or AP_FP_Greater(AbsReal(VT-1),Threshold);
                    end
                    else
                    begin
                        RErr := RErr or AP_FP_Greater(AbsReal(VT),Threshold);
                    end;
                    VT := APVDotProduct(@B[I][0], 0, N-1, @B[J][0], 0, N-1);
                    if I=J then
                    begin
                        RErr := RErr or AP_FP_Greater(AbsReal(VT-1),Threshold);
                    end
                    else
                    begin
                        RErr := RErr or AP_FP_Greater(AbsReal(VT),Threshold);
                    end;
                    
                    //
                    // test for difference in A and B
                    //
                    if N>=2 then
                    begin
                        RErr := RErr or AP_FP_Eq(A[I,J],B[I,J]);
                    end;
                    Inc(J);
                end;
                Inc(I);
            end;
            
            //
            // Random orthogonal, complex
            //
            Unset2DC(CA);
            Unset2DC(CB);
            CMatrixRndOrthogonal(N, CA);
            CMatrixRndOrthogonal(N, CB);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    
                    //
                    // orthogonality test
                    //
                    CT := C_Complex(0.0);
                    for i_ := 0 to N-1 do
                    begin
                        CT := C_Add(CT,C_Mul(CA[I,i_],Conj(CA[J,i_])));
                    end;
                    if I=J then
                    begin
                        CErr := CErr or AP_FP_Greater(AbsComplex(C_SubR(CT,1)),Threshold);
                    end
                    else
                    begin
                        CErr := CErr or AP_FP_Greater(AbsComplex(CT),Threshold);
                    end;
                    CT := C_Complex(0.0);
                    for i_ := 0 to N-1 do
                    begin
                        CT := C_Add(CT,C_Mul(CB[I,i_],Conj(CB[J,i_])));
                    end;
                    if I=J then
                    begin
                        CErr := CErr or AP_FP_Greater(AbsComplex(C_SubR(CT,1)),Threshold);
                    end
                    else
                    begin
                        CErr := CErr or AP_FP_Greater(AbsComplex(CT),Threshold);
                    end;
                    
                    //
                    // test for difference in A and B
                    //
                    if N>=2 then
                    begin
                        CErr := CErr or C_Equal(CA[I,J],CB[I,J]);
                    end;
                    Inc(J);
                end;
                Inc(I);
            end;
            
            //
            // From the right real tests:
            // 1. E*Q is orthogonal
            // 2. Q1<>Q2 (routine result is changing)
            // 3. (E E)'*Q = (Q' Q')' (correct handling of non-square matrices)
            //
            Unset2D(A);
            Unset2D(B);
            SetLength(A, N-1+1, N-1+1);
            SetLength(B, N-1+1, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := 0;
                    B[I,J] := 0;
                    Inc(J);
                end;
                A[I,I] := 1;
                B[I,I] := 1;
                Inc(I);
            end;
            RMatrixRndOrthogonalFromTheRight(A, N, N);
            RMatrixRndOrthogonalFromTheRight(B, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    
                    //
                    // orthogonality test
                    //
                    VT := APVDotProduct(@A[I][0], 0, N-1, @A[J][0], 0, N-1);
                    if I=J then
                    begin
                        RErr := RErr or AP_FP_Greater(AbsReal(VT-1),Threshold);
                    end
                    else
                    begin
                        RErr := RErr or AP_FP_Greater(AbsReal(VT),Threshold);
                    end;
                    VT := APVDotProduct(@B[I][0], 0, N-1, @B[J][0], 0, N-1);
                    if I=J then
                    begin
                        RErr := RErr or AP_FP_Greater(AbsReal(VT-1),Threshold);
                    end
                    else
                    begin
                        RErr := RErr or AP_FP_Greater(AbsReal(VT),Threshold);
                    end;
                    
                    //
                    // test for difference in A and B
                    //
                    if N>=2 then
                    begin
                        RErr := RErr or AP_FP_Eq(A[I,J],B[I,J]);
                    end;
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
                    R2[I,J] := 2*RandomReal-1;
                    R2[I+N,J] := R2[I,J];
                    Inc(J);
                end;
                Inc(I);
            end;
            RMatrixRndOrthogonalFromTheRight(R2, 2*N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    RErr := RErr or AP_FP_Greater(AbsReal(R2[I+N,J]-R2[I,J]),Threshold);
                    Inc(J);
                end;
                Inc(I);
            end;
            
            //
            // From the left real tests:
            // 1. Q*E is orthogonal
            // 2. Q1<>Q2 (routine result is changing)
            // 3. Q*(E E) = (Q Q) (correct handling of non-square matrices)
            //
            Unset2D(A);
            Unset2D(B);
            SetLength(A, N-1+1, N-1+1);
            SetLength(B, N-1+1, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    A[I,J] := 0;
                    B[I,J] := 0;
                    Inc(J);
                end;
                A[I,I] := 1;
                B[I,I] := 1;
                Inc(I);
            end;
            RMatrixRndOrthogonalFromTheLeft(A, N, N);
            RMatrixRndOrthogonalFromTheLeft(B, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    
                    //
                    // orthogonality test
                    //
                    VT := APVDotProduct(@A[I][0], 0, N-1, @A[J][0], 0, N-1);
                    if I=J then
                    begin
                        RErr := RErr or AP_FP_Greater(AbsReal(VT-1),Threshold);
                    end
                    else
                    begin
                        RErr := RErr or AP_FP_Greater(AbsReal(VT),Threshold);
                    end;
                    VT := APVDotProduct(@B[I][0], 0, N-1, @B[J][0], 0, N-1);
                    if I=J then
                    begin
                        RErr := RErr or AP_FP_Greater(AbsReal(VT-1),Threshold);
                    end
                    else
                    begin
                        RErr := RErr or AP_FP_Greater(AbsReal(VT),Threshold);
                    end;
                    
                    //
                    // test for difference in A and B
                    //
                    if N>=2 then
                    begin
                        RErr := RErr or AP_FP_Eq(A[I,J],B[I,J]);
                    end;
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
                    R1[I,J] := 2*RandomReal-1;
                    R1[I,J+N] := R1[I,J];
                    Inc(J);
                end;
                Inc(I);
            end;
            RMatrixRndOrthogonalFromTheLeft(R1, N, 2*N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    RErr := RErr or AP_FP_Greater(AbsReal(R1[I,J]-R1[I,J+N]),Threshold);
                    Inc(J);
                end;
                Inc(I);
            end;
            
            //
            // From the right complex tests:
            // 1. E*Q is orthogonal
            // 2. Q1<>Q2 (routine result is changing)
            // 3. (E E)'*Q = (Q' Q')' (correct handling of non-square matrices)
            //
            Unset2DC(CA);
            Unset2DC(CB);
            SetLength(CA, N-1+1, N-1+1);
            SetLength(CB, N-1+1, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    CA[I,J] := C_Complex(0);
                    CB[I,J] := C_Complex(0);
                    Inc(J);
                end;
                CA[I,I] := C_Complex(1);
                CB[I,I] := C_Complex(1);
                Inc(I);
            end;
            CMatrixRndOrthogonalFromTheRight(CA, N, N);
            CMatrixRndOrthogonalFromTheRight(CB, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    
                    //
                    // orthogonality test
                    //
                    CT := C_Complex(0.0);
                    for i_ := 0 to N-1 do
                    begin
                        CT := C_Add(CT,C_Mul(CA[I,i_],Conj(CA[J,i_])));
                    end;
                    if I=J then
                    begin
                        CErr := CErr or AP_FP_Greater(AbsComplex(C_SubR(CT,1)),Threshold);
                    end
                    else
                    begin
                        CErr := CErr or AP_FP_Greater(AbsComplex(CT),Threshold);
                    end;
                    CT := C_Complex(0.0);
                    for i_ := 0 to N-1 do
                    begin
                        CT := C_Add(CT,C_Mul(CB[I,i_],Conj(CB[J,i_])));
                    end;
                    if I=J then
                    begin
                        CErr := CErr or AP_FP_Greater(AbsComplex(C_SubR(CT,1)),Threshold);
                    end
                    else
                    begin
                        CErr := CErr or AP_FP_Greater(AbsComplex(CT),Threshold);
                    end;
                    
                    //
                    // test for difference in A and B
                    //
                    CErr := CErr or C_Equal(CA[I,J],CB[I,J]);
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
                    C2[I,J] := C_Complex(2*RandomReal-1);
                    C2[I+N,J] := C2[I,J];
                    Inc(J);
                end;
                Inc(I);
            end;
            CMatrixRndOrthogonalFromTheRight(C2, 2*N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    CErr := CErr or AP_FP_Greater(AbsComplex(C_Sub(C2[I+N,J],C2[I,J])),Threshold);
                    Inc(J);
                end;
                Inc(I);
            end;
            
            //
            // From the left complex tests:
            // 1. Q*E is orthogonal
            // 2. Q1<>Q2 (routine result is changing)
            // 3. Q*(E E) = (Q Q) (correct handling of non-square matrices)
            //
            Unset2DC(CA);
            Unset2DC(CB);
            SetLength(CA, N-1+1, N-1+1);
            SetLength(CB, N-1+1, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    CA[I,J] := C_Complex(0);
                    CB[I,J] := C_Complex(0);
                    Inc(J);
                end;
                CA[I,I] := C_Complex(1);
                CB[I,I] := C_Complex(1);
                Inc(I);
            end;
            CMatrixRndOrthogonalFromTheLeft(CA, N, N);
            CMatrixRndOrthogonalFromTheLeft(CB, N, N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    
                    //
                    // orthogonality test
                    //
                    CT := C_Complex(0.0);
                    for i_ := 0 to N-1 do
                    begin
                        CT := C_Add(CT,C_Mul(CA[I,i_],Conj(CA[J,i_])));
                    end;
                    if I=J then
                    begin
                        CErr := CErr or AP_FP_Greater(AbsComplex(C_SubR(CT,1)),Threshold);
                    end
                    else
                    begin
                        CErr := CErr or AP_FP_Greater(AbsComplex(CT),Threshold);
                    end;
                    CT := C_Complex(0.0);
                    for i_ := 0 to N-1 do
                    begin
                        CT := C_Add(CT,C_Mul(CB[I,i_],Conj(CB[J,i_])));
                    end;
                    if I=J then
                    begin
                        CErr := CErr or AP_FP_Greater(AbsComplex(C_SubR(CT,1)),Threshold);
                    end
                    else
                    begin
                        CErr := CErr or AP_FP_Greater(AbsComplex(CT),Threshold);
                    end;
                    
                    //
                    // test for difference in A and B
                    //
                    CErr := CErr or C_Equal(CA[I,J],CB[I,J]);
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
                    C1[I,J] := C_Complex(2*RandomReal-1);
                    C1[I,J+N] := C1[I,J];
                    Inc(J);
                end;
                Inc(I);
            end;
            CMatrixRndOrthogonalFromTheLeft(C1, N, 2*N);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    CErr := CErr or AP_FP_Greater(AbsComplex(C_Sub(C1[I,J],C1[I,J+N])),Threshold);
                    Inc(J);
                end;
                Inc(I);
            end;
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // Testing GCond
    //
    N:=2;
    while N<=MaxN do
    begin
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // real test
            //
            Unset2D(A);
            Cond := Exp(Ln(1000)*RandomReal);
            RMatrixRndCond(N, Cond, A);
            SetLength(B, N+1, N+1);
            I:=1;
            while I<=N do
            begin
                J:=1;
                while J<=N do
                begin
                    B[I,J] := A[I-1,J-1];
                    Inc(J);
                end;
                Inc(I);
            end;
            if ObsoleteSVDDecomposition(B, N, N, W, V) then
            begin
                MaxW := W[1];
                MinW := W[1];
                I:=2;
                while I<=N do
                begin
                    if AP_FP_Greater(W[I],MaxW) then
                    begin
                        MaxW := W[I];
                    end;
                    if AP_FP_Less(W[I],MinW) then
                    begin
                        MinW := W[I];
                    end;
                    Inc(I);
                end;
                VT := MaxW/MinW/Cond;
                if AP_FP_Greater(AbsReal(Ln(VT)),Ln(1+Threshold)) then
                begin
                    RErr := True;
                end;
            end;
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // Symmetric/SPD
    // N = 2 .. 30
    //
    N:=2;
    while N<=MaxN do
    begin
        
        //
        // SPD matrices
        //
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // Generate A
            //
            Unset2D(A);
            Cond := Exp(Ln(1000)*RandomReal);
            SPDMatrixRndCond(N, Cond, A);
            
            //
            // test condition number
            //
            SPDErr := SPDErr or AP_FP_Greater(SVDCond(A, N)/Cond-1,Threshold);
            
            //
            // test SPD
            //
            SPDErr := SPDErr or  not IsSPD(A, N, True);
            
            //
            // test that A is symmetic
            //
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    SPDErr := SPDErr or AP_FP_Greater(AbsReal(A[I,J]-A[J,I]),Threshold);
                    Inc(J);
                end;
                Inc(I);
            end;
            
            //
            // test for difference between A and B (subsequent matrix)
            //
            Unset2D(B);
            SPDMatrixRndCond(N, Cond, B);
            if N>=2 then
            begin
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        SPDErr := SPDErr or AP_FP_Eq(A[I,J],B[I,J]);
                        Inc(J);
                    end;
                    Inc(I);
                end;
            end;
            Inc(Pass);
        end;
        
        //
        // HPD matrices
        //
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // Generate A
            //
            Unset2DC(CA);
            Cond := Exp(Ln(1000)*RandomReal);
            HPDMatrixRndCond(N, Cond, CA);
            
            //
            // test HPD
            //
            HPDErr := HPDErr or  not IsHPD(CA, N);
            
            //
            // test that A is Hermitian
            //
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    HPDErr := HPDErr or AP_FP_Greater(AbsComplex(C_Sub(CA[I,J],Conj(CA[J,I]))),Threshold);
                    Inc(J);
                end;
                Inc(I);
            end;
            
            //
            // test for difference between A and B (subsequent matrix)
            //
            Unset2DC(CB);
            HPDMatrixRndCond(N, Cond, CB);
            if N>=2 then
            begin
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        HPDErr := HPDErr or C_Equal(CA[I,J],CB[I,J]);
                        Inc(J);
                    end;
                    Inc(I);
                end;
            end;
            Inc(Pass);
        end;
        
        //
        // Symmetric matrices
        //
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // test condition number
            //
            Unset2D(A);
            Cond := Exp(Ln(1000)*RandomReal);
            SMatrixRndCond(N, Cond, A);
            SErr := SErr or AP_FP_Greater(SVDCond(A, N)/Cond-1,Threshold);
            
            //
            // test for difference between A and B
            //
            Unset2D(B);
            SMatrixRndCond(N, Cond, B);
            if N>=2 then
            begin
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        SErr := SErr or AP_FP_Eq(A[I,J],B[I,J]);
                        Inc(J);
                    end;
                    Inc(I);
                end;
            end;
            Inc(Pass);
        end;
        
        //
        // Hermitian matrices
        //
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // Generate A
            //
            Unset2DC(CA);
            Cond := Exp(Ln(1000)*RandomReal);
            HMatrixRndCond(N, Cond, CA);
            
            //
            // test that A is Hermitian
            //
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    HErr := HErr or AP_FP_Greater(AbsComplex(C_Sub(CA[I,J],Conj(CA[J,I]))),Threshold);
                    Inc(J);
                end;
                Inc(I);
            end;
            
            //
            // test for difference between A and B (subsequent matrix)
            //
            Unset2DC(CB);
            HMatrixRndCond(N, Cond, CB);
            if N>=2 then
            begin
                I:=0;
                while I<=N-1 do
                begin
                    J:=0;
                    while J<=N-1 do
                    begin
                        HErr := HErr or C_Equal(CA[I,J],CB[I,J]);
                        Inc(J);
                    end;
                    Inc(I);
                end;
            end;
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    WasErrors := RErr or CErr or SErr or SPDErr or HErr or HPDErr;
    if  not Silent then
    begin
        Write(Format('TESTING MATRIX GENERATOR'#13#10'',[]));
        Write(Format('REAL TEST:                               ',[]));
        if  not RErr then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('COMPLEX TEST:                            ',[]));
        if  not CErr then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('SYMMETRIC TEST:                          ',[]));
        if  not SErr then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('HERMITIAN TEST:                          ',[]));
        if  not HErr then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('SPD TEST:                                ',[]));
        if  not SPDErr then
        begin
            Write(Format('OK'#13#10'',[]));
        end
        else
        begin
            Write(Format('FAILED'#13#10'',[]));
        end;
        Write(Format('HPD TEST:                                ',[]));
        if  not HPDErr then
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
Test whether matrix is SPD
*************************************************************************)
function IsSPD(A : TReal2DArray; N : AlglibInteger; IsUpper : Boolean):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    AJJ : Double;
    v : Double;
    i_ : AlglibInteger;
begin
    A := DynamicArrayCopy(A);
    
    //
    //     Test the input parameters.
    //
    Assert(N>=0, 'Error in SMatrixCholesky: incorrect function arguments');
    
    //
    //     Quick return if possible
    //
    Result := True;
    if N<=0 then
    begin
        Exit;
    end;
    if IsUpper then
    begin
        
        //
        // Compute the Cholesky factorization A = U'*U.
        //
        J:=0;
        while J<=N-1 do
        begin
            
            //
            // Compute U(J,J) and test for non-positive-definiteness.
            //
            V := 0.0;
            for i_ := 0 to J-1 do
            begin
                V := V + A[i_,J]*A[i_,J];
            end;
            AJJ := A[J,J]-V;
            if AP_FP_Less_Eq(AJJ,0) then
            begin
                Result := False;
                Exit;
            end;
            AJJ := SQRT(AJJ);
            A[J,J] := AJJ;
            
            //
            // Compute elements J+1:N of row J.
            //
            if J<N-1 then
            begin
                I:=J+1;
                while I<=N-1 do
                begin
                    V := 0.0;
                    for i_ := 0 to J-1 do
                    begin
                        V := V + A[i_,I]*A[i_,J];
                    end;
                    A[J,I] := A[J,I]-V;
                    Inc(I);
                end;
                V := 1/AJJ;
                APVMul(@A[J][0], J+1, N-1, V);
            end;
            Inc(J);
        end;
    end
    else
    begin
        
        //
        // Compute the Cholesky factorization A = L*L'.
        //
        J:=0;
        while J<=N-1 do
        begin
            
            //
            // Compute L(J,J) and test for non-positive-definiteness.
            //
            V := APVDotProduct(@A[J][0], 0, J-1, @A[J][0], 0, J-1);
            AJJ := A[J,J]-V;
            if AP_FP_Less_Eq(AJJ,0) then
            begin
                Result := False;
                Exit;
            end;
            AJJ := SQRT(AJJ);
            A[J,J] := AJJ;
            
            //
            // Compute elements J+1:N of column J.
            //
            if J<N-1 then
            begin
                I:=J+1;
                while I<=N-1 do
                begin
                    V := APVDotProduct(@A[I][0], 0, J-1, @A[J][0], 0, J-1);
                    A[I,J] := A[I,J]-V;
                    Inc(I);
                end;
                V := 1/AJJ;
                for i_ := J+1 to N-1 do
                begin
                    A[i_,J] := V*A[i_,J];
                end;
            end;
            Inc(J);
        end;
    end;
end;


function ObsoleteSVDDecomposition(var a : TReal2DArray;
     m : AlglibInteger;
     n : AlglibInteger;
     var w : TReal1DArray;
     var v : TReal2DArray):Boolean;
var
    nm : AlglibInteger;
    MinMN : AlglibInteger;
    l : AlglibInteger;
    k : AlglibInteger;
    j : AlglibInteger;
    jj : AlglibInteger;
    its : AlglibInteger;
    i : AlglibInteger;
    z : Double;
    y : Double;
    x : Double;
    vscale : Double;
    s : Double;
    h : Double;
    g : Double;
    f : Double;
    c : Double;
    anorm : Double;
    rv1 : TReal1DArray;
    Flag : Boolean;
begin
    SetLength(rv1, N+1);
    SetLength(W, N+1);
    SetLength(V, N+1, N+1);
    Result := True;
    if M<N then
    begin
        MinMN := M;
    end
    else
    begin
        MinMN := N;
    end;
    g := 0.0;
    vscale := 0.0;
    anorm := 0.0;
    i:=1;
    while i<=n do
    begin
        l := i+1;
        rv1[i] := vscale*g;
        g := 0;
        s := 0;
        vscale := 0;
        if i<=m then
        begin
            k:=i;
            while k<=m do
            begin
                vscale := vscale+AbsReal(a[k,i]);
                Inc(k);
            end;
            if AP_FP_Neq(vscale,0.0) then
            begin
                k:=i;
                while k<=m do
                begin
                    a[k,i] := a[k,i]/vscale;
                    s := s+a[k,i]*a[k,i];
                    Inc(k);
                end;
                f := a[i,i];
                g := -ExtSign(sqrt(s), f);
                h := f*g-s;
                a[i,i] := f-g;
                if i<>n then
                begin
                    j:=l;
                    while j<=n do
                    begin
                        s := 0.0;
                        k:=i;
                        while k<=m do
                        begin
                            s := s+a[k,i]*a[k,j];
                            Inc(k);
                        end;
                        f := s/h;
                        k:=i;
                        while k<=m do
                        begin
                            a[k,j] := a[k,j]+f*a[k,i];
                            Inc(k);
                        end;
                        Inc(j);
                    end;
                end;
                k:=i;
                while k<=m do
                begin
                    a[k,i] := vscale*a[k,i];
                    Inc(k);
                end;
            end;
        end;
        w[i] := vscale*g;
        g := 0.0;
        s := 0.0;
        vscale := 0.0;
        if (i<=m) and (i<>n) then
        begin
            k:=l;
            while k<=n do
            begin
                vscale := vscale+AbsReal(a[i,k]);
                Inc(k);
            end;
            if AP_FP_Neq(vscale,0.0) then
            begin
                k:=l;
                while k<=n do
                begin
                    a[i,k] := a[i,k]/vscale;
                    s := s+a[i,k]*a[i,k];
                    Inc(k);
                end;
                f := a[i,l];
                g := -ExtSign(sqrt(s), f);
                h := f*g-s;
                a[i,l] := f-g;
                k:=l;
                while k<=n do
                begin
                    rv1[k] := a[i,k]/h;
                    Inc(k);
                end;
                if i<>m then
                begin
                    j:=l;
                    while j<=m do
                    begin
                        s := 0.0;
                        k:=l;
                        while k<=n do
                        begin
                            s := s+a[j,k]*a[i,k];
                            Inc(k);
                        end;
                        k:=l;
                        while k<=n do
                        begin
                            a[j,k] := a[j,k]+s*rv1[k];
                            Inc(k);
                        end;
                        Inc(j);
                    end;
                end;
                k:=l;
                while k<=n do
                begin
                    a[i,k] := vscale*a[i,k];
                    Inc(k);
                end;
            end;
        end;
        anorm := MyMax(anorm, AbsReal(w[i])+AbsReal(rv1[i]));
        Inc(i);
    end;
    i:=n;
    while i>=1 do
    begin
        if i<n then
        begin
            if AP_FP_Neq(g,0.0) then
            begin
                j:=l;
                while j<=n do
                begin
                    v[j,i] := a[i,j]/a[i,l]/g;
                    Inc(j);
                end;
                j:=l;
                while j<=n do
                begin
                    s := 0.0;
                    k:=l;
                    while k<=n do
                    begin
                        s := s+a[i,k]*v[k,j];
                        Inc(k);
                    end;
                    k:=l;
                    while k<=n do
                    begin
                        v[k,j] := v[k,j]+s*v[k,i];
                        Inc(k);
                    end;
                    Inc(j);
                end;
            end;
            j:=l;
            while j<=n do
            begin
                v[i,j] := 0.0;
                v[j,i] := 0.0;
                Inc(j);
            end;
        end;
        v[i,i] := 1.0;
        g := rv1[i];
        l := i;
        Dec(i);
    end;
    i:=MinMN;
    while i>=1 do
    begin
        l := i+1;
        g := w[i];
        if i<n then
        begin
            j:=l;
            while j<=n do
            begin
                a[i,j] := 0.0;
                Inc(j);
            end;
        end;
        if AP_FP_Neq(g,0.0) then
        begin
            g := 1.0/g;
            if i<>n then
            begin
                j:=l;
                while j<=n do
                begin
                    s := 0.0;
                    k:=l;
                    while k<=m do
                    begin
                        s := s+a[k,i]*a[k,j];
                        Inc(k);
                    end;
                    f := s/a[i,i]*g;
                    k:=i;
                    while k<=m do
                    begin
                        a[k,j] := a[k,j]+f*a[k,i];
                        Inc(k);
                    end;
                    Inc(j);
                end;
            end;
            j:=i;
            while j<=m do
            begin
                a[j,i] := a[j,i]*g;
                Inc(j);
            end;
        end
        else
        begin
            j:=i;
            while j<=m do
            begin
                a[j,i] := 0.0;
                Inc(j);
            end;
        end;
        a[i,i] := a[i,i]+1.0;
        Dec(i);
    end;
    k:=n;
    while k>=1 do
    begin
        its:=1;
        while its<=MaxSVDIterations do
        begin
            Flag := True;
            l:=k;
            while l>=1 do
            begin
                nm := l-1;
                if AP_FP_Eq(AbsReal(rv1[l])+anorm,anorm) then
                begin
                    Flag := False;
                    Break;
                end;
                if AP_FP_Eq(AbsReal(w[nm])+anorm,anorm) then
                begin
                    Break;
                end;
                Dec(l);
            end;
            if Flag then
            begin
                c := 0.0;
                s := 1.0;
                i:=l;
                while i<=k do
                begin
                    f := s*rv1[i];
                    if AP_FP_Neq(AbsReal(f)+anorm,anorm) then
                    begin
                        g := w[i];
                        h := Pythag(f, g);
                        w[i] := h;
                        h := 1.0/h;
                        c := g*h;
                        s := -f*h;
                        j:=1;
                        while j<=m do
                        begin
                            y := a[j,nm];
                            z := a[j,i];
                            a[j,nm] := y*c+z*s;
                            a[j,i] := -y*s+z*c;
                            Inc(j);
                        end;
                    end;
                    Inc(i);
                end;
            end;
            z := w[k];
            if l=k then
            begin
                if AP_FP_Less(z,0.0) then
                begin
                    w[k] := -z;
                    j:=1;
                    while j<=n do
                    begin
                        v[j,k] := -v[j,k];
                        Inc(j);
                    end;
                end;
                Break;
            end;
            if its=MaxSVDIterations then
            begin
                Result := False;
                Exit;
            end;
            x := w[l];
            nm := k-1;
            y := w[nm];
            g := rv1[nm];
            h := rv1[k];
            f := ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
            g := Pythag(f, 1);
            f := ((x-z)*(x+z)+h*(y/(f+ExtSign(g, f))-h))/x;
            c := 1.0;
            s := 1.0;
            j:=l;
            while j<=nm do
            begin
                i := j+1;
                g := rv1[i];
                y := w[i];
                h := s*g;
                g := c*g;
                z := Pythag(f, h);
                rv1[j] := z;
                c := f/z;
                s := h/z;
                f := x*c+g*s;
                g := -x*s+g*c;
                h := y*s;
                y := y*c;
                jj:=1;
                while jj<=n do
                begin
                    x := v[jj,j];
                    z := v[jj,i];
                    v[jj,j] := x*c+z*s;
                    v[jj,i] := -x*s+z*c;
                    Inc(jj);
                end;
                z := Pythag(f, h);
                w[j] := z;
                if AP_FP_Neq(z,0.0) then
                begin
                    z := 1.0/z;
                    c := f*z;
                    s := h*z;
                end;
                f := c*g+s*y;
                x := -s*g+c*y;
                jj:=1;
                while jj<=m do
                begin
                    y := a[jj,j];
                    z := a[jj,i];
                    a[jj,j] := y*c+z*s;
                    a[jj,i] := -y*s+z*c;
                    Inc(jj);
                end;
                Inc(j);
            end;
            rv1[l] := 0.0;
            rv1[k] := f;
            w[k] := x;
            Inc(its);
        end;
        Dec(k);
    end;
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
Unsets 2D array.
*************************************************************************)
procedure Unset2DC(var A : TComplex2DArray);
begin
    SetLength(A, 0+1, 0+1);
    A[0,0] := C_Complex(2*RandomReal-1);
end;


(*************************************************************************
Tests whether A is HPD
*************************************************************************)
function IsHPD(A : TComplex2DArray; N : AlglibInteger):Boolean;
var
    J : AlglibInteger;
    AJJ : Double;
    V : Complex;
    R : Double;
    T : TComplex1DArray;
    T2 : TComplex1DArray;
    T3 : TComplex1DArray;
    I : AlglibInteger;
    A1 : TComplex2DArray;
    i_ : AlglibInteger;
begin
    A := DynamicArrayCopy(A);
    SetLength(T, N-1+1);
    SetLength(T2, N-1+1);
    SetLength(T3, N-1+1);
    Result := True;
    
    //
    // Compute the Cholesky factorization A = U'*U.
    //
    J:=0;
    while J<=N-1 do
    begin
        
        //
        // Compute U(J,J) and test for non-positive-definiteness.
        //
        V := C_Complex(0.0);
        for i_ := 0 to J-1 do
        begin
            V := C_Add(V,C_Mul(Conj(A[i_,J]),A[i_,J]));
        end;
        AJJ := C_Sub(A[J,J],V).X;
        if AP_FP_Less_Eq(AJJ,0) then
        begin
            A[J,J] := C_Complex(AJJ);
            Result := False;
            Exit;
        end;
        AJJ := SQRT(AJJ);
        A[J,J] := C_Complex(AJJ);
        
        //
        // Compute elements J+1:N-1 of row J.
        //
        if J<N-1 then
        begin
            for i_ := 0 to J-1 do
            begin
                T2[i_] := Conj(A[i_,J]);
            end;
            for i_ := J+1 to N-1 do
            begin
                T3[i_] := A[J,i_];
            end;
            I:=J+1;
            while I<=N-1 do
            begin
                V := C_Complex(0.0);
                for i_ := 0 to J-1 do
                begin
                    V := C_Add(V,C_Mul(A[i_,I],T2[i_]));
                end;
                T3[I] := C_Sub(T3[I],V);
                Inc(I);
            end;
            for i_ := J+1 to N-1 do
            begin
                A[J,i_] := T3[i_];
            end;
            R := 1/AJJ;
            for i_ := J+1 to N-1 do
            begin
                A[J,i_] := C_MulR(A[J,i_],R);
            end;
        end;
        Inc(J);
    end;
end;


(*************************************************************************
SVD condition number
*************************************************************************)
function SVDCond(const a : TReal2DArray; N : AlglibInteger):Double;
var
    a1 : TReal2DArray;
    v : TReal2DArray;
    w : TReal1DArray;
    i : AlglibInteger;
    j : AlglibInteger;
    MinW : Double;
    MaxW : Double;
begin
    SetLength(a1, n+1, n+1);
    i:=1;
    while i<=n do
    begin
        j:=1;
        while j<=n do
        begin
            a1[i,j] := a[i-1,j-1];
            Inc(j);
        end;
        Inc(i);
    end;
    if  not ObsoleteSVDDecomposition(a1, n, n, w, v) then
    begin
        result := 0;
        Exit;
    end;
    MinW := W[1];
    MaxW := W[1];
    I:=2;
    while I<=N do
    begin
        if AP_FP_Less(W[I],MinW) then
        begin
            MinW := W[I];
        end;
        if AP_FP_Greater(W[I],MaxW) then
        begin
            MaxW := W[I];
        end;
        Inc(I);
    end;
    Result := MaxW/MinW;
end;


function ExtSign(a : Double; b : Double):Double;
begin
    if AP_FP_Greater_Eq(b,0) then
    begin
        Result := AbsReal(a);
    end
    else
    begin
        Result := -AbsReal(a);
    end;
end;


function MyMax(a : Double; b : Double):Double;
begin
    if AP_FP_Greater(a,b) then
    begin
        Result := a;
    end
    else
    begin
        Result := b;
    end;
end;


function Pythag(A : Double; B : Double):Double;
begin
    if AP_FP_Less(AbsReal(A),AbsReal(B)) then
    begin
        Result := AbsReal(B)*Sqrt(1+AP_Sqr(A/B));
    end
    else
    begin
        Result := AbsReal(A)*Sqrt(1+AP_Sqr(B/A));
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testmatgenunit_test_silent():Boolean;
begin
    Result := TestMatGen(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testmatgenunit_test():Boolean;
begin
    Result := TestMatGen(False);
end;


end.