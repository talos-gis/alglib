unit testablasunit;
interface
uses Math, Sysutils, Ap, ablasf, ablas;

function TestABLAS(Silent : Boolean):Boolean;
procedure RefCMatrixRightTRSM(M : AlglibInteger;
     N : AlglibInteger;
     const A : TComplex2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TComplex2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger);
procedure RefCMatrixLeftTRSM(M : AlglibInteger;
     N : AlglibInteger;
     const A : TComplex2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TComplex2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger);
procedure RefRMatrixRightTRSM(M : AlglibInteger;
     N : AlglibInteger;
     const A : TReal2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TReal2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger);
procedure RefRMatrixLeftTRSM(M : AlglibInteger;
     N : AlglibInteger;
     const A : TReal2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TReal2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger);
function InternalCMatrixTRInverse(var A : TComplex2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsunitTriangular : Boolean):Boolean;
function InternalRMatrixTRInverse(var A : TReal2DArray;
     N : AlglibInteger;
     IsUpper : Boolean;
     IsunitTriangular : Boolean):Boolean;
procedure RefCMatrixSYRK(N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Double;
     const A : TComplex2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     Beta : Double;
     var C : TComplex2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger;
     IsUpper : Boolean);
procedure RefRMatrixSYRK(N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Double;
     const A : TReal2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     Beta : Double;
     var C : TReal2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger;
     IsUpper : Boolean);
procedure RefCMatrixGEMM(M : AlglibInteger;
     N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Complex;
     const A : TComplex2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     const B : TComplex2DArray;
     IB : AlglibInteger;
     JB : AlglibInteger;
     OpTypeB : AlglibInteger;
     Beta : Complex;
     var C : TComplex2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger);
procedure RefRMatrixGEMM(M : AlglibInteger;
     N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Double;
     const A : TReal2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     const B : TReal2DArray;
     IB : AlglibInteger;
     JB : AlglibInteger;
     OpTypeB : AlglibInteger;
     Beta : Double;
     var C : TReal2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger);
function testablasunit_test_silent():Boolean;
function testablasunit_test():Boolean;

implementation

procedure NaiveMatrixMatrixMultiply(const A : TReal2DArray;
     AI1 : AlglibInteger;
     AI2 : AlglibInteger;
     AJ1 : AlglibInteger;
     AJ2 : AlglibInteger;
     TransA : Boolean;
     const B : TReal2DArray;
     BI1 : AlglibInteger;
     BI2 : AlglibInteger;
     BJ1 : AlglibInteger;
     BJ2 : AlglibInteger;
     TransB : Boolean;
     Alpha : Double;
     var C : TReal2DArray;
     CI1 : AlglibInteger;
     CI2 : AlglibInteger;
     CJ1 : AlglibInteger;
     CJ2 : AlglibInteger;
     Beta : Double);forward;
function TestTRSM(MinN : AlglibInteger; MaxN : AlglibInteger):Boolean;forward;
function TestSYRK(MinN : AlglibInteger; MaxN : AlglibInteger):Boolean;forward;
function TestGEMM(MinN : AlglibInteger; MaxN : AlglibInteger):Boolean;forward;
function TestTrans(MinN : AlglibInteger; MaxN : AlglibInteger):Boolean;forward;
function TestRANK1(MinN : AlglibInteger; MaxN : AlglibInteger):Boolean;forward;
function TestMV(MinN : AlglibInteger; MaxN : AlglibInteger):Boolean;forward;
function TestCopy(MinN : AlglibInteger; MaxN : AlglibInteger):Boolean;forward;


function TestABLAS(Silent : Boolean):Boolean;
var
    Threshold : Double;
    TRSMErrors : Boolean;
    SYRKErrors : Boolean;
    GEMMErrors : Boolean;
    TRANSErrors : Boolean;
    RANK1Errors : Boolean;
    MVErrors : Boolean;
    CopyErrors : Boolean;
    WasErrors : Boolean;
    RA : TReal2DArray;
begin
    TRSMErrors := False;
    SYRKErrors := False;
    GEMMErrors := False;
    TRANSErrors := False;
    RANK1Errors := False;
    MVErrors := False;
    CopyErrors := False;
    WasErrors := False;
    Threshold := 10000*MachineEpsilon;
    TRSMErrors := TRSMErrors or TestTRSM(1, 3*ABLASBlockSize(RA)+1);
    SYRKErrors := SYRKErrors or TestSYRK(1, 3*ABLASBlockSize(RA)+1);
    GEMMErrors := GEMMErrors or TestGEMM(1, 3*ABLASBlockSize(RA)+1);
    TRANSErrors := TRANSErrors or TestTRANS(1, 3*ABLASBlockSize(RA)+1);
    RANK1Errors := RANK1Errors or TestRANK1(1, 3*ABLASBlockSize(RA)+1);
    MVErrors := MVErrors or TestMV(1, 3*ABLASBlockSize(RA)+1);
    CopyErrors := CopyErrors or TestCopy(1, 3*ABLASBlockSize(RA)+1);
    
    //
    // report
    //
    WasErrors := TRSMErrors or SYRKErrors or GEMMErrors or TRANSErrors or RANK1Errors or MVErrors or CopyErrors;
    if  not Silent then
    begin
        Write(Format('TESTING ABLAS'#13#10'',[]));
        Write(Format('* TRSM:                                  ',[]));
        if TRSMErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* SYRK:                                  ',[]));
        if SYRKErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* GEMM:                                  ',[]));
        if GEMMErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* TRANS:                                 ',[]));
        if TRANSErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* RANK1:                                 ',[]));
        if RANK1Errors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* MV:                                    ',[]));
        if MVErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* COPY:                                  ',[]));
        if CopyErrors then
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
Reference implementation

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure RefCMatrixRightTRSM(M : AlglibInteger;
     N : AlglibInteger;
     const A : TComplex2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TComplex2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger);
var
    A1 : TComplex2DArray;
    A2 : TComplex2DArray;
    TX : TComplex1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    VC : Complex;
    RUpper : Boolean;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if N*M=0 then
    begin
        Exit;
    end;
    SetLength(A1, N, N);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            A1[I,J] := C_Complex(0);
            Inc(J);
        end;
        Inc(I);
    end;
    if IsUpper then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=I;
            while J<=N-1 do
            begin
                A1[I,J] := A[I1+I,J1+J];
                Inc(J);
            end;
            Inc(I);
        end;
    end
    else
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=I do
            begin
                A1[I,J] := A[I1+I,J1+J];
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    RUpper := IsUpper;
    if IsUnit then
    begin
        I:=0;
        while I<=N-1 do
        begin
            A1[I,I] := C_Complex(1);
            Inc(I);
        end;
    end;
    SetLength(A2, N, N);
    if OpType=0 then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                A2[I,J] := A1[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    if OpType=1 then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                A2[I,J] := A1[J,I];
                Inc(J);
            end;
            Inc(I);
        end;
        RUpper :=  not RUpper;
    end;
    if OpType=2 then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                A2[I,J] := Conj(A1[J,I]);
                Inc(J);
            end;
            Inc(I);
        end;
        RUpper :=  not RUpper;
    end;
    InternalCMatrixTRInverse(A2, N, RUpper, False);
    SetLength(TX, N);
    I:=0;
    while I<=M-1 do
    begin
        i1_ := (J2) - (0);
        for i_ := 0 to N-1 do
        begin
            TX[i_] := X[I2+I,i_+i1_];
        end;
        J:=0;
        while J<=N-1 do
        begin
            VC := C_Complex(0.0);
            for i_ := 0 to N-1 do
            begin
                VC := C_Add(VC,C_Mul(TX[i_],A2[i_,J]));
            end;
            X[I2+I,J2+J] := VC;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Reference implementation

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure RefCMatrixLeftTRSM(M : AlglibInteger;
     N : AlglibInteger;
     const A : TComplex2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TComplex2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger);
var
    A1 : TComplex2DArray;
    A2 : TComplex2DArray;
    TX : TComplex1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    VC : Complex;
    RUpper : Boolean;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if N*M=0 then
    begin
        Exit;
    end;
    SetLength(A1, M, M);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=M-1 do
        begin
            A1[I,J] := C_Complex(0);
            Inc(J);
        end;
        Inc(I);
    end;
    if IsUpper then
    begin
        I:=0;
        while I<=M-1 do
        begin
            J:=I;
            while J<=M-1 do
            begin
                A1[I,J] := A[I1+I,J1+J];
                Inc(J);
            end;
            Inc(I);
        end;
    end
    else
    begin
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=I do
            begin
                A1[I,J] := A[I1+I,J1+J];
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    RUpper := IsUpper;
    if IsUnit then
    begin
        I:=0;
        while I<=M-1 do
        begin
            A1[I,I] := C_Complex(1);
            Inc(I);
        end;
    end;
    SetLength(A2, M, M);
    if OpType=0 then
    begin
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                A2[I,J] := A1[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    if OpType=1 then
    begin
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                A2[I,J] := A1[J,I];
                Inc(J);
            end;
            Inc(I);
        end;
        RUpper :=  not RUpper;
    end;
    if OpType=2 then
    begin
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                A2[I,J] := Conj(A1[J,I]);
                Inc(J);
            end;
            Inc(I);
        end;
        RUpper :=  not RUpper;
    end;
    InternalCMatrixTRInverse(A2, M, RUpper, False);
    SetLength(TX, M);
    J:=0;
    while J<=N-1 do
    begin
        i1_ := (I2) - (0);
        for i_ := 0 to M-1 do
        begin
            TX[i_] := X[i_+i1_,J2+J];
        end;
        I:=0;
        while I<=M-1 do
        begin
            VC := C_Complex(0.0);
            for i_ := 0 to M-1 do
            begin
                VC := C_Add(VC,C_Mul(A2[I,i_],TX[i_]));
            end;
            X[I2+I,J2+J] := VC;
            Inc(I);
        end;
        Inc(J);
    end;
end;


(*************************************************************************
Reference implementation

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure RefRMatrixRightTRSM(M : AlglibInteger;
     N : AlglibInteger;
     const A : TReal2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TReal2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger);
var
    A1 : TReal2DArray;
    A2 : TReal2DArray;
    TX : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    VR : Double;
    RUpper : Boolean;
    i_ : AlglibInteger;
begin
    if N*M=0 then
    begin
        Exit;
    end;
    SetLength(A1, N, N);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            A1[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    if IsUpper then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=I;
            while J<=N-1 do
            begin
                A1[I,J] := A[I1+I,J1+J];
                Inc(J);
            end;
            Inc(I);
        end;
    end
    else
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=I do
            begin
                A1[I,J] := A[I1+I,J1+J];
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    RUpper := IsUpper;
    if IsUnit then
    begin
        I:=0;
        while I<=N-1 do
        begin
            A1[I,I] := 1;
            Inc(I);
        end;
    end;
    SetLength(A2, N, N);
    if OpType=0 then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                A2[I,J] := A1[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    if OpType=1 then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                A2[I,J] := A1[J,I];
                Inc(J);
            end;
            Inc(I);
        end;
        RUpper :=  not RUpper;
    end;
    InternalRMatrixTRInverse(A2, N, RUpper, False);
    SetLength(TX, N);
    I:=0;
    while I<=M-1 do
    begin
        APVMove(@TX[0], 0, N-1, @X[I2+I][0], J2, J2+N-1);
        J:=0;
        while J<=N-1 do
        begin
            VR := 0.0;
            for i_ := 0 to N-1 do
            begin
                VR := VR + TX[i_]*A2[i_,J];
            end;
            X[I2+I,J2+J] := VR;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Reference implementation

  -- ALGLIB routine --
     15.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure RefRMatrixLeftTRSM(M : AlglibInteger;
     N : AlglibInteger;
     const A : TReal2DArray;
     I1 : AlglibInteger;
     J1 : AlglibInteger;
     IsUpper : Boolean;
     IsUnit : Boolean;
     OpType : AlglibInteger;
     var X : TReal2DArray;
     I2 : AlglibInteger;
     J2 : AlglibInteger);
var
    A1 : TReal2DArray;
    A2 : TReal2DArray;
    TX : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    VR : Double;
    RUpper : Boolean;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    if N*M=0 then
    begin
        Exit;
    end;
    SetLength(A1, M, M);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=M-1 do
        begin
            A1[I,J] := 0;
            Inc(J);
        end;
        Inc(I);
    end;
    if IsUpper then
    begin
        I:=0;
        while I<=M-1 do
        begin
            J:=I;
            while J<=M-1 do
            begin
                A1[I,J] := A[I1+I,J1+J];
                Inc(J);
            end;
            Inc(I);
        end;
    end
    else
    begin
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=I do
            begin
                A1[I,J] := A[I1+I,J1+J];
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    RUpper := IsUpper;
    if IsUnit then
    begin
        I:=0;
        while I<=M-1 do
        begin
            A1[I,I] := 1;
            Inc(I);
        end;
    end;
    SetLength(A2, M, M);
    if OpType=0 then
    begin
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                A2[I,J] := A1[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    if OpType=1 then
    begin
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                A2[I,J] := A1[J,I];
                Inc(J);
            end;
            Inc(I);
        end;
        RUpper :=  not RUpper;
    end;
    InternalRMatrixTRInverse(A2, M, RUpper, False);
    SetLength(TX, M);
    J:=0;
    while J<=N-1 do
    begin
        i1_ := (I2) - (0);
        for i_ := 0 to M-1 do
        begin
            TX[i_] := X[i_+i1_,J2+J];
        end;
        I:=0;
        while I<=M-1 do
        begin
            VR := APVDotProduct(@A2[I][0], 0, M-1, @TX[0], 0, M-1);
            X[I2+I,J2+J] := VR;
            Inc(I);
        end;
        Inc(J);
    end;
end;


(*************************************************************************
Internal subroutine.
Triangular matrix inversion

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************)
function InternalCMatrixTRInverse(var A : TComplex2DArray;
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
                    if I+1<J then
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
            if J+1<N then
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
Internal subroutine.
Triangular matrix inversion

  -- LAPACK routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     February 29, 1992
*************************************************************************)
function InternalRMatrixTRInverse(var A : TReal2DArray;
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
Reference SYRK subroutine.

  -- ALGLIB routine --
     16.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure RefCMatrixSYRK(N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Double;
     const A : TComplex2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     Beta : Double;
     var C : TComplex2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger;
     IsUpper : Boolean);
var
    AE : TComplex2DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    VC : Complex;
    i_ : AlglibInteger;
begin
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if IsUpper and (J>=I) or  not IsUpper and (J<=I) then
            begin
                if AP_FP_Eq(Beta,0) then
                begin
                    C[I+IC,J+JC] := C_Complex(0);
                end
                else
                begin
                    C[I+IC,J+JC] := C_MulR(C[I+IC,J+JC],Beta);
                end;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(Alpha,0) then
    begin
        Exit;
    end;
    if N*K>0 then
    begin
        SetLength(AE, N, K);
    end;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=K-1 do
        begin
            if OpTypeA=0 then
            begin
                AE[I,J] := A[IA+I,JA+J];
            end;
            if OpTypeA=2 then
            begin
                AE[I,J] := Conj(A[IA+J,JA+I]);
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
            VC := C_Complex(0);
            if K>0 then
            begin
                VC := C_Complex(0.0);
                for i_ := 0 to K-1 do
                begin
                    VC := C_Add(VC,C_Mul(AE[I,i_],Conj(AE[J,i_])));
                end;
            end;
            VC := C_MulR(VC,Alpha);
            if IsUpper and (J>=I) then
            begin
                C[IC+I,JC+J] := C_Add(VC,C[IC+I,JC+J]);
            end;
            if  not IsUpper and (J<=I) then
            begin
                C[IC+I,JC+J] := C_Add(VC,C[IC+I,JC+J]);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Reference SYRK subroutine.

  -- ALGLIB routine --
     16.12.2009
     Bochkanov Sergey
*************************************************************************)
procedure RefRMatrixSYRK(N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Double;
     const A : TReal2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     Beta : Double;
     var C : TReal2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger;
     IsUpper : Boolean);
var
    AE : TReal2DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    VR : Double;
begin
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if IsUpper and (J>=I) or  not IsUpper and (J<=I) then
            begin
                if AP_FP_Eq(Beta,0) then
                begin
                    C[I+IC,J+JC] := 0;
                end
                else
                begin
                    C[I+IC,J+JC] := C[I+IC,J+JC]*Beta;
                end;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    if AP_FP_Eq(Alpha,0) then
    begin
        Exit;
    end;
    if N*K>0 then
    begin
        SetLength(AE, N, K);
    end;
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=K-1 do
        begin
            if OpTypeA=0 then
            begin
                AE[I,J] := A[IA+I,JA+J];
            end;
            if OpTypeA=1 then
            begin
                AE[I,J] := A[IA+J,JA+I];
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
            VR := 0;
            if K>0 then
            begin
                VR := APVDotProduct(@AE[I][0], 0, K-1, @AE[J][0], 0, K-1);
            end;
            VR := Alpha*VR;
            if IsUpper and (J>=I) then
            begin
                C[IC+I,JC+J] := VR+C[IC+I,JC+J];
            end;
            if  not IsUpper and (J<=I) then
            begin
                C[IC+I,JC+J] := VR+C[IC+I,JC+J];
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Reference GEMM,
ALGLIB subroutine
*************************************************************************)
procedure RefCMatrixGEMM(M : AlglibInteger;
     N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Complex;
     const A : TComplex2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     const B : TComplex2DArray;
     IB : AlglibInteger;
     JB : AlglibInteger;
     OpTypeB : AlglibInteger;
     Beta : Complex;
     var C : TComplex2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger);
var
    AE : TComplex2DArray;
    BE : TComplex2DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    VC : Complex;
    i_ : AlglibInteger;
begin
    SetLength(AE, M, K);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=K-1 do
        begin
            if OpTypeA=0 then
            begin
                AE[I,J] := A[IA+I,JA+J];
            end;
            if OpTypeA=1 then
            begin
                AE[I,J] := A[IA+J,JA+I];
            end;
            if OpTypeA=2 then
            begin
                AE[I,J] := Conj(A[IA+J,JA+I]);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    SetLength(BE, K, N);
    I:=0;
    while I<=K-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if OpTypeB=0 then
            begin
                BE[I,J] := B[IB+I,JB+J];
            end;
            if OpTypeB=1 then
            begin
                BE[I,J] := B[IB+J,JB+I];
            end;
            if OpTypeB=2 then
            begin
                BE[I,J] := Conj(B[IB+J,JB+I]);
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            VC := C_Complex(0.0);
            for i_ := 0 to K-1 do
            begin
                VC := C_Add(VC,C_Mul(AE[I,i_],BE[i_,J]));
            end;
            VC := C_Mul(Alpha,VC);
            if C_NotEqualR(Beta,0) then
            begin
                VC := C_Add(VC,C_Mul(Beta,C[IC+I,JC+J]));
            end;
            C[IC+I,JC+J] := VC;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Reference GEMM,
ALGLIB subroutine
*************************************************************************)
procedure RefRMatrixGEMM(M : AlglibInteger;
     N : AlglibInteger;
     K : AlglibInteger;
     Alpha : Double;
     const A : TReal2DArray;
     IA : AlglibInteger;
     JA : AlglibInteger;
     OpTypeA : AlglibInteger;
     const B : TReal2DArray;
     IB : AlglibInteger;
     JB : AlglibInteger;
     OpTypeB : AlglibInteger;
     Beta : Double;
     var C : TReal2DArray;
     IC : AlglibInteger;
     JC : AlglibInteger);
var
    AE : TReal2DArray;
    BE : TReal2DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    VC : Double;
    i_ : AlglibInteger;
begin
    SetLength(AE, M, K);
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=K-1 do
        begin
            if OpTypeA=0 then
            begin
                AE[I,J] := A[IA+I,JA+J];
            end;
            if OpTypeA=1 then
            begin
                AE[I,J] := A[IA+J,JA+I];
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    SetLength(BE, K, N);
    I:=0;
    while I<=K-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            if OpTypeB=0 then
            begin
                BE[I,J] := B[IB+I,JB+J];
            end;
            if OpTypeB=1 then
            begin
                BE[I,J] := B[IB+J,JB+I];
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            VC := 0.0;
            for i_ := 0 to K-1 do
            begin
                VC := VC + AE[I,i_]*BE[i_,J];
            end;
            VC := Alpha*VC;
            if AP_FP_Neq(Beta,0) then
            begin
                VC := VC+Beta*C[IC+I,JC+J];
            end;
            C[IC+I,JC+J] := VC;
            Inc(J);
        end;
        Inc(I);
    end;
end;


procedure NaiveMatrixMatrixMultiply(const A : TReal2DArray;
     AI1 : AlglibInteger;
     AI2 : AlglibInteger;
     AJ1 : AlglibInteger;
     AJ2 : AlglibInteger;
     TransA : Boolean;
     const B : TReal2DArray;
     BI1 : AlglibInteger;
     BI2 : AlglibInteger;
     BJ1 : AlglibInteger;
     BJ2 : AlglibInteger;
     TransB : Boolean;
     Alpha : Double;
     var C : TReal2DArray;
     CI1 : AlglibInteger;
     CI2 : AlglibInteger;
     CJ1 : AlglibInteger;
     CJ2 : AlglibInteger;
     Beta : Double);
var
    ARows : AlglibInteger;
    ACols : AlglibInteger;
    BRows : AlglibInteger;
    BCols : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    L : AlglibInteger;
    R : AlglibInteger;
    V : Double;
    X1 : TReal1DArray;
    X2 : TReal1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    
    //
    // Setup
    //
    if  not TransA then
    begin
        ARows := AI2-AI1+1;
        ACols := AJ2-AJ1+1;
    end
    else
    begin
        ARows := AJ2-AJ1+1;
        ACols := AI2-AI1+1;
    end;
    if  not TransB then
    begin
        BRows := BI2-BI1+1;
        BCols := BJ2-BJ1+1;
    end
    else
    begin
        BRows := BJ2-BJ1+1;
        BCols := BI2-BI1+1;
    end;
    Assert(ACols=BRows, 'NaiveMatrixMatrixMultiply: incorrect matrix sizes!');
    if (ARows<=0) or (ACols<=0) or (BRows<=0) or (BCols<=0) then
    begin
        Exit;
    end;
    L := ARows;
    R := BCols;
    K := ACols;
    SetLength(X1, K+1);
    SetLength(X2, K+1);
    I:=1;
    while I<=L do
    begin
        J:=1;
        while J<=R do
        begin
            if  not TransA then
            begin
                if  not TransB then
                begin
                    i1_ := (AJ1)-(BI1);
                    V := 0.0;
                    for i_ := BI1 to BI2 do
                    begin
                        V := V + B[i_,BJ1+J-1]*A[AI1+I-1,i_+i1_];
                    end;
                end
                else
                begin
                    V := APVDotProduct(@B[BI1+J-1][0], BJ1, BJ2, @A[AI1+I-1][0], AJ1, AJ2);
                end;
            end
            else
            begin
                if  not TransB then
                begin
                    i1_ := (AI1)-(BI1);
                    V := 0.0;
                    for i_ := BI1 to BI2 do
                    begin
                        V := V + B[i_,BJ1+J-1]*A[i_+i1_,AJ1+I-1];
                    end;
                end
                else
                begin
                    i1_ := (AI1)-(BJ1);
                    V := 0.0;
                    for i_ := BJ1 to BJ2 do
                    begin
                        V := V + B[BI1+J-1,i_]*A[i_+i1_,AJ1+I-1];
                    end;
                end;
            end;
            if AP_FP_Eq(Beta,0) then
            begin
                C[CI1+I-1,CJ1+J-1] := Alpha*V;
            end
            else
            begin
                C[CI1+I-1,CJ1+J-1] := Beta*C[CI1+I-1,CJ1+J-1]+Alpha*V;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
?Matrix????TRSM tests

Returns False for passed test, True - for failed
*************************************************************************)
function TestTRSM(MinN : AlglibInteger; MaxN : AlglibInteger):Boolean;
var
    N : AlglibInteger;
    M : AlglibInteger;
    MX : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    OpType : AlglibInteger;
    UpperType : AlglibInteger;
    UnitType : AlglibInteger;
    XOffsI : AlglibInteger;
    XOffsJ : AlglibInteger;
    AOffsIType : AlglibInteger;
    AOffsJType : AlglibInteger;
    AOffsI : AlglibInteger;
    AOffsJ : AlglibInteger;
    RefRA : TReal2DArray;
    RefRXL : TReal2DArray;
    RefRXR : TReal2DArray;
    RefCA : TComplex2DArray;
    RefCXL : TComplex2DArray;
    RefCXR : TComplex2DArray;
    RA : TReal2DArray;
    CA : TComplex2DArray;
    RXR1 : TReal2DArray;
    RXL1 : TReal2DArray;
    CXR1 : TComplex2DArray;
    CXL1 : TComplex2DArray;
    RXR2 : TReal2DArray;
    RXL2 : TReal2DArray;
    CXR2 : TComplex2DArray;
    CXL2 : TComplex2DArray;
    Threshold : Double;
begin
    Threshold := AP_Sqr(MaxN)*100*MachineEpsilon;
    Result := False;
    MX:=MinN;
    while MX<=MaxN do
    begin
        
        //
        // Select random M/N in [1,MX] such that max(M,N)=MX
        //
        M := 1+RandomInteger(MX);
        N := 1+RandomInteger(MX);
        if AP_FP_Greater(RandomReal,0.5) then
        begin
            M := MX;
        end
        else
        begin
            N := MX;
        end;
        
        //
        // Initialize RefRA/RefCA by random matrices whose upper
        // and lower triangle submatrices are non-degenerate
        // well-conditioned matrices.
        //
        // Matrix size is 2Mx2M (four copies of same MxM matrix
        // to test different offsets)
        //
        SetLength(RefRA, 2*M, 2*M);
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                RefRA[I,J] := 0.2*RandomReal-0.1;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=0;
        while I<=M-1 do
        begin
            RefRA[I,I] := (2*RandomInteger(1)-1)*(2*M+RandomReal);
            Inc(I);
        end;
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                RefRA[I+M,J] := RefRA[I,J];
                RefRA[I,J+M] := RefRA[I,J];
                RefRA[I+M,J+M] := RefRA[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
        SetLength(RefCA, 2*M, 2*M);
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                RefCA[I,J].X := 0.2*RandomReal-0.1;
                RefCA[I,J].Y := 0.2*RandomReal-0.1;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=0;
        while I<=M-1 do
        begin
            RefCA[I,I].X := (2*RandomInteger(2)-1)*(2*M+RandomReal);
            RefCA[I,I].Y := (2*RandomInteger(2)-1)*(2*M+RandomReal);
            Inc(I);
        end;
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                RefCA[I+M,J] := RefCA[I,J];
                RefCA[I,J+M] := RefCA[I,J];
                RefCA[I+M,J+M] := RefCA[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Generate random XL/XR.
        //
        // XR is NxM matrix (matrix for 'Right' subroutines)
        // XL is MxN matrix (matrix for 'Left' subroutines)
        //
        SetLength(RefRXR, N, M);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                RefRXR[I,J] := 2*RandomReal-1;
                Inc(J);
            end;
            Inc(I);
        end;
        SetLength(RefRXL, M, N);
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                RefRXL[I,J] := 2*RandomReal-1;
                Inc(J);
            end;
            Inc(I);
        end;
        SetLength(RefCXR, N, M);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                RefCXR[I,J].X := 2*RandomReal-1;
                RefCXR[I,J].Y := 2*RandomReal-1;
                Inc(J);
            end;
            Inc(I);
        end;
        SetLength(RefCXL, M, N);
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                RefCXL[I,J].X := 2*RandomReal-1;
                RefCXL[I,J].Y := 2*RandomReal-1;
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // test different types of operations, offsets, and so on...
        //
        // to avoid unnecessary slowdown we don't test ALL possible
        // combinations of operation types. We just generate one random
        // set of parameters and test it.
        //
        SetLength(RA, 2*M, 2*M);
        SetLength(RXR1, N, M);
        SetLength(RXR2, N, M);
        SetLength(RXL1, M, N);
        SetLength(RXL2, M, N);
        SetLength(CA, 2*M, 2*M);
        SetLength(CXR1, N, M);
        SetLength(CXR2, N, M);
        SetLength(CXL1, M, N);
        SetLength(CXL2, M, N);
        OpType := RandomInteger(3);
        UpperType := RandomInteger(2);
        UnitType := RandomInteger(2);
        XOffsI := RandomInteger(2);
        XOffsJ := RandomInteger(2);
        AOffsIType := RandomInteger(2);
        AOffsJType := RandomInteger(2);
        AOffsI := M*AOffsIType;
        AOffsJ := M*AOffsJType;
        
        //
        // copy A, XR, XL (fill unused parts with random garbage)
        //
        I:=0;
        while I<=2*M-1 do
        begin
            J:=0;
            while J<=2*M-1 do
            begin
                if (I>=AOffsI) and (I<AOffsI+M) and (J>=AOffsJ) and (J<AOffsJ+M) then
                begin
                    CA[I,J] := RefCA[I,J];
                    RA[I,J] := RefRA[I,J];
                end
                else
                begin
                    CA[I,J] := C_Complex(RandomReal);
                    RA[I,J] := RandomReal;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                if (I>=XOffsI) and (J>=XOffsJ) then
                begin
                    CXR1[I,J] := RefCXR[I,J];
                    CXR2[I,J] := RefCXR[I,J];
                    RXR1[I,J] := RefRXR[I,J];
                    RXR2[I,J] := RefRXR[I,J];
                end
                else
                begin
                    CXR1[I,J] := C_Complex(RandomReal);
                    CXR2[I,J] := CXR1[I,J];
                    RXR1[I,J] := RandomReal;
                    RXR2[I,J] := RXR1[I,J];
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                if (I>=XOffsI) and (J>=XOffsJ) then
                begin
                    CXL1[I,J] := RefCXL[I,J];
                    CXL2[I,J] := RefCXL[I,J];
                    RXL1[I,J] := RefRXL[I,J];
                    RXL2[I,J] := RefRXL[I,J];
                end
                else
                begin
                    CXL1[I,J] := C_Complex(RandomReal);
                    CXL2[I,J] := CXL1[I,J];
                    RXL1[I,J] := RandomReal;
                    RXL2[I,J] := RXL1[I,J];
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Test CXR
        //
        CMatrixRightTRSM(N-XOffsI, M-XOffsJ, CA, AOffsI, AOffsJ, UpperType=0, UnitType=0, OpType, CXR1, XOffsI, XOffsJ);
        RefCMatrixRightTRSM(N-XOffsI, M-XOffsJ, CA, AOffsI, AOffsJ, UpperType=0, UnitType=0, OpType, CXR2, XOffsI, XOffsJ);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                Result := Result or AP_FP_Greater(AbsComplex(C_Sub(CXR1[I,J],CXR2[I,J])),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Test CXL
        //
        CMatrixLeftTRSM(M-XOffsI, N-XOffsJ, CA, AOffsI, AOffsJ, UpperType=0, UnitType=0, OpType, CXL1, XOffsI, XOffsJ);
        RefCMatrixLeftTRSM(M-XOffsI, N-XOffsJ, CA, AOffsI, AOffsJ, UpperType=0, UnitType=0, OpType, CXL2, XOffsI, XOffsJ);
        I:=0;
        while I<=M-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                Result := Result or AP_FP_Greater(AbsComplex(C_Sub(CXL1[I,J],CXL2[I,J])),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
        if OpType<2 then
        begin
            
            //
            // Test RXR
            //
            RMatrixRightTRSM(N-XOffsI, M-XOffsJ, RA, AOffsI, AOffsJ, UpperType=0, UnitType=0, OpType, RXR1, XOffsI, XOffsJ);
            RefRMatrixRightTRSM(N-XOffsI, M-XOffsJ, RA, AOffsI, AOffsJ, UpperType=0, UnitType=0, OpType, RXR2, XOffsI, XOffsJ);
            I:=0;
            while I<=N-1 do
            begin
                J:=0;
                while J<=M-1 do
                begin
                    Result := Result or AP_FP_Greater(AbsReal(RXR1[I,J]-RXR2[I,J]),Threshold);
                    Inc(J);
                end;
                Inc(I);
            end;
            
            //
            // Test RXL
            //
            RMatrixLeftTRSM(M-XOffsI, N-XOffsJ, RA, AOffsI, AOffsJ, UpperType=0, UnitType=0, OpType, RXL1, XOffsI, XOffsJ);
            RefRMatrixLeftTRSM(M-XOffsI, N-XOffsJ, RA, AOffsI, AOffsJ, UpperType=0, UnitType=0, OpType, RXL2, XOffsI, XOffsJ);
            I:=0;
            while I<=M-1 do
            begin
                J:=0;
                while J<=N-1 do
                begin
                    Result := Result or AP_FP_Greater(AbsReal(RXL1[I,J]-RXL2[I,J]),Threshold);
                    Inc(J);
                end;
                Inc(I);
            end;
        end;
        Inc(MX);
    end;
end;


(*************************************************************************
SYRK tests

Returns False for passed test, True - for failed
*************************************************************************)
function TestSYRK(MinN : AlglibInteger; MaxN : AlglibInteger):Boolean;
var
    N : AlglibInteger;
    K : AlglibInteger;
    MX : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    UpperType : AlglibInteger;
    XOffsI : AlglibInteger;
    XOffsJ : AlglibInteger;
    AOffsIType : AlglibInteger;
    AOffsJType : AlglibInteger;
    AOffsI : AlglibInteger;
    AOffsJ : AlglibInteger;
    AlphaType : AlglibInteger;
    BetaType : AlglibInteger;
    RefRA : TReal2DArray;
    RefRC : TReal2DArray;
    RefCA : TComplex2DArray;
    RefCC : TComplex2DArray;
    Alpha : Double;
    Beta : Double;
    RA1 : TReal2DArray;
    RA2 : TReal2DArray;
    CA1 : TComplex2DArray;
    CA2 : TComplex2DArray;
    RC : TReal2DArray;
    RCT : TReal2DArray;
    CC : TComplex2DArray;
    CCT : TComplex2DArray;
    Threshold : Double;
begin
    Threshold := MaxN*100*MachineEpsilon;
    Result := False;
    MX:=MinN;
    while MX<=MaxN do
    begin
        
        //
        // Select random M/N in [1,MX] such that max(M,N)=MX
        //
        K := 1+RandomInteger(MX);
        N := 1+RandomInteger(MX);
        if AP_FP_Greater(RandomReal,0.5) then
        begin
            K := MX;
        end
        else
        begin
            N := MX;
        end;
        
        //
        // Initialize RefRA/RefCA by random Hermitian matrices,
        // RefRC/RefCC by random matrices
        //
        // RA/CA size is 2Nx2N (four copies of same NxN matrix
        // to test different offsets)
        //
        SetLength(RefRA, 2*N, 2*N);
        SetLength(RefCA, 2*N, 2*N);
        I:=0;
        while I<=N-1 do
        begin
            RefRA[I,I] := 2*RandomReal-1;
            RefCA[I,I] := C_Complex(2*RandomReal-1);
            J:=I+1;
            while J<=N-1 do
            begin
                RefRA[I,J] := 2*RandomReal-1;
                RefCA[I,J].X := 2*RandomReal-1;
                RefCA[I,J].Y := 2*RandomReal-1;
                RefRA[J,I] := RefRA[I,J];
                RefCA[J,I] := Conj(RefCA[I,J]);
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
                RefRA[I+N,J] := RefRA[I,J];
                RefRA[I,J+N] := RefRA[I,J];
                RefRA[I+N,J+N] := RefRA[I,J];
                RefCA[I+N,J] := RefCA[I,J];
                RefCA[I,J+N] := RefCA[I,J];
                RefCA[I+N,J+N] := RefCA[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
        SetLength(RefRC, N, K);
        SetLength(RefCC, N, K);
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=K-1 do
            begin
                RefRC[I,J] := 2*RandomReal-1;
                RefCC[I,J].X := 2*RandomReal-1;
                RefCC[I,J].Y := 2*RandomReal-1;
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // test different types of operations, offsets, and so on...
        //
        // to avoid unnecessary slowdown we don't test ALL possible
        // combinations of operation types. We just generate one random
        // set of parameters and test it.
        //
        SetLength(RA1, 2*N, 2*N);
        SetLength(RA2, 2*N, 2*N);
        SetLength(CA1, 2*N, 2*N);
        SetLength(CA2, 2*N, 2*N);
        SetLength(RC, N, K);
        SetLength(RCT, K, N);
        SetLength(CC, N, K);
        SetLength(CCT, K, N);
        UpperType := RandomInteger(2);
        XOffsI := RandomInteger(2);
        XOffsJ := RandomInteger(2);
        AOffsIType := RandomInteger(2);
        AOffsJType := RandomInteger(2);
        AlphaType := RandomInteger(2);
        BetaType := RandomInteger(2);
        AOffsI := N*AOffsIType;
        AOffsJ := N*AOffsJType;
        Alpha := AlphaType*(2*RandomReal-1);
        Beta := BetaType*(2*RandomReal-1);
        
        //
        // copy A, C (fill unused parts with random garbage)
        //
        I:=0;
        while I<=2*N-1 do
        begin
            J:=0;
            while J<=2*N-1 do
            begin
                if (I>=AOffsI) and (I<AOffsI+N) and (J>=AOffsJ) and (J<AOffsJ+N) then
                begin
                    CA1[I,J] := RefCA[I,J];
                    CA2[I,J] := RefCA[I,J];
                    RA1[I,J] := RefRA[I,J];
                    RA2[I,J] := RefRA[I,J];
                end
                else
                begin
                    CA1[I,J] := C_Complex(RandomReal);
                    CA2[I,J] := CA1[I,J];
                    RA1[I,J] := RandomReal;
                    RA2[I,J] := RA1[I,J];
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=K-1 do
            begin
                if (I>=XOffsI) and (J>=XOffsJ) then
                begin
                    RC[I,J] := RefRC[I,J];
                    RCT[J,I] := RefRC[I,J];
                    CC[I,J] := RefCC[I,J];
                    CCT[J,I] := RefCC[I,J];
                end
                else
                begin
                    RC[I,J] := RandomReal;
                    RCT[J,I] := RC[I,J];
                    CC[I,J] := C_Complex(RandomReal);
                    CCT[J,I] := CCT[J,I];
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Test complex
        // Only one of transform types is selected and tested
        //
        if AP_FP_Greater(RandomReal,0.5) then
        begin
            CMatrixSYRK(N-XOffsI, K-XOffsJ, Alpha, CC, XOffsI, XOffsJ, 0, Beta, CA1, AOffsI, AOffsJ, UpperType=0);
            RefCMatrixSYRK(N-XOffsI, K-XOffsJ, Alpha, CC, XOffsI, XOffsJ, 0, Beta, CA2, AOffsI, AOffsJ, UpperType=0);
        end
        else
        begin
            CMatrixSYRK(N-XOffsI, K-XOffsJ, Alpha, CCT, XOffsJ, XOffsI, 2, Beta, CA1, AOffsI, AOffsJ, UpperType=0);
            RefCMatrixSYRK(N-XOffsI, K-XOffsJ, Alpha, CCT, XOffsJ, XOffsI, 2, Beta, CA2, AOffsI, AOffsJ, UpperType=0);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                Result := Result or AP_FP_Greater(AbsComplex(C_Sub(CA1[I,J],CA2[I,J])),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Test real
        // Only one of transform types is selected and tested
        //
        if AP_FP_Greater(RandomReal,0.5) then
        begin
            RMatrixSYRK(N-XOffsI, K-XOffsJ, Alpha, RC, XOffsI, XOffsJ, 0, Beta, RA1, AOffsI, AOffsJ, UpperType=0);
            RefRMatrixSYRK(N-XOffsI, K-XOffsJ, Alpha, RC, XOffsI, XOffsJ, 0, Beta, RA2, AOffsI, AOffsJ, UpperType=0);
        end
        else
        begin
            RMatrixSYRK(N-XOffsI, K-XOffsJ, Alpha, RCT, XOffsJ, XOffsI, 1, Beta, RA1, AOffsI, AOffsJ, UpperType=0);
            RefRMatrixSYRK(N-XOffsI, K-XOffsJ, Alpha, RCT, XOffsJ, XOffsI, 1, Beta, RA2, AOffsI, AOffsJ, UpperType=0);
        end;
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                Result := Result or AP_FP_Greater(AbsReal(RA1[I,J]-RA2[I,J]),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
        Inc(MX);
    end;
end;


(*************************************************************************
GEMM tests

Returns False for passed test, True - for failed
*************************************************************************)
function TestGEMM(MinN : AlglibInteger; MaxN : AlglibInteger):Boolean;
var
    M : AlglibInteger;
    N : AlglibInteger;
    K : AlglibInteger;
    MX : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    AOffsI : AlglibInteger;
    AOffsJ : AlglibInteger;
    AOpType : AlglibInteger;
    AOpTypeR : AlglibInteger;
    BOffsI : AlglibInteger;
    BOffsJ : AlglibInteger;
    BOpType : AlglibInteger;
    BOpTypeR : AlglibInteger;
    COffsI : AlglibInteger;
    COffsJ : AlglibInteger;
    RefRA : TReal2DArray;
    RefRB : TReal2DArray;
    RefRC : TReal2DArray;
    RefCA : TComplex2DArray;
    RefCB : TComplex2DArray;
    RefCC : TComplex2DArray;
    AlphaR : Double;
    BetaR : Double;
    AlphaC : Complex;
    BetaC : Complex;
    RC1 : TReal2DArray;
    RC2 : TReal2DArray;
    CC1 : TComplex2DArray;
    CC2 : TComplex2DArray;
    Threshold : Double;
begin
    Threshold := MaxN*100*MachineEpsilon;
    Result := False;
    MX:=MinN;
    while MX<=MaxN do
    begin
        
        //
        // Select random M/N/K in [1,MX] such that max(M,N,K)=MX
        //
        M := 1+RandomInteger(MX);
        N := 1+RandomInteger(MX);
        K := 1+RandomInteger(MX);
        I := RandomInteger(3);
        if I=0 then
        begin
            M := MX;
        end;
        if I=1 then
        begin
            N := MX;
        end;
        if I=2 then
        begin
            K := MX;
        end;
        
        //
        // Initialize A/B/C by random matrices with size (MaxN+1)*(MaxN+1)
        //
        SetLength(RefRA, MaxN+1, MaxN+1);
        SetLength(RefRB, MaxN+1, MaxN+1);
        SetLength(RefRC, MaxN+1, MaxN+1);
        SetLength(RefCA, MaxN+1, MaxN+1);
        SetLength(RefCB, MaxN+1, MaxN+1);
        SetLength(RefCC, MaxN+1, MaxN+1);
        I:=0;
        while I<=MaxN do
        begin
            J:=0;
            while J<=MaxN do
            begin
                RefRA[I,J] := 2*RandomReal-1;
                RefRB[I,J] := 2*RandomReal-1;
                RefRC[I,J] := 2*RandomReal-1;
                RefCA[I,J].X := 2*RandomReal-1;
                RefCA[I,J].Y := 2*RandomReal-1;
                RefCB[I,J].X := 2*RandomReal-1;
                RefCB[I,J].Y := 2*RandomReal-1;
                RefCC[I,J].X := 2*RandomReal-1;
                RefCC[I,J].Y := 2*RandomReal-1;
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // test different types of operations, offsets, and so on...
        //
        // to avoid unnecessary slowdown we don't test ALL possible
        // combinations of operation types. We just generate one random
        // set of parameters and test it.
        //
        SetLength(RC1, MaxN+1, MaxN+1);
        SetLength(RC2, MaxN+1, MaxN+1);
        SetLength(CC1, MaxN+1, MaxN+1);
        SetLength(CC2, MaxN+1, MaxN+1);
        AOffsI := RandomInteger(2);
        AOffsJ := RandomInteger(2);
        AOpType := RandomInteger(3);
        AOpTypeR := RandomInteger(2);
        BOffsI := RandomInteger(2);
        BOffsJ := RandomInteger(2);
        BOpType := RandomInteger(3);
        BOpTypeR := RandomInteger(2);
        COffsI := RandomInteger(2);
        COffsJ := RandomInteger(2);
        AlphaR := RandomInteger(2)*(2*RandomReal-1);
        BetaR := RandomInteger(2)*(2*RandomReal-1);
        if AP_FP_Greater(RandomReal,0.5) then
        begin
            AlphaC.X := 2*RandomReal-1;
            AlphaC.Y := 2*RandomReal-1;
        end
        else
        begin
            AlphaC := C_Complex(0);
        end;
        if AP_FP_Greater(RandomReal,0.5) then
        begin
            BetaC.X := 2*RandomReal-1;
            BetaC.Y := 2*RandomReal-1;
        end
        else
        begin
            BetaC := C_Complex(0);
        end;
        
        //
        // copy C
        //
        I:=0;
        while I<=MaxN do
        begin
            J:=0;
            while J<=MaxN do
            begin
                RC1[I,J] := RefRC[I,J];
                RC2[I,J] := RefRC[I,J];
                CC1[I,J] := RefCC[I,J];
                CC2[I,J] := RefCC[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Test complex
        //
        CMatrixGEMM(M, N, K, AlphaC, RefCA, AOffsI, AOffsJ, AOpType, RefCB, BOffsI, BOffsJ, BOpType, BetaC, CC1, COffsI, COffsJ);
        RefCMatrixGEMM(M, N, K, AlphaC, RefCA, AOffsI, AOffsJ, AOpType, RefCB, BOffsI, BOffsJ, BOpType, BetaC, CC2, COffsI, COffsJ);
        I:=0;
        while I<=MaxN do
        begin
            J:=0;
            while J<=MaxN do
            begin
                Result := Result or AP_FP_Greater(AbsComplex(C_Sub(CC1[I,J],CC2[I,J])),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // Test real
        //
        RMatrixGEMM(M, N, K, AlphaR, RefRA, AOffsI, AOffsJ, AOpTypeR, RefRB, BOffsI, BOffsJ, BOpTypeR, BetaR, RC1, COffsI, COffsJ);
        RefRMatrixGEMM(M, N, K, AlphaR, RefRA, AOffsI, AOffsJ, AOpTypeR, RefRB, BOffsI, BOffsJ, BOpTypeR, BetaR, RC2, COffsI, COffsJ);
        I:=0;
        while I<=MaxN do
        begin
            J:=0;
            while J<=MaxN do
            begin
                Result := Result or AP_FP_Greater(AbsReal(RC1[I,J]-RC2[I,J]),Threshold);
                Inc(J);
            end;
            Inc(I);
        end;
        Inc(MX);
    end;
end;


(*************************************************************************
transpose tests

Returns False for passed test, True - for failed
*************************************************************************)
function TestTrans(MinN : AlglibInteger; MaxN : AlglibInteger):Boolean;
var
    M : AlglibInteger;
    N : AlglibInteger;
    MX : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    AOffsI : AlglibInteger;
    AOffsJ : AlglibInteger;
    BOffsI : AlglibInteger;
    BOffsJ : AlglibInteger;
    V1 : Double;
    V2 : Double;
    Threshold : Double;
    RefRA : TReal2DArray;
    RefRB : TReal2DArray;
    RefCA : TComplex2DArray;
    RefCB : TComplex2DArray;
begin
    Result := False;
    Threshold := 1000*MachineEpsilon;
    MX:=MinN;
    while MX<=MaxN do
    begin
        
        //
        // Select random M/N in [1,MX] such that max(M,N)=MX
        // Generate random V1 and V2 which are used to fill
        // RefRB/RefCB with control values.
        //
        M := 1+RandomInteger(MX);
        N := 1+RandomInteger(MX);
        if RandomInteger(2)=0 then
        begin
            M := MX;
        end
        else
        begin
            N := MX;
        end;
        V1 := RandomReal;
        V2 := RandomReal;
        
        //
        // Initialize A by random matrix with size (MaxN+1)*(MaxN+1)
        // Fill B with control values
        //
        SetLength(RefRA, MaxN+1, MaxN+1);
        SetLength(RefRB, MaxN+1, MaxN+1);
        SetLength(RefCA, MaxN+1, MaxN+1);
        SetLength(RefCB, MaxN+1, MaxN+1);
        I:=0;
        while I<=MaxN do
        begin
            J:=0;
            while J<=MaxN do
            begin
                RefRA[I,J] := 2*RandomReal-1;
                RefCA[I,J].X := 2*RandomReal-1;
                RefCA[I,J].Y := 2*RandomReal-1;
                RefRB[I,J] := I*V1+J*V2;
                RefCB[I,J] := C_Complex(I*V1+J*V2);
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // test different offsets (zero or one)
        //
        // to avoid unnecessary slowdown we don't test ALL possible
        // combinations of operation types. We just generate one random
        // set of parameters and test it.
        //
        AOffsI := RandomInteger(2);
        AOffsJ := RandomInteger(2);
        BOffsI := RandomInteger(2);
        BOffsJ := RandomInteger(2);
        RMatrixTranspose(M, N, RefRA, AOffsI, AOffsJ, RefRB, BOffsI, BOffsJ);
        I:=0;
        while I<=MaxN do
        begin
            J:=0;
            while J<=MaxN do
            begin
                if (I<BOffsI) or (I>=BOffsI+N) or (J<BOffsJ) or (J>=BOffsJ+M) then
                begin
                    Result := Result or AP_FP_Greater(AbsReal(RefRB[I,J]-(V1*I+V2*J)),Threshold);
                end
                else
                begin
                    Result := Result or AP_FP_Greater(AbsReal(RefRB[I,J]-RefRA[AOffsI+J-BOffsJ,AOffsJ+I-BOffsI]),Threshold);
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        CMatrixTranspose(M, N, RefCA, AOffsI, AOffsJ, RefCB, BOffsI, BOffsJ);
        I:=0;
        while I<=MaxN do
        begin
            J:=0;
            while J<=MaxN do
            begin
                if (I<BOffsI) or (I>=BOffsI+N) or (J<BOffsJ) or (J>=BOffsJ+M) then
                begin
                    Result := Result or AP_FP_Greater(AbsComplex(C_SubR(RefCB[I,J],V1*I+V2*J)),Threshold);
                end
                else
                begin
                    Result := Result or AP_FP_Greater(AbsComplex(C_Sub(RefCB[I,J],RefCA[AOffsI+J-BOffsJ,AOffsJ+I-BOffsI])),Threshold);
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        Inc(MX);
    end;
end;


(*************************************************************************
rank-1tests

Returns False for passed test, True - for failed
*************************************************************************)
function TestRANK1(MinN : AlglibInteger; MaxN : AlglibInteger):Boolean;
var
    M : AlglibInteger;
    N : AlglibInteger;
    MX : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    AOffsI : AlglibInteger;
    AOffsJ : AlglibInteger;
    UOffs : AlglibInteger;
    VOffs : AlglibInteger;
    Threshold : Double;
    RefRA : TReal2DArray;
    RefRB : TReal2DArray;
    RefCA : TComplex2DArray;
    RefCB : TComplex2DArray;
    RU : TReal1DArray;
    RV : TReal1DArray;
    CU : TComplex1DArray;
    CV : TComplex1DArray;
begin
    Result := False;
    Threshold := 1000*MachineEpsilon;
    MX:=MinN;
    while MX<=MaxN do
    begin
        
        //
        // Select random M/N in [1,MX] such that max(M,N)=MX
        //
        M := 1+RandomInteger(MX);
        N := 1+RandomInteger(MX);
        if RandomInteger(2)=0 then
        begin
            M := MX;
        end
        else
        begin
            N := MX;
        end;
        
        //
        // Initialize A by random matrix with size (MaxN+1)*(MaxN+1)
        // Fill B with control values
        //
        SetLength(RefRA, MaxN+MaxN, MaxN+MaxN);
        SetLength(RefRB, MaxN+MaxN, MaxN+MaxN);
        SetLength(RefCA, MaxN+MaxN, MaxN+MaxN);
        SetLength(RefCB, MaxN+MaxN, MaxN+MaxN);
        I:=0;
        while I<=2*MaxN-1 do
        begin
            J:=0;
            while J<=2*MaxN-1 do
            begin
                RefRA[I,J] := 2*RandomReal-1;
                RefCA[I,J].X := 2*RandomReal-1;
                RefCA[I,J].Y := 2*RandomReal-1;
                RefRB[I,J] := RefRA[I,J];
                RefCB[I,J] := RefCA[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
        SetLength(RU, 2*M);
        SetLength(CU, 2*M);
        I:=0;
        while I<=2*M-1 do
        begin
            RU[I] := 2*RandomReal-1;
            CU[I].X := 2*RandomReal-1;
            CU[I].Y := 2*RandomReal-1;
            Inc(I);
        end;
        SetLength(RV, 2*N);
        SetLength(CV, 2*N);
        I:=0;
        while I<=2*N-1 do
        begin
            RV[I] := 2*RandomReal-1;
            CV[I].X := 2*RandomReal-1;
            CV[I].Y := 2*RandomReal-1;
            Inc(I);
        end;
        
        //
        // test different offsets (zero or one)
        //
        // to avoid unnecessary slowdown we don't test ALL possible
        // combinations of operation types. We just generate one random
        // set of parameters and test it.
        //
        AOffsI := RandomInteger(MaxN);
        AOffsJ := RandomInteger(MaxN);
        UOffs := RandomInteger(M);
        VOffs := RandomInteger(N);
        CMatrixRank1(M, N, RefCA, AOffsI, AOffsJ, CU, UOffs, CV, VOffs);
        I:=0;
        while I<=2*MaxN-1 do
        begin
            J:=0;
            while J<=2*MaxN-1 do
            begin
                if (I<AOffsI) or (I>=AOffsI+M) or (J<AOffsJ) or (J>=AOffsJ+N) then
                begin
                    Result := Result or AP_FP_Greater(AbsComplex(C_Sub(RefCA[I,J],RefCB[I,J])),Threshold);
                end
                else
                begin
                    Result := Result or AP_FP_Greater(AbsComplex(C_Sub(RefCA[I,J],C_Add(RefCB[I,J],C_Mul(CU[I-AOffsI+UOffs],CV[J-AOffsJ+VOffs])))),Threshold);
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        RMatrixRank1(M, N, RefRA, AOffsI, AOffsJ, RU, UOffs, RV, VOffs);
        I:=0;
        while I<=2*MaxN-1 do
        begin
            J:=0;
            while J<=2*MaxN-1 do
            begin
                if (I<AOffsI) or (I>=AOffsI+M) or (J<AOffsJ) or (J>=AOffsJ+N) then
                begin
                    Result := Result or AP_FP_Greater(AbsReal(RefRA[I,J]-RefRB[I,J]),Threshold);
                end
                else
                begin
                    Result := Result or AP_FP_Greater(AbsReal(RefRA[I,J]-(RefRB[I,J]+RU[I-AOffsI+UOffs]*RV[J-AOffsJ+VOffs])),Threshold);
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        Inc(MX);
    end;
end;


(*************************************************************************
MV tests

Returns False for passed test, True - for failed
*************************************************************************)
function TestMV(MinN : AlglibInteger; MaxN : AlglibInteger):Boolean;
var
    M : AlglibInteger;
    N : AlglibInteger;
    MX : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    AOffsI : AlglibInteger;
    AOffsJ : AlglibInteger;
    XOffs : AlglibInteger;
    YOffs : AlglibInteger;
    OpCA : AlglibInteger;
    OpRA : AlglibInteger;
    Threshold : Double;
    RV1 : Double;
    RV2 : Double;
    CV1 : Complex;
    CV2 : Complex;
    RefRA : TReal2DArray;
    RefCA : TComplex2DArray;
    RX : TReal1DArray;
    RY : TReal1DArray;
    CX : TComplex1DArray;
    CY : TComplex1DArray;
    i_ : AlglibInteger;
    i1_ : AlglibInteger;
begin
    Result := False;
    Threshold := 1000*MachineEpsilon;
    MX:=MinN;
    while MX<=MaxN do
    begin
        
        //
        // Select random M/N in [1,MX] such that max(M,N)=MX
        //
        M := 1+RandomInteger(MX);
        N := 1+RandomInteger(MX);
        if RandomInteger(2)=0 then
        begin
            M := MX;
        end
        else
        begin
            N := MX;
        end;
        
        //
        // Initialize A by random matrix with size (MaxN+MaxN)*(MaxN+MaxN)
        // Initialize X by random vector with size (MaxN+MaxN)
        // Fill Y by control values
        //
        SetLength(RefRA, MaxN+MaxN, MaxN+MaxN);
        SetLength(RefCA, MaxN+MaxN, MaxN+MaxN);
        I:=0;
        while I<=2*MaxN-1 do
        begin
            J:=0;
            while J<=2*MaxN-1 do
            begin
                RefRA[I,J] := 2*RandomReal-1;
                RefCA[I,J].X := 2*RandomReal-1;
                RefCA[I,J].Y := 2*RandomReal-1;
                Inc(J);
            end;
            Inc(I);
        end;
        SetLength(RX, 2*MaxN);
        SetLength(CX, 2*MaxN);
        SetLength(RY, 2*MaxN);
        SetLength(CY, 2*MaxN);
        I:=0;
        while I<=2*MaxN-1 do
        begin
            RX[I] := 2*RandomReal-1;
            CX[I].X := 2*RandomReal-1;
            CX[I].Y := 2*RandomReal-1;
            RY[I] := I;
            CY[I] := C_Complex(I);
            Inc(I);
        end;
        
        //
        // test different offsets (zero or one)
        //
        // to avoid unnecessary slowdown we don't test ALL possible
        // combinations of operation types. We just generate one random
        // set of parameters and test it.
        //
        AOffsI := RandomInteger(MaxN);
        AOffsJ := RandomInteger(MaxN);
        XOffs := RandomInteger(MaxN);
        YOffs := RandomInteger(MaxN);
        OpCA := RandomInteger(3);
        OpRA := RandomInteger(2);
        CMatrixMV(M, N, RefCA, AOffsI, AOffsJ, OpCA, CX, XOffs, CY, YOffs);
        I:=0;
        while I<=2*MaxN-1 do
        begin
            if (I<YOffs) or (I>=YOffs+M) then
            begin
                Result := Result or C_NotEqualR(CY[I],I);
            end
            else
            begin
                CV1 := CY[I];
                if OpCA=0 then
                begin
                    i1_ := (XOffs)-(AOffsJ);
                    CV2 := C_Complex(0.0);
                    for i_ := AOffsJ to AOffsJ+N-1 do
                    begin
                        CV2 := C_Add(CV2,C_Mul(RefCA[AOffsI+I-YOffs,i_],CX[i_+i1_]));
                    end;
                end;
                if OpCA=1 then
                begin
                    i1_ := (XOffs)-(AOffsI);
                    CV2 := C_Complex(0.0);
                    for i_ := AOffsI to AOffsI+N-1 do
                    begin
                        CV2 := C_Add(CV2,C_Mul(RefCA[i_,AOffsJ+I-YOffs],CX[i_+i1_]));
                    end;
                end;
                if OpCA=2 then
                begin
                    i1_ := (XOffs)-(AOffsI);
                    CV2 := C_Complex(0.0);
                    for i_ := AOffsI to AOffsI+N-1 do
                    begin
                        CV2 := C_Add(CV2,C_Mul(Conj(RefCA[i_,AOffsJ+I-YOffs]),CX[i_+i1_]));
                    end;
                end;
                Result := Result or AP_FP_Greater(AbsComplex(C_Sub(CV1,CV2)),Threshold);
            end;
            Inc(I);
        end;
        RMatrixMV(M, N, RefRA, AOffsI, AOffsJ, OpRA, RX, XOffs, RY, YOffs);
        I:=0;
        while I<=2*MaxN-1 do
        begin
            if (I<YOffs) or (I>=YOffs+M) then
            begin
                Result := Result or AP_FP_Neq(RY[I],I);
            end
            else
            begin
                RV1 := RY[I];
                if OpRA=0 then
                begin
                    RV2 := APVDotProduct(@RefRA[AOffsI+I-YOffs][0], AOffsJ, AOffsJ+N-1, @RX[0], XOffs, XOffs+N-1);
                end;
                if OpRA=1 then
                begin
                    i1_ := (XOffs)-(AOffsI);
                    RV2 := 0.0;
                    for i_ := AOffsI to AOffsI+N-1 do
                    begin
                        RV2 := RV2 + RefRA[i_,AOffsJ+I-YOffs]*RX[i_+i1_];
                    end;
                end;
                Result := Result or AP_FP_Greater(AbsReal(RV1-RV2),Threshold);
            end;
            Inc(I);
        end;
        Inc(MX);
    end;
end;


(*************************************************************************
COPY tests

Returns False for passed test, True - for failed
*************************************************************************)
function TestCopy(MinN : AlglibInteger; MaxN : AlglibInteger):Boolean;
var
    M : AlglibInteger;
    N : AlglibInteger;
    MX : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    AOffsI : AlglibInteger;
    AOffsJ : AlglibInteger;
    BOffsI : AlglibInteger;
    BOffsJ : AlglibInteger;
    Threshold : Double;
    RV1 : Double;
    RV2 : Double;
    CV1 : Complex;
    CV2 : Complex;
    RA : TReal2DArray;
    RB : TReal2DArray;
    CA : TComplex2DArray;
    CB : TComplex2DArray;
begin
    Result := False;
    Threshold := 1000*MachineEpsilon;
    MX:=MinN;
    while MX<=MaxN do
    begin
        
        //
        // Select random M/N in [1,MX] such that max(M,N)=MX
        //
        M := 1+RandomInteger(MX);
        N := 1+RandomInteger(MX);
        if RandomInteger(2)=0 then
        begin
            M := MX;
        end
        else
        begin
            N := MX;
        end;
        
        //
        // Initialize A by random matrix with size (MaxN+MaxN)*(MaxN+MaxN)
        // Initialize X by random vector with size (MaxN+MaxN)
        // Fill Y by control values
        //
        SetLength(RA, MaxN+MaxN, MaxN+MaxN);
        SetLength(CA, MaxN+MaxN, MaxN+MaxN);
        SetLength(RB, MaxN+MaxN, MaxN+MaxN);
        SetLength(CB, MaxN+MaxN, MaxN+MaxN);
        I:=0;
        while I<=2*MaxN-1 do
        begin
            J:=0;
            while J<=2*MaxN-1 do
            begin
                RA[I,J] := 2*RandomReal-1;
                CA[I,J].X := 2*RandomReal-1;
                CA[I,J].Y := 2*RandomReal-1;
                RB[I,J] := 1+2*I+3*J;
                CB[I,J] := C_Complex(1+2*I+3*J);
                Inc(J);
            end;
            Inc(I);
        end;
        
        //
        // test different offsets (zero or one)
        //
        // to avoid unnecessary slowdown we don't test ALL possible
        // combinations of operation types. We just generate one random
        // set of parameters and test it.
        //
        AOffsI := RandomInteger(MaxN);
        AOffsJ := RandomInteger(MaxN);
        BOffsI := RandomInteger(MaxN);
        BOffsJ := RandomInteger(MaxN);
        CMatrixCopy(M, N, CA, AOffsI, AOffsJ, CB, BOffsI, BOffsJ);
        I:=0;
        while I<=2*MaxN-1 do
        begin
            J:=0;
            while J<=2*MaxN-1 do
            begin
                if (I<BOffsI) or (I>=BOffsI+M) or (J<BOffsJ) or (J>=BOffsJ+N) then
                begin
                    Result := Result or C_NotEqualR(CB[I,J],1+2*I+3*J);
                end
                else
                begin
                    Result := Result or AP_FP_Greater(AbsComplex(C_Sub(CA[AOffsI+I-BOffsI,AOffsJ+J-BOffsJ],CB[I,J])),Threshold);
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        RMatrixCopy(M, N, RA, AOffsI, AOffsJ, RB, BOffsI, BOffsJ);
        I:=0;
        while I<=2*MaxN-1 do
        begin
            J:=0;
            while J<=2*MaxN-1 do
            begin
                if (I<BOffsI) or (I>=BOffsI+M) or (J<BOffsJ) or (J>=BOffsJ+N) then
                begin
                    Result := Result or AP_FP_Neq(RB[I,J],1+2*I+3*J);
                end
                else
                begin
                    Result := Result or AP_FP_Greater(AbsReal(RA[AOffsI+I-BOffsI,AOffsJ+J-BOffsJ]-RB[I,J]),Threshold);
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        Inc(MX);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testablasunit_test_silent():Boolean;
begin
    Result := TestABLAS(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testablasunit_test():Boolean;
begin
    Result := TestABLAS(False);
end;


end.