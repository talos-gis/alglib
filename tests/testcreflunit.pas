unit testcreflunit;
interface
uses Math, Sysutils, Ap, creflections;

function TestCRefl(Silent : Boolean):Boolean;
function testcreflunit_test_silent():Boolean;
function testcreflunit_test():Boolean;

implementation

function TestCRefl(Silent : Boolean):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    N : AlglibInteger;
    M : AlglibInteger;
    MaxMN : AlglibInteger;
    X : TComplex1DArray;
    V : TComplex1DArray;
    WORK : TComplex1DArray;
    H : TComplex2DArray;
    A : TComplex2DArray;
    B : TComplex2DArray;
    C : TComplex2DArray;
    Tmp : Complex;
    Beta : Complex;
    Tau : Complex;
    Err : Double;
    MER : Double;
    MEL : Double;
    MEG : Double;
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    WasErrors : Boolean;
    Threshold : Double;
    i_ : AlglibInteger;
begin
    Threshold := 1000*MachineEpsilon;
    PassCount := 1000;
    MER := 0;
    MEL := 0;
    MEG := 0;
    Pass:=1;
    while Pass<=PassCount do
    begin
        
        //
        // Task
        //
        N := 1+RandomInteger(10);
        M := 1+RandomInteger(10);
        MaxMN := Max(M, N);
        
        //
        // Initialize
        //
        SetLength(X, MaxMN+1);
        SetLength(V, MaxMN+1);
        SetLength(WORK, MaxMN+1);
        SetLength(H, MaxMN+1, MaxMN+1);
        SetLength(A, MaxMN+1, MaxMN+1);
        SetLength(B, MaxMN+1, MaxMN+1);
        SetLength(C, MaxMN+1, MaxMN+1);
        
        //
        // GenerateReflection
        //
        I:=1;
        while I<=N do
        begin
            X[I].X := 2*RandomReal-1;
            X[I].Y := 2*RandomReal-1;
            V[I] := X[I];
            Inc(I);
        end;
        ComplexGenerateReflection(V, N, Tau);
        Beta := V[1];
        V[1] := C_Complex(1);
        I:=1;
        while I<=N do
        begin
            J:=1;
            while J<=N do
            begin
                if I=J then
                begin
                    H[I,J] := C_RSub(1,C_Mul(C_Mul(Tau,V[I]),Conj(V[J])));
                end
                else
                begin
                    H[I,J] := C_Opposite(C_Mul(C_Mul(Tau,V[I]),Conj(V[J])));
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        Err := 0;
        I:=1;
        while I<=N do
        begin
            Tmp := C_Complex(0.0);
            for i_ := 1 to N do
            begin
                Tmp := C_Add(Tmp,C_Mul(Conj(H[i_,I]),X[i_]));
            end;
            if I=1 then
            begin
                Err := Max(Err, AbsComplex(C_Sub(Tmp,Beta)));
            end
            else
            begin
                Err := Max(Err, AbsComplex(Tmp));
            end;
            Inc(I);
        end;
        Err := Max(Err, AbsReal(Beta.Y));
        MEG := Max(MEG, Err);
        
        //
        // ApplyReflectionFromTheLeft
        //
        I:=1;
        while I<=M do
        begin
            X[I].X := 2*RandomReal-1;
            X[I].Y := 2*RandomReal-1;
            V[I] := X[I];
            Inc(I);
        end;
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=N do
            begin
                A[I,J].X := 2*RandomReal-1;
                A[I,J].Y := 2*RandomReal-1;
                B[I,J] := A[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
        ComplexGenerateReflection(V, M, Tau);
        Beta := V[1];
        V[1] := C_Complex(1);
        ComplexApplyReflectionFromTheLeft(B, Tau, V, 1, M, 1, N, WORK);
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=M do
            begin
                if I=J then
                begin
                    H[I,J] := C_RSub(1,C_Mul(C_Mul(Tau,V[I]),Conj(V[J])));
                end
                else
                begin
                    H[I,J] := C_Opposite(C_Mul(C_Mul(Tau,V[I]),Conj(V[J])));
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=N do
            begin
                Tmp := C_Complex(0.0);
                for i_ := 1 to M do
                begin
                    Tmp := C_Add(Tmp,C_Mul(H[I,i_],A[i_,J]));
                end;
                C[I,J] := Tmp;
                Inc(J);
            end;
            Inc(I);
        end;
        Err := 0;
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=N do
            begin
                Err := Max(Err, AbsComplex(C_Sub(B[I,J],C[I,J])));
                Inc(J);
            end;
            Inc(I);
        end;
        MEL := Max(MEL, Err);
        
        //
        // ApplyReflectionFromTheRight
        //
        I:=1;
        while I<=N do
        begin
            X[I] := C_Complex(2*RandomReal-1);
            V[I] := X[I];
            Inc(I);
        end;
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=N do
            begin
                A[I,J] := C_Complex(2*RandomReal-1);
                B[I,J] := A[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
        ComplexGenerateReflection(V, N, Tau);
        Beta := V[1];
        V[1] := C_Complex(1);
        ComplexApplyReflectionFromTheRight(B, Tau, V, 1, M, 1, N, WORK);
        I:=1;
        while I<=N do
        begin
            J:=1;
            while J<=N do
            begin
                if I=J then
                begin
                    H[I,J] := C_RSub(1,C_Mul(C_Mul(Tau,V[I]),Conj(V[J])));
                end
                else
                begin
                    H[I,J] := C_Opposite(C_Mul(C_Mul(Tau,V[I]),Conj(V[J])));
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=N do
            begin
                Tmp := C_Complex(0.0);
                for i_ := 1 to N do
                begin
                    Tmp := C_Add(Tmp,C_Mul(A[I,i_],H[i_,J]));
                end;
                C[I,J] := Tmp;
                Inc(J);
            end;
            Inc(I);
        end;
        Err := 0;
        I:=1;
        while I<=M do
        begin
            J:=1;
            while J<=N do
            begin
                Err := Max(Err, AbsComplex(C_Sub(B[I,J],C[I,J])));
                Inc(J);
            end;
            Inc(I);
        end;
        MER := Max(MER, Err);
        Inc(Pass);
    end;
    
    //
    // Overflow crash test
    //
    SetLength(X, 10+1);
    SetLength(V, 10+1);
    I:=1;
    while I<=10 do
    begin
        V[I] := C_Complex(MaxRealNumber*0.01*(2*RandomReal-1));
        Inc(I);
    end;
    ComplexGenerateReflection(V, 10, Tau);
    
    //
    // report
    //
    WasErrors := AP_FP_Greater(MEG,Threshold) or AP_FP_Greater(MEL,Threshold) or AP_FP_Greater(MER,Threshold);
    if  not Silent then
    begin
        Write(Format('TESTING COMPLEX REFLECTIONS'#13#10'',[]));
        Write(Format('Generate error:                          %5.4e'#13#10'',[
            MEG]));
        Write(Format('Apply(L) error:                          %5.4e'#13#10'',[
            MEL]));
        Write(Format('Apply(R) error:                          %5.4e'#13#10'',[
            MER]));
        Write(Format('Threshold:                               %5.4e'#13#10'',[
            Threshold]));
        Write(Format('Overflow crash test:                     PASSED'#13#10'',[]));
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
Silent unit test
*************************************************************************)
function testcreflunit_test_silent():Boolean;
begin
    Result := TestCRefl(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testcreflunit_test():Boolean;
begin
    Result := TestCRefl(False);
end;


end.