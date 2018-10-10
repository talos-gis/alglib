unit testreflectionsunit;
interface
uses Math, Sysutils, Ap, reflections;

function TestReflections(Silent : Boolean):Boolean;
function testreflectionsunit_test_silent():Boolean;
function testreflectionsunit_test():Boolean;

implementation

function TestReflections(Silent : Boolean):Boolean;
var
    I : AlglibInteger;
    J : AlglibInteger;
    N : AlglibInteger;
    M : AlglibInteger;
    MaxMN : AlglibInteger;
    X : TReal1DArray;
    V : TReal1DArray;
    WORK : TReal1DArray;
    H : TReal2DArray;
    A : TReal2DArray;
    B : TReal2DArray;
    C : TReal2DArray;
    Tmp : Double;
    Beta : Double;
    Tau : Double;
    Err : Double;
    MER : Double;
    MEL : Double;
    MEG : Double;
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    Threshold : Double;
    TaskType : AlglibInteger;
    XScale : Double;
    i_ : AlglibInteger;
begin
    PassCount := 10;
    Threshold := 100*MachineEpsilon;
    MER := 0;
    MEL := 0;
    MEG := 0;
    Pass:=1;
    while Pass<=PassCount do
    begin
        N:=1;
        while N<=10 do
        begin
            M:=1;
            while M<=10 do
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
                // GenerateReflection, three tasks are possible:
                // * random X
                // * zero X
                // * non-zero X[1], all other are zeros
                // * random X, near underflow scale
                // * random X, near overflow scale
                //
                TaskType:=0;
                while TaskType<=4 do
                begin
                    XScale := 1;
                    if TaskType=0 then
                    begin
                        I:=1;
                        while I<=N do
                        begin
                            X[I] := 2*RandomReal-1;
                            Inc(I);
                        end;
                    end;
                    if TaskType=1 then
                    begin
                        I:=1;
                        while I<=N do
                        begin
                            X[I] := 0;
                            Inc(I);
                        end;
                    end;
                    if TaskType=2 then
                    begin
                        X[1] := 2*RandomReal-1;
                        I:=2;
                        while I<=N do
                        begin
                            X[I] := 0;
                            Inc(I);
                        end;
                    end;
                    if TaskType=3 then
                    begin
                        I:=1;
                        while I<=N do
                        begin
                            X[I] := (RandomInteger(21)-10)*MinRealNumber;
                            Inc(I);
                        end;
                        XScale := 10*MinRealNumber;
                    end;
                    if TaskType=4 then
                    begin
                        I:=1;
                        while I<=N do
                        begin
                            X[I] := (2*RandomReal-1)*MaxRealNumber;
                            Inc(I);
                        end;
                        XScale := MaxRealNumber;
                    end;
                    APVMove(@V[0], 1, N, @X[0], 1, N);
                    GenerateReflection(V, N, Tau);
                    Beta := V[1];
                    V[1] := 1;
                    I:=1;
                    while I<=N do
                    begin
                        J:=1;
                        while J<=N do
                        begin
                            if I=J then
                            begin
                                H[I,J] := 1-Tau*V[I]*V[J];
                            end
                            else
                            begin
                                H[I,J] := -Tau*V[I]*V[J];
                            end;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    Err := 0;
                    I:=1;
                    while I<=N do
                    begin
                        Tmp := APVDotProduct(@H[I][0], 1, N, @X[0], 1, N);
                        if I=1 then
                        begin
                            Err := Max(Err, AbsReal(Tmp-Beta));
                        end
                        else
                        begin
                            Err := Max(Err, AbsReal(Tmp));
                        end;
                        Inc(I);
                    end;
                    MEG := Max(MEG, Err/XScale);
                    Inc(TaskType);
                end;
                
                //
                // ApplyReflectionFromTheLeft
                //
                I:=1;
                while I<=M do
                begin
                    X[I] := 2*RandomReal-1;
                    V[I] := X[I];
                    Inc(I);
                end;
                I:=1;
                while I<=M do
                begin
                    J:=1;
                    while J<=N do
                    begin
                        A[I,J] := 2*RandomReal-1;
                        B[I,J] := A[I,J];
                        Inc(J);
                    end;
                    Inc(I);
                end;
                GenerateReflection(V, M, Tau);
                Beta := V[1];
                V[1] := 1;
                ApplyReflectionFromTheLeft(B, Tau, V, 1, M, 1, N, WORK);
                I:=1;
                while I<=M do
                begin
                    J:=1;
                    while J<=M do
                    begin
                        if I=J then
                        begin
                            H[I,J] := 1-Tau*V[I]*V[J];
                        end
                        else
                        begin
                            H[I,J] := -Tau*V[I]*V[J];
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
                        Tmp := 0.0;
                        for i_ := 1 to M do
                        begin
                            Tmp := Tmp + H[I,i_]*A[i_,J];
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
                        Err := Max(Err, AbsReal(B[I,J]-C[I,J]));
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
                    X[I] := 2*RandomReal-1;
                    V[I] := X[I];
                    Inc(I);
                end;
                I:=1;
                while I<=M do
                begin
                    J:=1;
                    while J<=N do
                    begin
                        A[I,J] := 2*RandomReal-1;
                        B[I,J] := A[I,J];
                        Inc(J);
                    end;
                    Inc(I);
                end;
                GenerateReflection(V, N, Tau);
                Beta := V[1];
                V[1] := 1;
                ApplyReflectionFromTheRight(B, Tau, V, 1, M, 1, N, WORK);
                I:=1;
                while I<=N do
                begin
                    J:=1;
                    while J<=N do
                    begin
                        if I=J then
                        begin
                            H[I,J] := 1-Tau*V[I]*V[J];
                        end
                        else
                        begin
                            H[I,J] := -Tau*V[I]*V[J];
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
                        Tmp := 0.0;
                        for i_ := 1 to N do
                        begin
                            Tmp := Tmp + A[I,i_]*H[i_,J];
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
                        Err := Max(Err, AbsReal(B[I,J]-C[I,J]));
                        Inc(J);
                    end;
                    Inc(I);
                end;
                MER := Max(MER, Err);
                Inc(M);
            end;
            Inc(N);
        end;
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
        V[I] := MaxRealNumber*0.01*(2*RandomReal-1);
        Inc(I);
    end;
    GenerateReflection(V, 10, Tau);
    Result := AP_FP_Less_Eq(MEG,Threshold) and AP_FP_Less_Eq(MEL,Threshold) and AP_FP_Less_Eq(MER,Threshold);
    if  not Silent then
    begin
        Write(Format('TESTING REFLECTIONS'#13#10'',[]));
        Write(Format('Pass count is %0d'#13#10'',[
            PassCount]));
        Write(Format('Generate     absolute error is       %5.4e'#13#10'',[
            MEG]));
        Write(Format('Apply(Left)  absolute error is       %5.4e'#13#10'',[
            MEL]));
        Write(Format('Apply(Right) absolute error is       %5.4e'#13#10'',[
            MER]));
        Write(Format('Overflow crash test passed'#13#10'',[]));
        if Result then
        begin
            Write(Format('TEST PASSED'#13#10'',[]));
        end
        else
        begin
            Write(Format('TEST FAILED'#13#10'',[]));
        end;
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testreflectionsunit_test_silent():Boolean;
begin
    Result := TestReflections(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testreflectionsunit_test():Boolean;
begin
    Result := TestReflections(False);
end;


end.