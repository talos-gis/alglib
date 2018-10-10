unit testblasunit;
interface
uses Math, Sysutils, Ap, blas;

function TestBLAS(Silent : Boolean):Boolean;
function testblasunit_test_silent():Boolean;
function testblasunit_test():Boolean;

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


function TestBLAS(Silent : Boolean):Boolean;
var
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    N : AlglibInteger;
    I : AlglibInteger;
    I1 : AlglibInteger;
    I2 : AlglibInteger;
    J : AlglibInteger;
    J1 : AlglibInteger;
    J2 : AlglibInteger;
    L : AlglibInteger;
    K : AlglibInteger;
    R : AlglibInteger;
    I3 : AlglibInteger;
    J3 : AlglibInteger;
    Col1 : AlglibInteger;
    Col2 : AlglibInteger;
    Row1 : AlglibInteger;
    Row2 : AlglibInteger;
    X1 : TReal1DArray;
    X2 : TReal1DArray;
    A : TReal2DArray;
    B : TReal2DArray;
    C1 : TReal2DArray;
    C2 : TReal2DArray;
    Err : Double;
    E1 : Double;
    E2 : Double;
    E3 : Double;
    V : Double;
    Scl1 : Double;
    Scl2 : Double;
    Scl3 : Double;
    Was1 : Boolean;
    Was2 : Boolean;
    Trans1 : Boolean;
    Trans2 : Boolean;
    Threshold : Double;
    N2Errors : Boolean;
    HSNErrors : Boolean;
    AMaxErrors : Boolean;
    MVErrors : Boolean;
    ITErrors : Boolean;
    CTErrors : Boolean;
    MMErrors : Boolean;
    WasErrors : Boolean;
    i_ : AlglibInteger;
begin
    N2Errors := False;
    AMaxErrors := False;
    HSNErrors := False;
    MVErrors := False;
    ITErrors := False;
    CTErrors := False;
    MMErrors := False;
    WasErrors := False;
    Threshold := 10000*MachineEpsilon;
    
    //
    // Test Norm2
    //
    PassCount := 1000;
    E1 := 0;
    E2 := 0;
    E3 := 0;
    Scl2 := 0.5*MaxRealNumber;
    Scl3 := 2*MinRealNumber;
    Pass:=1;
    while Pass<=PassCount do
    begin
        N := 1+RandomInteger(1000);
        I1 := RandomInteger(10);
        I2 := N+I1-1;
        SetLength(X1, I2+1);
        SetLength(X2, I2+1);
        I:=I1;
        while I<=I2 do
        begin
            X1[I] := 2*RandomReal-1;
            Inc(I);
        end;
        V := 0;
        I:=I1;
        while I<=I2 do
        begin
            V := V+AP_Sqr(X1[I]);
            Inc(I);
        end;
        V := Sqrt(V);
        E1 := Max(E1, AbsReal(V-VectorNorm2(X1, I1, I2)));
        I:=I1;
        while I<=I2 do
        begin
            X2[I] := Scl2*X1[I];
            Inc(I);
        end;
        E2 := Max(E2, AbsReal(V*Scl2-VectorNorm2(X2, I1, I2)));
        I:=I1;
        while I<=I2 do
        begin
            X2[I] := Scl3*X1[I];
            Inc(I);
        end;
        E3 := Max(E3, AbsReal(V*Scl3-VectorNorm2(X2, I1, I2)));
        Inc(Pass);
    end;
    E2 := E2/Scl2;
    E3 := E3/Scl3;
    N2Errors := AP_FP_Greater_Eq(E1,Threshold) or AP_FP_Greater_Eq(E2,Threshold) or AP_FP_Greater_Eq(E3,Threshold);
    
    //
    // Testing VectorAbsMax, Column/Row AbsMax
    //
    SetLength(X1, 5+1);
    X1[1] := 2.0;
    X1[2] := 0.2;
    X1[3] := -1.3;
    X1[4] := 0.7;
    X1[5] := -3.0;
    AMaxErrors := (VectorIdxAbsMax(X1, 1, 5)<>5) or (VectorIdxAbsMax(X1, 1, 4)<>1) or (VectorIdxAbsMax(X1, 2, 4)<>3);
    N := 30;
    SetLength(X1, N+1);
    SetLength(A, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        J:=1;
        while J<=N do
        begin
            A[I,J] := 2*RandomReal-1;
            Inc(J);
        end;
        Inc(I);
    end;
    Was1 := False;
    Was2 := False;
    Pass:=1;
    while Pass<=1000 do
    begin
        J := 1+RandomInteger(N);
        I1 := 1+RandomInteger(N);
        I2 := I1+RandomInteger(N+1-I1);
        for i_ := I1 to I2 do
        begin
            X1[i_] := A[i_,J];
        end;
        if VectorIdxAbsMax(X1, I1, I2)<>ColumnIdxAbsMax(A, I1, I2, J) then
        begin
            Was1 := True;
        end;
        I := 1+RandomInteger(N);
        J1 := 1+RandomInteger(N);
        J2 := J1+RandomInteger(N+1-J1);
        APVMove(@X1[0], J1, J2, @A[I][0], J1, J2);
        if VectorIdxAbsMax(X1, J1, J2)<>RowIdxAbsMax(A, J1, J2, I) then
        begin
            Was2 := True;
        end;
        Inc(Pass);
    end;
    AMaxErrors := AMaxErrors or Was1 or Was2;
    
    //
    // Testing upper Hessenberg 1-norm
    //
    SetLength(A, 3+1, 3+1);
    SetLength(X1, 3+1);
    A[1,1] := 2;
    A[1,2] := 3;
    A[1,3] := 1;
    A[2,1] := 4;
    A[2,2] := -5;
    A[2,3] := 8;
    A[3,1] := 99;
    A[3,2] := 3;
    A[3,3] := 1;
    HSNErrors := AP_FP_Greater(AbsReal(UpperHessenberg1Norm(A, 1, 3, 1, 3, X1)-11),Threshold);
    
    //
    // Testing MatrixVectorMultiply
    //
    SetLength(A, 3+1, 5+1);
    SetLength(X1, 3+1);
    SetLength(X2, 2+1);
    A[2,3] := 2;
    A[2,4] := -1;
    A[2,5] := -1;
    A[3,3] := 1;
    A[3,4] := -2;
    A[3,5] := 2;
    X1[1] := 1;
    X1[2] := 2;
    X1[3] := 1;
    X2[1] := -1;
    X2[2] := -1;
    MatrixVectorMultiply(A, 2, 3, 3, 5, False, X1, 1, 3, 1.0, X2, 1, 2, 1.0);
    MatrixVectorMultiply(A, 2, 3, 3, 5, True, X2, 1, 2, 1.0, X1, 1, 3, 1.0);
    E1 := AbsReal(X1[1]+5)+AbsReal(X1[2]-8)+AbsReal(X1[3]+1)+AbsReal(X2[1]+2)+AbsReal(X2[2]+2);
    X1[1] := 1;
    X1[2] := 2;
    X1[3] := 1;
    X2[1] := -1;
    X2[2] := -1;
    MatrixVectorMultiply(A, 2, 3, 3, 5, False, X1, 1, 3, 1.0, X2, 1, 2, 0.0);
    MatrixVectorMultiply(A, 2, 3, 3, 5, True, X2, 1, 2, 1.0, X1, 1, 3, 0.0);
    E2 := AbsReal(X1[1]+3)+AbsReal(X1[2]-3)+AbsReal(X1[3]+1)+AbsReal(X2[1]+1)+AbsReal(X2[2]+1);
    MVErrors := AP_FP_Greater_Eq(E1+E2,Threshold);
    
    //
    // testing inplace transpose
    //
    N := 10;
    SetLength(A, N+1, N+1);
    SetLength(B, N+1, N+1);
    SetLength(X1, N-1+1);
    I:=1;
    while I<=N do
    begin
        J:=1;
        while J<=N do
        begin
            A[I,J] := RandomReal;
            Inc(J);
        end;
        Inc(I);
    end;
    PassCount := 10000;
    Was1 := False;
    Pass:=1;
    while Pass<=PassCount do
    begin
        I1 := 1+RandomInteger(N);
        I2 := I1+RandomInteger(N-I1+1);
        J1 := 1+RandomInteger(N-(I2-I1));
        J2 := J1+(I2-I1);
        CopyMatrix(A, I1, I2, J1, J2, B, I1, I2, J1, J2);
        InplaceTranspose(B, I1, I2, J1, J2, X1);
        I:=I1;
        while I<=I2 do
        begin
            J:=J1;
            while J<=J2 do
            begin
                if AP_FP_Neq(A[I,J],B[I1+(J-J1),J1+(I-I1)]) then
                begin
                    Was1 := True;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        Inc(Pass);
    end;
    ITErrors := Was1;
    
    //
    // testing copy and transpose
    //
    N := 10;
    SetLength(A, N+1, N+1);
    SetLength(B, N+1, N+1);
    I:=1;
    while I<=N do
    begin
        J:=1;
        while J<=N do
        begin
            A[I,J] := RandomReal;
            Inc(J);
        end;
        Inc(I);
    end;
    PassCount := 10000;
    Was1 := False;
    Pass:=1;
    while Pass<=PassCount do
    begin
        I1 := 1+RandomInteger(N);
        I2 := I1+RandomInteger(N-I1+1);
        J1 := 1+RandomInteger(N);
        J2 := J1+RandomInteger(N-J1+1);
        CopyAndTranspose(A, I1, I2, J1, J2, B, J1, J2, I1, I2);
        I:=I1;
        while I<=I2 do
        begin
            J:=J1;
            while J<=J2 do
            begin
                if AP_FP_Neq(A[I,J],B[J,I]) then
                begin
                    Was1 := True;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        Inc(Pass);
    end;
    CTErrors := Was1;
    
    //
    // Testing MatrixMatrixMultiply
    //
    N := 10;
    SetLength(A, 2*N+1, 2*N+1);
    SetLength(B, 2*N+1, 2*N+1);
    SetLength(C1, 2*N+1, 2*N+1);
    SetLength(C2, 2*N+1, 2*N+1);
    SetLength(X1, N+1);
    SetLength(X2, N+1);
    I:=1;
    while I<=2*N do
    begin
        J:=1;
        while J<=2*N do
        begin
            A[I,J] := RandomReal;
            B[I,J] := RandomReal;
            Inc(J);
        end;
        Inc(I);
    end;
    PassCount := 1000;
    Was1 := False;
    Pass:=1;
    while Pass<=PassCount do
    begin
        I:=1;
        while I<=2*N do
        begin
            J:=1;
            while J<=2*N do
            begin
                C1[I,J] := 2.1*I+3.1*J;
                C2[I,J] := C1[I,J];
                Inc(J);
            end;
            Inc(I);
        end;
        L := 1+RandomInteger(N);
        K := 1+RandomInteger(N);
        R := 1+RandomInteger(N);
        I1 := 1+RandomInteger(N);
        J1 := 1+RandomInteger(N);
        I2 := 1+RandomInteger(N);
        J2 := 1+RandomInteger(N);
        I3 := 1+RandomInteger(N);
        J3 := 1+RandomInteger(N);
        Trans1 := AP_FP_Greater(RandomReal,0.5);
        Trans2 := AP_FP_Greater(RandomReal,0.5);
        if Trans1 then
        begin
            Col1 := L;
            Row1 := K;
        end
        else
        begin
            Col1 := K;
            Row1 := L;
        end;
        if Trans2 then
        begin
            Col2 := K;
            Row2 := R;
        end
        else
        begin
            Col2 := R;
            Row2 := K;
        end;
        Scl1 := RandomReal;
        Scl2 := RandomReal;
        MatrixMatrixMultiply(A, I1, I1+Row1-1, J1, J1+Col1-1, Trans1, B, I2, I2+Row2-1, J2, J2+Col2-1, Trans2, Scl1, C1, I3, I3+L-1, J3, J3+R-1, Scl2, X1);
        NaiveMatrixMatrixMultiply(A, I1, I1+Row1-1, J1, J1+Col1-1, Trans1, B, I2, I2+Row2-1, J2, J2+Col2-1, Trans2, Scl1, C2, I3, I3+L-1, J3, J3+R-1, Scl2);
        Err := 0;
        I:=1;
        while I<=L do
        begin
            J:=1;
            while J<=R do
            begin
                Err := Max(Err, AbsReal(C1[I3+I-1,J3+J-1]-C2[I3+I-1,J3+J-1]));
                Inc(J);
            end;
            Inc(I);
        end;
        if AP_FP_Greater(Err,Threshold) then
        begin
            Was1 := True;
            Break;
        end;
        Inc(Pass);
    end;
    MMErrors := Was1;
    
    //
    // report
    //
    WasErrors := N2Errors or AMaxErrors or HSNErrors or MVErrors or ITErrors or CTErrors or MMErrors;
    if  not Silent then
    begin
        Write(Format('TESTING BLAS'#13#10'',[]));
        Write(Format('VectorNorm2:                             ',[]));
        if N2Errors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('AbsMax (vector/row/column):              ',[]));
        if AMaxErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('UpperHessenberg1Norm:                    ',[]));
        if HSNErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('MatrixVectorMultiply:                    ',[]));
        if MVErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('InplaceTranspose:                        ',[]));
        if ITErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('CopyAndTranspose:                        ',[]));
        if CTErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('MatrixMatrixMultiply:                    ',[]));
        if MMErrors then
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
Silent unit test
*************************************************************************)
function testblasunit_test_silent():Boolean;
begin
    Result := TestBLAS(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testblasunit_test():Boolean;
begin
    Result := TestBLAS(False);
end;


end.