unit testinverseupdateunit;
interface
uses Math, Sysutils, Ap, inverseupdate;

function TestInverseUpdate(Silent : Boolean):Boolean;
function testinverseupdateunit_test_silent():Boolean;
function testinverseupdateunit_test():Boolean;

implementation

procedure MakeACopy(const A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var B : TReal2DArray);forward;
procedure MatLU(var A : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger;
     var Pivots : TInteger1DArray);forward;
procedure GenerateRandomOrthogonalMatrix(var A0 : TReal2DArray;
     N : AlglibInteger);forward;
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
function MatrixDiff(const A : TReal2DArray;
     const B : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger):Double;forward;
function UpdAndInv(var A : TReal2DArray;
     const U : TReal1DArray;
     const V : TReal1DArray;
     N : AlglibInteger):Boolean;forward;


function TestInverseUpdate(Silent : Boolean):Boolean;
var
    A : TReal2DArray;
    InvA : TReal2DArray;
    B1 : TReal2DArray;
    B2 : TReal2DArray;
    U : TReal1DArray;
    V : TReal1DArray;
    N : AlglibInteger;
    MaxN : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    UpdRow : AlglibInteger;
    UpdCol : AlglibInteger;
    Val : Double;
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    WasErrors : Boolean;
    Threshold : Double;
    C : Double;
begin
    WasErrors := False;
    MaxN := 10;
    PassCount := 100;
    Threshold := 1.0E-6;
    
    //
    // process
    //
    N:=1;
    while N<=MaxN do
    begin
        SetLength(A, N-1+1, N-1+1);
        SetLength(B1, N-1+1, N-1+1);
        SetLength(B2, N-1+1, N-1+1);
        SetLength(U, N-1+1);
        SetLength(V, N-1+1);
        Pass:=1;
        while Pass<=PassCount do
        begin
            C := Exp(RandomReal*Ln(10));
            GenerateRandomMatrixCond(A, N, C);
            MakeACopy(A, N, N, InvA);
            if  not InvMat(InvA, N) then
            begin
                WasErrors := True;
                Break;
            end;
            
            //
            // Test simple update
            //
            UpdRow := RandomInteger(N);
            UpdCol := RandomInteger(N);
            Val := 0.1*(2*RandomReal-1);
            I:=0;
            while I<=N-1 do
            begin
                if I=UpdRow then
                begin
                    U[I] := Val;
                end
                else
                begin
                    U[I] := 0;
                end;
                if I=UpdCol then
                begin
                    V[I] := 1;
                end
                else
                begin
                    V[I] := 0;
                end;
                Inc(I);
            end;
            MakeACopy(A, N, N, B1);
            if  not UpdAndInv(B1, U, V, N) then
            begin
                WasErrors := True;
                Break;
            end;
            MakeACopy(InvA, N, N, B2);
            RMatrixInvUpdateSimple(B2, N, UpdRow, UpdCol, Val);
            WasErrors := WasErrors or AP_FP_Greater(MatrixDiff(B1, B2, N, N),Threshold);
            
            //
            // Test row update
            //
            UpdRow := RandomInteger(N);
            I:=0;
            while I<=N-1 do
            begin
                if I=UpdRow then
                begin
                    U[I] := 1;
                end
                else
                begin
                    U[I] := 0;
                end;
                V[I] := 0.1*(2*RandomReal-1);
                Inc(I);
            end;
            MakeACopy(A, N, N, B1);
            if  not UpdAndInv(B1, U, V, N) then
            begin
                WasErrors := True;
                Break;
            end;
            MakeACopy(InvA, N, N, B2);
            RMatrixInvUpdateRow(B2, N, UpdRow, V);
            WasErrors := WasErrors or AP_FP_Greater(MatrixDiff(B1, B2, N, N),Threshold);
            
            //
            // Test column update
            //
            UpdCol := RandomInteger(N);
            I:=0;
            while I<=N-1 do
            begin
                if I=UpdCol then
                begin
                    V[I] := 1;
                end
                else
                begin
                    V[I] := 0;
                end;
                U[I] := 0.1*(2*RandomReal-1);
                Inc(I);
            end;
            MakeACopy(A, N, N, B1);
            if  not UpdAndInv(B1, U, V, N) then
            begin
                WasErrors := True;
                Break;
            end;
            MakeACopy(InvA, N, N, B2);
            RMatrixInvUpdateColumn(B2, N, UpdCol, U);
            WasErrors := WasErrors or AP_FP_Greater(MatrixDiff(B1, B2, N, N),Threshold);
            
            //
            // Test full update
            //
            I:=0;
            while I<=N-1 do
            begin
                V[I] := 0.1*(2*RandomReal-1);
                U[I] := 0.1*(2*RandomReal-1);
                Inc(I);
            end;
            MakeACopy(A, N, N, B1);
            if  not UpdAndInv(B1, U, V, N) then
            begin
                WasErrors := True;
                Break;
            end;
            MakeACopy(InvA, N, N, B2);
            RMatrixInvUpdateUV(B2, N, U, V);
            WasErrors := WasErrors or AP_FP_Greater(MatrixDiff(B1, B2, N, N),Threshold);
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    if  not Silent then
    begin
        Write(Format('TESTING INVERSE UPDATE (REAL)'#13#10'',[]));
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
Generate matrix with given condition number C (2-norm)
*************************************************************************)
procedure GenerateRandomOrthogonalMatrix(var A0 : TReal2DArray;
     N : AlglibInteger);
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
    A : TReal2DArray;
    SM : Double;
begin
    if N<=0 then
    begin
        Exit;
    end;
    SetLength(W, N+1);
    SetLength(V, N+1);
    SetLength(A, N+1, N+1);
    SetLength(A0, N-1+1, N-1+1);
    
    //
    // Prepare A
    //
    I:=1;
    while I<=N do
    begin
        J:=1;
        while J<=N do
        begin
            if I=J then
            begin
                A[I,J] := 1;
            end
            else
            begin
                A[I,J] := 0;
            end;
            Inc(J);
        end;
        Inc(I);
    end;
    
    //
    // Calculate A using Stewart algorithm
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


procedure GenerateRandomMatrixCond(var A0 : TReal2DArray;
     N : AlglibInteger;
     C : Double);
var
    L1 : Double;
    L2 : Double;
    Q1 : TReal2DArray;
    Q2 : TReal2DArray;
    CC : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
begin
    GenerateRandomOrthogonalMatrix(Q1, N);
    GenerateRandomOrthogonalMatrix(Q2, N);
    SetLength(CC, N-1+1);
    L1 := 0;
    L2 := Ln(1/C);
    CC[0] := Exp(L1);
    I:=1;
    while I<=N-2 do
    begin
        CC[I] := Exp(RandomReal*(L2-L1)+L1);
        Inc(I);
    end;
    CC[N-1] := Exp(L2);
    SetLength(A0, N-1+1, N-1+1);
    I:=0;
    while I<=N-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            A0[I,J] := 0;
            K:=0;
            while K<=N-1 do
            begin
                A0[I,J] := A0[I,J]+Q1[I,K]*CC[K]*Q2[J,K];
                Inc(K);
            end;
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
    MatLU(A, N, N, Pivots);
    Result := InvMatLU(A, Pivots, N);
end;


(*************************************************************************
Diff
*************************************************************************)
function MatrixDiff(const A : TReal2DArray;
     const B : TReal2DArray;
     M : AlglibInteger;
     N : AlglibInteger):Double;
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    Result := 0;
    I:=0;
    while I<=M-1 do
    begin
        J:=0;
        while J<=N-1 do
        begin
            Result := Max(Result, AbsReal(B[I,J]-A[I,J]));
            Inc(J);
        end;
        Inc(I);
    end;
end;


(*************************************************************************
Update and inverse
*************************************************************************)
function UpdAndInv(var A : TReal2DArray;
     const U : TReal1DArray;
     const V : TReal1DArray;
     N : AlglibInteger):Boolean;
var
    Pivots : TInteger1DArray;
    I : AlglibInteger;
    R : Double;
begin
    I:=0;
    while I<=N-1 do
    begin
        R := U[I];
        APVAdd(@A[I][0], 0, N-1, @V[0], 0, N-1, R);
        Inc(I);
    end;
    MatLU(A, N, N, Pivots);
    Result := InvMatLU(A, Pivots, N);
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testinverseupdateunit_test_silent():Boolean;
begin
    Result := TestInverseUpdate(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testinverseupdateunit_test():Boolean;
begin
    Result := TestInverseUpdate(False);
end;


end.