unit testldltunit;
interface
uses Math, Sysutils, Ap, ldlt;

function TestLDLT(Silent : Boolean):Boolean;
function testldltunit_test_silent():Boolean;
function testldltunit_test():Boolean;

implementation

procedure GenerateMatrix(var A : TReal2DArray;
     N : AlglibInteger;
     Task : AlglibInteger);forward;


function TestLDLT(Silent : Boolean):Boolean;
var
    A : TReal2DArray;
    A2 : TReal2DArray;
    L : TReal2DArray;
    D : TReal2DArray;
    U : TReal2DArray;
    T : TReal2DArray;
    T2 : TReal2DArray;
    P : TInteger1DArray;
    N : AlglibInteger;
    Pass : AlglibInteger;
    MTask : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    MinIJ : AlglibInteger;
    UpperIn : Boolean;
    CR : Boolean;
    V : Double;
    Err : Double;
    WasErrors : Boolean;
    PassCount : AlglibInteger;
    MaxN : AlglibInteger;
    HTask : AlglibInteger;
    Threshold : Double;
    i_ : AlglibInteger;
begin
    Err := 0;
    PassCount := 10;
    MaxN := 20;
    Threshold := 1000*MachineEpsilon;
    WasErrors := False;
    
    //
    // Test
    //
    N:=1;
    while N<=MaxN do
    begin
        SetLength(A, N-1+1, N-1+1);
        SetLength(A2, N-1+1, N-1+1);
        SetLength(L, N-1+1, N-1+1);
        SetLength(U, N-1+1, N-1+1);
        SetLength(D, N-1+1, N-1+1);
        SetLength(T, N-1+1, N-1+1);
        SetLength(T2, N-1+1, N-1+1);
        MTask:=0;
        while MTask<=2 do
        begin
            HTask:=0;
            while HTask<=1 do
            begin
                Pass:=1;
                while Pass<=PassCount do
                begin
                    UpperIn := HTask=0;
                    
                    //
                    // Prepare task:
                    // * A contains symmetric matrix
                    // * A2 contains its upper (or lower) half
                    //
                    GenerateMatrix(A, N, MTask);
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            A2[I,J] := A[I,J];
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
                            if UpperIn then
                            begin
                                if J<I then
                                begin
                                    A2[I,J] := 0;
                                end;
                            end
                            else
                            begin
                                if I<J then
                                begin
                                    A2[I,J] := 0;
                                end;
                            end;
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    
                    //
                    // LDLt
                    //
                    SMatrixLDLT(A2, N, UpperIn, P);
                    
                    //
                    // Test (upper or lower)
                    //
                    if UpperIn then
                    begin
                        
                        //
                        // Unpack D
                        //
                        I:=0;
                        while I<=N-1 do
                        begin
                            J:=0;
                            while J<=N-1 do
                            begin
                                D[I,J] := 0;
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        K := 0;
                        while K<N do
                        begin
                            if P[K]>=0 then
                            begin
                                D[K,K] := A2[K,K];
                                K := K+1;
                            end
                            else
                            begin
                                D[K,K] := A2[K,K];
                                D[K,K+1] := A2[K,K+1];
                                D[K+1,K] := A2[K,K+1];
                                D[K+1,K+1] := A2[K+1,K+1];
                                K := K+2;
                            end;
                        end;
                        
                        //
                        // Unpack U
                        //
                        I:=0;
                        while I<=N-1 do
                        begin
                            J:=0;
                            while J<=N-1 do
                            begin
                                U[I,J] := 0;
                                Inc(J);
                            end;
                            U[I,I] := 1;
                            Inc(I);
                        end;
                        K := 0;
                        while K<N do
                        begin
                            
                            //
                            // unpack Uk
                            //
                            I:=0;
                            while I<=N-1 do
                            begin
                                J:=0;
                                while J<=N-1 do
                                begin
                                    T[I,J] := 0;
                                    Inc(J);
                                end;
                                T[I,I] := 1;
                                Inc(I);
                            end;
                            if P[K]>=0 then
                            begin
                                I:=0;
                                while I<=K-1 do
                                begin
                                    T[I,K] := A2[I,K];
                                    Inc(I);
                                end;
                            end
                            else
                            begin
                                I:=0;
                                while I<=K-1 do
                                begin
                                    T[I,K] := A2[I,K];
                                    T[I,K+1] := A2[I,K+1];
                                    Inc(I);
                                end;
                            end;
                            
                            //
                            // multiply U
                            //
                            I:=0;
                            while I<=N-1 do
                            begin
                                J:=0;
                                while J<=N-1 do
                                begin
                                    V := 0.0;
                                    for i_ := 0 to N-1 do
                                    begin
                                        V := V + T[I,i_]*U[i_,J];
                                    end;
                                    T2[I,J] := V;
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
                                    U[I,J] := T2[I,J];
                                    Inc(J);
                                end;
                                Inc(I);
                            end;
                            
                            //
                            // permutations
                            //
                            if P[K]>=0 then
                            begin
                                J:=0;
                                while J<=N-1 do
                                begin
                                    V := U[K,J];
                                    U[K,J] := U[P[K],J];
                                    U[P[K],J] := V;
                                    Inc(J);
                                end;
                                K := K+1;
                            end
                            else
                            begin
                                J:=0;
                                while J<=N-1 do
                                begin
                                    V := U[K,J];
                                    U[K,J] := U[N+P[K],J];
                                    U[N+P[K],J] := V;
                                    Inc(J);
                                end;
                                K := K+2;
                            end;
                        end;
                        
                        //
                        // Calculate U*D*U', store result in T2
                        //
                        I:=0;
                        while I<=N-1 do
                        begin
                            J:=0;
                            while J<=N-1 do
                            begin
                                V := 0.0;
                                for i_ := 0 to N-1 do
                                begin
                                    V := V + U[I,i_]*D[i_,J];
                                end;
                                T[I,J] := V;
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
                                V := APVDotProduct(@T[I][0], 0, N-1, @U[J][0], 0, N-1);
                                T2[I,J] := V;
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                    end
                    else
                    begin
                        
                        //
                        // Unpack D
                        //
                        I:=0;
                        while I<=N-1 do
                        begin
                            J:=0;
                            while J<=N-1 do
                            begin
                                D[I,J] := 0;
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                        K := 0;
                        while K<N do
                        begin
                            if P[K]>=0 then
                            begin
                                D[K,K] := A2[K,K];
                                K := K+1;
                            end
                            else
                            begin
                                D[K,K] := A2[K,K];
                                D[K,K+1] := A2[K+1,K];
                                D[K+1,K] := A2[K+1,K];
                                D[K+1,K+1] := A2[K+1,K+1];
                                K := K+2;
                            end;
                        end;
                        
                        //
                        // Unpack L
                        //
                        I:=0;
                        while I<=N-1 do
                        begin
                            J:=0;
                            while J<=N-1 do
                            begin
                                L[I,J] := 0;
                                Inc(J);
                            end;
                            L[I,I] := 1;
                            Inc(I);
                        end;
                        K := 0;
                        while K<N do
                        begin
                            
                            //
                            // permutations
                            //
                            if P[K]>=0 then
                            begin
                                I:=0;
                                while I<=N-1 do
                                begin
                                    V := L[I,K];
                                    L[I,K] := L[I,P[K]];
                                    L[I,P[K]] := V;
                                    Inc(I);
                                end;
                            end
                            else
                            begin
                                I:=0;
                                while I<=N-1 do
                                begin
                                    V := L[I,K+1];
                                    L[I,K+1] := L[I,N+P[K+1]];
                                    L[I,N+P[K+1]] := V;
                                    Inc(I);
                                end;
                            end;
                            
                            //
                            // unpack Lk
                            //
                            I:=0;
                            while I<=N-1 do
                            begin
                                J:=0;
                                while J<=N-1 do
                                begin
                                    T[I,J] := 0;
                                    Inc(J);
                                end;
                                T[I,I] := 1;
                                Inc(I);
                            end;
                            if P[K]>=0 then
                            begin
                                I:=K+1;
                                while I<=N-1 do
                                begin
                                    T[I,K] := A2[I,K];
                                    Inc(I);
                                end;
                            end
                            else
                            begin
                                I:=K+2;
                                while I<=N-1 do
                                begin
                                    T[I,K] := A2[I,K];
                                    T[I,K+1] := A2[I,K+1];
                                    Inc(I);
                                end;
                            end;
                            
                            //
                            // multiply L
                            //
                            I:=0;
                            while I<=N-1 do
                            begin
                                J:=0;
                                while J<=N-1 do
                                begin
                                    V := 0.0;
                                    for i_ := 0 to N-1 do
                                    begin
                                        V := V + L[I,i_]*T[i_,J];
                                    end;
                                    T2[I,J] := V;
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
                                    L[I,J] := T2[I,J];
                                    Inc(J);
                                end;
                                Inc(I);
                            end;
                            
                            //
                            // Next K
                            //
                            if P[K]>=0 then
                            begin
                                K := K+1;
                            end
                            else
                            begin
                                K := K+2;
                            end;
                        end;
                        
                        //
                        // Calculate L*D*L', store result in T2
                        //
                        I:=0;
                        while I<=N-1 do
                        begin
                            J:=0;
                            while J<=N-1 do
                            begin
                                V := 0.0;
                                for i_ := 0 to N-1 do
                                begin
                                    V := V + L[I,i_]*D[i_,J];
                                end;
                                T[I,J] := V;
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
                                V := APVDotProduct(@T[I][0], 0, N-1, @L[J][0], 0, N-1);
                                T2[I,J] := V;
                                Inc(J);
                            end;
                            Inc(I);
                        end;
                    end;
                    
                    //
                    // Test
                    //
                    I:=0;
                    while I<=N-1 do
                    begin
                        J:=0;
                        while J<=N-1 do
                        begin
                            Err := Max(Err, AbsReal(A[I,J]-T2[I,J]));
                            Inc(J);
                        end;
                        Inc(I);
                    end;
                    Inc(Pass);
                end;
                Inc(HTask);
            end;
            Inc(MTask);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    WasErrors := AP_FP_Greater(Err,Threshold);
    if  not Silent then
    begin
        Write(Format('TESTING LDLT DECOMPOSITION'#13#10'',[]));
        Write(Format('ERROR:                                   %5.4e'#13#10'',[
            Err]));
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


procedure GenerateMatrix(var A : TReal2DArray;
     N : AlglibInteger;
     Task : AlglibInteger);
var
    I : AlglibInteger;
    J : AlglibInteger;
begin
    if Task=0 then
    begin
        
        //
        // Zero matrix
        //
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=N-1 do
            begin
                A[I,J] := 0;
                Inc(J);
            end;
            Inc(I);
        end;
    end;
    if Task=1 then
    begin
        
        //
        // Sparse matrix
        //
        I:=0;
        while I<=N-1 do
        begin
            J:=I+1;
            while J<=N-1 do
            begin
                if AP_FP_Greater(RandomReal,0.95) then
                begin
                    A[I,J] := 2*RandomReal-1;
                end
                else
                begin
                    A[I,J] := 0;
                end;
                A[J,I] := A[I,J];
                Inc(J);
            end;
            if AP_FP_Greater(RandomReal,0.95) then
            begin
                A[I,I] := (2*RandomInteger(2)-1)*(0.8+RandomReal);
            end
            else
            begin
                A[I,I] := 0;
            end;
            Inc(I);
        end;
    end;
    if Task=2 then
    begin
        
        //
        // Dense matrix
        //
        I:=0;
        while I<=N-1 do
        begin
            J:=I+1;
            while J<=N-1 do
            begin
                A[I,J] := 2*RandomReal-1;
                A[J,I] := A[I,J];
                Inc(J);
            end;
            A[I,I] := (2*RandomInteger(2)-1)*(0.8+RandomReal);
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testldltunit_test_silent():Boolean;
begin
    Result := TestLDLT(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testldltunit_test():Boolean;
begin
    Result := TestLDLT(False);
end;


end.