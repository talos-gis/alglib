
program _demo;
Array[0]
var
    N : AlglibInteger;
    I : AlglibInteger;
    J : AlglibInteger;
    C1 : Double;
    X : Double;
    A : TReal2DArray;
begin
    Write(Format('                 CONDITION NUMBERS'#13#10'',[]));
    Write(Format('OF VANDERMONDE AND CHEBYSHEV INTERPOLATION MATRICES'#13#10''#13#10'',[]));
    Write(Format('    VANDERMONDE   CHEBYSHEV'#13#10'',[]));
    Write(Format('  N      1-norm      1-norm'#13#10'',[]));
    N:=2;
    while N<=14 do
    begin
        SetLength(A, N, N);
        Write(Format('%3d',[
            N]));
        
        //
        // Vandermone matrix
        //
        I:=0;
        while I<=N-1 do
        begin
            X := AP_Double(2*I)/(N-1)-1;
            A[I,0] := 1;
            J:=1;
            while J<=N-1 do
            begin
                A[I,J] := A[I,J-1]*X;
                Inc(J);
            end;
            Inc(I);
        end;
        C1 := 1/RMatrixRCond1(A, N);
        Write(Format(' %11.1f',[
            C1]));
        
        //
        // Chebyshev interpolation matrix
        //
        I:=0;
        while I<=N-1 do
        begin
            X := AP_Double(2*I)/(N-1)-1;
            A[I,0] := 1;
            if N>=2 then
            begin
                A[I,1] := X;
            end;
            J:=2;
            while J<=N-1 do
            begin
                A[I,J] := 2*X*A[I,J-1]-A[I,J-2];
                Inc(J);
            end;
            Inc(I);
        end;
        C1 := 1/RMatrixRCond1(A, N);
        Write(Format(' %11.1f'#13#10'',[
            C1]));
        Inc(N);
    end;
end.