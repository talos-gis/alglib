unit testgq;
interface
uses Math, Sysutils, Ap, blas, rotations, tdevd, gammafunc, gq;

function TestGQunit(Silent : Boolean):Boolean;
procedure BuildGaussHermiteQuadrature(n : AlglibInteger;
     var x : TReal1DArray;
     var w : TReal1DArray);
function testgq_test_silent():Boolean;
function testgq_test():Boolean;

implementation

function MapKind(K : AlglibInteger):Double;forward;
procedure BuildGaussLegendreQuadrature(n : AlglibInteger;
     var x : TReal1DArray;
     var w : TReal1DArray);forward;
procedure BuildGaussJacobiQuadrature(n : AlglibInteger;
     Alpha : Double;
     Beta : Double;
     var x : TReal1DArray;
     var w : TReal1DArray);forward;
procedure BuildGaussLaguerreQuadrature(n : AlglibInteger;
     Alpha : Double;
     var x : TReal1DArray;
     var w : TReal1DArray);forward;


(*************************************************************************
Test
*************************************************************************)
function TestGQunit(Silent : Boolean):Boolean;
var
    Alpha : TReal1DArray;
    Beta : TReal1DArray;
    X : TReal1DArray;
    W : TReal1DArray;
    X2 : TReal1DArray;
    W2 : TReal1DArray;
    Err : Double;
    N : AlglibInteger;
    I : AlglibInteger;
    Info : AlglibInteger;
    AKind : AlglibInteger;
    BKind : AlglibInteger;
    AlphaC : Double;
    BetaC : Double;
    ErrTol : Double;
    NonStrictErrTol : Double;
    StrictErrTol : Double;
    RecErrors : Boolean;
    SpecErrors : Boolean;
    WasErrors : Boolean;
begin
    RecErrors := False;
    SpecErrors := False;
    WasErrors := False;
    ErrTol := 1.0E-12;
    NonStrictErrTol := 1.0E-6;
    StrictErrTol := 1000*MachineEpsilon;
    
    //
    // Three tests for rec-based Gauss quadratures with known weights/nodes:
    // 1. Gauss-Legendre with N=2
    // 2. Gauss-Legendre with N=5
    // 3. Gauss-Chebyshev with N=1, 2, 4, 8, ..., 512
    //
    Err := 0;
    SetLength(Alpha, 2);
    SetLength(Beta, 2);
    Alpha[0] := 0;
    Alpha[1] := 0;
    Beta[1] := AP_Double(1)/(4*1*1-1);
    GQGenerateRec(Alpha, Beta, 2.0, 2, Info, X, W);
    if Info>0 then
    begin
        Err := Max(Err, AbsReal(X[0]+Sqrt(3)/3));
        Err := Max(Err, AbsReal(X[1]-Sqrt(3)/3));
        Err := Max(Err, AbsReal(W[0]-1));
        Err := Max(Err, AbsReal(W[1]-1));
        I:=0;
        while I<=0 do
        begin
            RecErrors := RecErrors or AP_FP_Greater_Eq(X[I],X[I+1]);
            Inc(I);
        end;
    end
    else
    begin
        RecErrors := True;
    end;
    SetLength(Alpha, 5);
    SetLength(Beta, 5);
    Alpha[0] := 0;
    I:=1;
    while I<=4 do
    begin
        Alpha[I] := 0;
        Beta[I] := AP_Sqr(I)/(4*AP_Sqr(I)-1);
        Inc(I);
    end;
    GQGenerateRec(Alpha, Beta, 2.0, 5, Info, X, W);
    if Info>0 then
    begin
        Err := Max(Err, AbsReal(X[0]+Sqrt(245+14*Sqrt(70))/21));
        Err := Max(Err, AbsReal(X[0]+X[4]));
        Err := Max(Err, AbsReal(X[1]+Sqrt(245-14*Sqrt(70))/21));
        Err := Max(Err, AbsReal(X[1]+X[3]));
        Err := Max(Err, AbsReal(X[2]));
        Err := Max(Err, AbsReal(W[0]-(322-13*Sqrt(70))/900));
        Err := Max(Err, AbsReal(W[0]-W[4]));
        Err := Max(Err, AbsReal(W[1]-(322+13*Sqrt(70))/900));
        Err := Max(Err, AbsReal(W[1]-W[3]));
        Err := Max(Err, AbsReal(W[2]-AP_Double(128)/225));
        I:=0;
        while I<=3 do
        begin
            RecErrors := RecErrors or AP_FP_Greater_Eq(X[I],X[I+1]);
            Inc(I);
        end;
    end
    else
    begin
        RecErrors := True;
    end;
    N := 1;
    while N<=512 do
    begin
        SetLength(Alpha, N);
        SetLength(Beta, N);
        I:=0;
        while I<=N-1 do
        begin
            Alpha[I] := 0;
            if I=0 then
            begin
                Beta[I] := 0;
            end;
            if I=1 then
            begin
                Beta[I] := AP_Double(1)/2;
            end;
            if I>1 then
            begin
                Beta[I] := AP_Double(1)/4;
            end;
            Inc(I);
        end;
        GQGenerateRec(Alpha, Beta, Pi, N, Info, X, W);
        if Info>0 then
        begin
            I:=0;
            while I<=N-1 do
            begin
                Err := Max(Err, AbsReal(X[I]-Cos(Pi*(N-I-0.5)/N)));
                Err := Max(Err, AbsReal(W[I]-Pi/N));
                Inc(I);
            end;
            I:=0;
            while I<=N-2 do
            begin
                RecErrors := RecErrors or AP_FP_Greater_Eq(X[I],X[I+1]);
                Inc(I);
            end;
        end
        else
        begin
            RecErrors := True;
        end;
        N := N*2;
    end;
    RecErrors := RecErrors or AP_FP_Greater(Err,ErrTol);
    
    //
    // Three tests for rec-based Gauss-Lobatto quadratures with known weights/nodes:
    // 1. Gauss-Lobatto with N=3
    // 2. Gauss-Lobatto with N=4
    // 3. Gauss-Lobatto with N=6
    //
    Err := 0;
    SetLength(Alpha, 2);
    SetLength(Beta, 2);
    Alpha[0] := 0;
    Alpha[1] := 0;
    Beta[0] := 0;
    Beta[1] := AP_Double(1*1)/(4*1*1-1);
    GQGenerateGaussLobattoRec(Alpha, Beta, 2.0, -1, +1, 3, Info, X, W);
    if Info>0 then
    begin
        Err := Max(Err, AbsReal(X[0]+1));
        Err := Max(Err, AbsReal(X[1]));
        Err := Max(Err, AbsReal(X[2]-1));
        Err := Max(Err, AbsReal(W[0]-AP_Double(1)/3));
        Err := Max(Err, AbsReal(W[1]-AP_Double(4)/3));
        Err := Max(Err, AbsReal(W[2]-AP_Double(1)/3));
        I:=0;
        while I<=1 do
        begin
            RecErrors := RecErrors or AP_FP_Greater_Eq(X[I],X[I+1]);
            Inc(I);
        end;
    end
    else
    begin
        RecErrors := True;
    end;
    SetLength(Alpha, 3);
    SetLength(Beta, 3);
    Alpha[0] := 0;
    Alpha[1] := 0;
    Alpha[2] := 0;
    Beta[0] := 0;
    Beta[1] := AP_Double(1*1)/(4*1*1-1);
    Beta[2] := AP_Double(2*2)/(4*2*2-1);
    GQGenerateGaussLobattoRec(Alpha, Beta, 2.0, -1, +1, 4, Info, X, W);
    if Info>0 then
    begin
        Err := Max(Err, AbsReal(X[0]+1));
        Err := Max(Err, AbsReal(X[1]+Sqrt(5)/5));
        Err := Max(Err, AbsReal(X[2]-Sqrt(5)/5));
        Err := Max(Err, AbsReal(X[3]-1));
        Err := Max(Err, AbsReal(W[0]-AP_Double(1)/6));
        Err := Max(Err, AbsReal(W[1]-AP_Double(5)/6));
        Err := Max(Err, AbsReal(W[2]-AP_Double(5)/6));
        Err := Max(Err, AbsReal(W[3]-AP_Double(1)/6));
        I:=0;
        while I<=2 do
        begin
            RecErrors := RecErrors or AP_FP_Greater_Eq(X[I],X[I+1]);
            Inc(I);
        end;
    end
    else
    begin
        RecErrors := True;
    end;
    SetLength(Alpha, 5);
    SetLength(Beta, 5);
    Alpha[0] := 0;
    Alpha[1] := 0;
    Alpha[2] := 0;
    Alpha[3] := 0;
    Alpha[4] := 0;
    Beta[0] := 0;
    Beta[1] := AP_Double(1*1)/(4*1*1-1);
    Beta[2] := AP_Double(2*2)/(4*2*2-1);
    Beta[3] := AP_Double(3*3)/(4*3*3-1);
    Beta[4] := AP_Double(4*4)/(4*4*4-1);
    GQGenerateGaussLobattoRec(Alpha, Beta, 2.0, -1, +1, 6, Info, X, W);
    if Info>0 then
    begin
        Err := Max(Err, AbsReal(X[0]+1));
        Err := Max(Err, AbsReal(X[1]+Sqrt((7+2*Sqrt(7))/21)));
        Err := Max(Err, AbsReal(X[2]+Sqrt((7-2*Sqrt(7))/21)));
        Err := Max(Err, AbsReal(X[3]-Sqrt((7-2*Sqrt(7))/21)));
        Err := Max(Err, AbsReal(X[4]-Sqrt((7+2*Sqrt(7))/21)));
        Err := Max(Err, AbsReal(X[5]-1));
        Err := Max(Err, AbsReal(W[0]-AP_Double(1)/15));
        Err := Max(Err, AbsReal(W[1]-(14-Sqrt(7))/30));
        Err := Max(Err, AbsReal(W[2]-(14+Sqrt(7))/30));
        Err := Max(Err, AbsReal(W[3]-(14+Sqrt(7))/30));
        Err := Max(Err, AbsReal(W[4]-(14-Sqrt(7))/30));
        Err := Max(Err, AbsReal(W[5]-AP_Double(1)/15));
        I:=0;
        while I<=4 do
        begin
            RecErrors := RecErrors or AP_FP_Greater_Eq(X[I],X[I+1]);
            Inc(I);
        end;
    end
    else
    begin
        RecErrors := True;
    end;
    RecErrors := RecErrors or AP_FP_Greater(Err,ErrTol);
    
    //
    // Three tests for rec-based Gauss-Radau quadratures with known weights/nodes:
    // 1. Gauss-Radau with N=2
    // 2. Gauss-Radau with N=3
    // 3. Gauss-Radau with N=3 (another case)
    //
    Err := 0;
    SetLength(Alpha, 1);
    SetLength(Beta, 2);
    Alpha[0] := 0;
    Beta[0] := 0;
    Beta[1] := AP_Double(1*1)/(4*1*1-1);
    GQGenerateGaussRadauRec(Alpha, Beta, 2.0, -1, 2, Info, X, W);
    if Info>0 then
    begin
        Err := Max(Err, AbsReal(X[0]+1));
        Err := Max(Err, AbsReal(X[1]-AP_Double(1)/3));
        Err := Max(Err, AbsReal(W[0]-0.5));
        Err := Max(Err, AbsReal(W[1]-1.5));
        I:=0;
        while I<=0 do
        begin
            RecErrors := RecErrors or AP_FP_Greater_Eq(X[I],X[I+1]);
            Inc(I);
        end;
    end
    else
    begin
        RecErrors := True;
    end;
    SetLength(Alpha, 2);
    SetLength(Beta, 3);
    Alpha[0] := 0;
    Alpha[1] := 0;
    I:=0;
    while I<=2 do
    begin
        Beta[I] := AP_Sqr(I)/(4*AP_Sqr(I)-1);
        Inc(I);
    end;
    GQGenerateGaussRadauRec(Alpha, Beta, 2.0, -1, 3, Info, X, W);
    if Info>0 then
    begin
        Err := Max(Err, AbsReal(X[0]+1));
        Err := Max(Err, AbsReal(X[1]-(1-Sqrt(6))/5));
        Err := Max(Err, AbsReal(X[2]-(1+Sqrt(6))/5));
        Err := Max(Err, AbsReal(W[0]-AP_Double(2)/9));
        Err := Max(Err, AbsReal(W[1]-(16+Sqrt(6))/18));
        Err := Max(Err, AbsReal(W[2]-(16-Sqrt(6))/18));
        I:=0;
        while I<=1 do
        begin
            RecErrors := RecErrors or AP_FP_Greater_Eq(X[I],X[I+1]);
            Inc(I);
        end;
    end
    else
    begin
        RecErrors := True;
    end;
    SetLength(Alpha, 2);
    SetLength(Beta, 3);
    Alpha[0] := 0;
    Alpha[1] := 0;
    I:=0;
    while I<=2 do
    begin
        Beta[I] := AP_Sqr(I)/(4*AP_Sqr(I)-1);
        Inc(I);
    end;
    GQGenerateGaussRadauRec(Alpha, Beta, 2.0, +1, 3, Info, X, W);
    if Info>0 then
    begin
        Err := Max(Err, AbsReal(X[2]-1));
        Err := Max(Err, AbsReal(X[1]+(1-Sqrt(6))/5));
        Err := Max(Err, AbsReal(X[0]+(1+Sqrt(6))/5));
        Err := Max(Err, AbsReal(W[2]-AP_Double(2)/9));
        Err := Max(Err, AbsReal(W[1]-(16+Sqrt(6))/18));
        Err := Max(Err, AbsReal(W[0]-(16-Sqrt(6))/18));
        I:=0;
        while I<=1 do
        begin
            RecErrors := RecErrors or AP_FP_Greater_Eq(X[I],X[I+1]);
            Inc(I);
        end;
    end
    else
    begin
        RecErrors := True;
    end;
    RecErrors := RecErrors or AP_FP_Greater(Err,ErrTol);
    
    //
    // test recurrence-based special cases (Legendre, Jacobi, Hermite, ...)
    // against another implementation (polynomial root-finder)
    //
    N:=1;
    while N<=20 do
    begin
        
        //
        // test gauss-legendre
        //
        Err := 0;
        GQGenerateGaussLegendre(N, Info, X, W);
        if Info>0 then
        begin
            BuildGaussLegendreQuadrature(N, X2, W2);
            I:=0;
            while I<=N-1 do
            begin
                Err := Max(Err, AbsReal(X[I]-X2[I]));
                Err := Max(Err, AbsReal(W[I]-W2[I]));
                Inc(I);
            end;
        end
        else
        begin
            SpecErrors := True;
        end;
        SpecErrors := SpecErrors or AP_FP_Greater(Err,ErrTol);
        
        //
        // Test Gauss-Jacobi.
        // Since task is much more difficult we will use less strict
        // threshold.
        //
        Err := 0;
        AKind:=0;
        while AKind<=9 do
        begin
            BKind:=0;
            while BKind<=9 do
            begin
                AlphaC := MapKind(AKind);
                BetaC := MapKind(BKind);
                GQGenerateGaussJacobi(N, AlphaC, BetaC, Info, X, W);
                if Info>0 then
                begin
                    BuildGaussJacobiQuadrature(N, AlphaC, BetaC, X2, W2);
                    I:=0;
                    while I<=N-1 do
                    begin
                        Err := Max(Err, AbsReal(X[I]-X2[I]));
                        Err := Max(Err, AbsReal(W[I]-W2[I]));
                        Inc(I);
                    end;
                end
                else
                begin
                    SpecErrors := True;
                end;
                Inc(BKind);
            end;
            Inc(AKind);
        end;
        SpecErrors := SpecErrors or AP_FP_Greater(Err,NonStrictErrTol);
        
        //
        // special test for Gauss-Jacobi (Chebyshev weight
        // function with analytically known nodes/weights)
        //
        Err := 0;
        GQGenerateGaussJacobi(N, -0.5, -0.5, Info, X, W);
        if Info>0 then
        begin
            I:=0;
            while I<=N-1 do
            begin
                Err := Max(Err, AbsReal(X[I]+Cos(Pi*(I+0.5)/N)));
                Err := Max(Err, AbsReal(W[I]-Pi/N));
                Inc(I);
            end;
        end
        else
        begin
            SpecErrors := True;
        end;
        SpecErrors := SpecErrors or AP_FP_Greater(Err,StrictErrTol);
        
        //
        // Test Gauss-Laguerre
        //
        Err := 0;
        AKind:=0;
        while AKind<=9 do
        begin
            AlphaC := MapKind(AKind);
            GQGenerateGaussLaguerre(N, AlphaC, Info, X, W);
            if Info>0 then
            begin
                BuildGaussLaguerreQuadrature(N, AlphaC, X2, W2);
                I:=0;
                while I<=N-1 do
                begin
                    Err := Max(Err, AbsReal(X[I]-X2[I]));
                    Err := Max(Err, AbsReal(W[I]-W2[I]));
                    Inc(I);
                end;
            end
            else
            begin
                SpecErrors := True;
            end;
            Inc(AKind);
        end;
        SpecErrors := SpecErrors or AP_FP_Greater(Err,NonStrictErrTol);
        
        //
        // Test Gauss-Hermite
        //
        Err := 0;
        GQGenerateGaussHermite(N, Info, X, W);
        if Info>0 then
        begin
            BuildGaussHermiteQuadrature(N, X2, W2);
            I:=0;
            while I<=N-1 do
            begin
                Err := Max(Err, AbsReal(X[I]-X2[I]));
                Err := Max(Err, AbsReal(W[I]-W2[I]));
                Inc(I);
            end;
        end
        else
        begin
            SpecErrors := True;
        end;
        SpecErrors := SpecErrors or AP_FP_Greater(Err,NonStrictErrTol);
        Inc(N);
    end;
    
    //
    // end
    //
    WasErrors := RecErrors or SpecErrors;
    if  not Silent then
    begin
        Write(Format('TESTING GAUSS QUADRATURES'#13#10'',[]));
        Write(Format('FINAL RESULT:                             ',[]));
        if WasErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* SPECIAL CASES (LEGENDRE/JACOBI/..)      ',[]));
        if SpecErrors then
        begin
            Write(Format('FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('OK'#13#10'',[]));
        end;
        Write(Format('* RECURRENCE-BASED:                       ',[]));
        if RecErrors then
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
Gauss-Hermite, another variant
*************************************************************************)
procedure BuildGaussHermiteQuadrature(n : AlglibInteger;
     var x : TReal1DArray;
     var w : TReal1DArray);
var
    i : AlglibInteger;
    j : AlglibInteger;
    r : Double;
    r1 : Double;
    p1 : Double;
    p2 : Double;
    p3 : Double;
    dp3 : Double;
    pipm4 : Double;
    tmp : Double;
begin
    SetLength(x, n-1+1);
    SetLength(w, n-1+1);
    pipm4 := Power(Pi, -0.25);
    i:=0;
    while i<=(n+1) div 2-1 do
    begin
        if i=0 then
        begin
            r := Sqrt(2*n+1)-1.85575*Power(2*n+1, -AP_Double(1)/6);
        end
        else
        begin
            if i=1 then
            begin
                r := r-1.14*Power(n, 0.426)/r;
            end
            else
            begin
                if i=2 then
                begin
                    r := 1.86*r-0.86*x[0];
                end
                else
                begin
                    if i=3 then
                    begin
                        r := 1.91*r-0.91*x[1];
                    end
                    else
                    begin
                        r := 2*r-x[i-2];
                    end;
                end;
            end;
        end;
        repeat
            p2 := 0;
            p3 := pipm4;
            j:=0;
            while j<=n-1 do
            begin
                p1 := p2;
                p2 := p3;
                p3 := p2*r*Sqrt(AP_Double(2)/(j+1))-p1*Sqrt(AP_Double(j)/(j+1));
                Inc(j);
            end;
            dp3 := Sqrt(2*j)*p2;
            r1 := r;
            r := r-p3/dp3;
        until AP_FP_Less(AbsReal(r-r1),MachineEpsilon*(1+AbsReal(r))*100);
        x[i] := r;
        w[i] := 2/(dp3*dp3);
        x[n-1-i] := -x[i];
        w[n-1-i] := w[i];
        Inc(i);
    end;
    i:=0;
    while i<=N-1 do
    begin
        j:=0;
        while j<=n-2-i do
        begin
            if AP_FP_Greater_Eq(x[j],x[j+1]) then
            begin
                Tmp := x[j];
                x[j] := x[j+1];
                x[j+1] := Tmp;
                Tmp := w[j];
                w[j] := w[j+1];
                w[j+1] := Tmp;
            end;
            Inc(j);
        end;
        Inc(i);
    end;
end;


(*************************************************************************
Maps:
    0   =>  -0.9
    1   =>  -0.5
    2   =>  -0.1
    3   =>   0.0
    4   =>  +0.1
    5   =>  +0.5
    6   =>  +0.9
    7   =>  +1.0
    8   =>  +1.5
    9   =>  +2.0
*************************************************************************)
function MapKind(K : AlglibInteger):Double;
begin
    Result := 0;
    if K=0 then
    begin
        Result := -0.9;
    end;
    if K=1 then
    begin
        Result := -0.5;
    end;
    if K=2 then
    begin
        Result := -0.1;
    end;
    if K=3 then
    begin
        Result := 0.0;
    end;
    if K=4 then
    begin
        Result := +0.1;
    end;
    if K=5 then
    begin
        Result := +0.5;
    end;
    if K=6 then
    begin
        Result := +0.9;
    end;
    if K=7 then
    begin
        Result := +1.0;
    end;
    if K=8 then
    begin
        Result := +1.5;
    end;
    if K=9 then
    begin
        Result := +2.0;
    end;
end;


(*************************************************************************
Gauss-Legendre, another variant
*************************************************************************)
procedure BuildGaussLegendreQuadrature(n : AlglibInteger;
     var x : TReal1DArray;
     var w : TReal1DArray);
var
    i : AlglibInteger;
    j : AlglibInteger;
    r : Double;
    r1 : Double;
    p1 : Double;
    p2 : Double;
    p3 : Double;
    dp3 : Double;
    Tmp : Double;
begin
    SetLength(x, n-1+1);
    SetLength(w, n-1+1);
    i:=0;
    while i<=(n+1) div 2-1 do
    begin
        r := cos(Pi*(4*i+3)/(4*n+2));
        repeat
            p2 := 0;
            p3 := 1;
            j:=0;
            while j<=n-1 do
            begin
                p1 := p2;
                p2 := p3;
                p3 := ((2*j+1)*r*p2-j*p1)/(j+1);
                Inc(j);
            end;
            dp3 := n*(r*p3-p2)/(r*r-1);
            r1 := r;
            r := r-p3/dp3;
        until AP_FP_Less(AbsReal(r-r1),MachineEpsilon*(1+AbsReal(r))*100);
        x[i] := r;
        x[n-1-i] := -r;
        w[i] := 2/((1-r*r)*dp3*dp3);
        w[n-1-i] := 2/((1-r*r)*dp3*dp3);
        Inc(i);
    end;
    i:=0;
    while i<=N-1 do
    begin
        j:=0;
        while j<=n-2-i do
        begin
            if AP_FP_Greater_Eq(x[j],x[j+1]) then
            begin
                Tmp := x[j];
                x[j] := x[j+1];
                x[j+1] := Tmp;
                Tmp := w[j];
                w[j] := w[j+1];
                w[j+1] := Tmp;
            end;
            Inc(j);
        end;
        Inc(i);
    end;
end;


(*************************************************************************
Gauss-Jacobi, another variant
*************************************************************************)
procedure BuildGaussJacobiQuadrature(n : AlglibInteger;
     Alpha : Double;
     Beta : Double;
     var x : TReal1DArray;
     var w : TReal1DArray);
var
    i : AlglibInteger;
    j : AlglibInteger;
    r : Double;
    r1 : Double;
    t1 : Double;
    t2 : Double;
    t3 : Double;
    p1 : Double;
    p2 : Double;
    p3 : Double;
    pp : Double;
    an : Double;
    bn : Double;
    a : Double;
    b : Double;
    c : Double;
    tmpsgn : Double;
    Tmp : Double;
    alfbet : Double;
    temp : Double;
    its : AlglibInteger;
begin
    SetLength(x, n-1+1);
    SetLength(w, n-1+1);
    i:=0;
    while i<=n-1 do
    begin
        if i=0 then
        begin
            an := alpha/n;
            bn := beta/n;
            t1 := (1+alpha)*(2.78/(4+n*n)+0.768*an/n);
            t2 := 1+1.48*an+0.96*bn+0.452*an*an+0.83*an*bn;
            r := (t2-t1)/t2;
        end
        else
        begin
            if i=1 then
            begin
                t1 := (4.1+alpha)/((1+alpha)*(1+0.156*alpha));
                t2 := 1+0.06*(n-8)*(1+0.12*alpha)/n;
                t3 := 1+0.012*beta*(1+0.25*AbsReal(alpha))/n;
                r := r-t1*t2*t3*(1-r);
            end
            else
            begin
                if i=2 then
                begin
                    t1 := (1.67+0.28*alpha)/(1+0.37*alpha);
                    t2 := 1+0.22*(n-8)/n;
                    t3 := 1+8*beta/((6.28+beta)*n*n);
                    r := r-t1*t2*t3*(x[0]-r);
                end
                else
                begin
                    if i<n-2 then
                    begin
                        r := 3*x[i-1]-3*x[i-2]+x[i-3];
                    end
                    else
                    begin
                        if i=n-2 then
                        begin
                            t1 := (1+0.235*beta)/(0.766+0.119*beta);
                            t2 := 1/(1+0.639*(n-4)/(1+0.71*(n-4)));
                            t3 := 1/(1+20*alpha/((7.5+alpha)*n*n));
                            r := r+t1*t2*t3*(r-x[i-2]);
                        end
                        else
                        begin
                            if i=n-1 then
                            begin
                                t1 := (1+0.37*beta)/(1.67+0.28*beta);
                                t2 := 1/(1+0.22*(n-8)/n);
                                t3 := 1/(1+8*alpha/((6.28+alpha)*n*n));
                                r := r+t1*t2*t3*(r-x[i-2]);
                            end;
                        end;
                    end;
                end;
            end;
        end;
        alfbet := alpha+beta;
        repeat
            temp := 2+alfbet;
            p1 := (alpha-beta+temp*r)*0.5;
            p2 := 1;
            j:=2;
            while j<=n do
            begin
                p3 := p2;
                p2 := p1;
                temp := 2*j+alfbet;
                a := 2*j*(j+alfbet)*(temp-2);
                b := (temp-1)*(alpha*alpha-beta*beta+temp*(temp-2)*r);
                c := 2*(j-1+alpha)*(j-1+beta)*temp;
                p1 := (b*p2-c*p3)/a;
                Inc(j);
            end;
            pp := (n*(alpha-beta-temp*r)*p1+2*(n+alpha)*(n+beta)*p2)/(temp*(1-r*r));
            r1 := r;
            r := r1-p1/pp;
        until AP_FP_Less(AbsReal(r-r1),MachineEpsilon*(1+AbsReal(r))*100);
        x[i] := r;
        w[i] := exp(LnGamma(alpha+n, tmpsgn)+LnGamma(beta+n, tmpsgn)-LnGamma(n+1, tmpsgn)-LnGamma(n+alfbet+1, tmpsgn))*temp*Power(2, alfbet)/(pp*p2);
        Inc(i);
    end;
    i:=0;
    while i<=N-1 do
    begin
        j:=0;
        while j<=n-2-i do
        begin
            if AP_FP_Greater_Eq(x[j],x[j+1]) then
            begin
                Tmp := x[j];
                x[j] := x[j+1];
                x[j+1] := Tmp;
                Tmp := w[j];
                w[j] := w[j+1];
                w[j+1] := Tmp;
            end;
            Inc(j);
        end;
        Inc(i);
    end;
end;


(*************************************************************************
Gauss-Laguerre, another variant
*************************************************************************)
procedure BuildGaussLaguerreQuadrature(n : AlglibInteger;
     Alpha : Double;
     var x : TReal1DArray;
     var w : TReal1DArray);
var
    i : AlglibInteger;
    j : AlglibInteger;
    r : Double;
    r1 : Double;
    p1 : Double;
    p2 : Double;
    p3 : Double;
    dp3 : Double;
    tsg : Double;
    tmp : Double;
begin
    SetLength(x, n-1+1);
    SetLength(w, n-1+1);
    i:=0;
    while i<=n-1 do
    begin
        if i=0 then
        begin
            r := (1+Alpha)*(3+0.92*Alpha)/(1+2.4*n+1.8*Alpha);
        end
        else
        begin
            if i=1 then
            begin
                r := r+(15+6.25*Alpha)/(1+0.9*Alpha+2.5*n);
            end
            else
            begin
                r := r+((1+2.55*(i-1))/(1.9*(i-1))+1.26*(i-1)*Alpha/(1+3.5*(i-1)))/(1+0.3*Alpha)*(r-x[i-2]);
            end;
        end;
        repeat
            p2 := 0;
            p3 := 1;
            j:=0;
            while j<=n-1 do
            begin
                p1 := p2;
                p2 := p3;
                p3 := ((-r+2*j+alpha+1)*p2-(j+Alpha)*p1)/(j+1);
                Inc(j);
            end;
            dp3 := (n*p3-(n+Alpha)*p2)/r;
            r1 := r;
            r := r-p3/dp3;
        until AP_FP_Less(AbsReal(r-r1),MachineEpsilon*(1+AbsReal(r))*100);
        x[i] := r;
        w[i] := -Exp(LnGamma(Alpha+n, tsg)-LnGamma(n, tsg))/(dp3*n*p2);
        Inc(i);
    end;
    i:=0;
    while i<=N-1 do
    begin
        j:=0;
        while j<=n-2-i do
        begin
            if AP_FP_Greater_Eq(x[j],x[j+1]) then
            begin
                Tmp := x[j];
                x[j] := x[j+1];
                x[j+1] := Tmp;
                Tmp := w[j];
                w[j] := w[j+1];
                w[j+1] := Tmp;
            end;
            Inc(j);
        end;
        Inc(i);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testgq_test_silent():Boolean;
begin
    Result := TestGQunit(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testgq_test():Boolean;
begin
    Result := TestGQunit(False);
end;


end.