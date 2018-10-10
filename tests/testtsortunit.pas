unit testtsortunit;
interface
uses Math, Sysutils, Ap, tsort;

function TestSort(Silent : Boolean):Boolean;
function testtsortunit_test_silent():Boolean;
function testtsortunit_test():Boolean;

implementation

procedure Unset2D(var A : TComplex2DArray);forward;
procedure Unset1D(var A : TReal1DArray);forward;
procedure Unset1DI(var A : TInteger1DArray);forward;
procedure TestSortResults(const ASorted : TReal1DArray;
     const P1 : TInteger1DArray;
     const P2 : TInteger1DArray;
     const AOriginal : TReal1DArray;
     N : AlglibInteger;
     var WasErrors : Boolean);forward;


(*************************************************************************
Testing tag sort
*************************************************************************)
function TestSort(Silent : Boolean):Boolean;
var
    WasErrors : Boolean;
    N : AlglibInteger;
    I : AlglibInteger;
    Pass : AlglibInteger;
    PassCount : AlglibInteger;
    MaxN : AlglibInteger;
    A : TReal1DArray;
    A0 : TReal1DArray;
    A2 : TReal1DArray;
    P1 : TInteger1DArray;
    P2 : TInteger1DArray;
begin
    WasErrors := False;
    MaxN := 100;
    PassCount := 10;
    
    //
    // Test
    //
    N:=1;
    while N<=MaxN do
    begin
        Pass:=1;
        while Pass<=PassCount do
        begin
            
            //
            // (probably) distinct sort
            //
            Unset1DI(P1);
            Unset1DI(P2);
            SetLength(A, N-1+1);
            SetLength(A0, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                A[I] := 2*RandomReal-1;
                A0[I] := A[I];
                Inc(I);
            end;
            TagSort(A0, N, P1, P2);
            TestSortResults(A0, P1, P2, A, N, WasErrors);
            
            //
            // non-distinct sort
            //
            Unset1DI(P1);
            Unset1DI(P2);
            SetLength(A, N-1+1);
            SetLength(A0, N-1+1);
            I:=0;
            while I<=N-1 do
            begin
                A[I] := I div 2;
                A0[I] := A[I];
                Inc(I);
            end;
            TagSort(A0, N, P1, P2);
            TestSortResults(A0, P1, P2, A, N, WasErrors);
            Inc(Pass);
        end;
        Inc(N);
    end;
    
    //
    // report
    //
    if  not Silent then
    begin
        Write(Format('TESTING TAGSORT'#13#10'',[]));
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
Unsets 2D array.
*************************************************************************)
procedure Unset2D(var A : TComplex2DArray);
begin
    SetLength(A, 0+1, 0+1);
    A[0,0] := C_Complex(2*RandomReal-1);
end;


(*************************************************************************
Unsets 1D array.
*************************************************************************)
procedure Unset1D(var A : TReal1DArray);
begin
    SetLength(A, 0+1);
    A[0] := 2*RandomReal-1;
end;


(*************************************************************************
Unsets 1D array.
*************************************************************************)
procedure Unset1DI(var A : TInteger1DArray);
begin
    SetLength(A, 0+1);
    A[0] := RandomInteger(3)-1;
end;


procedure TestSortResults(const ASorted : TReal1DArray;
     const P1 : TInteger1DArray;
     const P2 : TInteger1DArray;
     const AOriginal : TReal1DArray;
     N : AlglibInteger;
     var WasErrors : Boolean);
var
    I : AlglibInteger;
    A2 : TReal1DArray;
    T : Double;
    F : TInteger1DArray;
begin
    SetLength(A2, N-1+1);
    SetLength(F, N-1+1);
    
    //
    // is set ordered?
    //
    I:=0;
    while I<=N-2 do
    begin
        WasErrors := WasErrors or AP_FP_Greater(ASorted[I],ASorted[I+1]);
        Inc(I);
    end;
    
    //
    // P1 correctness
    //
    I:=0;
    while I<=N-1 do
    begin
        WasErrors := WasErrors or AP_FP_Neq(ASorted[I],AOriginal[P1[I]]);
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        F[I] := 0;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        F[P1[I]] := F[P1[I]]+1;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        WasErrors := WasErrors or (F[I]<>1);
        Inc(I);
    end;
    
    //
    // P2 correctness
    //
    I:=0;
    while I<=N-1 do
    begin
        A2[I] := AOriginal[I];
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        if P2[I]<>I then
        begin
            T := A2[I];
            A2[I] := A2[P2[I]];
            A2[P2[I]] := T;
        end;
        Inc(I);
    end;
    I:=0;
    while I<=N-1 do
    begin
        WasErrors := WasErrors or AP_FP_Neq(ASorted[I],A2[I]);
        Inc(I);
    end;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testtsortunit_test_silent():Boolean;
begin
    Result := TestSort(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testtsortunit_test():Boolean;
begin
    Result := TestSort(False);
end;


end.