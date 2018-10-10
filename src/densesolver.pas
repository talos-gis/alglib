(*************************************************************************
Copyright (c) 2007-2008, Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************)
unit densesolver;
interface
uses Math, Sysutils, Ap, reflections, bidiagonal, qr, lq, blas, rotations, bdsvd, svd, lu, trlinsolve, rcond, tsort, xblas;

type
DenseSolverReport = record
    R1 : Double;
    RInf : Double;
end;


DenseSolverLSReport = record
    R2 : Double;
    CX : TReal2DArray;
    N : AlglibInteger;
    K : AlglibInteger;
end;



procedure RMatrixSolveM(const A : TReal2DArray;
     N : AlglibInteger;
     const B : TReal2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);
procedure RMatrixSolveLS(const A : TReal2DArray;
     NRows : AlglibInteger;
     NCols : AlglibInteger;
     const B : TReal1DArray;
     Threshold : Double;
     var Info : AlglibInteger;
     var Rep : DenseSolverLSReport;
     var X : TReal1DArray);
procedure RMatrixSolve(const A : TReal2DArray;
     N : AlglibInteger;
     const B : TReal1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal1DArray);

implementation

function DenseSolverRFSMax(N : AlglibInteger;
     R1 : Double;
     RInf : Double):AlglibInteger;forward;
function DenseSolverRFSMaxV2(N : AlglibInteger;
     R2 : Double):AlglibInteger;forward;


(*************************************************************************
Dense solver.

This  subroutine  solves  a  system  A*X=B,  where A is NxN non-denegerate
real matrix, X and B are NxM real matrices.

Additional features include:
* automatic detection of degenerate cases
* iterative improvement

INPUT PARAMETERS
    A       -   array[0..N-1,0..N-1], system matrix
    N       -   size of A
    B       -   array[0..N-1,0..M-1], right part
    M       -   size of right part
    
OUTPUT PARAMETERS
    Info    -   return code:
                * -3    if A is singular, or VERY close to singular.
                        X is filled by zeros in such cases.
                * -1    if N<=0 or M<=0 was passed
                *  1    if task is solved (matrix A may be near  singular,
                        check R1/RInf parameters for condition numbers).
    Rep     -   solver report, see below for more info
    X       -   array[0..N-1,0..M-1], it contains:
                * solution of A*X=B if A is non-singular (well-conditioned
                  or ill-conditioned, but not very close to singular)
                * zeros,  if  A  is  singular  or  VERY  close to singular
                  (in this case Info=-3).

SOLVER REPORT

Subroutine sets following fields of the Rep structure:
* R1        reciprocal of condition number: 1/cond(A), 1-norm.
* RInf      reciprocal of condition number: 1/cond(A), inf-norm.

SEE ALSO:
    DenseSolverR() - solves A*x = b, where x and b are Nx1 matrices.

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixSolveM(const A : TReal2DArray;
     N : AlglibInteger;
     const B : TReal2DArray;
     M : AlglibInteger;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal2DArray);
var
    I : AlglibInteger;
    J : AlglibInteger;
    K : AlglibInteger;
    RFS : AlglibInteger;
    NRFS : AlglibInteger;
    P : TInteger1DArray;
    XC : TReal1DArray;
    Y : TReal1DArray;
    BC : TReal1DArray;
    XA : TReal1DArray;
    XB : TReal1DArray;
    TX : TReal1DArray;
    DA : TReal2DArray;
    V : Double;
    VErr : Double;
    SmallErr : Boolean;
    TerminateNextTime : Boolean;
    i_ : AlglibInteger;
begin
    
    //
    // prepare: check inputs, allocate space...
    //
    if (N<=0) or (M<=0) then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(DA, N, N);
    SetLength(X, N, M);
    SetLength(Y, N);
    SetLength(XC, N);
    SetLength(BC, N);
    SetLength(TX, N+1);
    SetLength(XA, N+1);
    SetLength(XB, N+1);
    
    //
    // factorize matrix, test for exact/near singularity
    //
    I:=0;
    while I<=N-1 do
    begin
        APVMove(@DA[I][0], 0, N-1, @A[I][0], 0, N-1);
        Inc(I);
    end;
    RMatrixLU(DA, N, N, P);
    Rep.R1 := RMatrixLURCond1(DA, N);
    Rep.RInf := RMatrixLURCondInf(DA, N);
    if AP_FP_Less(Rep.R1,10*MachineEpsilon) or AP_FP_Less(Rep.RInf,10*MachineEpsilon) then
    begin
        I:=0;
        while I<=N-1 do
        begin
            J:=0;
            while J<=M-1 do
            begin
                X[I,J] := 0;
                Inc(J);
            end;
            Inc(I);
        end;
        Rep.R1 := 0;
        Rep.RInf := 0;
        Info := -3;
        Exit;
    end;
    Info := 1;
    
    //
    // solve
    //
    K:=0;
    while K<=M-1 do
    begin
        
        //
        // First, non-iterative part of solution process:
        // * pivots
        // * L*y = b
        // * U*x = y
        //
        for i_ := 0 to N-1 do
        begin
            BC[i_] := B[i_,K];
        end;
        I:=0;
        while I<=N-1 do
        begin
            if P[I]<>I then
            begin
                V := BC[I];
                BC[I] := BC[P[I]];
                BC[P[I]] := V;
            end;
            Inc(I);
        end;
        Y[0] := BC[0];
        I:=1;
        while I<=N-1 do
        begin
            V := APVDotProduct(@DA[I][0], 0, I-1, @Y[0], 0, I-1);
            Y[I] := BC[I]-V;
            Inc(I);
        end;
        XC[N-1] := Y[N-1]/DA[N-1,N-1];
        I:=N-2;
        while I>=0 do
        begin
            V := APVDotProduct(@DA[I][0], I+1, N-1, @XC[0], I+1, N-1);
            XC[I] := (Y[I]-V)/DA[I,I];
            Dec(I);
        end;
        
        //
        // Iterative improvement of xc:
        // * calculate r = bc-A*xc using extra-precise dot product
        // * solve A*y = r
        // * update x:=x+r
        //
        // This cycle is executed until one of two things happens:
        // 1. maximum number of iterations reached
        // 2. last iteration decreased error to the lower limit
        //
        NRFS := DenseSolverRFSMax(N, Rep.R1, Rep.RInf);
        TerminateNextTime := False;
        RFS:=0;
        while RFS<=NRFS-1 do
        begin
            if TerminateNextTime then
            begin
                Break;
            end;
            
            //
            // generate right part
            //
            SmallErr := True;
            I:=0;
            while I<=N-1 do
            begin
                APVMove(@XA[0], 0, N-1, @A[I][0], 0, N-1);
                XA[N] := -1;
                APVMove(@XB[0], 0, N-1, @XC[0], 0, N-1);
                XB[N] := B[I,K];
                XDot(XA, XB, N+1, TX, V, VErr);
                BC[I] := -V;
                SmallErr := SmallErr and AP_FP_Less(AbsReal(V),4*VErr);
                Inc(I);
            end;
            if SmallErr then
            begin
                TerminateNextTime := True;
            end;
            
            //
            // solve
            //
            I:=0;
            while I<=N-1 do
            begin
                if P[I]<>I then
                begin
                    V := BC[I];
                    BC[I] := BC[P[I]];
                    BC[P[I]] := V;
                end;
                Inc(I);
            end;
            Y[0] := BC[0];
            I:=1;
            while I<=N-1 do
            begin
                V := APVDotProduct(@DA[I][0], 0, I-1, @Y[0], 0, I-1);
                Y[I] := BC[I]-V;
                Inc(I);
            end;
            TX[N-1] := Y[N-1]/DA[N-1,N-1];
            I:=N-2;
            while I>=0 do
            begin
                V := APVDotProduct(@DA[I][0], I+1, N-1, @TX[0], I+1, N-1);
                TX[I] := (Y[I]-V)/DA[I,I];
                Dec(I);
            end;
            
            //
            // update
            //
            APVAdd(@XC[0], 0, N-1, @TX[0], 0, N-1);
            Inc(RFS);
        end;
        
        //
        // Store xc
        //
        for i_ := 0 to N-1 do
        begin
            X[i_,K] := XC[i_];
        end;
        Inc(K);
    end;
end;


(*************************************************************************
Dense solver.

This subroutine finds solution of the linear system A*X=B with non-square,
possibly degenerate A.  System  is  solved in the least squares sense, and
general least squares solution  X = X0 + CX*y  which  minimizes |A*X-B| is
returned. If A is non-degenerate, solution in the  usual sense is returned

Additional features include:
* iterative improvement

INPUT PARAMETERS
    A       -   array[0..NRows-1,0..NCols-1], system matrix
    NRows   -   vertical size of A
    NCols   -   horizontal size of A
    B       -   array[0..NCols-1], right part
    Threshold-  a number in [0,1]. Singular values  beyond  Threshold  are
                considered  zero.  Set  it to 0.0, if you don't understand
                what it means, so the solver will choose good value on its
                own.
                
OUTPUT PARAMETERS
    Info    -   return code:
                * -4    SVD subroutine failed
                * -1    if NRows<=0 or NCols<=0 or Threshold<0 was passed
                *  1    if task is solved
    Rep     -   solver report, see below for more info
    X       -   array[0..N-1,0..M-1], it contains:
                * solution of A*X=B if A is non-singular (well-conditioned
                  or ill-conditioned, but not very close to singular)
                * zeros,  if  A  is  singular  or  VERY  close to singular
                  (in this case Info=-3).

SOLVER REPORT

Subroutine sets following fields of the Rep structure:
* R2        reciprocal of condition number: 1/cond(A), 2-norm.
* N         = NCols
* K         dim(Null(A))
* CX        array[0..N-1,0..K-1], kernel of A.
            Columns of CX store such vectors that A*CX[i]=0.

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixSolveLS(const A : TReal2DArray;
     NRows : AlglibInteger;
     NCols : AlglibInteger;
     const B : TReal1DArray;
     Threshold : Double;
     var Info : AlglibInteger;
     var Rep : DenseSolverLSReport;
     var X : TReal1DArray);
var
    SV : TReal1DArray;
    U : TReal2DArray;
    VT : TReal2DArray;
    RP : TReal1DArray;
    UTB : TReal1DArray;
    SUTB : TReal1DArray;
    Tmp : TReal1DArray;
    TA : TReal1DArray;
    TX : TReal1DArray;
    Buf : TReal1DArray;
    W : TReal1DArray;
    I : AlglibInteger;
    J : AlglibInteger;
    NSV : AlglibInteger;
    KernelIdx : AlglibInteger;
    V : Double;
    VErr : Double;
    SVDFailed : Boolean;
    ZeroA : Boolean;
    RFS : AlglibInteger;
    NRFS : AlglibInteger;
    TerminateNextTime : Boolean;
    SmallErr : Boolean;
    i_ : AlglibInteger;
begin
    if (NRows<=0) or (NCols<=0) or AP_FP_Less(Threshold,0) then
    begin
        Info := -1;
        Exit;
    end;
    if AP_FP_Eq(Threshold,0) then
    begin
        Threshold := 1000*MachineEpsilon;
    end;
    
    //
    // Factorize A first
    //
    SVDFailed :=  not RMatrixSVD(A, NRows, NCols, 1, 2, 2, SV, U, VT);
    ZeroA := AP_FP_Eq(SV[0],0);
    if SVDFailed or ZeroA then
    begin
        if SVDFailed then
        begin
            Info := -4;
        end
        else
        begin
            Info := 1;
        end;
        SetLength(X, NCols);
        I:=0;
        while I<=NCols-1 do
        begin
            X[I] := 0;
            Inc(I);
        end;
        Rep.N := NCols;
        Rep.K := NCols;
        SetLength(Rep.CX, NCols, NCols);
        I:=0;
        while I<=NCols-1 do
        begin
            J:=0;
            while J<=NCols-1 do
            begin
                if I=J then
                begin
                    Rep.CX[I,J] := 1;
                end
                else
                begin
                    Rep.CX[I,J] := 0;
                end;
                Inc(J);
            end;
            Inc(I);
        end;
        Rep.R2 := 0;
        Exit;
    end;
    NSV := Min(NCols, NRows);
    if NSV=NCols then
    begin
        Rep.R2 := SV[NSV-1]/SV[0];
    end
    else
    begin
        Rep.R2 := 0;
    end;
    Rep.N := NCols;
    Info := 1;
    
    //
    // Iterative improvement of xc combined with solution:
    // 1. xc = 0
    // 2. calculate r = bc-A*xc using extra-precise dot product
    // 3. solve A*y = r
    // 4. update x:=x+r
    // 5. goto 2
    //
    // This cycle is executed until one of two things happens:
    // 1. maximum number of iterations reached
    // 2. last iteration decreased error to the lower limit
    //
    SetLength(UTB, NSV);
    SetLength(SUTB, NSV);
    SetLength(X, NCols);
    SetLength(Tmp, NCols);
    SetLength(TA, NCols+1);
    SetLength(TX, NCols+1);
    SetLength(Buf, NCols+1);
    I:=0;
    while I<=NCols-1 do
    begin
        X[I] := 0;
        Inc(I);
    end;
    KernelIdx := NSV;
    I:=0;
    while I<=NSV-1 do
    begin
        if AP_FP_Less_Eq(SV[I],Threshold*SV[0]) then
        begin
            KernelIdx := I;
            Break;
        end;
        Inc(I);
    end;
    Rep.K := NCols-KernelIdx;
    NRFS := DenseSolverRFSMaxV2(NCols, Rep.R2);
    TerminateNextTime := False;
    SetLength(RP, NRows);
    RFS:=0;
    while RFS<=NRFS do
    begin
        if TerminateNextTime then
        begin
            Break;
        end;
        
        //
        // calculate right part
        //
        if RFS=0 then
        begin
            APVMove(@RP[0], 0, NRows-1, @B[0], 0, NRows-1);
        end
        else
        begin
            SmallErr := True;
            I:=0;
            while I<=NRows-1 do
            begin
                APVMove(@TA[0], 0, NCols-1, @A[I][0], 0, NCols-1);
                TA[NCols] := -1;
                APVMove(@TX[0], 0, NCols-1, @X[0], 0, NCols-1);
                TX[NCols] := B[I];
                XDot(TA, TX, NCols+1, Buf, V, VErr);
                RP[I] := -V;
                SmallErr := SmallErr and AP_FP_Less(AbsReal(V),4*VErr);
                Inc(I);
            end;
            if SmallErr then
            begin
                TerminateNextTime := True;
            end;
        end;
        
        //
        // solve A*dx = rp
        //
        I:=0;
        while I<=NCols-1 do
        begin
            Tmp[I] := 0;
            Inc(I);
        end;
        I:=0;
        while I<=NSV-1 do
        begin
            UTB[I] := 0;
            Inc(I);
        end;
        I:=0;
        while I<=NRows-1 do
        begin
            V := RP[I];
            APVAdd(@UTB[0], 0, NSV-1, @U[I][0], 0, NSV-1, V);
            Inc(I);
        end;
        I:=0;
        while I<=NSV-1 do
        begin
            if I<KernelIdx then
            begin
                SUTB[I] := UTB[I]/SV[I];
            end
            else
            begin
                SUTB[I] := 0;
            end;
            Inc(I);
        end;
        I:=0;
        while I<=NSV-1 do
        begin
            V := SUTB[I];
            APVAdd(@Tmp[0], 0, NCols-1, @VT[I][0], 0, NCols-1, V);
            Inc(I);
        end;
        
        //
        // update x:  x:=x+dx
        //
        APVAdd(@X[0], 0, NCols-1, @Tmp[0], 0, NCols-1);
        Inc(RFS);
    end;
    
    //
    // fill CX
    //
    if Rep.K>0 then
    begin
        SetLength(Rep.CX, NCols, Rep.K);
        I:=0;
        while I<=Rep.K-1 do
        begin
            for i_ := 0 to NCols-1 do
            begin
                Rep.CX[i_,I] := VT[KernelIdx+I,i_];
            end;
            Inc(I);
        end;
    end;
end;


(*************************************************************************
Dense solver.

Similar to RMatrixSolveM() but solves task with one right part  (where b/x
are vectors, not matrices).

See RMatrixSolveM()  description  for  more  information  about subroutine
parameters.

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************)
procedure RMatrixSolve(const A : TReal2DArray;
     N : AlglibInteger;
     const B : TReal1DArray;
     var Info : AlglibInteger;
     var Rep : DenseSolverReport;
     var X : TReal1DArray);
var
    BM : TReal2DArray;
    XM : TReal2DArray;
    i_ : AlglibInteger;
begin
    if N<=0 then
    begin
        Info := -1;
        Exit;
    end;
    SetLength(BM, N, 1);
    for i_ := 0 to N-1 do
    begin
        BM[i_,0] := B[i_];
    end;
    RMatrixSolveM(A, N, BM, 1, Info, Rep, XM);
    SetLength(X, N);
    for i_ := 0 to N-1 do
    begin
        X[i_] := XM[i_,0];
    end;
end;


(*************************************************************************
Internal subroutine.
Returns maximum count of RFS iterations as function of:
1. machine epsilon
2. task size.
3. condition number
*************************************************************************)
function DenseSolverRFSMax(N : AlglibInteger;
     R1 : Double;
     RInf : Double):AlglibInteger;
begin
    Result := 2;
end;


(*************************************************************************
Internal subroutine.
Returns maximum count of RFS iterations as function of:
1. machine epsilon
2. task size.
3. norm-2 condition number
*************************************************************************)
function DenseSolverRFSMaxV2(N : AlglibInteger; R2 : Double):AlglibInteger;
begin
    Result := DenseSolverRFSMax(N, 0, 0);
end;


end.