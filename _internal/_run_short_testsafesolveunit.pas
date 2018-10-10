
program _test;
uses Sysutils, testsafesolveunit;

var
    MySeed: Cardinal;
begin
    if StrToIntDef(ParamStr(1),0)=0 then
    begin
        Randomize();
        MySeed:=RandSeed;
    end
    else
        MySeed:=StrToIntDef(ParamStr(1),0);
    RandSeed:=MySeed;
    try 
        if not testsafesolveunit_test_silent() then
            raise Exception.Create('');
    except on E: Exception do 
        begin
            WriteLn('safesolve                        FAILED(seed=',MySeed,')');
            Halt(1);
        end;
    end;
    WriteLn('safesolve                        OK');
    Halt(0);
end.
