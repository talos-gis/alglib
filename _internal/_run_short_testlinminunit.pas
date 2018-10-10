
program _test;
uses Sysutils, testlinminunit;

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
        if not testlinminunit_test_silent() then
            raise Exception.Create('');
    except on E: Exception do 
        begin
            WriteLn('linmin                           FAILED(seed=',MySeed,')');
            Halt(1);
        end;
    end;
    WriteLn('linmin                           OK');
    Halt(0);
end.
