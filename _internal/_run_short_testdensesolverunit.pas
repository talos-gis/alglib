
program _test;
uses Sysutils, testdensesolverunit;

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
        if not testdensesolverunit_test_silent() then
            raise Exception.Create('');
    except on E: Exception do 
        begin
            WriteLn('densesolver                      FAILED(seed=',MySeed,')');
            Halt(1);
        end;
    end;
    WriteLn('densesolver                      OK');
    Halt(0);
end.
