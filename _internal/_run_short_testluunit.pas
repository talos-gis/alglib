
program _test;
uses Sysutils, testluunit;

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
        if not testluunit_test_silent() then
            raise Exception.Create('');
    except on E: Exception do 
        begin
            WriteLn('lu                               FAILED(seed=',MySeed,')');
            Halt(1);
        end;
    end;
    WriteLn('lu                               OK');
    Halt(0);
end.