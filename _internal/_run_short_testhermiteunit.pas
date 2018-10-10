
program _test;
uses Sysutils, testhermiteunit;

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
        if not testhermiteunit_test_silent() then
            raise Exception.Create('');
    except on E: Exception do 
        begin
            WriteLn('hermite                          FAILED(seed=',MySeed,')');
            Halt(1);
        end;
    end;
    WriteLn('hermite                          OK');
    Halt(0);
end.
