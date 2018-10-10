
program _test;
uses Sysutils, testcdetunit;

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
        if not testcdetunit_test_silent() then
            raise Exception.Create('');
    except on E: Exception do 
        begin
            WriteLn('cdet                             FAILED(seed=',MySeed,')');
            Halt(1);
        end;
    end;
    WriteLn('cdet                             OK');
    Halt(0);
end.
