
program _test;
uses Sysutils, testconvunit;

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
        if not testconvunit_test_silent() then
            raise Exception.Create('');
    except on E: Exception do 
        begin
            WriteLn('conv                             FAILED(seed=',MySeed,')');
            Halt(1);
        end;
    end;
    WriteLn('conv                             OK');
    Halt(0);
end.
