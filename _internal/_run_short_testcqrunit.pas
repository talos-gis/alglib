
program _test;
uses Sysutils, testcqrunit;

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
        if not testcqrunit_test_silent() then
            raise Exception.Create('');
    except on E: Exception do 
        begin
            WriteLn('cqr                              FAILED(seed=',MySeed,')');
            Halt(1);
        end;
    end;
    WriteLn('cqr                              OK');
    Halt(0);
end.
