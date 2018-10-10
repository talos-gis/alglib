
program _test;
uses Sysutils, testmlpeunit;

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
        if not testmlpeunit_test_silent() then
            raise Exception.Create('');
    except on E: Exception do 
        begin
            WriteLn('mlpe                             FAILED(seed=',MySeed,')');
            Halt(1);
        end;
    end;
    WriteLn('mlpe                             OK');
    Halt(0);
end.
