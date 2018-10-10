
program _test;
uses Sysutils, testsvdunit;

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
        if not testsvdunit_test_silent() then
            raise Exception.Create('');
    except on E: Exception do 
        begin
            WriteLn('svd                              FAILED(seed=',MySeed,')');
            Halt(1);
        end;
    end;
    WriteLn('svd                              OK');
    Halt(0);
end.
