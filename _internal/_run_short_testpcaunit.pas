
program _test;
uses Sysutils, testpcaunit;

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
        if not testpcaunit_test_silent() then
            raise Exception.Create('');
    except on E: Exception do 
        begin
            WriteLn('pca                              FAILED(seed=',MySeed,')');
            Halt(1);
        end;
    end;
    WriteLn('pca                              OK');
    Halt(0);
end.
