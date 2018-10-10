
program _test;
uses Sysutils, testldaunit;

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
        if not testldaunit_test_silent() then
            raise Exception.Create('');
    except on E: Exception do 
        begin
            WriteLn('lda                              FAILED(seed=',MySeed,')');
            Halt(1);
        end;
    end;
    WriteLn('lda                              OK');
    Halt(0);
end.
