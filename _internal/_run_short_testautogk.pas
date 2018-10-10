
program _test;
uses Sysutils, testautogk;

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
        if not testautogk_test_silent() then
            raise Exception.Create('');
    except on E: Exception do 
        begin
            WriteLn('autogk                           FAILED(seed=',MySeed,')');
            Halt(1);
        end;
    end;
    WriteLn('autogk                           OK');
    Halt(0);
end.
