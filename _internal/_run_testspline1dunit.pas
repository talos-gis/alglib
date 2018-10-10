
program _test;
uses Sysutils, testspline1dunit;

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
        if not testspline1dunit_test() then
        begin
            WriteLn('SEED: ', MySeed);
            Halt(1);
        end;
    except on E: Exception do 
        begin
            WriteLn('SEED: ', MySeed);
            WriteLn('AP exception generated!');
            WriteLn('\"',E.Message,'\"');
            Halt(1);
        end;
    end;
    Halt(0);
end.
