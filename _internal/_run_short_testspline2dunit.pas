
program _test;
uses Sysutils, testspline2dunit;

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
        if not testspline2dunit_test_silent() then
            raise Exception.Create('');
    except on E: Exception do 
        begin
            WriteLn('spline2d                         FAILED(seed=',MySeed,')');
            Halt(1);
        end;
    end;
    WriteLn('spline2d                         OK');
    Halt(0);
end.
