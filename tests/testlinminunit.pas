unit testlinminunit;
interface
uses Math, Sysutils, Ap, linmin;

function TestLinMin(Silent : Boolean):Boolean;
function testlinminunit_test_silent():Boolean;
function testlinminunit_test():Boolean;

implementation

function TestLinMin(Silent : Boolean):Boolean;
var
    WasErrors : Boolean;
begin
    WasErrors := False;
    if  not Silent then
    begin
        Write(Format('TESTING LINMIN'#13#10'',[]));
        if WasErrors then
        begin
            Write(Format('TEST FAILED'#13#10'',[]));
        end
        else
        begin
            Write(Format('TEST PASSED'#13#10'',[]));
        end;
        Write(Format(''#13#10''#13#10'',[]));
    end;
    Result :=  not WasErrors;
end;


(*************************************************************************
Silent unit test
*************************************************************************)
function testlinminunit_test_silent():Boolean;
begin
    Result := TestLinMin(True);
end;


(*************************************************************************
Unit test
*************************************************************************)
function testlinminunit_test():Boolean;
begin
    Result := TestLinMin(False);
end;


end.