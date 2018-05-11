% Run Test Suite

suiteClass = matlab.unittest.TestSuite.fromClass(?Test_Aerostruct);
result = run(suiteClass);
disp(result)