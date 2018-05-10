classdef Test_Aerostruct < matlab.unittest.TestCase
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        origPath
        fortran_flag
    end
    
    methods (TestClassSetup)
        function LoadOpenAeroStruct(testCase)
            % load Python and import OpenAeroStruct
            try
                % for Unix systems
                py.sys.setdlopenflags(int32(10));  % Set RTLD_NOW and RTLD_DEEPBIND
            catch
            end
            % load Python from virtual environment with OpenMDAO 1.7.3 installed
            [~,~,isloaded] = pyversion;
            % import OpenAeroStruct python module
            OAS_PATH = py.os.path.abspath('../..');
            P = py.sys.path;
            if count(P,OAS_PATH) == 0
                insert(P,int32(0),OAS_PATH);
            end
            py.importlib.import_module('OpenAeroStruct');
            % set fortran flag
            testCase.fortran_flag = py.OpenAeroStruct.run_classes.fortran_flag;
            % seed random number generator
            rng(42);
        end
    end
    
    methods (Test)
        function test_aerostruct_analysis(testCase)
            prob_dict = struct;
            prob_dict.type = 'aerostruct';
            prob_dict.optimize = false;
            prob_dict.record_db = false;
            
            OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
            
            surf_dict = struct;
            surf_dict.num_y = int64(13);
            surf_dict.num_x = int64(2);
            surf_dict.wing_type = 'CRM';
            surf_dict.CD0 = 0.015;
            surf_dict.symmetry = false;
            
            OAS_prob.add_surface(surf_dict);
            OAS_prob.setup();
            out = struct(OAS_prob.run(pyargs('matlab',true)));
            CM = np2mat(out.CM);
            testCase.verifyEqual(out.wing_CL,0.6587983088529461,'AbsTol',1e-5);
            testCase.verifyEqual(out.wing_failure,0.13716279310143381,'AbsTol',1e-5);
            testCase.verifyEqual(out.fuelburn,55565.087226705218,'AbsTol',1e-2);
            testCase.verifyEqual(CM(2),-0.18836163204083048,'AbsTol',1e-2);
        end
        
        function test_aerostruct_symmetry(testCase)
            prob_dict = struct;
            prob_dict.type = 'aerostruct';
            prob_dict.optimize = false;
            prob_dict.record_db = false;
            
            OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
            
            surf_dict = struct;
            surf_dict.symmetry = true;
            surf_dict.num_y = int64(13);
            surf_dict.num_x = int64(2);
            surf_dict.wing_type = 'CRM';
            surf_dict.CD0 = 0.015;
            
            OAS_prob.add_surface(surf_dict);
            OAS_prob.setup();
            out = struct(OAS_prob.run(pyargs('matlab',true)));
            CM = np2mat(out.CM);
            testCase.verifyEqual(out.wing_CL,0.69060502679333224,'AbsTol',1e-5);
            testCase.verifyEqual(out.wing_failure,0.064759950520982532,'AbsTol',1e-5);
            testCase.verifyEqual(out.fuelburn,57109.065516474155,'AbsTol',1e-1);
            testCase.verifyEqual(CM(2),-0.19380236992046351,'AbsTol',1e-2);
        end
    end
    
    methods (Test, TestTags = {'Fortran'})
        function test_aero_optimization_chord_monotonic(testCase)
            prob_dict = struct;
            prob_dict.type = 'aero';
            prob_dict.optimize = true;
            prob_dict.record_db = false;
            prob_dict.with_viscous = false;
            
            OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
            
            surf_dict = struct;
            surf_dict.symmetry = true;
            surf_dict.num_y = int64(13);
            surf_dict.num_x = int64(2);
            surf_dict.wing_type = 'CRM';
            surf_dict.CD0 = 0.015;
            
            OAS_prob.add_surface(surf_dict);
            OAS_prob.setup();
            out = struct(OAS_prob.run(pyargs('matlab',true)));
            CM = np2mat(out.CM);
            testCase.verifyEqual(out.wing_CL,0.69060502679333224,'AbsTol',1e-5);
            testCase.verifyEqual(out.wing_failure,0.064759950520982532,'AbsTol',1e-5);
            testCase.verifyEqual(out.fuelburn,57109.065516474155,'AbsTol',1e-1);
            testCase.verifyEqual(CM(2),-0.19380236992046351,'AbsTol',1e-2);
        end
    end
    
end

