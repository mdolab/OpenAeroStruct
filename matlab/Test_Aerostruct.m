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
            % load Python from virtual environment with OpenMDAO 1.7.4 installed
            pyversion;
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
        function test_input_validation(testCase)
            prob_dict = struct;
            prob_dict.type = 'aerostruct';
            prob_dict.optimize = false;
            prob_dict.record_db = false;
            prob_dict.with_viscous = true;
            prob_dict.cg = [30, 0, 5];
            
            OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
            
            surf_dict = struct;
            surf_dict.name = 'wing';
            surf_dict.num_y = 11;
            surf_dict.num_x = 3;
            surf_dict.wing_type = 'CRM';
            surf_dict.CD0 = 0.015;
            surf_dict.symmetry = true;
            surf_dict.num_twist_cp = 2;
            surf_dict.num_thickness_cp = 2;
            OAS_prob.add_surface(surf_dict);
            
            surf_dict = struct;
            surf_dict.name = 'tail';
            surf_dict.num_y = 7;
            surf_dict.num_x = 2;
            surf_dict.span = 20;
            surf_dict.root_chord = 5;
            surf_dict.wing_type = 'rect';
            surf_dict.offset = [50, 0, 5];
            surf_dict.twist_cp = -9.5;
            OAS_prob.add_surface(surf_dict);
            
            
            OAS_prob.setup();
            out = struct(OAS_prob.run(pyargs('matlab',true)));
            CM = np2mat(out.CM);
            testCase.verifyEqual(out.wing_CL,0.836561056638825,'AbsTol',1e-5);
            testCase.verifyEqual(out.wing_failure,0.425301457001547,'AbsTol',1e-5);
            testCase.verifyEqual(out.fuelburn,97377.16195092892,'AbsTol',1e-2);
            testCase.verifyEqual(CM(2),-0.005158355060936,'AbsTol',1e-2);
            
        end
    end
    
    methods (Test, TestTags={'Aerostruct'})
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
        function test_aerostruct_analysis_symmetry(testCase)
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
        function test_aerostruct_analysis_set_variables(testCase)
            prob_dict = struct;
            prob_dict.type = 'aerostruct';
            prob_dict.with_viscous = true;
            prob_dict.optimize = false;
            prob_dict.record_db = false;  % using sqlitedict locks a process
            prob_dict.print_level = 0;
            prob_dict.alpha = 0.;
            
            OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
            
            % Create a dictionary to store options about the surface
            surf_dict = struct;
            surf_dict.name = 'wing';
            surf_dict.num_y = 7;
            surf_dict.num_x = 2;
            surf_dict.wing_type = 'CRM';
            surf_dict.CD0 = 0.015;
            surf_dict.symmetry = true;
            surf_dict.num_twist_cp = 2;
            surf_dict.num_thickness_cp = 2;
            surf_dict.num_chord_cp = 1;
            surf_dict.exact_failure_constraint = true;
            surf_dict.span_cos_spacing = 0.5;
            % Add the specified wing surface to the problem
            OAS_prob.add_surface(surf_dict);

            % Multiple lifting surfaces
            surf_dict = struct;
            surf_dict.name = 'tail';
            surf_dict.num_y = 7;
            surf_dict.num_x = 2;
            surf_dict.span = 20.;
            surf_dict.root_chord = 5.;
            surf_dict.wing_type = 'rect';
            surf_dict.offset = [50., 0., 5.];
            surf_dict.twist_cp = -9.5;
            surf_dict.exact_failure_constraint = true;
            OAS_prob.add_surface(surf_dict)
            
            OAS_prob.add_desvar('wing.twist_cp');
            OAS_prob.add_desvar('wing.thickness_cp');
            OAS_prob.add_desvar('wing.taper');
            OAS_prob.add_desvar('wing.chord_cp');

            OAS_prob.setup()
            % Actually run the problem
            input = {'wing.twist_cp',[12.803738284992180 14.737846154728121],...
                'wing.thickness_cp',[0.037776846454264, 0.071832717954386],...
                'wing.taper',0.2,'wing.chord_cp',0.9,...
                'matlab',true};
            output = struct(OAS_prob.run(pyargs(input{:})));
            testCase.verifyEqual(output.fuelburn,101898.5636 ,'AbsTol',1e-1); 
        end
    end
%         function test_aerostruct_symmetry_deriv(testCase)
%             prob_dict = struct;
%             prob_dict.type = 'aerostruct';
%             prob_dict.optimize = false;
%             prob_dict.record_db = false;
%             
%             OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
%             
%             surf_dict = struct;
%             surf_dict.symmetry = true;
%             surf_dict.num_y = int64(13);
%             surf_dict.num_x = int64(2);
%             surf_dict.wing_type = 'CRM';
%             surf_dict.CD0 = 0.015;
%             
%             OAS_prob.add_surface(surf_dict);
%             OAS_prob.setup();
%             out = struct(OAS_prob.run(pyargs('matlab',true)));
%     end
    
    methods (Test, TestTags = {'Aerostruct','Fortran'})
        function test_aerostruct_optimization(testCase)
            prob_dict = struct;
            prob_dict.type = 'aerostruct';
            prob_dict.optimize = true;
            prob_dict.with_viscous = true;
            prob_dict.record_db = false;
            
            OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
            
            surf_dict = struct;
            surf_dict.num_y = int64(7);
            surf_dict.num_x = int64(2);
            surf_dict.wing_type = 'CRM';
            surf_dict.CD0 = 0.015;
            surf_dict.symmetry = false;
            surf_dict.num_twist_cp = int64(2);
            surf_dict.num_thickness_cp = int64(2);
            
            OAS_prob.add_surface(surf_dict);
            
            OAS_prob.add_desvar('wing.twist_cp',pyargs('lower',-15.,'upper',15.));
            OAS_prob.add_desvar('wing.thickness_cp',pyargs('lower',0.01,'upper',0.5,'scaler',1e2));
            OAS_prob.add_constraint('wing_perf.failure',pyargs('upper',0.));
            OAS_prob.add_constraint('wing_perf.thickness_intersects',pyargs('upper',0.));
            OAS_prob.add_desvar('alpha',pyargs('lower',-10.,'upper',10.));
            OAS_prob.add_constraint('L_equals_W',pyargs('equals',0.));
            OAS_prob.add_objective('fuelburn',pyargs('scaler',1e-5));
            
            OAS_prob.setup();
            out = struct(OAS_prob.run(pyargs('matlab',true)));
            CM = np2mat(out.CM);
            testCase.verifyEqual(out.fuelburn,96889.255792361335,'AbsTol',1e0);
            testCase.verifyEqual(out.wing_failure, 0., 'AbsTol',1e-4);
            testCase.verifyEqual(CM(2), -0.14194155955058388, 'AbsTol',1e-2);
        end
        function test_aerostruct_optimization_symmetry(testCase)
            prob_dict = struct;
            prob_dict.type = 'aerostruct';
            prob_dict.optimize = true;
            prob_dict.with_viscous = true;
            prob_dict.record_db = false;
            
            OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
            
            surf_dict = struct;
            surf_dict.num_y = int64(7);
            surf_dict.num_x = int64(3);
            surf_dict.wing_type = 'CRM';
            surf_dict.CD0 = 0.015;
            surf_dict.symmetry = true;
            surf_dict.num_twist_cp = int64(2);
            surf_dict.num_thickness_cp = int64(2);
            
            OAS_prob.add_surface(surf_dict);
            
            OAS_prob.add_desvar('wing.twist_cp',pyargs('lower',-15.,'upper',15.));
            OAS_prob.add_desvar('wing.thickness_cp',pyargs('lower',0.01,'upper',0.5,'scaler',1e2));
            OAS_prob.add_constraint('wing_perf.failure',pyargs('upper',0.));
            OAS_prob.add_constraint('wing_perf.thickness_intersects',pyargs('upper',0.));
            OAS_prob.add_desvar('alpha',pyargs('lower',-10.,'upper',10.));
            OAS_prob.add_constraint('L_equals_W',pyargs('equals',0.));
            OAS_prob.add_objective('fuelburn',pyargs('scaler',1e-4));
            
            OAS_prob.setup();
            out = struct(OAS_prob.run(pyargs('matlab',true)));
%             CM = np2mat(out.CM);
            testCase.verifyEqual(out.fuelburn,96077.224922178371,'AbsTol',1e0);
            testCase.verifyEqual(out.wing_failure, 0., 'AbsTol',1e-5);
        end
    end
    
    methods (Test, TestTags={'Aero'})
        function test_aero_analysis_flat_viscous_full(testCase)
           prob_dict = struct;
           prob_dict.type = 'aero';
           prob_dict.optimize = false;
           prob_dict.record_db = false;
           prob_dict.with_viscous = true;
           
           surf_dict = struct;
           surf_dict.symmetry = false;
           
           OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
           
           OAS_prob.add_surface(surf_dict);
           OAS_prob.setup()
           out = struct(OAS_prob.run(pyargs('matlab',true)));
           
           testCase.verifyEqual(out.wing_CL, .45655138, 'AbsTol', 1e-5);
           testCase.verifyEqual(out.wing_CD, 0.018942466133780547, 'AbsTol', 1e-5);
        end
        function test_aero_analysis_flat_side_by_side(testCase)
           prob_dict = struct;
           prob_dict.type = 'aero';
           prob_dict.optimize = false;
           prob_dict.record_db = false;
           
           OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
           
           surf_dict = struct;
           surf_dict.name = 'wing';
           surf_dict.span = 5.;
           surf_dict.num_y = 3;
           surf_dict.span_cos_spacing = 0;
           surf_dict.symmetry = false;
           surf_dict.offset = [0, -2.5, 0];
           OAS_prob.add_surface(surf_dict);
           
           surf_dict = struct;
           surf_dict.name = 'tail';
           surf_dict.span = 5.;
           surf_dict.num_y = 3;
           surf_dict.span_cos_spacing = 0;
           surf_dict.symmetry = false;
           surf_dict.offset = [0, 2.5, 0];
           OAS_prob.add_surface(surf_dict);
           
           OAS_prob.setup()
           out = struct(OAS_prob.run(pyargs('matlab',true)));
           
           testCase.verifyEqual(out.wing_CL, 0.46173591841167183, 'AbsTol', 1e-5);
           testCase.verifyEqual(out.tail_CL, 0.46173591841167183, 'AbsTol', 1e-5);
           testCase.verifyEqual(out.wing_CD, .005524603647, 'AbsTol', 1e-5);
           testCase.verifyEqual(out.tail_CD, .005524603647, 'AbsTol', 1e-5);
        end
    end
    
    methods (Test, TestTags={'Aero','Fortran'})
        function test_aero_optimization(testCase)
            % Need to use SLSQP here because SNOPT finds a different optimum
            prob_dict = struct;
            prob_dict.type = 'aero';
            prob_dict.optimize = true;
            prob_dict.record_db = false;
            prob_dict.optimizer = 'SLSQP';
            
            OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
            
            OAS_prob.add_surface();
            
            OAS_prob.add_desvar('wing.twist_cp',pyargs('lower',-10.,'upper',15.));
            OAS_prob.add_desvar('wing.sweep',pyargs('lower',10.,'upper',30.));
            OAS_prob.add_desvar('wing.dihedral',pyargs('lower',-10.,'upper',20.));
            OAS_prob.add_constraint('wing_perf.CL',pyargs('equals',0.5));
            OAS_prob.add_objective('wing_perf.CD', pyargs('scaler',1e4));
            
            OAS_prob.setup();
            out = struct(OAS_prob.run(pyargs('matlab',true)));
            testCase.verifyEqual(out.wing_CD,0.0049392534859265614,'AbsTol',1e-5);
        end    
    end
    
    methods (Test, TestTags={'Struct'})
        function test_struct_analysis(testCase)
            prob_dict = struct;
            prob_dict.type = 'struct';
            prob_dict.optimize = false;
            prob_dict.record_db = false;
            
            OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
            
            surf_dict = struct;
            surf_dict.symmetry = false;
            surf_dict.t_over_c = 0.15;
            
            OAS_prob.add_surface(surf_dict);

            OAS_prob.setup();
            out = struct(OAS_prob.run(pyargs('matlab',true)));
            testCase.verifyEqual(out.wing_structural_weight,988.13495481064024,'AbsTol',1e-3);
        end 
        function test_struct_analysis_symmetry(testCase)
            prob_dict = struct;
            prob_dict.type = 'struct';
            prob_dict.optimize = false;
            prob_dict.record_db = false;
            
            OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
            
            surf_dict = struct;
            surf_dict.symmetry = true;
            surf_dict.t_over_c = 0.15;
            
            OAS_prob.add_surface(surf_dict);

            OAS_prob.setup();
            out = struct(OAS_prob.run(pyargs('matlab',true)));
            testCase.verifyEqual(out.wing_structural_weight,988.13495481063956,'AbsTol',1e-3);
        end
        function test_struct_optimization_symmetry(testCase)
            prob_dict = struct;
            prob_dict.type = 'struct';
            prob_dict.optimize = true;
            prob_dict.record_db = false;
            
            OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
            
            surf_dict = struct;
            surf_dict.symmetry = true;
            surf_dict.t_over_c = 0.15;
            OAS_prob.add_surface(surf_dict);

            OAS_prob.add_desvar('wing.thickness_cp', pyargs('lower',0.001,'upper',0.25,'scaler',1e2));
            OAS_prob.add_constraint('wing.failure',pyargs('upper',0.));
            OAS_prob.add_constraint('wing.thickness_intersects',pyargs('upper',0.));
            OAS_prob.add_objective('wing.structural_weight',pyargs('scaler',1e-3));
            
            OAS_prob.setup();
            out = struct(OAS_prob.run(pyargs('matlab',true)));
            
            testCase.verifyEqual(out.wing_structural_weight,1144.8503583047038,'AbsTol',1e-2);
        end
        
    end
    
    methods (Test, TestTags={'Struct','Fortran'})
        function test_struct_optimization(testCase)
            prob_dict = struct;
            prob_dict.type = 'struct';
            prob_dict.optimize = true;
            prob_dict.record_db = false;
            
            OAS_prob = py.OpenAeroStruct.run_classes.OASProblem(prob_dict);
            
            surf_dict = struct;
            surf_dict.symmetry = false;
            surf_dict.t_over_c = 0.15;
            OAS_prob.add_surface(surf_dict);

            OAS_prob.add_desvar('wing.thickness_cp', pyargs('lower',0.001,'upper',0.25,'scaler',1e2));
            OAS_prob.add_constraint('wing.failure',pyargs('upper',0.));
            OAS_prob.add_constraint('wing.thickness_intersects',pyargs('upper',0.));
            OAS_prob.add_objective('wing.structural_weight',pyargs('scaler',1e-3));
            
            OAS_prob.setup();
            out = struct(OAS_prob.run(pyargs('matlab',true)));
            
            testCase.verifyEqual(out.wing_structural_weight,1154.4491377169238,'AbsTol',1e-2);
        end
    end
end

