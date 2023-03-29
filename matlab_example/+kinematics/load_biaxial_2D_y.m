classdef load_biaxial_2D_y
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lambda11
        lambda12
        lambda21
        lambda22
        p
        ft
        name
        free = 0
    end
    
    methods
        function obj = load_biaxial_2D_y(wave,period,amp,varargin)
            %UNTITLED4 Construct an instance of this class
            %   Detailed explanation goes here
            obj.name     = '2D Biaxial Stretch test';
            obj.lambda11 = amp(1);
            obj.lambda12 = amp(2);
            obj.lambda21 = amp(3);
            obj.lambda22 = amp(4);
            obj.p        = period;
            if not(isempty(varargin))
                obj.free     = varargin{1};
            end
            switch wave
                case 'SawCos'
                    obj.ft       = @saw_cos_wave;
                otherwise
                    error(['Unknown wave form ' wave ' requested'])
            end
            
        end
        
        function F = arg_func(obj,t)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            x = obj.ft(obj.p, t);
            F = eye(2);
            F(1,1) = obj.lambda11*x + F(1,1);
            F(1,2) = obj.lambda12*x;
            F(2,1) = obj.lambda21*x;
            F(2,2) = obj.lambda22*x + F(2,2);
            switch obj.free
                case 1
                    F(1,1) = 2.0 - F(2,2);
                case 2
                    F(2,2) = 2.0 - F(1,1);
            end            
        end
    end
end

