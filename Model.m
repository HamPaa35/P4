classdef Model
    %MODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        p_start
        A
        phi
    end
    
    methods
        function obj = Model(p_start, A, phi)
            obj.p_start = p_start;
            obj.A = A;
            obj.phi = phi;
        end
    end
end

