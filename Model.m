classdef Model
    %MODEL This class is for holding data for GMM and HMM model
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

