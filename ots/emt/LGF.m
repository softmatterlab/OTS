classdef LGF < handle
    % LGF < handle : Logarithm of factorial
    %	LGF is a class to precalcualte the logarithms of factorials.
    %
    % LGF properties:
    %   J       -   maximum value for computing factorial
    %   d       -   digits
    %   fact    -   precalculated symbolic log(k!)
    %
    % LGF methods:
    %   LGF     -   constructor
    %   getfact	-   getting precalculated symbolic log(k!)
    %
    % See also CG.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01

    properties
        J       % maximum J
        d       % digits
        fact    % precalculated symbolic log(k!)
    end    
    methods
        function obj = LGF(varargin)
            % LGF() precalculates the logarithms of factorials log(J!)
            %   for J up to 1000. 
            %
            % LGF('PropertyName',PropertyValue) sets the property PropertyName to PropertyValue. 
            %   The properties listed below can be used:
            %       J	-   maximum J [default J=1000]
            %       d	-   digits [default d=100]
            %
            %  See also LGF, CG.
            
            % maximum J
            J = 1e+3;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'J')
                    J = varargin{n+1};
                end
            end

            % digits
            d = 1e+2;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'d')
                    d = varargin{n+1};
                end
            end
            
            obj.J = J;

            obj.d = d;

            dtmp = digits;
            digits(obj.d);
            log_k_fac = sym('log(k!)');
            obj.fact = zeros(1,J+1);
            for n = 0:1:J
                obj.fact(n+1) = vpa(subs(log_k_fac,'k',n));
            end
            digits(dtmp);
        end
        function fact = getfact(lgf,n)
            % GETFACT log(N!)
            %
            % FACT = GETFACT(LGF,N) retrieves the log(N!).
            %
            %  See also LGF, CG.

            fact = lgf.fact(n+1);
        end
    end
end