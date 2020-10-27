classdef Coefficients
    % Coefficients : Multipole amplitudes
    %   Coefficients C are multipole amplitude for incident and scattered field.
    %   They are used in, e.g., Coefficients, TMatrix, IncidentField,
    %   IncidentFieldPlaneWave.
    %
    % Coefficients properties:
    %   C                  -   multipole amplitudes
    %
    % Coefficients methods:
    %   Coefficients       -   constructor
    %   lmax               -   maximum number of coefficients 
    %   c                  -   coefficients corresponding to indices p,l,m
    %   c1                 -   coefficients corresponding to indices p=1,l,m
    %   c2                 -   coefficients corresponding to indices p=2,l,m
    %   uplus              -   unitary plus (+c = c)
    %   uminus             -   unitary mines (cm = -c)
    %   plus               -   summation of coefficients (c = c1+c2)
    %   minus              -   subtraction of coefficients (c = c1-c2)
    %   mtimes             -   scalar product (c = c1*c2)
    %   times              -   scalar product (c = c1.*c2)
    %   rdivide            -   division (c = c1/b)
    %   real               -   real part of coefficients
    %   imag               -   imaginary part of coefficients
    %   conj               -   connjugate of coefficients
    %
    % Coefficients static methods:
    %   index              -   index i corresponding to p, l, m
    %   plm                -   p, l, m corresponding to index i
    %
    % See also Coefficients, TMatrix.
    
    %   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01
    
    properties
        C  % multipole amplitudes
    end
    methods
        function obj = Coefficients(C)
            % COEFFICIENTS(C) constructs multipole coeficients.
            %
            % See also Coefficients, TMatrix.
            
            obj.C = sparse(reshape(C,numel(C),1));
        end
        function L = lmax(coeff)
            % LMAX Maximum number of coefficients 
            %
            % L = LMAX(COEFF) returns maximum number of coefficients. 
            %
            % See also Coefficients.

            L = floor(sqrt((length(coeff.C)-1)/2));
        end
        function coeff_plm = c(coeff,p,l,m)
            % C Coefficients for indices p, l, m 
            %
            % COEFF_PLM = C(COEFF,P,L,M) gives coefficients for indices P, L, M.
            %
            % See also Coefficients.

            coeff_plm = coeff.C(Coefficients.index(p,l,m));
        end
        function coeff_plm = c1(coeff,l,m)
            % C1 Coefficients for index p=1, l, m 
            %
            % COEFF_PLM = C1(COEFF,L,M) gives coefficients for indices P = 1 L, M.
            %
            % See also Coefficients.

            coeff_plm = coeff.C(Coefficients.index(1,l,m));
        end
        function coeff_plm = c2(coeff,l,m)
            % C2 Coefficients for index p=2, l, m 
            %
            % COEFF_PLM = C2(COEFF,L,M) gives coefficients for indices P = 2 L, M.
            %
            % See also Coefficients.

            coeff_plm = coeff.C(Coefficients.index(2,l,m));
        end
        function coeff = uplus(coeff)            
            % UPLUS Unitary plus (coefficients)
            %
            % C = UPLUS(C) Unitary plus (+C = C)
            %
            % See also Coefficients.
        end
        function coeff_m = uminus(coeff)
            % UMINUS Unitary minus (coefficients)
            %
            % COEFF_M = UMINUS(COEFF) Unitary minus (-COEFF_M)
            %
            % See also Coefficients.

            coeff_m = coeff;
            coeff_m.C = -coeff_m.C;
        end
        function coeff = plus(coeff1,coeff2)
            % PLUS Binary addition (coefficients)
            %
            % C = PLUS(C1,C2) Binary addition (C = C1+C2).
            %   The coefficient C is the sum of the coefficients C1 and C2.
            %
            % See also Coefficients.

            coeff = coeff1;
            coeff.C = coeff1.C+coeff2.C;
        end
        function coeff = minus(coeff1,coeff2)
            % MINUS Binary subtraction (coefficients)
            %
            % C = MINUS(C1,C2) Binary subtraction (C = C1-C2).
            %   The coefficient C is the difference of the coefficients C1 and C2.
            %
            % See also Coefficients.

            coeff = coeff1;
            coeff.C = coeff1.C-coeff2.C;
        end
        function m = mtimes(a,b)
            % MTIMES Scalar product (coefficients)
            %
            % C = MTIMES(A,B) scaler product of coefficients (C = C1*C2).
            %   The coefficients of C are the coefficients C1 multiplied by
            %   the coefficients C2.
            %
            % See also Coefficients.

            if isa(a,'Coefficients') && isa(b,'Coefficients')
                coeff1 = a;
                coeff2 = b;
                m = coeff1;
                m.C = coeff1.C.*coeff2.C;
            elseif isa(a,'Coefficients')
                coeff1 = a;
                m = coeff1;
                m.C = coeff1.C.*b;
            elseif isa(b,'Coefficients')
                coeff2 = b;
                m = coeff2;
                m.C = a.*coeff2.C;
            else
                m = a.*b;
            end
        end
        function m = times(a,b)
            % MTIMES Scalar product (coefficients)
            %
            % C = TIMES(C1,C2) scaler product of coefficients (C = C1.*C2).
            %   The coefficients of C are the coefficients C1 multiplied by
            %   the coefficients C2.
            %
            % See also Coefficients.

            m = mtimes(a,b)
        end
        function coeff_d = rdivide(coeff,b)
            % RDIVIDE Devision (coefficients)
            %
            % C_D = RDIVIDE(C,B) division of coefficients by a constnat (C_D = C/B).
            %   The coefficients of C_D are the the coefficients of  C divided by B.
            %
            % See also Coefficients.

            coeff_d = coeff;
            coeff_d.C = coeff.C./b;
        end
        function coeff2 = real(coeff)
            % REAL Real part of coefficients
            %
            % C2 = REAL(C) determines the real part of C.
            %
            % See also Coefficients.

            coeff2 = coeff;
            coeff2.C = real(coeff.C);
        end
        function coeff2 = imag(coeff)
            % IMAG Imaginary part of coefficients
            %
            % C2 = IMAG(C) determines the imaginary part of C.
            %
            % See also Coefficients.

            coeff2 = coeff;
            coeff2.C = imag(coeff.C);
        end
        function coeff2 = conj(coeff)
            % CONJ Conjugate of coefficients
            %
            % C2 = CONJ(C) determines the conjugate of COEFF.
            %
            % See also Coefficients.

            coeff2 = coeff;
            coeff2.C = conj(coeff.C);
        end
    end
    methods (Static)
        function i = index(p,l,m)
            % INDEX Index corresponding to p, l, m
            %
            % I = INDEX(P,L,M) calculates the index corresponding to P, L, M. 
            %   P,L,M are the indices of multipole fields. 
            %
            % See also Coefficients, TMatrix, Multipole.
            
            Check.isinteger('p must be equal to 1 or 2',p,'>=',1,'<=',2)
            Check.isinteger('l must be a non-nengative integer',l,'>=',0)
            Check.isinteger('m must be an integer between -l and +l',m,'>=',-l,'<=',l)
            
            i = 2*(l.^2+l+m+1)+(p-2);
        end
        function [p,l,m] = plm(i)
            % PLM Indices p, l, m corresponding to i
            %
            % [P,L,M] = PLM(I,Ip) calculates multipole indices
            %   P, L, M corresponding to index I.
            %
            % See also Coefficients, TMatrix, Multipole.
            
            Check.isinteger('i must be a positive integer',i,'>',0)
            
            % calculates p
            p = mod(i,2);
            p(p==0) = 2;
            % if mod(i,2)~=0
            %     p = 1;
            % else
            %     p = 2;
            % end
            
            % calculates l
            i = i - (p-2); % i = 2*(l^2+l+m+1)
            i = i/2; % i = l^2+l+m+1
            i = i - 1; % i = l^2+l+m
            l = floor(sqrt(i));
            
            % calculates m
            m = i-(l.^2+l);
        end
    end
end