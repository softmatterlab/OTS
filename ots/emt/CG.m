classdef CG < handle
    % CG < handle : Clebsch-Gordan coefficients
    %   CG precalcualtes the Clebsch-Gordan coefficients.
    %   CG can also store the precalcualted values and retrieve them from a file.
    %
    % CG properties:
    %   FileName        - file name
    %   FilePath        - file path
    %   J               - maximum value for computing factorial
    %   JHp1            - maximum J/2+1 [CG.J/2+1]
    %   J2              - maximum J^2 [CG.J^2]
    %   d               - digits
    %   Jmax_created	- maximum J created
    %   Jmax_completed  - maximum J completed 
    %   fact            - precalculated symbolic log(k!)
    %   CM              - sparse matrix with Clebsch-Gordan coefficients
    %
    % CG methods:
    %   CG          -  constructor
    %   indices     -  indices i1,12 corresponding to j1,j2,j,m1,m2,m
    %   jms         -  j1,j2,j,m1,m2,m corresponding to indices i1,12
    %   exact       -  calculate Clebsch-Gordan coefficients using exact formula
    %   create      -  create Clebsch-Gordan coefficients
    %   complete	-  complete set of Clebsch-Gordan coefficients
    %   Cpos        -  Clebsch-Gordan coefficients for positive value of m1,m2 and m
    %   C           -  Clebsch-Gordan coefficients
    %   save        -  save file
    %   wigner3j	-  Wigner 3j symbol
    %
    % See also LGF.

    %   Author: Giovanni Volpe
    %   Revision: 1.0.0
    %   Date: 2015/01/01
    
    properties
        FileName
        FilePath
        J       % maximum J
        JHp1    % maximum J/2+1 [CG.J/2+1]
        J2      % maximum J^2 [CG.J^2]
        d       % digits
        Jmax_created    % maximum J created
        Jmax_completed  % maximum J completed
        fact    % precalculated symbolic log(k!)
        CM      % sparse matrix with Clebsch-Gordan coefficients
    end
    methods
        function obj = CG(varargin)
            % CG() loads the Clebsch-Gordan coefficients from a file
            %   'cg.mat' stored in the same directory as the class CG.
            %   This command is equivalent to
            %   CG('regenerate',false, ...
            %       'filename','cg', ...
            %       'filepath',which('CG'), ...
            %       'J',1000, ...
            %       'digits',100)
            %
            % CG('regenerate',true,'J',Jmax) regenerates the Clebsh-Gordan
            %   coefficients up to Jmax.
            %
            % CG('filename',FileName) loads the Clebsh-Gordan coefficients
            %   from FileName (in the directory of the class CG).
            %
            % CG('pathname',PathName) loads the Clebsh-Gordan coefficients
            %   from the file cg.mat in the directory PathName.
            %
            % CG('PropertyName',PropertyValue) sets the property PropertyName to PropertyValue.
            %   The properties listed below can be used:
            %       regenerate  -  regenarate [defalut regenerate=false]
            %       filename    -  data file name [defalut filename='cg']
            %       filepath    -  data path name [defalut pathname=which('CG')]
            %       J           -  maximum of J [defalut J=1000]
            %       d           -  digits [defalut d=100]
            %   Note that the last two properties are only effective is
            %   regenerate=false.
            %
            % See also CG.

            % data file name
            regenerate = false;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'regenerate')
                    regenerate = varargin{n+1};
                end
            end

            % data file name
            obj.FileName = 'cg';
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'filename')
                    obj.FileName = varargin{n+1};
                end
            end
            
            % data path name
            obj.FilePath = which('CG');
            obj.FilePath = obj.FilePath(1:end-4);
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'filepath')
                    obj.FilePath = varargin{n+1};
                end
            end
            
            if regenerate
                
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
                for n = 0:1:J
                    obj.fact{n+1} = vpa(subs(log_k_fac,'k',n));
                end
                digits(dtmp);

                obj.CM = sparse(J^3,J^2);
                obj.Jmax_created = -1;
                obj.Jmax_completed = -1;            
            else
                
                load([obj.FilePath filesep() obj.FileName '.mat']) % CM J Jmax_created Jmax_completed d fact
                
                obj.J = J;

                obj.d = d;
                obj.fact = fact;

                obj.CM = CM;
                obj.Jmax_created = Jmax_created;
                obj.Jmax_completed = Jmax_completed;
            end
            
            obj.JHp1 = J/2+1;
            obj.J2 = J^2;
        end
        function [i1,i2] = indices(cg,j1,j2,j,m1,m2,m)
            % INDICES Scaler indices
            %
            % [I1,I2] = INDICES(Cg,J1,J2,J,M1,M2,M) determine the indices
            %   corresponding to J1, J2, J, M1, M2, M. 
            %
            %  See also CG.
            
            i1 = (j1+1) ...
                + (j2+1)*cg.J ...
                + (j+1)*cg.J2;
            i2 = (m1+cg.JHp1) ...
                + (m2+cg.JHp1)*cg.J;
        end
        function [j1,j2,j,m1,m2,m] = jms(cg,i1,i2)
            % JMS Indices J1, J2, J, M1, M2, M
            %
            % [J1,J2,J,M1,M2,M] = JMS(CG,I1,I2) determine the complete set
            %   of indices J1, J2, J, M1, M2, M corresponding to indices I1 and I2. 
            %
            %  See also CG.
            
            j1 = rem(i1,cg.J)-1;
            
            i1 = (i1-(j1+1))/cg.J;
            j2 = rem(i1,cg.J)-1;
            
            i1 = (i1-(j2+1))/cg.J;
            j = rem(i1,cg.J)-1;
            
            m1 = rem(i2,cg.J)-cg.JHp1;
            
            i2 = (i2 - (m1+cg.JHp1))/cg.J;
            m2 = rem(i2,cg.J)-cg.JHp1;
            
            m = m1+m2;
        end
        function c = exact(cg,j1,j2,j,m1,m2,m)
            % EXACT Calculate Clebsch-Gordan coefficients using exact formula
            %
            % EXACT(CG,J1,J2,J,M1,M2,M) calculate the set of Clebsch-Gordan coefficients 
            %   for J1,J2,J,M1,M2,M using the exact formula.
            %
            %  See also CG.

            if (abs(m1)>j1)||(abs(m2)>j2)||(abs(m)>j)||(m1+m2~=m)||(j<abs(j1-j2))||(j>(j1+j2))

                c = 0;
                
            else

                dtmp = digits;
                digits(cg.d);
                
                % calculate denominator
                kmin = max([0,j2-j-m1,j1-j+m2]);
                kmax = min([j1+j2-j,j1-m1,j2+m2]);
                
                Sk = 0;
                for k = kmin:1:kmax
                    Sk = Sk + (-1)^(k) ...
                        * exp( ...
                        -cg.fact{k + 1} ...
                        -cg.fact{j1+j2-j-k + 1} ...
                        -cg.fact{j1-m1-k + 1}...
                        -cg.fact{j2+m2-k + 1} ...
                        -cg.fact{j-j2+m1+k + 1} ...
                        -cg.fact{j-j1-m2+k + 1} ...
                        );
                end                

                % Calculate Clebsch-Gordan coefficient
                c = sqrt(2*j+1) ...
                    * exp( 0.5 * ...
                    ( ...
                    cg.fact{j1+j2-j + 1} ...
                    +cg.fact{j+j1-j2 + 1} ...
                    +cg.fact{j+j2-j1 + 1} ...
                    -cg.fact{j1+j2+j+1 + 1} ...
                    +cg.fact{j1+m1 + 1} ...
                    +cg.fact{j1-m1 + 1} ...
                    +cg.fact{j2-m2 + 1} ...
                    +cg.fact{j2+m2 + 1} ...
                    +cg.fact{j+m + 1} ...
                    +cg.fact{j-m + 1} ...
                    ) ...
                    ) * Sk;

                c = eval(c);

                if abs(c)<1e-15
                    c = 0;
                end
                
                digits(dtmp);
            end
        end
        function create(cg,Jmax,p)
            % CREATE Create set of Clebsch-Gordan coefficients 
            %
            % COMPLETE(CG,JMAX,P) calculate set of Clebsch-Gordan coefficients 
            %   untill JMAX (maximum number of J) and display them with
            %   probability P (defalt value is 0.01). P = 1 will show all of them..
            %
            %  See also CG.
            
            if nargin<3
                p = 0.01;
            end
            
            tic
            for j1 = 0:1:Jmax
                for j2 = 0:1:Jmax
                    for j = abs(j2-j1):1:min(j1+j2,Jmax)  % j = j2:1:j1+j2
                        for m1 = 0:1:j1 % m1 = -j1:1:j1  
                            for m2 = m1:1:j2 % m2 = -j2:1:j2
                                m = m1+m2;
                                
                                [i1,i2] = cg.indices(j1,j2,j,m1,m2,m);
                                if cg.CM(i1,i2)==0
                                    cg.CM(i1,i2) = cg.exact(j1,j2,j,m1,m2,m);
                                
                                    if rand>1-p
                                        disp([num2str(cg.CM(i1,i2)) ' for j1,j2,m1,m2,j,m = ' int2str(j1) ' ' int2str(j2) ' ' int2str(m1) ' ' int2str(m2) ' ' int2str(j)  ' ' int2str(m) ' [' num2str(toc) 's]'])
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            cg.Jmax_created = max(Jmax,cg.Jmax_created);
        end
        function complete(cg,Jmax)
            % COMPLETE Complete set of Clebsch-Gordan coefficients 
            %
            % COMPLETE(CG,JMAX) complete set of Clebsch-Gordan coefficients 
            %   until JMAX (maximum number of J). Must be used after
            %   create.
            %
            %  See also CG, CG.create.
        
            Jmax = min(Jmax,cg.Jmax_created);
        
            % C(j2,j1,j;m2,m1,m) - permutation 213
            for j1 = 0:1:Jmax
                for j2 = 0:1:Jmax
                    for j = abs(j2-j1):1:min(j1+j2,Jmax)
        
                        s = (-1)^(j-j1-j2);
        
                        for m1 = 0:1:j1
                            for m2 = m1+1:1:j2
                                m = m1+m2;
                                
                                if (abs(m1)<=j1)&&(abs(m2)<=j2)&&(abs(m)<=j)&&(m1+m2==m)&&(j>=abs(j1-j2))&&(j<=(j1+j2))
                                
                                    [i1,i2] = cg.indices(j1,j2,j,m1,m2,m);

                                    [k1,k2] = cg.indices(j2,j1,j,m2,m1,m);
                                    cg.CM(k1,k2) = s * cg.CM(i1,i2);
                                    % disp(['GENERATED ' int2str(j2) ' ' int2str(j1) ' ' int2str(j) ' ' int2str(m2) ' ' int2str(m1) ' ' int2str(m)])
                                else 
                                    % disp(['NOT GENERATED ' int2str(j2) ' ' int2str(j1) ' ' int2str(j) ' ' int2str(m2) ' ' int2str(m1) ' ' int2str(m)])
                                end
                            end
                        end
                        
                        disp(['1st phase: j1=' int2str(j1) ' j2=' int2str(j2) ' j=' int2str(j) ' completed'])
                    end
                end
            end
        
            for j1 = 0:1:Jmax
                for j2 = 0:1:Jmax
                    for j = abs(j2-j1):1:min(j1+j2,Jmax)
                        
                        s = (-1)^(j-j1-j2);
                        s1 = (-1)^(-j1) * sqrt((2*j2+1)/(2*j+1));
                        s2 = (-1)^(-j2) * sqrt((2*j1+1)/(2*j+1));
        
                        for m1 = 0:1:j1
                            for m2 = 0:1:j2
                                m = m1+m2;
        
                                [t1,t2] = cg.indices( ...
                                    [j1 j1  j1  j   j2  j], ...
                                    [j2 j2  j   j1  j   j2], ...
                                    [j  j   j2  j2  j1  j1], ...
                                    [m1 -m1 m1  m   -m2 -m], ...
                                    [m2 -m2 -m  -m1 m   m2], ...
                                    [m  -m  -m2 m2  m1  -m1] ...
                                    );

                                % [i1,i2] = cg.indices(j1,j2,j,m1,m2,m);
                                c = cg.CM(t1(1),t2(1));

                                % C(j1,j2,j;-m1,-m2,-m) - change of sign
                                if (abs(m1)<=j1)&&(abs(m2)<=j2)&&(abs(m)<=j)&&(m1+m2==m)&&(j>=abs(j1-j2))&&(j<=(j1+j2))
                                    % [k1,k2] = cg.indices(j1,j2,j,-m1,-m2,-m);
                                    cg.CM(t1(2),t2(2)) = s * c;
                                    % disp(['GENERATED ' int2str(j1) ' ' int2str(j2) ' ' int2str(j) ' ' int2str(-m1) ' ' int2str(-m2) ' ' int2str(-m)])
                                else 
                                    % disp(['NOT GENERATED ' int2str(j1) ' ' int2str(j2) ' ' int2str(j) ' ' int2str(-m1) ' ' int2str(-m2) ' ' int2str(-m)])
                                end
                                    

                                % C(j1,j,j2;m1,-m,-m2) - permutation 132
                                if (abs(m1)<=j1)&&(abs(m)<=j)&&(abs(m2)<=j2)&&(m1-m==-m2)&&(j2>=abs(j-j1))&&(j2<=(j+j1))
                                    % [k1,k2] = cg.indices(j1,j,j2,m1,-m,-m2);
                                    cg.CM(t1(3),t2(3)) = (-1)^m1 * s1* c;
                                    % disp(['GENERATED ' int2str(j1) ' ' int2str(j) ' ' int2str(j2) ' ' int2str(m1) ' ' int2str(-m) ' ' int2str(-m2)])
                                else 
                                    % disp(['NOT GENERATED ' int2str(j1) ' ' int2str(j) ' ' int2str(j2) ' ' int2str(m1) ' ' int2str(-m) ' ' int2str(-m2)])
                                end

                                % C(j,j1,j2;m,-m1,m2) - permutation 312
                                if (abs(m)<=j)&&(abs(m1)<=j1)&&(abs(m2)<=j2)&&(m-m1==m2)&&(j2>=abs(j-j1))&&(j2<=(j+j1))
                                    % [k1,k2] = cg.indices(j,j1,j2,m,-m1,m2);
                                    cg.CM(t1(4),t2(4)) = (-1)^m1 * s1 * c;
                                    % disp(['GENERATED ' int2str(j) ' ' int2str(j1) ' ' int2str(j2) ' ' int2str(m) ' ' int2str(-m1) ' ' int2str(m2)])
                                else 
                                    % disp(['NOT GENERATED ' int2str(j) ' ' int2str(j1) ' ' int2str(j2) ' ' int2str(m) ' ' int2str(-m1) ' ' int2str(m2)])
                                end

                                % C(j2,j,j1;-m2,m,m1) - permutation 231
                                if (abs(m2)<=j2)&&(abs(m)<=j)&&(abs(m1)<=j1)&&(-m2+m==m1)&&(j1>=abs(j2-j))&&(j1<=(j2+j))
                                    % [k1,k2] = cg.indices(j2,j,j1,-m2,m,m1);
                                    cg.CM(t1(5),t2(5)) = (-1)^(-m2) * s2 * c;
                                    % disp(['GENERATED ' int2str(j2) ' ' int2str(j) ' ' int2str(j1) ' ' int2str(-m2) ' ' int2str(m) ' ' int2str(m1)])
                                else 
                                    % disp(['NOT GENERATED ' int2str(j2) ' ' int2str(j) ' ' int2str(j1) ' ' int2str(-m2) ' ' int2str(m) ' ' int2str(m1)])
                                end

                                % C(j,j2,j1;-m,m2,-m1) - permutation 321
                                if (abs(m)<=j)&&(abs(m2)<=j2)&&(abs(m1)<=j1)&&(-m+m2==-m1)&&(j1>=abs(j2-j))&&(j1<=(j2+j))
                                    % [k1,k2] = cg.indices(j,j2,j1,-m,m2,-m1);
                                    cg.CM(t1(6),t2(6)) = (-1)^(-m2) * s2 * c;
                                    % disp(['GENERATED ' int2str(j) ' ' int2str(j2) ' ' int2str(j1) ' ' int2str(-m) ' ' int2str(m2) ' ' int2str(-m1)])
                                else 
                                    % disp(['NOT GENERATED ' int2str(j) ' ' int2str(j2) ' ' int2str(j1) ' ' int2str(-m) ' ' int2str(m2) ' ' int2str(-m1)])
                                end
                            end
                        end
                        
                        disp(['2nd phase: j1=' int2str(j1) ' j2=' int2str(j2) ' j=' int2str(j) ' completed'])
                    end
                end
            end
        
            cg.Jmax_completed = max(Jmax,cg.Jmax_completed);
            
            % function complete(cg,Jmax,p)
            % 
            %     if nargin<3
            %         p = 0.001;
            %     end
            % 
            %     Jmax = min(Jmax,cg.Jmax_created);
            % 
            %     % C(j2,j1,j;m2,m1,m) - permutation 213
            %     for j1 = 0:1:Jmax
            %         for j2 = 0:1:Jmax
            %             for j = abs(j2-j1):1:min(j1+j2,Jmax)
            %                 for m1 = 0:1:j1
            %                     for m2 = m1+1:1:j2
            %                         m = m1+m2;
            % 
            %                         [i1,i2] = cg.indices(j1,j2,j,m1,m2,m);
            % 
            %                         [k1,k2] = cg.indices(j2,j1,j,m2,m1,m);
            %                         cg.CM(k1,k2) = (-1)^(j-j1-j2) * cg.CM(i1,i2);
            % 
            %                         if rand>1-p
            %                             disp([num2str(cg.CM(k1,k2)) ' for j1,j2,m1,m2,j,m = ' int2str(j2) ' ' int2str(j1) ' ' int2str(m2) ' ' int2str(m1) ' ' int2str(j)  ' ' int2str(m)])
            %                         end
            %                     end
            %                 end
            %             end
            %         end
            %     end
            % 
            %     % C(j1,j2,j;-m1,-m2,-m) - change of sign
            %     for j1 = 0:1:Jmax
            %         for j2 = 0:1:Jmax
            %             for j = abs(j2-j1):1:min(j1+j2,Jmax)
            %                 for m1 = 0:1:j1
            %                     for m2 = 0:1:j2
            %                         m = m1+m2;
            % 
            %                         [i1,i2] = cg.indices(j1,j2,j,m1,m2,m);
            % 
            %                         [k1,k2] = cg.indices(j1,j2,j,-m1,-m2,-m);
            %                         cg.CM(k1,k2) = (-1)^(j-j1-j2) * cg.CM(i1,i2);
            % 
            %                         if rand>1-p
            %                             disp([num2str(cg.CM(k1,k2)) ' for j1,j2,m1,m2,j,m = ' int2str(j1) ' ' int2str(j2) ' ' int2str(-m1) ' ' int2str(-m2) ' ' int2str(j)  ' ' int2str(-m)])
            %                         end
            %                     end
            %                 end
            %             end
            %         end
            %     end
            % 
            %     % C(j,j2,j1;-m,m2,-m1) - permutation 321
            %     for j1 = 0:1:Jmax
            %         for j2 = 0:1:Jmax
            %             for j = abs(j2-j1):1:min(j1+j2,Jmax)
            %                 for m1 = 0:1:j1
            %                     for m2 = 0:1:j2
            %                         m = m1+m2;
            % 
            %                         [i1,i2] = cg.indices(j1,j2,j,m1,m2,m);
            % 
            %                         [k1,k2] = cg.indices(j,j2,j1,-m,m2,-m1);
            %                         cg.CM(k1,k2) = (-1)^(-m2-j2) * sqrt((2*j1+1)/(2*j+1)) * cg.CM(i1,i2);
            % 
            %                         if rand>1-p
            %                             disp([num2str(cg.CM(k1,k2)) ' for j1,j2,m1,m2,j,m = ' int2str(j) ' ' int2str(j2) ' ' int2str(-m) ' ' int2str(m2) ' ' int2str(j1)  ' ' int2str(-m1)])
            %                         end
            %                     end
            %                 end
            %             end
            %         end
            %     end
            % 
            %     % C(j1,j,j2;m1,-m,-m2) - permutation 132
            %     for j1 = 0:1:Jmax
            %         for j2 = 0:1:Jmax
            %             for j = abs(j2-j1):1:min(j1+j2,Jmax)
            %                 for m1 = 0:1:j1
            %                     for m2 = 0:1:j2
            %                         m = m1+m2;
            % 
            %                         [i1,i2] = cg.indices(j1,j2,j,m1,m2,m);
            % 
            %                         [k1,k2] = cg.indices(j1,j,j2,m1,-m,-m2);
            %                         cg.CM(k1,k2) = (-1)^(m1-j1) * sqrt((2*j2+1)/(2*j+1)) * cg.CM(i1,i2);
            % 
            %                         if rand>1-p
            %                             disp([num2str(cg.CM(k1,k2)) ' for j1,j2,m1,m2,j,m = ' int2str(j1) ' ' int2str(j) ' ' int2str(m1) ' ' int2str(-m) ' ' int2str(j2)  ' ' int2str(-m2)])
            %                         end
            %                     end
            %                 end
            %             end
            %         end
            %     end
            % 
            %     % C(j,j1,j2;m,-m1,m2) - permutation 312
            %     for j1 = 0:1:Jmax
            %         for j2 = 0:1:Jmax
            %             for j = abs(j2-j1):1:min(j1+j2,Jmax)
            %                 for m1 = 0:1:j1
            %                     for m2 = 0:1:j2
            %                         m = m1+m2;
            % 
            %                         [i1,i2] = cg.indices(j1,j2,j,m1,m2,m);
            % 
            %                         [k1,k2] = cg.indices(j,j1,j2,m,-m1,m2);
            %                         cg.CM(k1,k2) = (-1)^(m1-j1) * sqrt((2*j2+1)/(2*j+1)) * cg.CM(i1,i2);
            % 
            %                         if rand>1-p
            %                             disp([num2str(cg.CM(k1,k2)) ' for j1,j2,m1,m2,j,m = ' int2str(j) ' ' int2str(j1) ' ' int2str(m) ' ' int2str(-m1) ' ' int2str(j2)  ' ' int2str(m2)])
            %                         end
            %                     end
            %                 end
            %             end
            %         end
            %     end
            % 
            %     % C(j2,j,j1;-m2,m,m1) - permutation 231
            %     for j1 = 0:1:Jmax
            %         for j2 = 0:1:Jmax
            %             for j = abs(j2-j1):1:min(j1+j2,Jmax)
            %                 for m1 = 0:1:j1
            %                     for m2 = 0:1:j2
            %                         m = m1+m2;
            % 
            %                         [i1,i2] = cg.indices(j1,j2,j,m1,m2,m);
            % 
            %                         [k1,k2] = cg.indices(j2,j,j1,-m2,m,m1);
            %                         cg.CM(k1,k2) = (-1)^(-m2-j2) * sqrt((2*j1+1)/(2*j+1)) * cg.CM(i1,i2);
            % 
            %                         if rand>1-p
            %                             disp([num2str(cg.CM(k1,k2)) ' for j1,j2,m1,m2,j,m = ' int2str(j2) ' ' int2str(j) ' ' int2str(-m2) ' ' int2str(m) ' ' int2str(j1)  ' ' int2str(m1)])
            %                         end
            %                     end
            %                 end
            %             end
            %         end
            %     end
            % 
            %     cg.Jmax_completed = max(Jmax,cg.Jmax_completed);
            % end
            
            % function complete(cg,Jmax,p)
            % 
            %     if nargin<3
            %         p = 0.001;
            %     end
            % 
            %     Jmax = min(Jmax,cg.Jmax_created);
            % 
            %     % C(j2,j1,j;m2,m1,m) - permutation 213
            %     [I1,I2] = find(cg.CM~=0);
            %     [J1,J2,J,M1,M2,M] = cg.jms(I1,I2);
            %     Ci = diag(cg.CM(I1,I2));
            % 
            %     [K1,K2] = cg.indices(J2,J1,J,M2,M1,M);
            % 
            %     Ck = (-1).^(J-J1-J2) .* Ci;
            % 
            %     for k = 1:1:length(K1)
            %         cg.CM(K1(k),K2(k)) = Ck(k);
            %     end
            % 
            %     % other premutations
            %     [I1,I2] = find(cg.CM~=0);
            %     [J1,J2,J,M1,M2,M] = cg.jms(I1,I2);
            %     Ci = diag(cg.CM(I1,I2));
            % 
            %     % C(j1,j2,j;-m1,-m2,-m) - change of sign
            %     [K1,K2] = cg.indices(J1,J2,J,-M1,-M2,-M);
            % 
            %     Ck = (-1).^(J-J1-J2) .* Ci;
            % 
            %     for k = 1:1:length(K1)
            %         cg.CM(K1(k),K2(k)) = Ck(k);
            %     end
            % 
            %     % C(j,j2,j1;-m,m2,-m1) - permutation 321
            %     [K1,K2] = cg.indices(J,J2,J1,-M,M2,-M1);
            % 
            %     Ck = (-1).^(-M2-J2) .* sqrt((2*J1+1)./(2*J+1)) .* Ci;
            % 
            %     for k = 1:1:length(K1)
            %         cg.CM(K1(k),K2(k)) = Ck(k);
            %     end
            % 
            %     % C(j1,j,j2;m1,-m,-m2) - permutation 132
            %     [K1,K2] = cg.indices(J1,J,J2,M1,-M,-M2);
            % 
            %     Ck = (-1).^(M1-J1) .* sqrt((2*J2+1)./(2*J+1)) .* Ci;
            % 
            %     for k = 1:1:length(K1)
            %         cg.CM(K1(k),K2(k)) = Ck(k);
            %     end
            % 
            %     % C(j,j1,j2;m,-m1,m2) - permutation 312
            %     [K1,K2] = cg.indices(J,J1,J2,M,-M1,M2);
            % 
            %     Ck = (-1).^(M1-J1) .* sqrt((2*J2+1)./(2*J+1)) .* Ci;
            % 
            %     for k = 1:1:length(K1)
            %         cg.CM(K1(k),K2(k)) = Ck(k);
            %     end
            % 
            %     % C(j2,j,j1;-m2,m,m1) - permutation 231
            %     [K1,K2] = cg.indices(J2,J,J1,-M2,M,M1);
            % 
            %     Ck = (-1).^(-M2-J2) .* sqrt((2*J1+1)./(2*J+1)) .* Ci;
            % 
            %     for k = 1:1:length(K1)
            %         cg.CM(K1(k),K2(k)) = Ck(k);
            %     end
            % 
            % 
            %     cg.Jmax_completed = max(Jmax,cg.Jmax_completed);
            % end            
        end
        function c = Cpos(cg,j1,j2,j,m1,m2,m)
            % CPOS Clebsch-Gordan coefficients for positive value of m1 and m2
            %
            % C = CPOS(CG,J1,J2,J,M1,M2,M) determine the Clebsch-Gordan coefficients 
            %   using precalculation CG for J1, J2, J and positive values
            %   of M1, M2, M multipole indices.
            %
            %  See also CG.
            
            c = zeros(size(j1));

            j1 = reshape(j1,1,numel(j1));
            j2 = reshape(j2,1,numel(j1));
            j = reshape(j,1,numel(j1));
            m1 = reshape(m1,1,numel(j1));
            m2 = reshape(m2,1,numel(j1));
            m = reshape(m,1,numel(j1));
            
            % Pre-calculated
            ind = find(m1>=0 & m2>=0 & m2>=m1);
            if ~isempty(ind)
                [i1,i2] = cg.indices(j1(ind),j2(ind),j(ind),m1(ind),m2(ind),m(ind));
                c(ind) = diag(cg.CM(i1,i2));
            end
            
            % C(j2,j1,j;m2,m1,m) - permutation 213
            ind = find(m1>=0 & m2>=0 & m2<m1);
            if ~isempty(ind)
                [i1,i2] = cg.indices(j2(ind),j1(ind),j(ind),m2(ind),m1(ind),m(ind));
                c(ind) = (-1).^(j1(ind)+j2(ind)-j(ind)) .* diag(cg.CM(i1,i2))';
            end            
        end
        function c = C(cg,j1,j2,j,m1,m2,m)
            % C Clebsch-Gordan coefficients
            %
            % C = C(CG,J1,J2,J,M1,M2,M) return the precalcualted Clebsch-Gordan coefficients 
            %   for J1, J2, J, M1, M2, M.
            %
            %  See also CG.
            
            % [i1,i2] = cg.indices(j1,j2,j,m1,m2,m);
            % c1 = diag(cg.CM(i1,i2));
            % c1 = full(reshape(c1,size(j1)));
            % c1(m1+m2~=m) = 0;
            
            c = zeros(size(j1));

            j1 = reshape(j1,1,numel(j1));
            j2 = reshape(j2,1,numel(j1));
            j = reshape(j,1,numel(j1));
            m1 = reshape(m1,1,numel(j1));
            m2 = reshape(m2,1,numel(j1));
            m = reshape(m,1,numel(j1));
            
            % m1,m2>=0
            ind = find(m1>=0 & m2>=0 & m1+m2==m);
            if ~isempty(ind)
                c(ind) = cg.Cpos(j1(ind),j2(ind),j(ind),m1(ind),m2(ind),m(ind));
            end
                        
            % C(j1,j2,j;-m1,-m2,-m) - change of sign
            ind = find(m1<=0 & m2<=0 & m1+m2==m);
            if ~isempty(ind)
                c(ind) = (-1).^(j1(ind)+j2(ind)-j(ind)) .* cg.Cpos(j1(ind),j2(ind),j(ind),-m1(ind),-m2(ind),-m(ind));
            end
                        
            % C(j,j2,j1;-m,m2,-m1) - permutation 321
            ind = find(m1<0 & m2>=0 & m<=0 & m1+m2==m);
            if ~isempty(ind)
                c(ind) = (-1).^(j2(ind)+m2(ind)) .* sqrt((2*j(ind)+1)./(2*j1(ind)+1)) .* cg.Cpos(j(ind),j2(ind),j1(ind),-m(ind),m2(ind),-m1(ind));
            end
            
            % C(j1,j,j2;m1,-m,-m2) - permutation 132
            ind = find(m1>=0 & m2<0 & m<=0 & m1+m2==m);
            if ~isempty(ind)
                c(ind) = (-1).^(j1(ind)-m1(ind)) .* sqrt((2*j(ind)+1)./(2*j2(ind)+1)) .* cg.Cpos(j1(ind),j(ind),j2(ind),m1(ind),-m(ind),-m2(ind));
            end
                        
            % C(j,j1,j2;m,-m1,m2) - permutation 312
            ind = find(m1<0 & m2>=0 & m>=0 & m1+m2==m);
            if ~isempty(ind)
                c(ind) = (-1).^(j1(ind)-m1(ind)) .* sqrt((2*j(ind)+1)./(2*j2(ind)+1)) .* cg.Cpos(j(ind),j1(ind),j2(ind),m(ind),-m1(ind),m2(ind));
            end
            
            % C(j2,j,j1;-m2,m,m1) - permutation 231
            ind = find(m1>=0 & m2<0 & m>=0 & m1+m2==m);
            if ~isempty(ind)
                c(ind) = (-1).^(j2(ind)+m2(ind)) .* sqrt((2*j(ind)+1)./(2*j1(ind)+1)) .* cg.Cpos(j2(ind),j(ind),j1(ind),-m2(ind),m(ind),m1(ind));
            end
        end
        function save(cg,varargin)
            % SAVE Save file
            %
            % SAVE(CG) save the Clebsch-Gordan coefficients in a file.
            %
            % SAVE('PropertyName',PropertyValue) sets the property PropertyName to PropertyValue. 
            %   The properties listed below can be used:
            %       filename	-   data file name [default filename='cg']
            %       filepath	-   data file path [default filepath=which('CG')]
            %
            % See also CG.
            
            % data file name
            FileName = cg.FileName;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'filename')
                    FileName = varargin{n+1};
                end
            end
            
            % data path name
            FilePath = cg.FilePath;
            for n = 1:2:length(varargin)
                if strcmpi(varargin{n},'filepath')
                    FilePath = varargin{n+1};
                end
            end
            
            J = cg.J;
            d = cg.d;
            Jmax_created = cg.Jmax_created;
            Jmax_completed = cg.Jmax_completed;
            fact = cg.fact;
            CM = cg.CM;          
            
            save([FilePath filesep() FileName '.mat'],'J','d','Jmax_created','Jmax_completed','fact','CM')
        end
        function w = wigner3j(cg,j1,j2,j,m1,m2,mm)
            % WIGNER3J Wigner 3j symbol
            %
            % WIGNER3J(CG,J1,J2,J,M1,M2,M) determine the Wigner 3j symbol
            %   using the precalculated CG coefficients 
            %   for J1, J2, J, M1, M2, M.
            %
            % See also CG.

            w = (-1).^(j1-j2+m1+m2)./sqrt(2*j+1).*cg.C(j1,j2,j,m1,m2,-mm);
        end
    end
end