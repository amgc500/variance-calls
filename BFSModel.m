classdef BFSModel < handle & PriceModel
    %BFSMODEL Implementation of a model specified by Beiglboeck,
    %           Friz and Sturm. NB: Vols are divided by 10 for practical
    %           reasons.
    %   Detailed explanation goes here    
    
    properties (SetAccess = immutable)
        s % initial price
        r % interest rate
        T % Time to maturity
        scale % scaling factor for volatility relative to BFS paper
    end
            
    properties (Access = private)
        computed_input
        computed_output
        computed_type
        computed_num
    end
    
    methods (Access = public)
        function obj = BFSModel(s,r,T,scale)
            if nargin == 0
                obj.s = 1;
                obj.r = 0; 
                obj.T = 3;
                obj.scale = 0.1;
            elseif nargin == 3
                obj.s = s;
                obj.r = r; 
                obj.T = T;
                obj.scale = 0.1;
            else
                obj.s = s;
                obj.r = r; 
                obj.T = T;
                obj.scale = scale;
            end
            obj.computed_input = {};
            obj.computed_output = {};
            obj.computed_type = {};
            obj.computed_num=0;
        end
        
        function P = Call(obj,K)
            if isvector(K)
                % P = zeros(size(K));
                z = find(cellfun(@(q) (isequal(q,K)),obj.computed_input)==1,1);
                if ~isempty(z)
                    display('Already computed!')
                    if strcmp(obj.computed_type{z},'Call')
                        P = obj.computed_output{z};
                        return;
                    elseif strcmp(obj.computed_type{z},'Put')
                        P = obj.computed_output{z} + obj.s - exp(-obj.r * obj.T)*K;
                        return;
                    else
                        error('BFSM:MErr1','Unknown type');
                    end
                end   
                
                if obj.T<=1
                    P = BlackScholesModel(obj.s,obj.scale*sqrt(2),obj.r,obj.T).Call(K);
                elseif obj.T <= 2
                    P = (BlackScholesModel(obj.s,obj.scale*sqrt((2+3*(obj.T-1))/obj.T),obj.r,obj.T).Call(K) + BlackScholesModel(obj.s,obj.scale*sqrt((2+1*(obj.T-1))/obj.T),obj.r,obj.T).Call(K))/2;
                else
                    P = (BlackScholesModel(obj.s,obj.scale*sqrt((5+1*(obj.T-2))/obj.T),obj.r,obj.T).Call(K) + BlackScholesModel(obj.s,obj.scale*sqrt((3+3*(obj.T-2))/obj.T),obj.r,obj.T).Call(K))/2;
                end
                obj.computed_num = obj.computed_num+1;
                obj.computed_input{obj.computed_num} = K;
                obj.computed_output{obj.computed_num} = P;
                obj.computed_type{obj.computed_num} = 'Call';
                
            else
                warning('BFSM:CLS1','Strikes must be a vector');
                P = 0;
            end
            
        end
        
        function P = Put(obj,K)
            if isvector(K)
                % P = zeros(size(K));
                z = find(cellfun(@(q) (isequal(q,K)),obj.computed_input)==1,1);
                if ~isempty(z)
                    if strcmp(obj.computed_type{z},'Put')
                        P = obj.computed_output{z};
                        return;
                    elseif strcmp(obj.computed_type{z},'Call')
                        P = obj.computed_output{z} - obj.s + exp(-obj.r * obj.T)*K;
                        return;
                    else
                        error('BFSM:Err2','Unknown type');
                    end
                end
                if obj.T<=1
                    P = BlackScholesModel(obj.s,obj.scale*sqrt(2),obj.r,obj.T).Put(K);
                elseif obj.T <= 2
                    P = (BlackScholesModel(obj.s,obj.scale*sqrt((2+3*(obj.T-1))/obj.T),obj.r,obj.T).Put(K) + BlackScholesModel(obj.s,obj.scale*sqrt((2+1*(obj.T-1))/obj.T),obj.r,obj.T).Put(K))/2;
                else
                    P = (BlackScholesModel(obj.s,obj.scale*sqrt((5+1*(obj.T-2))/obj.T),obj.r,obj.T).Put(K) + BlackScholesModel(obj.s,obj.scale*sqrt((3+3*(obj.T-2))/obj.T),obj.r,obj.T).Put(K))/2;
                end

                obj.computed_num = obj.computed_num+1;
                obj.computed_input{obj.computed_num} = K;
                obj.computed_output{obj.computed_num} = P;
                obj.computed_type{obj.computed_num} = 'Put';
            else
                warning('BFSM:PTS1','Strikes must be a vector');
                P = 0;
            end
        end

        function [p,iv] = simulate(obj,t)
            n = length(t);
            
            omega = 1*(randn(1,1)<0);
            % p = obj.s + zeros(size(t));
            
            iv = obj.intVarOmega(t,omega);
            
            rn = randn(1,(n-1));
            p = obj.s*exp(cumsum([0 (sqrt(diff(iv)).*rn)]) + obj.r*t-iv/2);
            
        end
            
        function iv = intVarOmega(obj,t,omega)
            ind1 = 1*(t<=1);
            ind2 = 1*(t<=2)-ind1;
            ind3 = ones(size(t)) - ind2-ind1;
            if omega==1
                iv = obj.scale^2*(2*t.*ind1 + (3*(t-1)+2).*ind2 + ((t-2)+5).*ind3);
            else
                iv = obj.scale^2*(2*t.*ind1 + (t+1).*ind2 + (3*(t-2)+3).*ind3);
            end
        end
        
        function [p] = VaroptionPrice(obj,F)
            % Price a function of integrated variance.
            
            p = (F(obj.intVarOmega(obj.T,1)) + F(obj.intVarOmega(obj.T,0)))/2;
        end

        function S = saveobj(obj)
            obj.computed_input = {};
            obj.computed_output = {};
            obj.computed_type = {};
            obj.computed_num=0;
            S = obj;
        end            
        
    end

    
    
    
end

