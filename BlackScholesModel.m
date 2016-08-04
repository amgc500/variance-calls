classdef BlackScholesModel < handle & PriceModel
    %BLACKSCHOLESMODEL Representation of Black Scholes
    %   Detailed explanation goes here
    
    properties (SetAccess = immutable)
        s % initial price
        sigma % initial sigma
        r % interest rate
        T % Time to maturity
        
    end

    
    properties (Access = private)
        computed_input
        computed_output
        computed_type
        computed_num
    end
        
    
    methods (Access = public)
        function obj = BlackScholesModel(s,sigma,r,T)
            if nargin == 0
                obj.s = 1;
                obj.sigma = .04; 
                obj.r = 0; 
                obj.T = 1;
                
            else
                obj.s = s;
                obj.sigma = sigma; 
                obj.r = r; 
                obj.T = T;
            end
            obj.computed_input = {};
            obj.computed_output = {};
            obj.computed_type = {};
            obj.computed_num=0;
        end

        
        function P = Call(obj,K)
            if isvector(K)
                P = zeros(size(K));
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
                        error('BS:Err1','Unknown type');
                    end
                end   
                
                if (obj.T == 0)
                    P = max(obj.s-K,0);
                elseif (obj.sigma == 0)
                    P = max(obj.s-exp(-obj.r*obj.T)*K,0);
                else
                    for j = 1:length(K)
                        if(K(j) == 0)
                            P(j) = obj.s;
                        else
                            d_1 = (log(obj.s/K(j)) + (obj.r+obj.sigma^2/2)*obj.T)/obj.sigma/sqrt(obj.T);
                            d_2 = d_1 - obj.sigma*sqrt(obj.T);
                            P(j) = obj.s * normcdf(d_1) - K(j)*exp(-obj.r*obj.T) * normcdf(d_2);
                        end
                    end
                end
                  
                obj.computed_num = obj.computed_num+1;
                obj.computed_input{obj.computed_num} = K;
                obj.computed_output{obj.computed_num} = P;
                obj.computed_type{obj.computed_num} = 'Call';
                
            else
                warning('BS:CLS1','Strikes must be a vector');
                P = 0;
            end
            
        end
        
        function P = Put(obj,K)
            if isvector(K)
                P = zeros(size(K));
                z = find(cellfun(@(q) (isequal(q,K)),obj.computed_input)==1,1);
                if ~isempty(z)
                    if strcmp(obj.computed_type{z},'Put')
                        P = obj.computed_output{z};
                        return;
                    elseif strcmp(obj.computed_type{z},'Call')
                        P = obj.computed_output{z} - obj.s + exp(-obj.r * obj.T)*K;
                        return;
                    else
                        error('BS:Err2','Unknown type');
                    end
                end
                
                if (obj.T == 0)
                    P = max(K-obj.s,0);
                elseif (obj.sigma == 0)
                    P = max(exp(-obj.r*obj.T)*K-obj.s,0);
                else
                    for j = 1:length(K)
                        if(K(j) == 0)
                            P(j) = 0;
                        else
                            d_1 = (log(obj.s/K(j)) + (obj.r+obj.sigma^2/2)*obj.T)/obj.sigma/sqrt(obj.T);
                            d_2 = d_1 - obj.sigma*sqrt(obj.T);
                            P(j) = -obj.s * normcdf(-d_1) + K(j)*exp(-obj.r*obj.T) * normcdf(-d_2);
                        end
                    end
                end
                
                obj.computed_num = obj.computed_num+1;
                obj.computed_input{obj.computed_num} = K;
                obj.computed_output{obj.computed_num} = P;
                obj.computed_type{obj.computed_num} = 'Put';
            else
                warning('BS:PTS1','Strikes must be a vector');
                P = 0;
            end
        end

        function [p,iv] = simulate(obj,t)
            n = length(t);
            iv = t*obj.sigma^2;
            
            rn = randn(1,(n-1));
            p = obj.s*exp(cumsum([0 obj.sigma.*rn]) + obj.r*t-iv/2);
            
        end
        
        function [p] = VaroptionPrice(obj,F)
            % Price a function of integrated variance.
            
            p = F(obj.T*obj.sigma^2);
        end
        
        function S = saveobj(obj)
            obj.computed_input = {};
            obj.computed_output = {};
            obj.computed_type = {};
            obj.computed_num=0;
            S = obj;
        end            

    end
        
    methods
    end
    
end

