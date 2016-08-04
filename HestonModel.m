%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file contains code to accompany the paper 
% "Optimal robust bounds for variance options" 
% by Alexander M G Cox and Jiajie Wang
%
% Copyright (c) 2013 Alexander Cox and Jiajie Wang
%
% Please send comments, corrections etc. to
% a.m.g.cox@bath.ac.uk
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as
% published by the Free Software Foundation, either version 3 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
%
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef HestonModel < handle & PriceModel
    %HESTONMODEL Representation of Heston
    %   Stores the parameters of a Heston Model, and computes prices as
    %   required.
    %   
    %   To compute the Call price, we use the algorithm as described by
    %   Janek, Kluge, Weron and Wystup (2010) - FX smile in the Heston
    %   model.
    %
    %   We also apply some smoothing to the call prices at high strikes,
    %   where the prices may be relatively very accurate, but terms such as
    %   second derivatives can become very unstable.

    
    
    properties (SetAccess = immutable)
        s % initial price
        V0 % initial vol
        xi % vol of vol
        r % interest rate
        T % Time to maturity
        kappa % speed of mean reversion
        theta % long-run mean variance
        rho % correlation of asset and vol bm        
    end
    
    properties (Access = private)
        computed_input
        computed_output
        computed_type
        computed_num
    end
        
    
    methods (Access = public)
        function obj = HestonModel(s,V0,xi,r,T,kappa,theta,rho)
            if nargin == 0
                obj.s = 1;
                obj.V0 = .047; 
                obj.xi = .39; 
                obj.r = 0; 
                obj.kappa = 1.3;
                obj.theta = 0.035;
                obj.rho = -0.7;
                obj.T = 1;
                
            else
                obj.s = s;
                obj.V0 = V0; 
                obj.xi = xi; 
                obj.r = r; 
                obj.kappa = kappa;
                obj.theta = theta;
                obj.rho = rho;
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
                        error('HPM:Err1','Unknown type');
                    end
                end   
                for i = 1:length(K)
                    I1 = quadgk(@(z) obj.HVI1(z,K(i)),0,inf,'RelTol',1e-10,'AbsTol',0);
                    I2 = quadgk(@(z) obj.HVI2(z,K(i)),0,inf,'RelTol',1e-10,'AbsTol',0);
            
                    % Calculate option price
                    P(i) = obj.s.*exp(-obj.r.*obj.T).*(0.5+1/pi.*I1) - K(i).*exp(-obj.r.*obj.T).*(0.5+1/pi.*I2);
                    P(i) = max(max(P(i),0),obj.s-K(i).*exp(-obj.r.*obj.T));
                end
                
                % Find 'feasible' Call prices that satisfy usual conditions
                % Find nearby solution that is convex, decreasing, above
                % lower bounds
                
                n = length(K);

                if n >=2
                    A = diag([0 1./diff(K)'] + [1./diff(K)' 0],0) -diag(1./diff(K),1) -diag(1./diff(K),-1);
                    
                    
                    A(1,1) = 1;
                    A(1,2) = -1;
                    A(n,n-1) = -1;
                    A(n,n) = 1;
                    
                    b = zeros(n,1);
                    
                    b(1) = (K(2)-K(1))*exp(-obj.r.*obj.T);
                    
                    % First derivative:
                    H = diag(1./diff(K)',1)-diag([1./diff(K)' 0]);
                    
                    d = H*P;
                    
                    d = min(d,zeros(size(d)));
                    
                    
                    k = round(0.02*n);
                    
                    med = zeros(n,1);
                    
                    for i = 1:(n-k)
                        med(i) = median(d(i+1:i+k));
                    end
                    
                    % Smooth the calls in a sensible manner.
                    
                    d(1) = max(-exp(-obj.r.*obj.T),min(d(1),min(med(i),0)));
                    
                    for i = 2:(n-1)
                        d(i) = max(d(i-1),min(d(i),min(med(i),min(0,min((P(i+1)-P(i))/(K(i+1)-K(i)),(P(n)-P(i))/(K(n)-K(i)))))));
                        P(i+1) =  P(i) + d(i)*(K(i+1)-K(i));
                    end
                    
                    if (max(A*P-b)>eps)
                        warning('HM2','Problem with finding Feasible Call Prices');
                        disp(max(A*P-b));
                    end
                    
                end
                
                
                
                obj.computed_num = obj.computed_num+1;
                obj.computed_input{obj.computed_num} = K;
                obj.computed_output{obj.computed_num} = P;
                obj.computed_type{obj.computed_num} = 'Call';
                
            else
                warning('HPM:CLS1','Strikes must be a vector');
                P = 0;
            end
            
        end
        
        function P = Put(obj,K)
            if isvector(K)
                z = find(cellfun(@(q) (isequal(q,K)),obj.computed_input)==1,1);
                if ~isempty(z)
                    if strcmp(obj.computed_type{z},'Put')
                        P = obj.computed_output{z};
                        return;
                    elseif strcmp(obj.computed_type{z},'Call')
                        P = obj.computed_output{z} - obj.s + exp(-obj.r * obj.T)*K;
                        return;
                    else
                        error('HPM:Err2','Unknown type');
                    end
                end
                
                C = obj.Call(K);
                P = C - obj.s + K*exp(-obj.r*obj.T);
                
                obj.computed_num = obj.computed_num+1;
                obj.computed_input{obj.computed_num} = K;
                obj.computed_output{obj.computed_num} = P;
                obj.computed_type{obj.computed_num} = 'Put';
            else
                warning('HPM:PTS1','Strikes must be a vector');
                P = 0;
            end
        end

        function [p,iv] = simulate(obj,t)
            n = length(t);
            dt = diff(t);
            iv = zeros(size(t));

            v = obj.V0;
            rn1 = randn(1,n);
            rn2 = obj.rho*rn1 +sqrt(1-obj.rho^2)*randn(1,n); 
            for i=2:n
                v = max((v+obj.kappa*(obj.theta-v)*dt(i-1) + rn1(i-1)*sqrt(dt(i-1)*v)*obj.xi),0); % stop at zero.
                iv(i) = iv(i-1)+v*dt(i-1);
            end
            
            p = obj.s*exp(cumsum([0 (sqrt(diff(iv)).*rn2(2:n))]) + obj.r*t-iv/2);
            
        end

        function [p] = VaroptionPrice(obj,F)
            % Price a function of integrated variance. Method used is
            % inversion of MGF using the Abate and Whitt algorithm as
            % described in Leblanc and Scaillet (1998).
            
            ep = 0.00001; % Something small
            m = (1-obj.VarMGF(ep))/ep; % Mean of integrated variance
            
           
            % Compute a sensible grid
            M = 400; % Number of points between 0 and m in grid
            N = 14; % Multiples of m for grid length
            x = linspace(0,N*m,ceil(N*M));
            
            % Abate and Whitt algorithm parameters: Leblanc and Scaillet
            % reccomend A = 40, n = 1000
            A = 40;
            n = 4000; % 22/5/12: 1000 too small; 10000 works, but slow. Try 4000.
            
            f = zeros(size(x));
            k = 1:n;
            for i = 2:length(x)
                f(i) = exp(A/2)/x(i)*(real(obj.VarMGF(A/(2*x(i))))/2 + (-1).^k *real(obj.VarMGF((A+2*pi*1i*k')/(2*x(i)))));
            end
            
            % Check that the density is plausible
            tol = 0.025;
            test  = trapz(x,f);
            if (abs(test-1) > tol)
                display('Variance option price may be inaccurate. Density sums to:')
                display(test)
            end
                        
            p = trapz(x,f.*F(x));
                       
        end
            
        function psi = VarCharFn(obj,q)
            
            z = -(sqrt(obj.kappa^2-2*(obj.xi^2)*1i*q));
            
            psi = ((exp(obj.kappa*obj.T/2)./(cosh(z*obj.T/2) + (obj.kappa./z).*sinh(z*obj.T/2))) .* exp(obj.xi^2/(2*obj.kappa*obj.theta)*1i*q*obj.V0 * 2 .* sinh(z*obj.T/2)./(z .* (cosh(z*obj.T/2)+obj.kappa.*sinh(z*obj.T/2)./z)))).^(2*obj.kappa*obj.theta/(obj.xi^2));
        end

        
        function psi = VarMGF(obj,x)

            % See Leblanc & Scaillet for details:
            phi = obj.kappa*obj.theta;
            lambda = obj.kappa;
            alpha = obj.xi^2;
            
            mu = lambda * sqrt(1+2*x*alpha/lambda^2);
            k = (mu - lambda)./(mu+lambda);
            
            A = (1+k)./mu .* (1-exp(-mu*obj.T))./(1+k.*exp(-mu*obj.T));
            b = 2*phi./alpha.*log((1+k.*exp(-mu*obj.T))./(1+k)) + obj.T.*phi.*(mu-lambda)./alpha;
            
            psi = exp(-A.*x.*obj.V0 - b);
        
        end

        function S = saveobj(obj)
            obj.computed_input = {};
            obj.computed_output = {};
            obj.computed_type = {};
            obj.computed_num=0;
            S = obj;
        end

        
    end


    methods (Access=protected)
        function int = HVI1(obj,z,K)
            
            % Implement functions in Janek et al.
            b = obj.kappa - obj.xi*obj.rho;
            
            d = sqrt((obj.rho*obj.xi*z*1i - b).^2 - obj.xi^2.*(1i*z-z.^2));
            g = (b-obj.rho*obj.xi*z*1i+d)./(b-obj.rho*obj.xi*z*1i-d);
            
            C = obj.r*z*1i*obj.T+obj.kappa*obj.theta*((b-obj.rho*obj.xi*z*1i-d)*obj.T - 2*log((exp(-d*obj.T)-g)./(1-g)))/obj.xi^2;
            D = ((b-obj.rho*obj.xi*z*1i+d)/obj.xi^2).*((exp(-d*obj.T)-1)./(exp(-d*obj.T)-g));
            
            f = exp(C+D*obj.V0+1i*z*log(obj.s));
            
            int = real(exp(-1i*z.*log(K)).*f./(1i*z));
        
        end
    
        function int = HVI2(obj,z,K)
            
            % Implement functions in Janek et al.
            b = obj.kappa;
            
            d = sqrt((obj.rho*obj.xi*z*1i - b).^2 + obj.xi^2*z.*(1i+z));
            g = (b-obj.rho*obj.xi*z*1i+d)./(b-obj.rho*obj.xi*z*1i-d);
            
            C = obj.r*z*1i*obj.T+obj.kappa*obj.theta*((b-obj.rho*obj.xi*z*1i-d)*obj.T - 2*log((exp(-d*obj.T)-g)./(1-g)))/obj.xi^2;
            D = ((b-obj.rho*obj.xi*z*1i+d)/obj.xi^2).*((exp(-d*obj.T)-1)./(exp(-d*obj.T)-g));
            
            f = exp(C+D*obj.V0+1i*z*log(obj.s));
            
            int = real(exp(-1i*z.*log(K)).*f./(1i*z));
        
        end       
            
            
    end
end
