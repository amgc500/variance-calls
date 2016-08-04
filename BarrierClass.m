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

classdef BarrierClass < handle
    %BARRIERCLASS Holder for general barrier functions
    %   This contains functions which are used by both Rost and Root
    %   solutions
    properties (Abstract, SetAccess = immutable)
        PM % Pricing Model
        x % log-asset price
        z % asset price
        Cl % Compute Calls
        r % Computed barrier
        time_vec % vector of times
        n_low % Lower index for sensible barriers/prices
        n_high % Upper index for sensible barriers/prices
        
    end
    
    properties (Abstract, SetAccess = protected)
        n_samp; % Number of samples computed so far
        XTau_samp; % X values of stopped process
        TauD_samp; % time values of stopped process
        EmpPotential % Empirical potential computed using samples
        EmpCall % Empirical call prices computed using samples
        G % Surface of hedge function
        H % Static holding
        CurrentFn % Function used for last call to construct hedge
        PlotType; % 0 == time on x-axis for barriers, 1 == space on x-axis
    end
        
    
    methods (Abstract)
        
        [XTau,TauD] = simulate(n)
        [G,H] = ComputeHedge(f)
        g = interG(z,t)
        g = interGDash(z,t)
        h = interH(z)
        bool = IsInBarrier(z,t)
    end
    
    methods
        function obj = PlotBarrier(obj,col)
            
            if nargin ==1
                col = 12;
            end
            X = jet(16);
            
            if obj.PlotType == 0
                plot(obj.r(obj.n_low:obj.n_high),obj.z(obj.n_low:obj.n_high),'Color',X(col,:))
            else
                plot(obj.z(obj.n_low:obj.n_high),obj.r(obj.n_low:obj.n_high),'Color',X(col,:))
            end               
        end

        function obj = PlotBarrierLogScale(obj,col)
            if nargin ==1
                col = 12;
            end
            X = jet(16);
            
            if obj.PlotType == 0
                plot(obj.r(obj.n_low:obj.n_high),obj.x(obj.n_low:obj.n_high),'Color',X(col,:))
            else
                plot(obj.x(obj.n_low:obj.n_high),obj.r(obj.n_low:obj.n_high),'Color',X(col,:))
            end
        end

        function obj = PlotEmpCall(obj,col)
            if nargin ==1
                col = 12;
            end
            X = jet(16);
            
            plot(obj.z(obj.n_low:obj.n_high),obj.EmpCall(obj.n_low:obj.n_high),'Color',X(col,:))
            hold on
            plot(obj.z(obj.n_low:obj.n_high),obj.Cl(obj.n_low:obj.n_high),'Color',X(8,:))
        end        

        
        function obj = PlotImpliedDensity(obj)
            if isempty(obj.EmpPotential)
                obj.simulate(100);
                display('Run 100 simulations.')
            end
            finish = -2*obj.Cl + obj.PM.s - obj.z;
            hold off;
            plot(obj.x(2:length(obj.x)-1),-diff(diff(obj.EmpPotential)./diff(obj.z))./(diff(obj.z(1:length(obj.z)-1)+obj.z(2:length(obj.z)))))
            hold on;
            plot(obj.x(2:length(obj.x)-1),diff(diff(-finish)./diff(obj.z))./(diff(obj.z(1:length(obj.z)-1)+obj.z(2:length(obj.z)))))
            hold off;
        end
        
        function obj = CompareImpliedVol(obj,col)
            if nargin ==1
                col = 12;
            end
            X = jet(16);

            Cl_sim = zeros(size(obj.x)); % NB: more accurate to adjust mean to fit true mean.
            m = obj.PM.s-mean(obj.XTau_samp);
            for j = 1:length(obj.z)
                Cl_sim(j) = sum(max(obj.XTau_samp+m-obj.z(j),0))/length(obj.XTau_samp);
            end
    
            iv_Cl_sim = zeros(length(obj.z),1);
            for i = 1:length(obj.z)
                iv_Cl_sim(i) = blsimpv(obj.PM.s,obj.z(i)*exp(-obj.PM.r*obj.PM.T),0,obj.PM.T,Cl_sim(i));
            end

            iv_Cl = zeros(length(obj.z),1);
            for i = 1:length(obj.z)
                iv_Cl(i) = blsimpv(obj.PM.s,obj.z(i)*exp(-obj.PM.r*obj.PM.T),0,obj.PM.T,obj.Cl(i));
            end

            plot(obj.z(obj.n_low:obj.n_high),iv_Cl_sim(obj.n_low:obj.n_high),'Color',X(col,:))
            hold on
            plot(obj.z(obj.n_low:obj.n_high),iv_Cl(obj.n_low:obj.n_high),'Color',X(8,:))
            hold off
        end
        
        function p = price(obj,f)
            % Compute the price of the option with payoff F, where F' = f
            % and F(0) = 0.
            % NB: assumes that price is linear beyond end of time-range,
            % and in x-range, go to zero sufficiently quickly that payoff
            % outside the range is approximately zero.
            obj.ComputeHedge(f);
            
            df = exp(-obj.PM.r*obj.PM.T);
            n = obj.space_steps;
                                    
            % Compute price of a wavelet-type option which pays one at
            % Z/df)
            
            h_wavelet = obj.Cl(1:(n-2))./diff(obj.z(1:n-1))/df...
                - (obj.z(3:n)-obj.z(1:(n-2)))./diff(obj.z(2:n))./diff(obj.z(1:(n-1)))/df.*obj.Cl(2:(n-1)) ...
                + obj.Cl(3:n)./diff(obj.z(2:n))/df;
                        
            p = obj.H(2:n-1)'*h_wavelet + obj.interG(obj.PM.s,0);
        end
        
        function [p,t,iv,h_a,h_c,delta] = SimulateHedge(obj,f,n)
            
            % h_a is the traded portfolio plus initial case, h_c is the
            % intrinsic value of the payoff of calls.
            
            [p,t,iv,h_a,h_c,delta] = SimulateHedgeLong(obj,f,n,1);
            
        end
        
        function [Pyff,h_f,iv_f] = SimulateHedgeFinal(obj,f,n,N)
            % Return the payoff, final hedge value, final int var for N
            % runs of length n.
            
            Pyff = zeros(N,1);
            h_f = zeros(N,1);
            iv_f = zeros(N,1);
            
            for i = 1:N
                [~,~,iv,h_a,h_c,~] = obj.SimulateHedgeLong(f,n,1);
                Pyff(i) = trapz(iv,f(iv));
                h_f(i) = h_a(n) + h_c(n);
                iv_f(i) = iv(n);
            end
        end
        
        function [p,t,iv,h_a,h_c,delta] = SimulateHedgeLong(obj,f,n,T_mult)
            
            % As above, allow longer simulations (beyond T) tau_mult specifies the factor.
            
            % h_a is the traded portfolio plus initial cash, h_c is the
            % intrinsic value of the payoff of calls.
            
            
            p0 = obj.price(f);
            
            display('Price is:')
            display(p0)
            
            t = linspace(0,obj.PM.T*T_mult,n);
            [p,iv] = obj.PM.simulate(t);
            [h_a,h_c,delta] = obj.PathwiseSimulateHedgeLong(f,p,t,iv);
            
        end

        function [h_a,h_c,delta] = PathwiseSimulateHedgeLong(obj,f,p,t,iv)
            obj.price(f);
            
            h_a = zeros(size(t)) + obj.interG(p(1),0);
            h_c = zeros(size(t));
            
            delta = zeros(size(t));
            h_c(1) = obj.interH(exp(-obj.PM.r*t(1))*p(1));
            
            for i=2:length(t)
                delta(i-1) = obj.interGDash(p(i-1),iv(i-1));
                h_a(i) = h_a(i-1) + delta(i-1)*(p(i)-p(i-1));
                h_c(i) = obj.interH(exp(-obj.PM.r*t(i))*p(i));
            end
            
            delta(length(delta)) = obj.interGDash(p(length(delta)),iv(length(delta)));
        end            
                  
        function TestHedgeLong(obj,f,n,T_mult)
            [p,t,iv,h_a,h_c,delta] = obj.SimulateHedgeLong(f,n,T_mult);
            
            g = zeros(1,n);
            for i = 1:n
               g(i) = obj.interG(p(i),iv(i));
            end
            
            F = cumtrapz(iv,f(iv));

            bar = zeros(length(t));
            for i = 1:length(t)
                bar(i) = obj.IsInBarrier(p(i)*exp(-obj.PM.r*t(i)),iv(i))*1.0;
            end
            scale = h_a(length(h_a))-g(length(g));
            

            display('Computed Hedge (in time)')
            plot(t,p,t,iv,t,h_a+h_c,t,F)
            pause
            display('Perfect Super/Subhedge (in time)')
            plot(t,h_a+h_c,t,g+h_c,t,F)
            pause
            display('Perfect Super/Subhedge (in int var)')
            plot(iv,h_a+h_c,iv,g+h_c,iv,F)
            pause
            display('Diff between hedge and G. Should be >=0 (Rost)/ <=0 (Root)')
            plot(t,h_a-g)
            pause
            display('Diff between hedge and G. Should be >=0 (Rost)/ <=0 (Root) and IsInBarrier')
            plot(iv,h_a-g,iv,bar*scale)
            pause
            display('Price and Barrier')
            plot(obj.r,obj.z,iv,p)
            pause
            display('Delta')
            plot(t,delta)
            
            
        end
        
        function SetPlotType(obj,val)
            obj.PlotType = val;
        end

        function S = saveobj(obj)
            EP = obj.EmpPotential;
            EC = obj.EmpCall;
            obj.EmpPotential = [];
            obj.EmpCall = [];
            obj.G = [];
            obj.H = [];
            obj.CurrentFn = [];
            S = obj;
            obj.EmpPotential = EP;
            obj.EmpCall = EC;
        end
        
%         function obj = loadobj(obj)
%             obj.n_samp = n;
%             obj.XTau_samp = XTau;
%             obj.TauD_samp = TauD;
%             obj.EmpPotential = zeros(size(obj.x));
%             for j = 1:length(obj.z)
%                 obj.EmpPotential(j) = -sum(abs(XTau-obj.z(j)))/length(XTau);
%             end
%             obj.EmpCall = zeros(size(obj.x));
%             for j = 1:length(obj.z)
%                 obj.EmpCall(j) = sum(max(XTau-exp(-obj.PM.r*obj.PM.T)*obj.z(j),0))/length(XTau);
%             end
%         end 
            
    end
end

