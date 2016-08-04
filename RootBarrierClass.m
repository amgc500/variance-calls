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

classdef RootBarrierClass < handle & BarrierClass
    %ROOTBARRIERCLASS Find Root Barrier of given Call prices
    %   Implementation of the Root Barrier, and associated functions.
    
    properties (SetAccess = immutable)
        PM % Pricing Model
        x % log-asset price
        z % asset price
        space_steps % Number of space steps
        time_steps % Number of time steps
        time_vec % vector of times
        Cl % Compute Calls
        r % Computed barrier
        Sind % Index of initial value (so PM.s = z(Sind))
        atmVol % At the money volatility computed from option prices. Used for computing domains etc
        n_low % Lower index for sensible barriers/prices
        n_high % Upper index for sensible barriers/prices
        Kushner % ==1 if Kushner approximation to diffusion used. Also truncates the boundary of the call prices.
        u_plot_points % Number of points at which u is saved for plotting
        u_save % value of u at saved points
    end

    properties (SetAccess = protected)
        n_samp; % Number of samples computed so far
        XTau_samp; % X values of stopped process
        TauD_samp; % time values of stopped process
        EmpPotential % Empirical potential computed using samples
        EmpCall % Empirical call prices computed using samples
        G % Surface of hedge function
        H % Static holding
        CurrentFn % Function used for last call to construct hedge
        PlotType % Which way round to plot barriers
    end

    methods
        function obj = RootBarrierClass(PM,space_size,time_size,Kushner,upp)
            if nargin ==0
                obj.PM = HestonModel();
                obj.space_steps = 1000;
                obj.time_steps = 1000;
                obj.Kushner = 1;
                obj.u_plot_points = 8;
            elseif nargin ==1
                obj.PM = PM;
                obj.space_steps = 1000;
                obj.time_steps = 1000;
                obj.Kushner = 1;
                obj.u_plot_points = 8;
            elseif nargin ==2
                obj.PM = PM;
                obj.space_steps = space_size;
                obj.time_steps = 1000;
                obj.Kushner = 1;
                obj.u_plot_points = 8;
            elseif nargin ==3
                obj.PM = PM;
                obj.space_steps = space_size;
                obj.time_steps = time_size;
                obj.Kushner = 1;
                obj.u_plot_points = 8;
            elseif nargin ==4
                obj.PM = PM;
                obj.space_steps = space_size;
                obj.time_steps = time_size;
                obj.Kushner = Kushner;
                obj.u_plot_points = 8;                            
            else
                obj.PM = PM;
                obj.space_steps = space_size;
                obj.time_steps = time_size;
                obj.Kushner = Kushner;
                obj.u_plot_points = upp;                            
            end
            obj.n_samp = 0;
            obj.EmpPotential = [];
            obj.EmpCall = [];
            obj.G = [];
            obj.H = [];
            obj.CurrentFn = [];
            obj.PlotType = 0;
            
            S = obj.PM.s;

            atmCall = obj.PM.Call(S*exp(obj.PM.r*obj.PM.T));
            obj.atmVol = blsimpv(S,S*exp(obj.PM.r*obj.PM.T),obj.PM.r,obj.PM.T,atmCall);
            
            S_factor_up = (1.5*exp((4*obj.atmVol-(obj.atmVol^2)/2)*obj.PM.T));
            S_factor_down = (1.5*exp((4*obj.atmVol+(obj.atmVol^2))*obj.PM.T));
            S_min = S/S_factor_down;
            S_max = S*S_factor_up;
            n = obj.space_steps;
            m = obj.time_steps;
            
            
            n0 = floor((log(S)-log(S_min))/(log(S_max)-log(S_min))*(n-1));
            S_max = exp(log(S_min) + (log(S)-log(S_min))*(n-1)/(n0-1));
            obj.x = linspace(log(S_min),log(S_max),n)'; % log asset price
            obj.z = exp(obj.x); % Asset prices
            obj.Sind = n0;

            
            obj.n_low = find((obj.z>S*exp(-4*obj.atmVol-(obj.atmVol^2)/2)),1);
            obj.n_high = find((obj.z>S*exp(3*obj.atmVol-(obj.atmVol^2)/2)),1);
            
            
            % Compute possible initial and final potentials. Vary these for different
            % measures.
            start = -abs(obj.z-S);
            
            display('Computing Prices...')
            obj.Cl = obj.PM.Call(obj.z*exp(-obj.PM.r*obj.PM.T));
            display('...Done')

            % Use call prices to estimate a suitable upper time limit
            T = 10*(obj.atmVol^2)*obj.PM.T ; % Upper bound on integrated variance

            dt = T/obj.time_steps; % Size of timestep
            t = linspace(0,T,obj.time_steps+1);
            obj.time_vec = t;
            

            finish = -2*obj.Cl + S - obj.z;
  
            if (obj.Kushner == 1) %Make the measure have bounded support
                k = 2;
                while ((start(1) + (finish(k)-start(1))*(obj.z(k+1)-obj.z(1))/(obj.z(k)-obj.z(1)) - finish(k+1) <0)&& (k < n))
                    k = k+1;
                end
                for j = 1:k
                    finish(j) = start(1) + (finish(k)-start(1))*(obj.z(j)-obj.z(1))/(obj.z(k)-obj.z(1));
                end
                
                k = n-1;
                while ((start(n) + (finish(k)-start(n))*(obj.z(k-1)-obj.z(n))/(obj.z(k)-obj.z(n)) - finish(k-1) <0)&& (k > 1))
                    k = k-1;
                end
                for j = n:-1:k
                    finish(j) = start(n) + (finish(k)-start(n))*(obj.z(j)-obj.z(n))/(obj.z(k)-obj.z(n));
                end
            end
           
            % Control numerical solution: theta = 0, explicit, =1 implicit = 1/2
            % Crank-Nicolson.
            theta = 0.5; %Note: will be redefined in Kushner method!
            
            % Compute the boundary corrections to the differentiation matrices
            D1 = [-1 1 0 0 0];
            D2 = [1 -2 1 0 0];
            D3 = [-1 3 -3 1 0];
            D4 = [1 -4 6 -4 1];
            Hmat = [D1;D2;D3;D4];

            k = 0:4;
            A = [k.^1;k.^2/2;k.^3/3;k.^4/4]';

            v1 = [1 0 0 0]/(Hmat*A); % For computing boundary estimate of first derivative
            v2 = [0 1 0 0]/(Hmat*A); % For computing boundary estimate of second derivative

            % Diffusion and drift coefficients:
            c = spdiags(ones(size(obj.x))/2,0,n,n); %Diffusion

            dx = obj.x(2)-obj.x(1);

            if (obj.Kushner == 1) %Constant on boundary, and Kushner approximation

                % Drift is minus 1/2 always.
                % First and second derivatives as sparse matrices:
                Pminus = spdiags([ones(n,1) -ones(n,1)]/(2*dx),[-1 0],n,n);
                Pminus(n,n) = 0;
                Pminus(n,n-1) = 0;
                Pminus(1,1) = 0;
                
                
                
                Q = spdiags([ones(n,1) -2*ones(n,1) ones(n,1)]/(dx^2),[-1 0 1],n,n);
                c(1,1) = 0;
                c(n,n) = 0;
                
                
                D = (c*Q + Pminus);

                % CFL Condition:
                %theta = max(1-1/(2*max(diff(obj.time_vec))*max(-diag(D))),theta)
                theta  = 1; % Use implicit time step for unconditional stability
                
                CFL = max(diff(obj.time_vec))*(1-theta)*max(-diag(D));
                if CFL>=1
                    warning('CFL condition fails');
                    disp(CFL);
                end

                
                
            else
                
                d = spdiags(-(ones(size(obj.x))/2),0,n,n); %Drift

                % First and second derivatives as sparse matrices:
                P = spdiags([-ones(n,1) zeros(n,1) ones(n,1)]/(2*dx),[-1 0 1],n,n);
                P(1,1:5) = v1*Hmat/dx;
                P(n,n:-1:n-4) = -v1*Hmat/dx;
                
                % Second derivative
                Q = spdiags([ones(n,1) -2*ones(n,1) ones(n,1)]/(dx^2),[-1 0 1],n,n);
                Q(1,1:5) = v2*Hmat/dx^2;
                Q(n,n:-1:n-4) = v2*Hmat/dx^2;
                
                D = c*Q + d*P;
            end
            
            % Set up value of potentials
            u0 = start;
            u = zeros(n,m+1);
            u(:,1) = u0;
            
            
            %Set up to solve du/dt = (c*Q + d*P)u
            if obj.Kushner == 1
                M = (speye(n)-dt*D);
                v = M*finish;
                z = zeros(size(start));
                
            else
                R = speye(n) + (1-theta)*dt*D;
                L = speye(n) - theta*dt*D;
                RL = L\R;
            end

            display('Computing Barrier...')            

            for i=1:m
                if obj.Kushner == 1
                    [u(:,i+1),z] = CPA(M,v-u(:,i),z);
                    u(:,i+1) = u(:,i+1) + finish;                    
                else
                    u(:,i+1) = RL*u(:,i);
                    u(:,i+1) = max(u(:,i+1),finish);
                    u(:,i+1) = min(u(:,i),u(:,i+1)); % enforce decreasing
                end
            end

            display('...done.')
 
            % Save some values for later plotting. Conentrate towards
            % smaller values of t.
            obj.u_save = zeros(n,obj.u_plot_points);
            
            u_t = [1 ceil((obj.time_steps-1)*((1:(obj.u_plot_points-1)).^2)/((obj.u_plot_points-1)^2))];
            
            for j = 1:obj.u_plot_points
                obj.u_save(:,j) = u(:,u_t(j));
            end
                       
            % compute the barrier.
            obj.r = zeros(1,n);
 
            for k = 1:n
                p = find(u(k,:) == finish(k));
                if  numel(p) == 0
                    obj.r(k) = T;
                else
                    q = p(1);
                    obj.r(k) = t(q);
                end
            end
        end
        
        function [XTau,TauD] = simulate(obj,n)
            if obj.n_samp >= n
                XTau = obj.XTau_samp(1:n);
                TauD = obj.XTau_samp(1:n);
            else
                TauD = (1:n)*0+obj.time_vec(length(obj.time_vec));
                TauD(1:obj.n_samp) = obj.TauD_samp;
                XTau = (1:n)*0;
                XTau(1:obj.n_samp) = obj.XTau_samp;
                
                m = length(obj.time_vec);
                Y = (1:m+1)*0;
                
                t_step = obj.time_vec;

                for k = (obj.n_samp+1):n
                    Y(1) = 0; %Sample log process
                    for i = 2:length(t_step)
                        dt = t_step(i)-t_step(i-1);
                        Y(i) = Y(i-1)+sqrt(dt)*randn(1)-0.5*dt;
                        if obj.IsInBarrier(obj.PM.s*exp(Y(i)),t_step(i))
                            TauD(k) = t_step(i);
                            XTau(k) = obj.PM.s*exp(Y(i));
                            break
                        else
                            continue
                        end
                    end
                end
                obj.n_samp = n;
                obj.XTau_samp = XTau;
                obj.TauD_samp = TauD;
                obj.EmpPotential = zeros(size(obj.x));
                for j = 1:length(obj.z)
                    obj.EmpPotential(j) = -sum(abs(XTau-obj.z(j)))/length(XTau);
                end
                obj.EmpCall = zeros(size(obj.x));
                for j = 1:length(obj.z)
                    obj.EmpCall(j) = sum(max(XTau-exp(-obj.PM.r*obj.PM.T)*obj.z(j),0))/length(XTau);
                end                
            end
        end

        function [G,H] = ComputeHedge(obj,f)
            % Compute the hedge of an option with payoff F = f', F(0) = 0
            % Assume that the second derivative of f is constant beyond
            % T.
            
            % e.g. ComputeHedge(@(t) t.^2)
            %  or  ComputeHedge(@(t) max(t-k,0))
            if(isempty(obj.CurrentFn)||~strcmp(func2str(f),func2str(obj.CurrentFn))||any(~(f(obj.z)==obj.CurrentFn(obj.z))))
                t = obj.time_vec;
                n = obj.space_steps;
                m = length(t);
                X = obj.x;
                f_t = 1*f(t);
                M = zeros(length(X),length(t));
                
                
                % Compute the boundary corrections to the differentiation matrices
                D1 = [-1 1 0 0 0];
                D2 = [1 -2 1 0 0];
                D3 = [-1 3 -3 1 0];
                D4 = [1 -4 6 -4 1];
                Hmat = [D1;D2;D3;D4];
                
                k = 0:4;
                A = [k.^1;k.^2/2;k.^3/3;k.^4/4]';
                
                v1 = [1 0 0 0]/(Hmat*A); % For computing boundary estimate of first derivative
                v2 = [0 1 0 0]/(Hmat*A); % For computing boundary estimate of second derivative


                % First and second derivatives as sparse matrices:
                dx = obj.x(2)-obj.x(1);
                P = spdiags([-ones(n,1) zeros(n,1) ones(n,1)]/(2*dx),[-1 0 1],n,n);
                P(1,1:5) = v1*Hmat/dx;
                P(n,n:-1:n-4) = -v1*Hmat/dx;
            
                % Second derivative
                Q = spdiags([ones(n,1) -2*ones(n,1) ones(n,1)]/(dx^2),[-1 0 1],n,n);
                Q(1,1:5) = v2*Hmat/dx^2;
                Q(n,n:-1:n-4) = v2*Hmat/dx^2;

                % Diffusion and drift coefficients:
                c = spdiags(-ones(size(obj.x))/2,0,n,n); %Diffusion
                d = spdiags((ones(size(obj.x))/2),0,n,n); %Drift
                
                M(:,m) = f_t(m);
            
                %Compute hedge:
                dt = t(m)-t(m-1);
                D = dt*(c*Q+d*P)/2; % Factor of 1/2 = Crank-Nicholson
                
                display('Computing Hedge...')
                
                for i = (m-1):-1:1
                    ind = spdiags((obj.r > t(i))',0,n,n);
                    M(:,i) = (speye(n) + ind*D)\(ind*(speye(n)-D)*M(:,i+1)-ind*D*(speye(n)-ind)*ones(n,1)*f_t(i)) +(speye(n)-ind)*ones(n,1)*f_t(i);
                end
                
                display('...done.')
                
                zd = obj.z;
                Z = zeros(n,1);
                
                Z(obj.Sind:n) = cumtrapz(zd(obj.Sind:n),cumtrapz(zd(obj.Sind:n),M(obj.Sind:n,1)./(zd(obj.Sind:n)).^2))*2;
                Z(obj.Sind:-1:1) = cumtrapz(zd(obj.Sind:-1:1),cumtrapz(zd(obj.Sind:-1:1),M(obj.Sind:-1:1,1)./(zd(obj.Sind:-1:1)).^2))*2;
                Mint = cumtrapz(t,M')';
                G = Mint - Z*ones(1,length(t));
                
                Ind = (obj.r'*ones(1,length(t)) > ones(length(obj.r),1)*t);
                
                H = trapz(t,((ones(n,1)*f_t-M).*Ind)')' +Z;
                
                obj.G = G;
                obj.H = H;
                obj.CurrentFn = f;
            else
                G = obj.G;
                H = obj.H;
            end
        end
        
        function g = interG(obj,z,t)
            % return an interpolated value of G at asset price z, time t
            i1 = find((obj.time_vec >=t),1);
            if i1 == 1
                i1 = 1;
                tfrac = 0;
            elseif isempty(i1)
                i1 = length(obj.time_vec)-1;
                tfrac = 1;
            else
                i1 = i1-1;
                tfrac = (t-obj.time_vec(i1))/(obj.time_vec(i1+1)-obj.time_vec(i1));
            end
            
            X = log(z);
            j1 = find((obj.x >=X),1);
            if j1 == 1
                j1 = 1;
                jfrac = 0;
            elseif isempty(j1)
                j1 = length(obj.x)-1;
                jfrac = 1;
            else
                j1 = j1-1;
                jfrac = (X-obj.x(j1))/(obj.x(j1+1)-obj.x(j1));
            end
            
            g = (1-tfrac)*((1-jfrac)*obj.G(j1,i1) + jfrac*obj.G(j1+1,i1)) + ...
                tfrac*((1-jfrac)*obj.G(j1,i1+1) + jfrac*obj.G(j1+1,i1+1));
        end

        function g = interGDash(obj,z,t)
            % return an interpolated value of GDash (derivative wrt discounted asset price) at asset price z, time t
            i1 = find((obj.time_vec >=t),1);
            if i1 == 1
                i1 = 1;
                tfrac = 0;
            elseif isempty(i1)
                i1 = length(obj.time_vec)-1;
                tfrac = 1;
            else
                i1 = i1-1;
                tfrac = (t-obj.time_vec(i1))/(obj.time_vec(i1+1)-obj.time_vec(i1));
            end
            
            X = log(z);
            
            j1 = find((obj.x >=X),1);
            if j1 == 1
                j1 = 1;
                g = ((1-tfrac)*(obj.G(j1+1,i1)-obj.G(j1,i1))/(obj.z(j1+1)-obj.z(j1)) + tfrac*(obj.G(j1+1,i1+1)-obj.G(j1,i1+1))/(obj.z(j1+1)-obj.z(j1)));
            elseif isempty(j1)
                j1 = length(obj.x)-1;
                g = ((1-tfrac)*(obj.G(j1+1,i1)-obj.G(j1,i1))/(obj.z(j1+1)-obj.z(j1)) + tfrac*(obj.G(j1+1,i1+1)-obj.G(j1,i1+1))/(obj.z(j1+1)-obj.z(j1)));
            else
                j1 = j1-1;
                jfrac = (X-obj.x(j1))/(obj.x(j1+1)-obj.x(j1));
       
                if ((j1 == 1)&&(jfrac < 0.5))
                    jfrac = 0.5;
                end

                if ((j1 == length(obj.x)-1)&&(jfrac >= 0.5))
                    jfrac = 0.499;
                end
                
                if jfrac < 0.5
                    g = (0.5-jfrac)*((1-tfrac)*(obj.G(j1,i1)-obj.G(j1-1,i1))/(obj.z(j1)-obj.z(j1-1)) + tfrac*(obj.G(j1,i1+1)-obj.G(j1-1,i1+1))/(obj.z(j1)-obj.z(j1-1)))...
                        + (0.5+jfrac)*((1-tfrac)*(obj.G(j1+1,i1)-obj.G(j1,i1))/(obj.z(j1+1)-obj.z(j1)) + tfrac*(obj.G(j1+1,i1+1)-obj.G(j1,i1+1))/(obj.z(j1+1)-obj.z(j1)));
                else
                    g = (jfrac-0.5)*((1-tfrac)*(obj.G(j1+2,i1)-obj.G(j1+1,i1))/(obj.z(j1+2)-obj.z(j1+1)) + tfrac*(obj.G(j1+2,i1+1)-obj.G(j1+1,i1+1))/(obj.z(j1+2)-obj.z(j1+1)))...
                        + (1.5-jfrac)*((1-tfrac)*(obj.G(j1+1,i1)-obj.G(j1,i1))/(obj.z(j1+1)-obj.z(j1)) + tfrac*(obj.G(j1+1,i1+1)-obj.G(j1,i1+1))/(obj.z(j1+1)-obj.z(j1)));
                end
            end
        end

        
        function h = interH(obj,z)
            % return an interpolated value of H at asset price z
            
            X = log(z);
            j1 = find((obj.x >=X),1);
            if j1 == 1
                j1 = 1;
                jfrac = 0;
            elseif isempty(j1)
                j1 = length(obj.x)-1;
                jfrac = 1;
            else
                j1 = j1-1;
                jfrac = (X-obj.x(j1))/(obj.x(j1+1)-obj.x(j1));
            end
            
            h = ((1-jfrac)*obj.H(j1) + jfrac*obj.H(j1+1));
        end
        
        
        
        function obj = u_plot(obj)
            X = jet(obj.u_plot_points);

            hold off
                        
            for i = 1:obj.u_plot_points
                plot(obj.z(obj.n_low:obj.n_high),obj.u_save(obj.n_low:obj.n_high,i),'Color',X(i,:))
                hold on
            end
            hold off          
        end        
        
        
        function bool = IsInBarrier(obj,z,t)
                
            dx = obj.x(2)-obj.x(1);
            int_part = min(max(floor((log(z)-obj.x(1))/dx)+1,1),obj.space_steps-1);
            frac = min(max((log(z)-obj.x(int_part))/dx,0),1);
            r_crit = (1-frac)*obj.r(int_part) + frac*obj.r(int_part+1);
            if r_crit <= t
                bool = true;
            else
                bool = false;
            end
        end
            
            
    end
    
end

