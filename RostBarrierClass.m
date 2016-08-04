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

classdef RostBarrierClass < handle & BarrierClass
    %ROSTBARRIERCLASS Find Rost's Barrier of given Call prices
    %   Given certain initial data, this will compute the Rost barrier.
    %   The class also includes certain methods for computing hedging
    %   strategies.
    
    properties (SetAccess = immutable)
        PM % Heston Parameters
        x % log-asset price
        y % Equally spaces steps -> arctan transform of asset scale
        z % asset price
        space_steps % Number of space steps
        time_s % number of time-steps
        time_vec % vector of times
        time_step_change % Index where time step changes
        Cl % Compute Calls
        r % Computed barrier
        atmVol % At the money volatility computed from option prices. Used for computing domains etc
        n_low % Lower index for sensible barriers/prices
        n_high % Upper index for sensible barriers/prices
        t_diff % number of time-step regimes
        Kushner % 1 = use Kushner approimations for derivatives (see Barles & Jacobsen). This will also truncate any given data.
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
        function obj = RostBarrierClass(PM,space_size,time_size,time_changes,Kushner,upp)
            if nargin ==0
                obj.PM = HestonModel();
                obj.space_steps = 1000;
                obj.time_s = 1000;
                obj.t_diff = 3;
                obj.Kushner=1;
                obj.u_plot_points =8;
            elseif nargin ==1
                obj.PM = PM;
                obj.space_steps = 1000;
                obj.time_s = 1000;
                obj.t_diff = 3;
                obj.Kushner=1;
                obj.u_plot_points =8;
            elseif nargin ==2
                obj.PM = PM;
                obj.space_steps = space_size;
                obj.time_s = 1000;
                obj.t_diff = 3;
                obj.Kushner=1;
                obj.u_plot_points =8;
            elseif nargin ==3
                obj.PM = PM;
                obj.space_steps = space_size;
                obj.time_s = time_size;
                obj.t_diff = 3;
                obj.Kushner=1;
                obj.u_plot_points =8;
            elseif nargin ==4
                obj.PM = PM;
                obj.space_steps = space_size;
                obj.time_s = time_size;
                obj.t_diff = time_changes;
                obj.Kushner=1;
                obj.u_plot_points =8;
            elseif nargin ==5
                obj.PM = PM;
                obj.space_steps = space_size;
                obj.time_s = time_size;
                obj.t_diff = time_changes;
                obj.Kushner=Kushner;
                obj.u_plot_points =8;
            else
                obj.PM = PM;
                obj.space_steps = space_size;
                obj.time_s = time_size;
                obj.t_diff = time_changes;
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

            S_factor_up = (1.5*exp((6*obj.atmVol-(obj.atmVol^2)/2)*obj.PM.T));
            S_factor_down = (1.5*exp((6*obj.atmVol+(obj.atmVol^2))*obj.PM.T));
            S_min = S/S_factor_down;
            S_max = S*S_factor_up;

            n = obj.space_steps;
            
            m = obj.time_s;
            
            % The change of variables has been chosen to cluster the space points near
            % the starting point.
            obj.y = linspace(atan(log(S_min/S)),atan(log(S_max/S)),n)'; % Space for solving vi
            obj.x = tan(obj.y); % log asset price
            obj.z = exp(obj.x)*S; % Asset prices

            obj.n_low = find((obj.z>S*exp(-5*obj.atmVol-(obj.atmVol^2))/1.3),1);
            obj.n_high = find((obj.z>S*exp(3*obj.atmVol-(obj.atmVol^2)/2)),1);

            
            % Compute possible initial and final potentials. Vary these for different
            % measures.
            start = -abs(obj.z-S);
            
            m_vec = ceil((0:(obj.t_diff))*(m-1)/obj.t_diff+1);
            
            obj.time_step_change = m_vec;
            
            t_crit = [0 obj.atmVol*exp(linspace(log(2^(-obj.t_diff+1)),log(100),obj.t_diff))];
            obj.time_vec = zeros(m,1);
            
            for i = 1:obj.t_diff
                obj.time_vec(m_vec(i):m_vec(i+1)) = linspace(t_crit(i),t_crit(i+1),m_vec(i+1)-m_vec(i)+1);
            end
            
                        
            display('Computing Prices...')
            obj.Cl = obj.PM.Call(obj.z*exp(-obj.PM.r*obj.PM.T));
            display('...Done')
            
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
            c = spdiags((cos(obj.y).^4)/2,0,n,n);

            dy = obj.y(2)-obj.y(1);

            if (obj.Kushner == 1) %Constant on boundary, and Kushner approximation

                % First and second derivatives as sparse matrices:
                Pplus = spdiags([-ones(n,1) ones(n,1)]/(2*dy),[0 1],n,n);                
                Pminus = spdiags([ones(n,1) -ones(n,1)]/(2*dy),[-1 0],n,n);
                
                
                Q = spdiags([ones(n,1) -2*ones(n,1) ones(n,1)]/(dy^2),[-1 0 1],n,n);
                c(1,1) = 0;
                c(n,n) = 0;
                
                dplus =spdiags(max(-(sin(obj.y).*(cos(obj.y).^3) + (cos(obj.y).^2)/2),0),0,n,n);
                dminus =spdiags(max(0,-(-(sin(obj.y).*(cos(obj.y).^3) + (cos(obj.y).^2)/2))),0,n,n);
                dplus(1,1) = 0;
                dplus(n,n) = 0;
                dminus(1,1) = 0;
                dminus(n,n) = 0;
                
                
                D = (c*Q + dplus*Pplus + dminus*Pminus);

                % CFL Condition:
                %theta = max(1-1/(2*max(diff(obj.time_vec))*max(-diag(D))),theta)
                theta  = 1; % Use implicit time step for unconditional stability
                
                CFL = max(diff(obj.time_vec))*(1-theta)*max(-diag(D));
                if CFL>=1
                    warning('CFL condition fails');
                    disp(CFL);
                end

                
                
            else
                
                % First and second derivatives as sparse matrices:
                P = spdiags([-ones(n,1) zeros(n,1) ones(n,1)]/(2*dy),[-1 0 1],n,n);
                P(1,1:5) = v1*Hmat/dy;
                P(n,n:-1:n-4) = -v1*Hmat/dy;
                
                Q = spdiags([ones(n,1) -2*ones(n,1) ones(n,1)]/(dy^2),[-1 0 1],n,n);
                Q(1,1:5) = v2*Hmat/dy^2;
                Q(n,n:-1:n-4) = v2*Hmat/dy^2;
                
                d = spdiags(-(sin(obj.y).*(cos(obj.y).^3) + (cos(obj.y).^2)/2),0,n,n);
                D = (c*Q + d*P);

            end
                
            % Set up value of potentials
            u0 = finish-start;
            u = zeros(n,m);
            u(:,1) = finish-start;
            z0 = zeros(n,1);
            z0(Q*u(:,1)>0) = 1;
            
            display('Computing Barrier...')
                
            for i = 1:obj.t_diff
                
                dt = obj.time_vec(m_vec(i)+1)-obj.time_vec(m_vec(i));

                % Use a Complementary Pivot Algorithm to solve
                M = (speye(n) - dt*theta*D);
                for j = m_vec(i):m_vec(i+1)-1
                    z0 = CPA(M,-dt*D*u(:,j),z0);
                    u(:,j+1) = u(:,j)+z0;
                end
                
            end
                
            % Save some values for later plotting
            j = 1;
            obj.u_save = zeros(n,obj.u_plot_points);
            for i = 1:ceil((obj.time_s-1)/obj.u_plot_points):obj.time_s
                obj.u_save(:,j) = u(:,i);
                j = j+1;
            end
            
            
            % compute the reversed barrier.
            obj.r = zeros(1,n);
            T = max(obj.time_vec);
 
            for k = 1:n
                p = find(u(k,:) ~= u0(k));
                if  numel(p) == 0
                    obj.r(k) = T;
                else
                    q = p(1);
                    obj.r(k) = obj.time_vec(q);
                end
            end

            % 'Smooth' the reversed barrier. NB: remove for general starting measures
            if (u(:,1) == finish-abs(obj.z-S))
                S_ind =  find(obj.z>=S,1,'first');
                curr_max = obj.r(S_ind);
                for i = (S_ind+1):length(obj.x)
                    obj.r(i) = max(curr_max,obj.r(i));
                    curr_max = obj.r(i);
                end
                
                S_ind =  find(obj.z<=S,1,'last');
                curr_max = obj.r(S_ind);
                for i = (S_ind-1):-1:1
                    obj.r(i) = max(curr_max,obj.r(i));
                    curr_max = obj.r(i);
                end
            end
            
            display('...done.')
            
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
                Y = (1:m)*0;
                
                dy = obj.y(2)-obj.y(1);
                
                % Slight offset to make start smoother
                %m1 = obj.time_step_change(2);
                %dt = obj.time_vec(2)-obj.time_vec(1);
                
                %t_step = obj.time_vec + [0 ones(1,m1)*dt/2 ones(1,length(obj.time_vec)-1-m1)*0]';

                
                t_step = obj.time_vec;

                for k = (obj.n_samp+1):n
                    
                    Y(1) = 0; %Sample log process

                    for i = 2:length(t_step)
                        dt = t_step(i)-t_step(i-1);
                        Y(i) = Y(i-1)+sqrt(dt)*randn(1)-0.5*dt;
                        int_part = min(max(floor(atan(Y(i))/dy-obj.y(1)/dy)+1,1),obj.space_steps-1);
                        frac = (atan(Y(i))-obj.y(int_part))/dy;
                        r_crit = (1-frac)*obj.r(int_part) + frac*obj.r(int_part+1);
                        if r_crit > t_step(i)
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
            % time_s2.
            
            % e.g. ComputeHedge(@(t) t.^2)
            %  or  ComputeHedge(@(t) max(t-k,0))
            if(isempty(obj.CurrentFn)||~strcmp(func2str(f),func2str(obj.CurrentFn))||any(~(f(obj.x)==obj.CurrentFn(obj.x))))
                t = obj.time_vec';
                n = obj.space_steps;
                m = length(t);
                m_vec = obj.time_step_change;
                x0 = obj.x;
                f_t = 1*f(t);
                F = cumtrapz(t,f_t);
                M = zeros(length(x0),length(t));
                
                
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
                dy = obj.y(2)-obj.y(1);
                P = spdiags([-ones(n,1) zeros(n,1) ones(n,1)]/(2*dy),[-1 0 1],n,n);
                P(1,1:5) = v1*Hmat/dy;
                P(n,n:-1:n-4) = -v1*Hmat/dy;
                
                Q = spdiags([ones(n,1) -2*ones(n,1) ones(n,1)]/(dy^2),[-1 0 1],n,n);
                Q(1,1:5) = v2*Hmat/dy^2;
                Q(n,n:-1:n-4) = v2*Hmat/dy^2;
                
                % Diffusion and drift coefficients:
                c = -spdiags((cos(obj.y).^4)/2,0,n,n);
                d = spdiags((sin(obj.y).*(cos(obj.y).^3) + (cos(obj.y).^2)/2),0,n,n);
                
                M(:,m) = f_t(m);
            
                display('Computing Hedge...')
                
                for i = obj.t_diff:-1:1
                    dt = obj.time_vec(m_vec(i)+1)-obj.time_vec(m_vec(i));
                    D = dt*(c*Q+d*P)/2; % Factor of 1/2 = Crank-Nicholson

                    for j = (m_vec(i+1)-1):-1:(m_vec(i))
                        ind = spdiags((obj.r < t(j))',0,n,n);
                        M(:,j) = (speye(n) + ind*D)\(ind*(speye(n)-D)*M(:,j+1)-ind*D*(speye(n)-ind)*ones(n,1)*f_t(j)) +(speye(n)-ind)*ones(n,1)*f_t(j);
                    end
                end
                
                
                display('...done.')
                
                Z = -2*f_t(m).*log(obj.z./obj.PM.s);
                Mint = cumtrapz(t,M')';
                G = F(m) + Mint -(Mint(:,m)+ Z)*ones(1,length(t));
                
                Ind = (obj.r'*ones(1,length(t)) < ones(length(obj.r),1)*t);
                
                H = trapz(t,((M-ones(n,1)*f_t).*Ind)')' +Z;
                
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
            
            y0 = atan(log(z/obj.PM.s));
            j1 = find((obj.y >=y0),1);
            if j1 == 1
                j1 = 1;
                jfrac = 0;
            elseif isempty(j1)
                j1 = length(obj.y)-1;
                jfrac = 1;
            else
                j1 = j1-1;
                jfrac = (y0-obj.y(j1))/(obj.y(j1+1)-obj.y(j1));
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
            
            y0 = atan(log(z/obj.PM.s));
            
            j1 = find((obj.y >=y0),1);
            if j1 == 1
                j1 = 1;
                g = ((1-tfrac)*(obj.G(j1+1,i1)-obj.G(j1,i1))/(obj.z(j1+1)-obj.z(j1)) + tfrac*(obj.G(j1+1,i1+1)-obj.G(j1,i1+1))/(obj.z(j1+1)-obj.z(j1)));
            elseif isempty(j1)
                j1 = length(obj.y)-1;
                g = ((1-tfrac)*(obj.G(j1+1,i1)-obj.G(j1,i1))/(obj.z(j1+1)-obj.z(j1)) + tfrac*(obj.G(j1+1,i1+1)-obj.G(j1,i1+1))/(obj.z(j1+1)-obj.z(j1)));
            else
                j1 = j1-1;
                jfrac = (y0-obj.y(j1))/(obj.y(j1+1)-obj.y(j1));
       
                if ((j1 == 1)&&(jfrac < 0.5))
                    jfrac = 0.5;
                end

                if ((j1 == length(obj.y)-1)&&(jfrac >= 0.5))
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
            
            y0 = atan(log(z/obj.PM.s));
            j1 = find((obj.y >=y0),1);
            if j1 == 1
                j1 = 1;
                jfrac = 0;
            elseif isempty(j1)
                j1 = length(obj.y)-1;
                jfrac = 1;
            else
                j1 = j1-1;
                jfrac = (y0-obj.y(j1))/(obj.y(j1+1)-obj.y(j1));
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

            dy = obj.y(2)-obj.y(1);
            int_part = min(max(floor(atan(log(z))/dy-obj.y(1)/dy)+1,1),obj.space_steps-1);
            frac = (atan(log(z))-obj.y(int_part))/dy;
            r_crit = (1-frac)*obj.r(int_part) + frac*obj.r(int_part+1);
            if r_crit > t
                bool = true;
            else
                bool = false;
            end
        end
    end
end

