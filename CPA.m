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

function [z,z0] = CPA(M,q,z0)
%CPA Complementary Pivot Algorithm
%   The Complementary Pivot Algorithm, as described by Murty, K.G., "Linear
%   Complementarity, Linear and Nonlinear Programming", available from:
%
%   http://ioe.engin.umich.edu/people/fac/books/murty/linear_complementarity_webbook/
%
%   The algorithm will use a pivoting method to solve the LCP: 
%      Mz + q = w,   z,w >= 0,   z.w = 0.
%
%   The required input data is the matrix M, and the vector q. An
%   additional optional argument, z0, represents a possible starting point
%   for the algorithm. Non-zero entries of z0 should correspond to the
%   expected non-zero entries of the vector z. 
%
%   The algorithm will return the value of z, and optionally the vector z0
%   which contains the non-zero values of z.


n = length(q); % Problem size

% If z0 not specified, assume z = 0 as initial solution
if nargin == 2
    z0 = zeros(size(q));
end

% Find basis co-ordinates
basicz = find(z0~=0);
basicw = find(z0 == 0);
b = length(basicz); % Number of basis elements

I = speye(n);

% Compute intial basic z,w:
M2 = [M(:,basicz) -I(:,basicw)];
y = -M2\q;
[mn,loc] = min(y);

z = zeros(size(q));
%w = zeros(size(q));
z0 = zeros(size(q));

if mn >= 0 % initial allocation is also feasible
    
    z(basicz) = y(1:b);
    % w(nonbasic) = y(b+1:n);
    z0(basicz) = ones(1,b);
else
    base = [basicz; -basicw]; % identity of basis variable in given column;    
    if loc <= b
        pivot = basicz(loc);
        last = 1; % basic z left
        % basicz = basicz(1:b ~= loc);
    else
        pivot = basicw(loc-b);
        last = -1; % basic w left
        % basicw = basicw(1:(n-b) ~= (loc-b));
    end
    base(loc) = 0;
    
    M2(:,loc) = -M2*ones(n,1);
    y = y - mn*ones(n,1);
    y(loc) = -mn;
    
    
    
    % Iterate, choosing new pivots.
    while last ~= 0 % Not yet chosen artificial
        if last == -1 % Last selected w basic to leave
            h = M(:,pivot);
        else
            h = zeros(size(q));
            h(pivot) = -1;
        end

        e = M2\h;
        if (max(e)<=0)
            error('Ray termination')
        end
        ind = find(e>0);
        theta = min(y(ind)./e(ind));
        ind = ind(y(ind)./e(ind) == theta);
        if ~isempty(find(base(ind)==0, 1)) % Choose the artificial basis vector if possible
            loc = ind(base(ind) == 0);
        else %Choose randomly
            loc = ind(randi(length(ind),1));
        end
        
        
        y = y - theta*e;
        y(loc) = theta;

        temp = last*pivot; % New basis vector
        last = sign(base(loc));
        pivot = base(loc)*last;
        base(loc) = -temp;
        M2(:,loc) = h;
        
    end
    
    z0(base(base>0))=1;
    z(base(base>0)) = y(base>0);
end

