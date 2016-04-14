function [varargout] = odess(varargin)
%ODESS Single-Step, Variable-Order, Runge-Kutta ODE Solution.
% ODESS advances the numerical solution of a set of differential equations
% one step from t to t + h. A variable order Runge-Kutta algorithm is used,
% where a 5th order step of length h is taken. If this step fails, the same
% derivative calculations are used to take a 3rd order step of length
% 3*h/5. If this step fails, the same derivative calculations are used to
% take a 2nd order step of length h/5. If all these attempts fail, a
% revised step length h is computed and the process repeats until a step
% succeeds. This mechanism is well suited for differential equations that
% have discontinuities or sharp changes appearing at unknown times.
% ODESS advances the solution a single time step from t to t + h, thereby
% providing intimate access to the solution and ODE parameters every step.
%
% [T,Y,YP]=ODESS(ODEFUN,t,y,yp) integrates the first order differential
% equations computed in the M-file ODEFUN(t,y) from the time point t where
% y = y(t) and yp = ODEFUN(t,y). ODEFUN must be a function handle and must
% return a column vector of derivatives given the time t and vector y.
% [T,Y,YP] are the results at the time T = t + h, where Y=Y(T) is a row
% vector of the solution, and YP=ODEFUN(T,Y) is a row vector of the slopes.
%
% ODESS('PName1',PValue1,'PName2',PValue2,...) sets ODESS parameters to the
% corresponding property values. ODESS(P) sets parameters using the
% structure P having fieldnames equal to parameter names and corresponding
% contents equal to property values. ODESS('Reset') resets all parameters
% to their default values including 'NextStep' so that a new problem may be
% started.
%
% ODESS('PName') returns the value associated with the parameter 'PName'.
% ODESS({'PName1' 'PName2' ...}) returns parameter values associated with
% the listed property values in a cell array.
% P=ODESS with no input argument returns a structure P of all settable
% parameter values as described above.
%
%  Parameter    Default
%  Name         Value    Limits       Description
% 'RelTol'      1e-3      >1e-10      Relative Error Tolerance
% 'AbsTol'      1e-6      >0          Absolute Error Tolerance
% 'MinStep'     1e-10     >1e-12      Minimum Stepsize, Hmin
% 'MaxStep'       1       >=Hmin      Maximum Stepsize
% 'NextStep'     []       >=Hmin      Next Stepsize, Hnext
% 'SafetyFactor' 0.9     (0.75,.95)   Stepsize Safety Factor
% 'GrowthLimit'    5     (2,20)       Stepsize Growth Limit Ratio
% 'ShrinkLimit'  0.1     (.05,.5)     Stepsize Shrink Limit Ratio
% 'FallBack'     'on'    ('on','off') Enable Fall Back on Full Step Failure
%
% The stepsize safety factor is the normalized amount of the predicted next
% successful step size h that is attempted on the next step.
% The stepsize growth limit is the maximum normalized amount that the step
% size h is permitted to grow from one successful step to the next.
% The stepsize shrink limit is the maximum normalized amount that the step
% size h is permitted to shrink from one successful step to the next.
% If fallback is 'off', only 5th order steps are attempted.
%
% [T,Y,YP,Stats]=ODESS(ODEFUN,t,y,yp) in addition returns statistics
% Stats = [Ierr Fail Order] where Ierr identifies the variable Y(Ierr)
% which dominated the error in the step, Fail is the number of failed steps
% encountered in this integration step and Order is the order
% of the accepted solution. Order is 2, 3, or 5.  
%
% Typical Usage:
%                 t=0; % initial time
%                 y=[y1;y2;...;yn];   % initial conditions
%                 fun=@odefun;        % create function handle to ODEFUN
%                 yp=feval(fun,t,y);  % compute initial slopes
%                 odess('Reset')      % reset parameters if new problem
%                 odess('PName',Pvalue,...) % optionally set parameters
%                 while test
%                    [t,y,yp]=odess(fun,t,y,yp);
%                    % process data as needed
%                    % modify ode parameters if needed
%                    % do not change t,y,yp from call to call!
%                    % break from loop when desired
%                 end
% Reference:
% J.R. Cash and A.H. Harp, "A Variable Order Runge-Kutta Method for Initial
% Value Problems with Rapidly Varying Right Hand Sides," ACM Trans. on
% Math. Software, vol. 16, no. 3, pp. 201-222, 1990.

% D.C. Hanselman, University of Maine, Orono, ME 04469
% Mastering MATLAB 7
% masteringmatlab@yahoo.com

persistent P Pc A B o5 Eo5 o3 Eo3 o2 Eo2

if isempty(P)
   [P,Pc] = local_getdefaults;
   [A,B,o5,Eo5,o3,Eo3,o2,Eo2] = local_getparameters;
end

if nargin==0                          % return structure of property values
   varargout{1} = P;
elseif nargin==1 && ischar(varargin{1})...                 % ODESS('Reset')
       && strcmpi(varargin{1},'reset')
   [P,Pc] = local_getdefaults;
elseif nargin==1 && ischar(varargin{1})                    % ODESS('PName')
   [fn,errmsg] = local_isfield(Pc,varargin{1});
   error(errmsg)
   varargout{1} = P.(fn);
elseif nargin==1 && iscellstr(varargin{1}) % ODESS({'PName1' 'PName2' ...})
   arg = varargin{1};
   N = length(arg);
   out = cell(size(arg));
   for k = 1:N
      [fn,errmsg] = local_isfield(Pc,arg{k});
      error(errmsg)
      out{k} = P.(fn);
   end
   varargout{1} = out;
elseif nargin==1 && isstruct(varargin{1})                        % ODESS(P)
   arg = varargin{1};
   fnarg = fieldnames(arg);
   N = length(fnarg);
   for k = 1:N
      [fn,errmsg] = local_isfield(Pc,fnarg{k});
      error(errmsg)
      P.(fn) = arg.(fnarg{k});
   end
   P = local_checklimits(P);
elseif rem(nargin,2)==0 && ischar(varargin{1}) %ODESS('PName1',PValue1,...) 
   N = nargin;
   for k = 1:2:N-1
      arg = varargin{k};
      if ischar(arg)
         [fn,errmsg] = local_isfield(Pc,arg);
         error(errmsg)
         P.(fn) = varargin{k+1};
      else
         error('String Property Name Expected.')
      end
   end
   P = local_checklimits(P);
elseif nargin==4 && isa(varargin{1},'function_handle')...
                 && nargout>=3                       % ODESS(ODEFUN,t,y,yp)
   odefun = varargin{1};
   t = varargin{2};
   yi = varargin{3}(:);
   yp = varargin{4}(:);
   f = zeros(length(yi),6);
   f(:,1) = yp;  % slopes at t

   if isempty(P.NextStep) % first call, guess a stepsize
      ypm = max(abs(yp),eps);
      P.NextStep = P.SafetyFactor*...
                 min(abs((P.RelTol*abs(yi)+P.AbsTol)./ypm)).^0.2;
   end
   h = P.NextStep;  % get step size to try
   
   fail = 0; % reset number of failed tries

   while 1  % advance one step
      h  = min(max([h; P.MinStep; 16*eps*t]), P.MaxStep); % limit stepsize
       disp([h, P.MaxStep])
      hA = A*h;
      hB = B*h;
      disp(hA)
      f(:,2) = feval(odefun,t+hA(1),yi+f*hB(:,1));
      f(:,3) = feval(odefun,t+hA(2),yi+f*hB(:,2));
      f(:,4) = feval(odefun,t+hA(3),yi+f*hB(:,3));
      disp([t, t+hA(4)])
      f(:,5) = feval(odefun,t+hA(4),yi+f*hB(:,4));
      f(:,6) = feval(odefun,t+hA(5),yi+f*hB(:,5));
      tn = t+h;  % prospective solution
      yn = yi+h*f*o5;

      yerr = abs(h*f*Eo5);                                  % compute errors
      ytol = P.RelTol*max(abs(yi),abs(yn))+P.AbsTol;        % tolerances
      [err,ierr] = max(yerr-ytol);                          % worst variable
      err = max(yerr(ierr),eps);                            % worst error
      tol = ytol(ierr);                                     % worst tolerance
      hratio = P.SafetyFactor*(tol/err)^0.2;

      if err<=tol                                   % 5th order step worked
         if fail==0 % find new step size only if no failure
            P.NextStep = h*min(P.GrowthLimit,hratio);
         else       % retain last successful step size
            P.NextStep = h;
         end
         stats = [ierr fail 5];
         break
      elseif strcmp(P.FallBack,'on')    % step failed, check 3rd order step
         if h<=max(P.MinStep,16*eps*t)
            error('Hmin reached: Step Failure at t = %.4g.',t)
         end
         tn = t+3*h/5;  % prospective solution
         yn = yi+h*f*o3;

         yerr = abs(h*f*Eo3);                               % compute errors
         ytol = P.RelTol*max(abs(yi),abs(yn))+P.AbsTol;     % tolerances
         [err,ierr] = max(yerr-ytol);                       % worst variable
         err = max(yerr(ierr),eps);                         % worst error
         tol = ytol(ierr);                                  % worst tolerance
         if err<=tol  % 3rd order solution works
            P.NextStep = 3*h/5;
            stats = [ierr fail 3];
            break
         else                       % step failed, check 2nd order solution
            tn = t+h/5;  % prospective solution
            yn = yi+h*f*o2;

            yerr = abs(h*f*Eo2);                            % compute errors
            ytol = P.RelTol*max(abs(yi),abs(yn))+P.AbsTol;  % tolerances
            [err,ierr] = max(yerr-ytol);                    % worst variable
            err = max(yerr(ierr),eps);                      % worst error
            tol = ytol(ierr);                               % worst tolerance
            if err<=tol  % 2nd order solution works
               P.NextStep = 2*h/5;
               stats = [ierr fail 2];
               break
            else                   % all fallbacks fail, go back to order 5
               if fail==0 % first failure, get a new step size
                  fail = 1;
                  h = h*max(P.ShrinkLimit,hratio);
               else       % multiple failure, so just cut step in half
                  fail = fail+1;
                  h = h/2;
               end
            end
         end
      else                     % no fallback, try a new stepsize at order 5
         if fail==0 % first failure, get a new step size
            fail = 1;
            h = h*max(P.ShrinkLimit,hratio);
         else       % multiple failure, so just cut step in half
            fail = fail+1;
            h = h/2;
         end
      end
   end
   varargout = {tn,yn.',feval(odefun,tn,yn).'};
   if nargout>3
      varargout{4} = stats;
   end
else
   error('Unknown Input Syntax or Not Enough Output Arguments.')
end
%--------------------------------------------------------------------------
function [out,errmsg] = local_isfield(fnames,str)             % local_isfield
% compare str to fnames, if found, return complete fieldname
% otherwise return error and empty string.
% fnames is cell array, str is a string
% outputs are strings

idx = find(strcmpi(fnames,str)); % check for exact match first

if isempty(idx) % no exact match, so look for more general match
   idx = find(strncmpi(str,fnames,length(str)));
end
if numel(idx)==1 % unique match found
   out = fnames{idx};
   errmsg = '';
else             % trouble
   out = '';
   errmsg = sprintf('Unknown or Not Unique Property: %s',str);
end
%--------------------------------------------------------------------------
function [P,fn] = local_getdefaults                       % local_getdefaults
tmp = {'RelTol' 'AbsTol' 'MinStep' 'MaxStep' 'NextStep' ...
     'SafetyFactor' 'GrowthLimit' 'ShrinkLimit' 'FallBack'
       1e-3     1e-6      1e-10      1         nan    ...
           0.9            5              0.1        'on'};
P = struct(tmp{:});
fn = fieldnames(P);
%--------------------------------------------------------------------------
function P = local_checklimits(P)                         % local_checklimits
P.RelTol = max(P.RelTol(1),1e-10);
P.AbsTol = max(P.AbsTol(1),0);
P.MinStep = max(P.MinStep(1),1e-12);
P.MaxStep = max(P.MaxStep(1),P.MinStep);
P.NextStep = [];
P.SafetyFactor = min(0.95,max(0.75,P.SafetyFactor(1)));
P.GrowthLimit = min(20,max(2,P.GrowthLimit(1)));
P.ShrinkLimit = min(0.5,max(0.05,P.ShrinkLimit(1)));
if strncmpi(P.FallBack,'of',2)
   P.FallBack = 'off';
else % default is 'on'
   P.FallBack = 'on';
end
%--------------------------------------------------------------------------
function [A,B,o5,Eo5,o3,Eo3,o2,Eo2] = local_getparameters
A = [1/5 3/10 3/5 1 7/8];  % 1 by 5, time increments
B = [1/5        0       0         0            0        0
   3/40       9/40    0         0            0        0
   3/10       -9/10   6/5       0            0        0
   -11/54     5/2     -70/27    35/27        0        0
   1631/55296 175/512 575/13824 44275/110592 253/4096 0]';  % 5 by 6

o5 = [37/378;0;250/621;125/594;0;512/1771];  % 6 by 1
Eo5 = o5-[2825/27648;0;18575/48384;13525/55296;277/14336;1/4];
o3 = [1/10;0;2/5;1/10;0;0];
Eo3 = [1;0;-2;1;0;0];
o2 = [1;1;0;0;0;0]/10;
Eo2 = [-1;1;0;0;0;0];
%--------------------------------------------------------------------------
