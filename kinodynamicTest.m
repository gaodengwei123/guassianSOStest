%% ================Kinodynamic method========================start
% [A,B] = tv_poly_linearize(f0, @(t) xT, @(t) uT);
%
% % [K,S] = lqr(A(T)+eye(4),B(T),Q(T),R(T));
%  % check the controlablity of system
% if rank(ctrb(A(T),B(T))) ~= size(A(T),1)
%     disp('System is not controllable - aborting');
%     return;
% end
%
% state_dims = size(A(T),1);
% % set the state variable
% syms t x t_s
% obj.x0 = sym('x0',[state_dims,1]);
% obj.x0 = sym(obj.x0, 'real');
% obj.x1 = sym('x1',[state_dims,1]);
% obj.x1 = sym(obj.x1, 'real');
%
% x_bar = expm(A(T)*t)*obj.x0;
% G = int(expm(A(T)*(t-x))*B(T)/R(T)*B(T)'*expm(A(T)'*(t-x)), x, 0, t);
%
% d = G\(obj.x1-x_bar);
% tau_star = 1-2*(A(T)*obj.x1)'*d-d'*B(T)/R(T)*B(T)'*d;
%
% solution = expm([A(T), B(T)/R(T)*B(T)';zeros(state_dims), -A(T)']*(t-t_s))*[obj.x1;subs(d,t,t_s)];
%
% control = R(T)\B(T)'*solution(state_dims+1:2*state_dims,:);
% states = solution(1:state_dims);
%
% % cost = int(x'*Q*x+control'*R*control, t, 0, t_s);
%
% %% calculate the factors for the resulting polynomial explicitly
% p1 = [];
% it = 0;
% while it < 20 && length(p1) <= 1 %just so mupad knows it's a polynomial
%     p1 = feval(symengine,'coeff',simplifyFraction(tau_star*t^it), t, 'All');
%     it = it+1;
% end
% if it > 20
%     disp('either the result is not a polynomial or the degree is too high');
% end
%
% p([obj.x0', obj.x1']) = fliplr(p1);
% end_time = matlabFunction(p);
% eval_states_internal = matlabFunction(states);
% eval_control_internal = matlabFunction(control);
%
% in = num2cell([x0', xT']);
% time = roots(end_time(in{:}));
%
% time = time(imag(time)==0);
% time = min(time(time>=0));
%
% inj = num2cell([time, x0', xT']);
% states = @(t)eval_states_internal(t, inj{:});
% inputs = @(t)eval_control_internal(t, inj{:});
%
% pd = states(1:2);
% dpd = states(4:5);
% psi = states(3);
% dpsi = states(6);
%% ================Kinodynamic method==========================end