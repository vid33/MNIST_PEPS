function [R Q]=rq_m_greater_n(A, varargin)
%For now, use this when m>n

% function [R Q]=rq(A)
%
% Perform RQ factorization
% A [m x n]
% returns R [m x n] triangular superior and
% Q [n x n] unitary (i.e., Q'*Q=I)
% such that A = R*Q
%
% Note: Standard RQ requires m<=n.
%
% For n ~=m, default output R has non-zero elements populate
% at the m-first rows
% Use [R Q] = rq(A, 'last') to populate non-sero R-rows at other end.
%
% For "non standard" RQ, i.e. for A with m > n, zero are padded in the
% outputs R and Q. However Q is no longer unitary. More precisely
% Q'*Q = eye(n) but Q*Q' is
% diag([zeros(1,m-n) ones(1,n)]) without 'last' option, and
% diag([ones(1,n) zeros(1,m-n)]) with 'last' and R becomes triangular.
%
% Author: Bruno Luong
% last Update: 05/Oct/2008

[m n]=size(A);

[Q R]=qr(flipud(A).');
R=flipud(R.'); % m x n
Q=Q.'; % n x n
if m>n
   % warning('RQ:DimensionBizarre',...
   %         'RQ: number of rows is larger the number of columns'); 
    % Padding zeros
    R(end,m)=0; % m x m
    Q(m,end)=0; % m x n
end

R(:,1:m)=R(:,m:-1:1);
Q(1:m,:)=Q(m:-1:1,:);

last=strcmpi(getoption('',varargin{:}),'last');
% Populate R at last rows
if xor(m>n,last)
    R=circshift(R,[0 n-m]);
    Q=circshift(Q,[n-m 0]);
end

%my addition - get rid of the zero padding
if(m>n)
    R = R(1:m, 1:n);
    Q = Q(1:n, 1:n);
    
end

end

% Get option if provided
function res=getoption(default, option)
  if nargin<2 || isempty(option)
      res=default;
  else
      res=option;
  end
end

