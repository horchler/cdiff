function [d,err]=cdiff(f,x,n)
%CDIFF  First- and second-order complex step derivative approximation.
%   CDIFF(F,X) returns the first order derivative approximation of the function
%   F evaluated at X. F is a function handle with a single input argument that
%   returns an output of the same dimensions as the floating-point array X.
%
%   CDIFF(F,X,N) specifies the order of complex step derivative approximation.
%   The parameter N is a scalar integer, either 0, 1, or 2. If N is 0, the
%   function F evaluated at X is returned.
%
%   Note:
%       If N is 2, specifying a second derivative, the approximation used only
%       achieves a numerical accuracy on the order of SQRT(EPS(CLASS(X))),
%       rather than an accuracy of about EPS(CLASS(X)) for first derivatives.
%
%   See also: DIFF, EPS

%   Inspired by:
%   http://blogs.mathworks.com/cleve/2013/10/14/complex-step-differentiation/

%   Andrew D. Horchler, horchler @ gmail . com, Created 7-7-13
%   Revision: 1.0, 10-30-13


if ~isa(f,'function_handle') || nargin(f) ~= 1 || nargout(f) ~= -1
    error('cdiff:NotFunctionHandle',...
         ['The function to be integrated must be a function handle with '...
          'one input and one output.']);
end
if isempty(x) || ~isfloat(x)
    error('cdiff:InvalidX','X must be a non-empty floating-point array.');
end
if nargin < 3
    n = 1;
elseif ~isscalar(n) || ~isnumeric(n) || ~any(n == [0 1 2])
    error('cdiff:InvalidOrder',...
          'The optional derivative order must be a scalar integer: 0, 1 or 2.');
end

if n == 1
    h = 2^-28;
    d = imag(f(x+1i*h))/h;
elseif n == 2
    fx = f(x);
    if isa(fx,'single') || isa(x,'single')
        h = 2^-7;
    else
        h = 2^-13;
    end
    d = (fx-real(f(x+1i*h)))*2/h^2;
else
    d = f(x);
end

if nargout > 1
    if n == 0
        err = cast(0,class(d));
    else
        err = cast(abs(d-double(subs(diff(sym(f),n),x))),class(d));
    end
end