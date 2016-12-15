function [ y ] = sinc( x )
%SINC cardinal sine function as defined for digital signal processing
%   sinc(x) = sin(pi * x)/pi * x
%   sinc(0) = 1

if x == 0
    y = 1;
    return;
end

y = sin(pi .* x)/pi.*x;

end

