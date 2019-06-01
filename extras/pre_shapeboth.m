% Find the Pre_shape for any arbitary shape
function [X0]=pre_shapeboth(X)
% s is the shape sample in the n*k, n points and each point lies in R^k


[n, m]   = size(X);
% [ny, my] = size(Y);

if ny ~= n
    error(message('stats:procrustes:InputSizeMismatch'));
% elseif my > m
%     error(message('stats:procrustes:TooManyColumns'));
% end

% Center at the origin.
muX = mean(X,1);
% muY = mean(Y,1);
X0 = X - repmat(muX, n, 1);
% Y0 = Y - repmat(muY, n, 1);

ssqX = sum(X0.^2,1);
% ssqY = sum(Y0.^2,1);
constX = all(ssqX <= abs(eps(class(X))*n*muX).^2);
% constY = all(ssqY <= abs(eps(class(X))*n*muY).^2);
ssqX = sum(ssqX);
% ssqY = sum(ssqY);

if ~constX %&& ~constY
    % The "centered" Frobenius norm.
    normX = sqrt(ssqX); % == sqrt(trace(X0*X0'))
%     normY = sqrt(ssqY); % == sqrt(trace(Y0*Y0'))

    % Scale to equal (unit) norm.
    X0 = X0 / normX;
%     Y0 = Y0 / normY;

    % Make sure they're in the same dimension space.
%     if my < m
%         Y0 = [Y0 zeros(n, m-my)];
%     end


end
