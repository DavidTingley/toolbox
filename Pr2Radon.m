function [slope integral] = Pr2Radon(Pr)
% computes max projection line, using radon transform, on Pr matrix

plotting = 0;

Pr(isnan(Pr)) = 0;
theta = 0:.1:180;
[R xp] = radon(Pr,theta);
y(1) = round(size(Pr,1)./2);
x(1) = round(size(Pr,2)./2);


[Y,I] = max(R);
[a b]  = max(Y);
angle = theta(b);
offset = x(1) - (I(b));



y(2) = y(1) + offset*sin(deg2rad(-angle));
x(2) = x(1) + offset*cos(deg2rad(-angle));
coeffs = polyfit(x, y, 1);
xx = 1:size(Pr,2);
yy = (-1/coeffs(1))*(xx - x(1)) + y(1);
if plotting
    imagesc(Pr)
    hold on
    line([x(1) x(2)],[y(1) y(2)])
    plot(xx,yy)
end
coeffs = polyfit(xx, yy, 1);
slope = coeffs(1);

if abs(slope) < size(Pr,1) ./ size(Pr,2) % rise/run limit to calc integral (must be in the frame)
    for i=1:length(xx)
        curve(i) = Pr(round(yy(i)),xx(i));
    end
    integral = nansum(curve);
else
    integral = NaN;
    slope = NaN;
end
integral = double(integral); % weird typecasting fix

