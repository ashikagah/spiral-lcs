
function o_color = mppink(m);

if nargin==0;
  m=64;
end

s = linspace(0,1,m);

for i=1:length(s)
    red(i) = (1000.^s(i)-1)/(1000^1-1);
    green(i) = 0.78 * (1000.^s(i)-1)/(1000^1-1);
    blue(i) = 0.80 * (1000.^s(i)-1)/(1000^1-1);
end

o_color=[red; green; blue;]';