
function o_color = mporange(m);

if nargin==0;
  m=64;
end

s = linspace(0,1,m);

for i=1:length(s)
    red(i) = (1000.^s(i)-1)/(1000^1-1);
    green(i) = 0.5 * (1000.^s(i)-1)/(1000^1-1);
    blue(i) = 0;
end

o_color=[red; green; blue;]';