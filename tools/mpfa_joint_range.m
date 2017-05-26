
x0 = 20;
x = [5; 0; -5; 0];
%x = [-4; 0; 2; 3];
%x = [-4; 1; 2; 3];

y0 = 10;
y = [5; 0; 5; 0];
%y = [-2; 1; 0; -1];
%y = [-2; 0; 0; 1];


u0 = 0;
u = [6; 0; 4; 0];

v0 = 0;
v = [0; 3; 0; 7];

x0 = x0 + u0 + v0;
x = x + u + v;

y0 = y0 + u0 - v0;
y = y + u - v;


% Get all permutations of noise symbol extremities
terms = size(x, 1);
x = repmat(x, 1, 2^terms);
y = repmat(y, 1, 2^terms);
e = zeros(terms, 2^terms);

for i = 1:2^terms
    e(:, i) = bitget(i - 1, 1:terms);
end
e(e == 0) = -1;

xx = x0 + sum(x .* e);
yy = y0 + sum(y .* e);

k = convhull(xx' ,yy');
plot(xx(k), yy(k));
