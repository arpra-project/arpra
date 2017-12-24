
function arpra_affine_approx(f)

% Chebyshev approximation of f in [a, b]
disp(f);

a = sym('a');
b = sym('b');
u = sym('u');
alpha = sym('aplha');
gamma = sym('gamma');
delta = sym('delta');

assume(a < b);

% handle singularities


alpha = (f(b) - f(a)) / (b - a);
pretty(simplify(alpha, 'Steps', 100));
u = solve(diff(f(u)) == alpha, u, 'PrincipalValue', true);
pretty(simplify(u, 'Steps', 100));

da = f(a) - alpha * a;
pretty(simplify(da, 'Steps', 100));
db = f(b) - alpha * b;
pretty(simplify(db, 'Steps', 100));
du = f(u) - alpha * u;
pretty(simplify(du, 'Steps', 100));


%dlo = min([da, db, du])
%dhi = max([da, db, du])

% 
% gamma = (d2 + du) ./ 2;
% r1 = d2 - gamma;
% r2 = gamma - du;
% delta = max(r1, r2);
% 
% 
% figure();
% hold on;
% %grid on;
% axis([-10, 10, -10, 10]);
% 
% line([xl0, xl0], [-10, 10], 'color', 'r');
% line([xh0, xh0], [-10, 10], 'color', 'r');
% %line([xc0, xc0], [-10, 10], 'color', 'r');
% 
% plot(x, y, 'b');
% plot(x, ax, 'r');
% plot(x, ax + gamma, 'g');
% plot(x, ax + gamma + delta, 'g');
% plot(x, ax + gamma - delta, 'g');
% 
% hold off;
