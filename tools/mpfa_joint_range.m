
function mpfa_joint_range (x, y, t_start, t_stop)

xc_data = fopen([x, '_c.dat']);
xr_data = fopen([x, '_r.dat']);
xs_data = fopen([x, '_s.dat']);
xd_data = fopen([x, '_d.dat']);
yc_data = fopen([y, '_c.dat']);
yr_data = fopen([y, '_r.dat']);
ys_data = fopen([y, '_s.dat']);
yd_data = fopen([y, '_d.dat']);

for i = 1:(t_start - 1)
    fgetl(xc_data);
    fgetl(xr_data);
    fgetl(xs_data);
    fgetl(xd_data);
    fgetl(yc_data);
    fgetl(yr_data);
    fgetl(ys_data);
    fgetl(yd_data);
end

% maxterms = 15;
% e = zeros(maxterms, 2^maxterms);
% for i = 1:2^maxterms
%     e(:, i) = double(bitget(i - 1, 1:maxterms))';
% end
% e(e == 0) = -1;

figure;
hold on;

% Plot interval regions
for i = t_start:t_stop
    disp(num2str(i));

    xc = sscanf(fgetl(xc_data), '%f');
    xr = sscanf(fgetl(xr_data), '%f');
    yc = sscanf(fgetl(yc_data), '%f');
    yr = sscanf(fgetl(yr_data), '%f');

    pos = [(xc - xr), (yc - yr), (2 .* xr), (2 .* yr)];
    %pos = [-xr, -yr, (2 .* xr), (2 .* yr)];
    rectangle('Position', pos, 'EdgeColor', 'b');
    drawnow;
end

% % Plot affine regions
% for i = t_start:t_stop
%     disp(num2str(i));
% 
%     %xc = sscanf(fgetl(xc_data), '%f');
%     %xr = sscanf(fgetl(xr_data), '%f');
%     xs = sscanf(fgetl(xs_data), '%u');
%     xd = sscanf(fgetl(xd_data), '%f');
%     %yc = sscanf(fgetl(yc_data), '%f');
%     %yr = sscanf(fgetl(yr_data), '%f');
%     ys = sscanf(fgetl(ys_data), '%u');
%     yd = sscanf(fgetl(yd_data), '%f');
% 
%     us = union(xs, ys);
%     if isrow(us)
%         us = us';
%     end
%     terms = size(us, 1);
% 
%     ix = ismember(us, xs);
%     xxd = zeros(1, terms);
%     xxd(ix) = xd;
%     iy = ismember(us, ys);
%     yyd = zeros(1, terms);
%     yyd(iy) = yd;
% 
%     % Get all permutations of noise symbol extremities
%     xx = zeros(2^terms, 1);
%     yy = zeros(2^terms, 1);
%     for j = 1:2^terms
%         e = double(bitget(j - 1, 1:terms))';
%         e(e == 0) = -1;
%         %xx(j) = xc + xxd * e;
%         %yy(j) = yc + yyd * e;
%         xx(j) = xxd * e;
%         yy(j) = yyd * e;
%     end
% 
% %     % Get all permutations of noise symbol extremities
% %     %xx = xc + xxd * e(1:terms, 1:2^terms);
% %     %yy = yc + yyd * e(1:terms, 1:2^terms);
% %     xx = xxd * e(1:terms, 1:2^terms);
% %     yy = yyd * e(1:terms, 1:2^terms);
% 
%     k = convhull(xx, yy);
%     plot(xx(k), yy(k));
%     drawnow;
% end

xlabel(x); ylabel(y);
hold off;

fclose(xc_data);
fclose(xr_data);
fclose(xs_data);
fclose(xd_data);
fclose(yc_data);
fclose(yr_data);
fclose(ys_data);
fclose(yd_data);

end
