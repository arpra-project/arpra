
xc_data = fopen('../out_V_centre');
xn_data = fopen('../out_V_nterms');
xs_data = fopen('../out_V_symbols');
xd_data = fopen('../out_V_deviations');

yc_data = fopen('../out_N_centre');
yn_data = fopen('../out_N_nterms');
ys_data = fopen('../out_N_symbols');
yd_data = fopen('../out_N_deviations');

%figure;
figure(2); clf;
%hold on;
%axis([-20, 20, 0, 1]);

t = 20;

for i = 1:t
    %sprintf('%.80f', xc)

    xc = sscanf(fgetl(xc_data), '%f');
    xn = sscanf(fgetl(xn_data), '%u');
    xs = sscanf(fgetl(xs_data), '%u');
    xd = sscanf(fgetl(xd_data), '%f');

    yc = sscanf(fgetl(yc_data), '%f');
    yn = sscanf(fgetl(yn_data), '%u');
    ys = sscanf(fgetl(ys_data), '%u');
    yd = sscanf(fgetl(yd_data), '%f');

    us = union(xs, ys);
    if isrow(us)
        us = us';
    end;
    
    ix = find(ismember(us, xs));
    xxd = zeros(size(us));
    xxd(ix) = xd;

    iy = find(ismember(us, ys));
    yyd = zeros(size(us));
    yyd(iy) = yd;
    
    % Get all permutations of noise symbol extremities
    terms = size(us, 1);
    xxd = repmat(xxd, 1, 2^terms);
    yyd = repmat(yyd, 1, 2^terms);
    e = zeros(terms, 2^terms);
    for j = 1:2^terms
        e(:, j) = bitget(j - 1, 1:terms);
    end
    e(e == 0) = -1;

    xxd = xc + sum(xxd .* e);
    yyd = yc + sum(yyd .* e);

    k = convhull(xxd, yyd);
    plot(xxd(k), yyd(k));
    drawnow;
end

hold off;

fclose(xc_data);
fclose(xn_data);
fclose(xs_data);
fclose(xd_data);

fclose(yc_data);
fclose(yn_data);
fclose(ys_data);
fclose(yd_data);
