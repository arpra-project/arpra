
function mpfa_joint_range (x, y, t)

    xc_data = fopen([x, '_c.dat']);
    xn_data = fopen([x, '_n.dat']);
    xs_data = fopen([x, '_s.dat']);
    xd_data = fopen([x, '_d.dat']);

    yc_data = fopen([y, '_c.dat']);
    yn_data = fopen([y, '_n.dat']);
    ys_data = fopen([y, '_s.dat']);
    yd_data = fopen([y, '_d.dat']);

    maxterms = 10;
    e = zeros(maxterms, 2^maxterms);
    for j = 1:2^maxterms
        e(:, j) = bitget(j - 1, 1:maxterms);
    end
    e(e == 0) = -1;

    figure;
    hold on;

    for i = 1:t
        disp(num2str(i));

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
        end

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
        xxd = sum(xxd .* e(1:terms, 1:2^terms));
        yyd = sum(yyd .* e(1:terms, 1:2^terms));
        %xxd = xc + sum(xxd .* e(1:terms, 1:2^terms));
        %yyd = yc + sum(yyd .* e(1:terms, 1:2^terms));

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
end
