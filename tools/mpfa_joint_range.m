
function mpfa_joint_range (x, y, t)

    %xc_data = fopen([x, '_c.dat']);
    xs_data = fopen([x, '_s.dat']);
    xd_data = fopen([x, '_d.dat']);
    %yc_data = fopen([y, '_c.dat']);
    ys_data = fopen([y, '_s.dat']);
    yd_data = fopen([y, '_d.dat']);

    for i = 1:(t(1) - 1)
        %sscanf(fgetl(xc_data), '%f');
        sscanf(fgetl(xs_data), '%u');
        sscanf(fgetl(xd_data), '%f');
        %sscanf(fgetl(yc_data), '%f');
        sscanf(fgetl(ys_data), '%u');
        sscanf(fgetl(yd_data), '%f');
    end

%    maxterms = 22;
%    e = zeros(maxterms, 2^maxterms);
%    for j = 1:2^maxterms
%        e(:, j) = bitget(j - 1, 1:maxterms);
%    end
%    e(e == 0) = -1;

    figure;
    hold on;

    for i = t
        disp(num2str(i));

        %xc = sscanf(fgetl(xc_data), '%f');
        xs = sscanf(fgetl(xs_data), '%u');
        xd = sscanf(fgetl(xd_data), '%f');
        %yc = sscanf(fgetl(yc_data), '%f');
        ys = sscanf(fgetl(ys_data), '%u');
        yd = sscanf(fgetl(yd_data), '%f');

        us = union(xs, ys);
        if isrow(us)
            us = us';
        end

        ix = ismember(us, xs);
        xxd = zeros(size(us));
        xxd(ix) = xd;
        iy = ismember(us, ys);
        yyd = zeros(size(us));
        yyd(iy) = yd;

        % Get all permutations of noise symbol extremities
        terms = size(us, 1);
        %xx = ones(2^terms, 1) .* xc;
        %yy = ones(2^terms, 1) .* yc;
        xx = zeros(2^terms, 1);
        yy = zeros(2^terms, 1);

        for j = 1:2^terms
            e = double(bitget(j - 1, 1:terms))';
            e(e == 0) = -1;
            xx(j) = sum(xxd .* e);
            yy(j) = sum(yyd .* e);
        end

%        % Get all permutations of noise symbol extremities
%        terms = size(us, 1);
%        xxd = repmat(xxd, 1, 2^terms);
%        yyd = repmat(yyd, 1, 2^terms);
%        %xxd = xc + sum(xxd .* e(1:terms, 1:2^terms));
%        %yyd = yc + sum(yyd .* e(1:terms, 1:2^terms));
%        xxd = sum(xxd .* e(1:terms, 1:2^terms));
%        yyd = sum(yyd .* e(1:terms, 1:2^terms));

        k = convhull(xx, yy);
        plot(xx(k), yy(k));
        drawnow;
    end

    xlabel(x); ylabel(y);
    hold off;

    %fclose(xc_data);
    fclose(xs_data);
    fclose(xd_data);
    %fclose(yc_data);
    fclose(ys_data);
    fclose(yd_data);
end
