function points = fibonacci_sphere(N)

    points = zeros(N, 3); % Preallocate points array
    phi = pi * (sqrt(5) - 1); % Golden angle in radians

    for i = 0:(N - 1)
        z = 1 - (i / (N - 1)) * 2; % z goes from 1 to -1
        radius = sqrt(1 - z^2); % Radius at z

        theta = phi * i; % Golden angle increment

        x = cos(theta) * radius;
        y = sin(theta) * radius;

        points(i + 1, :) = [x, y, z]; % Store the point
    end
end 