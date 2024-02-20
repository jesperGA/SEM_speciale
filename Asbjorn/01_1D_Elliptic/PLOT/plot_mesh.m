function plot_mesh(mesh)
    % plot_mesh: Plots the mesh structure.

    % Extract element and node data
    X = mesh.X;
    IX = mesh.IX;

    % Plot each element in the mesh
    for e = 1:size(IX, 1)
        edof = [(e - 1) * mesh.N + 1, e * mesh.N + 1];
        x = [X([edof], 2); X(flip(edof), 2)];
        y = [0; 0; X(end, 2) / 4; X(end, 2) / 4];

        % Plot element with color based on material property
        patch(x, y, 1 - mesh.mu(e) .* [1, 1, 1], 'EdgeColor', 'none')
    end

    % Plot the bounding box
    v = [0, 0; X(end, 2), 0; X(end, 2), X(end, 2) / 4; 0, X(end, 2) / 4];
    f = [1, 2, 3, 4];
    patch('Faces', f, 'Vertices', v, 'EdgeColor', 'k', 'FaceColor', 'none', 'LineWidth', 0.4);

    % Customize plot appearance
    set(gca, 'ytick', [])
    xlim([0 - X(end, 2) * 0.01, X(end, 2) * 1.01])
    xlabel('x coordinate')
    grid on

    enhance_plot(0, 0, 0, 0, 0);
end
