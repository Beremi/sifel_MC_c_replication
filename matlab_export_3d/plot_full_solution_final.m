function plot_full_solution_final(root, image_dir)

  if (nargin < 1) || isempty(root)
    root = fullfile(fileparts(pwd), 'replication_output_3d', 'full');
  end
  if (nargin < 2) || isempty(image_dir)
    image_dir = fullfile(pwd, 'images');
  end

  ensure_directory(image_dir);

  coord = readmatrix(fullfile(root, 'mesh', 'coord.txt'));
  elem = round(readmatrix(fullfile(root, 'mesh', 'elem.txt')));
  surf = round(readmatrix(fullfile(root, 'mesh', 'surf.txt')));
  U = readmatrix(fullfile(root, 'matlab', 'U_final.txt'));
  E_sifel = readmatrix(fullfile(root, 'matlab', 'E_final.txt'));
  return_type = readmatrix(fullfile(root, 'matlab', 'return_type_matlab.txt'));
  plastic_norm = readmatrix(fullfile(root, 'matlab', 'plastic_strain_norm.txt'));
  zeta_hist = readmatrix(fullfile(root, 'matlab', 'zeta_hist.txt'));
  settle_hist = readmatrix(fullfile(root, 'matlab', 'settle_hist.txt'));
  max_displacement_hist = readmatrix(fullfile(root, 'matlab', 'max_displacement_hist.txt'));
  plastic_count_hist = readmatrix(fullfile(root, 'matlab', 'plastic_count_hist.txt'));
  meta = read_metadata(fullfile(root, 'meta.txt'));

  return_type = return_type(:)';
  plastic_norm = plastic_norm(:)';
  zeta_hist = zeta_hist(:)';
  settle_hist = settle_hist(:)';
  max_displacement_hist = max_displacement_hist(:)';
  plastic_count_hist = plastic_count_hist(:)';

  centroids = element_centroids(coord, elem);
  displacement_norm = sqrt(sum(U.^2, 1));
  deviatoric_norm = deviatoric_strain_norm(E_sifel);

  plot_mesh(coord, surf, meta, ...
            fullfile(image_dir, 'full_solution_mesh.png'));
  plot_displacement_norm(coord, surf, displacement_norm, meta, ...
                         fullfile(image_dir, 'full_solution_displacement_norm.png'));
  plot_deviatoric_strain_norm(coord, elem, surf, deviatoric_norm, meta, ...
                              fullfile(image_dir, 'full_solution_deviatoric_strain_norm.png'));
  plot_load_history(zeta_hist, max_displacement_hist, plastic_count_hist, meta, ...
                    fullfile(image_dir, 'full_solution_load_history.png'));
  plot_return_type_surface(coord, surf, centroids, return_type, meta, ...
                           fullfile(image_dir, 'full_solution_return_type.png'));
  plot_plastic_strain_norm(coord, surf, centroids, plastic_norm, return_type, meta, ...
                           fullfile(image_dir, 'full_solution_plastic_strain_norm.png'));

end

function plot_mesh(coord, surf, meta, path)

  plot_coord = coord([1, 3, 2], :);

  fig = new_figure();
  patch('Faces', surf', 'Vertices', plot_coord', ...
        'FaceColor', [0.79, 0.82, 0.85], 'FaceAlpha', 1.0, ...
        'EdgeColor', [0.18, 0.22, 0.27], 'LineWidth', 0.50);
  hold on;
  scatter3(plot_coord(1,:), plot_coord(2,:), plot_coord(3,:), 20, [0.94, 0.94, 0.94], ...
           'filled', 'MarkerEdgeColor', [0.10, 0.12, 0.15], 'LineWidth', 0.45);

  title(sprintf('Undeformed 3D validation mesh: %.0f nodes, %.0f tetrahedra', ...
                meta.nv, meta.nelem));
  subtitle('P1 tetrahedral slope mesh, displayed with y/z axes switched');
  format_3d_axes('x', 'z', 'y');
  save_png(fig, path);

end

function plot_displacement_norm(coord, surf, displacement_norm, meta, path)

  plot_coord = coord([1, 3, 2], :);

  fig = new_figure();
  patch('Faces', surf', 'Vertices', plot_coord', ...
        'FaceVertexCData', displacement_norm(:), ...
        'FaceColor', 'interp', 'FaceAlpha', 1.0, ...
        'EdgeColor', [0.18, 0.22, 0.27], 'LineWidth', 0.35);
  hold on;
  colormap(parula(256));
  cb = colorbar;
  cb.Label.String = '||u||';
  title('3D displacement norm on undeformed mesh');
  subtitle(sprintf('max ||u|| %.6g, final zeta %.3g', max(displacement_norm), meta.final_zeta));
  format_3d_axes('x', 'z', 'y');
  save_png(fig, path);

end

function plot_deviatoric_strain_norm(coord, elem, surf, deviatoric_norm, meta, path)

  plot_coord = coord([1, 3, 2], :);
  nodal_dev = element_values_to_nodes(elem, size(coord, 2), deviatoric_norm);

  fig = new_figure();
  patch('Faces', surf', 'Vertices', plot_coord', ...
        'FaceVertexCData', nodal_dev(:), ...
        'FaceColor', 'interp', 'FaceAlpha', 1.0, ...
        'EdgeColor', [0.18, 0.22, 0.27], 'LineWidth', 0.35);
  hold on;
  colormap(parula(256));
  cb = colorbar;
  cb.Label.String = '||dev(eps)||';
  title('3D deviatoric strain norm on undeformed mesh');
  subtitle(sprintf('max ||dev(eps)|| %.6g, final zeta %.3g', max(deviatoric_norm), meta.final_zeta));
  format_3d_axes('x', 'z', 'y');
  save_png(fig, path);

end

function nodal_values = element_values_to_nodes(elem, nv, elem_values)

  nodal_values = zeros(1, nv);
  counts = zeros(1, nv);

  for e = 1:size(elem, 2)
    nodes = elem(:, e);
    nodal_values(nodes) = nodal_values(nodes) + elem_values(e);
    counts(nodes) = counts(nodes) + 1;
  end

  active = counts > 0;
  nodal_values(active) = nodal_values(active) ./ counts(active);

end

function dev_norm = deviatoric_strain_norm(E_sifel)

  exx = E_sifel(1, :);
  eyy = E_sifel(2, :);
  ezz = E_sifel(3, :);
  eyz = 0.5 * E_sifel(4, :);
  exz = 0.5 * E_sifel(5, :);
  exy = 0.5 * E_sifel(6, :);

  mean_normal = (exx + eyy + ezz) / 3;
  dxx = exx - mean_normal;
  dyy = eyy - mean_normal;
  dzz = ezz - mean_normal;

  dev_norm = sqrt(dxx.^2 + dyy.^2 + dzz.^2 + 2*(exy.^2 + exz.^2 + eyz.^2));

end

function centroids = element_centroids(coord, elem)

  nt = size(elem, 2);
  centroids = zeros(3, nt);
  for i = 1:nt
    centroids(:, i) = mean(coord(:, elem(:, i)), 2);
  end

end

function plot_load_history(zeta_hist, max_displacement_hist, plastic_count_hist, meta, path)

  fig = new_figure();
  yyaxis left
  plot(zeta_hist, max_displacement_hist, '-o', 'LineWidth', 1.8, 'MarkerSize', 4);
  ylabel('max ||u||');
  if isfield(meta, 'target_displacement')
    yline(meta.target_displacement, '--', 'target max ||u||', 'LineWidth', 1.1);
  end
  yyaxis right
  stairs(zeta_hist, plastic_count_hist, 'LineWidth', 1.8);
  ylabel('plastic integration points');
  xlabel('load factor zeta');
  title(sprintf('3D full-path load history: final zeta %.3g, plastic points %.0f', ...
                meta.final_zeta, meta.plastic_points));
  grid on;
  save_png(fig, path);

end

function plot_return_type_surface(coord, surf, centroids, return_type, meta, path)

  fig = new_figure();
  plot_coord = coord([1, 3, 2], :);
  plot_centroids = centroids([1, 3, 2], :);

  patch('Faces', surf', 'Vertices', plot_coord', ...
        'FaceColor', [0.84, 0.86, 0.88], 'FaceAlpha', 0.28, ...
        'EdgeColor', [0.55, 0.58, 0.62], 'LineWidth', 0.35);
  hold on;

  elastic = return_type == 0;
  plastic = return_type > 0;
  scatter3(plot_centroids(1, elastic), plot_centroids(2, elastic), plot_centroids(3, elastic), ...
           18, [0.15, 0.25, 0.35], 'filled', 'MarkerFaceAlpha', 0.22);
  scatter3(plot_centroids(1, plastic), plot_centroids(2, plastic), plot_centroids(3, plastic), ...
           90, return_type(plastic), 'filled', 'MarkerEdgeColor', [0.05, 0.05, 0.05], ...
           'LineWidth', 1.1);

  colormap(return_type_colormap());
  caxis([-0.5, 4.5]);
  cb = colorbar;
  cb.Ticks = 0:4;
  cb.TickLabels = {'elastic', 'smooth', 'left', 'right', 'apex'};
  title('3D return type on undeformed mesh');
  subtitle(sprintf('final zeta %.3g, settlement %.4g, smooth plastic points %.0f', ...
                   meta.final_zeta, meta.final_settlement, meta.smooth_points));
  format_3d_axes('x', 'z', 'y');
  save_png(fig, path);

end

function plot_plastic_strain_norm(coord, surf, centroids, plastic_norm, return_type, meta, path)

  fig = new_figure();
  plot_coord = coord([1, 3, 2], :);
  plot_centroids = centroids([1, 3, 2], :);

  patch('Faces', surf', 'Vertices', plot_coord', ...
        'FaceColor', [0.88, 0.88, 0.86], 'FaceAlpha', 0.24, ...
        'EdgeColor', [0.58, 0.58, 0.55], 'LineWidth', 0.35);
  hold on;

  plastic = return_type > 0;
  scatter3(plot_centroids(1, ~plastic), plot_centroids(2, ~plastic), plot_centroids(3, ~plastic), ...
           16, [0.50, 0.54, 0.58], 'filled', 'MarkerFaceAlpha', 0.18);
  scatter3(plot_centroids(1, plastic), plot_centroids(2, plastic), plot_centroids(3, plastic), ...
           120, plastic_norm(plastic), 'filled', ...
           'MarkerEdgeColor', [0.05, 0.05, 0.05], 'LineWidth', 1.2);

  colormap(parula(256));
  cb = colorbar;
  cb.Label.String = '||eps_p||';
  title('3D plastic strain norm on undeformed mesh');
  subtitle(sprintf('max ||eps_p|| %.6g at final zeta %.3g', max(plastic_norm), meta.final_zeta));
  format_3d_axes('x', 'z', 'y');
  save_png(fig, path);

end

function cmap = return_type_colormap()
  cmap = [
    0.60, 0.64, 0.68
    0.84, 0.20, 0.18
    0.18, 0.48, 0.76
    0.96, 0.64, 0.14
    0.45, 0.22, 0.62
  ];
end

function fig = new_figure()
  fig = figure('Visible', 'off', 'Color', 'w', ...
               'Renderer', 'painters', ...
               'Position', [100, 100, 1100, 760]);
end

function format_3d_axes(xlabel_text, ylabel_text, zlabel_text)
  axis equal;
  axis tight;
  grid on;
  box on;
  xlabel(xlabel_text);
  ylabel(ylabel_text);
  zlabel(zlabel_text);
  view(42, 24);
end

function save_png(fig, path)
  set(fig, 'PaperPositionMode', 'auto');
  print(fig, path, '-dpng', '-r180');
  close(fig);
  drawnow;
end

function ensure_directory(path)
  if exist(path, 'dir') ~= 7
    mkdir(path);
  end
end

function meta = read_metadata(path)

  fid = fopen(path, 'r');
  if fid < 0
    error('Cannot open metadata file: %s', path);
  end
  cleanup = onCleanup(@() fclose(fid));

  meta = struct();
  while true
    line = fgetl(fid);
    if ~ischar(line)
      break
    end
    parts = textscan(line, '%s %f');
    if ~isempty(parts{1})
      meta.(parts{1}{1}) = parts{2};
    end
  end

end
