function results = export_refinement_sweep_final(N_h_values)

  if (nargin < 1) || isempty(N_h_values)
    N_h_values = [3, 4];
  end

  output_root = fullfile(fileparts(pwd), 'replication_output_3d');
  image_root = fullfile(pwd, 'images');
  results = struct([]);

  global N_h_override

  for i = 1:numel(N_h_values)
    N_h_override = N_h_values(i);
    outdir = fullfile(output_root, sprintf('full_Nh%d', N_h_override));
    image_dir = fullfile(image_root, sprintf('full_Nh%d', N_h_override));

    export_loading_process_final(outdir, false);
    plot_full_solution_final(outdir, image_dir);

    meta = read_metadata(fullfile(outdir, 'meta.txt'));
    results(i).N_h = N_h_override;
    results(i).nt = meta.nt;
    results(i).nv = meta.nv;
    results(i).final_zeta = meta.final_zeta;
    results(i).final_max_displacement = meta.final_max_displacement;
    results(i).plastic_points = meta.plastic_points;
  end

  clear global N_h_override

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
