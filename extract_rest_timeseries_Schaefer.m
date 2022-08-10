function extract_rest_timeseries_Schaefer(resolution, subj_ls, outdir)
    % extract_rest_timeseries_Schaefer(resolution, subj_ls, outdir)
    % Given a resolution of the Schaefer parcellation (e.g. 200), extract and average the resting-state
    % timeseries within each parcel for every run (fsLR32k space), save into `outdir` in BIDS format.

    start_dir = pwd;
    proj_dir = '/data/project/parcellate_ABCD_preprocessed';
    addpath(genpath(fullfile(proj_dir, 'scripts', 'external', 'cifti-matlab')))

    if(~ischar(resolution))
        resolution = num2str(resolution);
    end

    indir = fullfile(proj_dir, 'data', 'inm7-superds', 'original', 'abcd', 'derivatives', 'abcd-hcp-pipeline');
    if(~exist('outdir', 'var'))
        outdir = fullfile(proj_dir, 'data', 'parcellated_timeseries');
    end

    % read Schaefer parcellation
    parcellation = fullfile(proj_dir, 'data', 'SchaeferParcellations', 'HCP', 'fslr32k', 'cifti', ...
        ['Schaefer2018_' resolution 'Parcels_17Networks_order.dlabel.nii']);
    x = ft_read_cifti(parcellation, 'mapname', 'array');

    % loop through subjects
    subjects = text2cell(subj_ls);
    ses = 'ses-baselineYear1Arm1';
    for i = 1:length(subjects)
        s = subjects{i};
        cd(indir)
        system(sprintf('datalad get -n %s', s));
        system(sprintf('git -C %s config --local --add remote.datalad.annex-ignore true', s));

        cd(fullfile(indir, s, ses, 'func'))
        [flag, msg] = system(sprintf('ls -d %s', [s '_' ses '_task-rest_run-*_bold_timeseries.dtseries.nii']))
        runs = strsplit(msg);
        for j = 1:length(runs)
            if(~isempty(runs{j}))
                runnum = strsplit(runs{j}, '_'); runnum = runnum{4};
                system(sprintf('datalad get -s inm7-storage %s', runs{j}));

                ts = ft_read_cifti(runs{j});
                pts = zeros(str2num(resolution), size(ts.dtseries, 2));
                for roi = 1:str2num(resolution)
                    pts(roi,:) = mean(ts.dtseries(x.dlabel==roi, :), 1);
                end
                
                mkdir(fullfile(outdir, s, ses, 'func'))
                save(fullfile(outdir, s, ses, 'func', [s '_' ses '_task-rest_' runnum ...
                    '_bold_atlas-Schaefer' resolution '_timeseries']), 'pts', '-v7.3')

                system(sprintf('datalad drop %s', runs{j}));
            end
        end

        cd(indir)
        system(sprintf('datalad uninstall %s', s));
        
    end

    rmpath(genpath(fullfile(proj_dir, 'scripts', 'external', 'cifti-matlab')))
    cd(start_dir)
end

function cell_array = text2cell(text_file)
    num_lines = 0;
    fid = fopen(text_file);
    while (~feof(fid))
        num_lines = num_lines + 1;
        cell_array{num_lines} = fgetl(fid);
    end
    fclose(fid);

end