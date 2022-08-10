function extract_rest_timeseries_Tian(scale, subj_ls, outdir)
    % extract_rest_timeseries_Tian(scale, subj_ls, outdir)
    % Given a scale of Tian's subcortical parcellation (e.g. 2, from 1 to 4), extract 
    % the average resting-state fMRI signal from each parcel for every run (MNI 2mm 
    % space), save into `outdir` in BIDS format.

    start_dir = pwd;
    proj_dir = '/data/project/parcellate_ABCD_preprocessed';

    if(~ischar(scale))
        scale = num2str(scale);
    end

    indir = fullfile(proj_dir, 'data', 'inm7-superds', 'original', 'abcd', 'derivatives', 'abcd-hcp-pipeline');
    if(~exist('outdir', 'var'))
        outdir = fullfile(proj_dir, 'data', 'parcellated_timeseries');
    end

    parcellation = fullfile(proj_dir, 'data', 'Tian2020MSA', '3T', 'Subcortex-Only', ...
        ['Tian_Subcortex_S' scale '_3T.nii.gz']);
    x = MRIread(parcellation);
    % loop through subjects
    subjects = text2cell(subj_ls);
    ses = 'ses-baselineYear1Arm1';
    for i = 1:length(subjects)
        s = subjects{i};
        cd(indir)
        system(sprintf('datalad get -n %s', s));
        system(sprintf('git -C %s config --local --add remote.datalad.annex-ignore true', s));

        cd(fullfile(indir, s, ses, 'func'))
        [flag, msg] = system(sprintf('ls -d %s', [s '_' ses '_task-rest_run-*_space-MNI_bold.nii.gz']))
        if(flag==0)
            runs = strsplit(msg);
            for j = 1:length(runs)
                if(~isempty(runs{j}))
                    runnum = strsplit(runs{j}, '_'); runnum = runnum{4};
                    out_name = fullfile(outdir, s, ses, 'func', [s '_' ses ...
                        '_task-rest_' runnum '_space-MNI_bold_atlas-TianS' scale '.mat']);
                    has_out = system(sprintf('ls -d %s', out_name));
                    if(has_out == 0)
                        fprintf('Output for run %s exists. Skip.\n', runs{j});
                        continue
                    end

                    system(sprintf('datalad get -s inm7-storage %s', runs{j}));

                    ts = MRIread(runs{j});
                    ts.vol = reshape(ts.vol, size(ts.vol,1)*size(ts.vol,2)*size(ts.vol,3), size(ts.vol,4));
                    pts = single(zeros(max(unique(x.vol(:))), size(ts.vol, 2)));
                    for roi = 1:max(unique(x.vol(:)))
                        idx = find(x.vol == roi);
                        pts(roi, :) = mean(ts.vol(idx, :), 1);
                    end
                    mkdir(fullfile(outdir, s, ses, 'func'))
                    save(out_name, 'pts', '-v7.3')
                    system(sprintf('datalad drop %s', runs{j}));
                end
            end
        end

        cd(indir)
        system(sprintf('datalad uninstall %s', s));
        
    end
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