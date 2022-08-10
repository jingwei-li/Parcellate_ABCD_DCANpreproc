function compute_RSFC_with_censor(scale, subj_ls, roi_ts_dir, out_censor_mat)
% compute_RSFC_with_censor(scale, subj_ls, roi_ts_dir, out_censor_mat)
%
% Compute ROI-to-ROI functional connectivity. The ROIs are defined by a combination of
% cortical Schaefer's parcellation and subcortical Tian's parcellation at a certain scale.
% For each subject, the computed functional connectivity will be saved in each subject folder
% under `roi_ts_dir`, following BIDS format.
%
% Input:
% - scale: choose from 1 to 10. 
%          Scale 1: Schaefer parcellation with 100 areas + Tian parcellation at scale 1 (16 ROIs)
%          Scale 2: Schaefer parcellation with 200 areas + Tian parcellation at scale 2 (32 ROIs)
%          Scale 3: Schaefer parcellation with 300 areas + Tian parcellation at scale 3 (50 ROIs)
%          From 4 to 10: Schaefer parcellation with (scale * 100) areas + Tian parcellation at scale 4 (54 ROIs)
%
% - subj_ls
%   Subject list who have preprocessed resting-state fMRI data.
%
% - roi_ts_dir
%   Full-path directory containing the parcellated timeseries (output folder of `extract_rest_timeseries_*.m`).
%
% - out_censor_mat
%   Output .mat filename containing which subjects, which runs passed motion censoring, which are not, and which
%   runs only have fsLR32k space files but not MNI space files.
%

start_dir = pwd;
proj_dir = '/data/project/parcellate_ABCD_preprocessed';
data_dir = fullfile(proj_dir, 'data', 'inm7-superds', 'original', 'abcd', 'derivatives', 'abcd-hcp-pipeline');

subjects = text2cell(subj_ls);
ses = 'ses-baselineYear1Arm1';
FD_threshold = 0.3;

Schaefer_res = 100*scale;
if(scale<4)
    Tian_res = scale;
else
    Tian_res = 4;
end

subjects_pass = {};
pass_runs = cell(1, length(subjects)); 
unpass_runs = cell(1, length(subjects)); 
noMNI_runs = cell(1, length(subjects));
for i = 1:length(subjects)
    s = subjects{i};
    fprintf('%s\n', s)
    pass_runs{i} = {};    unpass_runs{i} = {};   noMNI_runs{i} = {};

    cd(fullfile(roi_ts_dir, s, ses, 'func'))
    [flag, msg] = system(sprintf('ls -d %s', [s '_' ses '_task-rest_run-*_bold_atlas-Schaefer' ...
        num2str(Schaefer_res) '_timeseries.mat']));
    if(flag==0)
        cd(data_dir)
        system(sprintf('datalad get -n %s', s));
        system(sprintf('git -C %s config --local --add remote.datalad.annex-ignore true', s));

        runs = strsplit(msg);
        runs = runs(~cellfun(@isempty, runs));
        censor = nan(1, length(runs)); % for each run of current subject, does it pass the censoring threshold?
        
        currFC = cell(1, length(runs));

        for j = 1:length(runs)
            runnum = strsplit(runs{j}, '_'); runnum = runnum{4};
            run_MNI = fullfile(roi_ts_dir, s, ses, 'func', [s '_' ses '_task-rest_' ...
                runnum '_space-MNI_bold_atlas-TianS' num2str(Tian_res) '.mat']);
            if(~exist(run_MNI, 'file'))
                noMNI_runs{i} = [noMNI_runs{i} {runnum}];
                continue;
            end
            out_name = fullfile(roi_ts_dir, s, ses, 'func', [s '_' ses '_task-rest_' runnum ...
                '_RSFC_Schaefer' num2str(Schaefer_res) '_Tian' num2str(Tian_res) '.mat']);
            
            cd(fullfile(data_dir, s, ses, 'func'))
            mt_tsv = [s '_' ses '_task-rest_' runnum '_desc-includingFD_motion.tsv'];
            system(sprintf('datalad get -s inm7-storage %s', mt_tsv));
            % for some run, the "desc-includingFD_motion.tsv" file doesn't exist
            % calculate FD from 6 motion parameters
            if(~exist(mt_tsv))
                mt_tsv = [s '_' ses '_task-rest_' runnum '_motion.tsv'];
                system(sprintf('datalad get -s inm7-storage %s', mt_tsv));
                system(sprintf('cat %s | tr -s ''([\t]+)'' '','' > tmp.tsv', mt_tsv)) % replace multiple \t to a single comma
                % add comma to the beginning of the first line (because there are extra tabs from line 2 in the original file); 
                % remove the first comma of each line (remove the extra tab); 
                % remove the last comma of each line (because there are extra tabs at the end of each line in the original file)
                system('echo ",$(cat tmp.tsv)" > tmp2.tsv; cut -c 2- < tmp2.tsv > tmp3.tsv; sed -i ''s/.$//'' tmp3.tsv')
                mt = tdfread('tmp3.tsv', ',')
                % for some run, there isn't extra tab at the end of first line. Therefore the previous step would remove the 't'
                if(~isfield(mt, 'RotZDt'))
                    [mt.RotZDt] = mt.RotZD;
                    mt = rmfield(mt, 'RotZD');
                end
                mt.framewise_displacement = abs(mt.XDt) + abs(mt.YDt) + abs(mt.ZDt) + 50*pi/360 * (abs(mt.RotXDt) + abs(mt.RotYDt) + abs(mt.RotZDt));
            else
                mt = tdfread(mt_tsv, ' ');
            end
            
            cd(fullfile(roi_ts_dir, s, ses, 'func'))

            FD_outlier = mt.framewise_displacement>FD_threshold;
            if(length(find(FD_outlier)) > length(mt.framewise_displacement)/2)
                censor(j) = 1;
                unpass_runs{i} = [unpass_runs{i} {runnum}];
                % write a subfunction to calculate RSFC for this run
            else
                censor(j) = 0;
                pass_runs{i} = [pass_runs{i} {runnum}];

                % compute FC for current run with censoring
                ts_cort = load(runs{j});
                ts_subcort = load(run_MNI);
                currFC{j} = FC_per_run(ts_cort.pts, ts_subcort.pts, FD_outlier);
            end
        end

        if(any(censor==0))
            subjects_pass = [subjects_pass {s}];
            corr_mat = mean(cat(3, currFC{~cellfun(@isempty, currFC)}), 3);
            save(out_name, 'corr_mat')
        end
    else
        censor = [];
    end

end

save(out_censor_mat, 'subjects', 'subjects_pass', 'pass_runs', 'unpass_runs', 'noMNI_runs')
[subj_ls_dir, subj_ls_base] = fileparts(out_censor_mat);
cell2text(subjects_pass, fullfile(subj_ls_dir, [subj_ls_base '.txt']))


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

function cell2text(cell_var, filename)

    fid = fopen(filename, 'w');
    formatSpec = '%s\n';
    
    for row = 1:length(cell_var);
        fprintf(fid, formatSpec, cell_var{row});
    end
    
    fclose(fid);
end

function corr_mat = my_corr(X, Y)

    % Calculate correlation matrix between each column of two matrix.
    % 
    % 	corr_mat = my_corr(X, Y)
    % 	Input:
    % 		X: D x N1 matrix
    % 		Y: D x N2 matrix
    % 	Output:
    % 		corr_mat: N1 x N2 matrix    
    
    X = bsxfun(@minus, X, mean(X, 1));
    X = bsxfun(@times, X, 1./sqrt(sum(X.^2, 1)));
    
    Y = bsxfun(@minus, Y, mean(Y, 1));
    Y = bsxfun(@times, Y, 1./sqrt(sum(Y.^2, 1)));
    
    corr_mat = X' * Y;
end

function FC = FC_per_run(ts_cort, ts_subcort, FD_outlier)

    % ts_cort: Schaefer-parcellated timeseries from cortex. #ROI_1 x T matrix
    % ts_subcort: Tian-parcellated timeseries from subcortex. #ROI_2 x T matrix
    % FD_outlier: 0/1 vector with length = T

    % cortical to cortical
    cort = my_corr(ts_cort(:,FD_outlier==0)', ts_cort(:,FD_outlier==0)');
    subcort = my_corr(ts_subcort(:,FD_outlier==0)', ts_subcort(:,FD_outlier==0)');
    cort_subcort = my_corr(ts_cort(:,FD_outlier==0)', ts_subcort(:, FD_outlier==0)');
    FC = [[cort cort_subcort]; [cort_subcort' subcort]];
end