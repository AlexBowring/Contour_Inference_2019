nSubj=500;
nRlz=5;

cd /storage/maullz/Contour_Inference_2019/
addpath('/storage/essicd/spm8/')
addpath(genpath('/storage/maullz/Contour_Inference_2019/'))

rng('shuffle')
tID=str2num(getenv('SGE_TASK_ID'))
Sim_08(nSubj,['Sim_08_' num2str(nSubj) '_subjects_' sprintf('%03d',tID)],nRlz)
