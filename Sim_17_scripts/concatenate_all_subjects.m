function ConcatenateSims
Sim = 'Sim_17'; 
Nfiles = 100;
mag = [0.2, 0.7, 1, 2, 5];
smo = [1, 3, 5, 7];
home_dir = '/storage/maullz/Contour_Inference_2019/';

for j = 1:5
    for k = 1:4
    Base = fullfile(home_dir, 'Sim_17_results', [num2str(erase(num2str(mag(j)),'.')) '_effect_size'], [num2str(smo(k)) '_smoothing/']);
	cd(Base)

	a = dir([Sim,'_*.mat']);

	x = load(a(1).name);
    filename = strsplit(a(1).name,'_');
    filename = filename(1:end-1);
    filename = strjoin(filename,'_');
	filename = sprintf([filename,'.mat']);
	save(filename,'-struct','x');

        for l=2:Nfiles
        x = load(filename);
        y = load(a(l).name);

        vrs = fieldnames(x);

        x.(vrs{1}) = x.(vrs{1}) + y.(vrs{1}); % Adding number of realizations

            for m=5:58
                x.(vrs{m}) = [x.(vrs{m}); y.(vrs{m})]; % Concatenates the Monte Carlo max values
            end

        save(filename,'-struct','x');
        end
    end
end

end
