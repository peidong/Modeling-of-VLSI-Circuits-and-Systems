clear all;
clc;

% %% mc
fp = [];
fom = [];
sample_n=0;
counter = 0;
sample_unit = 100;
total_error_counter = 0;
threshold = 0.202e-9;
stop_fom = 0.05;
hspice_path = 'C:/synopsys/Hspice_C-2009.03-SP1/BIN\hspice.exe';

while(sample_n(end) < 20000)
    
	counter = counter+1;
	sample_n = [sample_n, sample_n(end)+sample_unit];

	% generate MC or QMC samples
	% call MC_gene or QMC_gene
	
	% write MC and QMC sample into sweep file (Not for 6-T SRAM application)
	fid1 = fopen('sweep_data_mc','w');
	fprintf(fid1, '.DATA data\n');
	fprintf(fid1, 'nvth1 nvth2 nvth3 nvth4 pvth1 pvth2 \n');
	for i = 1:sample_unit
			for j = 1:6
					fprintf(fid1, '%e\t', samples(i,j));
			end
			fprintf(fid1, '\n');    
	end
	fprintf(fid1, '.ENDDATA\n');
	fclose(fid1);

	% call HSPICE simulation
	dos([hspice_path, ' -i sram.sp']);
	x = loadsig('chain.tr0');

	%data processing
	error_counter = 0;
	for i = 1:sample_unit
	% get the performance of one simulation
	% Hereadd the code to get the performance a

		if (a > threshold)
			error_counter = error_counter+1;
		end

	end
	total_error_counter = total_error_counter + error_counter;
	fp = [fp total_error_counter/sample_n(end)];

	fom = [fom std(fp)/mean(fp)];

	str = sprintf('%d out of %d samples failed(%d/%d), failure rate = %e, FOM = %e', total_error_counter, sample_n(end), error_counter, sample_unit, fp(end), fom(end));
	disp(str);
	
	if(counter>=5)	% first 5 FOM may be not accurate
		if(fom(end)<=stop_fom)
			break;
		end
	end
    
end

%%
figure(1)
semilogx(sample_n(2:end),fp, '-*');
figure(2)
semilogx(sample_n(2:end),fom, '-*');




