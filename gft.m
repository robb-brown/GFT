1; % lol. Tells Matlab this is not a single function file

function partitions = octavePartitions(N)
	widths = 2.^(0:round(log(N)/log(2)-2));
	widths = [[1],widths,flip(widths)];
	partitions = [[0],cumsum(widths)]+1;
endfunction;

function widths = partitionWidths(partitions)
	widths = [shift(partitions,-1) - partitions];
	widths(length(partitions)) += max(partitions);
endfunction

function windows = boxcarWindows(partitions)
	windows = ones(1,max(partitions));
endfunction;
	
function SIG = GFT(sig,partitions,windows)
	SIG = fft(complex(sig));
	SIG .*= windows;
	for p = 1:(length(partitions)-1)
		SIG(partitions(p):partitions(p+1)-1) = ifft(SIG(partitions(p):partitions(p+1)-1));
	endfor;
endfunction;

function spectrogram = interpolateGFT(SIG,partitions,tAxis,fAxis,method='linear')
	% Interpolate a 1D GFT onto a grid. If axes is specified it should be a
	% list or tuple consisting of two arrays, the sampling points on the time and frequency
	% axes, respectively. Alternatively, M can be specified, which gives the number
	% of points along each axis.
	N = length(SIG);
	widths = partitionWidths(partitions);
	spectrogram = complex(length(partitions),zeros(length(tAxis)));
	% interpolate each frequency band in time
	for p = 1:length(partitions)
		% indices of sample points, plus 3 extra on each side in case of cubic interpolation
		indices = (-3:widths(p)+2);
		% time coordinates of samples
		t = indices .* (N/widths(p));
		% values at sample points
		if (p < length(partitions))
			f = SIG(partitions(p):partitions(p+1)-1)(mod(indices,widths(p))+1);
		else
			f = SIG(partitions(p):N)(mod(indices,widths(p))+1);
		endif
		if (length(f) > 1)
			spectrogram(p,:) = interp1(t,f,tAxis,method=method);
		else
			spectrogram(p,:) = f;
		endif
	endfor

	% Interpolate in frequency
	indices = mod(-3:length(partitions)+2,length(partitions));
	f = partitions(indices+1) + widths(indices+1)/2; 
	f(1:3) -= N; 
	f(length(f)-2:length(f)) += N;
	t = spectrogram(indices+1,:);
	spectrogram = interp1(f,t,fAxis,method=method);
endfunction


% function demo()
	% Create a fake signal
	N = 256;
	x = linspace(0,1,N);
	sig = zeros(1,length(x));
	
	% signal example 1 (a single delta)
	sig(N/2) = 1.0;
	
	% signal example 2 (a mixture of sinusoids and a delta)
	% sig(1:N/2) += (sin((N/16)*2*pi*x)*1.0)(1:N/2);
	% sig(N/2+1:N) += (cos((N/8)*2*pi*x)*1.0)(N/2+1:N);
	% sig(2*N/16+1:3*N/16) += (sin((N/4)*2*pi*x)*1.0)(2*N/16+1:3*N/16);
	% sig(N/2+N/4+1) = 2.0;

	% Do the transform
	partitions = octavePartitions(N);
	windows = boxcarWindows(partitions);
	SIG = GFT(sig,partitions,windows);
	
	% Interpolate to get a spectrogram
	% The third and fourth parameters set the time and frequency axes respectively,
	% and can be changed to raise or lower the resolution, or zoom in on
	% a feature of interest
	spectrogram = interpolateGFT(SIG,partitions,1:N,1:N);
	
	% Display
	figure(); 
	subplot(3,1,1);
	plot(x,sig,';signal;');
	ax = subplot(3,1,2);
	hold on;
	for p = partitions
		line([x(p),x(p)],[0,max(abs(SIG))],'Color',[1 0 0],'linestyle','--');
	endfor;
	plot(x,abs(SIG),';SIGNAL;');
	subplot(3,1,3);
	imagesc(abs(spectrogram));
% endfunction;