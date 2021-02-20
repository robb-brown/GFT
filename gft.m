1; % lol. Tells Matlab this is not a single function file

function ker = fourierKernel(signal,k)
	x = complex(linspace(0,1,length(signal)),0);
	ker = e.^(1j*2*pi*k*x);
endfunction;

function win = window(signal,t,k)
	l = length(signal);
	x = complex(linspace(0,1,l),0);
	win = 1*abs(k)/sqrt(2*pi)*e.^(-(x-0.5).^2*abs(k).^2/2);
	win = win./(sum(win));
	win = shift(win,-l/2);
	win = shift(win,t);
endfunction;
	
% @ operator creates a pointer to a function
function gft = GFT(signal,kernel=@fourierKernel,window=@window,windowmult=0.8,windowcut=16)
	transform = {};
	transformn = {};
	sig = signal;
	N = length(sig);
	k = length(sig); 			% will get divided by 2 before we start
	while k > 1
		k = k / 2;
		M = length(sig)
		sigma = k.*self.windowmult;
		ker = kernel(sig,k);
		kern = kernel(sig,-k);
		%printf("k= %d  M= %d" % (k,M))
		sk = sig.*ker;
		skn = sig.*kern;
		line = {}
		linen = {}
		for i = 1:M
			win = window(sig,i,sigma,windowcut);			
			line[i] = sum(sk.*win);		
			linen[i] = sum(skn.*win);
		endfor;
		transform[length(transform)+1] = line;
		transformn[length(transformn)+1]] = linen;
		if (k <= N/4 & M > 2)
			sig = downsampleft(sig,2);
		endif;
	endwhile;
	transform.append([sum(sig)])
	return (transform,transformn)


endfunction;


function demo()
	% Create a fake signal
	N = 256;
	x = linspace(0,1,N);
	sig = zeros(1,length(x),'double');
	sig(1:N/2) += (sin((N/16)*2*pi*x)*1.0)(1:N/2);
	sig(N/2:N) += (cos((N/8)*2*pi*x)*1.0)(N/2:N);
	sig(2*N/16:3*N/16) += (sin((N/4)*2*pi*x)*1.0)(2*N/16:3*N/16);
	sig(N/2+N/4) = 2.0;
	
	figure(); hold on;
	plot(x,sig);
	plot(x,window(sig,N/2,2));
	plot(x,fourierKernel(sig,2))
	
	ST = GFT(sig)
		
endfunction;