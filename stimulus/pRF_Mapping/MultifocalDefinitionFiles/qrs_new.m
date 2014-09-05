%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% qrs_new
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [moo,serind,red]=qrs_new(Nreg,prime)

% Gives the matrix on onsets as logical indices and serial indices for a given number of visual field areas
% Nreg is region count, N is the block count
% Nreg = 48; % Region total is 48
% N = 59; % First prime number greater than region count which full fills N = 4m -1; where m is an integer
% red is the 5thlevel CS sqeuence used for region switching

if nargin<2
	prime=17; % Shift prime
end

Kini=ceil((Nreg)/4); %
Nini=4*(1:(Kini))-1;
N=max(Nini(isprime(Nini)));

%Find first prime in the series 4K-1 above amount of regions
while N<Nreg
	Kini=Kini+1;
	Nini=4*(1:(Kini))-1;
	N=max(Nini(isprime(Nini)));
end

% Index active times for t
j=[0:N-1]; %Index
t=zeros(1,N); %t stimuli to all zero
t((mod(j.^2,N))+1)=1; %Stimuli impossed; if after j^2+1 impossed, a remainder exists, than 1, else zero

% Rotate for other regions
moo=zeros(Nreg,N);
serind=zeros(Nreg,length(find(t==1)));
t1=t;

for i=1:Nreg
	moo(i,:)=t1;
	numind=ind2sub([1,N],find(t1==1));
	serind(i,:)=numind;
	
	shift=mod(i*prime,Nreg); % Use N instead of Nreg if wrapping at the end of time instead of end of regions
	t1(1:shift)=t(end-shift+1:end);
	t1(shift+1:end)=t(1:end-shift);
	
end

% Rudy to define odd/even sector reversal pattern with maximum
% counterbalancing for given size and a base two
N2=63;
rudy = (((N2+1)*2)-1); %ie 128 or k = 7;
j2=[0:rudy-1]; %Index
t2=zeros(1,rudy); %t stimuli to all zero
t2((mod(j2.^2,rudy))+1)=1;
t2(1,rudy+1)=0; %complete 50/50 on off


red(1,:)=t2;


% Transpose to get time in rows
moo=moo';
serind=serind';

% Test identical time sequences and report
for i=2:Nreg;tf_regdat(i)=isequal(moo(:,1),moo(:,i));end;

if max(tf_regdat) == 0
	disp('No correlating time sequences');
else
	error('Correlating time sequences found! Check qrs_new input parameters.')
end
