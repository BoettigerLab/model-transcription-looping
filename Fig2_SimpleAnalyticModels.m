%% Analytic calculations from Markov Models
% compute the stationary distribution of some simple Markov chains and
% small versions of the condensate model analytically.
% 
% For additional analytic approaches, see Durrett 2012: Essentials of
% Stochastic Processes.
% For analytic approaches to dynamic properties including first passage
% times, see Boettiger 2013 and Boettiger et al. 2011. 

% SetFigureSavePath('C:\Shared\Documents\Jordan Looping Model\Revision1\Images\');
%% 
k_on = sym('k_on',{'real','positive'});
k_off = sym('k_off',{'real','positive'});
A = sym('A',{'real','positive'});


% 2 state
M2 = [-k_on*A, k_on*A;
    k_off,   -k_off  ];

p2 = null(M2'); % not normalized!
p2 = p2./(sum(p2)) % now it is. 

k_on1 = sym('k_on1',{'real','positive'});
k_off1 = sym('k_off1',{'real','positive'});
k_on2 = sym('k_on2',{'real','positive'});
k_off2 = sym('k_off2',{'real','positive'});
A = sym('A',{'real','positive'});

% 3 state
M3 = [-k_on1*A, k_on1*A, 0;
      k_off1, -k_off1 - k_on2*A, k_on2*A;
       0, k_off2,   -k_off2  ];

pSym = null(M3'); % warning, not normalized!
pSym = pSym./(sum(pSym)) % now it is. 

%%
sc = 1e-2;
on = .005*sc^2;
of1 = .05*sc;
of2 = .005*sc;
sqrt(of1.*of2./on^2)
f = subs(pSym(3),[k_on1,k_on2,k_off1,k_off2],[on,on,of1,of2]);
f2 = subs(p2(end),[k_on,k_off],[on,mean([of1,of2])]);
h = A^2./( (of1.*of2./on^2) + A^2);

f1 = figure(1); clf;
fplot(h,[0,10/sc]); hold on;
fplot(f,[0,10/sc]);
fplot(f2,[0,10/sc]);
legend('hill n=2','3-state','2-state','Location','Best');
ylim([0,1]);
xlabel('total A (molecules)'); ylabel('probability bound');
set(gcf,'color','w');
% set(gca,'xscale','log','yscale','log');

eval(subs(f,A,200))./eval(subs(f,A,100))
eval(subs(f,A,100))./eval(subs(f,A,50))

x1 = 200; x2 = 100;
y1 = eval(subs(f,A,x1));
y2 = eval(subs(f,A,x2));
r1 =1+(y1-y2)./y2 
r2 =1+(x1-x2)./x2
r1/r2

x1 = 200; x2 = 160;
y1 = eval(subs(f,A,x1));
y2 = eval(subs(f,A,x2));
r1 =1+(y1-y2)./y2 
r2 =1+(x1-x2)./x2
r1/r2

figure(2); clf; fplot(diff(f),[0,10/sc]);

% SaveFigure(f1,'name','analytic_3state','formats',{'eps','png'},'overwrite',1);

%% Simple condensate model
% [-a-e              a+e              0              0	0		0		0		0...
% [ r           -1+(1-a)2-e-r  1-(1-a)2+e     0 	0		0		0		0…
% [ …
% [ … 				 		 r	-1+(1-a)n -e-r    1-(1-a)n +e    		0 ... ]
clear M
syms a e r M
cMax = 6;  % Warning, analytic calculations with a large cMax are expensive and slow
% Build Matrix
% boundaries
M(1,1) = -a-e;
M(1,2) = a+e;
M(cMax,cMax-1) = r;
M(cMax,cMax) = -r;
for n=2:cMax-1 
    M(n,n) = -1+(1-a)^n -e-r; % main diag
    M(n,n-1) = r;  % lower diag
    M(n,n+1) = 1-(1-a)^n +e ; % upper diag
end
disp(M)
M2 = subs(M,[a,e,r],[addPol,enPr,losePol])
p2 = eval(null(M2')); p2 = p2./sum(p2(:))  % perfect this checks out. 
% in symbolic form
pSym = null(M'); % 
pSym = pSym./(sum(pSym)) % symbolic closed form solution
p3 = eval(subs(pSym,[a,e,r],[addPol,enPr,losePol])); % probability of all 6 states

N = 100; % number of contact frequencies to test
contactFreqs = linspace(0,1,N);
probActive = zeros(N,1);
for n=1:N
    enPr = contactFreqs(n);
    probActive(n) = eval(subs(pSym(end),[a,e,r],[addPol,enPr,losePol])); % probability of all 6 states
end
figure(2); clf; plot(contactFreqs,probActive);
xlabel('E-P contact frequency');
ylabel('probability of final (active) state');
