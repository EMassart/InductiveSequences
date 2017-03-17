function [new_eta, step_size, xnew, fnew, gfnew, lnf, lng, lnV, lnR,ls,lnftime,lngtime,lnRtime,lnVtime] = choose_step_size(fns, eta, x1, f1, gradf1, params)

% Code extracted from http://www.math.fsu.edu/~whuang2/papers/ARLBACMGM.htm
% References: 'A Riemannian Limited-Memory BFGS Algorithm for
% Computing the Matrix Geometric Mean', X. Yuan, W. Huang, P.-A. Absil, K.
% Gallivan.
% (Slight modifications between this version and the original are possible)

lnftime = 0;
lngtime = 0;
lnRtime = 0;
lnVtime = 0;


C = max(x1.c);
L = 1 + 0.5 * log(C);
step_size = 2/(1+L);


new_eta = fns.STV(x1, step_size, eta);

lnRtimestart = toc;
xnew = fns.R(x1,new_eta);
lnRtime = lnRtime + toc - lnRtimestart;

lnR = 1;
lnftimestart = toc;
[fnew, xnew] = fns.f(xnew);
lnftime = lnftime + toc - lnftimestart;
lnf = 1;

lngtimestart = toc;
[gfnew, xnew] = fns.gf(xnew);
lngtime = lngtime + toc - lngtimestart;
lng = 1;

lnV = 0;
ls = 1;
end