
snr_y = -5:1:20;
snr_de = 10.^(snr_y/10);
global p;global m;
p =5;m = 1;
global q
q = p^m;
global delta
delta = (0:1:q-1)-(q-1)/2;

% delta = (0:q-1)-(q-1)/2;
delta = delta/sqrt(mean(delta.^2));
del_delta = (delta(end)-delta(1))*(q)/(q-1);
Es = sum(delta.^2);
sigma_n = sqrt((Es/q)./(10.^(snr_y/10)));
y = -5*q:0.01:5*q;

% MI = zeros(1,length(snr_y));
%% PAM
C = zeros(1,length(snr_y));

for j = 1:1:length(snr_y)
    py = zeros(size(y));
    for i = 1:1:q
        py = 1/(q*sqrt(2*pi*sigma_n(j)^2))*exp(-(y-delta(i)).^2/(2*sigma_n(j)^2))+py;
    end
    C(j) = -sum(py.*log2(py+eps))*0.01-0.5*log2(2*pi*exp(1)*sigma_n(j)^2);
end

l = find(C==max(C),1,'first')-1;
c_snr =  griddedInterpolant(snr_y,C);
snr_c = griddedInterpolant(C(1:l),snr_y(1:l));



%% int s
C2 = zeros(1,length(snr_y));
y = -10*q:0.01:10*q;
delta = delta/sqrt(mean(delta.^2));
s = [-1 0 1];
s = delta;
n_s = length(s);
% c-s_q
sub_mat = delta' - s;
% mod delta
k_max = floor( (max(sub_mat,[],'all')-delta(1))/del_delta);
k_min = floor( (min(sub_mat,[],'all')-delta(1))/del_delta);
ec = delta'+del_delta*(k_min:k_max);
[~,index1] = min(abs(sub_mat(:)'-ec(:)));
index = mod(index1,q);
index(index==0)=q;
[~,isq]=min(abs(s-ec(:)));
r =  reshape(delta(index),q,[])+ec(isq);
% iout = find(abs(sub_mat)>(delta(q)+del_delta)/2);
% x= sub_mat;
% x(iout)= sub_mat(iout)' - (ec(index1(iout)));
% r = x+s;

r_set = unique(r)';
n_r = length(r_set);
mat_rc = zeros(q,n_r);
for i = 1:1:q
    mat_rc(i,:) = histcounts(r(i,:),[r_set,10]);    
end

% counts
p_r = sum(mat_rc,1)/(sum(mat_rc,'all')+eps);
p_rc = mat_rc./(sum(mat_rc,2)+eps);

for j = 1:1:length(snr_y)
    p_yr=1/(sqrt(2*pi*sigma_n(j)^2))*exp(-(y-r_set').^2/(2*sigma_n(j)^2));
    p_y = p_r*p_yr;
    p_yc = p_rc*p_yr;
    C2(j) = -sum(p_y.*log2(p_y+eps))*0.01+sum(sum(p_yc.*log2(p_yc+eps)))/q*0.01;
end
l = find(C2==max(C2),1,'first')-1;
c_snr2 =  griddedInterpolant(snr_y,C2);
snr_c2 = griddedInterpolant(C2(1:l),snr_y(1:l));

%% 
C3 = zeros(1,length(snr_y));
y = -10*q:0.01:10*q;
delta = [  -1.5748   -0.1417         0    0.1417    1.5748];
del_delta =   3.4495;
delta = delta/sqrt(mean(delta.^2));
s = [-2 0 2];
n_s = length(s);
% c-s_q
sub_mat = delta' - s;
% mod delta
k_max = floor( (max(sub_mat,[],'all')-delta(1))/del_delta);
k_min = floor( (min(sub_mat,[],'all')-delta(1))/del_delta);
ec = delta'+del_delta*((k_min:k_max));
% ec = unique(sort(sub_mat));
[~,index1] = min(abs(sub_mat(:)'-ec(:)));
index = mod(index1,q);
index(index==0)=5;
% r = sub_mat-reshape(delta(index),q,[])+s_q;
[~,isq]=min(abs(s-ec(:)));
% err = ec(isq)-s
r =  reshape(delta(index),q,[])+ec(isq);
% iout = find(abs(sub_mat)>(delta(q)+del_delta)/2);
% x= sub_mat;
% x(iout)= sub_mat(iout)' - (ec(index1(iout)));
% r = x+s;

r_set = unique(r)';
n_r = length(r_set);

mat_rc = zeros(q,n_r);
for i = 1:1:q
    mat_rc(i,:) = histcounts(r(i,:),[r_set,10]);    
end

% counts
p_r = sum(mat_rc,1)/(sum(mat_rc,'all')+eps);
p_rc = mat_rc./(sum(mat_rc,2)+eps);

for j = 1:1:length(snr_y)
    p_yr=1/(sqrt(2*pi*sigma_n(j)^2))*exp(-(y-r_set').^2/(2*sigma_n(j)^2));
    p_y = p_r*p_yr;
    p_yc = p_rc*p_yr;
    C3(j) = -sum(p_y.*log2(p_y+eps))*0.01+sum(sum(p_yc.*log2(p_yc+eps)))/q*0.01;
end
l = find(C3==max(C3),1,'first')-1;
c_snr3 =  griddedInterpolant(snr_y,C3);
snr_c3 = griddedInterpolant(C3(1:l),snr_y(1:l));
% c3 = plot(snr_y,C3);

%% NUPAM
delta =  [  -1.5102   -0.4682         0    0.4682    1.5102]
C4 = zeros(1,length(snr_y));
for j = 1:1:length(snr_y)
        py = zeros(size(y));
        for i = 1:1:q
            py = 1/(q*sqrt(2*pi*sigma_n(j)^2))*exp(-(y-delta(i)).^2/(2*sigma_n(j)^2))+py;
        end
        C4(j) = -sum(py.*log2(py+eps))*0.01-0.5*log2(2*pi*exp(1)*sigma_n(j)^2);    
end


%%
figure()
c1 = plot(snr_y,C,'--');hold on
c2 = plot(snr_y,C2);
c3 = plot(snr_y,C3);
c4 = plot(snr_y,C4);
% legend('PAM','C2');
xlabel('SNR/dB');
ylabel('MI');
title([num2str(q),'PAM']);
% axis([-5 10 0 2]);
%%
% figure()
% pq2c1 =  plot(snr_y,q2c1,'--');hold on
% pq2c2 =  plot(snr_y,q2c2);
% pq5c1 =  plot(snr_y,q5c1,'--');
% pq5c2 =  plot(snr_y,q5c2);
% pq5c3 =  plot(snr_y,q5c3);
% pq5c4 =  plot(snr_y,q5c4);
% % legend([pq2c1,pq2c2,pq5c1,pq5c2,pq5c3],'q=2','q=2','q=5','q=5','q=5,NP')
% % c3 = plot(snr_y,C3);
% % legend('PAM','C2');
% xlabel('SNR/dB');
% ylabel('MI');
% % title([num2str(q),'PAM']);


% %%
% global reJfunc
% load(['J',num2str(q),'.mat'])
% mat_prod = mod([0:q-1].*[0:q-1]',q);

