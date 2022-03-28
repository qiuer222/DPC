q = 5;
% n = floor((q-1)/2);
di = 0.01;
delta = zeros(1,q)
snr_y = 4:0.1:15;
y = -5*q:0.01:5*q;
sigma_n = sqrt((1)./(10.^(snr_y/10)));
C = zeros(1,length(snr_y));
MI_goal = 1;
figure;
[~,now] = size(snr_y);
for delta1 = -1:di:0
    
    %     for delta2 = delta1:di:-di
    %         for delta3 = -1:di:0
    %             for delta4 = -1:di:0
    %             for delta4 = -1:di:1
    %                 for delta5 = -1:di:1
    %                     delta = [delta1 delta2 delta3 delta4 delta5]
    %                     delta = [delta1 delta2 delta3  delta4 -delta4 -delta3 -delta2 -delta1]
    delta = [-1 delta1 0 -delta1 1]
    %         delta = delta/sqrt(mean(delta.^2)+eps)*5
    [~,isu] = size(unique(delta))
    if(isu~=q)
        continue;
    end
    %%
    delta = delta/sqrt(mean(delta.^2))
    for del_delta = (delta(end)-delta(1)):0.1:(delta(end)-delta(1))*(q+2)/(q-1)
        
        s_q = [-2  2];
        % s_q = delta(1:4);
        n_s = length(s_q);
        % c-s_q
        sub_mat = delta' - s_q;
        % mod delta
        k_max = floor( (max(sub_mat,[],'all')-delta(1))/del_delta);
        k_min = floor( (min(sub_mat,[],'all')-delta(1))/del_delta);
        ec = delta'+del_delta*(k_min:k_max);
        % ec = unique(sort(sub_mat));
        [~,index1] = min(abs(sub_mat(:)'-ec(:)));
        index = mod(index1,q);
        index(index==0)=5;
        % r = sub_mat-reshape(delta(index),q,[])+s_q;
        [~,isq]=min(abs(s_q-ec(:)));
        r =  reshape(delta(index),q,[])+ec(isq);
        r_set = unique(r)';
        % ec = ec(min(index):max(index));
        
        n_r = length(r_set);
        % r = reshape(delta(index),q,[])+s_q;
        mat_rc = zeros(q,n_r);
        for i = 1:1:q
            mat_rc(i,:) = histcounts(r(i,:),[r_set,10]);
        end
        
        % counts
        p_r = sum(mat_rc,1)/(sum(mat_rc,'all')+eps);
        p_rc = mat_rc./(sum(mat_rc,2)+eps);
        C2 = zeros(1,length(snr_y));
        for j = 1:1:length(snr_y)
            p_yr=1/(sqrt(2*pi*sigma_n(j)^2))*exp(-(y-r_set').^2/(2*sigma_n(j)^2));
            p_y = p_r*p_yr;
            p_yc = p_rc*p_yr;
            C2(j) = -sum(p_y.*log2(p_y+eps))*0.01+sum(sum(p_yc.*log2(p_yc+eps)))/q*0.01;
            if(C2(j)>MI_goal)
                if j<now
                    now = j;
                    delta_now = delta;
                    del_delta_now = del_delta
                    plot(snr_y(1:j),C2(1:j));hold on;         
                end
                break
            end
        end
%         plot(snr_y(1:j),C2(1:j));hold on;      
        %%
        %     temp = find(C2>0.5*log2(q),1);
        %     if temp<now
        %         now = temp
        %         delta_now = delta;
        %         plot(snr_y,C2);hold on;
        %     end
%         plot(snr_y(1:j),C2(1:j));hold on
    end
    
end
%         end
%     end
%         end
%     end
% end
%%
delta = delta_now
del_delta = del_delta_now
snr_y = -5:0.1:20;
sigma_n = sqrt((1)./(10.^(snr_y/10)));
        % s_q = delta(1:4);
        n_s = length(s_q);
        % c-s_q
        sub_mat = delta' - s_q;
        % mod delta
        k_max = floor( (max(sub_mat,[],'all')-delta(1))/del_delta);
        k_min = floor( (min(sub_mat,[],'all')-delta(1))/del_delta);
        ec = delta'+del_delta*(k_min:k_max);
        % ec = unique(sort(sub_mat));
        [~,index1] = min(abs(sub_mat(:)'-ec(:)));
        index = mod(index1,q);
        index(index==0)=5;
        % r = sub_mat-reshape(delta(index),q,[])+s_q;
        [~,isq]=min(abs(s_q-ec(:)));
        r =  reshape(delta(index),q,[])+ec(isq);
        r_set = unique(r)';
        % ec = ec(min(index):max(index));
        
        n_r = length(r_set);
        % r = reshape(delta(index),q,[])+s_q;
        mat_rc = zeros(q,n_r);
        for i = 1:1:q
            mat_rc(i,:) = histcounts(r(i,:),[r_set,10]);
        end
        
        % counts
        p_r = sum(mat_rc,1)/(sum(mat_rc,'all')+eps);
        p_rc = mat_rc./(sum(mat_rc,2)+eps);
        C2 = zeros(1,length(snr_y));
        for j = 1:1:length(snr_y)
            p_yr=1/(sqrt(2*pi*sigma_n(j)^2))*exp(-(y-r_set').^2/(2*sigma_n(j)^2));
            p_y = p_r*p_yr;
            p_yc = p_rc*p_yr;
            C2(j) = -sum(p_y.*log2(p_y+eps))*0.01+sum(sum(p_yc.*log2(p_yc+eps)))/q*0.01;
        end