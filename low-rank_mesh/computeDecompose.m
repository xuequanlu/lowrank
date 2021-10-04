
integer = 1:40000;
integer = integer';
num = numel(integer);

out = zeros(num,2);
for i=1:num
    tmp = i*3;
    facs = factor(tmp);
    num_facs = numel(facs);
    %dis = 1e10;
    
    if num_facs==1
        out_tmp = zeros(1,2);
        out_tmp(1,1) = facs(1);
        out_tmp(1,2) = 1;
    
    elseif num_facs==2
        out_tmp = zeros(1,2);
        out_tmp(1,1) = facs(2);%max
        out_tmp(1,2) = facs(1);
        
    elseif num_facs==3
        num_tmp = 3*2/(2*1*1);
        out_tmp = zeros(num_tmp, 2);
        out_tmp(:,1) = 0; out_tmp(:,2) = 1e10;
        count = 0;
        % C_3_1
        for j=1:num_facs
            mult = 1;
            for k=1:num_facs
                if j~=k % the remaining factors
                    mult = mult*facs(k);
                end
            end
            count = count + 1;
            out_tmp(j,:) = [facs(j), mult];
        end
        
    elseif num_facs==4
        num_tmp = 4*3*2*1/(1*3*2*1) + 4*3*2*1/(2*1*2*1);
        out_tmp = zeros(num_tmp, 2);
        out_tmp(:,1) = 0; out_tmp(:,2) = 1e10;
        count = 0;
        % C_4_1
         for j=1:num_facs
            mult = 1;
            for k=1:num_facs
                if j~=k % the remaining factors
                    mult = mult*facs(k);
                end
            end
            count = count + 1;
            out_tmp(count,:) = [facs(j), mult];
         end
        % C_4_2
        for j=1:num_facs-1
            for k=j+1:num_facs
                mult1 = facs(j)*facs(k);
                mult = 1;
                for t=1:num_facs
                    if t~=j && t~=k
                        mult = mult*facs(t);
                    end
                end
                count = count + 1;
                out_tmp(count,:) = [mult1, mult];
            end
        end
        
    elseif num_facs==5
        num_tmp = 5*4*3*2*1/(1*4*3*2*1) + 5*4*3*2*1/(2*1*3*2*1);
        out_tmp = zeros(num_tmp, 2);
        out_tmp(:,1) = 0; out_tmp(:,2) = 1e10;
        count = 0;
        % C_5_1
         for j=1:num_facs
            mult = 1;
            for k=1:num_facs
                if j~=k % the remaining factors
                    mult = mult*facs(k);
                end
            end
            count = count + 1;
            out_tmp(count,:) = [facs(j), mult];
         end
        % C_5_2
        for j=1:num_facs-1
            for k=j+1:num_facs
                mult1 = facs(j)*facs(k);
                mult = 1;
                for t=1:num_facs
                    if t~=j && t~=k
                        mult = mult*facs(t);
                    end
                end
                count = count + 1;
                out_tmp(count,:) = [mult1, mult];
            end
        end    
        
    elseif num_facs==6
        num_tmp = 6*5*4*3*2*1/(1*5*4*3*2*1) + 6*5*4*3*2*1/(2*1*4*3*2*1) + 6*5*4*3*2*1/(3*2*1*3*2*1);
        out_tmp = zeros(num_tmp, 2);
        out_tmp(:,1) = 0; out_tmp(:,2) = 1e10;
        count = 0;
        % C_6_1
         for j=1:num_facs
            mult = 1;
            for k=1:num_facs
                if j~=k % the remaining factors
                    mult = mult*facs(k);
                end
            end
            count = count + 1;
            out_tmp(count,:) = [facs(j), mult];
         end
        % C_6_2
        for j=1:num_facs-1
            for k=j+1:num_facs
                mult1 = facs(j)*facs(k);
                mult = 1;
                for t=1:num_facs
                    if t~=j && t~=k
                        mult = mult*facs(t);
                    end
                end
                count = count + 1;
                out_tmp(count,:) = [mult1, mult];
            end
        end
        % C_6_3
        for j=1:num_facs-2
            for k=j+1:num_facs-1
                for q=k+1:num_facs
                    mult1 = facs(j)*facs(k)*facs(q);
                    mult = 1;
                    for t=1:num_facs
                        if t~=j && t~=k && t~=q
                            mult = mult*facs(t);
                        end
                    end
                    count = count + 1;
                    out_tmp(count,:) = [mult1, mult];
                end
            end
        end

    elseif num_facs==7
        num_tmp = 7*6*5*4*3*2*1/(1*6*5*4*3*2*1) + 7*6*5*4*3*2*1/(2*1*5*4*3*2*1) +...
                            7*6*5*4*3*2*1/(3*2*1*4*3*2*1);
        out_tmp = zeros(num_tmp, 2);
        out_tmp(:,1) = 0; out_tmp(:,2) = 1e10;
        count = 0;
        % C_7_1
         for j=1:num_facs
            mult = 1;
            for k=1:num_facs
                if j~=k % the remaining factors
                    mult = mult*facs(k);
                end
            end
            count = count + 1;
            out_tmp(count,:) = [facs(j), mult];
         end
        % C_7_2
        for j=1:num_facs-1
            for k=j+1:num_facs
                mult1 = facs(j)*facs(k);
                mult = 1;
                for t=1:num_facs
                    if t~=j && t~=k
                        mult = mult*facs(t);
                    end
                end
                count = count + 1;
                out_tmp(count,:) = [mult1, mult];
            end
        end
        % C_7_3
        for j=1:num_facs-2
            for k=j+1:num_facs-1
                for q=k+1:num_facs
                    mult1 = facs(j)*facs(k)*facs(q);
                    mult = 1;
                    for t=1:num_facs
                        if t~=j && t~=k && t~=q
                            mult = mult*facs(t);
                        end
                    end
                    count = count + 1;
                    out_tmp(count,:) = [mult1, mult];
                end
            end
        end
        
    elseif num_facs==8
        num_tmp = 8*7*6*5*4*3*2*1/(1*7*6*5*4*3*2*1) + 8*7*6*5*4*3*2*1/(2*1*6*5*4*3*2*1) + ...
                                8*7*6*5*4*3*2*1/(3*2*1*5*4*3*2*1) + 8*7*6*5*4*3*2*1/(4*3*2*1*4*3*2*1);
        out_tmp = zeros(num_tmp, 2);
        out_tmp(:,1) = 0; out_tmp(:,2) = 1e10;
        count = 0;
        % C_8_1
         for j=1:num_facs
            mult = 1;
            for k=1:num_facs
                if j~=k % the remaining factors
                    mult = mult*facs(k);
                end
            end
            count = count + 1;
            out_tmp(count,:) = [facs(j), mult];
         end
        % C_8_2
        for j=1:num_facs-1
            for k=j+1:num_facs
                mult1 = facs(j)*facs(k);
                mult = 1;
                for t=1:num_facs
                    if t~=j && t~=k
                        mult = mult*facs(t);
                    end
                end
                count = count + 1;
                out_tmp(count,:) = [mult1, mult];
            end
        end
        % C_8_3
        for j=1:num_facs-2
            for k=j+1:num_facs-1
                for q=k+1:num_facs
                    mult1 = facs(j)*facs(k)*facs(q);
                    mult = 1;
                    for t=1:num_facs
                        if t~=j && t~=k && t~=q
                            mult = mult*facs(t);
                        end
                    end
                    count = count + 1;
                    out_tmp(count,:) = [mult1, mult];
                end
            end
        end
        % C_8_4
        for j=1:num_facs-3
            for k=j+1:num_facs-2
                for q=k+1:num_facs-1
                    for w=q+1:num_facs
                        mult1 = facs(j)*facs(k)*facs(q)*facs(w);
                        mult = 1;
                        for t=1:num_facs
                            if t~=j && t~=k && t~=q && t~=w
                                mult = mult*facs(t);
                            end
                        end
                        count = count + 1;
                        out_tmp(count,:) = [mult1, mult];
                    end
                end
            end
        end
    
    elseif num_facs==9
        num_tmp = 9*8*7*6*5*4*3*2*1/(1*8*7*6*5*4*3*2*1) + 9*8*7*6*5*4*3*2*1/(2*1*7*6*5*4*3*2*1) + ...
                                9*8*7*6*5*4*3*2*1/(3*2*1*6*5*4*3*2*1) + 9*8*7*6*5*4*3*2*1/(4*3*2*1*5*4*3*2*1);
        out_tmp = zeros(num_tmp, 2);
        out_tmp(:,1) = 0; out_tmp(:,2) = 1e10;
        count = 0;
        % C_9_1
         for j=1:num_facs
            mult = 1;
            for k=1:num_facs
                if j~=k % the remaining factors
                    mult = mult*facs(k);
                end
            end
            count = count + 1;
            out_tmp(count,:) = [facs(j), mult];
         end
        % C_9_2
        for j=1:num_facs-1
            for k=j+1:num_facs
                mult1 = facs(j)*facs(k);
                mult = 1;
                for t=1:num_facs
                    if t~=j && t~=k
                        mult = mult*facs(t);
                    end
                end
                count = count + 1;
                out_tmp(count,:) = [mult1, mult];
            end
        end
        % C_9_3
        for j=1:num_facs-2
            for k=j+1:num_facs-1
                for q=k+1:num_facs
                    mult1 = facs(j)*facs(k)*facs(q);
                    mult = 1;
                    for t=1:num_facs
                        if t~=j && t~=k && t~=q
                            mult = mult*facs(t);
                        end
                    end
                    count = count + 1;
                    out_tmp(count,:) = [mult1, mult];
                end
            end
        end
        % C_9_4
        for j=1:num_facs-3
            for k=j+1:num_facs-2
                for q=k+1:num_facs-1
                    for w=q+1:num_facs
                        mult1 = facs(j)*facs(k)*facs(q)*facs(w);
                        mult = 1;
                        for t=1:num_facs
                            if t~=j && t~=k && t~=q && t~=w
                                mult = mult*facs(t);
                            end
                        end
                        count = count + 1;
                        out_tmp(count,:) = [mult1, mult];
                    end
                end
            end
        end
    
    elseif num_facs==10
        num_tmp = 10*9*8*7*6*5*4*3*2*1/(1*9*8*7*6*5*4*3*2*1) + 10*9*8*7*6*5*4*3*2*1/(2*1*8*7*6*5*4*3*2*1) + ...
                           10*9*8*7*6*5*4*3*2*1/(3*2*1*7*6*5*4*3*2*1) + 10*9*8*7*6*5*4*3*2*1/(4*3*2*1*6*5*4*3*2*1) + ...
                           10*9*8*7*6*5*4*3*2*1/(5*4*3*2*1*5*4*3*2*1);
        out_tmp = zeros(num_tmp, 2);
        out_tmp(:,1) = 0; out_tmp(:,2) = 1e10;
        count = 0;
        % C_10_1
         for j=1:num_facs
            mult = 1;
            for k=1:num_facs
                if j~=k % the remaining factors
                    mult = mult*facs(k);
                end
            end
            count = count + 1;
            out_tmp(count,:) = [facs(j), mult];
         end
        % C_10_2
        for j=1:num_facs-1
            for k=j+1:num_facs
                mult1 = facs(j)*facs(k);
                mult = 1;
                for t=1:num_facs
                    if t~=j && t~=k
                        mult = mult*facs(t);
                    end
                end
                count = count + 1;
                out_tmp(count,:) = [mult1, mult];
            end
        end
        % C_10_3
        for j=1:num_facs-2
            for k=j+1:num_facs-1
                for q=k+1:num_facs
                    mult1 = facs(j)*facs(k)*facs(q);
                    mult = 1;
                    for t=1:num_facs
                        if t~=j && t~=k && t~=q
                            mult = mult*facs(t);
                        end
                    end
                    count = count + 1;
                    out_tmp(count,:) = [mult1, mult];
                end
            end
        end
        % C_10_4
        for j=1:num_facs-3
            for k=j+1:num_facs-2
                for q=k+1:num_facs-1
                    for w=q+1:num_facs
                        mult1 = facs(j)*facs(k)*facs(q)*facs(w);
                        mult = 1;
                        for t=1:num_facs
                            if t~=j && t~=k && t~=q && t~=w
                                mult = mult*facs(t);
                            end
                        end
                        count = count + 1;
                        out_tmp(count,:) = [mult1, mult];
                    end
                end
            end
        end
        % C_10_5
        for j=1:num_facs-4
            for k=j+1:num_facs-3
                for q=k+1:num_facs-2
                    for w=q+1:num_facs-1
                        for p=w+1:num_facs
                            mult1 = facs(j)*facs(k)*facs(q)*facs(w)*facs(p);
                            mult = 1;
                            for t=1:num_facs
                                if t~=j && t~=k && t~=q && t~=w && t~=p
                                    mult = mult*facs(t);
                                end
                            end
                            count = count + 1;
                            out_tmp(count,:) = [mult1, mult];
                        end
                    end
                end
            end
        end
        
    elseif num_facs==11
        num_tmp = 11*10*9*8*7*6*5*4*3*2*1/(1*10*9*8*7*6*5*4*3*2*1) +...
                           11*10*9*8*7*6*5*4*3*2*1/(2*1*9*8*7*6*5*4*3*2*1) + ...
                           11*10*9*8*7*6*5*4*3*2*1/(3*2*1*8*7*6*5*4*3*2*1) + ...
                           11*10*9*8*7*6*5*4*3*2*1/(4*3*2*1*7*6*5*4*3*2*1) + ...
                           11*10*9*8*7*6*5*4*3*2*1/(5*4*3*2*1*6*5*4*3*2*1);
        out_tmp = zeros(num_tmp, 2);
        out_tmp(:,1) = 0; out_tmp(:,2) = 1e10;
        count = 0;
        % C_11_1
         for j=1:num_facs
            mult = 1;
            for k=1:num_facs
                if j~=k % the remaining factors
                    mult = mult*facs(k);
                end
            end
            count = count + 1;
            out_tmp(count,:) = [facs(j), mult];
         end
        % C_11_2
        for j=1:num_facs-1
            for k=j+1:num_facs
                mult1 = facs(j)*facs(k);
                mult = 1;
                for t=1:num_facs
                    if t~=j && t~=k
                        mult = mult*facs(t);
                    end
                end
                count = count + 1;
                out_tmp(count,:) = [mult1, mult];
            end
        end
        % C_11_3
        for j=1:num_facs-2
            for k=j+1:num_facs-1
                for q=k+1:num_facs
                    mult1 = facs(j)*facs(k)*facs(q);
                    mult = 1;
                    for t=1:num_facs
                        if t~=j && t~=k && t~=q
                            mult = mult*facs(t);
                        end
                    end
                    count = count + 1;
                    out_tmp(count,:) = [mult1, mult];
                end
            end
        end
        % C_11_4
        for j=1:num_facs-3
            for k=j+1:num_facs-2
                for q=k+1:num_facs-1
                    for w=q+1:num_facs
                        mult1 = facs(j)*facs(k)*facs(q)*facs(w);
                        mult = 1;
                        for t=1:num_facs
                            if t~=j && t~=k && t~=q && t~=w
                                mult = mult*facs(t);
                            end
                        end
                        count = count + 1;
                        out_tmp(count,:) = [mult1, mult];
                    end
                end
            end
        end
        % C_11_5
        for j=1:num_facs-4
            for k=j+1:num_facs-3
                for q=k+1:num_facs-2
                    for w=q+1:num_facs-1
                        for p=w+1:num_facs
                            mult1 = facs(j)*facs(k)*facs(q)*facs(w)*facs(p);
                            mult = 1;
                            for t=1:num_facs
                                if t~=j && t~=k && t~=q && t~=w && t~=p
                                    mult = mult*facs(t);
                                end
                            end
                            count = count + 1;
                            out_tmp(count,:) = [mult1, mult];
                        end
                    end
                end
            end
        end
        
    elseif num_facs==12
        num_tmp = 12*11*10*9*8*7*6*5*4*3*2*1/(1*11*10*9*8*7*6*5*4*3*2*1) + ...
                           12*11*10*9*8*7*6*5*4*3*2*1/(2*1*10*9*8*7*6*5*4*3*2*1) + ...
                           12*11*10*9*8*7*6*5*4*3*2*1/(3*2*1*9*8*7*6*5*4*3*2*1) + ...
                           12*11*10*9*8*7*6*5*4*3*2*1/(4*3*2*1*8*7*6*5*4*3*2*1) + ...
                           12*11*10*9*8*7*6*5*4*3*2*1/(5*4*3*2*1*7*6*5*4*3*2*1) +...
                           12*11*10*9*8*7*6*5*4*3*2*1/(6*5*4*3*2*1*6*5*4*3*2*1);
        out_tmp = zeros(num_tmp, 2);
        out_tmp(:,1) = 0; out_tmp(:,2) = 1e10;
        count = 0;
        % C_12_1
         for j=1:num_facs
            mult = 1;
            for k=1:num_facs
                if j~=k % the remaining factors
                    mult = mult*facs(k);
                end
            end
            count = count + 1;
            out_tmp(count,:) = [facs(j), mult];
         end
        % C_12_2
        for j=1:num_facs-1
            for k=j+1:num_facs
                mult1 = facs(j)*facs(k);
                mult = 1;
                for t=1:num_facs
                    if t~=j && t~=k
                        mult = mult*facs(t);
                    end
                end
                count = count + 1;
                out_tmp(count,:) = [mult1, mult];
            end
        end
        % C_12_3
        for j=1:num_facs-2
            for k=j+1:num_facs-1
                for q=k+1:num_facs
                    mult1 = facs(j)*facs(k)*facs(q);
                    mult = 1;
                    for t=1:num_facs
                        if t~=j && t~=k && t~=q
                            mult = mult*facs(t);
                        end
                    end
                    count = count + 1;
                    out_tmp(count,:) = [mult1, mult];
                end
            end
        end
        % C_12_4
        for j=1:num_facs-3
            for k=j+1:num_facs-2
                for q=k+1:num_facs-1
                    for w=q+1:num_facs
                        mult1 = facs(j)*facs(k)*facs(q)*facs(w);
                        mult = 1;
                        for t=1:num_facs
                            if t~=j && t~=k && t~=q && t~=w
                                mult = mult*facs(t);
                            end
                        end
                        count = count + 1;
                        out_tmp(count,:) = [mult1, mult];
                    end
                end
            end
        end
        % C_12_5
        for j=1:num_facs-4
            for k=j+1:num_facs-3
                for q=k+1:num_facs-2
                    for w=q+1:num_facs-1
                        for p=w+1:num_facs
                            mult1 = facs(j)*facs(k)*facs(q)*facs(w)*facs(p);
                            mult = 1;
                            for t=1:num_facs
                                if t~=j && t~=k && t~=q && t~=w && t~=p
                                    mult = mult*facs(t);
                                end
                            end
                            count = count + 1;
                            out_tmp(count,:) = [mult1, mult];
                        end
                    end
                end
            end
        end
        % C_12_6
        for j=1:num_facs-5
            for k=j+1:num_facs-4
                for q=k+1:num_facs-3
                    for w=q+1:num_facs-2
                        for p=w+1:num_facs-1
                            for u=p+1:num_facs
                                mult1 = facs(j)*facs(k)*facs(q)*facs(w)*facs(p)*facs(u);
                                mult = 1;
                                for t=1:num_facs
                                    if t~=j && t~=k && t~=q && t~=w && t~=p && t~=u
                                        mult = mult*facs(t);
                                    end
                                end
                                count = count + 1;
                                out_tmp(count,:) = [mult1, mult];
                            end
                        end
                    end
                end
            end
        end
        
    elseif num_facs==13
        num_tmp = 13*12*11*10*9*8*7*6*5*4*3*2*1/(1*12*11*10*9*8*7*6*5*4*3*2*1) + ...
                           13*12*11*10*9*8*7*6*5*4*3*2*1/(2*1*11*10*9*8*7*6*5*4*3*2*1) + ...
                           13*12*11*10*9*8*7*6*5*4*3*2*1/(3*2*1*10*9*8*7*6*5*4*3*2*1) + ...
                           13*12*11*10*9*8*7*6*5*4*3*2*1/(4*3*2*1*9*8*7*6*5*4*3*2*1) + ...
                           13*12*11*10*9*8*7*6*5*4*3*2*1/(5*4*3*2*1*8*7*6*5*4*3*2*1) +...
                           13*12*11*10*9*8*7*6*5*4*3*2*1/(6*5*4*3*2*1*7*6*5*4*3*2*1);
        out_tmp = zeros(num_tmp, 2);
        out_tmp(:,1) = 0; out_tmp(:,2) = 1e10;
        count = 0;
        % C_13_1
         for j=1:num_facs
            mult = 1;
            for k=1:num_facs
                if j~=k % the remaining factors
                    mult = mult*facs(k);
                end
            end
            count = count + 1;
            out_tmp(count,:) = [facs(j), mult];
         end
        % C_13_2
        for j=1:num_facs-1
            for k=j+1:num_facs
                mult1 = facs(j)*facs(k);
                mult = 1;
                for t=1:num_facs
                    if t~=j && t~=k
                        mult = mult*facs(t);
                    end
                end
                count = count + 1;
                out_tmp(count,:) = [mult1, mult];
            end
        end
        % C_13_3
        for j=1:num_facs-2
            for k=j+1:num_facs-1
                for q=k+1:num_facs
                    mult1 = facs(j)*facs(k)*facs(q);
                    mult = 1;
                    for t=1:num_facs
                        if t~=j && t~=k && t~=q
                            mult = mult*facs(t);
                        end
                    end
                    count = count + 1;
                    out_tmp(count,:) = [mult1, mult];
                end
            end
        end
        % C_13_4
        for j=1:num_facs-3
            for k=j+1:num_facs-2
                for q=k+1:num_facs-1
                    for w=q+1:num_facs
                        mult1 = facs(j)*facs(k)*facs(q)*facs(w);
                        mult = 1;
                        for t=1:num_facs
                            if t~=j && t~=k && t~=q && t~=w
                                mult = mult*facs(t);
                            end
                        end
                        count = count + 1;
                        out_tmp(count,:) = [mult1, mult];
                    end
                end
            end
        end
        % C_13_5
        for j=1:num_facs-4
            for k=j+1:num_facs-3
                for q=k+1:num_facs-2
                    for w=q+1:num_facs-1
                        for p=w+1:num_facs
                            mult1 = facs(j)*facs(k)*facs(q)*facs(w)*facs(p);
                            mult = 1;
                            for t=1:num_facs
                                if t~=j && t~=k && t~=q && t~=w && t~=p
                                    mult = mult*facs(t);
                                end
                            end
                            count = count + 1;
                            out_tmp(count,:) = [mult1, mult];
                        end
                    end
                end
            end
        end
        % C_13_6
        for j=1:num_facs-5
            for k=j+1:num_facs-4
                for q=k+1:num_facs-3
                    for w=q+1:num_facs-2
                        for p=w+1:num_facs-1
                            for u=p+1:num_facs
                                mult1 = facs(j)*facs(k)*facs(q)*facs(w)*facs(p)*facs(u);
                                mult = 1;
                                for t=1:num_facs
                                    if t~=j && t~=k && t~=q && t~=w && t~=p && t~=u
                                        mult = mult*facs(t);
                                    end
                                end
                                count = count + 1;
                                out_tmp(count,:) = [mult1, mult];
                            end
                        end
                    end
                end
            end
        end
        
    elseif num_facs==14
        num_tmp = 14*13*12*11*10*9*8*7*6*5*4*3*2*1/(1*13*12*11*10*9*8*7*6*5*4*3*2*1) + ...
                           14*13*12*11*10*9*8*7*6*5*4*3*2*1/(2*1*12*11*10*9*8*7*6*5*4*3*2*1) + ...
                           14*13*12*11*10*9*8*7*6*5*4*3*2*1/(3*2*1*11*10*9*8*7*6*5*4*3*2*1) + ...
                           14*13*12*11*10*9*8*7*6*5*4*3*2*1/(4*3*2*1*10*9*8*7*6*5*4*3*2*1) + ...
                           14*13*12*11*10*9*8*7*6*5*4*3*2*1/(5*4*3*2*1*9*8*7*6*5*4*3*2*1) +...
                           14*13*12*11*10*9*8*7*6*5*4*3*2*1/(6*5*4*3*2*1*8*7*6*5*4*3*2*1) +...
                           14*13*12*11*10*9*8*7*6*5*4*3*2*1/(7*6*5*4*3*2*1*7*6*5*4*3*2*1);
        out_tmp = zeros(num_tmp, 2);
        out_tmp(:,1) = 0; out_tmp(:,2) = 1e10;
        count = 0;
        % C_14_1
         for j=1:num_facs
            mult = 1;
            for k=1:num_facs
                if j~=k % the remaining factors
                    mult = mult*facs(k);
                end
            end
            count = count + 1;
            out_tmp(count,:) = [facs(j), mult];
         end
        % C_14_2
        for j=1:num_facs-1
            for k=j+1:num_facs
                mult1 = facs(j)*facs(k);
                mult = 1;
                for t=1:num_facs
                    if t~=j && t~=k
                        mult = mult*facs(t);
                    end
                end
                count = count + 1;
                out_tmp(count,:) = [mult1, mult];
            end
        end
        % C_14_3
        for j=1:num_facs-2
            for k=j+1:num_facs-1
                for q=k+1:num_facs
                    mult1 = facs(j)*facs(k)*facs(q);
                    mult = 1;
                    for t=1:num_facs
                        if t~=j && t~=k && t~=q
                            mult = mult*facs(t);
                        end
                    end
                    count = count + 1;
                    out_tmp(count,:) = [mult1, mult];
                end
            end
        end
        % C_14_4
        for j=1:num_facs-3
            for k=j+1:num_facs-2
                for q=k+1:num_facs-1
                    for w=q+1:num_facs
                        mult1 = facs(j)*facs(k)*facs(q)*facs(w);
                        mult = 1;
                        for t=1:num_facs
                            if t~=j && t~=k && t~=q && t~=w
                                mult = mult*facs(t);
                            end
                        end
                        count = count + 1;
                        out_tmp(count,:) = [mult1, mult];
                    end
                end
            end
        end
        % C_14_5
        for j=1:num_facs-4
            for k=j+1:num_facs-3
                for q=k+1:num_facs-2
                    for w=q+1:num_facs-1
                        for p=w+1:num_facs
                            mult1 = facs(j)*facs(k)*facs(q)*facs(w)*facs(p);
                            mult = 1;
                            for t=1:num_facs
                                if t~=j && t~=k && t~=q && t~=w && t~=p
                                    mult = mult*facs(t);
                                end
                            end
                            count = count + 1;
                            out_tmp(count,:) = [mult1, mult];
                        end
                    end
                end
            end
        end
        % C_14_6
        for j=1:num_facs-5
            for k=j+1:num_facs-4
                for q=k+1:num_facs-3
                    for w=q+1:num_facs-2
                        for p=w+1:num_facs-1
                            for u=p+1:num_facs
                                mult1 = facs(j)*facs(k)*facs(q)*facs(w)*facs(p)*facs(u);
                                mult = 1;
                                for t=1:num_facs
                                    if t~=j && t~=k && t~=q && t~=w && t~=p && t~=u
                                        mult = mult*facs(t);
                                    end
                                end
                                count = count + 1;
                                out_tmp(count,:) = [mult1, mult];
                            end
                        end
                    end
                end
            end
        end
        % C_14_7
        for j=1:num_facs-6
            for k=j+1:num_facs-5
                for q=k+1:num_facs-4
                    for w=q+1:num_facs-3
                        for p=w+1:num_facs-2
                            for u=p+1:num_facs-1
                                for b=u+1:num_facs
                                    mult1 = facs(j)*facs(k)*facs(q)*facs(w)*facs(p)*facs(u)*facs(b);
                                    mult = 1;
                                    for t=1:num_facs
                                        if t~=j && t~=k && t~=q && t~=w && t~=p && t~=u && t~=b
                                            mult = mult*facs(t);
                                        end
                                    end
                                    count = count + 1;
                                    out_tmp(count,:) = [mult1, mult];
                                end
                            end
                        end
                    end
                end
            end
        end
        
    elseif num_facs==15
        num_tmp = 15*14*13*12*11*10*9*8*7*6*5*4*3*2*1/(1*14*13*12*11*10*9*8*7*6*5*4*3*2*1) + ...
                           15*14*13*12*11*10*9*8*7*6*5*4*3*2*1/(2*1*13*12*11*10*9*8*7*6*5*4*3*2*1) + ...
                           15*14*13*12*11*10*9*8*7*6*5*4*3*2*1/(3*2*1*12*11*10*9*8*7*6*5*4*3*2*1) + ...
                           15*14*13*12*11*10*9*8*7*6*5*4*3*2*1/(4*3*2*1*11*10*9*8*7*6*5*4*3*2*1) + ...
                           15*14*13*12*11*10*9*8*7*6*5*4*3*2*1/(5*4*3*2*1*10*9*8*7*6*5*4*3*2*1) +...
                           15*14*13*12*11*10*9*8*7*6*5*4*3*2*1/(6*5*4*3*2*1*9*8*7*6*5*4*3*2*1) +...
                           15*14*13*12*11*10*9*8*7*6*5*4*3*2*1/(7*6*5*4*3*2*1*8*7*6*5*4*3*2*1);
        out_tmp = zeros(num_tmp, 2);
        out_tmp(:,1) = 0; out_tmp(:,2) = 1e10;
        count = 0;
        % C_15_1
         for j=1:num_facs
            mult = 1;
            for k=1:num_facs
                if j~=k % the remaining factors
                    mult = mult*facs(k);
                end
            end
            count = count + 1;
            out_tmp(count,:) = [facs(j), mult];
         end
        % C_15_2
        for j=1:num_facs-1
            for k=j+1:num_facs
                mult1 = facs(j)*facs(k);
                mult = 1;
                for t=1:num_facs
                    if t~=j && t~=k
                        mult = mult*facs(t);
                    end
                end
                count = count + 1;
                out_tmp(count,:) = [mult1, mult];
            end
        end
        % C_15_3
        for j=1:num_facs-2
            for k=j+1:num_facs-1
                for q=k+1:num_facs
                    mult1 = facs(j)*facs(k)*facs(q);
                    mult = 1;
                    for t=1:num_facs
                        if t~=j && t~=k && t~=q
                            mult = mult*facs(t);
                        end
                    end
                    count = count + 1;
                    out_tmp(count,:) = [mult1, mult];
                end
            end
        end
        % C_15_4
        for j=1:num_facs-3
            for k=j+1:num_facs-2
                for q=k+1:num_facs-1
                    for w=q+1:num_facs
                        mult1 = facs(j)*facs(k)*facs(q)*facs(w);
                        mult = 1;
                        for t=1:num_facs
                            if t~=j && t~=k && t~=q && t~=w
                                mult = mult*facs(t);
                            end
                        end
                        count = count + 1;
                        out_tmp(count,:) = [mult1, mult];
                    end
                end
            end
        end
        % C_15_5
        for j=1:num_facs-4
            for k=j+1:num_facs-3
                for q=k+1:num_facs-2
                    for w=q+1:num_facs-1
                        for p=w+1:num_facs
                            mult1 = facs(j)*facs(k)*facs(q)*facs(w)*facs(p);
                            mult = 1;
                            for t=1:num_facs
                                if t~=j && t~=k && t~=q && t~=w && t~=p
                                    mult = mult*facs(t);
                                end
                            end
                            count = count + 1;
                            out_tmp(count,:) = [mult1, mult];
                        end
                    end
                end
            end
        end
        % C_15_6
        for j=1:num_facs-5
            for k=j+1:num_facs-4
                for q=k+1:num_facs-3
                    for w=q+1:num_facs-2
                        for p=w+1:num_facs-1
                            for u=p+1:num_facs
                                mult1 = facs(j)*facs(k)*facs(q)*facs(w)*facs(p)*facs(u);
                                mult = 1;
                                for t=1:num_facs
                                    if t~=j && t~=k && t~=q && t~=w && t~=p && t~=u
                                        mult = mult*facs(t);
                                    end
                                end
                                count = count + 1;
                                out_tmp(count,:) = [mult1, mult];
                            end
                        end
                    end
                end
            end
        end
        % C_15_7
        for j=1:num_facs-6
            for k=j+1:num_facs-5
                for q=k+1:num_facs-4
                    for w=q+1:num_facs-3
                        for p=w+1:num_facs-2
                            for u=p+1:num_facs-1
                                for b=u+1:num_facs
                                    mult1 = facs(j)*facs(k)*facs(q)*facs(w)*facs(p)*facs(u)*facs(b);
                                    mult = 1;
                                    for t=1:num_facs
                                        if t~=j && t~=k && t~=q && t~=w && t~=p && t~=u && t~=b
                                            mult = mult*facs(t);
                                        end
                                    end
                                    count = count + 1;
                                    out_tmp(count,:) = [mult1, mult];
                                end
                            end
                        end
                    end
                end
            end
        end
        
    end
    
    
    % compare
    %num_nn = size(out_tmp,1);
    out_dis = abs( out_tmp(:,1)-out_tmp(:,2) );
    [error,idx] = sort(out_dis);
    tempv = out_tmp( idx(1), : );
    out(i,:) = [max(tempv), min(tempv)];
    
end

