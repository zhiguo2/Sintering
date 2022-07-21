clc;
clear all;
% str = ["110.0" "115.0" "120.0" "125.0" "130.0" "135.0" "140.0" "142.0" "144.0" "146.0" "148.0" "150.0"];
str = ["170.0"];
number_file = length(str);
dmsh = 0.001*0.001;
file_write_1 = strcat('por_sinter_boun_',num2str(dmsh),'_',str(1),'-',str(length(str)),'.txt');
fid_1 = fopen(file_write_1,'w');
for loopFile = 1:number_file
    filename1 = strcat('interactions',str(loopFile),'.txt');  % 内部球
    filename2 = strcat('boun_wall_',str(loopFile),'.txt');    % 与墙接触的球
    % 读入结构信息
    sp = load(filename1);
    roll = length(sp);
    sp2 = load(filename2);
    % 找出y轴范围
    y_up = zeros(2*roll,1);
    for i=1:2:2*roll
        y_up(i) = sp((i+1)/2,4) + sp((i+1)/2,6);
        y_up(i+1) = sp((i+1)/2,9) + sp((i+1)/2,11);
    end
    y_min = 0;         % wall_id = 2
    y_max = max(y_up); % wall_id = 3

    % 画网格
    % 真实结构尺寸
    x_min = 0;     % wall_id = 0
    x_max = 0.3*0.001;   % wall_id = 1
    z_min = 0;     % wall_id = 4
    z_max = 0.3*0.001;   % wall_id = 5

    % 网格尺寸，单位网格长度为dmsh
    LX = ceil((x_max - x_min) / dmsh);
    LY = floor((y_max - y_min) / dmsh);
    LZ = ceil((z_max - z_min) / dmsh);

    rr  = zeros(LX,LY,LZ);
    pos = zeros(LX,LY,LZ);
    %% 内部球
    for i=1:length(sp)
        x_w1 = sp(i,3) - sp(i,6) - x_min;
        x_w2 = sp(i,8) - sp(i,11) - x_min;
        while (x_w1 <= x_min)
            x_w1 = x_w1 + dmsh;
            if (x_w1 <= x_min)
                continue
            else
                break
            end
        end
         while (x_w2 <= x_min)
            x_w2 = x_w2 + dmsh;
            if (x_w2 <= x_min)
                continue
            else
                break
            end
        end
        x_e1 = sp(i,3) + sp(i,6) - x_min;
        x_e2 = sp(i,8) + sp(i,11) - x_min;
        while (x_e1 >= x_max)
            x_e1 = x_e1 - dmsh;
            if (x_e1 >= x_max)
                continue
            else
                break
            end
        end
        while (x_e2 >= x_max)
            x_e2 = x_e2 - dmsh;
            if (x_e2 >= x_max)
                continue
            else
                break
            end
        end

        y_s1 =sp(i,4) - sp(i,6) - y_min;
        y_s2 =sp(i,9) - sp(i,11) - y_min;
        while (y_s1 <= y_min)
            y_s1 = y_s1 + dmsh;
            if (y_s1 <= y_min)
                continue
            else
                break
            end
        end
        while (y_s2 <= y_min)
            y_s2 = y_s2 + dmsh;
            if (y_s2 <= y_min)
                continue
            else
                break
            end
        end
        y_n1 = sp(i,4) + sp(i,6) - y_min;
        y_n2 = sp(i,9) + sp(i,11) - y_min;
        while (y_n1 >= y_max)
            y_n1 = y_n1 - dmsh;
            if (y_n1 >= y_max)
                continue
            else
                break
            end
        end
        while (y_n2 >= y_max)
            y_n2 = y_n2 - dmsh;
            if (y_n2 >= y_max)
                continue
            else
                break
            end
        end

        z_b1 = sp(i,5) - sp(i,6) - z_min;
        z_b2 = sp(i,10) - sp(i,11) - z_min;
        while (z_b1 <= z_min)
            z_b1 = z_b1 + dmsh;
            if (z_b1 <= z_min)
                continue
            else
                break
            end
        end
        while (z_b2 <= z_min)
            z_b2 = z_b2 + dmsh;
            if (z_b2 <= z_min)
                continue
            else
                break
            end
        end
        z_f1=sp(i,5) + sp(i,6) - z_min;
        z_f2=sp(i,10) + sp(i,11) - z_min;
        while (z_f2 >= z_max)
            z_f2 = z_f2 - dmsh;
            if (z_f2 >= z_max)
                continue
            else
                break
            end
        end    

        node_w1 = ceil(x_w1 / dmsh);
        node_e1 = floor(x_e1 / dmsh);
        node_s1 = ceil(y_s1 / dmsh);
        node_n1 = floor(y_n1 / dmsh);
        node_b1 = ceil(z_b1 / dmsh);
        node_f1 = floor(z_f1 / dmsh);

        node_w2 = ceil(x_w2 / dmsh);
        node_e2 = floor(x_e2 / dmsh);
        node_s2 = ceil(y_s2 / dmsh);
        node_n2 = floor(y_n2 / dmsh);
        node_b2 = ceil(z_b2 / dmsh);
        node_f2 = floor(z_f2 / dmsh);

        pos_sinter_range_x = [];
        pos_sinter_range_y = [];
        pos_sinter_range_z = [];
        count_overlap = 0;
        for ii = min([node_w1 node_w2]) : max([node_e1 node_e2])
            for jj = min([node_s1 node_s2]) : max([node_n1 node_n2])
                for kk = min([node_b1 node_b2]) : max([node_f1 node_f2])
                    pos_x = dmsh * ii;
                    pos_y = dmsh * jj;
                    pos_z = dmsh * kk;
                    length_from_sp1_center = ((pos_x-sp(i,3))^2+(pos_y-sp(i,4))^2+(pos_z-sp(i,5))^2)^0.5;
                    length_from_sp2_center = ((pos_x-sp(i,8))^2+(pos_y-sp(i,9))^2+(pos_z-sp(i,10))^2)^0.5;
                    if (length_from_sp1_center <= sp(i,6))
                        rr(ii,jj,kk)=1;
                    end
                    if (length_from_sp2_center <= sp(i,11))
                        rr(ii,jj,kk)=1;
                    end
                    if (length_from_sp1_center <= sp(i,6)) && (length_from_sp2_center <= sp(i,11))
                        count_overlap = count_overlap + 1;
                        pos_sinter_range_x(count_overlap) = ii;
                        pos_sinter_range_y(count_overlap) = jj;
                        pos_sinter_range_z(count_overlap) = kk;
                    end
                end
            end
        end
        %% sintering: volume conservation
        expand = 1;
        count = 0;
        pos_sinter_por = [];
        while count_overlap > count
            count = 0;
            expand = expand + 2;
            pos_sinter_por = [];
            min_ii = min(pos_sinter_range_x) - expand;
            min_jj = min(pos_sinter_range_y) - expand;
            min_kk = min(pos_sinter_range_z) - expand;
            max_ii = max(pos_sinter_range_x) + expand;
            max_jj = max(pos_sinter_range_y) + expand;
            max_kk = max(pos_sinter_range_z) + expand;
            while min(pos_sinter_range_x) <= expand
                min_ii = min(pos_sinter_range_x);
                break
            end
            while min(pos_sinter_range_y) <= expand
                min_jj = min(pos_sinter_range_y);
                break
            end
            while min(pos_sinter_range_z) <= expand
                min_kk = min(pos_sinter_range_z);
                break
            end
            while max(pos_sinter_range_x) >= (LX - expand)
                max_ii = LX;
                break
            end
            while max(pos_sinter_range_y) >= (LY - expand)
                max_jj = LY;
                break
            end
            while max(pos_sinter_range_z) >= (LZ - expand)
                max_kk = LZ;
                break
            end
            for ii = min_ii:max_ii
                for jj = min_jj:max_jj
                    for kk = min_kk:max_kk
                        if rr(ii,jj,kk) == 0 && ii~=1 && jj~=1 && kk~=1 && ii~=LX && jj~=LY && kk~=LZ
                            count = count + 1;
                            pos_sinter_por(count,1) = ii;
                            pos_sinter_por(count,2) = jj;
                            pos_sinter_por(count,3) = kk;
                            pos_sinter_por(count,4) = rr(ii-1,jj-1,kk-1) + rr(ii-1,jj-1,kk) + rr(ii-1,jj-1,kk+1) +...
                                rr(ii-1,jj,kk-1) + rr(ii-1,jj,kk) + rr(ii-1,jj,kk+1) +...
                                rr(ii-1,jj+1,kk-1) + rr(ii-1,jj+1,kk) + rr(ii-1,jj+1,kk+1) +...
                                rr(ii,jj-1,kk-1) + rr(ii,jj-1,kk) + rr(ii,jj-1,kk+1) +...
                                rr(ii,jj,kk-1)  + rr(ii,jj,kk+1) +...
                                rr(ii,jj+1,kk-1) + rr(ii,jj+1,kk) + rr(ii,jj+1,kk+1) +...
                                rr(ii+1,jj-1,kk-1) + rr(ii+1,jj-1,kk) + rr(ii+1,jj-1,kk+1) +...
                                rr(ii+1,jj,kk-1) + rr(ii+1,jj,kk) + rr(ii+1,jj,kk+1) +...
                                rr(ii+1,jj+1,kk-1) + rr(ii+1,jj+1,kk) + rr(ii+1,jj+1,kk+1);
                        end
                    end
                end
            end
        end
        if ~isempty(pos_sinter_por)
            new_pos_sinter_por = sortrows(pos_sinter_por,4,'descend');
        end
        if count_overlap ~= 0
            for ii = 1:count_overlap
                rr(new_pos_sinter_por(ii,1), new_pos_sinter_por(ii,2), new_pos_sinter_por(ii,3)) = 1;
            end
        end
    end
    %% 与墙接触的球
    for i=1:length(sp2)
        if sp2(i,2) == 0  % wall_id = 0, x_min
            h = sp2(i,6) - sp2(i,3);
            vol_overlap = 1/3*pi * h^2 * ( 3*sp2(i,6) - h );
            count_overlap = floor(vol_overlap / dmsh^3);
            
            expand = 1;
            count = 0;
            pos_sinter_por = [];
            while count_overlap > count
                count = 0;
                expand = expand + 2;
                pos_sinter_por = [];
                min_ii = 1;
                max_ii = 1 + expand;
                if floor(sp2(i,4)/dmsh) > expand
                    min_jj = ceil(sp2(i,4)/dmsh) - expand;
                else
                    min_jj = ceil(sp2(i,4)/dmsh);
                end
                if (LY - ceil(sp2(i,4)/dmsh)) > expand
                    max_jj = ceil(sp2(i,4)/dmsh) + expand;
                else
                    max_jj = ceil(sp2(i,4)/dmsh);
                end
                if floor(sp2(i,5)/dmsh) > expand
                    min_kk = ceil(sp2(i,5)/dmsh) - expand;
                else
                    min_kk = ceil(sp2(i,5)/dmsh);
                end
                if (LZ - ceil(sp2(i,5)/dmsh)) > expand
                    max_kk = ceil(sp2(i,5)/dmsh) + expand;
                else
                    max_kk = ceil(sp2(i,5)/dmsh);
                end
                for ii = min_ii:max_ii
                    for jj = min_jj:max_jj
                        for kk = min_kk:max_kk
                            if rr(ii,jj,kk) == 0 && ii~=1 && jj~=1 && kk~=1 && ii~=LX && jj~=LY && kk~=LZ
                                count = count + 1;
                                pos_sinter_por(count,1) = ii;
                                pos_sinter_por(count,2) = jj;
                                pos_sinter_por(count,3) = kk;
                                pos_sinter_por(count,4) = rr(ii-1,jj-1,kk-1) + rr(ii-1,jj-1,kk) + rr(ii-1,jj-1,kk+1) +...
                                    rr(ii-1,jj,kk-1) + rr(ii-1,jj,kk) + rr(ii-1,jj,kk+1) +...
                                    rr(ii-1,jj+1,kk-1) + rr(ii-1,jj+1,kk) + rr(ii-1,jj+1,kk+1) +...
                                    rr(ii,jj-1,kk-1) + rr(ii,jj-1,kk) + rr(ii,jj-1,kk+1) +...
                                    rr(ii,jj,kk-1)  + rr(ii,jj,kk+1) +...
                                    rr(ii,jj+1,kk-1) + rr(ii,jj+1,kk) + rr(ii,jj+1,kk+1) +...
                                    rr(ii+1,jj-1,kk-1) + rr(ii+1,jj-1,kk) + rr(ii+1,jj-1,kk+1) +...
                                    rr(ii+1,jj,kk-1) + rr(ii+1,jj,kk) + rr(ii+1,jj,kk+1) +...
                                    rr(ii+1,jj+1,kk-1) + rr(ii+1,jj+1,kk) + rr(ii+1,jj+1,kk+1);
                            end
                        end
                    end
                end
            end
            if ~isempty(pos_sinter_por)
                new_pos_sinter_por = sortrows(pos_sinter_por,4,'descend');
            end
            if count_overlap ~= 0
                for ii = 1:count_overlap
                    rr(new_pos_sinter_por(ii,1), new_pos_sinter_por(ii,2), new_pos_sinter_por(ii,3)) = 1;
                end
            end
        end
        if sp2(i,2) == 1  % wall_id = 1, x_max
            h = sp2(i,6) - (x_max - sp2(i,3));
            vol_overlap = 1/3*pi * h^2 * ( 3*sp2(i,6) - h );
            count_overlap = floor(vol_overlap / dmsh^3);
            
            expand = 1;
            count = 0;
            pos_sinter_por = [];
            while count_overlap > count
                count = 0;
                expand = expand + 2;
                pos_sinter_por = [];
                min_ii = floor(sp2(i,3)/dmsh) - expand;
                max_ii = floor(sp2(i,3)/dmsh);
                if floor(sp2(i,4)/dmsh) > expand
                    min_jj = ceil(sp2(i,4)/dmsh) - expand;
                else
                    min_jj = ceil(sp2(i,4)/dmsh);
                end
                if (LY - ceil(sp2(i,4)/dmsh)) > expand
                    max_jj = ceil(sp2(i,4)/dmsh) + expand;
                else
                    max_jj = ceil(sp2(i,4)/dmsh);
                end
                if floor(sp2(i,5)/dmsh) > expand
                    min_kk = ceil(sp2(i,5)/dmsh) - expand;
                else
                    min_kk = ceil(sp2(i,5)/dmsh);
                end
                if (LZ - ceil(sp2(i,5)/dmsh)) > expand
                    max_kk = ceil(sp2(i,5)/dmsh) + expand;
                else
                    max_kk = ceil(sp2(i,5)/dmsh);
                end
                for ii = min_ii:max_ii
                    for jj = min_jj:max_jj
                        for kk = min_kk:max_kk
                            if rr(ii,jj,kk) == 0 && ii~=1 && jj~=1 && kk~=1 && ii~=LX && jj~=LY && kk~=LZ
                                count = count + 1;
                                pos_sinter_por(count,1) = ii;
                                pos_sinter_por(count,2) = jj;
                                pos_sinter_por(count,3) = kk;
                                pos_sinter_por(count,4) = rr(ii-1,jj-1,kk-1) + rr(ii-1,jj-1,kk) + rr(ii-1,jj-1,kk+1) +...
                                    rr(ii-1,jj,kk-1) + rr(ii-1,jj,kk) + rr(ii-1,jj,kk+1) +...
                                    rr(ii-1,jj+1,kk-1) + rr(ii-1,jj+1,kk) + rr(ii-1,jj+1,kk+1) +...
                                    rr(ii,jj-1,kk-1) + rr(ii,jj-1,kk) + rr(ii,jj-1,kk+1) +...
                                    rr(ii,jj,kk-1)  + rr(ii,jj,kk+1) +...
                                    rr(ii,jj+1,kk-1) + rr(ii,jj+1,kk) + rr(ii,jj+1,kk+1) +...
                                    rr(ii+1,jj-1,kk-1) + rr(ii+1,jj-1,kk) + rr(ii+1,jj-1,kk+1) +...
                                    rr(ii+1,jj,kk-1) + rr(ii+1,jj,kk) + rr(ii+1,jj,kk+1) +...
                                    rr(ii+1,jj+1,kk-1) + rr(ii+1,jj+1,kk) + rr(ii+1,jj+1,kk+1);
                            end
                        end
                    end
                end
            end
            if ~isempty(pos_sinter_por)
                new_pos_sinter_por = sortrows(pos_sinter_por,4,'descend');
            end
            if count_overlap ~= 0
                for ii = 1:count_overlap
                    rr(new_pos_sinter_por(ii,1), new_pos_sinter_por(ii,2), new_pos_sinter_por(ii,3)) = 1;
                end
            end
        end
        if sp2(i,2) == 2  % wall_id = 2, y_min
            h = sp2(i,6) - sp2(i,4);
            vol_overlap = 1/3*pi * h^2 * ( 3*sp2(i,6) - h );
            count_overlap = floor(vol_overlap / dmsh^3);
            
            expand = 1;
            count = 0;
            pos_sinter_por = [];
            while count_overlap > count
                count = 0;
                expand = expand + 2;
                pos_sinter_por = [];
                min_jj = 1;
                max_jj = 1 + expand;
                if floor(sp2(i,3)/dmsh) > expand
                    min_ii = ceil(sp2(i,3)/dmsh) - expand;
                else
                    min_ii = ceil(sp2(i,3)/dmsh);
                end
                if (LX - ceil(sp2(i,3)/dmsh)) > expand
                    max_ii = ceil(sp2(i,3)/dmsh) + expand;
                else
                    max_ii = ceil(sp2(i,3)/dmsh);
                end
                if floor(sp2(i,5)/dmsh) > expand
                    min_kk = ceil(sp2(i,5)/dmsh) - expand;
                else
                    min_kk = ceil(sp2(i,5)/dmsh);
                end
                if (LZ - ceil(sp2(i,5)/dmsh)) > expand
                    max_kk = ceil(sp2(i,5)/dmsh) + expand;
                else
                    max_kk = ceil(sp2(i,5)/dmsh);
                end
                for ii = min_ii:max_ii
                    for jj = min_jj:max_jj
                        for kk = min_kk:max_kk
                            if rr(ii,jj,kk) == 0 && ii~=1 && jj~=1 && kk~=1 && ii~=LX && jj~=LY && kk~=LZ
                                count = count + 1;
                                pos_sinter_por(count,1) = ii;
                                pos_sinter_por(count,2) = jj;
                                pos_sinter_por(count,3) = kk;
                                pos_sinter_por(count,4) = rr(ii-1,jj-1,kk-1) + rr(ii-1,jj-1,kk) + rr(ii-1,jj-1,kk+1) +...
                                    rr(ii-1,jj,kk-1) + rr(ii-1,jj,kk) + rr(ii-1,jj,kk+1) +...
                                    rr(ii-1,jj+1,kk-1) + rr(ii-1,jj+1,kk) + rr(ii-1,jj+1,kk+1) +...
                                    rr(ii,jj-1,kk-1) + rr(ii,jj-1,kk) + rr(ii,jj-1,kk+1) +...
                                    rr(ii,jj,kk-1)  + rr(ii,jj,kk+1) +...
                                    rr(ii,jj+1,kk-1) + rr(ii,jj+1,kk) + rr(ii,jj+1,kk+1) +...
                                    rr(ii+1,jj-1,kk-1) + rr(ii+1,jj-1,kk) + rr(ii+1,jj-1,kk+1) +...
                                    rr(ii+1,jj,kk-1) + rr(ii+1,jj,kk) + rr(ii+1,jj,kk+1) +...
                                    rr(ii+1,jj+1,kk-1) + rr(ii+1,jj+1,kk) + rr(ii+1,jj+1,kk+1);
                            end
                        end
                    end
                end
            end
            if ~isempty(pos_sinter_por)
                new_pos_sinter_por = sortrows(pos_sinter_por,4,'descend');
            end
            if count_overlap ~= 0
                for ii = 1:count_overlap
                    rr(new_pos_sinter_por(ii,1), new_pos_sinter_por(ii,2), new_pos_sinter_por(ii,3)) = 1;
                end
            end
        end
        if sp2(i,2) == 3  % wall_id = 3, y_max
            h = sp2(i,6) - (y_max - sp2(i,4));
            vol_overlap = 1/3*pi * h^2 * ( 3*sp2(i,6) - h );
            count_overlap = floor(vol_overlap / dmsh^3);
            
            expand = 1;
            count = 0;
            pos_sinter_por = [];
            while count_overlap > count
                count = 0;
                expand = expand + 2;
                pos_sinter_por = [];
                min_jj = floor(sp2(i,4)/dmsh) - expand;
                max_jj = floor(sp2(i,4)/dmsh);
                if floor(sp2(i,3)/dmsh) > expand
                    min_ii = ceil(sp2(i,3)/dmsh) - expand;
                else
                    min_ii = ceil(sp2(i,3)/dmsh);
                end
                if (LX - ceil(sp2(i,3)/dmsh)) > expand
                    max_ii = ceil(sp2(i,3)/dmsh) + expand;
                else
                    max_ii = ceil(sp2(i,3)/dmsh);
                end
                if floor(sp2(i,5)/dmsh) > expand
                    min_kk = ceil(sp2(i,5)/dmsh) - expand;
                else
                    min_kk = ceil(sp2(i,5)/dmsh);
                end
                if (LZ - ceil(sp2(i,5)/dmsh)) > expand
                    max_kk = ceil(sp2(i,5)/dmsh) + expand;
                else
                    max_kk = ceil(sp2(i,5)/dmsh);
                end
                for ii = min_ii:max_ii
                    for jj = min_jj:max_jj
                        for kk = min_kk:max_kk
                            if rr(ii,jj,kk) == 0 && ii~=1 && jj~=1 && kk~=1 && ii~=LX && jj~=LY && kk~=LZ
                                count = count + 1;
                                pos_sinter_por(count,1) = ii;
                                pos_sinter_por(count,2) = jj;
                                pos_sinter_por(count,3) = kk;
                                pos_sinter_por(count,4) = rr(ii-1,jj-1,kk-1) + rr(ii-1,jj-1,kk) + rr(ii-1,jj-1,kk+1) +...
                                    rr(ii-1,jj,kk-1) + rr(ii-1,jj,kk) + rr(ii-1,jj,kk+1) +...
                                    rr(ii-1,jj+1,kk-1) + rr(ii-1,jj+1,kk) + rr(ii-1,jj+1,kk+1) +...
                                    rr(ii,jj-1,kk-1) + rr(ii,jj-1,kk) + rr(ii,jj-1,kk+1) +...
                                    rr(ii,jj,kk-1)  + rr(ii,jj,kk+1) +...
                                    rr(ii,jj+1,kk-1) + rr(ii,jj+1,kk) + rr(ii,jj+1,kk+1) +...
                                    rr(ii+1,jj-1,kk-1) + rr(ii+1,jj-1,kk) + rr(ii+1,jj-1,kk+1) +...
                                    rr(ii+1,jj,kk-1) + rr(ii+1,jj,kk) + rr(ii+1,jj,kk+1) +...
                                    rr(ii+1,jj+1,kk-1) + rr(ii+1,jj+1,kk) + rr(ii+1,jj+1,kk+1);
                            end
                        end
                    end
                end
            end
            if ~isempty(pos_sinter_por)
                new_pos_sinter_por = sortrows(pos_sinter_por,4,'descend');
            end
            if count_overlap ~= 0
                for ii = 1:count_overlap
                    rr(new_pos_sinter_por(ii,1), new_pos_sinter_por(ii,2), new_pos_sinter_por(ii,3)) = 1;
                end
            end
        end
        if sp2(i,2) == 4  % wall_id = 4, z_min
            h = sp2(i,6) - sp2(i,5);
            vol_overlap = 1/3*pi * h^2 * ( 3*sp2(i,6) - h );
            count_overlap = floor(vol_overlap / dmsh^3);
            
            expand = 1;
            count = 0;
            pos_sinter_por = [];
            while count_overlap > count
                count = 0;
                expand = expand + 2;
                pos_sinter_por = [];
                min_KK = 1;
                max_KK = 1 + expand;
                if floor(sp2(i,3)/dmsh) > expand
                    min_ii = ceil(sp2(i,3)/dmsh) - expand;
                else
                    min_ii = ceil(sp2(i,3)/dmsh);
                end
                if (LX - ceil(sp2(i,3)/dmsh)) > expand
                    max_ii = ceil(sp2(i,3)/dmsh) + expand;
                else
                    max_ii = ceil(sp2(i,3)/dmsh);
                end
                if floor(sp2(i,4)/dmsh) > expand
                    min_jj = ceil(sp2(i,4)/dmsh) - expand;
                else
                    min_jj= ceil(sp2(i,4)/dmsh);
                end
                if (LY - ceil(sp2(i,4)/dmsh)) > expand
                    max_jj = ceil(sp2(i,4)/dmsh) + expand;
                else
                    max_jj = ceil(sp2(i,4)/dmsh);
                end
                for ii = min_ii:max_ii
                    for jj = min_jj:max_jj
                        for kk = min_kk:max_kk
                            if rr(ii,jj,kk) == 0 && ii~=1 && jj~=1 && kk~=1 && ii~=LX && jj~=LY && kk~=LZ
                                count = count + 1;
                                pos_sinter_por(count,1) = ii;
                                pos_sinter_por(count,2) = jj;
                                pos_sinter_por(count,3) = kk;
                                pos_sinter_por(count,4) = rr(ii-1,jj-1,kk-1) + rr(ii-1,jj-1,kk) + rr(ii-1,jj-1,kk+1) +...
                                    rr(ii-1,jj,kk-1) + rr(ii-1,jj,kk) + rr(ii-1,jj,kk+1) +...
                                    rr(ii-1,jj+1,kk-1) + rr(ii-1,jj+1,kk) + rr(ii-1,jj+1,kk+1) +...
                                    rr(ii,jj-1,kk-1) + rr(ii,jj-1,kk) + rr(ii,jj-1,kk+1) +...
                                    rr(ii,jj,kk-1)  + rr(ii,jj,kk+1) +...
                                    rr(ii,jj+1,kk-1) + rr(ii,jj+1,kk) + rr(ii,jj+1,kk+1) +...
                                    rr(ii+1,jj-1,kk-1) + rr(ii+1,jj-1,kk) + rr(ii+1,jj-1,kk+1) +...
                                    rr(ii+1,jj,kk-1) + rr(ii+1,jj,kk) + rr(ii+1,jj,kk+1) +...
                                    rr(ii+1,jj+1,kk-1) + rr(ii+1,jj+1,kk) + rr(ii+1,jj+1,kk+1);
                            end
                        end
                    end
                end
            end
            if ~isempty(pos_sinter_por)
                new_pos_sinter_por = sortrows(pos_sinter_por,4,'descend');
            end
            if count_overlap ~= 0
                for ii = 1:count_overlap
                    rr(new_pos_sinter_por(ii,1), new_pos_sinter_por(ii,2), new_pos_sinter_por(ii,3)) = 1;
                end
            end
        end
        if sp2(i,2) == 5  % wall_id = 5, z_max
            h = sp2(i,6) - (x_max - sp2(i,5));
            vol_overlap = 1/3*pi * h^2 * ( 3*sp2(i,6) - h );
            count_overlap = floor(vol_overlap / dmsh^3);
            
            expand = 1;
            count = 0;
            pos_sinter_por = [];
            while count_overlap > count
                count = 0;
                expand = expand + 2;
                pos_sinter_por = [];
                min_kk = floor(sp2(i,5)/dmsh) - expand;
                max_kk = floor(sp2(i,5)/dmsh);
                if floor(sp2(i,3)/dmsh) > expand
                    min_ii = ceil(sp2(i,3)/dmsh) - expand;
                else
                    min_ii = ceil(sp2(i,3)/dmsh);
                end
                if (LX - ceil(sp2(i,3)/dmsh)) > expand
                    max_ii = ceil(sp2(i,3)/dmsh) + expand;
                else
                    max_ii = ceil(sp2(i,3)/dmsh);
                end
                if floor(sp2(i,4)/dmsh) > expand
                    min_jj = ceil(sp2(i,4)/dmsh) - expand;
                else
                    min_jj = ceil(sp2(i,4)/dmsh);
                end
                if (LY - ceil(sp2(i,4)/dmsh)) > expand
                    max_jj = ceil(sp2(i,4)/dmsh) + expand;
                else
                    max_jj = ceil(sp2(i,4)/dmsh);
                end
                for ii = min_ii:max_ii
                    for jj = min_jj:max_jj
                        for kk = min_kk:max_kk
                            if rr(ii,jj,kk) == 0 && ii~=1 && jj~=1 && kk~=1 && ii~=LX && jj~=LY && kk~=LZ
                                count = count + 1;
                                pos_sinter_por(count,1) = ii;
                                pos_sinter_por(count,2) = jj;
                                pos_sinter_por(count,3) = kk;
                                pos_sinter_por(count,4) = rr(ii-1,jj-1,kk-1) + rr(ii-1,jj-1,kk) + rr(ii-1,jj-1,kk+1) +...
                                    rr(ii-1,jj,kk-1) + rr(ii-1,jj,kk) + rr(ii-1,jj,kk+1) +...
                                    rr(ii-1,jj+1,kk-1) + rr(ii-1,jj+1,kk) + rr(ii-1,jj+1,kk+1) +...
                                    rr(ii,jj-1,kk-1) + rr(ii,jj-1,kk) + rr(ii,jj-1,kk+1) +...
                                    rr(ii,jj,kk-1)  + rr(ii,jj,kk+1) +...
                                    rr(ii,jj+1,kk-1) + rr(ii,jj+1,kk) + rr(ii,jj+1,kk+1) +...
                                    rr(ii+1,jj-1,kk-1) + rr(ii+1,jj-1,kk) + rr(ii+1,jj-1,kk+1) +...
                                    rr(ii+1,jj,kk-1) + rr(ii+1,jj,kk) + rr(ii+1,jj,kk+1) +...
                                    rr(ii+1,jj+1,kk-1) + rr(ii+1,jj+1,kk) + rr(ii+1,jj+1,kk+1);
                            end
                        end
                    end
                end
            end
            if ~isempty(pos_sinter_por)
                new_pos_sinter_por = sortrows(pos_sinter_por,4,'descend');
            end
            if count_overlap ~= 0
                for ii = 1:count_overlap
                    rr(new_pos_sinter_por(ii,1), new_pos_sinter_por(ii,2), new_pos_sinter_por(ii,3)) = 1;
                end
            end
        end
    end
    % 计算孔隙率
    for mm = 1:50
        sp_vol=0;
        tol_vol=0;
        tol_vol_1=0;
        surface_area=0;
        for ii=(mm+1):(LX-mm)
            for jj=1:LY
                for kk=(mm+1):(LZ-mm)
                    tol_vol=tol_vol+1;
                    if rr(ii,jj,kk)==1
                        sp_vol=sp_vol+1;
                        if ii~=1 && jj~=1 && kk~=1 && ii~=LX && jj~=LY && kk~=LZ
                            if rr(ii+1,jj,kk)==0 || rr(ii-1,jj,kk)==0 || rr(ii,jj-1,kk)==0 || rr(ii,jj+1,kk)==0 || rr(ii,jj,kk-1)==0 || rr(ii,jj,kk+1)==0 ...
                                surface_area=surface_area+1;
                            end
                        end
                    end
                end
            end
        end
        porosity(mm) = 1-sp_vol/tol_vol;
        surface_relative(mm) = surface_area/tol_vol;
        fprintf(fid_1,'%d %12.8f %12.8f\n',str2num(str(loopFile)), porosity(mm), surface_relative(mm));
    end
%     figure(1)
%     plot(porosity)
    %save 3d_structure rr
    % 写入Dat文件  
    boun = 12;
    file_write = strcat('geo_',str(loopFile),'_',num2str(LX-boun*2),'_',num2str(LY),'_',num2str(LZ-boun*2),'.dat');
    fid=fopen(file_write,'w');

    for i=(boun+1):(LX-boun)
        x_plane=rr(i,:,:);
        for j=1:LY
            for k=(boun+1):(LZ-boun)
                fprintf(fid,'%d ',x_plane(1,j,k));
            end

            fprintf(fid,'\r\n');
        end
    end

    fclose(fid);
end
fclose(fid_1);