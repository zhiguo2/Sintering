clc;
clear all;
% str = ["130.0"];
str = ["320.0" "340.0" "360.0" "380.0" "400.0"];
number_file = length(str);
dmsh = 0.001*0.001;
file_write_1 = strcat('por_no_sinter_boun_',num2str(dmsh),'_',str(1),'-',str(length(str)),'.txt');
fid_1 = fopen(file_write_1,'w');
for loopFile = 1:number_file
    filename = strcat('sphere_positions',str(loopFile),'.txt');
    % 读入结构信息
    sp=load(filename);
    roll=length(sp);
    % 找出x/y/z轴范围
    y_up = zeros(roll,1);
    for i=1:roll
        y_up(i) = sp(i,2) + sp(i,4);
    end
    y_min = 0;
    y_max = max(y_up);

    % 画网格
    % 真实结构尺寸
    x_min = 0;
    x_max = 0.3*0.001;
    z_min = 0;
    z_max = 0.3*0.001;
    % 网格尺寸，单位网格长度为dmsh
    LX=ceil((x_max-x_min)/dmsh);
    LY=ceil((y_max-y_min)/dmsh);
    LZ=ceil((z_max-z_min)/dmsh);

    rr=zeros(LX,LY,LZ);
    pos=zeros(LX,LY,LZ);

    for i=1:length(sp)
        x_w=sp(i,1)-sp(i,4)-x_min;
        while (x_w<=x_min)
            x_w=x_w+dmsh;
            if (x_w<=x_min)
                continue
            else
                break
            end
        end
        x_e=sp(i,1)+sp(i,4);
        while (x_e>=x_max)
            x_e=x_e-dmsh;
            if (x_e>=x_max)
                continue
            else
                break
            end
        end

        y_s=sp(i,2)-sp(i,4)-y_min;
        while (y_s<=y_min)
            y_s=y_s+dmsh;
            if (y_s<=y_min)
                continue
            else
                break
            end
        end
        y_n=sp(i,2)+sp(i,4)-y_min;
        while (y_n>=y_max)
            y_n=y_n-dmsh;
            if (y_n>=y_max)
                continue
            else
                break
            end
        end

        z_b=sp(i,3)-sp(i,4)-z_min;
        while (z_b<=z_min)
            z_b=z_b+dmsh;
            if (z_b<=z_min)
                continue
            else
                break
            end
        end
        z_f=sp(i,3)+sp(i,4)-z_min;
        while (z_f>=z_max)
            z_f=z_f-dmsh;
            if (z_f>=z_max)
                continue
            else
                break
            end
        end    

        node_w=ceil(x_w/dmsh);
        node_e=floor(x_e/dmsh);
        node_s=ceil(y_s/dmsh);
        node_n=floor(y_n/dmsh);
        node_b=ceil(z_b/dmsh);
        node_f=floor(z_f/dmsh);

        for ii=node_w:node_e
            for jj=node_s:node_n
                for kk=node_b:node_f
                    pos_x=dmsh*ii;
                    pos_y=dmsh*jj;
                    pos_z=dmsh*kk;
                    length_from_sp_center=((pos_x-sp(i,1))^2+(pos_y-sp(i,2)+y_min)^2+(pos_z-sp(i,3))^2)^0.5;
                    if (length_from_sp_center <= sp(i,4))
                        rr(ii,jj,kk)=1;
                    end
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

    % 写入Dat文件    
    boun = 6;
    file_write = strcat('geo_noSinter_',str(loopFile),'_',num2str(LX-boun*2),'_',num2str(LY),'_',num2str(LZ-boun*2),'.dat');
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
