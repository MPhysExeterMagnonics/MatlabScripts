% 
classdef FFTmaster < hgsetget % subclass hgsetget
   
properties

meshunit = 'm';
meshtype = 'rectangular';
xbase
ybase
zbase
xnodes
ynodes
znodes
xstepsize = 0.001;
ystepsize = 0.001;
zstepsize = 0.001;
xmin = 0;
ymin = 0;
zmin = 0;
xmax = 0.01;
ymax = 0.01;
zmax = 0.01;
dim = 3;

staticFile = 'static.stc';

% conversion factor for magnetization (A/m -> G)
convM = 1/1e3;
% conversion factor for magnetic field (T -> G)
%convH = 1e4; % muMax (T -> G)
convH = 4*pi/1e3; % OOMMF (A/m -> G)
    
% ground magnetization
M0;
% saturation magnetization
Ms;

header = '';   
format = '';

doubleTestVal = 123456789012345.0;
singleTestVal = 1234567.0;

end

methods
    
function run(obj)
  
currentpath=pwd;
path = input('Enter the path: ','s');
cd (path);

obj.M0 = obj.convM*obj.readFile(obj.staticFile); 
obj.Ms = sqrt(obj.M0(:,:,:,1).^2+obj.M0(:,:,:,2).^2+obj.M0(:,:,:,3).^2);

sq=zeros(size(obj.M0,1),size(obj.M0,2),size(obj.M0,3));

for x = 1:size(obj.M0,1)
    for y = 1:size(obj.M0,2)
        for z = 1:size(obj.M0,3)
            
            sq(x,y,z) = sqrt(obj.M0(x,y,z,1)^2+obj.M0(x,y,z,2)^2);
            
        end
    end
end



rot=zeros(size(obj.M0,1),size(obj.M0,2),size(obj.M0,3),3,3);

for x = 1:size(obj.M0,1)
    for y = 1:size(obj.M0,2)
        for z = 1:size(obj.M0,3)
            
            rot(x,y,z,:,:) =  [obj.M0(x,y,z,1)/obj.Ms(x,y,z), obj.M0(x,y,z,2)/obj.Ms(x,y,z), obj.M0(x,y,z,3)/obj.Ms(x,y,z);...
                              -obj.M0(x,y,z,2)/sq(x,y,z), obj.M0(x,y,z,1)/sq(x,y,z), 0;...
                              -obj.M0(x,y,z,1)*obj.M0(x,y,z,3)/(obj.Ms(x,y,z)*sq(x,y,z)), -obj.M0(x,y,z,2)*obj.M0(x,y,z,3)/(obj.Ms(x,y,z)*sq(x,y,z)), sq(x,y,z)/obj.Ms(x,y,z)];
        end
    end
end



for x = 1:size(obj.M0,1)
    for y = 1:size(obj.M0,2)
        for z = 1:size(obj.M0,3)
            for i = 1:3
            
                M0Loc(x,y,z,i) = rot(x,y,z,i,1)*obj.M0(x,y,z,1) + rot(x,y,z,i,2)*obj.M0(x,y,z,2) + rot(x,y,z,i,3)*obj.M0(x,y,z,3);        
                
            end 
        end   
    end
end



%%%%%%%%    Spatial FFT Calculation     %%%%%%%%
tic

nx	 = size(obj.M0,1);  
ny 	 = size(obj.M0,2);
dt   = 10e-12; 
comp = 2; % Out-of-plane component = 3
t_steps = 2000;

files = dir('*.omf');

data_box_all = ones(nx,ny,t_steps); 


for loop = 1:t_steps
	loop;
	
    DynamicFile = files(loop).name;
    M = obj.convM*obj.readFile(DynamicFile);
    
    for x = 1:size(obj.M0,1)
        for y = 1:size(obj.M0,2)
            for z = 1:size(obj.M0,3)
                for i = 1:3
            
                    MLoc(x,y,z,i) = rot(x,y,z,i,1)*M(x,y,z,1) + rot(x,y,z,i,2)*M(x,y,z,2) + rot(x,y,z,i,3)*M(x,y,z,3);        
                
                    end 
                end   
            end
    end 
    
   
    
    MLocDyn=MLoc - M0Loc;
    MDyn = M - obj.M0;
    data_box_all(:,:,loop) = MLocDyn(:,:,1,comp);
    data_box_phase(:,:,loop) = MDyn(:,:,1,comp);
    
end



data_box_phase = angle(fft(data_box_phase,[],3));
data_box_all = abs(fft(data_box_all,[],3));

save('FFT_data.mat','data_box_all','-v7.3');
save('FFT_phase_data.mat','data_box_phase','-v7.3');

%%%%%%%%      'Linear' FFT      %%%%%%%%

slice_sum=zeros(1,t_steps);
tf = isfinite(data_box_all);
 for k = 1:t_steps
    for i = 1:size(obj.M0,1)
        for j = 1:size(obj.M0,2)
            if  tf(i,j,k) == 1
            slice_sum(1,k) = slice_sum(1,k) + data_box_all(i,j,k); %seems to be returning zeros, not adding?
            else
            end
        end
     end
 end   
F = slice_sum; %%%% F is a matrix of Nan, data box all was working fine.
f = (1/(t_steps*dt))*[0:t_steps-1];
f_steps = f(2)-f(1);
f_low = (1e9)*input('What lower frequency band do you want the images for (GHz)? ');
f_high = (1e9)*input('What higher frequency band do you want the images for (GHz)? ');
f_low_bin = round(f_low/f_steps);
f_high_bin = round(f_high/f_steps);
r = ones(f_high_bin-f_low_bin,1);
for n = 0:1:t_steps-1
    if n<f_low_bin 
		continue
    elseif n>f_high_bin
		continue
    else
       r((n+1)-f_low_bin) = F(n);
    end
end
f1=figure;
clf(); set(gcf,'name','Memes','numbertitle','off') 
plot((f_low_bin:1:f_high_bin)*f_steps/1e9,r,'k');
title('Fourier Intensity');
xlabel('Frequency (GHz)');
ylabel('Power (arb. units)');
saveas(f1,strcat('Local_fft.fig'));
saveas(f1,strcat('Local_fft.png'));

%%%%%%%%    End of linear    %%%%%%%%

peaks = input('How many Frequency Bands do you want? ');
for i = 1:peaks

    time_intp=[dt:dt:t_steps*dt];
    frequencies = (1/(t_steps*dt))*[0:size(time_intp,2)-1].';

    frequency_step = frequencies(2)-frequencies(1); %1/(t_steps*dt)
    freq_low = 1e9*input('What lower frequency band do you want the images for (GHz)? ');
    freq_high = 1e9*input('What higher frequency band do you want the images for (GHz)? ');

    freq_low_bin = round(freq_low/frequency_step); 
    freq_high_bin = round(freq_high/frequency_step);

    G = fspecial('gaussian',[9 9], 3);
    
        

    for n = 0:1:t_steps-1 %n starts at zero as freq_low_bin was rounded down
	if n<freq_low_bin 
		continue
	elseif n>freq_high_bin
		continue %these two if's detail whether or not the specific slice you are looking at
        %is within the frequency range that you want.
    else
        
    x=[1,size(obj.M0,1)]*obj.xstepsize*1e9;
    y=[1,size(obj.M0,2)]*obj.ystepsize*1e9;
    f2 = figure;
    clf(); set(gcf,'name','Fourier Amplitude','numbertitle','off') 
    imagesc(x,y,transpose(data_box_all(:,:,(n+1))));
    axis xy;
	axis equal;
	axis tight;
    title(strcat('Fourier Amplitude (',num2str(frequencies(n+1)/1e9),'GHz)'));
    xlabel('x (nm)')
    ylabel('y (nm)')
    % c1 = colorbar('TickLabels',{'Min','','','','','','Max'}); 
    c1 = colorbar('Ticks',[min(min(data_box_all(:,:,(n+1)))),max(max(data_box_all(:,:,(n+1))))],'TickLabels',{'Min','Max'});
    c1.Label.String = 'Rel. Fourier Amplitude';
    colormap(jet);
	saveas(f2,strcat('Local_Spatial_fft_',num2str(frequencies(n+1)/1e9),'GHz.png'));
    matdata = transpose(data_box_all(:,:,(n+1)));
    save(strcat('FFT_intensity_',num2str(frequencies(n+1)/1e9),'GHz.mat'),'matdata','-v7.3');
    
    
    fftslice = data_box_all(:,round(size(data_box_all,2)/2),(n+1));         
    save('FFTslice.mat','fftslice');
    
    %Ip = imfilter(transpose(data_box_phase(:,:,(n+1))),G,'same','conv'); %have to add 1 due to down rounding
	%Ig = log10(Ig);
	f3=figure;
    clf(); set(gcf,'name','Phase','numbertitle','off') 
	imagesc(x,y,transpose(data_box_phase(:,:,(n+1))), [-pi pi]);
	colormap(hsv);
    c2 = colorbar('EastOutside','Ticks',[-pi,0,pi],'TickLabelInterpreter','latex', 'TickLabels', {'$-\pi$','0','$\pi$'}); 
    %set(gcf, 'Fontname','symbol');
    set(c2, 'Fontsize', 15);
    %c2.Label.String = 'Phase';
    %caxis([-1.5 4.5])
    axis xy;
	axis equal;
	axis tight;
    title(strcat('Fourier Phase (',num2str(frequencies(n+1)/1e9),'GHz)'));
    xlabel('x (nm)')
    ylabel('y (nm)')
	%axis off;
	saveas(f3,strcat('FFT_Phase_',num2str(frequencies(n+1)/1e9),'GHz.png'));
    
    %Ip = imfilter(transpose(data_box_phase(:,:,(n+1))),G,'same','conv'); %have to add 1 due to down rounding
	%Ig = log10(Ig);
	f4=figure;
    clf(); set(gcf,'name','Superpositionionio','numbertitle','off')
    conv = imagesc(x,y,transpose(data_box_phase(:,:,(n+1))), [-pi pi]);
    colormap(hsv);
    %set(gcf, 'Fontname','symbol');
    c3 = colorbar('EastOutside','Ticks',[-pi,0,pi],'TickLabelInterpreter','latex', 'TickLabels', {'$-\pi$','0','$\pi$'}); 
    %c3.Label.String = 'Phase';
    set(c3, 'Fontsize', 15);
    set(conv,'AlphaData',transpose(data_box_all(:,:,(n+1))),'AlphaDataMapping','scaled');
    %,'AlphaData',data_box_all(:,:,(n+1)));
    %set(data_box_all(:,:,(n+1)) , 'AlphaData', Ip )
    %caxis([-1.5 4.5])
    axis xy;
	axis equal;
	axis tight;
    title(strcat('Fourier Phase/Amplitude Superposition (',num2str(frequencies(n+1)/1e9),'GHz)'));
    xlabel('x (nm)')
    ylabel('y (nm)')
	%axis off;
	saveas(f4,strcat('FFT_Phase_conv_',num2str(frequencies(n+1)/1e9),'GHz.png'));
    end
    end

    end
    cd (currentpath);
end
%%%%%%%%     End of FFT     %%%%%%%%   
   
%%%%%%%%     Function to read .omf files     %%%%%%%
function M = readFile(obj,fName)     

        fid = fopen(fName);
        [IOmess, errnum] = ferror(fid);
        if (errnum ~= 0)
            disp(IOmess);
            return;
        end
        
        % read parameters file
        expr = '^#\s([\w\s:]+):\s([-.0-9e]+)';
        propertiesList = fieldnames(obj);
        obj.header = '';
        line = fgetl(fid);
        obj.header = strcat(obj.header,line,'\n');
        while (isempty(strfind(line,'Begin: Data Binary')))
            line = fgetl(fid);
            obj.header = strcat(obj.header,line,'\n');
            [~, ~, ~, ~, tokenStr, ~, splitStr] = regexp(line,expr);
            % read parameters
            if (size(tokenStr,1)>0)
                if (size(tokenStr{1,1},2)>1)
                    % seek properties
                    toks = tokenStr{1,1};
                    
                    if (strcmp(toks{1,1},'Desc:  Iteration'))
                        %obj.iteration = str2num(toks{1,2});
                    elseif (strcmp(toks{1,1},'Desc:  Total simulation time'))
                        %obj.totalSimTime = str2num(toks{1,2});
                    else
                        for i=1:size(propertiesList,1)
                            if(strcmp(propertiesList{i,1},toks{1,1}))
                                prop = toks{1,1};
                                val = toks{1,2};
                                
                                %  Is it numerical value?
                                [num,status] = str2num(val);
                                if (status) % yes, it's numerical
                                    set(obj,prop,num)
                                else % no, it's string
                                    set(obj,prop,val)
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % determine file format

        if (~isempty(strfind(line,'8')))
            obj.format = 'double';
            testVal = obj.doubleTestVal;
        elseif (~isempty(strfind(line,'4')))
            obj.format = 'single';
            testVal = obj.singleTestVal;
        else
            disp('Unknown format');
            return
        end
        
        % read first test value
        fTestVal = fread(fid, 1, obj.format, 0, 'ieee-le');
        if ((strcmp(obj.format,'single') && (fTestVal == obj.singleTestVal))...
                || (strcmp(obj.format,'double') && (fTestVal == obj.doubleTestVal)))
            %disp('Correct format')
        else
            disp('Wrong format');
            return;
        end
        
        
        data = fread(fid, obj.xnodes*obj.ynodes*obj.znodes*obj.dim,...
            obj.format, 0, 'ieee-le');
        
        line = fgetl(fid);
        if (isempty(strfind(line ,'# End: Data')) && isempty(strfind(line,'# End: Segment')))
            disp('End of file is incorrect. Something wrong');
            fclose(fid);
            % return;
        else
            fclose(fid);
        end

        Mx = data(1:3:size(data,1));
        My = data(2:3:size(data,1));
        Mz = data(3:3:size(data,1));
        
        M(:,:,:,1) = reshape(Mx, [obj.xnodes obj.ynodes obj.znodes]);
        M(:,:,:,2) = reshape(My, [obj.xnodes obj.ynodes obj.znodes]);
        M(:,:,:,3) = reshape(Mz, [obj.xnodes obj.ynodes obj.znodes]);
end

end
end