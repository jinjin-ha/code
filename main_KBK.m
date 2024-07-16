clear
%% INPUT
% ####################################################################
% S=[S11 S22 S33 S23 S31 S12];
% E=[E11 E22 E33 E23 E31 E12];
% YPSG=[PS11 PS22 0 0 0 0 PS12];    % for plane stress 
% Ang in [degree]; delta in [rad]

% TITLE
material='AA6022-T4(Same material with IJMS)';

NUT=7; Ang=[0:15:90];           % Input for UT # of Input & Angle
NBB=1;                          % Input for BB # of Input, [NBB=2] both, [NBB=1] RB, [NBB=0] SigB
NPS=3; PAng=[0 45 90];          % Input for PST # of Input & Angle

YF=3;   % Yield function choice: [YF=1] von-Mises, [YF=2] Hosford, [YF=3]=Yld2000-2d
MM=8;   % Yield function exponent
KK=2;   % Yield function constant
% a=[0.970233 1.054046 1.253166 1.12812 1.065 1.2534 0.9400 0.9089];   % Par-1
a=[0.95421 1.027886 1.053014 1.088182 1.008 0.9754 0.9068 1.0493];   % Par-2

% Input for uniaxial & biaxial tension
% Data [0 15 30 45 60 75 90 BB]; Put "1.0" if you do not have SigBB data
sr=[1.000 1.010 1.016 0.987 0.976 0.966 0.944 1.0]; % experimental stress ratio for uniaxial tension and biaxial
Er=[0.793 0.663 0.465 0.317 0.352 0.423 0.510 1.080]; % experimental r-value for uniaxial tension and biaxial

% Input for plane strain (By Yld2k-2d) in ##GLOBAL COORDINATE##
% EPSG Optimized by Par-2
EPSG=[1.081	0.597 0 0 0 0;        % PS00
      0.771	0.771 0 0 0 0.275;    % PS45
      0.647	1.004 0 0 0 0];       % PS90
% ####################################################################
%% 
% Isotropic
% sr=[1 1 1 1 1 1 1 1]; % experimental stress ratio for uniaxial tension and biaxial
% Er=[1 1 1 1 1 1 1 1]; % experimental r-value for uniaxial tension and biaxial

% Yr calculation
[~,Yr,YPSG,~]=func_yld(YF,MM,KK,a,Ang,PAng);

% YPSG for isotropic
% if YF==1 % if m=2
% YPSG1=[1.1547 0.5774 0 0 0 0]; % PS00: plane strain for 0
% YPSG2=[0.8661 0.8661 0 0 0 0.2886]; % PS45: plane strain for 45
% YPSG3=[0.5774 1.1547 0 0 0 0]; % PS90: plane strain for 90
% end
% if YF==2 && MM==6   % if m=6
% YPSG1=[1.1167 0.5584 0 0 0 0]; % PS00: plane strain for 0
% YPSG2=[0.8376 0.8376 0 0 0 0.2792]; % PS45: plane strain for 45
% YPSG3=[0.5584 1.1167 0 0 0 0]; % PS90: plane strain for 90
% end
% if YF==2 && MM==8   % if m=8
% YPSG1=[1.0894 0.5447 0 0 0 0]; % PS00: plane strain for 0
% YPSG2=[0.8171 0.8171 0 0 0 0.2724]; % PS45: plane strain for 45
% YPSG3=[0.5447 1.0894 0 0 0 0]; % PS90: plane strain for 90
% end

% if YF==3
% YPSG(1,:)=EPSG(1,:); % PS00: plane strain for 0
% YPSG(2,:)=EPSG(2,:); % PS45: plane strain for 45
% YPSG(3,:)=EPSG(3,:); % PS45: plane strain for 90
% end

%% KBK 
if YF==1;Ynm='vM';end
if YF==2;Ynm='HF';end
if YF==3;Ynm='Yld2k';end
for n=1:3   % [n=1]: UT, [n=2]: BT, [n=3]: PST
    
    if n==1 && NUT~=0 % Uniaxial tension
        for k=1:NUT
            % Stress and strain tensor rotation
            [S]=func_usg(Ang(k),sr(k));
            [EE]=func_ueg(Ang(k),Er(k));  % Experiment; Ang in [degree]
            [YE]=func_ueg(Ang(k),Yr(k));  % Yield function
            % cos, sigb, delta
            [cos]=func_cos(S);
            [SB]=func_phi(YF,MM,KK,a,S);
            [delt]=func_delt(EE,YE);

            cosU(k)=cos;
            deltU(k)=delt;
            phiU(k)=SB;
        end
    end

    if n==2 || n==1 || n==0 % Biaxial tension
            % Stress and strain tensor rotation
            S=[sr(NUT+1) sr(NUT+1) 0 0 0 0];   
            EE=[1 Er(NUT+1) -1-Er(NUT+1) 0 0 0];
            YE=[1 Yr(NUT+1) -1-Yr(NUT+1) 0 0 0];
            % cos, sigb, delta
            [cos]=func_cos(S);
            if n==2 || n==0; [SB]=func_phi(YF,MM,KK,a,S);end
            if n==2 || n==1; [delt]=func_delt(EE,YE);end
            cosB(1)=cos;
            deltB(1)=delt;
            phiB(1)=SB;
    end
    
    if n==3 && NPS~=0 % Plane strain tension
        for k=1:NPS
         	% cos, sigb, delta
            [cos]=func_cos(EPSG(k,:));
            [SB]=func_phi(YF,MM,KK,a,EPSG(k,:));
            cosP(k)=cos;
            deltP(k)=0;
            phiP(k)=SB;
        end
    end
    
end

%% DATA REPORT
% HEADER
fname0 = strcat('KBKrep-',Ynm,'.txt'); fid=fopen(fname0,'w');
fprintf(fid,'%s\n','## KBK representation ##');
fprintf(fid,'%s %s\n','Material:',material);
fprintf(fid,'%s\n','Experiment:');
% UT
if NUT==3; fprintf(fid,'%s\n','[UT] 3UT (0, 45, 90)');end
if NUT==7; fprintf(fid,'%s\n','[UT] 7UT (0, 15, 30, 45, 60, 75, 90)');end
% BB
if NBB==0; fprintf(fid,'%s\n','[BB] 1BB (SB)');end
if NBB==1; fprintf(fid,'%s\n','[BB] 1BB (RB)');end
if NBB==2; fprintf(fid,'%s\n','[BB] 1BB (SB, RB)');end
% PS
if NPS==1; fprintf(fid,'%s\n\n','[PS] 1PS (check angle)');end
if NPS==2; fprintf(fid,'%s\n\n','[PS] 2PS (check angle)');end
if NPS==3; fprintf(fid,'%s\n\n','[PS] 3PS (0, 45, 90)');end

fprintf(fid,'%s\n','Yield function:');
if YF==1; fprintf(fid,'%s\n\n','[von-Mises]');end
if YF==2; fprintf(fid,'%s\n','[Hosford]');fprintf(fid,'%2.1f  %s\n%2.1f  %s\n\n',MM,'M',KK,'K');end
if YF==3; fprintf(fid,'%s\n','[Yld2k-2d]');fprintf(fid,'%2.1f  %s\n%2.1f  %s\n',MM,'M',KK,'K');
    fprintf(fid,'\n%s\n','Parameter alpha');
    fprintf(fid,[repmat('%8.6s ',1,8),'\n'],'a1','a2','a3','a4','a5','a6','a7','a8');
    fprintf(fid,[repmat('%8.6f ',1,8),'\n\n'],a(1:8));
end

fprintf(fid,'%s\n','Input & Output:');
% UT
fprintf(fid,['%5s,',repmat('%10s,',1,6),'\n'],'[UT]','STRS-E','RVAL-E','RVAL-Y','COS','PHI','DELT');
form1='%5s,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,\n';
form2='%5s,%10.4f,%10.4f,%10.4f,%10.4f,%10s,%10.4f,\n';
form3='%5s,%10.4f,%10.4f,%10.4f,%10.4f,%10.4f,%10s,\n';
for i=1:NUT
    chrU = int2str(Ang(i));UTA= strcat('UT',chrU);
    fprintf(fid,form1,UTA,sr(i),Er(i),Yr(i),cosU(i),phiU(i),deltU(i));
end
fprintf(fid,'\n');
% BB
fprintf(fid,['%5s,',repmat('%10s,',1,6),'\n'],'[BB]','STRS-E','RVAL-E','RVAL-Y','COS','PHI','DEL');
if NBB==2; fprintf(fid,form1,'BB',sr(NUT+1),Er(NUT+1),Yr(NUT+1),cosB(1),phiB(1),deltB(1));end
if NBB==1; fprintf(fid,form2,'BB',sr(NUT+1),Er(NUT+1),Yr(NUT+1),cosB(1),'NONE',deltB(1));end
if NBB==0; fprintf(fid,form3,'BB',sr(NUT+1),Er(NUT+1),Yr(NUT+1),cosB(1),phiB(1),'NONE');end
fprintf(fid,'\n');
% PS
fprintf(fid,['%5s,',repmat('%10s,',1,6),'%15s\n'],'[PS]','SIGXX','SIGYY','SIGXY','COS','PHI','DEL','(Global Coord)');
for i=1:NPS
    chrP = int2str(PAng(i));PSA= strcat('PS',chrP);
    fprintf(fid,form1,PSA,EPSG(i,1),EPSG(i,2),EPSG(i,6),cosP(i),phiP(i),deltP(i));
end
fclose('all');
%% DATA PLOTTING
figure()
hold on
plot(cosU(:),phiU(:),'ko','displayName','UT');
if NBB==2 || NBB==0; plot(cosB(:),phiB(:),'bx','displayName','BB');end
plot(cosP(:),phiP(:),'rs','displayName','PS');

yline(1,'displayName','Yield condition');
xlim([-1 1])
% ylim([0.9 1.1])
xticks([-1 -0.5 0 0.5 1])
xlabel('cos(s_{R}:s)','fontsize',15)
ylabel('Normalized stress \sigma/\sigma_{0}','fontsize',15)
hold off
box on
legend show
legend('location','northwest')
fname1 = strcat('KBK-Phi-',Ynm,'.png');
print('-dpng', fname1);

figure()
hold on
plot(cosU(:),deltU(:),'ko','displayName','UT');
if NBB==2 || NBB==1; plot(cosB(:),deltB(:),'bx','displayName','BB');end
plot(cosP(:),deltP(:),'rs','displayName','PS');

yline(0,'displayName','Yield condition');
xlim([-1 1])
% ylim([-0.2 0.2])
xticks([-1 -0.5 0 0.5 1])
xlabel('cos(s:s)','fontsize',15)
ylabel('\delta (rad)','fontsize',15)
hold off
box on
legend show
legend('location','southwest')
fname1 = strcat('KBK-Delt-',Ynm,'.png');
print('-dpng', fname1);


fclose('all');