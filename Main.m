
clear,clc
close all


L = 4.0;
C = 1.0;
Alpha = 5.0;
BetaI = pi/25.0;
STime = 2.0;
Epsilon_e = 0.01;
Epsilon_i = 0.01;
NPoint = [128 ];
NP = 320;
BVal = [0.25*pi ];
CFL = [0.25 ];
INIT_Name = {'SIN','STEP','EXP'};
METHOD_Name = {'UPWIND','LEAPFROG','LAXWENDROFF','LAX','BEAMWARMING','DAMPED-BEAMWARMING'};


for IWave = 1:3
    
    Beta = BetaI;
    N = NP;
    Time = 0.0;
    dx = L/N;
    X = [0.0:dx:L];
    UI = WAVE( X,Time,IWave,Alpha,Beta,dx,N,C );
    
    fig = figure(1);
    fig.Units = 'normalized';
    fig.OuterPosition = [0 0 1 1];
    plot(X,UI,'k','Linewidth',2);grid on;
    xlabel('\fontsize{12}X','FontWeight','bold');
    ylabel('\fontsize{12}U','FontWeight','bold');
    title(['Initial Wave : ',INIT_Name{IWave}]);
    saveas(fig,['Initial_',INIT_Name{IWave},'.jpg'],'jpg');
    close all;
    
end

fprintf('================================================\n');
clearvars -except ICFL IWave IMethod INum L C Alpha Beta BetaI STime Epsilon_e Epsilon_i NPoint CFL INIT_Name METHOD_Name Xp Yp COLOR Leg IBeta BVal NP
COLOR = hsv(length(CFL));
for IWave = 1:3
    for IMethod = 1:5
        fig = figure(2);
        fig.Units = 'normalized';
        fig.OuterPosition = [0 0 1 1];
        for ICFL = 1:length(CFL)
            for INum = 1:length(NPoint)
                
                Beta = BetaI;
                Sigma = CFL(ICFL);
                N = NPoint(INum);
                n = 1;
                Time = 0.0;
                dx = L/N;
                dt = (Sigma*dx)/C;
                NTime = fix(STime/dt);
                X = [0.0:dx:L];
                UI = WAVE( X,Time,IWave,Alpha,Beta,dx,N,C );
                
                u(n,:) = UI;
                if (IMethod == 1)
                    [ u,G ] = UPWIND( u,C,n,N,dt,dx,NTime,Sigma,Beta );
                elseif (IMethod == 2)
                    [ u,G ] = LEAPFROG( u,C,n,N,dt,dx,NTime,Sigma,Beta );
                elseif (IMethod == 3)
                    [ u,G ] = LAXWENDROFF( u,C,n,N,dt,dx,NTime,Sigma,Beta );
                elseif (IMethod == 4)
                    [ u,G ] = LAX( u,C,n,N,dt,dx,NTime,Sigma,Beta );
                elseif (IMethod == 5)
                    [ u,G ] = BEAMWARMING( u,C,n,N,dt,dx,NTime,Sigma,Beta );
                end
                
                UE = WAVE( X,NTime*dt,IWave,Alpha,Beta,dx,N,C );
                
                Norm = sqrt(sum((u(NTime+1,:)-UE).^2))/(N+1);
                
                Xp(INum) = dx;
                Yp(INum) = Norm;
                
                fprintf('\t\t%d\t\t%d\t\t\t%d\t\t\t%d\n',IWave,IMethod,ICFL,INum);
                
                clearvars -except ICFL IWave IMethod INum L C Alpha Beta BetaI STime Epsilon_e Epsilon_i NPoint CFL INIT_Name METHOD_Name Xp Yp COLOR Leg IBeta BVal NP fig
            end
            [xData, yData] = prepareCurveData( log10(Xp), log10(Yp) );
            [fitresult, gof] = fit( xData, yData, 'poly1' );
            a = fitresult.p1;
            b = fitresult.p2;
            if(abs(a)<0.000001)
                a = 0.0;
            end
            if(abs(b)<0.000001)
                b = 0.0;
            end
            xt = linspace(min(Xp),max(Xp),200);
            yt = (xt.^a)*(10.0^b);
            
            loglog(Xp,Yp,'-s','linewidth',2.0,'color',COLOR(ICFL,:));hold on;grid on;
            loglog(xt,yt,'r--','linewidth',2.0,'color',COLOR(ICFL,:));
            Leg{2*ICFL-1} = ['CFL = ',num2str(CFL(ICFL))];
            Leg{2*ICFL} = ['SLOPE = ',num2str(a)];
            
            clearvars Xp Yp xData yData xt yt ft fitresult gof
        end
        title([METHOD_Name{IMethod},'-',INIT_Name{IWave}]);
        xlabel('\fontsize{12}Log(DX)','FontWeight','bold');
        ylabel('\fontsize{12}Log(Norm)','FontWeight','bold');
        legend(Leg);
        saveas(fig,['ORDER_',METHOD_Name{IMethod},'_',INIT_Name{IWave},'.jpg'],'jpg');
        close all;
    end
end

fprintf('================================================\n');
clearvars Leg
COLOR = hsv(length(CFL));
for IWave = 1:3
    for IMethod = 1:5
        fig = figure(3);
        fig.Units = 'normalized';
        fig.OuterPosition = [0 0 1 1];
        for ICFL = 1:length(CFL)
            for IBeta = 1:length(BVal)
                
                Sigma = CFL(ICFL);
                N = NP;
                Beta = BVal(IBeta);
                n = 1;
                Time = 0.0;
                dx = L/N;
                dt = (Sigma*dx)/C;
                NTime = fix(STime/dt);
                X = [0.0:dx:L];
                UI = WAVE( X,Time,IWave,Alpha,Beta,dx,N,C );
                
                u(n,:) = UI;
                if (IMethod == 1)
                    [ u,G ] = UPWIND( u,C,n,N,dt,dx,NTime,Sigma,Beta );
                elseif (IMethod == 2)
                    [ u,G ] = LEAPFROG( u,C,n,N,dt,dx,NTime,Sigma,Beta );
                elseif (IMethod == 3)
                    [ u,G ] = LAXWENDROFF( u,C,n,N,dt,dx,NTime,Sigma,Beta );
                elseif (IMethod == 4)
                    [ u,G ] = LAX( u,C,n,N,dt,dx,NTime,Sigma,Beta );
                elseif (IMethod == 5)
                    [ u,G ] = BEAMWARMING( u,C,n,N,dt,dx,NTime,Sigma,Beta );
                end
                
                Xp(IBeta) = Beta;
                Yp(IBeta) = abs(G);
                Zp(IBeta) = (-Sigma*Beta-angle(G));
                
                fprintf('\t\t%d\t\t%d\t\t\t%d\t\t\t%d\n',IWave,IMethod,ICFL,IBeta);
                
                clearvars -except ICFL IWave IMethod INum L C Alpha Beta BetaI STime Epsilon_e Epsilon_i NPoint CFL INIT_Name METHOD_Name Xp Yp Zp COLOR Leg IBeta BVal NP fig
            end
            Leg{ICFL} = ['CFL = ',num2str(CFL(ICFL))];
            
            subplot(2,1,1);plot(Xp,Yp,'-s','linewidth',2.0,'color',COLOR(ICFL,:));hold on;grid on;
            title([METHOD_Name{IMethod},'-',INIT_Name{IWave}]);
            xlabel('\fontsize{12}Beta','FontWeight','bold');
            ylabel('\fontsize{12}Altitude Error','FontWeight','bold');
            subplot(2,1,2);plot(Xp,Zp,'-s','linewidth',2.0,'color',COLOR(ICFL,:));hold on;grid on;
            xlabel('\fontsize{12}Beta','FontWeight','bold');
            ylabel('\fontsize{12}Phase Error','FontWeight','bold');
            
            clearvars Xp Yp Zp
        end
        legend(Leg);
        saveas(fig,['ALTITUDE_PHASE_',METHOD_Name{IMethod},'_',INIT_Name{IWave},'.jpg'],'jpg');
        close all;
    end
end

fprintf('================================================\n');
clearvars Leg
COLOR = hsv(length(BVal));
for IWave = 1:3
    for IMethod = 1:6
        for ICFL = 1:length(CFL)
            for IBeta = 1:length(BVal)
                
                Sigma = CFL(ICFL);
                N = NP;
                Beta = BVal(IBeta);
                n = 1;
                Time = 0.0;
                dx = L/N;
                dt = (Sigma*dx)/C;
                NTime = fix(STime/dt);
                X = [0.0:dx:L];
                UI = WAVE( X,Time,IWave,Alpha,Beta,dx,N,C );
                
                u(n,:) = UI;
                if (IMethod == 1)
                    [ u,G ] = UPWIND( u,C,n,N,dt,dx,NTime,Sigma,Beta );
                elseif (IMethod == 2)
                    [ u,G ] = LEAPFROG( u,C,n,N,dt,dx,NTime,Sigma,Beta );
                elseif (IMethod == 3)
                    [ u,G ] = LAXWENDROFF( u,C,n,N,dt,dx,NTime,Sigma,Beta );
                elseif (IMethod == 4)
                    [ u,G ] = LAX( u,C,n,N,dt,dx,NTime,Sigma,Beta );
                elseif (IMethod == 5)
                    [ u,G ] = BEAMWARMING( u,C,n,N,dt,dx,NTime,Sigma,Beta );
                elseif (IMethod == 6)
                    [ u ] = DAMPED_BEAMWARMING( u,C,n,N,dt,dx,NTime,Sigma,Beta,Epsilon_e,Epsilon_i );
                end
                fig = figure(4);
                fig.Units = 'normalized';
                fig.OuterPosition = [0 0 1 1];
                if (IBeta == 1)
                    UE = WAVE( X,NTime*dt,IWave,Alpha,Beta,dx,N,C );
                    plot(X,UE,'-k','linewidth',2.0);hold on;grid on;
                    Leg{IBeta} = 'EXACT';
                end
                plot(X,u(NTime+1,:),'-','linewidth',2.0,'color',COLOR(IBeta,:));
                Leg{IBeta+1} = ['Beta = ',num2str(CFL(IBeta))];
                
                fprintf('\t\t%d\t\t%d\t\t\t%d\t\t\t%d\n',IWave,IMethod,ICFL,IBeta);
                
                clearvars -except ICFL IWave IMethod INum L C Alpha Beta BetaI STime Epsilon_e Epsilon_i NPoint CFL INIT_Name METHOD_Name Xp Yp Zp COLOR Leg IBeta BVal NP fig
            end
            title([METHOD_Name{IMethod},'-',INIT_Name{IWave},'- CFL = ',num2str(CFL(ICFL))]);
            xlabel('\fontsize{12}X','FontWeight','bold');
            ylabel('\fontsize{12}U','FontWeight','bold');
            legend(Leg);
            saveas(fig,['COMPARE_',METHOD_Name{IMethod},'_',INIT_Name{IWave},'_',num2str(CFL(ICFL)),'.jpg'],'jpg');
            close all;
            
        end
    end
end





