clear all
clc
t = linspace(0,10,1001);
for k=1:length(t)
    tt=t(k);

    %Motion
    %case1
    syms r1 r2 z3 th1 th2 ;
    r1=200/sqrt(3);
    r2=100;
    th1=pi;
    th2=pi/6-1*tt;
    syms rd3 thd2 thd3 ;
    thd2=-1;
    syms rdd3 thdd2 thdd3 ;
    thdd2=0;
    eqx1 = r1*(cos(th1)+1i*sin(th1))+r2*(cos(th2)+1i*sin(th2))+abs(z3)*(cos(angle(z3))+1i*sin(angle(z3)))==0 ;
    sol=solve(eqx1,z3);
    r3=double(abs(sol));
    th3=double(angle(sol));
    eqv1 = 1i*thd2*r2*(cos(th2)+1i*sin(th2))+rd3*(cos(th3)+1i*sin(th3))+1i*thd3*r3*(cos(th3)+1i*sin(th3))==0;
    [rd3, thd3]=solve(eqv1,[rd3 thd3],"Real",true);
    rd3=double(rd3);
    thd3=double(thd3);
    eqa1 = (1i*r2*thdd2*(cos(th2)+1i*sin(th2))-r2*thd2^2*(cos(th2)+1i*sin(th2)))+(1i*r3*thdd3*(cos(th3)+1i*sin(th3))-r3*thd3^2*(cos(th3)+1i*sin(th3)))+(rdd3*(cos(th3)+1i*sin(th3))+2i*thd3*rd3*(cos(th3)+1i*sin(th3)))==0;
    [rdd3, thdd3]=solve(eqa1,[rdd3 thdd3],"Real",true);
    rdd3=double(rdd3);
    thdd3=double(thdd3);
    %case2
    syms r0 r4 r5 th0 th4 th5;
    r0=300/sqrt(3);
    th0=2*pi;
    th4=pi+th3;
    th5=-pi/2;
    syms rd4 rd5 thd4;
    thd4=thd3;
    syms rdd4 rdd5 thdd4;
    thdd4=thdd3;
    eqx2 = r0*(cos(th0)+1i*sin(th0))+r2*(cos(th2)+1i*sin(th2))+r4*(cos(th4)+1i*sin(th4))+r5*(cos(th5)+1i*sin(th5))==0;
    [r4, r5]=solve(eqx2,[r4 r5],"Real",true);
    r4=double(r4);
    r5=double(r5);
    eqv2 = 1i*thd2*r2*(cos(th2)+1i*sin(th2))+rd4*(cos(th4)+1i*sin(th4))+1i*thd4*r4*(cos(th4)+1i*sin(th4))+rd5*(cos(th5)+1i*sin(th5))==0;
    [rd4, rd5]=solve(eqv2,[rd4 rd5],"Real",true);
    rd4=double(rd4);
    rd5=double(rd5);
    eqa2 = (1i*r2*thdd2*(cos(th2)+1i*sin(th2))-r2*thd2^2*(cos(th2)+1i*sin(th2)))+(1i*r4*thdd4*(cos(th4)+1i*sin(th4))-r4*thd4^2*(cos(th4)+1i*sin(th4)))+(rdd4*(cos(th4)+1i*sin(th4))+2i*thd4*rd4*(cos(th4)+1i*sin(th4)))+(rdd5*(cos(th5)+1i*sin(th5)))==0;
    [rdd4, rdd5]=solve(eqa2,[rdd4 rdd5],"Real",true);
    rdd4=double(rdd4);
    rdd5=double(rdd5);

    %Forces
    %case 2
    syms m2 a2x a2y r2x r2y I2 alpha2 ;
    m2=1*(r2/1000);
    a2x=(1/2)*((-(r2/1000)*(thd2)^2)*cos(th2)+((r2/1000)*(thdd2))*sin(th2));
    a2y=(1/2)*((-(r2/1000)*(thd2)^2)*sin(th2)+(-(r2/1000)*(thdd2))*cos(th2));
    r2x=(1/2)*(1/1000)*(r2*cos(th2));
    r2y=(1/2)*(1/1000)*(r2*sin(th2));
    I2=(1/12)*m2*((r2/1000)^2);
    alpha2=thdd2;
    syms fb2x fb2y fo2x fo2y T2 ;
    eq1= fb2x+fo2x==m2*a2x;
    eq2= fb2y+fo2y==m2*a2y;
    eq3= T2+(r2x*fb2y-r2y*fb2x)-(r2x*fo2y-r2y*fo2x)==I2*alpha2;
    %case 3
    syms m3 a3x a3y rcx rcy rb3x rb3y rdx rdy I3 alpha3 ;
    m3=1*((1100/sqrt(3))/1000);
    r3cm=(1/1000)*((550/sqrt(3))-(1100/sqrt(3)-abs(r3)-abs(r4)));
    R=sqrt((r5^2)+(500/sqrt(3))^2);
    Rd=(r5*rd5)/(R);
    Rdd=(rd5^2+r5*rdd5-Rd^2)/(R);
    rd3cm=(1/1000)*Rd;
    rdd3cm=(1/1000)*Rdd;
    a3x0=(-0.550/sqrt(3))*((thd3.^2)*cos(-th3)+(-thdd3).*sin(-th3));
    %Ma3x=Mrdd3cm.*cos(pi-Mth4)-2.*Mrd3cm.*Mthd4.*sin(pi-Mth4)+Mr3cm.*Mthdd4.*sin(pi-Mth4)-Mr3cm.*((-Mthd4).^2).*cos(pi-Mth4);
    a3y0=rdd3cm*sin(pi-th4)-2*rd3cm*thd4*cos(pi-th4)-r3cm*thdd4*cos(pi-th4)-r3cm*((-thd4)^2)*sin(pi-th4);
    a3x=a3x0;
    a3y=a3y0;
    rcx=(1)*(r3cm*cos(-th3));
    rcy=(-1)*(r3cm*sin(-th3));
    rb3x=(1)*((r3cm-r3/1000)*cos(-th3));
    rb3y=(-1)*((r3cm-r3/1000)*sin(-th3));
    rdx=(-1/1000)*(550/sqrt(3))*cos(-th3);
    rdy=(1/1000)*(550/sqrt(3))*sin(-th3);
    I3=(1/12)*m3*((1.1/sqrt(3))^2);
    alpha3=thdd4;
    syms fcx fcy fb3x fb3y fdx fdy ;
    eq4= fb3x+fcx+fdx==m3*a3x;
    eq5= fb3y+fcy+fdy==m3*a3y;
    eq6= (rcx*fcy-rcy*fcx)+(rb3x*fb3y-rb3y*fb3x)+(rdx*fdy-rdy*fdx)==I3*alpha3;
    %case C
    syms mc aocx aocy theta3 ;
    mc=0.3;
    aocx=0;
    aocy=0;
    theta3=pi-th4;
    syms focx focy ;
    eq7= -fcx+focx==mc*aocx;
    eq8= -fcy+focy==mc*aocy;
    eq9= fcx*cos(theta3)-fcy*sin(theta3)==0;
    %case B
    syms mb abx aby ;
    mb=0.3;
    abx=2*a2x;
    aby=2*a2y;
    eq10= -fb3x-fb2x==mb*abx;
    eq11= -fb3y-fb2y==mb*aby;
    eq12= fb3x*cos(theta3)-fb3y*sin(theta3)==0;
    %case D
    syms md adx ady mu ;
    md=0.3;
    adx=0;
    ady=(1/1000)*rdd5;
    if (rd5 >=0)
        mu=+0.25;
    elseif (rd5 < 0)
        mu=-0.25;
    end
    syms fsx fsy ;
    eq13= -fdx+fsx==md*adx;
    eq14= -fdy+fsy==md*ady;
    eq15= fsy==mu*fsx;

    eqs =[eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10 eq11 eq12 eq13 eq14 eq15];
    vars = [fb2x fb2y fo2x fo2y T2 fcx fcy fb3x fb3y fdx fdy focx focy fsx fsy];
    sol = solve(eqs,vars);

    Mr2(k)=r2;
    Mr3(k)=r3;
    Mr4(k)=r4;
    Mr5(k)=r5;
    Mrd3(k)=rd3;
    Mrd4(k)=rd4;
    Mrd5(k)=rd5;
    Mrdd3(k)=rdd3;
    Mrdd4(k)=rdd4;
    Mrdd5(k)=rdd5;
    Mth2(k)=th2;
    Mth3(k)=th3;
    Mth4(k)=th4;
    Mthd2(k)=thd2;
    Mthd3(k)=thd3;
    Mthd4(k)=thd4;
    Mthdd2(k)=thdd2;
    Mthdd3(k)=thdd3;
    Mthdd4(k)=thdd4;

    Mr2x(k)=r2x;
    Mr2y(k)=r2y;
    Mrb3x(k)=rb3x;
    Mrb3y(k)=rb3y;
    Mrcx(k)=rcx;
    Mrcy(k)=rcy;
    Mrdx(k)=rdx;
    Mrdy(k)=rdy;
    Ma2x(k)=a2x;
    Ma2y(k)=a2y;
    Ma3x(k)=a3x;
    Ma3y(k)=a3y;
    Mabx(k)=abx;
    Maby(k)=aby;
    Mady(k)=ady;
    MRd(k)=Rd;
    MRdd(k)=Rdd;
    Mr3cm(k)=r3cm;
    Mrd3cm(k)=rd3cm;
    Mrdd3cm(k)=rdd3cm;

    MT2(k)=double(sol.T2);
    Mfb2x(k)=double(sol.fb2x);
    Mfb2y(k)=double(sol.fb2y);
    Mfo2x(k)=double(sol.fo2x);
    Mfo2y(k)=double(sol.fo2y);
    Mfcx(k)=double(sol.fcx);
    Mfcy(k)=double(sol.fcy);
    Mfb3x(k)=double(sol.fb3x);
    Mfb3y(k)=double(sol.fb3y);
    Mfdx(k)=double(sol.fdx);
    Mfdy(k)=double(sol.fdy);
    Mfocx(k)=double(sol.focx);
    Mfocy(k)=double(sol.focy);
    Mfsx(k)=double(sol.fsx);
    Mfsy(k)=double(sol.fsy);
    clc
    disp([num2str((k/length(t))*100),'%'])
end
clc

S=input('Do you want the animation? Y/N : ','s');
if(S=='y' || S=='Y')
    figure(1)
    pause(3);
    for k=1:length(t)
        clf
        tt=t(k);
        r2=Mr2(k);
        r4=Mr4(k);
        r5=Mr5(k);
        th2=Mth2(k);
        th4=Mth4(k);
        x0=[-510,510,510,-510,-510];
        y0=[-510,-510,510,510,-510];
        x01=[0,(500/sqrt(3)),(500/sqrt(3)),0,0];
        y01=[-500,-500,500,500,-500];
        x02=[0, 500/sqrt(3)];
        y02=[0 , 0];
        x2=[(300/sqrt(3)),((300/sqrt(3))+ (r2*cos(th2)))];
        y2=[ 0, (r2*sin(th2))];
        x3=[((300/sqrt(3))+(r2*cos(th2))),(500/sqrt(3))];
        y3=[(r2*sin(th2)),0];
        x4=[((300/sqrt(3))+(r2*cos(th2))),0];
        y4=[(r2*sin(th2)),((r2*sin(th2))+r4*sin(th4))];
        x5=((300/sqrt(3))+ (r2*cos(th2)));
        y5=(r2*sin(th2));
        x6=0;
        y6=r5;
        plot(x0,y0,'k',x01,y01,'k',x02,y02,'k')
        hold on
        axis equal
        title(['at t= ',num2str(tt),' seconds'])
        plot(300/sqrt(3),0,'Marker','o','MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319')
        plot(500/sqrt(3),0,'Marker','o','MarkerEdgeColor','#D95319','MarkerFaceColor','#D95319')
        plot(x2,y2,'color','#0072BD')
        plot(x3,y3,'color','#A2142F')
        plot(x4,y4,'color','#A2142F')
        plot(x5,y5,'Marker','o','MarkerEdgeColor','#77AC30','MarkerFaceColor','#77AC30')
        plot(x6,y6,'Marker','s','MarkerEdgeColor','#77AC30','MarkerFaceColor','#77AC30')

        drawnow
        ax = gca;
        ax.Units = 'pixels';
        pos = ax.Position;
        ti = ax.TightInset;
        rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
        F(k)=getframe(ax,rect);
    end
    mywriter = VideoWriter('SP1_810698349');
    mywriter.FrameRate = 100;
    open(mywriter);
    writeVideo(mywriter,F);
    close(mywriter)
elseif (S=='n' || S=='N' || S=='y' || S=='Y')
    figure(2)
    subplot(3,2,1)
    plot(t,MRd,'color','#0072BD')
    grid on
    title('C4 Velocity')
    ylabel('mm/s')
    subplot(3,2,2)
    plot(t,MRdd,'color','#A2142F')
    grid on
    title('C4 Acceleration')
    ylabel('mm/s^2')
    subplot(3,2,3)
    plot(t,Mrd5,'color','#0072BD')
    grid on
    title('D Velocity')
    ylabel('mm/s')
    subplot(3,2,4)
    plot(t,Mrdd5,'color','#A2142F')
    grid on
    title('D Acceleration')
    ylabel('mm/s^2')
    subplot(3,2,5)
    plot(t,Mthd3,'color','#0072BD')
    grid on
    title('C5 Angular Velocity')
    ylabel('rad/s')
    subplot(3,2,6)
    plot(t,Mthdd3,'color','#A2142F')
    grid on
    title('C5 Angular Acceleration')
    ylabel('rad/s^2')

    figure(3)
    subplot(4,2,1)
    plot(t,Mfcx,'color','#0072BD')
    grid on
    title('C5-X')
    ylabel('Newton')
    subplot(4,2,3)
    plot(t,Mfocy,'color','#A2142F')
    grid on
    title('C5-Y')
    ylabel('Newton')
    subplot(4,2,2)
    plot(t,Mfdx,'color','#0072BD')
    grid on
    title('D-X')
    ylabel('Newton')
    subplot(4,2,4)
    plot(t,Mfdy,'color','#A2142F')
    grid on
    title('D-Y')
    ylabel('Newton')
    subplot(4,2,5)
    plot(t,(sqrt(Mfcx.^2+Mfcy.^2)),'color','#0072BD')
    grid on
    title('C5-Mag')
    ylabel('Newton')
    subplot(4,2,6)
    plot(t,(sqrt(Mfdx.^2+Mfdy.^2)),'color','#A2142F')
    grid on
    title('D-Mag')
    ylabel('Newton')
    subplot(4,2,[7,8])
    plot(t,(sqrt(Mfdx.^2+Mfdy.^2)),'color','#77AC30')
    grid on
    title('T2')
    ylabel('Newton-meter')
end
