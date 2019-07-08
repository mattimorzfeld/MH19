function MakePlot(y,z,t,tObs,xAll,rmse,spread,rho,alg)
figure
subplot(511)
plot(t,y(1,:),'LineWidth',2)
subplot(512)
plot(t,y(2,:),'LineWidth',2)
subplot(513)
plot(t,y(3,:),'LineWidth',2)

subplot(511)
hold on, plot(tObs,z(1,:),'.','MarkerSize',20)
subplot(512)
hold on, plot(tObs,z(2,:),'.','MarkerSize',20)
subplot(513)
hold on, plot(tObs,z(3,:),'.','MarkerSize',20)

subplot(511)
hold on,plot(t,xAll(1,:),'LineWidth',2)
box off
title(alg)

subplot(512)
hold on,plot(t,xAll(2,:),'LineWidth',2)
box off
subplot(513)
hold on,plot(t,xAll(3,:),'LineWidth',2)
box off


subplot(514)
hold on, plot(tObs,rmse,'LineWidth',2)
hold on, plot(tObs,spread,'LineWidth',2)
ylabel('RMSE and Spread')

subplot(515)
hold on, plot(tObs,rho,'LineWidth',2)
ylabel('\rho')

xlabel('Time')
set(gcf,'Color','w')
set(gca,'FontSize',12)

drawnow
