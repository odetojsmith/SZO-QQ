x = 0:0.001:1;
x_cut = 0:0.0001:(sqrt(0.5)-0.5);


% 
x2 = [x, fliplr(x)];
inBetween = [func_1(x), fliplr(func_2(x))];

plot(x,func_1(x),'LineWidth',0.5);
hold on
plot(x,func_2(x),'LineWidth',0.5);
plot(0*x,0:0.001:1,'LineWidth',0.5)
set(gca,'Fontname','latex','linewidth',1.5,'FontSize',25,'TickLabelInterpreter','latex');
% 
fill(x2, inBetween, 1.3*[0.5660 0.6740 0.1880],'FaceAlpha',0.5);
plot(0,0,'-p','MarkerSize',40,'MarkerFaceColor','red')

hold on
xlabel('$x^{(1)}$','Interpreter','latex');
ylabel('$x^{(2)}$','Interpreter','latex');

xlim([0 1.1])
ylim([0 1.2])


function y = func_1(x)
y = x.^2;
end

function y = func_2(x)
y = 0*x+1;
end

function y = func_3(x)
y = 0.5-sqrt(0.5-(x+0.5).^2);
end

function y = func_4(x)
y = 0.5+sqrt(0.5-(x+0.5).^2);
end