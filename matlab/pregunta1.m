% PREGUNTA1 Actuador: RENC, RESC y respuesta total.
% Requiere MATLAB y Control System Toolbox. Ejecutar este script completo.
clear; close all; clc;
raiz = fileparts(fileparts(mfilename('fullpath')));
salida = fullfile(raiz,'figuras_matlab');
if ~exist(salida,'dir'), mkdir(salida); end

% Estados fisicos: x1 = y, x2 = dy/dt. No usar ss(tf(...)) con este x0.
A = [0 1; -5 -2]; B = [0;10]; C = [1 0]; D = 0;
sys = ss(A,B,C,D);
x0 = [-2;4]; t = (0:0.005:10)';
[yrenc,~] = initial(sys,x0,t);
[yresc,~] = step(sys,t);
[ytotal,~] = lsim(sys,ones(size(t)),t,x0);
yrenc = yrenc(:); yresc = yresc(:); ytotal = ytotal(:);

% Soluciones analiticas obtenidas mediante Laplace unilateral.
renc_teo = exp(-t).*(-2*cos(2*t)+sin(2*t));
resc_teo = 2-exp(-t).*(2*cos(2*t)+sin(2*t));
total_teo = 2-4*exp(-t).*cos(2*t);
errores = [max(abs(yrenc-renc_teo)), max(abs(yresc-resc_teo)), ...
    max(abs(ytotal-total_teo)), max(abs(ytotal-yrenc-yresc))];
assert(all(errores < 1e-8),'Fallo del contraste analitico o de superposicion.');
assert(all(real(eig(A)) < 0),'El modelo debe ser estable.');
fprintf('Polos del actuador:\n'); disp(eig(A));
fprintf('Errores maximos [RENC RESC Total Superposicion]:\n'); disp(errores);
fprintf('Prueba exacta: y(0)=-2; dy(0)=4; d2y(0)=12.\n');

figure('Theme','light','Color','w','Position',[100 100 900 800]);
teoricas = [renc_teo resc_teo total_teo];
numericas = [yrenc yresc ytotal];
titulos = {'RENC: initial','RESC: step','Respuesta total: lsim'};
for k=1:3
    subplot(3,1,k);
    plot(t,teoricas(:,k),'b-','LineWidth',1.5); hold on;
    plot(t(1:40:end),numericas(1:40:end,k),'ro','MarkerSize',4);
    grid on; xlabel('Tiempo [s]'); ylabel('Desviacion angular y');
    title(titulos{k}); legend('Analitica','MATLAB','Location','best');
end
print(gcf,fullfile(salida,'p1_respuestas.png'),'-dpng','-r200');
exportgraphics(gcf,fullfile(salida,'p1_respuestas.pdf'),'ContentType','vector');
figure('Theme','light','Color','w');
plot(real(eig(A)),imag(eig(A)),'rx','MarkerSize',12,'LineWidth',2);
grid on; xlim([-3 1]); ylim([-3 3]);
xlabel('Parte real [s^{-1}]'); ylabel('Parte imaginaria [s^{-1}]');
title('Polos del actuador: -1 +/- 2j');
print(gcf,fullfile(salida,'p1_polos.png'),'-dpng','-r200');
exportgraphics(gcf,fullfile(salida,'p1_polos.pdf'),'ContentType','vector');
writetable(table(t,yrenc,yresc,ytotal,renc_teo,resc_teo,total_teo), ...
    fullfile(salida,'p1_datos.csv'));
save(fullfile(salida,'p1_resultados.mat'),'A','B','C','D','x0','t', ...
    'yrenc','yresc','ytotal','errores');
