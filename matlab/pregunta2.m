% PREGUNTA2 Dos estanques: step, impulse y lsim.
% Requiere MATLAB y Control System Toolbox. Condiciones iniciales nulas.
clear; close all; clc;
raiz = fileparts(fileparts(mfilename('fullpath')));
salida = fullfile(raiz,'figuras_matlab');
if ~exist(salida,'dir'), mkdir(salida); end
A1=2; A2=1; R12=2; Rs1=2; Rs2=2;
A = [-(1/Rs1+1/R12)/A1, 1/(A1*R12); ...
    1/(A2*R12), -(1/R12+1/Rs2)/A2];
B = [1/A1;0]; C=eye(2); D=zeros(2,1);
sys=ss(A,B,C,D); t=(0:0.01:25)'; x0=[0;0];
[hstep,~]=step(sys,t); hstep=squeeze(hstep);
[himp,~]=impulse(sys,t); himp=squeeze(himp);
u=2*sin(0.4*t)+2.5;
[hvar,~]=lsim(sys,u,t,x0);

% Referencias exactas; z es la respuesta fasorial a frecuencia 0.4 rad/s.
K=-A\B; w=0.4; z=(1i*w*eye(2)-A)\B;
hstep_teo=zeros(length(t),2); himp_teo=hstep_teo; hvar_teo=hstep_teo;
for k=1:length(t)
    E=expm(A*t(k));
    hstep_teo(k,:)=(K-E*K)';
    himp_teo(k,:)=(E*B)';
    hvar_teo(k,:)=(2.5*K+2*imag(z*exp(1i*w*t(k))) ...
        -E*(2.5*K+2*imag(z)))';
end
err_step=max(max(abs(hstep-hstep_teo)));
err_imp=max(max(abs(himp-himp_teo)));
err_var=max(max(abs(hvar-hvar_teo)));
% lsim interpola la entrada muestreada: admitir ese error, no redondearlo a 0.
assert(err_step<1e-8 && err_imp<1e-8,'Fallo en step o impulse.');
assert(err_var<1e-4,'Error excesivo de lsim respecto al seno continuo.');
assert(all(real(eig(A))<0),'El modelo debe ser estable.');
assert(max(abs(K-[4/3;2/3]))<1e-12,'Equilibrio incorrecto.');
fprintf('Matrices del modelo:\n'); disp(A); disp(B); disp(C); disp(D);
fprintf('Polos:\n'); disp(eig(A));
fprintf('Equilibrio para u=1 [m]:\n'); disp(K);
fprintf('Errores maximos [step impulse lsim]: %.3e %.3e %.3e\n', ...
    err_step,err_imp,err_var);
fprintf('Niveles a los 25 s [m]:\n'); disp([hstep(end,:);hvar(end,:)]);

figure('Theme','light','Color','w','Position',[100 100 900 750]);
subplot(2,1,1); plot(t,hstep,'LineWidth',1.5); grid on;
xlabel('Tiempo [s]'); ylabel('Nivel [m]'); title('Escalon: u=1 m^3/s');
legend('h_1','h_2','Location','best');
subplot(2,1,2); plot(t,himp,'LineWidth',1.5); grid on;
xlabel('Tiempo [s]'); ylabel('Nivel [m]');
title('Impulso positivo: volumen de 1 m^3 en el estanque 1');
legend('h_1','h_2','Location','best');
print(gcf,fullfile(salida,'p2_escalon_impulso.png'),'-dpng','-r200');
exportgraphics(gcf,fullfile(salida,'p2_escalon_impulso.pdf'),'ContentType','vector');
figure('Theme','light','Color','w','Position',[100 100 900 700]);
subplot(2,1,1); plot(t,u,'k','LineWidth',1.4); grid on;
xlabel('Tiempo [s]'); ylabel('Caudal [m^3/s]'); title('u(t)=2 sin(0.4t)+2.5');
subplot(2,1,2); plot(t,hvar,'LineWidth',1.5); grid on;
xlabel('Tiempo [s]'); ylabel('Nivel [m]');
title('Niveles con condiciones iniciales nulas'); legend('h_1','h_2','Location','best');
print(gcf,fullfile(salida,'p2_variable.png'),'-dpng','-r200');
exportgraphics(gcf,fullfile(salida,'p2_variable.pdf'),'ContentType','vector');
% impulse(sys) agrega agua: para una purga usar una perturbacion negativa
% respecto de un equilibrio con suficiente agua; no niveles absolutos negativos.
writetable(table(t,u,hstep(:,1),hstep(:,2),himp(:,1),himp(:,2), ...
    hvar(:,1),hvar(:,2),'VariableNames',{'t','u','step_h1','step_h2', ...
    'impulse_h1','impulse_h2','variable_h1','variable_h2'}), ...
    fullfile(salida,'p2_datos.csv'));
save(fullfile(salida,'p2_resultados.mat'),'A','B','C','D','t','u', ...
    'hstep','himp','hvar','K','err_step','err_imp','err_var');
