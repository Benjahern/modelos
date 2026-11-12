% PREGUNTA2 Dos estanques: step, impulse y lsim.
% Requiere MATLAB R2025a o posterior y Control System Toolbox.
% Escalon y seno desde reposo; purga desde el equilibrio con caudal de 1 m^3/s.
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
% Purga ideal: extraer Vp=1 m^3 del estanque 1 en t=0.
% Se mantiene el caudal base ub=1 m^3/s durante la recuperacion.
ub=1; Vp=1; hbase=-(A\B)*ub;
sys_purga=ss(A,-B,C,D); % Entrada positiva = caudal extraido.
[dhpurga,~]=impulse(sys_purga,t); dhpurga=Vp*squeeze(dhpurga);
hpurga=dhpurga+hbase'; % Niveles absolutos, no desviaciones negativas.
u=2*sin(0.4*t)+2.5;
[hvar,~]=lsim(sys,u,t,x0);

% Referencias exactas; z es la respuesta fasorial a frecuencia 0.4 rad/s.
K=-A\B; w=0.4; z=(1i*w*eye(2)-A)\B;
hstep_teo=zeros(length(t),2); dhpurga_teo=hstep_teo; hvar_teo=hstep_teo;
for k=1:length(t)
    E=expm(A*t(k));
    hstep_teo(k,:)=(K-E*K)';
    dhpurga_teo(k,:)=(-E*B*Vp)';
    hvar_teo(k,:)=(2.5*K+2*imag(z*exp(1i*w*t(k))) ...
        -E*(2.5*K+2*imag(z)))';
end
err_step=max(max(abs(hstep-hstep_teo)));
err_imp=max(max(abs(dhpurga-dhpurga_teo)));
err_var=max(max(abs(hvar-hvar_teo)));
% lsim interpola la entrada muestreada: admitir ese error, no redondearlo a 0.
assert(err_step<1e-8 && err_imp<1e-8,'Fallo en step o impulse.');
assert(err_var<1e-4,'Error excesivo de lsim respecto al seno continuo.');
assert(all(real(eig(A))<0),'El modelo debe ser estable.');
assert(max(abs(K-[4/3;2/3]))<1e-12,'Equilibrio incorrecto.');
assert(min(hpurga(:))>=0,'La purga no debe producir niveles negativos.');
assert(abs([A1 A2]*dhpurga(1,:)'+Vp)<1e-12,'Volumen extraido incorrecto.');
assert(max(abs(hpurga(1,:)'-(hbase-B*Vp)))<1e-12,'Salto inicial incorrecto.');
fprintf('Matrices del modelo:\n'); disp(A); disp(B); disp(C); disp(D);
fprintf('Polos:\n'); disp(eig(A));
fprintf('Equilibrio para u=1 [m]:\n'); disp(K);
fprintf('Errores maximos [step purga lsim]: %.3e %.3e %.3e\n', ...
    err_step,err_imp,err_var);
fprintf('Niveles a los 25 s [m]:\n'); disp([hstep(end,:);hpurga(end,:);hvar(end,:)]);

figure('Theme','light','Color','w','Position',[100 100 900 750]);
subplot(2,1,1); plot(t,hstep,'LineWidth',1.5); grid on;
xlabel('Tiempo [s]'); ylabel('Nivel [m]'); title('Escalon: u=1 m^3/s');
legend('h_1','h_2','Location','best');
subplot(2,1,2); plot(t,hpurga,'LineWidth',1.5); grid on;
xlabel('Tiempo [s]'); ylabel('Nivel [m]');
title('Purga de 1 m^3 del estanque 1; caudal base de 1 m^3/s');
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
writetable(table(t,u,hstep(:,1),hstep(:,2),hpurga(:,1),hpurga(:,2), ...
    hvar(:,1),hvar(:,2),'VariableNames',{'t','u','step_h1','step_h2', ...
    'purga_h1','purga_h2','variable_h1','variable_h2'}), ...
    fullfile(salida,'p2_datos.csv'));
save(fullfile(salida,'p2_resultados.mat'),'A','B','C','D','t','u', ...
    'hstep','dhpurga','hpurga','hbase','ub','Vp','hvar','K','err_step','err_imp','err_var');
