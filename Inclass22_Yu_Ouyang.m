%Inclass 22

%1. Consider the case of the auto-activating gene that we discussed in class
%today. Make a bifurcation diagram for this system by varying the
%activated transcription rate for three cases - in which 1, 4, or 8 copies of the
%transcripton factor are necessary to activate transciption. Comment on your
%results. 

%%
kb = 3; ku = 0; x0 = 0.6;
rhs1 = @(t,x) (ku+kb*x)./(1+x)-x;
rhs4 = @(t,x) (ku+kb*x.^4)./(1+x.^4)-x;
rhs8 = @(t,x) (ku+kb*x.^8)./(1+x.^8)-x;
sol1 = ode23(rhs1,[0 10],x0);
sol4 = ode23(rhs4,[0 10],x0);
sol8 = ode23(rhs8,[0 10],x0);
plot(sol1.x,sol1.y);hold on;
plot(sol4.x,sol4.y);hold on;
plot(sol8.x,sol8.y);
xlabel('Time'); ylabel('Expression');
legend({'n=1','n=4','n=8'},'FontSize',12);

%% This is the bifuracation diagram codes
figure; hold on;
ku = 0;
for kb = 0:0.05:5
    % n = 1 case: x.^2 + (1-kb)x - ku = 0
    polycoeff1 = [1 1-kb -ku];
    rts1 = roots(polycoeff1);
    rts1 = rts1(imag(rts1) == 0);
   
    % n = 4 case: x.^5 - kb*x.^4 + x - ku = 0
    polycoeff4 = [1 -kb 0 0 1 -ku];
    rts4 = roots(polycoeff4);
    rts4 = rts4(imag(rts4) == 0);

    % n = 8 case: x.^9 - kb*x.^8 + x - ku = 0
    polycoeff8 = [1 -kb 0 0 0 0 0 0 1 -ku];
    rts8 = roots(polycoeff8);
    rts8 = rts8(imag(rts8) == 0);
    
    plot(kb*ones(length(rts1),1),rts1, 'r.');
    plot(kb*ones(length(rts4),1),rts4, 'k.');
    plot(kb*ones(length(rts8),1),rts8, 'b.');
end
hold off;
xlabel('k_b');ylabel('Fixed points');
legend({'n=1','n=4','n=8'},'FontSize',12);
set(gca,'FontSize',24);
    

% Yu Ouyang's Comment: There is a saddle node in n=4 and n=8 bifurcation
% diagram as expected.

% 2. Make a similar diagram for the case of an autorepressing gene in the
% case that 2 copies are need to turn off the gene. 
%% 
% In this case, we would flip the ku and kb in function, because the unbinding state
% would trigger the transcription. 
kb = 3; ku = 0; x0 = 0.3;
rhs1 = @(t,x) (kb+ku*x)./(1+x)-x;
rhs4 = @(t,x) (kb+ku*x.^4)./(1+x.^4)-x;
rhs8 = @(t,x) (kb+ku*x.^8)./(1+x.^8)-x;
sol1 = ode23(rhs1,[0 10],x0);
sol4 = ode23(rhs4,[0 10],x0);
sol8 = ode23(rhs8,[0 10],x0);
plot(sol1.x,sol1.y);hold on;
plot(sol4.x,sol4.y);hold on;
plot(sol8.x,sol8.y);
xlabel('Time'); ylabel('Expression');
legend({'n=1','n=4','n=8'},'FontSize',12);
%% This is the bifuracation diagram codes
figure; hold on;
ku = 0;
for kb = 0:0.05:5
    % n = 1 case: x.^3 - ku*x.^2 + x - kb = 0
    polycoeff2 = [1 1-ku 1 -kb];
    rts2 = roots(polycoeff2);
    rts2 = rts2(imag(rts2) == 0);
    plot(kb*ones(length(rts2),1),rts2, 'r.');
end
hold off;
xlabel('k_b');ylabel('Fixed points');
set(gca,'FontSize',24);
   


