clc;clear;tic;close all

%% Input initialization
% load the trained machine learning model
% ite is iteration number
% increment is weight increment for objective function
% n is concrete mixture input feature
% ub and lb is concrete mixture input upper and lower bound

% Price component is edited in fobj.m file

%% PF-SOS Initialization
load('MdlCRF.mat')
ite=40;
increment=0.015;
n=7;
lb=[0.2900    1.7000    10    14.0000   14.0000        0    0];
ub=[0.7200    6.4000   90.0000   32.0000   32.0000   10    3];

%% Initial Organism
nd=0:increment:1;
ratiolist=[nd' 1-(nd')];
ecosize=length(ratiolist);
eco=zeros(ecosize,n);
fitness=zeros(ecosize,1);
for i=1:ecosize
    ratio=ratiolist(i,:);
    eco(i,:)=rand(1,n).*(ub-lb)+lb;
    f=fobj(eco(i,:),MdlC);
    fitness(i,:)=sum(f.*ratio);

    if i==1
        PF(i,:)=f;
        PFeco(i,:)=eco(i,:);
    else
        if isempty(find(PF(:,1)>f(1)))==0 && isempty(find(PF(:,2)>f(2)))==0
            PF=vertcat(PF,f);
            PFeco=vertcat(PFeco,eco(i,:));
        end
    end
end

%% Start Iteration
for h=1:ite
    fprintf('Iteration %d \n',h)
    for i=1:ecosize
        ratio=ratiolist(i);
        % Update the best Organism
        if mod(h,2)==0 && h<0
            bestOrganism=eco(i,:);
        else
            bestOrganism=PFeco(ceil(unifrnd(0,size(PF,1))),:);
        end

        %% Commensialism Phase
        j=i;
        while i==j
            seed=randperm(ecosize);
            j=seed(1);
        end
        % Calculate new solution after Commensalism Phase
        ecoNew1=eco(i,:)+(rand(1,n)*2-1).*(bestOrganism-eco(j,:));
        ecoNew1=bound(ecoNew1,ub,lb);
        % Evaluate the fitness of the new solution
        f=fobj(ecoNew1,MdlC);

        if isempty(find(PF(:,1)>f(1)))==0 || isempty(find(PF(:,2)>f(2)))==0
            PF=vertcat(PF,f);
            PFeco=vertcat(PFeco,ecoNew1);
            DOMINATED  = checkDomination(PF);
            PF=PF(~DOMINATED,:);
            PFeco=PFeco(~DOMINATED,:);
            [~,b]=sort(PF(:,1));
            PF=PF(b,:);PFeco=PFeco(b,:);
            recPF{h}=[PF PFeco];
        end

        % Accept the new solution if the fitness is better
        for k=1:ecosize
            ratio=ratiolist(k,:);
            y=sum(f.*ratio);
            if y<fitness(k)
                fitness(k)=y;
                eco(k,:)=ecoNew1;
            end
        end
        % End of Commensalism Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Parasitism Phase
        j=i;
        while i==j
            seed=randperm(ecosize);
            j=seed(1);
        end
        % Determine Parasite Vector & Calculate the fitness
        parasiteVector=eco(i,:);
        seed=randperm(n);
        pick=seed(1:ceil(rand*n));  % select random dimension
        parasiteVector(:,pick)=rand(1,length(pick)).*(ub(pick)-lb(pick))+lb(pick);
        f=fobj(parasiteVector,MdlC);

        if isempty(find(PF(:,1)>f(1)))==0 || isempty(find(PF(:,2)>f(2)))==0
            PF=vertcat(PF,f);
            PFeco=vertcat(PFeco,parasiteVector);
            DOMINATED  = checkDomination(PF);
            PF=PF(~DOMINATED,:);
            PFeco=PFeco(~DOMINATED,:);
            [~,b]=sort(PF(:,1));
            PF=PF(b,:);PFeco=PFeco(b,:);
            recPF{h}=[PF PFeco];
        end

        % Kill organism j and replace it with the parasite
        % if the fitness is lower than the parasite
        for k=1:ecosize
            ratio=ratiolist(k,:);
            y=sum(f.*ratio);
            if y<fitness(k)
                fitness(k)=y;
                eco(k,:)=parasiteVector;
            end
        end
        % End of Parasitism Phase
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    DOMINATED  = checkDomination(PF);
    PF=PF(~DOMINATED,:);
    PFeco=PFeco(~DOMINATED,:);
    [~,b]=sort(PF(:,1));
    PF=PF(b,:);PFeco=PFeco(b,:);
    recPF{h}=[PF PFeco];

    clf
    hold on
    scatter(100./PF(:,1),PF(:,2),'filled')
    grid on
    hold off
    xlabel('Comp Strength')
    ylabel('price')
    pause(0.1)
end

function dom_vector = checkDomination(fitness)
Np = size(fitness,1);
dom_vector = zeros(Np,1);
all_perm = nchoosek(1:Np,2);    % Possible permutations
all_perm = [all_perm; [all_perm(:,2) all_perm(:,1)]];

d = dominates(fitness(all_perm(:,1),:),fitness(all_perm(:,2),:));
dominated_particles = unique(all_perm(d==1,2));
dom_vector(dominated_particles) = 1;
end

% Function that returns 1 if x dominates y and 0 otherwise
function d = dominates(x,y)
d = all(x<=y,2) & any(x<=y,2);
end