function ulimits(txtname,increment,range_l,range_u)

% txtname:
%       name of the text file that contains the input data
% increment:
%       resolution of the temperature range (e.g. 1K, 10K, 50K, 100K)
% range_l:
%       lower limit of the temperature range of evaluation
% range_u
%       upper limit of the temperature range of evaluation

useall='n';
while (isequal(useall,'n')==1)

    %%%Reading data
    [paper_type,paper_ID,considered,Tlower,Tupper,A,n,E_per_R,k,order,f_type,...
        y1,y2,y3,Tmid,multiplier,A_bi,n_bi,E_per_R_bi,T0,interval]=...
        textread(txtname,' %s %s  %f  %f %f %f  %f %f %f%f%f%f%f%f%f%f%f%f%f%f%f'...
        ,'emptyvalue',0);

    % paper_type:
    %       possible values:    Review,Experiment,Theory
    % paper_ID:
    %       NIST squib or similar identifier
    % considered:
    %       0       these data are not used
    %       1       these data are used, and using single Arrhenius expression
    %       2       these data are used, and using bi-Arrhenius expression
    % Tlower:
    %       lower limit of temperature range
    % Tupper:
    %       upper limit of tempreature range
    %       examples:
    %               300 2000   (temperature range 300-2000K)
    %               298 298    (data for T=298K only)
    %               298 0      (data for T=298K only)
    %
    % A,n,E_per_R
    %       Arrhenius parameters A, n, E/R
    % k:
    %       rate coefficient at 298K
    % order:
    %       reaction order (1, 2 or 3)
    % f_type:
    %       0       uncertainty parameter f is not stated in the paper
    %       1       temperature independent parameter f
    %       2       uncertainty parameter f is stated at Tlower and Tupper
    %       3       uncertainty parameter f is stated at Tlower, Tmid and Tupper
    % y1:
    %       parameter f at Tlower
    % y2:
    %       parameter f at Tmid (if applicable)
    % y3:
    %       parameter f at Tupper (if applicable)
    % Tmid:
    %       middle temperature value where parameter f is available
    %
    % multiplier:
    %       use when different third body efficiencies given
    %
    %
    % N:
    %       number of lines containing data (including the middle line)
    %       data for the middle line should be in the last row
    % A_bi,n_bi,E_per_R_bi
    %       bi-Arrhenius parameters A, n, E/R (if you would like to use single
    % Arrhenius expression, these values are zero)
    %

    [N,M]=size(Tlower);
    
    CompareCount = 0;
    for i=1:N
        if(strcmp(paper_type(i),'Compare') == 1)
            CompareCount = CompareCount + 1;
            continue;
        end

        if(CompareCount > 0)
            CompareCount = CompareCount + 1;
        end
    end

    N = N - CompareCount;%如果compare有3个就减3

    if(CompareCount > 0)
        paper_ID = paper_ID(1:N);
        paper_type = paper_type(1:N);
    end

    for i=1:N

        if (isequal(considered(i),0)==0 && isequal(considered(i),1)==0 && isequal(considered(i),2)==0)
            disp('Please check your text file, the flag is not correct in column 3 line')
            disp(i)
            return
        end

        if (isequal(order(i),1)==0 && isequal(order(i),2)==0 && isequal(order(i),3)==0)

            disp('Please check your text file, the order of reaction is not correct in column 10 line')
            disp(i)
            return
        end

        if (isequal(considered(i),1)==1 && A_bi(i)~=0 || isequal(considered(i),1)==1 && n_bi(i)~=0 || isequal(considered(i),1)==1 && E_per_R_bi(i)~=0)
            disp('Please check your text file,you use single Arrhenius expression but you defined second set of Arrhenius parameters in line')
            disp(i)
            return
        end


        if (isequal(considered(i),2)==1 && A_bi(i)==0 && n_bi(i)==0  && E_per_R_bi(i)==0)
            disp('Please check your text file,you use double Arrhenius expression but you did not define second set of Arrhenius parameters in line')
            disp(i)
            return
        end
        
        if (isequal(considered(i),2)==1 && T0(i)~=0)
            disp('Please check your text file, extended Arrhenius expression with T0 is not available in the case of duplicate Arrhenius expression')
            disp(i)
            return
        end
        

    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % definition of paper_type: 0 means identical type to the previous one

    v=ismember(paper_type,'0');

    [P,Q]=size(paper_type);

    for i=1:P-1

        if(v(i,1)==0 && v(i+1,1)==1)
            paper_type(i+1)=paper_type(i);
            v(i+1,1)=0;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Tlower(N)=range_l;
    Tupper(N)=range_u;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % A tempreature grid is defined by Tlower, Tupper, and increment. At each
    % temperatures values the rate coefficient is calculated from the various
    % Arrhenius expressions within there own temperature range


    %%%%%%%%%%%%%%%%

    for i=1:N
        if (Tupper(i)==0)
            Tupper(i)=Tlower(i);
        end
    end
    %%%%%%%%%%%%%%%


    gridpoints=zeros(1,N);
    for i=1:N

        if  (mod(Tlower(i)/increment,1)==0 || mod(Tupper(i)/increment,1)==0)

            gridpoints(i)=floor((Tupper(i)-Tlower(i))/increment)+1;

        end

        if((Tupper(i)-Tlower(i))/increment>floor((Tupper(i)-Tlower(i))/increment)&&...
                mod(Tlower(i)/increment,1)~=0 &&  mod(Tupper(i)/increment,1)~=0)

            gridpoints(i)=floor((Tupper(i))/increment)-floor((Tlower(i))/increment);
        end

        if( (Tupper(i)-Tlower(i))/increment<=floor((Tupper(i)-Tlower(i))/increment)&&...
                mod(Tlower(i)/increment,1)~=0 &&  mod(Tupper(i)/increment,1)~=0)

            gridpoints(i)=floor((Tupper(i)-Tlower(i))/increment);
        end
    end



    inside_inc=ceil((range_u-range_l)/increment);
    Xcoord=zeros(N,inside_inc+2);
    Xcoord_gp=zeros(N,max(gridpoints)+1);

    for i=1:N
        for j=1:inside_inc+2
            Xcoord(i,inside_inc+3-j)=1/(1/Tupper(i)+(j-1)*(1/Tlower(i)-1/Tupper(i))/(inside_inc+1));
        end
        if (gridpoints(i)~=0)
            for k=0:gridpoints(i)-1
                Xcoord_gp(i,k+1)=ceil(Tlower(i)/increment)*increment+k*increment;
            end
        end

        if(gridpoints(i)==0 && mod(Tlower(i)/increment,1)==0 || mod(Tupper(i)/increment,1)==0)
            Xcoord_gp(i,1)=ceil(Tlower(i)/increment)*increment;
        end


    end


    % Evaluation of the rate coefficients at the grid points

    % Xcoord:
    %           temperature values at the grid points for each data set
    % Ycoord:
    %           rate coefficient values at the grid points for each data set

    Ycoord=zeros(N,inside_inc+2);
    Ycoord_gp=zeros(N,max(gridpoints)+1);

    for i=1:N
        if considered(i)~=2
            for k=1:inside_inc+2
                if(Xcoord(i,k)~=0)
                    Ycoord(i,k)=(log(A(i))+n(i)*log(Xcoord(i,k))-...
                        E_per_R(i)*(Xcoord(i,k)+T0(i))*1/(Xcoord(i,k)^2+T0(i)^2)+log(1/multiplier(i)))./log(10);
                end
            end
        end
        if considered(i)==2
            for  k=1:inside_inc+2
                if(Xcoord(i,k)~=0)
                    E_min=min(E_per_R(i),E_per_R_bi(i));

                    Ycoord(i,k)=(-E_min*1/Xcoord(i,k)+log(A(i)*Xcoord(i,k)^(n(i))*...
                        exp((-E_per_R(i)+E_min)*1/Xcoord(i,k))+A_bi(i)*Xcoord(i,k)^(n_bi(i))*...
                        exp((-E_per_R_bi(i)+E_min)*1/Xcoord(i,k))))./log(10)+log(1/multiplier(i))./log(10);

                end
            end
        end
    end


    [Xko_size_1,Xko_size_2]=size(Xcoord_gp);

    for i=1:Xko_size_1
        if considered(i)~=2
            for k=1:Xko_size_2
                if(Xcoord_gp(i,k)~=0)
                    Ycoord_gp(i,k)=(log(A(i))+n(i)*log(Xcoord_gp(i,k))-...
                        E_per_R(i)*(Xcoord_gp(i,k)+T0(i))*1/(Xcoord_gp(i,k)^2+T0(i)^2)+log(1/multiplier(i)))./log(10);
                end
            end
        end
        if considered(i)==2
            for  k=1:Xko_size_2
                if(Xcoord_gp(i,k)~=0)
                    E_min=min(E_per_R(i),E_per_R_bi(i));

                    Ycoord_gp(i,k)=(-E_min*1/Xcoord_gp(i,k)+log(A(i)*Xcoord_gp(i,k)^(n(i))*...
                        exp((-E_per_R(i)+E_min)*1/Xcoord_gp(i,k))+A_bi(i)*Xcoord_gp(i,k)^(n_bi(i))*...
                        exp((-E_per_R_bi(i)+E_min)*1/Xcoord_gp(i,k))))./log(10)+log(1/multiplier(i))./log(10);

                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Adding uncertainty limits if needed


    additional=zeros(N,inside_inc+2);
    Xadded=zeros(N,inside_inc+2);

    additional_gp=zeros(N,max(gridpoints)+1);
    Xadded_gp=zeros(N,max(gridpoints)+1);

    for i=1:N
        if (f_type(i)==1 || f_type(i)==2||f_type(i)==3)
            Xadded(i,:)=Xcoord(i,:);
            for k=1:inside_inc+2
                if (f_type(i)==1)
                    additional(i,k)=y1(i);
                end

                if (f_type(i)==2)
                    additional(i,k)=(1000/Tlower(i)-1000/Xadded(i,k))/(1000./Tlower(i)-1000./Tupper(i))*...
                        (y3(i)-y1(i))+y1(i);
                end

                if (f_type(i)==3)
                    if (1000/Xadded(i,k)<=1000./Tmid(i))
                        additional(i,k)=(1000/Tupper(i)-1000/Xadded(i,k))/(1000./Tmid(i)-1000./Tupper(i))*...
                            (y3(i)-y2(i))+y3(i);
                    end

                    if (1000/Tmid(i)<=1000/Xadded(i,k))
                        additional(i,k)=(1000/Tmid(i)-1000/Xadded(i,k))/(1000/Tlower(i)-1000/Tmid(i))*...
                            (y2(i)-y1(i))+y2(i);
                    end
                end
            end
        end

    end


    for i=1:N

        if (f_type(i)==1 || f_type(i)==2||f_type(i)==3)


            Xadded_gp(i,:)=Xcoord_gp(i,:);


            for k=1:gridpoints(i)+1
                if (f_type(i)==1)
                    additional_gp(i,k)=y1(i);
                end

                if (f_type(i)==2)
                    additional_gp(i,k)=(1000/Tlower(i)-1000/Xadded_gp(i,k))/(1000./Tlower(i)-1000./Tupper(i))*...
                        (y3(i)-y1(i))+y1(i);
                end

                if (f_type(i)==3)

                    if (1000/Xadded_gp(i,k)<=1000./Tmid(i))
                        additional_gp(i,k)=(1000/Tupper(i)-1000/Xadded_gp(i,k))/(1000./Tmid(i)-1000./Tupper(i))*...
                            (y3(i)-y2(i))+y3(i);
                    end

                    if (1000/Tmid(i)<=1000/Xadded_gp(i,k))
                        additional_gp(i,k)=(1000/Tmid(i)-1000/Xadded_gp(i,k))/(1000/Tlower(i)-1000/Tmid(i))*...
                            (y2(i)-y1(i))+y2(i);
                    end
                end


            end
        end

    end





    for i=1:N

        Xaddeddouble(2*i,:)=Xadded(i,:);
        Xaddeddouble(2*i-1,:)=Xadded(i,:);
        Xaddeddouble_gp(2*i,:)=Xadded_gp(i,:);
        Xaddeddouble_gp(2*i-1,:)=Xadded_gp(i,:);

        paper_ID_2(2*i)=strcat(paper_ID(i),'+f');
        paper_ID_2(2*i-1)=strcat(paper_ID(i),'-f');
        paper_ID_2_gp(2*i)=strcat(paper_ID(i),'+f');
        paper_ID_2_gp(2*i-1)=strcat(paper_ID(i),'-f');

        paper_type_2(2*i)=paper_type(i);
        paper_type_2(2*i-1)=paper_type(i);
        paper_type_2_gp(2*i)=paper_type(i);
        paper_type_2_gp(2*i-1)=paper_type(i);

        Yaddeddouble(2*i,:)=Ycoord(i,:)+additional(i,:);
        Yaddeddouble(2*i-1,:)=Ycoord(i,:)-additional(i,:);
        Yaddeddouble_gp(2*i,:)=Ycoord_gp(i,:)+additional_gp(i,:);
        Yaddeddouble_gp(2*i-1,:)=Ycoord_gp(i,:)-additional_gp(i,:);

        numbering(2*i,1)=i;
        numbering(2*i-1,1)=i;
        numbering_gp(2*i,1)=i;
        numbering_gp(2*i-1,1)=i;

        A_2(2*i,:)=A(i,:);
        A_2(2*i-1,:)=A(i,:);
        A_2_gp(2*i,:)=A(i,:);
        A_2_gp(2*i-1,:)=A(i,:);

        n_2(2*i,:)=n(i,:);
        n_2(2*i-1,:)=n(i,:);
        n_2_gp(2*i,:)=n(i,:);
        n_2_gp(2*i-1,:)=n(i,:);

        E_per_R_2(2*i,:)=E_per_R(i,:);
        E_per_R_2(2*i-1,:)=E_per_R(i,:);
        E_per_R_2_gp(2*i,:)=E_per_R(i,:);
        E_per_R_2_gp(2*i-1,:)=E_per_R(i,:);

        Tlower_2(2*i,:)=Tlower(i,:);
        Tlower_2(2*i-1,:)=Tlower(i,:);
        Tlower_2_gp(2*i,:)=Tlower(i,:);
        Tlower_2_gp(2*i-1,:)=Tlower(i,:);

        Tupper_2(2*i,:)=Tupper(i,:);
        Tupper_2(2*i-1,:)=Tupper(i,:);
        Tupper_2_gp(2*i,:)=Tupper(i,:);
        Tupper_2_gp(2*i-1,:)=Tupper(i,:);



    end



    ZZ=ismember(f_type,[1,2,3]);



    cleaning=[];
    for i=1:N
        if (ZZ(i,1)==0)

            cleaning=horzcat(cleaning,2*i,2*i-1);
        end
    end

    Xaddeddouble(cleaning,:)=[];
    Xaddeddouble_gp(cleaning,:)=[];

    Yaddeddouble(cleaning,:)=[];
    Yaddeddouble_gp(cleaning,:)=[];

    numbering(cleaning,:)=[];
    numbering_gp(cleaning,:)=[];

    Tlower_2(cleaning,:)=[];
    Tlower_2_gp(cleaning,:)=[];

    Tupper_2(cleaning,:)=[];
    Tupper_2_gp(cleaning,:)=[];

    A_2(cleaning,:)=[];
    A_2_gp(cleaning,:)=[];

    n_2(cleaning,:)=[];
    n_2_gp(cleaning,:)=[];

    E_per_R_2(cleaning,:)=[];
    E_per_R_2_gp(cleaning,:)=[];

    paper_ID_2(cleaning)=[];
    paper_ID_2_gp(cleaning)=[];

    paper_type_2(cleaning)=[];
    paper_type_2_gp(cleaning)=[];

    numbering=num2str(numbering);
    numbering_gp=num2str(numbering_gp);


    % Elimination of data that are not considered


    zz=linspace(1,N,N);

    del=[];

    for i=1:N
        if (considered(i)==0)
            del=horzcat(del,i);
        end
    end


    Xcoord(del,:)=[];
    Xcoord_gp(del,:)=[];

    Ycoord(del,:)=[];
    Ycoord_gp(del,:)=[];
    Ycoord;
    Tlower(del,:)=[];
    Tupper(del,:)=[];
    paper_type(del)=[];
    paper_ID(del)=[];
    A(del,:)=[];
    A_bi(del,:)=[];
    n(del,:)=[];
    n_bi(del,:)=[];
    E_per_R(del,:)=[];
    E_per_R_bi(del,:)=[];
    zz(del)=[];
    zz=num2str(zz');
    considered(del,:)=[];

    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    size(Xaddeddouble_gp);
    size(Xcoord_gp);

    Xnew=vertcat(Xaddeddouble,Xcoord);
    Xnew_gp=vertcat(Xaddeddouble_gp,Xcoord_gp);

    Ynew=vertcat(Yaddeddouble,Ycoord);
    Ynew_gp=vertcat(Yaddeddouble_gp,Ycoord_gp);

    numbering_new=strvcat(numbering,zz);
    paper_type_new=vertcat(paper_type_2',paper_type);
    paper_ID_new=vertcat(paper_ID_2',paper_ID);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%% narrowing the temperature range

    Xrange=0.*Xnew;
    Xrange_gp=0.*Xnew_gp;

    Yrange=0.*Ynew;
    Yrange_gp=0.*Ynew_gp;

    [Xsize_1,Xsize_2]=size(Xnew);
    [Xsize_1_gp,Xsize_2_gp]=size(Xnew_gp);


    for i=1:Xsize_1
        for j=1:Xsize_2
            if (range_l<=Xnew(i,j) && Xnew(i,j)<=range_u)
                Xrange(i,j)=Xnew(i,j);
                Yrange(i,j)=Ynew(i,j);
            end
        end
    end


    for i=1:Xsize_1_gp
        for j=1:Xsize_2_gp
            if (range_l<=Xnew_gp(i,j) && Xnew_gp(i,j)<=range_u)
                Xrange_gp(i,j)=Xnew_gp(i,j);
                Yrange_gp(i,j)=Ynew_gp(i,j);
            end
        end
    end

    [tt_Xnew,tt_Ynew,tt_numbering_new,tt_paper_type_new,tt_paper_ID_new] = taotao(txtname,increment,range_l,range_u);
    full_paper_type_new = vertcat(paper_type_new(:,1),tt_paper_type_new);
    full_paper_ID_new = vertcat(paper_ID_new,tt_paper_ID_new);

    p1 = plot(1000./Xnew(1:Xsize_1-1,:)',Ynew(1:Xsize_1-1,:)','--.','LineWidth',0.2);
    xlabel('1000/T')
    ylabel('log_{10}(k)')
    grid on
    hold on

    p2 = plot(1000./Xnew(Xsize_1,:),Ynew(Xsize_1,:),'-r','Linewidth',1.5);
    title('Arrhenius plot')
    hold on
    [tt_Xnew,tt_Ynew,tt_numbering_new,tt_paper_type_new,tt_paper_ID_new] = taotao(txtname,increment,range_l,range_u);
    p3 = plot(1000./tt_Xnew(1,:),tt_Ynew(1,:),'Color',[0.1 0.8 0.9] ,'Linewidth',1.5);%'Color',[0.1 0.6 0.1] NPU;'Color',[0.1 0.8 0.9] Taotao; '-m' henry
    hold on
    p4 = plot(1000./tt_Xnew(2,:),tt_Ynew(2,:),'-m','Linewidth',1.5);
    hold on
    p5 = plot(1000./tt_Xnew(3,:),tt_Ynew(3,:),'Color',[0.1 0.6 0.1],'Linewidth',1.5);
    hold on

    %有4个compare的时候才画这个，一般都是3个
    if(CompareCount == 4)
        p6 = plot(1000./tt_Xnew(4,:),tt_Ynew(4,:),'Color',[0.6 0.3 0.4],'Linewidth',1.5);
        p = [p1' p2 p3 p4 p5 p6];
    else
        p = [p1' p2 p3 p4 p5];
    end
    
    title('Arrhenius plot')

    gridLegend(p,2,strcat(full_paper_type_new,';',full_paper_ID_new),'location','EastOutside');
    hold on

    %plot(1000./Xnew_gp(end,:),Ynew_gp(end,:),'-r','Linewidth',1.5)
    %title('Arrhenius plot')
    %hold on
    %legend(strcat(full_paper_type_new,';',full_paper_ID_new),'location','EastOutside')

    hold on

    first_answer=[];
    while(isequal(first_answer,'y')==0 && isequal(first_answer,'n')==0)
        [useall]=input('Do you want to use all plotted data for the determination of the uncertainty band? (y/n)','s');
        first_answer=useall;
    end
    if (isequal(useall,'n')==1)
        close(figure(1));


        disp('Press  ENTER if you finished the modification of the input data file!')
        pause
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 deteremination of f(T) points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Division of the temperature range according to the temperature increment



if  (mod(range_l/increment,1)==0 || mod(range_u/increment,1)==0)

    num_of_gridpoints=floor((range_u-range_l)/increment)+1;

end

if((range_u-range_l)/increment>floor((range_u-range_l)/increment)&&...
        mod(range_l/increment,1)~=0 &&  mod(range_u/increment,1)~=0)

    num_of_gridpoints=floor((range_u-range_l)/increment)+1;
end

if((range_u-range_l)<=floor((range_u-range_l)/increment)&&...
        mod(range_l/increment,1)~=0 &&  mod(range_u/increment,1)~=0)

    num_of_gridpoints=floor((range_u-range_l)/increment);
end


if (num_of_gridpoints~=0)
    for k=0:num_of_gridpoints-1
        T_grid(k+1)=ceil(range_l/increment)*increment+k*increment;
    end
end

if(num_of_gridpoints==0 && mod(range_l/increment,1)==0 || mod(range_u/increment,1)==0)
    T_grid(1)=ceil(range_l/increment)*increment;
end






[T_size1,T_size2]=size(T_grid);

distance=zeros(Xsize_1_gp,T_size2);

for i=1:T_size2

    for j=1:Xsize_1_gp
        for k=1:Xsize_2_gp
            if (Xrange_gp(j,k)==T_grid(i))
                distance(j,i)=Yrange_gp(j,k);
            end
        end
    end
end

recommended=Yrange_gp(Xsize_1_gp,1:T_size2);
F=ones(Xsize_1_gp,1)*Yrange_gp(Xsize_1_gp,1:T_size2);

for i=1:Xsize_1_gp
    for j=1:T_size2
        if(distance(i,j)==0)
            F(i,j)=0;
        else
        end

    end
end


C=(F-distance);


%取距离最远的的作为f_origianl
for j=1:T_size2
    f_values(j)=max(abs(C(:,j)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(1000./T_grid,recommended+f_values,'g*')
hold on
plot(1000./T_grid,recommended-f_values,'g*')
hold on
xlabel('1000/T')
ylabel('log_{10}(k)')
title('Arrhenius plot with main +/- f(T)_{original} points added')
axis([1000/range_u*0.9 1000/range_l*1.1 min(recommended-f_values)-1 max(recommended+f_values)+1])
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






% gridpoints without experimental/theoretical information are not considered
inside=[];

k=1;
for i=1:num_of_gridpoints

    if (f_values(i)~=0)
        inside(k)=i;
        k=k+1;
    end
end

[CC,DD]=size(inside);

% Presence of temperature range without experimental/theoretical information is
% indicated

if (num_of_gridpoints~=DD)
    disp('Temperature range without experimental/theoretical information')
end


figure(2)
plot(T_grid(inside),f_values(inside),...
    'om',...
    'LineWidth',1,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',6)


xlabel('T ')
ylabel('f(T)')
title('T-f_{original}(T) values (before modification)')
grid on
hold off


format long
horzcat(T_grid(inside)',f_values(inside)');




outp_1=fopen('T_f_points.txt','wt+');
fprintf(outp_1,'%d\t%d\n',horzcat(T_grid(inside)',f_values(inside)')');
fclose(outp_1);



modif='y';
while (isequal(modif,'y')==1)


    [Xcorr,Ycorr]=...
        textread('T_f_points.txt',' %f %f');


    [items,items2]=size(Xcorr);
    figure(3)
    plot(Xcorr,Ycorr,...
        'om',...
        'LineWidth',1,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','b',...
        'MarkerSize',6)


    xlabel('T ')
    ylabel('f(T)')
    title('T-f_{original}(T) values (after modification)')
    grid on
    hold on





    second_answer=[];
    while(isequal(second_answer,'y')==0 && isequal(second_answer,'n')==0)
        [modif]=input('Would you like to modify the T-f(T) values? (y/n)','s');
        second_answer=modif;
    end
    if (isequal(modif,'y')==1)
        close(figure(3));
        disp('Press ENTER if you have finished the modification of T_f_points.txt!')
        pause
    end

end


third_answer=[];
while(isequal(third_answer,'U')==0 && isequal(third_answer,'N')==0)
    [distribut]=input('Do you assume Uniform or Normal distribution? (U,N)','s');
    third_answer=distribut;
end
%Assuming uniform distribution
if  (ismember(distribut,'U')==1)
    close(figure(3));
    [constant]=input('Please give the assumed constant f value: ');

    sigmaalpha=constant/3;



    coma=fopen('CovMatrix.txt','wt+');
    fprintf(coma,'* Assumed distribution\n');
    fprintf(coma,'uniform\n');
    fprintf(coma,'* Uncertainty type\n');
    fprintf(coma,'3slog10k\n');
    fprintf(coma,'* covariance matrix [(a,n,e)x(a,n,e)], where a=lnA, e=E/R:\n');
    fprintf(coma,'%f\t%f\t%f\n',(sigmaalpha*log(10))^2,0,0);
    fprintf(coma,'%f\t%f\t%f\n',0,0,0);
    fprintf(coma,'%f\t%f\t%f\n',0,0,0);
    fclose(coma);


    sprintf(' Standard deviation of log_{10} k = %f',sigmaalpha)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(3)
    plot(Xcorr,ones(items,1).*constant,...
        'om',...
        'LineWidth',1,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','b',...
        'MarkerSize',6)
    hold on
    grid on
    plot(Xcorr,ones(items,1).*constant,...
        'om',...
        'LineWidth',1,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor','r',...
        'MarkerSize',5)
    hold on
    grid on
    plot(Xcorr,ones(items,1).*constant,'-m','LineWidth',2)
    legend('T-f_{original}(T)','T-f_{extreme}(T)','T-f_{prior}(T)','location','EastOutside')
    title('T-f_{original}(T), T-f_{extreme}(T) and T-f_{prior}(T) values')
    hold off
    xlabel('T ')
    ylabel('f(T)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(4)
    plot(1000./Xnew',Ynew','--.','LineWidth',0.2)
    xlabel('1000/T')
    ylabel('log_{10}(k)')
    grid on
    hold on
    plot(1000./Xrange(Xsize_1,:),Yrange(Xsize_1,:),'-r','Linewidth',1.5)
    hold on
    plot(1000./Xcorr,recommended'+ones(items,1).*constant,'b','LineWidth',2)
    plot(1000./Xcorr,recommended'-ones(items,1).*constant,'b','LineWidth',2)
    title('Arrhenius plot with main +/- f_{prior}(T) points added')
    axis([1000/range_u*0.9 1000/range_l*1.1 min(recommended'-ones(items,1).*constant)-1 max(recommended'+ones(items,1).*constant)+1])
    hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Assuming normal distribution
else
    
    if (ismember(distribut,'N')==1)

        fourth_answer=[];
        while(isequal(fourth_answer,'An')==0 && isequal(fourth_answer,'AE')==0 && isequal(fourth_answer,'AnE')==0)
            [type_of_pars]=input('Uncertainty of Arrhenius parameters to be determined (AnE,AE,An):','s');
            fourth_answer=type_of_pars;
        end


        if (isequal(type_of_pars,'An')==1 || isequal(type_of_pars,'AE')==1)
            dumny_lines=16;
        end
        if (isequal(type_of_pars,'AnE')==1)
            dumny_lines=17;
        end

%Create input file to program UBAC
        if (isequal(type_of_pars,'AE')==1)

            output=fopen('data.txt','wt+');
            fprintf(output,'* Uncertain Arrhenius parameters (An or AE or AnE)\n');
            fprintf(output,'AE\n');
            fprintf(output,'* Nominal values of all Arrhenius parameters (values for A,n,E/R(K) if present !)\n');
            fprintf(output,'%d\t%d\t%d\n',vertcat(A(Xsize_1)',n(Xsize_1),E_per_R(Xsize_1)'));
            fprintf(output,'* if n is uncertain, set limits for n: n_min, n_max\n%d\t%d \n*Number of data (1st row),data in rows: temperature uncertainty (on log10k scale) \n%d\n',n(Xsize_1)-2,n(Xsize_1)+2,items);
            fprintf(output,'%d\t%.3e\n',horzcat(Xcorr,Ycorr)');
            fclose(output);
        end


        if (isequal(type_of_pars,'An')==1)

            output=fopen('data.txt','wt+');
            fprintf(output,'* Uncertain Arrhenius parameters (An or AE or AnE)\n');
            fprintf(output,'An\n');
            fprintf(output,'* Nominal values of all Arrhenius parameters (values for A,n,E/R(K) if present !)\n');
            fprintf(output,'%d\t%d\t%d\n',vertcat(A(Xsize_1)',n(Xsize_1),E_per_R(Xsize_1)'));
            fprintf(output,'* if n is uncertain, set limits for n: n_min, n_max\n%d\t%d \n*Number of data (1st row),data in rows: temperature uncertainty (on log10k scale) \n%d\n',n(Xsize_1)-2,n(Xsize_1)+2,items);
            fprintf(output,'%d\t%.3e\n',horzcat(Xcorr,Ycorr)');
            fclose(output);

        end



        if (isequal(type_of_pars,'AnE')==1)
            output=fopen('data.txt','wt+');
            fprintf(output,'* Uncertain Arrhenius parameters (An or AE or AnE)\n');
            fprintf(output,'AnE\n');
            fprintf(output,'* Nominal values of all Arrhenius parameters (values for A,n,E/R(K) if present !)\n');
            fprintf(output,'%d\t%d\t%d\n',vertcat(A(Xsize_1)',n(Xsize_1),E_per_R(Xsize_1)'));
            fprintf(output,'* if n is uncertain, set limits for n: n_min, n_max\n%d\t%d \n*Number of data (1st row),data in rows: temperature uncertainty (on log10k scale) \n%d\n',n(Xsize_1)-2,n(Xsize_1)+2,items);
            fprintf(output,'%d\t%.3e\n',horzcat(Xcorr,Ycorr)');
            fclose(output);
        end

        %Calling UBAC
        if(isunix)
            !./ubac.x
        else
            !ubac.exe
        end
              
% Reading f_extreme values from the output file of UBAC

        chose=fopen('data.txt.input_for_jpdap.txt','r');
        for i=1:6
            fgetl(chose);
        end

        [read]=fscanf(chose,'%f',[2,items]);
        [read]=[read]';

        temp=read(:,1);
        f_consist=read(:,2);
% Rename the output file which created by UBAC as an input file to JPDAP
        fclose(chose);
        copyfile('data.txt.input_for_jpdap.txt','data.txt')



%Plot of the f_extreme values

        plot(temp',f_consist,...
            'om',...
            'LineWidth',1,...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor','r',...
            'MarkerSize',6)
        grid on
        xlabel('T')
        ylabel('f(T)')
        title('T-f_{original}(T) and T-f_{extreme}(T) values')
        hold on


%Calling program JPDAP

if(isunix)
    !./jpdap.x
else
    !jpdap.exe
end


        format long
        fid2=fopen('data.txt_fit_minRMSD.txt','r');
        for i=1:dumny_lines
            fgetl(fid2);
        end


        [second_data]=fscanf(fid2,'%f',[5,items]);
        [second_data]=[second_data]';

        fclose(fid2);
       
%Create CovMatrix.txt

        coma2=fopen('CovMatrix.txt','wt+');
        fprintf(coma2,'* Assumed distribution\n');
        fprintf(coma2,'normal\n');
        covFile = fopen('data.txt_fit_minRMSD.txt','r');
        while 1
            line = fgetl(covFile);
            if line == -1
                break
            end
            fprintf(coma2,'%s\n',line);
        end
        fclose(covFile);
        fclose(coma2);

% Create Figure 3, f_original, f_extre, f_prior in one figure

        plot(second_data(:,1),second_data(:,5),'-m','LineWidth',2)
        legend('T-f_{original}(T)','T-f_{extreme}(T)','T-f_{prior}(T)','location','EastOutside')
        title('T-f_{original}(T), T-f_{extreme}(T) and T-f_{prior}(T) values')
        hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(4)
        plot(1000./Xnew',Ynew','--.','LineWidth',0.2)
        hold on
        plot(1000./tt_Xnew',tt_Ynew','-','Linewidth',2)
        xlabel('1000/T')
        ylabel('log_{10}(k)')
        grid on
        hold on
        plot(1000./Xrange(Xsize_1,:),Yrange(Xsize_1,:),'-r','Linewidth',1.5)
        hold on          
        plot(1000./second_data(:,1),recommended'+second_data(:,5),'b','LineWidth',2)
        plot(1000./second_data(:,1),recommended'-second_data(:,5),'b','LineWidth',2)
        title('Arrhenius plot with main +/- f_{prior}(T) points added')
        axis([1000/range_u*0.9 1000/range_l*1.1 min(recommended'-second_data(:,5))-1 max(recommended'+second_data(:,5))+1])
        hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    else
        disp('Not valid distribution ID')
    end

end



