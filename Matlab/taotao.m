function[Xnew,Ynew,numbering_new, paper_type_new, paper_ID_new] = taotao(txtname,increment,range_l,range_u)

    [paper_type,paper_ID,considered,Tlower,Tupper,A,n,E_per_R,k,order,f_type,...
        y1,y2,y3,Tmid,multiplier,A_bi,n_bi,E_per_R_bi,T0,interval]=...
        textread(txtname,' %s %s  %f  %f %f %f  %f %f %f%f%f%f%f%f%f%f%f%f%f%f%f'...
        ,'emptyvalue',0);

    [N,M]=size(Tlower);
    
    JustSaveTaoTao = 0;
    for i=1:N
        if(strcmp(paper_type(i),'Compare') == 1)
            JustSaveTaoTao = 1;
        end
    end

    a = 3;

    if(JustSaveTaoTao == 1)
        paper_ID = paper_ID(N-a:N);
        paper_type = paper_type(N-a:N);
        considered = considered(N-a:N);
        Tlower = Tlower(N-a:N);
        Tupper = Tupper(N-a:N);
        A = A(N-a:N);
        n = n(N-a:N);
        E_per_R = E_per_R(N-a:N);
        k = k(N-a:N);
        order = order(N-a:N);
        f_type = f_type(N-a:N);
        y1 = y1(N-a:N);
        y2 = y2(N-a:N);
        y3 = y3(N-a:N);
        Tmid = Tmid(N-a:N);
        multiplier = multiplier(N-a:N);
        A_bi = A_bi(N-a:N);
        n_bi = n_bi(N-a:N);
        E_per_R_bi = E_per_R_bi(N-a:N);
        T0 = T0(N-a:N);
    end
    
    [N,M]=size(Tlower);

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

    % figure(1)
    % plot(1000./Xnew',Ynew','--.','LineWidth',0.2)
    % xlabel('1000/T')
    % ylabel('log_{10}(k)')
    % grid on
    % hold on


    % % plot(1000./Xnew(Xsize_1,:),Ynew(Xsize_1,:),'--','Linewidth',1.5)
    % % title('Arrhenius plot')
    % % hold on
    % legend(strcat(numbering_new,'.',paper_type_new,';',paper_ID_new),'location','EastOutside')
    % hold on
end