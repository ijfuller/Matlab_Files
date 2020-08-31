function [VF_Matrix] = View_Factor_Matrix_Function(N_el_domain,F_Panel_1_to_Panel,F_Panel_1_Floor,F_Panel_2_Floor,F_Roof_Floor,F_Aper_Floor,F_Aper_Roof,A_Roof,A_Panel,A_Aperature,A_Top_Lip,A_Bottom_Lip,F_Top_Lip_Floor,F_Bottom_Lip_Floor,F_Top_Lip_Roof,F_Bottom_Lip_Roof,Top_Lip,Bottom_Lip,Save_View_Factor_File,Save_View_Factor_Filename,panel_n,panel_m,cover_n,cover_m);
    VF_Matrix=zeros(N_el_domain,N_el_domain);
    A_Panel_El=A_Panel/(panel_n*panel_m);
    A_Floor_El=A_Roof/(cover_n*cover_m);
    
    [F_Panel_1_Roof,F_Panel_2_Roof] = VF_Reorder_Panel_Roof(panel_n,panel_m,cover_n,cover_m,F_Panel_1_Floor,F_Panel_2_Floor);
    F_Panel_1_Aper=F_Panel_1_to_Panel(1,:,4);
    F_Panel_2_Aper=F_Panel_1_to_Panel(1,:,5);
    F_Panel_3_Floor=fliplr(flipud((F_Panel_2_Roof)));
    F_Panel_3_Roof=fliplr(flipud((F_Panel_2_Floor)));
    F_Panel_3_Aper=fliplr(F_Panel_2_Aper);
    F_Panel_4_Floor=fliplr(flipud((F_Panel_1_Roof)));
    F_Panel_4_Roof=fliplr(flipud((F_Panel_1_Floor)));
    F_Panel_4_Aper=fliplr(F_Panel_1_Aper);

    Panel_to_Panel=zeros(panel_n*panel_m,panel_m*panel_n); %Initialize the panel to panel array
    for i=1:3 %Loop over number of panel to panel interactions
       Panel_to_Panel(:,:,i+1)=[F_Panel_1_to_Panel(:,:,i)]; %Adding all Panel to Panel VF to one array
    end

    %%%for panel 1%%%
    for ii=1:panel_m*panel_n %Loop over number of panel elements
        for jj=1:(4*panel_m*panel_n) %Loop over number of panel elements
            VF_Matrix(ii,jj)=Panel_to_Panel(ii,jj); %Add panel to panel VF for elements of panel 1 to other panels
        end
        for jj=1+(4*panel_m*panel_n):(cover_n*cover_m)+(4*panel_m*panel_n)
           VF_Matrix(ii,jj)=F_Panel_1_Floor(ii,(jj-(4*panel_m*panel_n))); %Add view factor from panel 1 to roof
        end
        for jj=1+(4*panel_m*panel_n)+(cover_n*cover_m):(4*panel_m*panel_n)+2*(cover_n*cover_m)
           VF_Matrix(ii,jj)=F_Panel_1_Roof(ii,(jj-(4*panel_m*panel_n)-(cover_n*cover_m))); %Add view factor from panel 1 to floor
        end
        for jj=1+(4*panel_m*panel_n)+(2*cover_n*cover_m)
            VF_Matrix(ii,jj)=F_Panel_1_Aper(ii);
        end
    end

    %%%for panel 2%%%
    Panel_to_Panel_Aug=[Panel_to_Panel(:,:,2) Panel_to_Panel(:,:,1) Panel_to_Panel(:,:,2) Panel_to_Panel(:,:,3)];
    for ii=(1+(panel_n*panel_m)):(2*panel_m*panel_n) %Loop over number of panel elements
        for jj=1:(4*panel_m*panel_n) %Loop over number of panel elements
            VF_Matrix(ii,jj)=Panel_to_Panel_Aug(ii-(panel_m*panel_n),jj); %Add panel to panel VF for elements of panel 2 to other panels
        end
        for jj=1+(4*panel_m*panel_n):(cover_n*cover_m)+(4*panel_m*panel_n)
           VF_Matrix(ii,jj)=F_Panel_2_Floor(ii-(panel_m*panel_n),(jj-(4*panel_m*panel_n))); %Add view factor from panel 1 to roof
        end
        for jj=1+(4*panel_m*panel_n)+(cover_n*cover_m):(4*panel_m*panel_n)+2*(cover_n*cover_m)
           VF_Matrix(ii,jj)=F_Panel_2_Roof(ii-(panel_m*panel_n),(jj-(4*panel_m*panel_n)-(cover_n*cover_m))); %Add view factor from panel 1 to floor
        end
        for jj=1+(4*panel_m*panel_n)+(2*cover_n*cover_m)
            VF_Matrix(ii,jj)=F_Panel_2_Aper(ii-(panel_m*panel_n));
        end
    end

    %%%for panel 3%%%
    Panel_to_Panel_Aug=flipud([Panel_to_Panel(:,:,3) Panel_to_Panel(:,:,2) Panel_to_Panel(:,:,1) Panel_to_Panel(:,:,2)]);
    for ii=(1+(panel_n*panel_m*2)):(3*panel_m*panel_n) %Loop over number of panel elements
        for jj=1:(4*panel_m*panel_n) %Loop over number of panel elements
            VF_Matrix(ii,jj)=Panel_to_Panel_Aug(ii-(panel_m*panel_n*2),jj); %Add panel to panel VF for elements of panel 2 to other panels
        end
        for jj=1+(4*panel_m*panel_n):(cover_n*cover_m)+(4*panel_m*panel_n)
           VF_Matrix(ii,jj)=F_Panel_3_Floor(ii-(panel_m*panel_n*2),(jj-(4*panel_m*panel_n))); %Add view factor from panel 1 to roof
        end
        for jj=1+(4*panel_m*panel_n)+(cover_n*cover_m):(4*panel_m*panel_n)+2*(cover_n*cover_m)
           VF_Matrix(ii,jj)=F_Panel_3_Roof(ii-(panel_m*panel_n*2),(jj-(4*panel_m*panel_n)-(cover_n*cover_m))); %Add view factor from panel 1 to floor
        end
        for jj=1+(4*panel_m*panel_n)+(2*cover_n*cover_m)
            VF_Matrix(ii,jj)=F_Panel_3_Aper(ii-(panel_m*panel_n*2));
        end
    end

    %%%for panel 4%%%
    Panel_to_Panel_Aug=flipud([Panel_to_Panel(:,:,4) Panel_to_Panel(:,:,3) Panel_to_Panel(:,:,2) Panel_to_Panel(:,:,1)]);
    for ii=(1+(panel_n*panel_m*3)):(4*panel_m*panel_n) %Loop over number of panel elements
        for jj=1:(4*panel_m*panel_n) %Loop over number of panel elements
            VF_Matrix(ii,jj)=Panel_to_Panel_Aug(ii-(panel_m*panel_n*3),jj); %Add panel to panel VF for elements of panel 2 to other panels
        end
        for jj=1+(4*panel_m*panel_n):(cover_n*cover_m)+(4*panel_m*panel_n)
           VF_Matrix(ii,jj)=F_Panel_4_Floor(ii-(panel_m*panel_n*3),(jj-(4*panel_m*panel_n))); %Add view factor from panel 1 to roof
        end
        for jj=1+(4*panel_m*panel_n)+(cover_n*cover_m):(4*panel_m*panel_n)+2*(cover_n*cover_m)
           VF_Matrix(ii,jj)=F_Panel_4_Roof(ii-(panel_m*panel_n*3),(jj-(4*panel_m*panel_n)-(cover_n*cover_m))); %Add view factor from panel 1 to floor
        end
        for jj=1+(4*panel_m*panel_n)+(2*cover_n*cover_m)
            VF_Matrix(ii,jj)=F_Panel_4_Aper(ii-(panel_m*panel_n*3));
        end
    end

    %%%For floor and roof elements%%%
    for ii=(1+(4*panel_m*panel_n)):((cover_n*cover_m*2)+(4*panel_m*panel_n))
        VF_Matrix(ii,:)=VF_Matrix(:,ii)*(A_Panel_El/A_Floor_El);
    end
    for ii=(1+(4*panel_m*panel_n)):((cover_n*cover_m)+(4*panel_m*panel_n))
        for jj=(1+(4*panel_m*panel_n)+(cover_n*cover_m)):((cover_n*cover_m*2)+(4*panel_m*panel_n))
            VF_Matrix(ii,jj)=F_Roof_Floor(ii-(4*panel_m*panel_n),jj-(4*panel_m*panel_n)-(cover_n*cover_m));
            VF_Matrix(jj,ii)=F_Roof_Floor(ii-(4*panel_m*panel_n),jj-(4*panel_m*panel_n)-(cover_n*cover_m));
        end
        for jj=1+(4*panel_n*panel_m)+(2*cover_n*cover_m)
           VF_Matrix(ii,jj)=(F_Aper_Floor(1,ii-(4*panel_m*panel_n)))*(A_Aperature/A_Floor_El); 
        end
        for jj=1+(4*panel_n*panel_m)+(2*cover_n*cover_m)
           VF_Matrix(ii+(cover_n*cover_m),jj)=(F_Aper_Floor(1,ii-(4*panel_m*panel_n)))*(A_Aperature/A_Floor_El); 
        end
    end

    %%%For Aperture%%%
    ii=1+(4*panel_m*panel_n)+(2*cover_n*cover_m);
    for jj=1:(4*panel_n*panel_m)
        VF_Matrix((N_el_domain),jj)=VF_Matrix(jj,(N_el_domain))*(A_Panel_El/A_Aperature);
    end
    for jj=1+(4*panel_n*panel_m):(4*panel_n*panel_m)+(2*cover_n*cover_m)
        VF_Matrix((N_el_domain),jj)=VF_Matrix(jj,(N_el_domain))*(A_Floor_El/A_Aperature);
    end

    if Save_View_Factor_File==1
        xlswrite(Save_View_Factor_Filename,VF_Matrix)
    end

    %%Check
%     for ii=1:N_el_domain
%        check(ii)=sum(VF_Matrix(ii,:));
%     end
%     check
end