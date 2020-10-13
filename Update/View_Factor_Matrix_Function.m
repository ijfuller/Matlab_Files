function [VF_Matrix] = View_Factor_Matrix_Function(N_el_domain,F_Panel_1_to_Panel,F_Panel_1_Floor,F_Panel_2_Floor,F_Roof_Floor,F_Aper_Floor,F_Aper_Roof,A_Roof,A_Panel,A_Aperature,A_Top_Lip,A_Bottom_Lip,F_Top_Lip_Floor,F_Bottom_Lip_Floor,F_Top_Lip_Roof,F_Bottom_Lip_Roof,Top_Lip,Bottom_Lip,Save_View_Factor_File,Save_View_Factor_Filename,panel_n,panel_m,cover_n,cover_m);
    VF_Matrix=zeros(N_el_domain,N_el_domain);
    A_Panel_El=A_Panel/(panel_n*panel_m);
    A_Floor_El=A_Roof;
    
F_Panel_1_Aper=F_Panel_1_to_Panel(1,:,4);
    F_Panel_2_Aper=F_Panel_1_to_Panel(1,:,5);
    
    F_Panel_3_Aper=fliplr(F_Panel_2_Aper);
   
    F_Panel_4_Aper=fliplr(F_Panel_1_Aper);
    
       for i=1:(panel_m*panel_n)
        for j=1:(cover_n*cover_m)
           if mod(i,2)==0
               F_Panel_3_Floor(i,j)=F_Panel_2_Floor(i-1,j);
               F_Panel_4_Floor(i,j)=F_Panel_1_Floor(i-1,j);
           else
               F_Panel_3_Floor(i,j)=F_Panel_2_Floor(i+1,j);
               F_Panel_4_Floor(i,j)=F_Panel_1_Floor(i+1,j);
           end
        end
   end
    F_Panel_3_Floor=fliplr(F_Panel_3_Floor);
    F_Panel_4_Floor=fliplr(F_Panel_4_Floor);
    
        for i=1:(panel_m*panel_n)
        for j=1:(cover_n*cover_m)
            if i+(panel_n*panel_m)/2 > (panel_n*panel_m)
                F_Panel_1_Roof(i,j)=F_Panel_1_Floor(i-(panel_n*panel_m)/2, j);
                F_Panel_2_Roof(i,j)=F_Panel_2_Floor(i-(panel_n*panel_m)/2, j);
                F_Panel_3_Roof(i,j)=F_Panel_3_Floor(i-(panel_n*panel_m)/2, j);
                F_Panel_4_Roof(i,j)=F_Panel_4_Floor(i-(panel_n*panel_m)/2, j);
            else
                F_Panel_1_Roof(i,j)=F_Panel_1_Floor(i+(panel_n*panel_m)/2, j);
                F_Panel_2_Roof(i,j)=F_Panel_2_Floor(i+(panel_n*panel_m)/2, j);
                F_Panel_3_Roof(i,j)=F_Panel_3_Floor(i+(panel_n*panel_m)/2, j);
                F_Panel_4_Roof(i,j)=F_Panel_4_Floor(i+(panel_n*panel_m)/2, j);
            end
        end
    end

    Panel_to_Panel=zeros(panel_n*panel_m,panel_m*panel_n); %Initialize the panel to panel array
    for i=1:3 %Loop over number of panel to panel interactions
       Panel_to_Panel(:,:,i+1)=[F_Panel_1_to_Panel(:,:,i)]; %Adding all Panel to Panel VF to one array
    end
    
    %F_Panel_1_Roof=F_Panel_1_Roof.*ones(panel_n*panel_m,2*cover_n*cover_m) %Expanding P1-Roof Matrix 

    %%%for panel 1%%%
    for ii=1:panel_m*panel_n %Loop over number of panel elements
        for jj=1:(4*panel_m*panel_n) %Loop over number of panel elements
            VF_Matrix(ii,jj)=Panel_to_Panel(ii,jj); %Add panel to panel VF for elements of panel 1 to other panels
        end
        for jj=1+(4*panel_m*panel_n):(cover_n*cover_m)+(4*panel_m*panel_n)
           VF_Matrix(ii,jj)=F_Panel_1_Floor(ii,(jj-(4*panel_m*panel_n))); %Add view factor from panel 1 to roof
        end
        for jj=1+(4*panel_m*panel_n)+(cover_n*cover_m):(4*panel_m*panel_n)+2*(cover_n*cover_m)
           VF_Matrix(ii,jj)=F_Panel_1_Roof(ii,(jj-((4*panel_m*panel_n)+(cover_n*cover_m)))); %Add view factor from panel 1 to floor
        end
        for jj=1+(4*panel_m*panel_n)+(2*cover_n*cover_m)
            VF_Matrix(ii,jj)=F_Panel_1_Aper(ii);
        end
    end

    %%%for panel 2%%%
    Panel_to_Panel_Aug=[fliplr(flipud(Panel_to_Panel(:,:,2))) Panel_to_Panel(:,:,1) Panel_to_Panel(:,:,2) Panel_to_Panel(:,:,3)];    
    for ii=(1+(panel_n*panel_m)):(2*panel_m*panel_n) %Loop over number of panel elements
        for jj=1:(4*panel_m*panel_n) %Loop over number of panel elements
            VF_Matrix(ii,jj)=Panel_to_Panel_Aug(ii-(panel_m*panel_n),jj); %Add panel to panel VF for elements of panel 2 to other panels
        end
        for jj=1+(4*panel_m*panel_n):(cover_n*cover_m)+(4*panel_m*panel_n)
           VF_Matrix(ii,jj)=F_Panel_2_Floor(ii-(panel_n*panel_m),(jj-(4*panel_m*panel_n))); %Add view factor from panel 1 to roof
        end
        for jj=1+(4*panel_m*panel_n)+(cover_n*cover_m):(4*panel_m*panel_n)+2*(cover_n*cover_m)
           VF_Matrix(ii,jj)=F_Panel_2_Roof(ii-(panel_m*panel_n),(jj-((4*panel_m*panel_n)+(cover_n*cover_m)))); %Add view factor from panel 1 to floor
        end
        for jj=1+(4*panel_m*panel_n)+(2*cover_n*cover_m)
            VF_Matrix(ii,jj)=F_Panel_2_Aper(ii-(panel_m*panel_n));
        end
    end

    %%%for panel 3%%%
    Panel_to_Panel_Aug=[fliplr(flipud(Panel_to_Panel(:,:,3))) fliplr(flipud(Panel_to_Panel(:,:,2))) Panel_to_Panel(:,:,1) Panel_to_Panel(:,:,2)];
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
    Panel_to_Panel_Aug=[fliplr(flipud(Panel_to_Panel(:,:,4))) fliplr(flipud(Panel_to_Panel(:,:,3))) fliplr(flipud(Panel_to_Panel(:,:,2))) Panel_to_Panel(:,:,1)];
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
 
    %%%For floor elements%%%
    for ii=(1+(4*panel_m*panel_n)):((cover_n*cover_m)+(4*panel_m*panel_n))
        VF_Matrix(ii,:)=VF_Matrix(:,ii)*(A_Panel_El/A_Floor_El(ii-(4*panel_m*panel_n)));
        for jj=(1+(4*panel_m*panel_n)+(cover_n*cover_m)):((cover_n*cover_m*2)+(4*panel_m*panel_n))
            VF_Matrix(ii,jj)=F_Roof_Floor(ii-(4*panel_m*panel_n),jj-(4*panel_m*panel_n)-(cover_n*cover_m));
        end
        for jj=1+((cover_n*cover_m*2)+(4*panel_m*panel_n))
            VF_Matrix(ii,jj)=F_Aper_Floor((jj-((cover_n*cover_m*2)+(4*panel_m*panel_n))),(ii-(4*panel_m*panel_n)))*(A_Aperature/A_Roof(ii-(4*panel_m*panel_n)));
        end
    end
    
        %%%For roof elements%%%
    for ii=(1+(4*panel_m*panel_n)+(cover_n*cover_m)):((2*cover_n*cover_m)+(4*panel_m*panel_n))
        VF_Matrix(ii,:)=VF_Matrix(:,ii)*(A_Panel_El/A_Floor_El(ii-(4*panel_m*panel_n)-(cover_n*cover_m)));
        for jj=(1+(4*panel_m*panel_n)):((cover_n*cover_m)+(4*panel_m*panel_n))
            VF_Matrix(ii,jj)=F_Roof_Floor(ii-(4*panel_m*panel_n)-(cover_n*cover_m),jj-(4*panel_m*panel_n));
        end
        for jj=1+((cover_n*cover_m*2)+(4*panel_m*panel_n))
            VF_Matrix(ii,jj)=F_Aper_Floor((jj-((cover_n*cover_m*2)+(4*panel_m*panel_n))),(ii-(4*panel_m*panel_n)-(cover_n*cover_m)))*(A_Aperature/A_Roof(ii-(4*panel_m*panel_n+cover_n*cover_m)));
        end
    end

    %%For Aperture%%%
    ii=1+(4*panel_m*panel_n)+(2*cover_n*cover_m);
    for jj=1:(4*panel_n*panel_m)
        VF_Matrix((N_el_domain),jj)=VF_Matrix(jj,(N_el_domain))*(A_Panel_El/A_Aperature);
    end
    for jj=1+(4*panel_n*panel_m):(4*panel_n*panel_m)+(cover_n*cover_m)
        VF_Matrix((N_el_domain),jj)=VF_Matrix(jj,(N_el_domain))*(A_Floor_El(jj-(4*panel_n*panel_m))/A_Aperature);
    end
    for jj=1+(4*panel_n*panel_m)+(cover_n*cover_m):(4*panel_n*panel_m)+(2*cover_n*cover_m)
        VF_Matrix((N_el_domain),jj)=VF_Matrix(jj,(N_el_domain))*(A_Floor_El(jj-(4*panel_n*panel_m+(cover_n*cover_m)))/A_Aperature);
    end

% %     if Save_View_Factor_File==1
% %         xlswrite(Save_View_Factor_Filename,VF_Matrix)
% %     end

    %Check
    for ii=1:N_el_domain
       check(ii)=sum(VF_Matrix(ii,:));
    end
    check
end