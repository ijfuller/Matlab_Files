function [VF_Matrix] = View_Factor_Matrix_Temp_Function_v3(N_el_domain,F_Panel_1_to_Panel,F_Panel_1_Floor,F_Panel_2_Floor,F_Roof_Floor,F_Aper_Floor,F_Aper_Roof,A_Roof,A_Panel,A_Aperature,A_Top_Lip,A_Bottom_Lip,F_Top_Lip_Floor,F_Bottom_Lip_Floor,F_Top_Lip_Roof,F_Bottom_Lip_Roof,Top_Lip,Bottom_Lip,Save_View_Factor_File,Save_View_Factor_Filename,panel_n,panel_m,cover_n,cover_m);
    
N_el_domain=(4*panel_n*panel_m)+(2*cover_n*cover_m)+1;
    if Top_Lip > 0
        N_el_domain=N_el_domain+1;
    else
        F_Aper_Roof=F_Aper_Floor;
    end
    if Bottom_Lip > 0
        N_el_domain=N_el_domain+1;
    end

    VF_Matrix=zeros(N_el_domain,N_el_domain);
    A_Panel_El=A_Panel/(panel_n*panel_m);
    A_Floor_El=A_Roof;

    
    F_Panel_1_Aper=F_Panel_1_to_Panel(1,:,4);
    F_Panel_2_Aper=F_Panel_1_to_Panel(1,:,5);

    F_Panel_3_Aper=reshape(flipud(reshape(F_Panel_2_Aper,[panel_n,panel_m])),[1,panel_n*panel_m]);
    F_Panel_4_Aper=reshape(flipud(reshape(F_Panel_1_Aper,[panel_n,panel_m])),[1,panel_n*panel_m]);


    F_Panel_3_Floor=reshape(reshape(fliplr(reshape(reshape(F_Panel_2_Floor',[panel_n*panel_m*cover_n*cover_m,1]),[panel_n*cover_n*cover_m,panel_m])')',[panel_n*panel_m*cover_n*cover_m,1]),[cover_n*cover_m,panel_n*panel_m])';
    F_Panel_4_Floor=reshape(reshape(fliplr(reshape(reshape(F_Panel_1_Floor',[panel_n*panel_m*cover_n*cover_m,1]),[panel_n*cover_n*cover_m,panel_m])')',[panel_n*panel_m*cover_n*cover_m,1]),[cover_n*cover_m,panel_n*panel_m])';
    F_Panel_1_Roof=reshape(reshape(flipud(reshape(reshape(F_Panel_1_Floor',[panel_n*panel_m*cover_n*cover_m,1]),[panel_n*cover_n*cover_m,panel_m])')',[panel_n*panel_m*cover_n*cover_m,1]),[cover_n*cover_m,panel_n*panel_m])';;
    F_Panel_2_Roof=reshape(reshape(flipud(reshape(reshape(F_Panel_2_Floor',[panel_n*panel_m*cover_n*cover_m,1]),[panel_n*cover_n*cover_m,panel_m])')',[panel_n*panel_m*cover_n*cover_m,1]),[cover_n*cover_m,panel_n*panel_m])';;
    F_Panel_3_Roof=reshape(reshape(flipud(reshape(reshape(F_Panel_3_Floor',[panel_n*panel_m*cover_n*cover_m,1]),[panel_n*cover_n*cover_m,panel_m])')',[panel_n*panel_m*cover_n*cover_m,1]),[cover_n*cover_m,panel_n*panel_m])';;
    F_Panel_4_Roof=reshape(reshape(flipud(reshape(reshape(F_Panel_4_Floor',[panel_n*panel_m*cover_n*cover_m,1]),[panel_n*cover_n*cover_m,panel_m])')',[panel_n*panel_m*cover_n*cover_m,1]),[cover_n*cover_m,panel_n*panel_m])';;


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
           VF_Matrix(ii,jj)=F_Panel_1_Roof(ii,(jj-((4*panel_m*panel_n)+(cover_n*cover_m)))); %Add view factor from panel 1 to floor
        end
        for jj=N_el_domain
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
        for jj=N_el_domain
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
        for jj=N_el_domain
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
        for jj=N_el_domain
            VF_Matrix(ii,jj)=F_Panel_4_Aper(ii-(panel_m*panel_n*3));
        end
    end

    %%%For floor elements%%%
    for ii=(1+(4*panel_m*panel_n)):((cover_n*cover_m)+(4*panel_m*panel_n))
        VF_Matrix(ii,:)=VF_Matrix(:,ii)*(A_Panel_El/A_Floor_El(ii-(4*panel_m*panel_n)));
        for jj=(1+(4*panel_m*panel_n)+(cover_n*cover_m)):((cover_n*cover_m*2)+(4*panel_m*panel_n))
            VF_Matrix(ii,jj)=F_Roof_Floor(ii-(4*panel_m*panel_n),jj-(4*panel_m*panel_n)-(cover_n*cover_m));
        end
        for jj=1+(4*panel_n*panel_m)+(2*cover_n*cover_m)
            VF_Matrix(ii,N_el_domain)=F_Aper_Floor((jj-((cover_n*cover_m*2)+(4*panel_m*panel_n))),(ii-(4*panel_m*panel_n)))*(A_Aperature/A_Roof(ii-(4*panel_m*panel_n)));
        end
    end

        %%%For roof elements%%%
    for ii=(1+(4*panel_m*panel_n)+(cover_n*cover_m)):((2*cover_n*cover_m)+(4*panel_m*panel_n))
        VF_Matrix(ii,:)=VF_Matrix(:,ii)*(A_Panel_El/A_Floor_El(ii-(4*panel_m*panel_n)-(cover_n*cover_m)));
        for jj=(1+(4*panel_m*panel_n)):((cover_n*cover_m)+(4*panel_m*panel_n))
            VF_Matrix(ii,jj)=F_Roof_Floor(ii-(4*panel_m*panel_n)-(cover_n*cover_m),jj-(4*panel_m*panel_n));
        end
        for jj=1+(4*panel_n*panel_m)+(2*cover_n*cover_m)
            VF_Matrix(ii,N_el_domain)=F_Aper_Roof((jj-((cover_n*cover_m*2)+(4*panel_m*panel_n))),(ii-(4*panel_m*panel_n)-(cover_n*cover_m)))*(A_Aperature/A_Roof(ii-(4*panel_m*panel_n+cover_n*cover_m)));
        end
    end

    %%%For Aperture%%%
    for jj=1:(4*panel_n*panel_m)
        VF_Matrix((N_el_domain),jj)=VF_Matrix(jj,(N_el_domain))*(A_Panel_El/A_Aperature);
    end
    for jj=1+(4*panel_n*panel_m):(4*panel_n*panel_m)+(cover_n*cover_m)
        VF_Matrix((N_el_domain),jj)=VF_Matrix(jj,(N_el_domain))*(A_Floor_El(jj-(4*panel_n*panel_m))/A_Aperature);
    end
    for jj=1+(4*panel_n*panel_m)+(cover_n*cover_m):(4*panel_n*panel_m)+(2*cover_n*cover_m)
        VF_Matrix((N_el_domain),jj)=VF_Matrix(jj,(N_el_domain))*(A_Floor_El(jj-(4*panel_n*panel_m+(cover_n*cover_m)))/A_Aperature);
    end

    %%%Top Lip%%%
    if Top_Lip > 0
        F_Panel_1_Top_Lip=F_Panel_1_to_Panel(1,:,6);
        F_Panel_2_Top_Lip=F_Panel_1_to_Panel(1,:,7);
        F_Panel_3_Top_Lip=fliplr(reshape(fliplr(reshape(F_Panel_2_Top_Lip,[panel_n,panel_m])),[1,panel_n*panel_m]));
        F_Panel_4_Top_Lip=fliplr(reshape(fliplr(reshape(F_Panel_1_Top_Lip,[panel_n,panel_m])),[1,panel_n*panel_m]));
        for jj=1:panel_m*panel_n %Loop over 1st panel elements
            VF_Matrix(jj,(4*panel_n*panel_m+2*cover_n*cover_m+1))=F_Panel_1_Top_Lip(jj);
        end
        for jj=1+panel_n*panel_m:2*panel_m*panel_n %Loop over 2nd panel elements
            VF_Matrix(jj,(4*panel_n*panel_m+2*cover_n*cover_m+1))=F_Panel_2_Top_Lip(jj-panel_n*panel_m);
        end
        for jj=1+panel_n*panel_m*2:3*panel_m*panel_n %Loop over 3rd panel elements
            VF_Matrix(jj,(4*panel_n*panel_m+2*cover_n*cover_m+1))=F_Panel_3_Top_Lip(jj-panel_n*panel_m*2);
        end
        for jj=1+panel_n*panel_m*3:4*panel_m*panel_n %Loop over 4th panel elements
            VF_Matrix(jj,(4*panel_n*panel_m+2*cover_n*cover_m+1))=F_Panel_4_Top_Lip(jj-panel_n*panel_m*3);
        end
        for jj=1+4*panel_m*panel_n:4*panel_n*panel_m+cover_n*cover_m %Add floor-top lip
            VF_Matrix(jj,(1+4*panel_n*panel_m+2*cover_n*cover_m))=(F_Top_Lip_Floor(1,jj-4*panel_m*panel_n))*(A_Top_Lip/A_Floor_El(jj-4*panel_m*panel_n));
        end
        for jj=1+4*panel_m*panel_n+cover_n*cover_m:4*panel_n*panel_m+2*cover_n*cover_m %Add roof-top lip
            VF_Matrix(jj,(1+4*panel_n*panel_m+2*cover_n*cover_m))=(F_Top_Lip_Roof(1,jj-(4*panel_m*panel_n+cover_n*cover_m)))*(A_Top_Lip/A_Floor_El(jj-(4*panel_m*panel_n+cover_n*cover_m)));
        end

        for jj=1:4*panel_n*panel_m %Top Lip-Panel Elements
           VF_Matrix((1+4*panel_n*panel_m+2*cover_n*cover_m),jj)=VF_Matrix(jj,(1+4*panel_n*panel_m+2*cover_n*cover_m))*(A_Panel_El/A_Top_Lip); 
        end
        for jj=1+4*panel_n*panel_m:4*panel_n*panel_m+cover_n*cover_m %Top Lip-Floor
           VF_Matrix((1+4*panel_n*panel_m+2*cover_n*cover_m),jj)=VF_Matrix(jj,(1+4*panel_n*panel_m+2*cover_n*cover_m))*(A_Floor_El(jj-4*panel_n*panel_m)/A_Top_Lip); 
        end
        for jj=1+4*panel_n*panel_m+cover_n*cover_m:4*panel_n*panel_m+2*cover_n*cover_m %Top Lip-Roof
           VF_Matrix((1+4*panel_n*panel_m+2*cover_n*cover_m),jj)=VF_Matrix(jj,(1+4*panel_n*panel_m+2*cover_n*cover_m))*(A_Floor_El(jj-(4*panel_n*panel_m+cover_n*cover_m))/A_Top_Lip); 
        end
    end

    if Bottom_Lip > 0
        F_Panel_1_Bottom_Lip=F_Panel_1_to_Panel(1,:,8);
        F_Panel_2_Bottom_Lip=F_Panel_1_to_Panel(1,:,9);
        F_Panel_3_Bottom_Lip=fliplr(reshape(fliplr(reshape(F_Panel_2_Bottom_Lip,[panel_n,panel_m])),[1,panel_n*panel_m]));
        F_Panel_4_Bottom_Lip=fliplr(reshape(fliplr(reshape(F_Panel_1_Bottom_Lip,[panel_n,panel_m])),[1,panel_n*panel_m]));
        for jj=1:panel_m*panel_n %Loop over 1st panel elements
            VF_Matrix(jj,(4*panel_n*panel_m+2*cover_n*cover_m+2))=F_Panel_1_Bottom_Lip(jj);
        end
        for jj=1+panel_n*panel_m:2*panel_m*panel_n %Loop over 2nd panel elements
            VF_Matrix(jj,(4*panel_n*panel_m+2*cover_n*cover_m+2))=F_Panel_2_Bottom_Lip(jj-panel_n*panel_m);
        end
        for jj=1+panel_n*panel_m*2:3*panel_m*panel_n %Loop over 3rd panel elements
            VF_Matrix(jj,(4*panel_n*panel_m+2*cover_n*cover_m+2))=F_Panel_3_Bottom_Lip(jj-panel_n*panel_m*2);
        end
        for jj=1+panel_n*panel_m*3:4*panel_m*panel_n %Loop over 4th panel elements
            VF_Matrix(jj,(4*panel_n*panel_m+2*cover_n*cover_m+2))=F_Panel_4_Bottom_Lip(jj-panel_n*panel_m*3);
        end
        for jj=1+4*panel_m*panel_n:4*panel_n*panel_m+cover_n*cover_m %Add floor-bottom lip
            VF_Matrix(jj,(2+4*panel_n*panel_m+2*cover_n*cover_m))=(F_Bottom_Lip_Floor(1,jj-4*panel_m*panel_n))*(A_Bottom_Lip/A_Floor_El(jj-4*panel_m*panel_n));
        end
        for jj=1+4*panel_m*panel_n+cover_n*cover_m:4*panel_n*panel_m+2*cover_n*cover_m %Add roof-bottom lip
            VF_Matrix(jj,(2+4*panel_n*panel_m+2*cover_n*cover_m))=(F_Bottom_Lip_Roof(1,jj-(4*panel_m*panel_n+cover_n*cover_m)))*(A_Bottom_Lip/A_Floor_El(jj-(4*panel_m*panel_n+cover_n*cover_m)));
        end

        for jj=1:4*panel_n*panel_m %Bottom Lip-Panel Elements
           VF_Matrix((2+4*panel_n*panel_m+2*cover_n*cover_m),jj)=VF_Matrix(jj,(2+4*panel_n*panel_m+2*cover_n*cover_m))*(A_Panel_El/A_Bottom_Lip); 
        end
        for jj=1+4*panel_n*panel_m:4*panel_n*panel_m+cover_n*cover_m %Bottom Lip-Floor
           VF_Matrix((2+4*panel_n*panel_m+2*cover_n*cover_m),jj)=VF_Matrix(jj,(2+4*panel_n*panel_m+2*cover_n*cover_m))*(A_Floor_El(jj-4*panel_n*panel_m)/A_Bottom_Lip); 
        end
        for jj=1+4*panel_n*panel_m+cover_n*cover_m:4*panel_n*panel_m+2*cover_n*cover_m %Bottom Lip-Floor
           VF_Matrix((2+4*panel_n*panel_m+2*cover_n*cover_m),jj)=VF_Matrix(jj,(2+4*panel_n*panel_m+2*cover_n*cover_m))*(A_Floor_El(jj-(4*panel_n*panel_m+cover_n*cover_m))/A_Bottom_Lip); 
        end 
    end


        if Save_View_Factor_File==1
            xlswrite(Save_View_Factor_Filename,VF_Matrix)
        end

%         VF_Matrix
% 
%         %Check
%         for ii=1:N_el_domain%1:N_el_domain
%            check(ii)=sum(VF_Matrix(ii,:));
%         end
%         check
display('-----------------------------------------------')
display('View Factor Matrix Completed')
display('-----------------------------------------------')

        
end