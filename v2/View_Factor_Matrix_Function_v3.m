function [VF_Matrix] = View_Factor_Matrix_Function_v3(N_el_domain,F_Panel_1_to_Panel,F_Panel_2_to_Panel,F_Panel_1_Floor,F_Panel_2_Floor,F_Roof_Floor,F_Aper_Floor,F_Aper_Roof,A_Roof,A_Panel,A_Aperature,A_Top_Lip,A_Bottom_Lip,F_Top_Lip_Floor,F_Bottom_Lip_Floor,F_Top_Lip_Roof,F_Bottom_Lip_Roof,Top_Lip,Bottom_Lip,Save_View_Factor_File,Save_View_Factor_Filename,panel_m)
    VF_Matrix=zeros(N_el_domain,N_el_domain);
    A_Panel=A_Panel/panel_m;
    Panel_to_Panel=zeros(panel_m,panel_m);
    for i=1:3
       Panel_to_Panel=[Panel_to_Panel F_Panel_1_to_Panel(:,:,i)]; %Adding all Panel to Panel VF to one array
    end
    F_P1_to_Floor=[flipud(F_Panel_1_Floor) F_Panel_1_Floor]; %Panel 1 to floor
    F_P2_to_Floor=[F_Panel_2_Floor flipud(F_Panel_2_Floor)]; %Panel 2 to floor
    for i=4:5
       Panel_to_Aperature((i-3),:)=[F_Panel_1_to_Panel(1,:,i)]; %Panel 1 and 2 to aperature
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%for panel 1%%%
    for ii=1:panel_m
        for jj=1:(4*panel_m)
            VF_Matrix(ii,jj)=Panel_to_Panel(ii,jj); %Add panel to panel VF for elements of panel 1
        end
        for jj=(1+4*panel_m):(2+4*panel_m)
           VF_Matrix(ii,jj)=F_P1_to_Floor(ii,jj-panel_m*4); %Add panel 1 to floor/roof
        end
        for jj=(4*panel_m)+3
           kk=1;
           VF_Matrix(ii,jj)=Panel_to_Aperature(kk,ii); %Add panel 1 to aperature
        end
    end


    %%%for panel 2%%%
    for ii=(1+panel_m):(panel_m*2)
        for jj=1:(panel_m)
            VF_Matrix(ii,jj)=Panel_to_Panel(ii-panel_m,jj+panel_m); %Add panel to adjacent for panel 2
        end
        for jj=(1+panel_m):(panel_m*4)
           VF_Matrix(ii,jj)=Panel_to_Panel(ii-panel_m,jj-panel_m); %Add rest of panel to panel for panel 2
        end
        for jj=(1+4*panel_m):(2+4*panel_m)
           VF_Matrix(ii,jj)=F_P2_to_Floor(ii-panel_m,jj-panel_m*4); %Add panel 2 to floor/roof
        end
        for jj=(4*panel_m)+3
           kk=2;
           VF_Matrix(ii,jj)=Panel_to_Aperature(kk,ii-panel_m); %Add panel 2 to aperature
        end
    end

    %%%for panel 3%%%
    for ii=(1+panel_m*2):(panel_m*3)
        for jj=1:(panel_m)
            VF_Matrix(ii,jj)=Panel_to_Panel(ii-panel_m*2,jj+panel_m*2); %2 panel offset
        end
        for jj=(1+panel_m):(panel_m*2)
           VF_Matrix(ii,jj)=Panel_to_Panel(ii-panel_m*2,jj); %1 panel offset
        end
        for jj=(1+panel_m*2):(panel_m*4)
           VF_Matrix(ii,jj)=Panel_to_Panel(ii-panel_m*2,jj-panel_m*2); %1 panel offset & zero panel offset
        end
        for jj=(1+4*panel_m):(2+4*panel_m)
           VF_Matrix(ii,jj)=F_P2_to_Floor(ii-panel_m*2,jj-panel_m*4); %Add panel 3 to floor/roof
        end
        for jj=(4*panel_m)+3
           kk=2;
           VF_Matrix(ii,jj)=Panel_to_Aperature(kk,ii-panel_m*2); %Add panel 2 to aperature
        end
    end

    %%%for panel 4%%%
    for ii=(1+panel_m*3):(panel_m*4)
        for jj=1:panel_m
           VF_Matrix(ii,jj)=Panel_to_Panel(ii-panel_m*3,jj+panel_m*3); %3 panel offset
        end
        for jj=(1+panel_m):(panel_m*2)
            VF_Matrix(ii,jj)=Panel_to_Panel(ii-panel_m*3,jj+panel_m*2-panel_m); %2 panel offset
        end
        for jj=(1+panel_m*2):(panel_m*3)
            VF_Matrix(ii,jj)=Panel_to_Panel(ii-panel_m*3,jj+panel_m-panel_m*2); %1 panel offset
        end
        for jj=(1+panel_m*3):(panel_m*4)
            VF_Matrix(ii,jj)=Panel_to_Panel(ii-panel_m*3,jj-panel_m*3); %0 panel offset
        end
        for jj=(1+4*panel_m):(2+4*panel_m)
           VF_Matrix(ii,jj)=F_P1_to_Floor(ii-panel_m*3,jj-panel_m*4); %Add panel 2 to floor/roof
        end
        for jj=(4*panel_m)+3
           kk=1;
           VF_Matrix(ii,jj)=Panel_to_Aperature(kk,ii-panel_m*3); %Add panel 2 to aperature
        end
    end

    %%%Roof/Floor%%%
    for ii=1+panel_m*4 %Roof
       for jj=1:panel_m*4
           VF_Matrix(ii,jj)=VF_Matrix(jj,ii)*(A_Panel/A_Roof); %Adding roof to panel 1-4
       end
       for jj=(2+panel_m*4)
          VF_Matrix(ii,jj)=F_Roof_Floor; %Roof to Floor
       end
       for jj=3+panel_m*4
          VF_Matrix(ii,jj)=F_Aper_Floor(1)*(A_Aperature/A_Roof); %Floor to aperature
       end
    end
    for ii=2+panel_m*4 %Floor
       for jj=1:panel_m*4
           VF_Matrix(ii,jj)=VF_Matrix(jj,ii)*(A_Panel/A_Roof); %Adding floor to panel 1-4
       end
       for jj=(1+panel_m*4)
          VF_Matrix(ii,jj)=F_Roof_Floor; %Floor to Roof
       end
       for jj=3+panel_m*4
          VF_Matrix(ii,jj)=F_Aper_Floor(1)*(A_Aperature/A_Roof); %Floor to aperature
       end
    end

    %%%Aperature%%%
    for ii=3+panel_m*4
       for jj=1:panel_m*4
          VF_Matrix(ii,jj)=VF_Matrix(jj,ii)*(A_Panel/A_Aperature); %Aperature to panel
       end
       for jj=1+panel_m*4:2+panel_m*4
          VF_Matrix(ii,jj)=VF_Matrix(jj,ii)*(A_Roof/A_Aperature); %Aperature to roof/floor
       end
    end
    
    %%%Code for top and bottom lip
    if Top_Lip > 0
        VF_Matrix(N_el_domain,:)=VF_Matrix(3+panel_m*4,:);
        VF_Matrix(3+panel_m*4,:)=0;
    end
    if Bottom_Lip > 0
        VF_Matrix(:,N_el_domain)=VF_Matrix(:,3+panel_m*4);
        VF_Matrix(:,3+panel_m*4)=0;
    end

    %%%If Top Lip is active
    if Top_Lip > 0
        F_P1_Top_Lip(:,:)=F_Panel_1_to_Panel(1,:,6);
        F_P2_Top_Lip(:,:)=F_Panel_1_to_Panel(1,:,7);
        for ii=[1,4]
            for jj=1:panel_m
                VF_Matrix(jj+(ii-1)*panel_m,panel_m*4+3)=F_P1_Top_Lip(jj); %Adds outer panel to top lip
            end
        end
        for ii=[2,3]
            for jj=1:panel_m
                VF_Matrix(jj+(ii-1)*panel_m,panel_m*4+3)=F_P2_Top_Lip(jj); %Adds inner panel to top lip
            end
        end
        VF_Matrix(panel_m*4+1,panel_m*4+3)=F_Bottom_Lip_Floor(1)*(A_Top_Lip/A_Roof); %Roof to top lip
        VF_Matrix(panel_m*4+2,panel_m*4+3)=F_Top_Lip_Floor(1)*(A_Top_Lip/A_Roof); %Floor to top lip

        for ii=1:panel_m*4
           VF_Matrix(panel_m*4+3,ii)=VF_Matrix(ii,panel_m*4+3)*(A_Panel/A_Top_Lip); %Top Lip to panels
        end
        for ii=1+panel_m*4:2+panel_m*4
           VF_Matrix(panel_m*4+3,ii)=VF_Matrix(ii,panel_m*4+3)*(A_Roof/A_Top_Lip);%Top Lip to roof/floor
        end
    end
    %%%If bottom lip is active
    if Bottom_Lip > 0
        F_P1_Bottom_Lip(:,:)=F_Panel_1_to_Panel(1,:,8);
        F_P2_Bottom_Lip(:,:)=F_Panel_1_to_Panel(1,:,9);
        for ii=[1,4]
            for jj=1:panel_m
                VF_Matrix(jj+(ii-1)*panel_m,panel_m*4+4)=F_P1_Bottom_Lip(jj); %Adds outer panel to bottom lip
            end
        end
        for ii=[2,3]
            for jj=1:panel_m
                VF_Matrix(jj+(ii-1)*panel_m,panel_m*4+4)=F_P2_Bottom_Lip(jj); %Adds inner panel to bottom lip
            end
        end
        VF_Matrix(panel_m*4+1,panel_m*4+4)=F_Top_Lip_Floor(1)*(A_Bottom_Lip/A_Roof); %Roof to bottom lip
        VF_Matrix(panel_m*4+2,panel_m*4+4)=F_Bottom_Lip_Floor(1)*(A_Bottom_Lip/A_Roof); %Floor to bottom lip

        for ii=1:panel_m*4
           VF_Matrix(panel_m*4+4,ii)=VF_Matrix(ii,panel_m*4+4)*(A_Panel/A_Bottom_Lip); %Top Lip to panels
        end
        for ii=1+panel_m*4:2+panel_m*4
           VF_Matrix(panel_m*4+4,ii)=VF_Matrix(ii,panel_m*4+4)*(A_Roof/A_Bottom_Lip);%Top Lip to roof/floor
        end
    end
    if Save_View_Factor_File==1
        xlswrite(Save_View_Factor_Filename,VF_Matrix)
    end
    VF_Matrix
end