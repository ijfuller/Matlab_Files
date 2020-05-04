function [VF_Matrix] = View_Factor_Matrix_Function_v2(N_el_domain,F_Panel_1_to_Panel,F_Panel_2_to_Panel,F_Panel_1_Floor,F_Panel_2_Floor,F_Roof_Floor,F_Aper_Floor,F_Aper_Roof,A_Roof,A_Panel,A_Aperature,A_Top_Lip,A_Bottom_Lip,F_Top_Lip_Floor,F_Bottom_Lip_Floor,F_Top_Lip_Roof,F_Bottom_Lip_Roof,Top_Lip,Bottom_Lip,Save_View_Factor_File,Save_View_Factor_Filename)
    VF_Matrix=zeros(N_el_domain,N_el_domain);%Row 1 is panel 1. Row 2 = panel 2. Roof=5. Floor=6. Aperature=7. Top Lip 8, Bottom Lip 9   
   
    Panel_to_Panel=[0 F_Panel_1_to_Panel(:,:,1) F_Panel_1_to_Panel(:,:,2) F_Panel_1_to_Panel(:,:,3)];
    Panel_to_Roof=[F_Panel_1_Floor F_Panel_2_Floor];
    Panel_to_Aperature=[F_Panel_1_to_Panel(:,:,4) F_Panel_1_to_Panel(:,:,5)];

    %%%For panels to panels 1-4%%%
    for ii=1:4
        for jj=1:4
            VF_Matrix(ii,jj)=Panel_to_Panel(abs(jj-ii)+1);
        end
    end
    %%%For panel to roof/floor%%%
    for ii=1:4
        for jj=5:6
            if ii==1 || ii==4
                VF_Matrix(ii,jj)=Panel_to_Roof(1);
            elseif ii==2 || ii==3
                VF_Matrix(ii,jj)=Panel_to_Roof(2);
            end
        end
    end
    %%%For panel to aperature 
    for ii=1:4
        for jj=7
            if ii==1 || ii==4
                VF_Matrix(ii,jj)=Panel_to_Aperature(1);
            elseif ii==2 || ii==3
                VF_Matrix(ii,jj)=Panel_to_Aperature(2);
            end
        end
    end
    %%%For roof/floor to panels and floor/roof
    for ii=5:6
        for jj=1:6
            if jj==5 || jj==6
                if jj==ii
                    VF_Matrix(ii,jj)=0;
                else
                    VF_Matrix(ii,jj)=F_Roof_Floor;
                end
            else
                VF_Matrix(ii,jj)=VF_Matrix(jj,ii)*(A_Panel/A_Roof);
            end
        end
    end
      %%%Roof to aperature 
      if Top_Lip > 0
            VF_Matrix(5,7)=F_Aper_Roof*(A_Aperature/A_Roof);  
      else
            VF_Matrix(5,7)=F_Aper_Floor*(A_Aperature/A_Roof);
      end
      %%%Floor to aperature
      VF_Matrix(6,7)=F_Aper_Floor*(A_Aperature/A_Roof);
      %%%Aperature to panels
      for ii=1:4
          VF_Matrix(7,ii)=VF_Matrix(ii,7)*(A_Panel/A_Aperature);
      end
      %%%Aperature to roof/floor
      for ii=5:6
          VF_Matrix(7,ii)=VF_Matrix(ii,7)*(A_Roof/A_Aperature);
      end
      
%%%%%%%For top and bottom lip%%%%%%%%%
      if Top_Lip > 0
      Panel_to_Top_Lip=[F_Panel_1_to_Panel(:,:,6), F_Panel_1_to_Panel(:,:,7)];
          %%%for panel 1 & 4 - top lip
          for ii=[1,4]
              for jj=8
                  VF_Matrix(ii,jj)=Panel_to_Top_Lip(1);
              end
          end
          %%%for panel 2 & 3 - top lip
          for ii=[2,3]
              for jj=8
                  VF_Matrix(ii,jj)=Panel_to_Top_Lip(2);
              end
          end 
          %%%roof/floor - top lip
          VF_Matrix(5,8)=F_Top_Lip_Roof*(A_Top_Lip/A_Roof);
          VF_Matrix(6,8)=F_Top_Lip_Floor*(A_Top_Lip/A_Roof);
          %%%top lip - panel 1-4
          for ii=8
              for jj=1:4
                  VF_Matrix(ii,jj)=VF_Matrix(jj,ii)*(A_Panel/A_Top_Lip);
              end
          end
          %%%top lip - roof/floor
          for ii=8
              for jj=5:6
                  VF_Matrix(ii,jj)=VF_Matrix(jj,ii)*(A_Roof/A_Top_Lip);
              end
          end
      end
      
      if Bottom_Lip > 0
      Panel_to_Bottom_Lip=[F_Panel_1_to_Panel(:,:,8), F_Panel_1_to_Panel(:,:,9)];
          %%%panel 1&4 to bottom lip
          for ii=[1,4]
              for jj=9
                  VF_Matrix(ii,jj)=Panel_to_Bottom_Lip(1);
              end
          end
          %%%panel 2&3 to bottom lip
          for ii=[2,3]
              for jj=9
                  VF_Matrix(ii,jj)=Panel_to_Bottom_Lip(2);
              end
          end
          %%%roof/floor to bottom lip
          VF_Matrix(5,9)=F_Bottom_Lip_Roof*(A_Bottom_Lip/A_Roof);
          VF_Matrix(6,9)=F_Bottom_Lip_Floor*(A_Bottom_Lip/A_Roof);
          %%%bottom lip - panel 1-4
          for ii=9
              for jj=1:4
                  VF_Matrix(ii,jj)=VF_Matrix(jj,ii)*(A_Panel/A_Bottom_Lip);
              end
          end
          %%%bottom lip - roof/floor
          for ii=9
              for jj=5:6
                  VF_Matrix(ii,jj)=VF_Matrix(jj,ii)*(A_Roof/A_Bottom_Lip);
              end
          end
      end
      VF_Matrix_Mod=VF_Matrix;
      if N_el_domain==9
          VF_Matrix(:,9)=0;
          VF_Matrix(9,:)=0;
          VF_Matrix(:,8)=0;
          VF_Matrix(8,:)=0;
          VF_Matrix(:,7)=0;
          VF_Matrix(7,:)=0;
          
          VF_Matrix(7,:)=VF_Matrix_Mod(8,:);
          VF_Matrix(:,7)=VF_Matrix_Mod(:,8);
          VF_Matrix(8,:)=VF_Matrix_Mod(9,:);
          VF_Matrix(:,8)=VF_Matrix_Mod(:,9);
          VF_Matrix(9,:)=VF_Matrix_Mod(7,:);
          VF_Matrix(:,9)=VF_Matrix_Mod(:,7);
      end
      if N_el_domain==8
          VF_Matrix(:,8)=0;
          VF_Matrix(8,:)=0;
          VF_Matrix(:,7)=0;
          VF_Matrix(7,:)=0;
          
          VF_Matrix(7,:)=VF_Matrix_Mod(8,:);
          VF_Matrix(:,7)=VF_Matrix_Mod(:,8);
          VF_Matrix(8,:)=VF_Matrix_Mod(7,:);
          VF_Matrix(:,8)=VF_Matrix_Mod(:,7);
      end
      xlswrite(Save_View_Factor_Filename,VF_Matrix)
end