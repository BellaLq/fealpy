function [K_bearing_system]=K_bearing_system_produce(order_load,component_order_node,order_section_shaft,order_bearing_vector,bearing_layout,bearing_stiffness_matrix_spectrum)
%%%�������
%��У�����������⣩ȫ������ϵ�����ḽ�ӵĸնȾ���K_bearing_system��[Fx(N) Fy(N) Fz(N) Mx(Nmm) My(Nmm) Tz(Nmm)]'=K_bearing_system*[delta_x(mm) delta_y(mm) delta_z(mm) theta_x(rad) theta_y(rad) theta_z(rad)]'
%%%��������
%������ţ�order_load
%���μ��Ľڵ������Ϣ����component_order_node��[1 ���μ���� 2 �������ͣ�1��ʾ�ӵأ�2��ʾ����ᣩ 3 ����״̬��0��ʾδ���1��ʾ��� 4 װ��������0��ʾ������1��ʾ�������� 5 �ڵ�������n�� 6��n+5 �ڵ����]
%��ϵ������ż�������Ϣ����order_section_shaft����1��Ϊ����ţ���2��Ϊ��Ӧ��Ľ�����
%��У�����������⣩���������order_bearing_vector
%��ϵ�����װ���������bearing_layout
%�����/������� [1 ������ 2 ������Ȧ���������� 3 �����Ȧ�ڵ���� 4 ������ľ���Ȧ��������˾��루mm�� 5 ������Ȧ���������� 6 �����Ȧ�ڵ���� 7 ������ľ���Ȧ��������˾��루mm�� 8 ��ж���غ�(N)��capacity_dynamic 9 ��ж���غ�(N)��capacity_static 10-12 �� 13 ������������bearing_element_column 14 ������ͣ�bearing_type 15 �������(kg)��bearing_mass 16 ��п��(mm)��bearing_width 17 �����Ȧ(mm)��bearing_inner 18 �����Ȧ(mm)��bearing_outer 19 ����г�ʼ�Ӵ���/������нӴ���(rad)��contact_angle_initial 20 �������Ȧ�����ʣ�bearing_fi 21 �������Ȧ�����ʣ�bearing_fe 22 �����������bearing_element_number 23 ������ֱ��/׶��й������о�(mm)��bearing_element_diameter 24 ׶��й�����׶��(rad)��taper_angle 25 ������й����峤��(mm)��bearing_element_length 26 �������Բ�ǰ뾶(mm)��bearing_element_fillet 27 ��о�����϶(mm)��bearing_radial_clearance 28 ���������϶(mm)��bearing_axial_clearance 29 ��Ȧ����λ�ã�0��ʾ��ߣ�1��ʾ�ұߣ���bearing_rib_position 30 Բ��������Ȧ����Ƿ�����ߣ�bearing_rib_inner_left 31 Բ��������Ȧ�ұ��Ƿ�����ߣ�bearing_rib_inner_right 32 Բ��������Ȧ����Ƿ�����ߣ�bearing_rib_outer_left 33 Բ��������Ȧ�ұ��Ƿ�����ߣ�bearing_rib_outer_right 34 �������������ͣ�type_element_modi 35 �����Ȧ����ģ��(Mpa)��E_outer 36 �����Ȧ���ɱȣ�v_outer 37 ��й����嵯��ģ��(Mpa),E_element 38 ��й����岴�ɱȣ�v_element 39 �����Ȧ����ģ��(Mpa)��E_inner 40 �����Ȧ���ɱȣ�v_inner 41 ��������Ƭ����]
%������� [1 ������ 2 ������Ȧ���������� 3 �����Ȧ�ڵ���� 4 ������ľ���Ȧ��������˾��루mm�� 5 ������Ȧ���������� 6 �����Ȧ�ڵ���� 7 ������ľ���Ȧ��������˾��루mm�� 8 ��ж���غ�(N)��capacity_dynamic 9 ��ж���غ�(N)��capacity_static 10-13 �� 14 ������ͣ�bearing_type 15 �������(kg)��bearing_mass 16 ��п��(mm)��bearing_width 17 �����Ȧ(mm)��bearing_inner 18 �����Ȧ(mm)��bearing_outer 19-41 ��]
%ϵͳȫ����������е�ȫ������ϵ�ĸնȾ���bearing_stiffness_matrix_spectrum��[������� �����ţ�5*1�� ��иնȾ���5*5��]��[Fx(N) Fy(N) Fz(N) Mx(Nmm) My(Nmm)]'=��иնȾ���*[delta_x(mm) delta_y(mm) delta_z(mm) theta_x(rad) theta_y(rad)]'
%%%����ʼ
%������ϵ�����������������sum_shaft_segment
sum_shaft_segment=sum(order_section_shaft(2,:));
%%%������У�����������⣩���ḽ�ӵĸնȾ���K_bearing_system��[Fx(N) Fy(N) Fz(N) Mx(Nmm) My(Nmm) Tz(Nmm)]'=K_bearing_system*[delta_x(mm) delta_y(mm) delta_z(mm) theta_x(rad) theta_y(rad) theta_z(rad)]'
%������У�����������⣩���ḽ�ӵĸնȾ���K_bearing_system��[Fx(N) Fy(N) Fz(N) Mx(Nmm) My(Nmm) Tz(Nmm)]'=K_bearing_system*[delta_x(mm) delta_y(mm) delta_z(mm) theta_x(rad) theta_y(rad) theta_z(rad)]'
K_bearing_system=zeros(6*sum_shaft_segment,6*sum_shaft_segment);
if isempty(bearing_stiffness_matrix_spectrum)~=1
    %��¼��order_load������������е�ȫ������ϵ�ĸնȾ���bearing_stiffness_matrix_load��[������� �����ţ�5*1�� ��иնȾ���5*5��]��[Fx(N) Fy(N) Fz(N) Mx(Nmm) My(Nmm)]'=��иնȾ���*[delta_x(mm) delta_y(mm) delta_z(mm) theta_x(rad) theta_y(rad)]'
    id_bearing_load=(bearing_stiffness_matrix_spectrum(:,1)==order_load);  %������ţ�order_load
    bearing_stiffness_matrix_load=bearing_stiffness_matrix_spectrum(id_bearing_load,:);
    %%%��У�����������⣩������number_bearing
    size_bearing_vector=size(order_bearing_vector);  %��У�����������⣩���������order_bearing_vector
    number_bearing=size_bearing_vector(1,1);
    %������У�����������⣩���ḽ�ӵĸնȾ���K_bearing_system��[Fx(N) Fy(N) Fz(N) Mx(Nmm) My(Nmm) Tz(Nmm)]'=K_bearing_system*[delta_x(mm) delta_y(mm) delta_z(mm) theta_x(rad) theta_y(rad) theta_z(rad)]'
    if number_bearing>0  %��У�����������⣩������number_bearing
        for ii=1:number_bearing  %��У�����������⣩������number_bearing
            %��¼��ii����У�����������⣩����Ĳ���
            %��¼��ii����У�����������⣩����ţ�order_bearing
            order_bearing=order_bearing_vector(ii,1);
            %��¼��ii����У�����������⣩װ���������bearing_layout_temp
            %�����/������� [1 ������ 2 ������Ȧ���������� 3 �����Ȧ�ڵ���� 4 ������ľ���Ȧ��������˾��루mm�� 5 ������Ȧ���������� 6 �����Ȧ�ڵ���� 7 ������ľ���Ȧ��������˾��루mm�� 8 ��ж���غ�(N)��capacity_dynamic 9 ��ж���غ�(N)��capacity_static 10-12 �� 13 ������������bearing_element_column 14 ������ͣ�bearing_type 15 �������(kg)��bearing_mass 16 ��п��(mm)��bearing_width 17 �����Ȧ(mm)��bearing_inner 18 �����Ȧ(mm)��bearing_outer 19 ����г�ʼ�Ӵ���/������нӴ���(rad)��contact_angle_initial 20 �������Ȧ�����ʣ�bearing_fi 21 �������Ȧ�����ʣ�bearing_fe 22 �����������bearing_element_number 23 ������ֱ��/׶��й������о�(mm)��bearing_element_diameter 24 ׶��й�����׶��(rad)��taper_angle 25 ������й����峤��(mm)��bearing_element_length 26 �������Բ�ǰ뾶(mm)��bearing_element_fillet 27 ��о�����϶(mm)��bearing_radial_clearance 28 ���������϶(mm)��bearing_axial_clearance 29 ��Ȧ����λ�ã�0��ʾ��ߣ�1��ʾ�ұߣ���bearing_rib_position 30 Բ��������Ȧ����Ƿ�����ߣ�bearing_rib_inner_left 31 Բ��������Ȧ�ұ��Ƿ�����ߣ�bearing_rib_inner_right 32 Բ��������Ȧ����Ƿ�����ߣ�bearing_rib_outer_left 33 Բ��������Ȧ�ұ��Ƿ�����ߣ�bearing_rib_outer_right 34 �������������ͣ�type_element_modi 35 �����Ȧ����ģ��(Mpa)��E_outer 36 �����Ȧ���ɱȣ�v_outer 37 ��й����嵯��ģ��(Mpa),E_element 38 ��й����岴�ɱȣ�v_element 39 �����Ȧ����ģ��(Mpa)��E_inner 40 �����Ȧ���ɱȣ�v_inner 41 ��������Ƭ����]
            id_bearing=(bearing_layout(:,1)==order_bearing);  %��ii����У�����������⣩����ţ�order_bearing
            %��¼��ii����У�����������⣩����Ȧ����������
            order_inner=bearing_layout(id_bearing,2);
            %��¼��ii����У�����������⣩����Ȧ�ڵ�����
            node_inner=bearing_layout(id_bearing,3);
            %��¼��ii����У�����������⣩����Ȧ����������
            order_outer=bearing_layout(id_bearing,5);
            %��¼��ii����У�����������⣩����Ȧ�ڵ�����
            node_outer=bearing_layout(id_bearing,6);
            %%%�����ii�������Ȧ�ڵ�������������λ�ã�Nbs_inner
            if order_inner>0  %�����Ȧ������
                position_inner=find(roundn(order_section_shaft(1,:),-8)==roundn(order_inner,-8));  %��ϵ������ż�������Ϣ����order_section_shaft����1��Ϊ����ţ���2��Ϊ��Ӧ��Ľ�����
                if position_inner==1
                    Nbs_inner=node_inner-1;  %��ii�������Ȧ�ڵ�������������λ��
                else
                    Nbs_inner=sum(order_section_shaft(2,1:position_inner-1))+node_inner-1;  %��ii�������Ȧ�ڵ�������������λ��
                end
            elseif order_inner<0  %�����Ȧ�������μ�
                %%%���������Ȧ�ڵ������μ���λ�ã�position_node_inner
                %��¼���Ϊorder_inner���μ��Ľڵ������Ϣ��component_order_node_inner
                id_component_inner=(roundn(component_order_node(:,1),-8)==roundn(order_inner,-8));  %���μ��Ľڵ������Ϣ����component_order_node��[1 ���μ���� 2 �������ͣ�1��ʾ�ӵأ�2��ʾ����ᣩ 3 ����״̬��0��ʾδ���1��ʾ��� 4 װ��������0��ʾ������1��ʾ�������� 5 �ڵ�������n�� 6��n+5 �ڵ����]
                %��¼���μ��Ľڵ�������number_node_inner
                number_node_inner=component_order_node(id_component_inner,5);
                %��¼���μ��Ľڵ������������component_order_node_inner
                component_order_node_inner=component_order_node(id_component_inner,6:number_node_inner+5);
                %���������Ȧ�ڵ������μ���λ�ã�position_node_inner
                position_node_inner=find(component_order_node_inner(1,:)==node_inner);  %��ii����У�����������⣩����Ȧ�ڵ�����
                %%%���������Ȧ�ڵ���������ϵ�ڵ��λ�ã�position_inner
                position_inner=find(roundn(order_section_shaft(1,:),-8)==roundn(order_inner,-8));  %��ϵ����/���μ���ż��ڵ���Ϣ����order_section_shaft����1��Ϊ��/���μ���ţ���2��Ϊ��Ӧ��/���μ��Ľڵ�����
                if position_inner==1
                    Nbs_inner=position_node_inner-1;
                else
                    Nbs_inner=sum(order_section_shaft(2,1:position_inner-1))+position_node_inner-1;  %��ϵ����/���μ���ż��ڵ���Ϣ����order_section_shaft����1��Ϊ��/���μ���ţ���2��Ϊ��Ӧ��/���μ��Ľڵ�����
                end
            end
            %%%
            %%%�����ii�������Ȧ�ڵ�������������λ�ã�Nbs_outer
            if order_outer>0  %�����Ȧ������
                position_outer=find(roundn(order_section_shaft(1,:),-8)==roundn(order_outer,-8));  %��ϵ������ż�������Ϣ����order_section_shaft����1��Ϊ����ţ���2��Ϊ��Ӧ��Ľ�����
                if position_outer==1
                    Nbs_outer=node_outer-1;  %��ii�������Ȧ�ڵ�������������λ��
                else
                    Nbs_outer=sum(order_section_shaft(2,1:position_outer-1))+node_outer-1;  %��ii�������Ȧ�ڵ�������������λ��
                end
            elseif order_outer<0  %�����Ȧ�������μ�
                %%%���������Ȧ�ڵ������μ���λ�ã�position_node_outer
                %��¼���Ϊorder_outer���μ��Ľڵ������Ϣ��component_order_node_outer
                id_component_outer=(roundn(component_order_node(:,1),-8)==roundn(order_outer,-8));  %���μ��Ľڵ������Ϣ����component_order_node��[1 ���μ���� 2 �������ͣ�1��ʾ�ӵأ�2��ʾ����ᣩ 3 ����״̬��0��ʾδ���1��ʾ��� 4 װ��������0��ʾ������1��ʾ�������� 5 �ڵ�������n�� 6��n+5 �ڵ����]
                %��¼���μ��Ľڵ�������number_node_outer
                number_node_outer=component_order_node(id_component_outer,5);
                %��¼���μ��Ľڵ������������component_order_node_outer
                component_order_node_outer=component_order_node(id_component_outer,6:number_node_outer+5);
                %���������Ȧ�ڵ������μ���λ�ã�position_node_outer
                position_node_outer=find(component_order_node_outer(1,:)==node_outer);  %��ii����У�����������⣩����Ȧ�ڵ�����
                %%%���������Ȧ�ڵ���������ϵ�ڵ��λ�ã�position_outer
                position_outer=find(roundn(order_section_shaft(1,:),-8)==roundn(order_outer,-8));  %��ϵ����/���μ���ż��ڵ���Ϣ����order_section_shaft����1��Ϊ��/���μ���ţ���2��Ϊ��Ӧ��/���μ��Ľڵ�����
                if position_outer==1
                    Nbs_outer=position_node_outer-1;
                else
                    Nbs_outer=sum(order_section_shaft(2,1:position_outer-1))+position_node_outer-1;  %��ϵ����/���μ���ż��ڵ���Ϣ����order_section_shaft����1��Ϊ��/���μ���ţ���2��Ϊ��Ӧ��/���μ��Ľڵ�����
                end
            end
            %��¼��ii����У�����������⣩�ĸնȾ���bearing_stiffness_matrix��[Fx(N) Fy(N) Fz(N) Mx(Nmm) My(Nmm)]'=bearing_stiffness_matrix*[delta_x(mm) delta_y(mm) delta_z(mm) theta_x(rad) theta_y(rad)]'
            %��order_load������������е�ȫ������ϵ�ĸնȾ���bearing_stiffness_matrix_load��[������� �����ţ�5*1�� ��иնȾ���5*5��]��[Fx(N) Fy(N) Fz(N) Mx(Nmm) My(Nmm)]'=��иնȾ���*[delta_x(mm) delta_y(mm) delta_z(mm) theta_x(rad) theta_y(rad)]'
            id_bearing_stiffness=(roundn(bearing_stiffness_matrix_load(:,2),-8)==roundn(order_bearing,-8));  %��ii����У�����������⣩����ţ�order_bearing
            bearing_stiffness_matrix=bearing_stiffness_matrix_load(id_bearing_stiffness,3:7);
            %%%�����ii����У�����������⣩���ḽ�ӵĸնȾ���K_bearing_system_temp��[Fx(N) Fy(N) Fz(N) Mx(Nmm) My(Nmm) Tz(Nmm)]'=K_bearing_system_temp*[delta_x(mm) delta_y(mm) delta_z(mm) theta_x(rad) theta_y(rad) theta_z(rad)]'
            if (order_inner~=0)&&(order_outer==0)  %�����Ȧ���ӵأ���Ȧ�ӵ�
                K_bearing_system(6*Nbs_inner+1:6*Nbs_inner+5,6*Nbs_inner+1:6*Nbs_inner+5)=K_bearing_system(6*Nbs_inner+1:6*Nbs_inner+5,6*Nbs_inner+1:6*Nbs_inner+5)+bearing_stiffness_matrix;
            elseif (order_inner==0)&&(order_outer~=0)  %�����Ȧ�ӵأ���Ȧ���ӵ�
                K_bearing_system(6*Nbs_outer+1:6*Nbs_outer+5,6*Nbs_outer+1:6*Nbs_outer+5)=K_bearing_system(6*Nbs_outer+1:6*Nbs_outer+5,6*Nbs_outer+1:6*Nbs_outer+5)+bearing_stiffness_matrix;
            else  %�����Ȧ���ӵأ���Ȧ���ӵ�
                K_bearing_system(6*Nbs_inner+1:6*Nbs_inner+5,6*Nbs_inner+1:6*Nbs_inner+5)=K_bearing_system(6*Nbs_inner+1:6*Nbs_inner+5,6*Nbs_inner+1:6*Nbs_inner+5)+bearing_stiffness_matrix;
                K_bearing_system(6*Nbs_outer+1:6*Nbs_outer+5,6*Nbs_outer+1:6*Nbs_outer+5)=K_bearing_system(6*Nbs_outer+1:6*Nbs_outer+5,6*Nbs_outer+1:6*Nbs_outer+5)+bearing_stiffness_matrix;
                K_bearing_system(6*Nbs_inner+1:6*Nbs_inner+5,6*Nbs_outer+1:6*Nbs_outer+5)=K_bearing_system(6*Nbs_inner+1:6*Nbs_inner+5,6*Nbs_outer+1:6*Nbs_outer+5)-bearing_stiffness_matrix;
                K_bearing_system(6*Nbs_outer+1:6*Nbs_outer+5,6*Nbs_inner+1:6*Nbs_inner+5)=K_bearing_system(6*Nbs_outer+1:6*Nbs_outer+5,6*Nbs_inner+1:6*Nbs_inner+5)-bearing_stiffness_matrix;
            end
        end
    end
    %%%ǿ�ƶԳƾ���
    K_bearing_system=0.5*(K_bearing_system+K_bearing_system');
end
%%%�������