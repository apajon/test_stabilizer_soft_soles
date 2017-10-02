%classdef soleFEM_newStiff<handle
classdef soleFEM_newStiff_nohandle
    %SOLEFEM Sole classe 3D
    
%     properties(GetAccess=private)
%         pname
%         fname 
%     end
%     properties(Dependent)
%         Noeudsinit
%     end       
    properties
        K
        Ks
        C
        Cs
        coor
        elements_surf
        elements_vol
        nEle
        trasl
        nTot        
        nodesDirichlet
        nodesDirichlet3
        nDirichlet
        nodesFreeSurf
        nodesFreeSurf3
        nFreeSurf
        nodesInt
        nodesInt3
        nInt
        nodesFree
        nodesFree3
        nFree
        IteOpt
        m_invKii_Kis
        dof
        E
        nu
    end     
    methods
        function obj=soleFEM_newStiff_nohandle(pname,fname)           
            [obj.elements_surf,obj.elements_vol,obj.coor] = input_mesh(pname,fname);     % reads nodes and elem. from GMSH file
            lx_foot = 0.23;
            lz_foot = 0.03;
            obj.nEle = size(obj.elements_surf,1) + size(obj.elements_vol,1);
            nodes = (1:1:size(obj.coor,1))';
            obj.nTot = length(nodes);
            obj.trasl = [lx_foot/2,0,lz_foot];
            % Dirichlet nodes
            obj.nodesDirichlet = find(obj.coor(:,3)==lz_foot); 
            obj.nDirichlet = length(obj.nodesDirichlet);
            % Free surface nodes
            % cont_surf_free = 0;
            % obj.elements_surf_free = obj.elements_surf;
            % ind_tot = [];
            % for i=1:obj.nDirichlet
            %     ind1 = [];
            %     ind2 = [];
            %     ind3 = [];
            %     ind1 = find(obj.elements_surf(:,1)==obj.nodesDirichlet(i));
            %     ind2 = find(obj.elements_surf(:,2)==obj.nodesDirichlet(i));
            %     ind3 = find(obj.elements_surf(:,3)==obj.nodesDirichlet(i));
            %     if ~isempty(ind1)
            %         ind_tot = [ind_tot; ind1];
            %     end
            %     if ~isempty(ind2)
            %         ind_tot = [ind_tot; ind2];
            %     end
            %     if ~isempty(ind3)
            %         ind_tot = [ind_tot; ind3];
            %     end                
            % end
            % ind_tot = unique(ind_tot,'first');
            % obj.elements_surf_free(ind_tot,:)=[];
            nodeSurf = unique(obj.elements_surf,'first');
            obj.nodesFreeSurf = setdiff(nodeSurf,obj.nodesDirichlet);
            obj.nFreeSurf = length(obj.nodesFreeSurf);
            % Nodes of volume - Internal nodes
            obj.nodesInt = nodes;
            obj.nodesInt([obj.nodesFreeSurf; obj.nodesDirichlet])=[];
            obj.nInt = length(obj.nodesInt);
            % Nodes free - Free surface nodes + Internal nodes
            obj.nodesFree = nodes;
            obj.nodesFree(obj.nodesDirichlet)=[];
            obj.nFree = length(obj.nodesFree);
            % Compute dof matrix
            cont = 0;
            obj.dof = zeros(size(obj.coor));
            contDir = 1;
            for i=1:obj.nTot
                if contDir <= obj.nDirichlet
                    if i~=obj.nodesDirichlet(contDir)
                        for j=1:3
                            cont = cont + 1;
                            obj.dof(i,j) = cont;
                        end
                    else
                        contDir = contDir + 1;
                    end
                else
                     for j=1:3
                        cont = cont + 1;
                        obj.dof(i,j) = cont;
                    end                   
                end
            end
            % Iteration
            obj.IteOpt = 1;
            % Nodes in x;y;z;
            obj.nodesFree3(1:3:3*obj.nFree-2) = 3*obj.nodesFree-2;
            obj.nodesFree3(2:3:3*obj.nFree-1) = 3*obj.nodesFree-1;
            obj.nodesFree3(3:3:3*obj.nFree-0) = 3*obj.nodesFree-0;            
            obj.nodesFreeSurf3(1:3:3*obj.nFreeSurf-2) = 3*obj.nodesFreeSurf-2;
            obj.nodesFreeSurf3(2:3:3*obj.nFreeSurf-1) = 3*obj.nodesFreeSurf-1;
            obj.nodesFreeSurf3(3:3:3*obj.nFreeSurf-0) = 3*obj.nodesFreeSurf-0;
            obj.nodesInt3(1:3:3*obj.nInt-2) = 3*obj.nodesInt-2;
            obj.nodesInt3(2:3:3*obj.nInt-1) = 3*obj.nodesInt-1;
            obj.nodesInt3(3:3:3*obj.nInt-0) = 3*obj.nodesInt-0;
            obj.nodesDirichlet3(1:3:3*obj.nDirichlet-2) = 3*obj.nodesDirichlet-2;
            obj.nodesDirichlet3(2:3:3*obj.nDirichlet-1) = 3*obj.nodesDirichlet-1;
            obj.nodesDirichlet3(3:3:3*obj.nDirichlet-0) = 3*obj.nodesDirichlet-0;          
        end
        function [obj]=setMaterial(obj,Young,Poisson)
            obj.E = Young;
            obj.nu = Poisson;
        end
        function [obj]=stiffness(obj)
            mu=obj.E/(2*(1+obj.nu));lambda=obj.E*obj.nu/((1+obj.nu)*(1-2*obj.nu));
            A = sparse(3*obj.nTot,3*obj.nTot);
            % obj.connec % element data file for tetraeder: 1-node/2-node/3-node/4-node
            for j = 1:size(obj.elements_vol,1)
              I = 3*obj.elements_vol(j,[1,1,1,2,2,2,3,3,3,4,4,4])-[2,1,0,2,1,0,2,1,0,2,1,0];
              A(I,I) = A(I,I) + stima(obj.coor(obj.elements_vol(j,:),:),lambda,mu); % this is the stiffness matrix
            end 
            obj.K = A(obj.nodesFree3,obj.nodesFree3);
            obj.IteOpt = obj.IteOpt + 1;        
            obj.C = inv(obj.K);
        end  
        function [obj]=stiffnessSurface(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Find the surface Stiffness Ks [Fi;Fs]=[Kii Kis;Ksi Kss]*[Xi;Xs]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            KInternalNodes = reshape(obj.dof(obj.nodesInt,:)',[1,3*obj.nInt]);
            KSurfNodes = reshape(obj.dof(obj.nodesFreeSurf,:)',[1,3*obj.nFreeSurf]);
            
            Kii = obj.K(KInternalNodes,KInternalNodes);
            Kis = obj.K(KInternalNodes,KSurfNodes);
            Ksi = obj.K(KSurfNodes,KInternalNodes);
            Kss = obj.K(KSurfNodes,KSurfNodes);
            
            obj.m_invKii_Kis = full(Kii\Kis);
            obj.Cs = full(obj.C(KSurfNodes,KSurfNodes));
            obj.Ks = full(obj.K(KSurfNodes,KSurfNodes));
        end
        function [dP,obj]=getdP(obj,dPSurf)
            dPSurf3 = reshape(dPSurf([1:obj.nFreeSurf],:)',[1, 3*obj.nFreeSurf]);
            dP = zeros(obj.nTot,3);
            dP(obj.nodesFreeSurf,:) = dPSurf;
            dPInt = -obj.m_invKii_Kis*dPSurf3';
            dP(obj.nodesInt,:) = reshape(dPInt,3,length(obj.nodesInt))';
        end        
    end
end