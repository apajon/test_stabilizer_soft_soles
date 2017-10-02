classdef soleFEM_newStiff<handle
    %SOLEFEM Sole classe 3D
       
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
        nodesSurf
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
        zpoles
        lx_foot
        ly_foot
        lz_foot
        Kid
        Kii
        Kis
        der_K
        der_A
        der_C
        mu
        lambda
    end     
    methods
        function obj=soleFEM_newStiff(pname,fname,l1,L1,e1)
            obj.lx_foot = L1;
            obj.ly_foot = l1;
            obj.lz_foot = e1;
            [obj.elements_surf,obj.elements_vol,obj.coor] = input_mesh(pname,fname);     % reads nodes and elem. from GMSH file
            obj.nEle = size(obj.elements_surf,1) + size(obj.elements_vol,1);
            nodes = (1:1:size(obj.coor,1))';
            obj.nTot = length(nodes);
            obj.trasl = [obj.lx_foot/2,0,obj.lz_foot];
            obj.zpoles = 0.0;
            % Dirichlet nodes
            obj.nodesDirichlet = find(obj.coor(:,3)==obj.lz_foot); 
            obj.nDirichlet = length(obj.nodesDirichlet);
            % Surface nodes
            obj.nodesSurf = unique(obj.elements_surf,'first');
            % Free surface nodes
            obj.nodesFreeSurf = setdiff(obj.nodesSurf,obj.nodesDirichlet);
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
        function [obj] = setMaterial(obj,Young,Poisson)
            obj.E = Young;
            obj.nu = Poisson;
        end
        function [obj] = stiffness(obj)
            obj.mu=obj.E/(2*(1+obj.nu));obj.lambda=obj.E*obj.nu/((1+obj.nu)*(1-2*obj.nu));
            A = sparse(3*obj.nTot,3*obj.nTot);
            % element data file for tetraeder: 1-node/2-node/3-node/4-node
            for j = 1:size(obj.elements_vol,1)
              I = 3*obj.elements_vol(j,[1,1,1,2,2,2,3,3,3,4,4,4])-[2,1,0,2,1,0,2,1,0,2,1,0];
              A(I,I) = A(I,I) + stima(obj.coor(obj.elements_vol(j,:),:),obj.lambda,obj.mu); % this is the stiffness matrix
            end 
            obj.K = A(obj.nodesFree3,obj.nodesFree3);
            obj.Kid = A(obj.nodesInt3,obj.nodesDirichlet3);
            obj.IteOpt = obj.IteOpt + 1;
            obj.C = inv(obj.K);
        end
        function [obj] = derStiff(obj)
            obj.der_K = cell(3*obj.nFree,1);
            obj.der_A = cell(3*obj.nTot,1);
            for i = 1:(3*obj.nTot)
                obj.der_A{i}=sparse(3*obj.nTot,3*obj.nTot);
            end
            for i = 1:(3*obj.nFree)
                obj.der_K{i}=sparse(3*obj.nFree,3*obj.nFree);
            end
            
            for j = 1:size(obj.elements_vol,1)
                I = 3*obj.elements_vol(j,[1,1,1,2,2,2,3,3,3,4,4,4])-[2,1,0,2,1,0,2,1,0,2,1,0];
                obj.der_A{I(1),1}(I,I) = obj.der_A{I(1),1}(I,I) + der_stima_dx1(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(2),1}(I,I) = obj.der_A{I(2),1}(I,I) + der_stima_dy1(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(3),1}(I,I) = obj.der_A{I(3),1}(I,I) + der_stima_dz1(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(4),1}(I,I) = obj.der_A{I(4),1}(I,I) + der_stima_dx2(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(5),1}(I,I) = obj.der_A{I(5),1}(I,I) + der_stima_dy2(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(6),1}(I,I) = obj.der_A{I(6),1}(I,I) + der_stima_dz2(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(7),1}(I,I) = obj.der_A{I(7),1}(I,I) + der_stima_dx3(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(8),1}(I,I) = obj.der_A{I(8),1}(I,I) + der_stima_dy3(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(9),1}(I,I) = obj.der_A{I(9),1}(I,I) + der_stima_dz3(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(10),1}(I,I) = obj.der_A{I(10),1}(I,I) + der_stima_dx4(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(11),1}(I,I) = obj.der_A{I(11),1}(I,I) + der_stima_dy4(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
                obj.der_A{I(12),1}(I,I) = obj.der_A{I(12),1}(I,I) + der_stima_dz4(obj.lambda,obj.mu,obj.coor(obj.elements_vol(j,:),:));
            end
            for i = 1:(3*obj.nFree)
                obj.der_K{i} = obj.der_A{obj.nodesFree3(i)}(obj.nodesFree3,obj.nodesFree3);
            end
            obj.der_C = cell(3*obj.nFree,1);
            for i = 1:(3*obj.nFree)
                obj.der_C{i} = sparse(3*obj.nFree,3*obj.nFree);
            end
            for i = 1:(3*obj.nFree)
                obj.der_C{i} = -obj.C * obj.der_K{i} * obj.C;
            end            
            
        end
        function [obj]=stiffnessSurface(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Find the surface Stiffness Ks [Fi;Fs]=[Kii Kis;Ksi Kss]*[Xi;Xs]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            KInternalNodes = reshape(obj.dof(obj.nodesInt,:)',[1,3*obj.nInt]);
            KNodesFreeSurf = reshape(obj.dof(obj.nodesFreeSurf,:)',[1,3*obj.nFreeSurf]);
            
            obj.Kii = obj.K(KInternalNodes,KInternalNodes);
            obj.Kis = obj.K(KInternalNodes,KNodesFreeSurf);
            Ksi = obj.K(KNodesFreeSurf,KInternalNodes);
            Kss = obj.K(KNodesFreeSurf,KNodesFreeSurf);
            
            obj.m_invKii_Kis = full(obj.Kii\obj.Kis);
            obj.Cs = full(obj.C(KNodesFreeSurf,KNodesFreeSurf));
            obj.Ks = full(obj.K(KNodesFreeSurf,KNodesFreeSurf));
        end
        function [dP,obj]=getdP(obj,dPFreeSurf,dPDir)
            dPSurf3 = reshape(dPFreeSurf([1:obj.nFreeSurf],:)',[1, 3*obj.nFreeSurf]);
            dP = zeros(obj.nTot,3);
            dP(obj.nodesFreeSurf,:) = dPFreeSurf;
            dP(obj.nodesDirichlet,:) = dPDir;
            if ~isempty(obj.m_invKii_Kis)
                dPInt = -obj.m_invKii_Kis*dPSurf3';
                dP(obj.nodesInt,:) = reshape(dPInt,3,length(obj.nodesInt))';
            end
        end        
        function [dP,obj]=getdPTot(obj,dPFreeSurf,dPDir)
            dPFreeSurf3 = reshape(dPFreeSurf([1:obj.nFreeSurf],:)',[1, 3*obj.nFreeSurf]);
            dPDir3 = reshape(dPDir([1:obj.nDirichlet],:)',[1, 3*obj.nDirichlet]);
            dP = zeros(obj.nTot,3);
            dP(obj.nodesFreeSurf,:) = dPFreeSurf;
            dP(obj.nodesDirichlet,:) = dPDir;
            dPInt3 = -inv(obj.Kii)*(obj.Kid*dPDir3' + obj.Kis*dPFreeSurf3');
            dP(obj.nodesInt,:) = reshape(dPInt3,3,length(obj.nodesInt))';
        end          
    end
end