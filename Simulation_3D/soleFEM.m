classdef soleFEM<handle
    %SOLEFEM Sole classe 3D
    
    properties(GetAccess=private)
        pname
        fname 
    end
    properties(Dependent)
        Noeudsinit
    end       
    properties
        K
        Ks
        C
        Cs
        coor
        dof
        dofsurf        
        connec
        analysis
        materials
        displ
        connecext
        TD
        trasl
        nTot
        nFree
        nFreeSurf
        nInt
        nodesDirichlet
        nodesDirichlet3
        nodesFreeSurf
        nodesFreeSurf3
        nodesInt
        nodesInt3
        IteOpt
        m_invKii_Kis  
        PSurfFree
        PSurfFree3
    end     
    methods
        function obj=soleFEM(pname,fname)
            obj.pname = pname;
            obj.fname = fname;
            [obj.analysis,obj.materials,obj.connec,obj.coor,...     % reads input file
             obj.dof,obj.displ,obj.TD]=read_input(obj.pname,obj.fname);
            Xmax = max(obj.coor(:,1));
            Xmin = min(obj.coor(:,1));
            Ymax = max(obj.coor(:,2));
            Ymin = min(obj.coor(:,2));
            Zmax = max(obj.coor(:,3));
            Zmin = min(obj.coor(:,3));
            obj.nTot = size(obj.coor,1);
            obj.trasl=[(Xmax-Xmin)/2,(Ymax-Ymin)/2,Zmax];
            obj.nodesDirichlet = find(obj.coor(:,3)==Zmax); % nodes of Dirichlet
            % Nodes of surface
            obj.nodesFreeSurf = find(obj.coor(:,1)==Xmin | obj.coor(:,1)==Xmax | obj.coor(:,2)==Ymin | obj.coor(:,2)==Ymax | obj.coor(:,3)==Zmin | obj.coor(:,3)==Zmax);
            obj.nodesFreeSurf(obj.nodesDirichlet)=[];
            % Nodes of volume - Internal nodes
            obj.nodesInt = find(obj.coor(:,1)~=Xmin & obj.coor(:,1)~=Xmax & obj.coor(:,2)~=Ymin & obj.coor(:,2)~=Ymax & obj.coor(:,3)~=Zmin & obj.coor(:,3)~=Zmax);
            obj.nInt = length(obj.nodesInt);
            obj.IteOpt = 1;
            obj.nFree = obj.analysis.NN - length(obj.nodesDirichlet);
            obj.nFreeSurf = length(obj.nodesFreeSurf);
            % Nodes in x;y;z;
            obj.nodesFreeSurf3(1:3:3*obj.nFreeSurf-2) = 3*obj.nodesFreeSurf-2;
            obj.nodesFreeSurf3(2:3:3*obj.nFreeSurf-1) = 3*obj.nodesFreeSurf-1;
            obj.nodesFreeSurf3(3:3:3*obj.nFreeSurf-0) = 3*obj.nodesFreeSurf-0;
            obj.nodesInt3(1:3:3*length(obj.nodesInt)-2) = 3*obj.nodesInt-2;
            obj.nodesInt3(2:3:3*length(obj.nodesInt)-1) = 3*obj.nodesInt-1;
            obj.nodesInt3(3:3:3*length(obj.nodesInt)-0) = 3*obj.nodesInt-0;
            obj.nodesDirichlet3(1:3:3*length(obj.nodesDirichlet)-2) = 3*obj.nodesDirichlet-2;
            obj.nodesDirichlet3(2:3:3*length(obj.nodesDirichlet)-1) = 3*obj.nodesDirichlet-1;
            obj.nodesDirichlet3(3:3:3*length(obj.nodesDirichlet)-0) = 3*obj.nodesDirichlet-0;           
        end
        function [obj]=stiffness(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % estimate of number of non-zero coefficients for Morse stocking
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            switch obj.analysis.Etag                         % coeffs of the relation between Nk and Ek
            %S.B modif des valeurs c1 et c2
            case{'4'}
              c1=3; c2=4;
            %S.B. fin de modif
            end
            Ek=c1*ones(obj.analysis.NN,1);               % initialisation of Ek
            for e=1:obj.analysis.NE,                     % loop over elements
              nodes=obj.connec(e,:);
              Ek(nodes)=Ek(nodes)+c2;                    % Ek is incremented on connectivity nodes
            end
            ncoeffs=9*sum(Ek);                           % sum of all the terms in Ek

            obj.K=spalloc(obj.analysis.neq,obj.analysis.neq,ncoeffs);% allocates sparse matrix
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % assemblage phase: system matrix and nodal forces due to imposed displ.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            D=3;                                           % problem dimensionality
            Dne=D*obj.analysis.ne;                         % number of nodal values in one surface el.
            for e=1:obj.analysis.NE,                       % stiffness matrix assemblage
              nodes = obj.connec(e,:);                     % element nodes
              T = obj.coor(nodes,:);                     % creates element
              Ke = stiff_linel_T4(T,obj.materials.A);
              dofe = reshape(obj.dof(nodes,:)',[1,Dne]);   % list of dof associated to element
              pe = find(dofe>0);                           % gets position of unknown displ. compon.
              Ie = dofe(pe);                               % gets value of associated DOFs 
              obj.K(Ie,Ie) = obj.K(Ie,Ie) + Ke(pe,pe);     % matrix assemblage
            end
            obj.K = -obj.K;
            obj.connecext = zeros(obj.analysis.nTD,3);
            for i = 1:obj.analysis.nTD
                obj.connecext(i,:) = obj.TD(i).nodes;
            end
            obj.IteOpt = obj.IteOpt + 1;        
            obj.C = inv(obj.K);
        end
        function [obj]=stiffnessSurface(obj)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Find the surface Stiffness Ks [Fi;Fs]=[Kii Kis;Ksi Kss]*[Xi;Xs]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            KInternalNodes = reshape(obj.dof(obj.nodesInt,:)',[1,3*length(obj.nodesInt)]);
            KSurfNodes = reshape(obj.dof(obj.nodesFreeSurf,:)',[1,3*length(obj.nodesFreeSurf)]);

            Kii = obj.K(KInternalNodes,KInternalNodes);
            Kis = obj.K(KInternalNodes,KSurfNodes);
            Ksi = obj.K(KSurfNodes,KInternalNodes);
            Kss = obj.K(KSurfNodes,KSurfNodes);

            %obj.Ks = (Kss-Ksi*(Kii\Kis);
            %obj.Ks = full(Kss);
            obj.Ks = full(obj.K(KSurfNodes,KSurfNodes));
            obj.dofsurf = obj.dof(obj.nodesFreeSurf,:);
            obj.m_invKii_Kis = full(Kii\Kis);
            obj.PSurfFree = obj.coor(obj.nodesFreeSurf,:);
            obj.PSurfFree3 = reshape(obj.PSurfFree([1:obj.nFreeSurf],:)',[1, 3*obj.nFreeSurf]);
            obj.Cs = full(obj.C(KSurfNodes,KSurfNodes));
        end
        function [dP,obj]=getdP(obj,dPSurf)
            dPSurf3 = reshape(dPSurf([1:obj.nFreeSurf],:)',[1, 3*obj.nFreeSurf]);
            dP = zeros(obj.nTot,3);
            dP(obj.nodesFreeSurf,:) = dPSurf;
            dPInt = -obj.m_invKii_Kis*dPSurf3';
            dP(obj.nodesInt,:) = reshape(dPInt,3,length(obj.nodesInt))';
        end        
%         function [RCR]=RCRSurface(obj,angle)
%             RBlock = [];
%             R = [1 0 0; 0 cos(angle) -sin(angle); 0 sin(angle) cos(angle)];              
%             for j=1:3:3*obj.nFreeSurf
%                 RBlock = blkdiag(RBlock,R);
%             end
%             RCR=RBlock*obj.Cs*RBlock';
%         end        
        
    end
end

