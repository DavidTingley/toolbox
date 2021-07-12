function [W]=make_graph(data,graph_type,param)
% function [W]=make_graph(data,graph_type,param)
%
% data - view to be translated into a similarity graph
% graph_type - Variable determines how similarity graphs of each view are constructed
% param - optional value to manually pick threshold for graph construction


[num_points,num_features]=size(data);




switch graph_type
    case 'knn'
        data_distances=squareform(pdist(data,@naneucdist));
        data_distances_squared=data_distances.^2;
        A=zeros(num_points);
        k=param;
        for example=1:num_points
            row=data_distances_squared(:,example);
            [sorted_distances,indices]=sort(row);
            neighbors=indices(2:k+1);
            for neighbor=1:k
                A(example,neighbors(neighbor))=1;
            end
        end
        W=A;
        
    case 'gaussianKernel'
        data_distances=squareform(pdist(data,@naneucdist));
        data_distances_squared=data_distances.^2;
        if nargin < 3
            param = nanmean(nanmean(sqrt(data_distances_squared)));
        end
        sigma=param;
        gaussianKernel=exp(-data_distances_squared./(2*sigma));
        W=gaussianKernel-eye(num_points);
        
    case 'epsilonBall'
        data_distances=squareform(pdist(data,@naneucdist));
        data_distances_squared=data_distances.^2;
        A=zeros(num_points);
        epsilon=param;
        for example=1:num_points
            row=data_distances_squared(:,example);
            neighbors=find(row<epsilon);
            for neighbor=1:length(neighbors)
                if neighbors(neighbor)~=example
                    A(example,neighbors(neighbor))=1;
                end
            end
        end
        W=A;
        
    case 'cosine'
        A=zeros(num_points);
        for i=1:num_points
            for j=1:num_points
                A(i,j) = getCosineSimilarity(data(:,i),data(:,j));
                if isnan(A(i,j))
                    disp
                end
            end
        end
        W=A;
        
end
