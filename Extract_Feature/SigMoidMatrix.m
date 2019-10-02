function Map_M = SigMoidMatrix(M_Oringin,window_size,IS_start,max_v,min_v)

MM_V = zeros(size(M_Oringin,1),size(M_Oringin,2));

for i=1:size(M_Oringin,1)
	for j=1:size(M_Oringin,2)
		x_v = M_Oringin(i,j);
		%y_v = 1/(1+exp(-x_v));
		%y_v = (x_v - min_v)/(max_v - min_v);
		y_v = x_v;
		MM_V(i,j) = y_v;
	end
end

X = zeros(window_size,size(M_Oringin,2));
if IS_start==1
	st_a = window_size - size(M_Oringin,1) + 1;
	en_b = window_size;
	X(st_a:en_b,:) = MM_V;
else
	st_a = 1;
	en_b = size(M_Oringin,1);
	X(st_a:en_b,:) = MM_V;
end
%X(end,:)=[];
Map_M = X(:);
Map_M = Map_M';