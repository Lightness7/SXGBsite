function feature_PSA = Read_PSA(file_name, Seq_s)


Seq_Length = length(Seq_s);
fin = fopen(file_name, 'r');

feature_PSA = [];

n_l = zeros(3, 1);
n_l = n_l';
c_str = '';
str = fgetl(fin);	
ret = 0;
row_index = 1;
while feof(fin) == 0;
	row_index = row_index + 1;
	str = fgetl(fin);
	if ((row_index >= 3)&(row_index <= Seq_Length + 2))
		vect_v = [];
		A = strread(str, '%s', 'delimiter', ' ');
		
		c_str = [c_str cell2mat(A(2))];
		n_l(1) = str2num(cell2mat(A(4)));
		n_l(2) = str2num(cell2mat(A(5)));
		n_l(3) = str2num(cell2mat(A(6)));
		
		
		vect_v = n_l;
		feature_PSA = [feature_PSA; vect_v];
		ret = ret + 1;
	end
		
	

end

TF1 = strcmp(c_str, Seq_s);
if TF1 == 1
else
	out_pri = 'NO!'
end


fclose(fin);