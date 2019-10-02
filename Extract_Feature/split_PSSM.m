function [SEQ,n_seq,pssm_vector] = split_PSSM(protein,window_radius,PSSM_M)


SEQ = [];pssm_vector = [];
n_seq = 0;
protein_s=cell2mat(protein);
n_seq = size(protein_s, 2);

X_str = 'XXXXXXXXXXX';
length_norm = 2*window_radius + 1;

max_k = max(PSSM_M);
max_v = max(max_k);

min_k = min(PSSM_M);
min_v = min(min_k);

for i = 1:n_seq
	fd = i;
	v_v_v = [];
	FF = [];
	if ((fd-window_radius) > 0)&((fd+window_radius) < n_seq)
		ppp = protein_s(fd-window_radius:fd+window_radius);
		v_v_v = PSSM_M(fd-window_radius:fd+window_radius, :);
		FF = SigMoidMatrix(v_v_v, length_norm, 0, max_v, min_v);
	elseif ((fd-window_radius) <= 0)
		N_need_X = window_radius - fd + 1;
		need_X_str = [];
		need_X_str = X_str(1:N_need_X);
		ppp = protein_s(1:fd+window_radius);
		ppp = [need_X_str ppp];
		
		v_v_v = PSSM_M(1:fd+window_radius, :);
		FF = SigMoidMatrix(v_v_v,length_norm, 1, max_v, min_v);
	else
		N_need_X = window_radius+fd-n_seq;
		need_X_str = [];
		need_X_str = X_str(1:N_need_X);
		
		ppp = protein_s(fd-window_radius:n_seq);
		ppp = [ppp need_X_str];
		
		v_v_v = PSSM_M(fd-window_radius:n_seq, :);
		FF = SigMoidMatrix(v_v_v,length_norm, 0, max_v, min_v);
	end
	kkk = cellstr(ppp);SEQ = [SEQ; kkk];
	pssm_vector = [pssm_vector; FF];
	
end






