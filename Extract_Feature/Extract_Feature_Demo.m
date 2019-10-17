

% feature extraction of training process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load label of data set
load('D:\Data\GTP\label_train.mat');

% load dataset in 'XXX.fasta' format
[HEADER, SEQ] = fastaread('D:\Data\GTP\GTP_TRAIN.fasta');
SEQ = SEQ';
HEADER=HEADER';
Protein_num = size(SEQ,1);

% Read PSSM files (for example, 'XXX.pssm') in the document
PSSM_path = 'D:\Data\GTP\PSSM_TRAIN\';
PSSM_file = '.pssm';
PSSM_path_file = [];
for i=1:Protein_num
path1 = [PSSM_path num2str(i-1) PSSM_file];
path2 = cellstr(path1);
PSSM_path_file = [PSSM_path_file; path2];
end

% The PSSM feature of contiguous residues is selected by 8 adjacent
% residues on both sides of the target amino acid residue with a window size of 8*2 + 1
PSSM_with_window=[];
window = 8;
window_size = window*2 + 1;
for i=1:Protein_num
Sequence_i = SEQ(i);
file_i = cell2mat(PSSM_path_file(i));
protein_i = cell2mat(Sequence_i);
Sequence_i_num = size(protein_i, 2);
PSSM_Matrix = Read_PSSM(file_i, Sequence_i_num);
[feature, sequence_num, PSSM_i] = split_PSSM(Sequence_i, window, PSSM_Matrix);
PSSM_with_window = [PSSM_with_window; PSSM_i];
end

% DCT compresses PSSM features of continuous residues and selects the first
% 9 rows of PSSM-DCT matrix as PSSM-DCT features
% DCT
PSSM_DCT_feature=[];
for i=1:size(PSSM_with_window,1)
	PSSM_i = [];
	PSSM_i = PSSM_with_window(i,:);
	PSSM_w_Matrix = reshape(PSSM_i,window_size,20);
	PSSM_DCT_i = dct2(PSSM_w_Matrix);
	PSSM_DCT_i = PSSM_DCT_i(1:9,:);i
	PSSM_DCT_feature(i, :) = PSSM_DCT_i(:);	
end

% DWT compresses PSSM features of continuous residues and obtains PSSM-DWT features
% DWT
% PSSM_DWT_feature=[];
% for i=1:size(PSSM_with_window,1)
% 	PSSM_ii = [];
% 	PSSM_ii = PSSM_with_window(i,:);
%     PSSM_DWT_ii = [];
% 	PSSM_DWT_ii = EstraggoConDWT(PSSM_ii,4);
% 	PSSM_DWT_feature = [PSSM_DWT_feature;PSSM_DWT_ii];
% end

% PSA, load the PSA data sets in the document, the format is 'XXX.psa'
% PSA
PSSM_path = 'D:\Data\GTP\GTP_Train_PSA';
PSSM_file = '.psa';
PSA_path_file = [];
for i=1:Protein_num 
use_head_str = cell2mat(HEADER(i,1));
use_head_str(5)=[];
use_head_str(6)=[];
path1 = [PSSM_path '\' use_head_str PSSM_file];
Path_p = cellstr(path1);
PSA_path_file = [PSA_path_file; Path_p];
end

PSA=[];
for i=1:Protein_num
Sequence_i = SEQ(i);
protein_i=cell2mat(Sequence_i);
PSA_path_file2 = cell2mat(PSA_path_file(i));
PSA_Matrix = Read_PSA(PSA_path_file2, protein_i);
PSA=[PSA;PSA_Matrix];i
end


[Train_PSSM_DCT_feature, ia, ic] = unique(PSSM_DCT_feature, 'rows');
Train_label_y = label_train(ia, :);
Train_PSA_feature = PSA(ia, :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% feature extraction of testing process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load label of data set
load('D:\Data\GTP\label_test.mat');

% load dataset in 'XXX.fasta' format
[HEADER, SEQ] = fastaread('D:\Data\GTP\GTP_TEST.fasta');
SEQ = SEQ';
HEADER=HEADER';
Protein_num = size(SEQ,1);

% Read PSSM files (for example, 'XXX.pssm') in the document
PSSM_path = 'D:\Data\GTP\PSSM_TEST\';
PSSM_file = '.pssm';
PSSM_path_file = [];
for i=1:Protein_num
path1 = [PSSM_path num2str(i-1) PSSM_file];
path2 = cellstr(path1);
PSSM_path_file = [PSSM_path_file; path2];
end

% The PSSM feature of contiguous residues is selected by 8 adjacent
% residues on both sides of the target amino acid residue with a window size of 8*2 + 1
PSSM_with_window=[];
window = 8;
window_size = window*2 + 1;
for i=1:Protein_num
Sequence_i = SEQ(i);
file_i = cell2mat(PSSM_path_file(i));
protein_i = cell2mat(Sequence_i);
Sequence_i_num = size(protein_i, 2);
PSSM_Matrix = Read_PSSM(file_i, Sequence_i_num);
[feature, sequence_num, PSSM_i] = split_PSSM(Sequence_i, window, PSSM_Matrix);
PSSM_with_window = [PSSM_with_window; PSSM_i];
end

% DCT compresses PSSM features of continuous residues and selects the first
% 9 rows of PSSM-DCT matrix as PSSM-DCT features
% DCT
PSSM_DCT_feature=[];
for i=1:size(PSSM_with_window,1)
	PSSM_i = [];
	PSSM_i = PSSM_with_window(i,:);
	PSSM_w_Matrix = reshape(PSSM_i,window_size,20);
	PSSM_DCT_i = dct2(PSSM_w_Matrix);
	PSSM_DCT_i = PSSM_DCT_i(1:9,:);i
	PSSM_DCT_feature(i, :) = PSSM_DCT_i(:);	
end

% DWT compresses PSSM features of continuous residues and obtains PSSM-DWT features
% DWT
% PSSM_DWT_feature=[];
% for i=1:size(PSSM_with_window,1)
% 	PSSM_ii = [];
% 	PSSM_ii = PSSM_with_window(i,:);
%     PSSM_DWT_ii = [];
% 	PSSM_DWT_ii = EstraggoConDWT(PSSM_ii,4);
% 	PSSM_DWT_feature = [PSSM_DWT_feature;PSSM_DWT_ii];
% end

% PSA, load the PSA data sets in the document, the format is 'XXX.psa'
% PSA
PSSM_path = 'D:\Data\GTP\GTP_Test_PSA';
PSSM_file = '.psa';
PSA_path_file = [];
for i=1:Protein_num 
use_head_str = cell2mat(HEADER(i,1));
use_head_str(5)=[];
use_head_str(6)=[];
path1 = [PSSM_path '\' use_head_str PSSM_file];
Path_p = cellstr(path1);
PSA_path_file = [PSA_path_file; Path_p];
end

PSA=[];
for i=1:Protein_num
Sequence_i = SEQ(i);
protein_i=cell2mat(Sequence_i);
PSA_path_file2 = cell2mat(PSA_path_file(i));
PSA_Matrix = Read_PSA(PSA_path_file2, protein_i);
PSA=[PSA;PSA_Matrix];i
end


[Test_PSSM_DCT_feature, ia, ic] = unique(PSSM_DCT_feature, 'rows');
Test_label_y = label_test(ia, :);
Test_PSA_feature = PSA(ia, :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sample normalization process
train_sample = [Train_PSSM_DCT_feature,Train_PSA_feature];
test_sample = [Test_PSSM_DCT_feature,Test_PSA_feature];

sample = [train_sample;test_sample];
sample = line_map(sample);

Train_X = sample(1:size(train_sample,1),:);
Test_X = sample(size(train_sample,1)+1:end,:);
train_y = Train_label_y;
test_y = Test_label_y;


save('D:\Data\GTP\GTP_DCT_PSA.mat', 'Train_X', 'train_y', 'Test_X', 'test_y')
