
import math
import numpy as np
from numpy import sort
import pandas as pd
import xgboost as xgb
import scipy.io as scio
import matplotlib.pyplot as plt
from sklearn import metrics
from imblearn.over_sampling import SMOTE
from imblearn.under_sampling import RandomUnderSampler
from collections import Counter


# load the data set in the 'XXX.mat' format after feature extraction
# the feature extraction part in SXGBsite is completed by MATLAB
data_path = "D:\Data\GTP\GTP_DCT_PSA.mat"
data_sample = scio.loadmat(data_path)

# load the dataset in the 'XXX.mat' format
# the training set 'train_X', the training set label 'train_y', the test set 'test_X', and the test set label 'test_y'
train_X = data_sample['Train_X']
train_y = np.ravel(data_sample['train_y'])
test_X = data_sample['Test_X']
test_y = np.ravel(data_sample['test_y'])
print(sorted(Counter(train_y).items()))

# the parameters of XGBoost in SXGBsite are as follows
params = {'booster': 'gbtree',
          'objective': 'binary:logistic',
          'eval_metric': 'auc',
          'learning_rate': 0.1,
          'min_child_weight': 1,
          'max_depth': 9,
          'gamma': 0.05,
          'lambda': 10,
          'silent': 1}

# SMOTE over-sampling process, where '19000' is the size of positive (or negative) samples
print(" SMOTE begin...")
print(" ...")
SMOTE_params = SMOTE(ratio={1: 19000}, random_state = 0)
train_X_SMOTED, train_y_SMOTED = SMOTE_params.fit_sample(train_X, train_y)
rus = RandomUnderSampler(ratio={0: 19000}, random_state = 0)
train_X_SMOTE, train_y_SMOTE = rus.fit_sample(train_X_SMOTED, train_y_SMOTED)
print(sorted(Counter(train_y_SMOTE).items()))
print(" SMOTE end.")

# build the prediction model by XGBoost
print(" Training Begin...")
dtrain = xgb.DMatrix(train_X_SMOTE, label=train_y_SMOTE)
dtest = xgb.DMatrix(test_X)
watchlist = [(dtrain, 'train')]
bst = xgb.train(params, dtrain, num_boost_round=200, evals = watchlist)
print(" Training End.")

# output the probability value label of the prediction results (range between 0 and 1) and the AUC of the prediction results
print(" Testing Begin...")
ypred = bst.predict(dtest)
thres = 0.5
y_pred = (ypred >= thres) * 1
print(" Testing End.")
print(' AUC: %.4f' % metrics.roc_auc_score(test_y,ypred))

# the function of Threshold Moving process, to get the best MCC
def cal_rate(pred_y, test_y, thres_min, thres_step, thres_max):
    mcc_list = []
    acc_list = []
    sn_list = []
    sp_list = []
    X = []

    p = np.copy(pred_y)
    P_number = len(p)

    print(" Threshold_find begin...")
    for k in range(thres_min, thres_max, thres_step):
        Threshold_t = k / 1000
        Predict_label = pred_y
        TP = 0
        FP = 0
        FN = 0
        TN = 0
        for i in range(P_number):
            PD = p[i]
            if PD >= Threshold_t:
                PD = 1
                Predict_label[i] = 1
            if PD == 1:
                if test_y[i] == 1:
                    TP += 1
                else:
                    FP += 1
            else:
                Predict_label[i] = 0
                if test_y[i] == 0:
                    TN += 1
                else:
                    FN += 1
        # print(TP+FP+TN+FN)
        accracy = float(TP + TN) / float(TP + FP + TN + FN)
        if TP + FP == 0:
            precision = 0
        else:
            precision = float(TP) / float(TP + FP)
        TPR = float(TP) / float(TP + FN)
        TNR = float(TN) / float(FP + TN)
        FNR = float(FN) / float(TP + FN)
        FPR = float(FP) / float(FP + TN)
        MCC = float(TP * TN - FP * FN) / float(math.sqrt((TP + FN) * (TN + FP) * (TP + FP) * (TN + FN)))
        Recall = TPR

        if Threshold_t == 0.500:
            print(" When the Threshold is 0.5 , the value of AUC , ACC , SN , Spec , Precesion , F1-score , Threshold is ...")
            print(' MCC: %.4f' % MCC)
            print(' ACC: %.4f' % accracy)
            print(' SN: %.4f' % TPR)
            print(' Spec: %.4f' % TNR)
            print(' Threshold: %.4f' % Threshold_t)

        mcc_list.append(MCC)
        acc_list.append(accracy)
        sn_list.append(TPR)
        sp_list.append(TNR)
        X.append(Threshold_t)

    print(" Threshold_find end.")
    print(" When the value of MCC is max , the value of AUC , ACC , SN , Spec , Precesion , F1-score , Threshold is ...")
    print(' max_MCC: %.4f' % max(mcc_list))
    maxmcc = mcc_list.index(max(mcc_list))
    print(' max_MCC_ACC: %.4f' % acc_list[maxmcc])
    print(' max_MCC_SN: %.4f' % sn_list[maxmcc])
    print(' max_MCC_Spec: %.4f' % sp_list[maxmcc])
    print(' max_MCC_Threshold: %.4f' % X[maxmcc])
    print(" Plot threshold_find.")

    # output SN, SP, ACC and MCC curves with Thresholds
    fig = plt.figure()
    plt.ylim(0, 1)
    plt.plot(X, mcc_list, c='b', label='MCC')
    plt.plot(X, acc_list, c='r', label='ACC')
    plt.plot(X, sn_list, c='y', label='SN')
    plt.plot(X, sp_list, c='g', label='Spec')
    plt.xlabel('Threshold')
    plt.ylabel('Values')
    plt.legend()
    plt.show()

    return


if __name__ == "__main__":
    # the 100 and 900 set here represent the Threshold Moving range of 100/1000 - 900/1000 (0.100 - 0.900) in steps of 1/1000 (0.001)
    print(cal_rate(ypred, test_y, 100, 1, 900))
