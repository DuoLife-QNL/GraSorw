from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
import numpy as np
import pandas as pd
import math

# Train
nBlocks = 17
# dynamic block thresthold training
onDemandLoad_d = pd.read_csv('conf/on-demand-th-1.csv')
fullLoad_d = pd.read_csv('conf/on-demand-th-0.csv')

def computeRatioAndProcessTime(df):
    df['process-block-time'] = df['load-block-time'] + df['exec-block-time']
    df['ratio'] = df['nwalks'] / df['nverts']
    return df

for df in [onDemandLoad_d, fullLoad_d]:
    df = computeRatioAndProcessTime(df)

TRAIN_MINIMUM_DATA = 8
def keepData(df, columnName, compareValue, *, lessThan = True):
    if lessThan:
        return df[df[columnName] < compareValue]
    else:
        return df[df[columnName] > compareValue]

def trainModel(df):
    model = LinearRegression()
    xs = df[['ratio']].to_numpy()
    ys = df[['process-block-time']].to_numpy()
    model.fit(xs, ys)
    return model

def getFactor(model):
    return model.coef_[0, 0], model.intercept_[0]

def trainBlock(df_FL, df_OL):
    model_FL = trainModel(df_FL)
    coef_FL, intercept_FL = getFactor(model_FL)
    df_OL = keepData(df_OL, 'process-block-time', intercept_FL, lessThan=True)
    model_OL = trainModel(df_OL)
    coef_OL, intercept_OL = getFactor(model_OL)
    ths = (intercept_OL - intercept_FL) / (coef_FL - coef_OL)
    return ths

def trainBlocks(df_log_FL, df_log_OL, nblocks, col_blockId):
    df_log_FL = keepData(df_log_FL, 'ratio', 1, lessThan=True)
    ths = [None] * nblocks
    for blockId in range(0, nblocks):
        df_FL = df_log_FL[df_log_FL[col_blockId] == blockId]
        df_OL = df_log_OL[df_log_OL[col_blockId] == blockId]
        # 如果数据点小于TRAIN_MINIMUM_DATA个，我们认为训练是不准确的，这个block就当作full-load（阈值 = 0）
        if df_FL.shape[0] < TRAIN_MINIMUM_DATA or df_OL.shape[0] < TRAIN_MINIMUM_DATA:
            ths[blockId] = 0
            continue
        ths[blockId] = trainBlock(df_FL, df_OL)
    return  ths

# train dynamic block ths
ths_d = trainBlocks(fullLoad_d, onDemandLoad_d, nBlocks, 'dynamic-block-id')

def writeFiles(ths, fileName):
    f = open(fileName, 'w')
    for blockId in range(0, len(ths)):
        f.write("{}\n".format(ths[blockId]))
    f.close()

ths = 'conf/ths.txt'
writeFiles(ths_d, ths_kron)