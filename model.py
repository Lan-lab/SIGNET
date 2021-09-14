import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import torch
import torch.nn as nn
import torch.nn.functional as F
import cv2
from torch.utils.data import Dataset, DataLoader, TensorDataset
from torch.autograd import Variable
from torchvision import transforms
from torchvision import models
from PIL import Image
from sklearn.model_selection import train_test_split
import os
device = torch.device('cpu')


def sigmoid(x):
    return 1.0/(1 + np.exp(-x))


def dsigmoid(x):
    return np.multiply(sigmoid(x), (1-sigmoid(x)))


def accuracy(x,y):
    x = x.data.numpy()
    index = 0
    for i in range(y.shape[0]):
        if int(y[i]) == int(x[i]):
            index = index + 1
    print(index/y.shape[0])
    return index/y.shape[0]

def top10(x):
    min_index = []
    i = 0
    y = x[0:len(x)]
    while i<10:
        index = x.index(min(x))
        min_index.append(index)
        x[index] = 1
        i = i+1
    return min_index


class MLP(nn.Module):
    def __init__(self,len):
        super(MLP, self).__init__()
        num = 1
        n = int(len)
        while len>2.5:
          len = len / 2
          num = num * 2
        num = int(num)
        self.fc1 = nn.Linear(n,num)
        self.fc2 = nn.Linear(num,int(num/4))
        self.fc3 = nn.Linear(int(num / 4), 2)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = F.relu(self.fc3(x))
        return x

def l1_regularizer(model, lambda_l1=0.01):
    lossl1 = 0
    for model_param_name, model_param_value in model.named_parameters():
            if model_param_name.endswith('weight'):
                lossl1 += lambda_l1 * model_param_value.abs().sum()
    return lossl1


if __name__ == "__main__":
    data = sc.read_10x_mtx(
        root,  # 包含有 `.mtx` 文件的目录
        var_names='gene_symbols',  # 用 gene symbols 作为变量名 (variables-axis index)
        cache=True)  # 使用缓存文件加快读取
    # filter gene expressed in less 3 cell & cell expressed less gene than 200
    data = sc.read_csv("~/data_cut_fc6000.csv")
    data = data.transpose()
    sc.pp.filter_genes(data, min_cells=3)
    # mouse data: mt-, human data: MT-
    mito_genes = data.var_names.str.startswith('mt-')
    data.var['mito'] = mito_genes
    qc = sc.pp.calculate_qc_metrics(data, qc_vars=['mito'])
    cell_qc_dataframe = qc[0]
    gene_qc_dataframe = qc[1]
    data.obs["n_genes"] = cell_qc_dataframe['n_genes_by_counts']
    data.obs["n_counts"] = cell_qc_dataframe['total_counts']
    data.obs["percent_mito"] = cell_qc_dataframe['pct_counts_mito']
    sc.pl.violin(data, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True)
    sc.pl.scatter(data, x='n_counts', y='percent_mito')
    sc.pl.scatter(data, x='n_counts', y='n_genes')
    data = data[data.obs.n_genes < 9000, :]
    data = data[data.obs.percent_mito < 5, :]
    raw_data = pd.DataFrame(data.X.transpose().copy())
    col_name = data.obs_names
    row_name = data.var_names
    raw_data.index = row_name
    raw_data.columns = col_name
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)
    sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=6, min_disp=0.125)
    sc.pl.highly_variable_genes(data)
    gene_list = data.var.highly_variable
    gene_list = gene_list.index[gene_list == True]
    gene_list = gene_list.to_list()
    raw_data_fe = raw_data.loc[gene_list]
    tf_gene = pd.read_table('~/tf_mm.txt',header=None)
    list1 = tf_gene[0].to_list()
    list2 = raw_data_fe.index
    list5 = raw_data.index
    list3 = [x for x in list1 if x in list2]
    list4 = [x for x in list2 if x not in list1]
    raw_data_tf = raw_data.loc[list3, ]
    raw_data_ntf = raw_data_fe.loc[list4, ]
    gene_ntf = raw_data_ntf.index
    gene_tf = raw_data_tf.index
    HL_estimator = []
    thresholds = []
    data_ntf_binary = raw_data_ntf
    threshold_index = raw_data_ntf.shape[1]*(raw_data_ntf.shape[1]-1)/4
    for i in range(raw_data_ntf.shape[0]):
        gene_expr = raw_data_ntf.loc[gene_ntf[i], ]
        marker = [0 for x in range(int(gene_expr.max())+1)]
        record = []
        value = []
        for j in range(gene_expr.__len__()):
            marker[int(gene_expr[j])] += 1
        for j in range(int(gene_expr.max())+1):
            if marker[j] != 0:
                value.append(j)
                record.append(marker[j])
        value_calculation = []
        record_calculation = []
        for j in range(int(value.__len__())):
            for k in range(j, value.__len__()):
                value_calculation.append((value[j]+value[k])//2)
                if j == k:
                    record_calculation.append(record[j]*(record[j]-1)/2)
                else:
                    record_calculation.append(record[j]*record[k])
        test = {'index': value_calculation, 'number': record_calculation}
        HL_estimator = pd.DataFrame(test)
        HL_estimator_new = HL_estimator.sort_values(by=['index'],ascending=True)
        thresin = 0
        for j in range(record_calculation.__len__()):
            thresin += HL_estimator_new.iat[j, 1]
            if thresin >= threshold_index:
                threshold = HL_estimator_new.iat[j, 0]
                for k in range(data_ntf_binary.shape[1]):
                    if data_ntf_binary.iat[i, k] > threshold:
                        data_ntf_binary.iat[i, k] = 1
                    else:
                        data_ntf_binary.iat[i, k] = 0
                thresholds.append(HL_estimator_new.iat[j, 0])
                break

    HL_estimator = []
    thresholds = []
    data_tf_binary = raw_data_tf
    threshold_index = raw_data_tf.shape[1]*(raw_data_tf.shape[1]-1)/4
    for i in range(raw_data_tf.shape[0]):
        gene_expr = raw_data_tf.loc[gene_tf[i], ]
        marker = [0 for x in range(int(gene_expr.max())+1)]
        record = []
        value = []
        for j in range(gene_expr.__len__()):
            marker[int(gene_expr[j])] += 1
        for j in range(int(gene_expr.max())+1):
            if marker[j] != 0:
                value.append(j)
                record.append(marker[j])
        value_calculation = []
        record_calculation = []
        for j in range(int(value.__len__())):
            for k in range(j, value.__len__()):
                value_calculation.append((value[j]+value[k])//2)
                if j == k:
                    record_calculation.append(record[j]*(record[j]-1)/2)
                else:
                    record_calculation.append(record[j]*record[k])
        test = {'index': value_calculation, 'number': record_calculation}
        HL_estimator = pd.DataFrame(test)
        HL_estimator_new = HL_estimator.sort_values(by=['index'],ascending=True)
        thresin = 0
        for j in range(record_calculation.__len__()):
            thresin += HL_estimator_new.iat[j, 1]
            if thresin >= threshold_index:
                threshold = HL_estimator_new.iat[j, 0]
                for k in range(data_tf_binary.shape[1]):
                    if data_tf_binary.iat[i, k] > threshold:
                        data_tf_binary.iat[i, k] = 1
                    else:
                        data_tf_binary.iat[i, k] = 0
                thresholds.append(HL_estimator_new.iat[j, 0])
                break
    del HL_estimator
    del HL_estimator_new
    del record_calculation
    del value_calculation
    del list1
    del list2
    del list4
    data_ntf_binary_train = np.asarray(data_ntf_binary.transpose())
    data_tf_binary_train = np.asarray(data_tf_binary.transpose())
    np.savetxt("~/data_tf_binary_fc.txt", data_tf_binary_train)
    np.savetxt("~/data_ntf_binary_fc.txt", data_ntf_binary_train)

    gene_tf_d = pd.DataFrame(gene_tf)
    gene_ntf_d = pd.DataFrame(gene_ntf)
    gene_tf_d.to_csv("~/gene_tf.csv")
    gene_ntf_d.to_csv("~/gene_ntf.csv")
    """""
    data_ntf_binary_train = np.loadtxt("/Users/mac/Desktop/data_ntf_binary_verified2_brain.txt")
    data_tf_binary_train = np.loadtxt("/Users/mac/Desktop/data_tf_binary_verified2_brain.txt")
    """""
    coexpressed_result = np.zeros((data_ntf_binary_train.shape[1], data_tf_binary_train.shape[1] + 1))
    for j in range(data_ntf_binary_train.shape[1]):
        iterations = 15
        data_test_train = data_ntf_binary_train[:, [j]].copy()
        data_X = data_tf_binary_train.copy()
        data_Y = data_test_train.copy()
        index = 0
        bootstrapping = []
        for i in range(data_test_train.shape[0]):
            if int(data_test_train[i]) > 0:
                index = index + 1
        if index <= np.floor(data_test_train.shape[0] / 6):
            X = list(range(data_test_train.shape[0]))
            pos_flag = 0
            neg_flag = 0
            pos = np.floor(data_test_train.shape[0] / 6) + 1
            neg = data_test_train.shape[0] - pos
            while neg_flag < neg or pos_flag < pos:
                sample = int(np.floor(np.random.random() * len(X)))
                flag = data_test_train[sample]
                if flag > 0:
                    if pos_flag < pos:
                        bootstrapping.append(sample)
                        pos_flag = pos_flag + 1
                else:
                    if neg_flag < neg:
                        bootstrapping.append(sample)
                        neg_flag = neg_flag + 1
            for i in range(data_test_train.shape[0]):
                data_X[i, :] = data_tf_binary_train[bootstrapping[i], :].copy()
                data_Y[i, :] = data_test_train[bootstrapping[i], :].copy()
        x_train, x_test, y_train, y_test = train_test_split(data_X, data_Y, test_size=0.33, random_state=42)
        x_train = torch.from_numpy(x_train)
        x_test = torch.from_numpy(x_test)
        y_train = torch.from_numpy(y_train)
        y_test = torch.from_numpy(y_test)
        train_dataset = TensorDataset(x_train, y_train)
        test_dataset = TensorDataset(x_test, y_test)
        train_data = DataLoader(dataset=train_dataset, batch_size=64, shuffle=True)
        test_data = DataLoader(dataset=test_dataset, batch_size=64)
        flag = True
        times = 0
        while flag and times < 10:
            times = times + 1
            model = MLP(data_tf_binary_train.shape[1])
            criterion = nn.CrossEntropyLoss()
            optimizer = torch.optim.SGD(model.parameters(), 1e-2, momentum=0.9)
            losses = []
            acces = []
            eval_losses = []
            eval_acces = []
            for e in range(iterations):
                train_loss = 0
                train_acc = 0  # 开始训练
                model.train()
                for i, data in enumerate(train_data):
                    optimizer.zero_grad()
                    (inputs, labels) = data
                    labels = torch.tensor(labels, dtype=torch.long)
                    inputs = Variable(inputs)
                    labels = Variable(labels)
                    labels = labels.squeeze(1)
                    # 前向传播
                    labels = torch.tensor(labels, dtype=torch.long)
                    inputs = torch.tensor(inputs, dtype=torch.float32)
                    out = model(inputs)
                    loss = criterion(out, labels) + l1_regularizer(model, 1e-3)
                    # 反向传播
                    loss.backward()
                    optimizer.step()
                    # 记录误差
                    train_loss += loss.item()
                    # 计算分类的准确率
                    _, pred = out.max(1)
                    num_correct = (pred == labels).sum().item()
                    acc = num_correct / inputs.shape[0]
                    train_acc += acc

                losses.append(train_loss / len(train_data))
                acces.append(train_acc / len(train_data))
                # 在测试集上检验效果
                eval_loss = 0
                eval_acc = 0
                model.eval()  # 将模型改为预测模式
                for im, label in test_data:
                    im = Variable(im)
                    label = Variable(label.long())
                    label = torch.argmax(label, -1)
                    im = torch.tensor(im, dtype=torch.float32)
                    out = model(im)
                    loss = criterion(out, label)
                    # 记录误差
                    eval_loss += loss.item()
                    # 记录准确率
                    _, pred = out.max(1)
                    num_correct = (pred == label).sum().item()
                    acc = num_correct / im.shape[0]
                    eval_acc += acc

                eval_losses.append(eval_loss / len(test_data))
                eval_acces.append(eval_acc / len(test_data))
                print('epoch: {}, Train Loss: {:.6f}, Train Acc: {:.6f}, Eval Loss: {:.6f}, Eval Acc: {:.6f}'
                      .format(e, train_loss / len(train_data), train_acc / len(train_data),
                              eval_loss / len(test_data), eval_acc / len(test_data)))
            if max(eval_acces) > 0.5:
                flag = False

        if times >= 10 or flag == True:
            print("the gene:", j, "has accuracy lower than 50% with 10 times trying")
        accu = []
        for i in range(data_tf_binary_train.shape[1]):
            test = data_tf_binary_train.copy()
            test[:, i] = 0
            test = torch.from_numpy(test)
            test = torch.tensor(test, dtype=torch.float32)
            result = model(test)
            result = torch.argmax(result, dim=1)
            accu.append(accuracy(result, data_test_train))
        accu.index(min(accu))
        plt.plot(accu)
        plt.savefig("~/Figure{}.png".format(j))
        plt.show()
        accu.append(
            accuracy(
                torch.argmax(model(torch.tensor(torch.from_numpy(data_tf_binary_train), dtype=torch.float32)), dim=1),
                data_test_train))
        coexpressed_result[j, :] = accu

    np.savetxt("~/co_fc_ntf.txt", coexpressed_result)

    data_ntf_binary_train = data_tf_binary_train
    coexpressed_result = np.zeros((data_ntf_binary_train.shape[1], data_tf_binary_train.shape[1] + 1))
    for j in range(data_ntf_binary_train.shape[1]):
        iterations = 50
        data_test_train = data_ntf_binary_train[:, [j]].copy()
        data_X = data_tf_binary_train.copy()
        data_X[:,j] = 0
        data_medium = data_tf_binary_train.copy()
        data_medium[:, j] = 0
        data_Y = data_test_train.copy()
        index = 0
        bootstrapping = []
        for i in range(data_test_train.shape[0]):
            if int(data_test_train[i]) > 0:
                index = index + 1
        if index <= np.floor(data_test_train.shape[0] / 6):
            X = list(range(data_test_train.shape[0]))
            pos_flag = 0
            neg_flag = 0
            pos = np.floor(data_test_train.shape[0] / 6) + 1
            neg = data_test_train.shape[0] - pos
            while neg_flag < neg or pos_flag < pos:
                sample = int(np.floor(np.random.random() * len(X)))
                flag = data_test_train[sample]
                if flag > 0:
                    if pos_flag < pos:
                        bootstrapping.append(sample)
                        pos_flag = pos_flag + 1
                else:
                    if neg_flag < neg:
                        bootstrapping.append(sample)
                        neg_flag = neg_flag + 1
            for i in range(data_test_train.shape[0]):
                data_X[i, :] = data_medium[bootstrapping[i], :].copy()
                data_Y[i, :] = data_test_train[bootstrapping[i], :].copy()
        x_train, x_test, y_train, y_test = train_test_split(data_X, data_Y, test_size=0.33, random_state=42)
        x_train = torch.from_numpy(x_train)
        x_test = torch.from_numpy(x_test)
        y_train = torch.from_numpy(y_train)
        y_test = torch.from_numpy(y_test)
        train_dataset = TensorDataset(x_train, y_train)
        test_dataset = TensorDataset(x_test, y_test)
        train_data = DataLoader(dataset=train_dataset, batch_size=64, shuffle=True)
        test_data = DataLoader(dataset=test_dataset, batch_size=64)
        flag = True
        times = 0
        while flag and times < 10:
            times = times + 1
            model = MLP(data_tf_binary_train.shape[1])
            criterion = nn.CrossEntropyLoss()
            optimizer = torch.optim.SGD(model.parameters(), 1e-2, momentum=0.9)
            losses = []
            acces = []
            eval_losses = []
            eval_acces = []
            for e in range(iterations):
                train_loss = 0
                train_acc = 0  # 开始训练
                model.train()
                for i, data in enumerate(train_data):
                    optimizer.zero_grad()
                    (inputs, labels) = data
                    labels = torch.tensor(labels, dtype=torch.long)
                    inputs = Variable(inputs)
                    labels = Variable(labels)
                    labels = labels.squeeze(1)
                    # 前向传播
                    labels = torch.tensor(labels, dtype=torch.long)
                    inputs = torch.tensor(inputs, dtype=torch.float32)
                    out = model(inputs)
                    loss = criterion(out, labels) + l1_regularizer(model, 1e-3)
                    # 反向传播
                    loss.backward()
                    optimizer.step()
                    # 记录误差
                    train_loss += loss.item()
                    # 计算分类的准确率
                    _, pred = out.max(1)
                    num_correct = (pred == labels).sum().item()
                    acc = num_correct / inputs.shape[0]
                    train_acc += acc

                losses.append(train_loss / len(train_data))
                acces.append(train_acc / len(train_data))
                # 在测试集上检验效果
                eval_loss = 0
                eval_acc = 0
                model.eval()  # 将模型改为预测模式
                for im, label in test_data:
                    im = Variable(im)
                    label = Variable(label.long())
                    label = torch.argmax(label, -1)
                    im = torch.tensor(im, dtype=torch.float32)
                    out = model(im)
                    loss = criterion(out, label)
                    # 记录误差
                    eval_loss += loss.item()
                    # 记录准确率
                    _, pred = out.max(1)
                    num_correct = (pred == label).sum().item()
                    acc = num_correct / im.shape[0]
                    eval_acc += acc

                eval_losses.append(eval_loss / len(test_data))
                eval_acces.append(eval_acc / len(test_data))
                print('epoch: {}, Train Loss: {:.6f}, Train Acc: {:.6f}, Eval Loss: {:.6f}, Eval Acc: {:.6f}'
                      .format(e, train_loss / len(train_data), train_acc / len(train_data),
                              eval_loss / len(test_data), eval_acc / len(test_data)))
            if max(eval_acces) > 0.5:
                flag = False

        if times >= 10 or flag == True:
            print("the gene:", j, "has accuracy lower than 50% with 10 times trying")
        accu = []
        for i in range(data_tf_binary_train.shape[1]):
            test = data_medium.copy()
            test[:, i] = 0
            test = torch.from_numpy(test)
            test = torch.tensor(test, dtype=torch.float32)
            result = model(test)
            result = torch.argmax(result, dim=1)
            accu.append(accuracy(result, data_test_train))
        accu.index(min(accu))
        plt.plot(accu)
        plt.savefig("~/Figure{}.png".format(j+data_ntf_binary_train.shape[1]))
        plt.show()
        accu.append(
            accuracy(
                torch.argmax(model(torch.tensor(torch.from_numpy(data_tf_binary_train), dtype=torch.float32)), dim=1),
                data_test_train))
        coexpressed_result[j, :] = accu

   np.savetxt("~/co_fc_tf.txt", coexpressed_result)
