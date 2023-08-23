'''
Introduction:
Simple Python Code for the Paper :
RaPlace：Place Recognition for Imaging Radar using Radon Transform and Mutable Threshold
'''

# Create Time: Aug 11th, 2023

import numpy as np
import math
import os
import cv2
import matplotlib.pyplot as plt
from skimage.transform import radon
from scipy import ndimage
import scipy.io as sio  # when the python file deal with .mat document

# Calculate the Euclidean distance between two two-dimensional postures
def dist_btn_pose(pose1, pose2):
    dist = math.sqrt((pose1[0] - pose2[0])**2 + (pose1[1] - pose2[1])**2)
    return dist

# Radar Similarity Appraisal using Discrete Fourier Transform
# resulting in a one-dimensional array of cross-correlation
def fast_dft(Mq, Mi):
    Fq = np.fft.fft(Mq, axis=0)  # fft along theta axis
    Fn = np.fft.fft(Mi, axis=0)
    corrmap_2d = np.fft.ifft(Fq * np.conj(Fn), axis=0)

    corrmap = np.sum(corrmap_2d, axis=-1)
    maxval = np.max(corrmap)
    return corrmap, maxval

# Another Radon Transform Method
# steps is a constant, usually steps = 180
def DiscreteRadonTransform(image, steps):
    channels = len(image[0])
    res = np.zeros((channels, steps), dtype='float64')
    for s in range(steps):
        rotation = ndimage.rotate(image, -s*180/steps, reshape=False).astype('float64')
        #print(sum(rotation).shape)
        res[:,s] = sum(rotation)
    return res

# Extract the middle row of the matrix
def rowkey(sino):
    row, col = sino.shape
    row_key = sino[(row + 1) // 2 - 1, :]
    # row_key = sino[math.floor((row + 1) / 2), :]
    return row_key

# Returns a list of file names under the specified path, except for special directories
def osdir(path):
    files = os.listdir(path)
    files = [file for file in files if not file.startswith('.')]

    return files

# Perform radon transformation
def generateRadon(data_dir, down_shape):
    # radar_data_dir = os.path.join(data_dir, 'sensor_data/radar/backward_cart/')

    radar_data_dir = 'D:/polar/polar'
    data_names = osdir(radar_data_dir)

    # gtpose = np.loadtxt(os.path.join(data_dir, 'global_pose.csv'), delimiter=',')

    gtpose = np.loadtxt('D:/polar/global_pose.csv', delimiter=',')
    # CSV format file storing datas are separated by commas
    
    gtpose_time = gtpose[:, 0]
    gtpose_xy = gtpose[:, [4, 8]]

    num_data = len(data_names)
    theta = np.arange(0, 180)
    sinoffts = []
    rowkeys = np.zeros((num_data, int(down_shape * 180)))
    xy_poses = np.zeros((num_data, 2))

    for data_idx in range(num_data):
        file_name = data_names[data_idx]
        data_time = float(file_name[:-4])
        data_path = os.path.join(radar_data_dir, file_name)

        tmp = cv2.imread(data_path, cv2.IMREAD_GRAYSCALE)
        R = radon(tmp, theta)
        xp = np.arange(-R.shape[1] // 2, R.shape[1] // 2)

        R = R / np.max(R)
        R = R.astype(np.float64)
        R = cv2.resize(R, (int(down_shape * R.shape[1]), int(down_shape * R.shape[0])))

        # sinofft = np.abs(np.fft.fft(R, axis=0)[:R.shape[0] // 2, :])
        sinofft = np.abs(np.fft.fft(R, axis=0))
        sinofft_rows = sinofft.shape[0]
        sinofft = sinofft[:sinofft_rows // 2, :]

        # nearest_idx = np.argmin(np.abs(data_time - gtpose_time))

        time_diff = np.abs(np.tile(data_time, len(gtpose_time)) - gtpose_time)
        nearest_time_gap = np.min(time_diff)
        nearest_idx = np.argmin(time_diff)

        xy_pose = gtpose_xy[nearest_idx]

        sinoffts.append(sinofft)
        rowkeys[data_idx, :] = rowkey(R).flatten()
        xy_poses[data_idx, :] = xy_pose

        if data_idx % 100 == 0:
            print(f"{data_idx} / {num_data} processed.")

    return sinoffts, rowkeys, xy_poses


data_path = 'D:/polar/polar'
print(f"Processing: {data_path[-12:-10]}{data_path[-4:-1]} within criteria m/")
down_shape = 0.1

data_sinofft, data_rowkeys, data_poses = generateRadon(data_path, down_shape)

revisit_criteria = 20
keyframe_gap = 1

num_candidates = 10
num_node_enough_apart = 50
num_top_n = 25

top_n = np.linspace(1, num_top_n, num_top_n)

# Entropy thresholds
middle_thres = 0.001
thresholds1 = np.linspace(0, middle_thres, 50)
thresholds2 = np.linspace(middle_thres, 0.01, 50)
thresholds = np.concatenate((thresholds1, thresholds2))
num_thresholds = len(thresholds)

# Main variables to store the result for drawing PR curve
num_hits = np.zeros((num_top_n, num_thresholds))
num_false_alarms = np.zeros((num_top_n, num_thresholds))
num_correct_rejections = np.zeros((num_top_n, num_thresholds))
num_misses = np.zeros((num_top_n, num_thresholds))

# Check the threshold for max value
case_for_hit = np.zeros((len(data_poses), 3))
case_for_fa = []

loop_log = []
exp_poses = []
exp_rowkeys = []
exp_sinofft = []
exp_sinograms = []

num_queries = len(data_poses)

for query_idx in range(num_queries - 1):
    query_sinofft = data_sinofft[query_idx]
    query_sinofft = (query_sinofft - np.mean(query_sinofft)) / np.std(query_sinofft)
    query_rowkey = data_rowkeys[query_idx, :]
    query_pose = data_poses[query_idx, :]

    exp_sinofft.append(query_sinofft)
    exp_rowkeys.append(query_rowkey)
    exp_poses.append(query_pose)

    if query_idx % keyframe_gap != 0:
        continue

    if query_idx % keyframe_gap != 0:
        continue

    can_sinofft = exp_sinofft[:-(num_node_enough_apart - 1)]

    tmpval = 0
    maxval = 0
    rotval = 0
    candnum = 0

    for cands in range(len(can_sinofft)):
        tmp_sinofft = can_sinofft[cands]
        fftresult, tmpval = fast_dft(query_sinofft, tmp_sinofft)
        if maxval < tmpval:
            maxval = tmpval
            candnum = cands

    nearest_idx = candnum
    fftresult, tmpval = fast_dft(query_sinofft, query_sinofft)
    min_dist = (tmpval - maxval) / 1000
    real_dist = dist_btn_pose(query_pose, exp_poses[nearest_idx])

    for topk in range(num_top_n):
        for thres_idx in range(num_thresholds):
            threshold = thresholds[thres_idx]
            reject = 0

            if min_dist > threshold:
                reject = 1

            if reject == 1:
                if real_dist < revisit_criteria:
                    num_correct_rejections[topk, thres_idx] += 1
                else:
                    num_misses[topk, thres_idx] += 1
            else:
                if real_dist < revisit_criteria:
                    num_hits[topk, thres_idx] += 1
                    case_for_hit[query_idx, 0] = 1
                    case_for_hit[query_idx, 1] = nearest_idx
                    case_for_hit[query_idx, 2] = dist_btn_pose(query_pose, exp_poses[nearest_idx])
                else:
                    num_false_alarms[topk, thres_idx] += 1

    if query_idx % 100 == 0:
        print(f"{query_idx/num_queries * 100} % processed")

# PR-result

ResultsDir = './pr_result/'
title_str = 'MulRan Sequence (radar polar)'
FigIdx = 2

# Create figure
plt.figure(FigIdx)
plt.clf()

TopNindexes = [0]  # Python indexing starts from 0
name = 'top1'

nTopNindexes = len(TopNindexes)

# Main
# SequenceNames = os.listdir(ResultsDir)
# SequenceNames = [name for name in SequenceNames if os.path.isdir(os.path.join(ResultsDir, name))]
# nSequences = len(SequenceNames)

nSequences = num_top_n

all_Precisions = []
all_Recalls = []

for ithTopN in range(nTopNindexes):
    TopNidx = TopNindexes[ithTopN]
    line_width = 4

    LineColors = plt.cm.rainbow(np.linspace(0, 1, nSequences))
    LineColors = np.flipud(LineColors)

    AUCs = np.zeros(nSequences)

    for ithSeq in range(nSequences):
        nCorrectRejectionsAll = num_correct_rejections
        nCorrectRejectionsForThisTopN = nCorrectRejectionsAll[TopNidx, :]

        nFalseAlarmsAll = num_false_alarms
        nFalseAlarmsForThisTopN = nFalseAlarmsAll[TopNidx, :]

        nHitsAll = num_hits
        nHitsForThisTopN = nHitsAll[TopNidx, :]

        nMissesAll = num_misses
        nMissesForThisTopN = nMissesAll[TopNidx, :]

        nThres = nCorrectRejectionsAll.shape[1]

        Precisions = np.empty(nThres)
        Recalls = np.empty(nThres)

        for ithThres in range(nThres):
            nCorrectRejections = nCorrectRejectionsForThisTopN[ithThres]
            nFalseAlarms = nFalseAlarmsForThisTopN[ithThres]
            nHits = nHitsForThisTopN[ithThres]
            nMisses = nMissesForThisTopN[ithThres]

            nTotalTestPlaces = nCorrectRejections + nFalseAlarms + nHits + nMisses

            Precision = nHits / (nHits + nFalseAlarms)
            Recall = nHits / (nHits + nMisses)

            Precisions[ithThres] = Precision
            Recalls[ithThres] = Recall

        all_Precisions.append(Precisions)
        all_Recalls.append(Recalls)

        # Draw
        # plot all the figures

        # plt.figure()
        plt.figure(figsize=(10, 6))
        plt.plot(Recalls, Precisions, linewidth=line_width, color=LineColors[ithSeq], label = 'MulRan Radar Dataset')
        plt.title(title_str, fontsize=12)
        plt.xlabel('Recall', fontsize=12)
        plt.ylabel('Precision', fontsize=12)
        plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.grid(True)
        plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
        # plt.legend(SequenceNames, loc='best', fontsize=9)
        # plt.legend(str(nSequences), loc='best', fontsize=9)
        plt.legend(loc='best')
        name = 'prcurve'
        # plt.savefig(name + '.pdf', format='pdf', bbox_inches='tight')

        plt.show()

# plot the last figure

# plt.figure(FigIdx)
# plt.figure(figsize=(10, 6))
# for ithSeq in range(nSequences):
#     plt.plot(all_Recalls[ithSeq], all_Precisions[ithSeq], linewidth=line_width, color=LineColors[ithSeq])
# plt.title(title_str, fontsize=12)
# plt.xlabel('Recall', fontsize=12)
# plt.ylabel('Precision', fontsize=12)
# plt.xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
# plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
# plt.xlim([0, 1])
# plt.ylim([0, 1])
# plt.grid(True)
# plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
# name = 'prcurve'
# plt.show()

