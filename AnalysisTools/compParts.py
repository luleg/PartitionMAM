import numpy as np
import os
from sklearn.metrics.cluster import contingency_matrix
from sklearn.metrics import ConfusionMatrixDisplay
import munkres as mk
import matplotlib.pyplot as plt
import argparse



def visualComp(assign1,assign2,saveFigs=False,showFigs=True,visu='both',path2Fig = '',figType='png'):

    # 1 - read the partition files and translate them into a contingency matrix

    n = len(assign1)
    assert n==len(assign2),"The two partitions must have been built on the same set"
    # Apply the HUngarian algorithm to find the "best" correspondance between the labels from the two classes
    munk = mk.Munkres()
    CostMat = contingency_matrix(assign1,assign2)
    CostMat = mk.make_cost_matrix(CostMat,lambda x:n-x)
    CostMat = munk.pad_matrix(CostMat,pad_value = n)
    PairsLabs = munk.compute(CostMat)
    # Sort labels, so that classes with highest number of joint elements are first
    nbCl = len(CostMat)
    CardClass = [0 for i in range(nbCl)]
    for cpt,(ind_r,ind_c) in enumerate(PairsLabs):
        CardClass[cpt] = (cpt,n-CostMat[ind_r][ind_c])

    CardClass.sort(reverse=True,key=lambda x:x[1])
    PairsLabs = [PairsLabs[i] for i,_ in CardClass]

    # dico of initial labels/new labels
    DicoAssign1 = {PairsLabs[u][0]:u for u in range(nbCl)}
    DicoAssign2 = {PairsLabs[u][1]:u for u in range(nbCl)}

    assign1 = [DicoAssign1[u] for u in assign1]
    assign2 = [DicoAssign2[u] for u in assign2]

    DicoAssign1 = list(set(assign1))
    DicoAssign1.sort()
    DicoAssign1 = {u:i for i,u in enumerate(DicoAssign1)}
    assign1 = [DicoAssign1[u] for u in assign1]

    DicoAssign2 = list(set(assign2))
    DicoAssign2.sort()
    DicoAssign2 = {u:i for i,u in enumerate(DicoAssign2)}
    assign2 = [DicoAssign2[u] for u in assign2]


    if (visu!='conf'):
        # To plot
        BreakPts2 = [-0.5]
        val2 = []
        BreakPts1 = [-0.5]
        Cls1 = list(set(assign1))
        Cls1.sort()
        for vv,oneClass in enumerate(Cls1):
            cc = [i for i,u in enumerate(assign1) if u == oneClass]
            BreakPts1.append(BreakPts1[-1]+len(cc))
            eltsInCl1 = [assign2[u] for u in cc]
            OtherClInCl1 = list(set(eltsInCl1))
            OtherClInCl1.sort()
            for anotherClass in OtherClInCl1:
                val2.append([anotherClass,vv])
                BreakPts2.append(BreakPts2[-1]+eltsInCl1.count(anotherClass))

        # BreakPts2 = BreakPts2[1:-1]
        BreakPts1 = BreakPts1[1:-1]

        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

        for cpt in range(len(BreakPts2)-1):
            if val2[cpt][1]%2:
                ax.plot([BreakPts2[cpt],BreakPts2[cpt+1]],[val2[cpt][0],val2[cpt][0]],'black')
            else:
                ax.plot([BreakPts2[cpt],BreakPts2[cpt+1]],[val2[cpt][0],val2[cpt][0]],'darkcyan')

        for bb in BreakPts1[:-1]:
            ax.plot([bb,bb],[-0.5,nbCl],'red',linewidth=0.5,linestyle='dashed')

        ax.plot([BreakPts1[-1],BreakPts1[-1]],[-0.5,nbCl],'red',linewidth=0.5,linestyle='dashed',label='Separators between parts from the 1st partition')


        ax.axis([0, n-0.5, -0.5, max(assign2)+0.5])
        ax.set_ylabel('New cluster labels of the 2nd partition')
        ax.legend()
        if (saveFigs):
            fig.set_size_inches(6, 4)
            fig.savefig(path2Fig+'.part.'+figType)

        if showFigs:
            fig.set_size_inches(6, 4)
            plt.show()




    if (visu!='part'):
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        ConfusionMatrixDisplay.from_predictions(assign1,assign2,normalize=None,include_values=True,ax=ax,xticks_rotation='vertical')
        ax.set_ylabel('1st Partition new labels')
        ax.set_xlabel('2nd Partition new labels')
        ax.set_ylim([max(assign1)+0.5,-0.5])
        ax.set_xlim([-0.5,max(assign2)+0.5])
        if (saveFigs):
            fig.set_size_inches(6, 4)
            fig.savefig(path2Fig+'.conf.'+figType)

        if showFigs:
            fig.set_size_inches(6, 4)
            plt.show()

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Relabel two partitions to ease their comparison, and display a visual comparison.')
    parser.add_argument('File2Part1', metavar='FP1.part',
                        help='Path to the file containing the 1st partition.')
    parser.add_argument('File2Part2', metavar='FP2.part',
                        help='Path to the file containing the 2nd partition.')
    parser.add_argument('-outFig', metavar='FoldToFig', default ='',
                        help='Path to the file to save the figures that compares the partitions. By default, figures are shown but not saved.')

    parser.add_argument('-ext', metavar='extension', default ='png',
                        help='Extension for the figures, if saved (default is png).')
    args = parser.parse_args()

    if args.outFig:
        saveFigs = True
    else:
        saveFigs = False


    tmp = np.loadtxt(open(args.File2Part1,'r'),dtype=int,delimiter=' ')
    n = tmp.shape[0]
    assign1 =[tmp[u,1] for u in range(n)]

    tmp = np.loadtxt(open(args.File2Part2,'r'),dtype=int,delimiter=' ')
    ntmp = tmp.shape[0]
    assert ntmp==n,"The two partitions must have been built on the same set"
    assign2 =[tmp[u,1] for u in range(n)]

    visualComp(assign1,assign2,saveFigs=saveFigs,showFigs=True,visu='both',path2Fig = args.outFig,figType=args.ext)
