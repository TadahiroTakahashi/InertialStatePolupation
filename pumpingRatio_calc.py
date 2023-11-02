import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import time
#not using
#from numba import jit

try:
    PATH = sys.argv[1] 
except IndexError:
    PATH = "./dataSet/0_test/"
    print("no argv[1],please set argv[1] as PATH")

#functions
def rk4(f, y0, t, h,W):
    """
    4次のRunge-Kutta法を使用して微分方程式を解く関数。
    f: 微分方程式の関数
    y0: 初期値
    t: 時間の配列
    h: 時間ステップ
    """
    y = np.zeros((len(t), len(y0)),dtype="float")
    y[0] = y0
    for i in range(0, len(t)-1):
        k1 = h * f(y[i], t[i],W)
        k2 = h * f(y[i] + 0.5*k1, t[i] + 0.5*h,W)
        k3 = h * f(y[i] + 0.5*k2, t[i] + 0.5*h,W)
        k4 = h * f(y[i] + k3, t[i] + h,W)
        y[i+1] = y[i] + (k1 + 2*k2 + 2*k3 + k4) / 6.0
    return y

def Nlevel_equations_for992nm(y, t, W,init=4):
    #return np.dot(W - np.diag(np.diag(W)), y) - np.dot(np.diag(np.sum(W, axis=1)), y)
    ans = np.dot(W, y) 
    ans[init]+= ϕ  
    ans[ans <= 0] = 0
    return ans


def laserExcite_add_to_W(W,s,δ,l,u):
    γ = matrix[l,u]
    Ω = np.sqrt(s*γ**2/2)
    γall = np.sum(matrix[:,u])-matrix[u,u]
    R = Ω**2/γall/(1+(2*δ/γall)**2)#2021Akatsuka(Metcalf) 
    W[l,l] -= R
    W[l,u] += R
    W[u,u] -= R
    W[u,l] += R
    return W

def atomLoss_add_to_W(W,outside):
    W[0,1:] = outside
    for i in range(W.shape[0]):
        if i != 0: 
            W[i,i] -= outside
        else:
            W[0,0] = 0
    return W

#install Grotorian chart
N = 19
try:
    input_book = pd.ExcelFile('./grotorian.xlsx')
except IOError:
    print("please set grotorian.xlsx in %s"%PATH)

input_sheet_name = input_book.sheet_names
#sheetName = "easiest"
sheetName = "referred"
matrix = np.array(input_book.parse(sheetName,header=2))
matrix = matrix[:N,3:3+N]
stateName = np.array(input_book.parse(sheetName))
colorName = stateName[3+N,3:]
stateName = stateName[0,3:]

#install paramSet
try:
    paramSet = pd.ExcelFile(PATH+'/paramSet.xlsx')
Except IOError:
    print("please set paramSet.xlsx in %s"%PATH)
input_sheet_name = input_book.sheet_names
sheetName = "params"

#install condSet
condDataRaw = paramSet.parse(sheetName,header=0)[:4].to_dict()
condDict = {condDataRaw['Unnamed: 0'][i]: condDataRaw['Unnamed: 2'][i] for i in condDataRaw['Unnamed: 0'].keys()}
print(condDict)
ϕ = condDict['flux']
h = condDict['calcStepTime']
step = condDict['NumOfStep']
outside = condDict['outside']

#install laser conditions
paramData = paramSet.parse(sheetName,header=6)
paramDictRaw = paramData.to_dict()
paramDict = {}
for key, value in paramDictRaw.items():
    if key != "WM":
        paramDict[key] = {paramDictRaw['WM'][k]: v for k, v in value.items()}
print(paramDict)
    

# maing matrix 
W = matrix*1.0
y0 = np.zeros(N)

for key,val in paramDict.items():
    W = laserExcite_add_to_W(W,s=val["int"],δ=val["detu"],l=val["lowerS"],u=val["upperS"])
W = atomLoss_add_to_W(W,outside)

print("calc start")

def hTuning():
    RKValid = False
    htest = h*1.0
    while RKValid == False:
        stepTest = 1000
        tTest = np.arange(0, htest*stepTest, htest)
        resultTest = rk4(Nlevel_equations_for992nm, y0, tTest, htest, W)
        check = np.sum(resultTest.T,axis=0)[-1]
        check_True = htest*stepTest*ϕ
        if abs(check - check_True) < 0.01:
            RKValid = True
            print("RK4 is valid")
        else:
            htest = htest*0.2
            print("RK4 is invalid, htest: %.2e"%htest)
    return htest

h = hTuning()
t = np.arange(0, h*step, h)
t0 = time.time()


result = rk4(Nlevel_equations_for992nm, y0, t, h, W)
#result = rk4_optimized(Nlevel_equations_for992nm_optimized, y0, t, h, W, tLim, ϕ)
t1 = time.time()
print("time: %.2f s"%(t1-t0))
result = result.T



MATPLOTLIB= True
#matplitlib
if MATPLOTLIB:
    plt.rcParams['lines.markeredgecolor'] = "black"
    plt.rcParams['lines.markerfacecolor'] = "white"
    plt.rcParams['lines.markersize'] = 10
    plt.rcParams['font.size'] = 14


    fig = plt.figure(figsize=(8,8))

    plt.subplots_adjust(left=0.08, right=0.98, bottom=0.15, top=0.95, wspace=0.45, hspace=0.3)
    axTime = plt.subplot2grid((2,1), (0,0))
    axTimeL = plt.subplot2grid((2,1), (1,0))
    #axLegend = plt.subplot2grid((1,2), (0,1),rowspan=1)
    #plot summation: should be linear-increasing function.
    axTime.plot(t*1000,np.sum(result,axis=0),"--",lw=2,c="black",label="total")
    axTime.plot(t*1000,np.sum(result[1:,:],axis=0),"--",lw=2,c="gray",label="alived atoms")
    for ax in [axTime,axTimeL]:#axTimeL is used just for enlarge near the x axis. 
        for (num,result_each) in enumerate(result):
            if num in [0,1,2,3,4,8,9,10,14,18]:
                ax.plot(t*1000, result_each,label=stateName[num],color=colorName[num],lw=1)
            elif num in [5,6,7]:
                ax.plot(t*1000, result_each,label=stateName[num],color=colorName[num],ls="dashed",lw=0.8)
        ax.set_xlim(-h*step*0.3*1000,h*step*1.1*1000)
        ax.set_xlabel('Time(ms)')
        ax.set_ylabel('Density')
        ax.legend(loc = "upper left",fontsize=10)
        ax.set_title('atom population')
        ax.grid(True)
    axTime.set_ylim(-0.05,1.05*ϕ*(h*step))
    axTimeL.set_ylim(-0.002,0.01*ϕ*(h*step))   
    plt.savefig(PATH+"/result.pdf")

    plt.show()


