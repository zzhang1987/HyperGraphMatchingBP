import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import cPickle as pickle
import pylab
Fname = 'MotorRes.pkl'
NofInstance = 20
f = open(Fname, "rb")
AllAcc = pickle.load(f)
AllRtime = pickle.load(f)
AllObj = pickle.load(f)
f.close()

SecondOrderMethods = ('Ours', 'FGM')
ThirdOrderMethods = ('Ours', 'OursBca', 'BCA',
                     'BCA-MP', 'BCA-IPFP', 'HGM',
                     'RRWHM', 'TM')

Type = 'pas'

WithEdge = False
NofOus = len(AllAcc) - 1

WithTriplet = True
MeanAcc = dict()


Colors = [[0., 0.5, 0.5, 1.],
          [1., 0.11692084, 0., 1.],
          [1, 0.00196078, 1., 1.],
          [0., 0.50392157, 1., 1.],
          [0.08538899, 1., 0.88235294, 1.],
          [0.49019608, 1., 0.47754586, 1.],
          [0.5, 1., 0.17273877, 1.],
          [1., 0.58169935, 0., 1.],
          [1., 0., 1., 1.],
          [0.1, 0.9, 0.1, 1]]
Markers = ['o', 'v', 'o', 'd', 'x', 'v', 'h', 'D', 'H', 'v']
LineStype = ['-', '--', '-.', '-', '-', '--', '-', '--', '-.', '-']

figData = pylab.figure(num=None, figsize=(4, 3), dpi=80,
                       facecolor='w', edgecolor='k')
ax = pylab.gca()
matplotlib.rc('font', family='Times New Roman')
cnt = 0
for methods in ThirdOrderMethods:
    FullMethods = (methods + 'WithEdge' + str(WithEdge) +
                   'WithTriplet' + str(True) + Type)
    CMeanAcc = np.zeros(NofOus + 1)
    for n in range(NofOus+1):
        CMeanAcc[n] = np.mean(list(AllAcc[n][FullMethods].values()))
    pylab.plot(range(0, NofOus + 1), CMeanAcc, label=methods,
               color=Colors[cnt], marker=Markers[cnt],
               linestyle=LineStype[cnt])
    cnt = cnt + 1

pylab.grid(True)

pylab.ylabel('Accuracy')
pylab.xlabel('#Outliers')
    
plt.figure(num=None, figsize=(4, 3), dpi=80, facecolor='w', edgecolor='k')
cnt = 0
MaxObj = np.zeros([NofOus + 1, NofInstance])
for n in range(NofOus + 1):
    cMaxObj = np.zeros(NofInstance)
    for methods in ThirdOrderMethods:
        FullMethods = (methods + 'WithEdge' + str(WithEdge) +
                       'WithTriplet' + str(True) + Type)
        CObjs = np.array(list(AllObj[n][FullMethods].values()))
        cMaxObj = np.array(np.max(np.array([cMaxObj.tolist(), CObjs.tolist()]),
                                  axis=0))
    MaxObj[n] = cMaxObj
for methods in ThirdOrderMethods:
    FullMethods = (methods + 'WithEdge' + str(WithEdge) +
                   'WithTriplet' + str(True) + Type)
    CMeanObj = np.zeros(NofOus + 1)
    for n in range(NofOus+1):
        CObj = np.array(list(AllObj[n][FullMethods].values()))
        CObj = CObj / MaxObj[n]
        CMeanObj[n] = CObj.mean()
    plt.plot(range(0, NofOus + 1), CMeanObj, label=methods,
             color=Colors[cnt], marker=Markers[cnt],
             linestyle=LineStype[cnt])
    cnt = cnt + 1
pylab.ylim([0.2, 1.05])

plt.ylabel('Obj. Value')
plt.xlabel('#Outliers')
plt.grid(True)

WithEdge = True

figData = pylab.figure(num=None, figsize=(4, 3), dpi=80,
                       facecolor='w', edgecolor='k')
ax = pylab.gca()
matplotlib.rc('font', family='Times New Roman')
cnt = 0
for methods in ThirdOrderMethods:
    FullMethods = (methods + 'WithEdge' + str(WithEdge) +
                   'WithTriplet' + str(True) + Type)
    CMeanAcc = np.zeros(NofOus + 1)
    for n in range(NofOus+1):
        CMeanAcc[n] = np.mean(list(AllAcc[n][FullMethods].values()))
    print((methods, CMeanAcc))

    pylab.plot(range(0, NofOus + 1), CMeanAcc, label=methods,
               color=Colors[cnt], marker=Markers[cnt],
               linestyle=LineStype[cnt])
    cnt = cnt + 1

for methods in SecondOrderMethods:
    FullMethods = (methods + 'WithEdge' + str(WithEdge) +
                   'WithTriplet' + str(False) + Type)
    CMeanAcc = np.zeros(NofOus + 1)
    for n in range(NofOus+1):
        CMeanAcc[n] = np.mean(list(AllAcc[n][FullMethods].values()))
    print((methods, CMeanAcc))
    
    pylab.plot(range(0, NofOus + 1), CMeanAcc, label=methods,
               color=Colors[cnt], marker=Markers[cnt],
               linestyle=LineStype[cnt])
    cnt = cnt + 1
    
pylab.grid(True)

pylab.ylabel('Accuracy')
pylab.xlabel('#Outliers')
figLegend = pylab.figure(figsize=(11.7, 0.4))
pylab.figlegend(*ax.get_legend_handles_labels(),
                loc='upper center',  ncol=10, shadow=True, fancybox=True)
figLegend.savefig('legend.pdf')
    
plt.figure(num=None, figsize=(4, 3), dpi=80, facecolor='w', edgecolor='k')
cnt = 0
MaxObj = np.zeros([NofOus + 1, NofInstance])
for n in range(NofOus + 1):
    cMaxObj = np.zeros(NofInstance)
    for methods in ThirdOrderMethods:
        FullMethods = (methods + 'WithEdge' + str(WithEdge) +
                       'WithTriplet' + str(True) + Type)
        CObjs = np.array(list(AllObj[n][FullMethods].values()))
        cMaxObj = np.array(np.max(np.array([cMaxObj.tolist(), CObjs.tolist()]),
                                  axis=0))
    MaxObj[n] = cMaxObj
for methods in ThirdOrderMethods:
    FullMethods = (methods + 'WithEdge' + str(WithEdge) +
                   'WithTriplet' + str(True) + Type)
    CMeanObj = np.zeros(NofOus + 1)
    for n in range(NofOus+1):
        CObj = np.array(list(AllObj[n][FullMethods].values()))
        CObj = CObj / MaxObj[n]
        CMeanObj[n] = CObj.mean()
    plt.plot(range(0, NofOus + 1), CMeanObj, label=methods,
             color=Colors[cnt], marker=Markers[cnt],
             linestyle=LineStype[cnt])
    cnt = cnt + 1

for methods in SecondOrderMethods:
    FullMethods = (methods + 'WithEdge' + str(WithEdge) +
                   'WithTriplet' + str(False) + Type)
    CMeanObj = np.zeros(NofOus + 1)
    for n in range(NofOus+1):
        CObj = np.array(list(AllObj[n][FullMethods].values()))
        CObj = CObj / MaxObj[n]
        CMeanObj[n] = CObj.mean()
    plt.plot(range(0, NofOus + 1), CMeanObj, label=methods,
             color=Colors[cnt], marker=Markers[cnt],
             linestyle=LineStype[cnt])
    cnt = cnt + 1
    
pylab.ylim([0.2, 1.05])

plt.ylabel('Obj. Value')
plt.xlabel('#Outliers')
plt.grid(True)

