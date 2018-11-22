import glob,sys,os
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display, HTML
import numpy as np
from astropy.io import fits

def interAssig(value,intervalDict,df_shutterCount, tolerance): 
    #auxiliary function for OverLapCorrection
    for keyInt in intervalDict:
        if (intervalDict[keyInt][0]-tolerance)<= value <=(intervalDict[keyInt][1]+tolerance):
            return keyInt,df_shutterCount['Counts'][keyInt]

def OverLapCorrection(folder_input, folder_output, filename_output, num_windows):
    
    # here fits and txt files are sorted, the last fits is excluded being the SumImg
    sorted_fits= sorted(glob.glob(folder_input+'/*.fits'))[:-1]
    sorted_TXT= sorted(glob.glob(folder_input+'/*.txt'))

    # the output folder is created if non-existing
    if not os.path.exists(folder_output):
        os.makedirs(folder_output)

    #     display(sorted_TXT[0]) #shutter counts
#     display(sorted_TXT[1]) # shutter times
#     display(sorted_TXT[2]) #spectra
#     display(sorted_TXT[3]) # status
    filename_spectra= sorted_TXT[2]
    filename_shuttercount = sorted_TXT[0]
    filename_shuttertime = sorted_TXT[1]
    df = pd.read_csv( filename_spectra,delim_whitespace=True,header=None,
                names=['Time','Spectra'])
    df['diff']=df['Time'].diff()
    df['Spectra']= df['Spectra']/df['Spectra'].max()
    df['diff']= df['diff']/df['diff'].max()

    df_shutterTime = pd.read_csv(filename_shuttertime,delim_whitespace=True,header=None,
                                 names=['t0','t1'], nrows=num_windows)

    df_shutterCount = pd.read_csv(filename_shuttercount,delim_whitespace=True,header=None,
                                 names=['Counts'], nrows=num_windows, index_col=0)

#     display(df_shutterCount)


    df_shutterTime=np.asarray(df_shutterTime.stack(level=[0]))

#     display(df_shutterTime)

    df.plot(x='Time',y=['Spectra','diff'],grid=True,figsize=(8,6))
    sumTime = 0
    TimeArray = np.zeros(num_windows*2)
    index=0
    for i in df_shutterTime:
        print(i)
        sumTime += i
        TimeArray[index] = sumTime
        index += 1
        plt.axvline(x=sumTime)


    plt.show()
#     display(TimeArray)
    interval = {int(key):(value1,value2) for key,value1,value2 in zip(range(len(TimeArray)),TimeArray[::2],TimeArray[1::2])}
    tolerance= (TimeArray[2]-TimeArray[1])/2 #to lower down
    
    dfName = pd.DataFrame(sorted_fits,columns=['name'])
    dfName['ToF'] = df['Time']
    dfName['ShutterWindow']= dfName['ToF'].apply(lambda i:interAssig(i,interval,df_shutterCount, tolerance))
    indexname=0

#     display(dfName.iloc[0,2][1])

# display(dfName[500:550])


    for i in dfName.groupby('ShutterWindow'):
        display(i[0])
    #     display(i[1]['name'])
        array = i[1]['name'].apply(lambda i:fits.open(i)[0].data).as_matrix()
    #     array = array.reshape(len(array),np.shape(array[0])[0],np.shape(array[0])[1])
    #     display(type(array),np.shape(array))
        sumim=np.zeros(np.shape(array[0]))

        for j in range(0,len(array)):
            sumim+=array[j]
    #         display(sumim)

            if j==0:
                P=0
            else:
                P=sumim/i[0][1]

            newim = array[j]/(1-P)
            newim= newim.astype(float)
            filename=filename_output+str(indexname).zfill(5)
            display(filename)
#           indexname+=1
#             filename_long =  i[1].name
#             filename_img = sfilename_long[len(folder_input)+1:]
            fits.writeto(folder_output+filename+'.fits',newim)
#            fits.writeto(folder_output+filename_output+str(i[0][0])+str(j).zfill(5)+'.fits',newim)

    