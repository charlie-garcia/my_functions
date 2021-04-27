from matplotlib import pyplot as plt

def clean_env():
    try:
        from IPython import get_ipython
        get_ipython().magic('clear')
        get_ipython().magic('reset -f')
        plt.close('all')
        
    except:
        pass

def GetMemory():
    import sys, os, psutil
    process = psutil.Process(os.getpid())
    print('Memory Usage : %.2f MB\n' %(process.memory_info()[0]/(1024**2)))

def VarSize(variable):
    import sys
    print('Variable size : %.2f MB\n' %(sys.getsizeof(variable)/(1024**2)))
