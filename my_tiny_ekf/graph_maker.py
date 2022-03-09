# -*- coding: utf-8 -*-
from calendar import c
from inspect import _void
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
# Data Load
import serial
import time
import signal
import threading
import ctypes
from ctypes import *
i = c_double(0)
pi = pointer(i)
x = []
y = []

z = []

port = 'COM6'
baud = 9600

exitThread = False


# %%
# ANIMATION FUNCTION
def func(num, dataSet, line, redDots):
    # NOTE: there is no .set_data() for 3 dim data...
    line.set_data(dataSet[0:2, :num])
    line.set_3d_properties(dataSet[2, :num])
    redDots.set_data(dataSet[0:2, :num])
    redDots.set_3d_properties(dataSet[2, :num])
    return line


# %%
def handler(signum, frame):
    exitThread = True


# 데이터 처리할 함수
def parsing_data(data):
    tmp = ''.join(data)
    print(tmp)



# 본 쓰레드
def readThread(ser):
    global line
    global exitThread
    global x
    global y
    global z

    while not exitThread:
        idx = 0
        for c in ser.read():
            if idx % 3 == 0:
                x.append(float(c))
            elif idx % 3 == 1:
                y.append(float(c))
            else:
                z.append(float(c))
            idx = idx + 1


# %% [markdown]
#

# %% [markdown]
#

# %%

# %%

if __name__ == "__main__":
    # %%
    # 종료 시그널 등록
    # signal.signal(signal.SIGINT, handler)
    #
    # ser = serial.Serial(port, baud, timeout=0)
    # if ser.readable():
    #     res = ser.readline()
    #     # print(res)
    #
    # thread = threading.Thread(target=readThread, args=(ser,))
    # thread.start()
    #
    # plot

    #pd array형태로  csv파일 읽어오기
    #####여기서 txt명과 데이터 개수를 적어주세요#####
    ##############################################
    file_name = "circle_test2"
    ##############################################


    #데이터 개수 바꾸려면 dll다시 설정해야함 기본 100개 데이터로 설정해둠.
    data_amount = 100
    #csv로 바꿔주고 9축 데이터 읽어옴
    df = pd.read_csv('./test_files/' + file_name + '.txt', sep = '\t')
    df.to_csv(r'./test_files/'+ file_name + '.csv')
    new_df = pd.read_csv('./test_files/'+ file_name + '.csv')
    m = new_df.values
    #print(m)
    data_matrix1 = m[0:data_amount, 3:6].astype(np.float64)
    data_matrix2 = m[0:data_amount, 6:9].astype(np.float64)
    data_matrix3 = m[0:data_amount, 13:16].astype(np.float64)
    
    print("입력값 A행렬:\n" , data_matrix1)
    print("입력값 W행렬:\n" , data_matrix2)
    print("입력값 H행렬:\n" , data_matrix3)

    #입력 배열 포인터에 할당하기
    filter1 = np.array(data_matrix1, dtype=np.float64)
    pointer_a = filter1.ctypes.data_as(ctypes.POINTER(ctypes.c_double*(data_amount*3)))
    filter2 = np.array(data_matrix2, dtype=np.float64)
    pointer_b = filter2.ctypes.data_as(ctypes.POINTER(ctypes.c_double*(data_amount*3)))
    filter3 = np.array(data_matrix3, dtype=np.float64)
    pointer_c = filter3.ctypes.data_as(ctypes.POINTER(ctypes.c_double*(data_amount*3)))
    
    #ctypes를 이용해서 dll 라이브러리의 함수에 9축 데이터를 3개의 array배열 입력 1개의 array배열 출력
    print("Dll function call")
    libc = ctypes.CDLL('./Dll_lib.dll')
    #함수 입력 형식 900개 double값
    libc.make_string.argtypes = {ctypes.POINTER(ctypes.c_double*(data_amount*3)), ctypes.POINTER(ctypes.c_double*(data_amount*3)), ctypes.POINTER(ctypes.c_double*(data_amount*3))}
    #함수 출력 형식 300개 double값
    libc.make_string.restype = ctypes.POINTER(ctypes.c_double*(data_amount*3))
    arrayptr = libc.make_string(pointer_a, pointer_b, pointer_c)
    c_array = [x for x in arrayptr.contents]
    print("S행렬 출력: ", len(c_array), "개 \n", c_array)
    
    # #여기는 실험영역...
    # # ctypes로 python배열에서 c++ 포인터 배열로 바꾸기. c++로 구현해야함.
    # filter = np.array([[1, 0, 1], [1, 0, 1], [1, -1, 0]], dtype=np.float64)
    # a = filter.ctypes.data_as(ctypes.POINTER(ctypes.c_double*9))
    # print([x for x in a.contents])


    #c_array 값 출력
    idx =0
    for c in c_array:
        if idx %3 ==0 :
            x.append(c)
        elif idx%3 ==1:
            y.append(c)
        elif idx%3 ==2:
            z.append(c)
        idx = idx + 1
    dataSet = np.array([x, y, z])
    #print(x)
    #print(y)
    #print(z)
    numDataPoints = 100

    # GET SOME MATPLOTLIB OBJECTS
    fig = plt.figure()
    ax = Axes3D(fig)
    redDots = plt.plot(dataSet[0], dataSet[1], dataSet[2], lw=2, c='r', marker='o')[0]  # For scatter plot
    # NOTE: Can't pass empty arrays into 3d version of plot()
    line = plt.plot(dataSet[0], dataSet[1], dataSet[2], lw=2, c='g')[0]  # For line plot

    # AXES PROPERTIES]
    ax.set_xlim3d([-10, 10])
    ax.set_ylim3d([-10, 10])
    ax.set_zlim3d([-10, 10])
    ax.set_xlabel('X(t)')
    ax.set_ylabel('Y(t)')
    ax.set_zlabel('Z(t)')
    ax.set_title('Trajectory of electron for E vector along [120]')

    # Creating the Animation object
    line_ani = animation.FuncAnimation(fig, func, frames=numDataPoints, fargs=(dataSet, line, redDots), interval=50,
                                       blit=False)
    # line_ani.save(r'Animation.mp4')

    plt.show()

