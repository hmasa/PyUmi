# -*- coding: utf-8 -*-
"""
Created on Fri May 15 20:49:09 2020

umi12_AddCircle.py

★N角形に内接する最大の円の半径を求める
☆(1)N角形の頂点のデータを準備する
☆(2)N角形をグラフに図示する
☆(3)N角形の内部に初期の検査エリアTest Area,TA=[Xmin,Xmax,Ymin,Ymax]を決める
☆(4)TAのX,Yを各M等分して（M+1)*（M+1)個の検査座標XYを準備する：M初期値=4
☆(5)[X[i],Y[j]]を中心としたN角形に内接する最大円の半径R[i,j]を求める
☆(6)R[i,j]の中から最大の値を持つR[k]＝R[i,j]max＝R[i_max,j_max]を選ぶ：ｋ初期値＝0
☆(7)Xmin=X[i_max - 1],Xmax=X[i_max + 1],Ymin=Y[j_max - 1],Ymax=X[j_max + 1]
☆(8) (4)～(6)を実行してR[k+1]を求める
☆(9)R[k+1] - R[k] > hantei なら(7)～(8)を繰り返す:kの上限＝50
   R[k+1] - R[k] <= hanteiになったら終了
☆(10)R[k+1]の円を描く
下記の（11）→(12)→（14）:90-緯度　→(15)：経度　→(16)→円をPLOT
    (11)円データの準備：北極を中心とした正N角形の頂点のデータを準備する
    (12)円データ（経度/緯度）を半径1の球上のｘｙｚ座標(SPx/SPy/SPz)に変換
    (13)ｘｙｚ座標をＸ軸を中心としてTHx(radian)回転
    (14)ｘｙｚ座標をＹ軸を中心としてTHy(radian)回転
    (15)ｘｙｚ座標をZ軸を中心としてTHz(radian)回転
    (16)-B:半径1の球上のｘｙｚ座標(SPx/SPy/SPz)を（経度/緯度）データに変換
@author: Hiraiwa
"""

import csv
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from datetime import datetime

# (1)N角形の頂点のデータを準備する
class OpenCsv:  # 2列N行のcsvファイルのデータを読み出す
    def op_csv(CF, newline='', encoding='utf8'):
        f = open(CF)
        dt = [i for i in csv.reader(f, delimiter=',', quotechar='"') ]
        return dt
        f.close
# (2)N角形をグラフに図示する
class PolygonData: # N角形の頂点座標のデータをリストにする
    x = []
    y = []
    def get_pgdata(x,y):    # (N+1)番目に始点データを加え、閉じたN角形にする
        N = len(x)
        x0 = float(x[0])
        y0 = float(y[0])
        x.append(x0)
        y.append(y0)
        
        for i in range(N+1):
            x[i] = float(x[i])  # 文字列データ X,Yを数値データ ｘ、ｙに変換
            y[i] = float(y[i])
        return x,y
    
class FixAxes:  # グラフ描画の時に円が楕円にならないようにX軸とY軸の長さを揃える
    def CalcXYlimit(x,y):
        x_min = min(x)
        x_max = max(x)
        y_min = min(y)
        y_max = max(y)
        x_len = x_max - x_min
        y_len = y_max - y_min     
        x_ctr = (x_max + x_min)/2
        y_ctr = (y_max + y_min)/2
        if x_len > y_len:
            x_left = x_ctr - x_len*0.55
            x_right = x_ctr + x_len*0.55
            y_left = y_ctr - x_len*0.55
            y_right = y_ctr + x_len*0.55
        else:
            x_left = x_ctr - y_len*0.55
            x_right = x_ctr + y_len*0.55
            y_left = y_ctr - y_len*0.55
            y_right = y_ctr + y_len*0.55
        xy_limit = [x_left, x_right, y_left,y_right]
        return xy_limit
    
# (3)N角形の内部に初期の検査エリアTest Area,TA=[Xmin,Xmax,Ymin,Ymax]を決める
class GetXY:
    def GetArea():
        AXmin,AXmax,AYmin,AYmax = (float(x) for x in \
                input('Test Area:Xmin Xmax,Ymin,Ymax>>').split())
        return AXmin,AXmax,AYmin,AYmax
    
# (4)TAのX,Yを各M等分して（M+1)*（M+1)個の検査座標XYを準備する：M初期値=4
class MxN():
    m = 4
    n = 4
    c_xy = []
    mn = []
    def CtrPoints(AXmin,AXmax,AYmin,AYmax):
        m = MxN.m
        n = MxN.n
        Cxy = []
        MN = []
        lx = (AXmax - AXmin)/m
        ly = (AYmax - AYmin)/n
        for i in range(m+1):
            for j in range(n+1):
                Cx = AXmin + i * lx
                Cy = AYmin + j * ly
                Cxy = [Cx, Cy]
                MN.append(Cxy)
        return MN

# (5)[X[i],Y[j]]を中心としたN角形に内接する最大円の半径R[i,j]を求める
class CalcRadius:   # N角形の頂点データ（ｘ[i]、ｙ[i]）と円中心（ｘｃ、ｙｃ）から最小半径を求める
#    x,y = PolygonData.get_pgdata(X,Y)
    def r_min(x,y,xc,yc):
        N = len(x)
        RCx = xc
        RCy = yc
        R = []
        # equatorRadius = 6378.137
        # Rdeg = np.pi*equatorRadius/360
        for i in range(N):
            # 1 Plane surface
            # Ri = np.sqrt((x[i]-RCx)**2 + (y[i]-RCy)**2)
            # 2 Spherical trigonometry
            la1 = np.pi * RCx/180   # dgree→rad変換。np.deg2rad(RCx)より計算速い
            la2 = np.pi * x[i]/180
            th1 = np.pi * RCy/180
            th2 = np.pi * y[i]/180
            Ri = np.arccos(np.sin(th1)*np.sin(th2)\
                           +np.cos(th1)*np.cos(th2)*np.cos(la2-la1))
            R.append(Ri)
        Rmin = min(R)
#        print(R)
        return Rmin * 180/np.pi
    
class DrawCircle:   # 円を描画する
     Xc = 0
     Yc = 0
     r = 1.0
     def DCircle(xc,yc,r0): # 中心（ｘｃ、ｙｃ）、半径ｒの円を描画する
         if r0 <= 0:
             return "r0 is minus value."
         else:
             c = patches.Circle(xy=(xc, yc), radius=r0, fc='w', ec='r')
             return c  

# (6)R[i,j]の中から最大の値を持つR[k]＝R[i,j]max＝R[i_max,j_max]を選ぶ：ｋ初期値＝0
class SelRmax():
    def GetXYR(ctrs):
        Count = len(ctrs)
        XYR = []
        for MN in range(Count):
            XC,YC = ctrs[MN]
            R = CalcRadius.r_min(x,y,XC,YC)
            xyr = [XC,YC,R]
            XYR.append(xyr)
        XYR_all = np.array(XYR)
        XYRmax = np.argmax(XYR_all, axis=0)
        NR = XYRmax[2]
        return XYR[NR]
    
# (11)円データの準備：北極を中心とした正N角形の頂点のデータを準備する
class SpCir:
    Nsp = 360
    # 円周をN等分する価を決める
    def MakeCData(r):
        Csp = []
        Eps = 0.5
        #　経度がゼロになるのを避ける係数
        for i in range (SpCir.Nsp):
           # Xsp = i/SpCir.Nsp * 360
           Xsp = (i - SpCir.Nsp/2 + Eps)/SpCir.Nsp  * 360
           Ysp = 90 - r
           Csp.append([Xsp,Ysp])
        # 最後の要素に０番目の要素を加えて散布図プロットで閉じた図形にする
        # Xsp0 = Csp[0][0]
        # Ysp0 = Csp[0][1]
        # Csp.append([Xsp0,Ysp0])
        # Csp.append(Csp[0])　
        return Csp
    # (12)円データ（経度/緯度）を半径1の球上のｘｙｚ座標(SPx/SPy/SPz)に変換
    def R2XYZ(long,lati,theta):
        SP = []
        count = len(long)
        # print('long =',long)
        # print('lati =',lati)
        for i in range(count):
            # Rx = np.pi*long[i]/180
            # Ry = np.pi*lati[i]/180
            Rx = np.deg2rad(long[i])
            Ry = np.deg2rad(lati[i])
            SPx = np.cos(Rx)*np.cos(Ry)
            SPy = np.sin(Rx)*np.cos(Ry)
            SPz = np.sin(Ry)
            SP.append([SPx,SPy,SPz])
            # print('Rx,Ry =',Rx,Ry)
            # print('K,I =',np.rad2deg(Rx),np.rad2deg(Ry))
        return SP        
    # (13)ｘｙｚ座標をＸ軸を中心としてTHx(radian)回転
    def RotateX(x,y,z,th):  # X軸の周りの回転：
    # X軸=(経度0,緯度0）と(経度180,緯度0）を結ぶ線。東経側はマイナス、西経側がプラス
        X = x
        Y = y*np.cos(th) - z*np.sin(th)
        Z = y*np.sin(th) + z*np.cos(th)
        return X,Y,Z
    # (14)ｘｙｚ座標をＹ軸を中心としてTHy(radian)回転+
    def RotateY(x,y,z,th):  # Y軸の周りの回転：
    # Ｙ軸=(東経90,緯度0）と(西経90,緯度0）を結ぶ線。アメリカ側はマイナス、インド側がプラス
        X = x*np.cos(th) + z*np.sin(th)
        Y = y
        Z = -1*x*np.sin(th) + z*np.cos(th)
        # print('i,th,np.cos,sin=',i,th,np.cos(th),np.sin(th))
        # print('i,X,Y,Z=',i,X,Y,Z)
        return X,Y,Z
    # (15)ｘｙｚ座標をZ軸を中心としてTHz(radian)回転
    def RotateZ(x,y,z,th):  # Z軸の周りの回転：
    # Ｙ軸=(北極,緯度90）と(南極,緯度90）を結ぶ線。西方向はマイナス、東方向がプラス
        X = x*np.cos(th) - y*np.sin(th)
        Y = x*np.sin(th) + y*np.cos(th)
        Z = z
        return X,Y,Z
    # (16)-B:半径1の球上のｘｙｚ座標(SPx/SPy/SPz)を（経度/緯度）データに変換
    def XYZ2R(rx,ry,rz):
        R = []
        count = len(rx)
        for i in range(count):
            Ly = np.rad2deg(np.arcsin(rz[i]))   # np.arcsin()戻り値は-1.57~1.57
            if ry[i] < 0:
                Lx = -90-np.rad2deg(np.arctan(rx[i]/ry[i]))   # np.arctan()戻り値は-∞~∞
            else:
                Lx = 90-np.rad2deg(np.arctan(rx[i]/ry[i]))   # np.arctan()戻り値は-∞~∞
            R.append([Lx,Ly])
        return R    

# 以下、main（）に含まれる記述    
#(1)N角形の頂点のデータを準備する
start_time = datetime.today()   # 計算開始時間を記録
xy = OpenCsv.op_csv('C:antarctica_full.csv')
#xy = OpenCsv.op_csv('C:hokkaido_full.csv')  # Cx0,Cy0=(142.5,43.5)
X = [row[0] for row in xy]    #　行X列の2次元リストから１列目データを抜き出す
Y = [row[1] for row in xy]    #　行X列の2次元リストから2列目データを抜き出す


#(2)N角形をグラフに図示する
x,y = PolygonData.get_pgdata(X,Y)
plt.figure(figsize=(10,10))   # figsize=(6,4)がdefault
# x_left, x_right, y_left,y_right0の値を決める
LR = FixAxes.CalcXYlimit(x,y)

ax = plt.axes() # X軸、Y軸.空枠を描画
ax.set_xlim( left=LR[0], right=LR[1])
ax.set_ylim( bottom=LR[2], top=LR[3])
ax.grid(which = "major", axis = "x", color = "blue", alpha = 0.8,\
        linestyle = "--", linewidth = 1)    # x軸に補助目盛線を設定
ax.grid(which = "major", axis = "y", color = "green", alpha = 0.8,\
        linestyle = "--", linewidth = 1)    # y軸に目盛線を設定

plt.plot(x,y) # Xn,Ynの頂点を持つ多角形をplot


#(3)N角形の内部に初期の検査エリアTest Area,TA=[Xmin,Xmax,Ymin,Ymax]を決める
Area= list(GetXY.GetArea())
# print(Area)
CXmin,CXmax,CYmin,CYmax = [Area[0],Area[1],Area[2],Area[3]]

# (8) (4)～(6)を実行してR[k+1]を求める
hantei = 1e-6
delta = 1
count = 0
R = []
R.append(0)

while delta > hantei:   #(9)R(k+1)-R(k)＜εになったら終了
    # (4)TAのX,Yを各M等分して（M+1)*（M+1)個の検査座標XYを準備する：M初期値=4
    Ctrs = MxN.CtrPoints(CXmin,CXmax,CYmin,CYmax)
    
    # (5)[X[i],Y[j]]を中心としたN角形に内接する最大円の半径R[i,j]を求める
    #    print('円中心のX座標Y座標',XC,YC)
    # MN = 10
    # XC,YC = Ctrs[MN]
    # R = CalcRadius.r_min(x,y,XC,YC)
    # c = DrawCircle.DCircle(XC,YC,R)
    # ax.add_patch(c)
    
    # (6)R[i,j]の中から最大の値を持つR[k]＝R[i,j]max＝R[i_max,j_max]を選ぶ：ｋ初期値＝0
    Rmax = np.array(SelRmax.GetXYR(Ctrs))
    print(Rmax)
    count += 1
    R.append(Rmax[2])
    delta = np.abs(R[count] - R[count-1])
    # Rmax = [CX, CY, R]
    if (delta <= 0):
        CXmin,CXmax,CYmin,CYmax = [Rmax[0]-(CXmax - CXmin)/(MxN.m*2),\
                                   Rmax[0]+(CXmax - CXmin)/(MxN.m*2),\
                                   Rmax[1]-(CYmax - CYmin)/(MxN.n*2),\
                                   Rmax[1]+(CYmax - CYmin)/(MxN.n*2)]
        delta = 1
    else:
    # (7)Xmin=X[i_max - 1],Xmax=X[i_max + 1],Ymin=Y[j_max - 1],Ymax=X[j_max + 1]
        CXmin,CXmax,CYmin,CYmax = [Rmax[0]-(CXmax - CXmin)/MxN.m,\
                                   Rmax[0]+(CXmax - CXmin)/MxN.m,\
                                   Rmax[1]-(CYmax - CYmin)/MxN.n,\
                                   Rmax[1]+(CYmax - CYmin)/MxN.n]
    if (count == 50):
        break
    
    Rmax = np.array(SelRmax.GetXYR(Ctrs))
#     c = DrawCircle.DCircle(Rmax[0],Rmax[1],Rmax[2])
# ax.add_patch(c)

# class data:
#     Rmax0 = 142.7777338
#     Rmax1= 43.42505741
#     Rmax2 = 0.9698363
    
# A = np.arcsin(np.deg2rad(data.Rmax1)) 
# A = 0.82  
# print('A =',A) 
Csp = SpCir.MakeCData(Rmax[2])
# Csp = SpCir.MakeCData(data.Rmax2)
# print('Csp =', Csp)  # 円データ（経度/緯度）の確認
Lx = [row[0] for row in Csp]    #　行X列の2次元リストから１列目データを抜き出す
Ly = [row[1] for row in Csp]    #　行X列の2次元リストから2列目データを抜き出す

# plt.figure(figsize=(10,10))
# plt.plot(Lx,Ly)

Sxyz = SpCir.R2XYZ(Lx,Ly,0)
# print('Sxyz=', Sxyz)  # 円データ（経度/緯度）の確認
Sx = [row[0] for row in Sxyz]    #　行X列の2次元リストから１列目データを抜き出す
Sy = [row[1] for row in Sxyz]    #　行X列の2次元リストから2列目データを抜き出す
Sz = [row[2] for row in Sxyz]    #　行X列の2次元リストから3列目データを抜き出す

Lxy0 = SpCir.XYZ2R(Sx,Sy,Sz)
# print('Lxy0 =',Lxy0)


# Tx = 90.0
# Ty = 90.0-Rmax[1] # 緯線に沿っての回転角(degree)
# Tz = Rmax[0]  # 経線に沿っての回転角(degree)
Tx = 90.0
Ty = 90.0-Rmax[1] # 緯線に沿っての回転角(degree)
Tz = Rmax[0]  # 経線に沿っての回転角(degree)
N = len(Sx)
RX = []
RY = []
RZ = []

# X軸を中心としてTHx(radian)回転
# for i in range(N):
#     RotX = SpCir.RotateX(Sx[i],Sy[i],Sz[i],np.deg2rad(Tx))
#     RX.append(RotX)
# print('RX =',RX)


# Ｙ軸を中心としてTHy(radian)回転
# 緯度を合わせる。緯度＝I(degree)に合わせるには
#　Ty = 90.0-I(北緯の場合)、Ty = 90.0+I(南緯の場合)
for i in range(N):
    RotY = SpCir.RotateY(Sx[i],Sy[i],Sz[i],np.deg2rad(Ty))
    RY.append(RotY)
# print('RY =',RY)

RYx = [row[0] for row in RY]    #　行X列の2次元リストから１列目データを抜き出す
RYy = [row[1] for row in RY]    #　行X列の2次元リストから2列目データを抜き出す
RYz = [row[2] for row in RY]    #　行X列の2次元リストから3列目データを抜き出す

Lxy1 = SpCir.XYZ2R(RYx,RYy,RYz)
# print('Lxy1 =',Lxy1)

# Z軸を中心としてTHz(radian)回転
# 経度を合わせる。経度＝K(degree)に合わせるには
#　Tz = +K(東経の場合)、Tz = -K(西経の場合)
for i in range(N):
    # RotZ = SpCir.RotateZ(Sx[i],Sy[i],Sz[i],thz)
    RotZ = SpCir.RotateZ(RYx[i],RYy[i],RYz[i],np.deg2rad(Tz)) # RY[]データを使う
    RZ.append(RotZ)
# print('RZ =',RZ)

RZx = [row[0] for row in RZ]    #　行X列の2次元リストから１列目データを抜き出す
RZy = [row[1] for row in RZ]    #　行X列の2次元リストから2列目データを抜き出す
RZz = [row[2] for row in RZ]    #　行X列の2次元リストから3列目データを抜き出す

Lxy2 = SpCir.XYZ2R(RZx,RZy,RZz)
# Lxy2 = SpCir.XYZ2R(RYx,RYy,RYz)
# print('Lxy2 =',Lxy2)
K = [row[0] for row in Lxy2]    #　行X列の2次元リストから１列目データを抜き出す
I = [row[1] for row in Lxy2]    #　行X列の2次元リストから2列目データを抜き出す
plt.plot(K,I,'-',color='red') 
plt.show()

end_time = datetime.today() # 計算終了時間を記録
time_delta = end_time - start_time# 計算かかった時間を記録

print("start_time = ",start_time)
print("end_time = ",end_time)
print("time delta = ",time_delta)
