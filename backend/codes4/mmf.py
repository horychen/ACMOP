import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import square
import os
import pickle

if True:
    # 输入参数
    Qs = 12                                                                                                     # 定子槽数
    p = 4                                                                                                       # 极对数
    m = 3                                                                                                       # 相数
    y = 1                                                                                                      # 节距
    alpha_u = 2 * np.pi / Qs                                                                                    # 槽距角 [rad]
    N = 30000                                                                                                   # 步数
    alpha = np.linspace(0, 2 * np.pi, N)                                                                        # 机械角度，机械角度划为N格
    alpha_ref = np.arange(0, 2 * np.pi, alpha_u)                                                                # 机械角度划分为Qs等份
                                                            
    # 电流幅值（瞬时值）
    i_u = 1
    i_v = 0*-0.5
    i_w = 0*-0.5

    # 标幺化的U相线电流密度[A/m]
    I_U = np.array([-1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1])                                             # U phase winding set
    ZQ_U = np.ones(Qs) * 2                                                                                      # 每槽导体数

    # U相电流密度[A/m]
    A_U = I_U * ZQ_U

    # 计算U相MMF
    MMF_U = np.zeros((N, Qs))
    for i in range(Qs):
        MMF_U[:, i] = A_U[i] * np.where(alpha >= alpha_ref[i], 1, -1)

    # 计算总的U相MMF
    MMF_U_pha = np.sum(MMF_U, axis=1)

    # 将U相平移得到V、W相MMF
    Delta_N = N // 3                                                                                            # 计算相位移动的点数
    MMF_V_pha = np.roll(MMF_U_pha, Delta_N)
    MMF_W_pha = np.roll(MMF_U_pha, -Delta_N)
    
    # 计算总的三相MMF
    MMF_tot = MMF_U_pha * i_u + MMF_V_pha * i_v + MMF_W_pha * i_w
    MMF_tot_norm = (MMF_tot - np.min(MMF_tot)) / (np.max(MMF_tot) - np.min(MMF_tot)) * 2 - 1                    # 归一化到[-1, 1]

    # 频谱分析
    MMF_y = np.fft.fft(MMF_tot_norm)                                                                            # 傅立叶变换
    # MMF_y = np.fft.fft(MMF_tot)                                                                            # 傅立叶变换（可选非归一化）
    P2 = np.abs(MMF_y / N)                                                                                      # 双边谱
    P1 = P2[:N//2+1]                                                                                            # 单边谱
    P1[1:-1] *= 2                                                                                               # 除了直流项，其他频率需要乘以2
    Order = np.arange(0, len(P1))               # 单边阶次
    
    # 存储列表到 pickle 文件
    data = np.vstack((Order, P1))
    with open('mmf.npy', 'wb') as f:
        np.save(f, data)
    # print(data)
    # exit()                                                                

    # 画MMF图
    plt.figure(figsize=(8, 6))
    # Fixed plotting: Proper subplots and grid for each phase, add labels and titles
    plt.subplot(4, 1, 1)
    plt.plot(alpha * 180 / np.pi, MMF_U_pha, color='crimson', linewidth=3, label='U phase')
    plt.grid(True)
    plt.legend(fontsize=12)
    plt.ylabel('MMF [A]', fontfamily='Times New Roman', fontsize=14)
    plt.title('U Phase MMF', fontsize=14, fontweight='bold')

    plt.subplot(4, 1, 2)
    plt.plot(alpha * 180 / np.pi, MMF_V_pha, color='royalblue', linewidth=2, label='V phase')
    plt.grid(True)
    plt.legend(fontsize=12)
    plt.ylabel('MMF [A]', fontfamily='Times New Roman', fontsize=14)
    plt.title('V Phase MMF', fontsize=14, fontweight='bold')

    plt.subplot(4, 1, 3)
    plt.plot(alpha * 180 / np.pi, MMF_W_pha, color='seagreen', linewidth=2, label='W phase')
    plt.grid(True)
    plt.legend(fontsize=12)
    plt.ylabel('MMF [A]', fontfamily='Times New Roman', fontsize=14)
    plt.title('W Phase MMF', fontsize=14, fontweight='bold')
        # plt.ylim([-1.1, 1.1])

    # 非归一化内容
    if False:
        plt.subplot(2, 1, 1)
        plt.plot(alpha * 180 / np.pi, MMF_tot, 'crimson', linewidth=3)
        plt.grid(True)
        plt.xlabel('Mechanical position [degree]', fontfamily='Times New Roman', fontsize=14)
        plt.ylabel('MMF [A]', labelpad=13, fontfamily='Times New Roman', fontsize=14)
        plt.xlim([0, 360])
        plt.ylim([-1.1, 1.1])

    # 画MMF频谱图
    plt.subplot(2, 1, 2)
    plt.bar(Order, P1, color='crimson')
    plt.grid(True)
    plt.xlabel('Harmonic Order', fontfamily='Times New Roman', fontsize=14)
    plt.ylabel('MMF [A]', labelpad=20, fontfamily='Times New Roman', fontsize=14)
    plt.xticks(np.arange(0, 41, 1))
    plt.xlim([0, 40])
    plt.ylim([0, 1])


    output_dir = "output_images"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    filename = f"MMF_plot_p{p}Qs{Qs}y{y}.png"

    plt.savefig(os.path.join(output_dir, filename), dpi=400, bbox_inches='tight', pad_inches=0)

    plt.tight_layout()
    plt.show()
    quit()

    max_harmonic_order = 40
    limited_order = Order[:max_harmonic_order]
    limited_P1 = P1[:max_harmonic_order]
    # 打印谐波分量
    for i, amplitude in enumerate(limited_P1):
        print(f"谐波次序 {i}: 幅值 = {amplitude}")

    # verification of the result: 2013 Pavel Ponomarev

else:     
# 输入参数
    Qs = 12                  # 定子槽数
    p = 4                   # 极对数
    m = 3                   # 相数
    y = 1                     # 节距
    alpha_u = 2 * np.pi / Qs  # 槽距角 [rad]
    N = 10000               # 步数
    alpha = np.linspace(0, 2 * np.pi, N)  # 机械角度划分为N格
    alpha_ref = np.linspace(0, 2 * np.pi, Qs+1)  # 机械角度划分为Qs等份

    # 给定U、V、W三相电流幅值
    i_u, i_v, i_w = 1, 0, 0

    # 标幺化的U相线电流密度[A/m]
    I_U = [-1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1]
    # 每槽导体数
    ZQ_U = [2] * Qs

    A_U = np.array(I_U) * np.array(ZQ_U)  # U相电流密度[A/m]

    # 计算U相MMF
    MMF_U = np.zeros((N, Qs))

    for i in range(Qs):
        MMF_U[:, i] = np.where(alpha > alpha_ref[i], A_U[i], -A_U[i])

    MMF_U_pha = np.sum(MMF_U, axis=1)  # 这里还未乘以电流幅值

    # 将U相平移得到V、W相MMF
    Detla_N = round(N / 3)  # 移相的点数
    MMF_V_pha = np.roll(MMF_U_pha, Detla_N)  # V相为在U相的基础上进行平移
    MMF_W_pha = np.roll(MMF_V_pha, Detla_N)  # W相为在V相的基础上进行平移

    MMF_tot = MMF_U_pha * i_u + MMF_V_pha * i_v + MMF_W_pha * i_w  # 三相总MMF

    # 归一化到[-1,1]
    MMF_tot_norm = 2 * (MMF_tot - np.min(MMF_tot)) / (np.max(MMF_tot) - np.min(MMF_tot)) - 1

    # 频谱分析
    MMF_y = np.fft.fft(MMF_tot_norm)  # 傅立叶变换
    P2 = np.abs(MMF_y / N)  # 双边谱
    P1 = P2[:N//2+1]  # 单边谱
    P1[1:-1] *= 2  # 除了直流项，其他频率需要乘以2
    Order = np.arange(len(P1))  # 单边阶次

    # 计算绕组系数
    kw = P1 * Order * np.pi * m / Qs / 2  # 参照论文 单相MMF

    # 画MMF
    plt.figure(figsize=(10, 8))
    plt.subplot(3, 1, 1)
    plt.plot(alpha, MMF_tot_norm, 'r', markerfacecolor=[1, 0, 0], linewidth=1)
    plt.grid(True, linestyle=':', color='k', alpha=1)
    plt.xlabel('Mechanical position [rad]', fontsize=14, fontname='Times New Roman')
    plt.ylabel('Current linkage [A]', fontsize=14, fontname='Times New Roman')
    plt.xlim([0, 2 * np.pi])
    plt.xticks(np.arange(0, 2 * np.pi + 1, 1))
    plt.ylim([-1.1, 1.1])
    plt.yticks(np.arange(-1, 1.1, 0.5))

    # 画MMF频谱
    plt.subplot(3, 1, 2)
    plt.bar(Order, P1, color='r')
    plt.grid(True, linestyle=':', color='k', alpha=1)
    plt.xlabel('Order', fontsize=14, fontname='Times New Roman')
    plt.ylabel('MMF [A]', fontsize=14, fontname='Times New Roman')
    plt.xlim([0, 30])
    plt.xticks(np.arange(0, 41, 1))
    plt.yticks(np.arange(min(P1), max(P1), 0.2))

    # 画绕组系数
    plt.subplot(3, 1, 3)
    plt.bar(Order, kw, color='r')
    plt.grid(True, linestyle=':', color='k', alpha=1)
    plt.xlabel('Order', fontsize=14, fontname='Times New Roman')
    plt.ylabel('Winding factor', fontsize=14, fontname='Times New Roman')
    plt.xlim([0, 30])
    plt.xticks(np.arange(0, 31, 1))
    plt.yticks(np.arange(min(kw), max(kw), 0.2))

    output_dir = "output_images"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    filename = f"MMF_plot_p{p}Qs{Qs}y{y}.png"

    plt.savefig(os.path.join(output_dir, filename), dpi=400, bbox_inches='tight', pad_inches=0)

    plt.tight_layout()
    plt.show()

    max_harmonic_order = 40
    limited_order = Order[:max_harmonic_order]
    limited_P1 = P1[:max_harmonic_order]
    # 打印谐波分量
    for i, amplitude in enumerate(limited_P1):
        print(f"Harmonic order {i}: Amplitude = {amplitude}")
    # print(min(kw), max(kw))