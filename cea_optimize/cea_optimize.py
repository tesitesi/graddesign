# -*- coding: utf-8 -*-

"""
ceaを用いて燃焼圧、混合比を最適化する。
cea-execに含まれるfcea2mが必須。windows用。
コードのいろんな場所にディレクトリ名とか書いていてメンテナンスしにくい(クソ)。
参考: https://qiita.com/ina111/items/f5c4eb35a848fdca04b8
Author :
    Yuki Kumon
Last Update :
    2019-01-08

追記
燃焼圧、混合比からペイロード比のグラフを出すところまで。排気圧も設定可。燃焼圧との比で設定。
冷却計算の準備のため、出口マッハ数や比熱を出力させた。
あとは、特性排気速度c*をCEAから出せるようにすれば冷却計算可能(これの出し方が分からない)。
グラフの凡例がなんか変だけど計算やグラフに支障はないはず。
Author :
    Ryo Sakamoto
Last Update :
    2019-01-07

自分のコードがうんこだったので全面的にパクっていくことにします。。。
Author :
    Yuki Kumon
Last Update :
    2019-01-11
"""

import numpy as np
import matplotlib.pyplot as plt
# from scipy import optimize
import csv
import os
import subprocess

from Pywork.Parameter import Parameter


def cea(f_ratio, p_c, p_ratio):
    file_name = '00temp'  # 一時的に作られるファイル

    # cea実行ファイルのディレクトリ(環境に合わせて書き換える)
    cea_dir = Parameter.cea_dir
    os.chdir(cea_dir)

    # inputファイルの中身(p, bar = 燃焼圧, o/f = 混合比)
    str = """problem    o/f={0},
    rocket  equilibrium  frozen  nfz=2
    p,bar = {1}
    pi/p= {2}
    subar = 3, supar = 2, 3, 5
react
    fuel  CH4 wt%=100 t,k=500
    oxid  O2(L) wt%=100. t,k=90.17
output
    plot o/f p pi/p t cp gam isp ivac ae/at cf machnumber C
end""".format(f_ratio, p_c, p_ratio)

    # inputファイルを書き込み
    input_file = file_name + '.inp'
    f = open(input_file, 'w')  # 書き込みモードで開く
    f.write(str)  # 引数の文字列をファイルに書き込む
    f.close()  # ファイルを閉じる

    # ceaを実行
    cmd = 'FCEA2m'
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    p.communicate((file_name + '\n').encode())  # 改行コード\nがないと動かない。エンコードしないと動かない。

    # csvファイルを読み込む(１行目: 燃焼室, 2行目: スローﾄ: 3行目: スロート出口)
    output_file = file_name + '.csv'
    with open(output_file) as f:
        reader = csv.reader(f, delimiter=",")
        header = next(reader)  # ヘッダーを読み飛ばす
        index = 0

        # 自分でいじるパラメータ
        """
        delta_v = 6433.45  # 速度増分[m/s]
        Mp = 5000  # ペイロード質量[kg]
        eta_s = 0.1  # 構造効率
        marge = 0.05  # 推進剤のマージン
        low_f = 425  # CH4密度(kg/m^3)
        low_o = 1140  # O2密度(kg/m^3)
        """
        delta_v = Parameter.delta_v  # 必要速度増分
        Mp = Parameter.Mp  # ペイロード質量[kg]
        eta_s = Parameter.eta_s  # 構造効率
        marge = Parameter.marge  # 推進剤のマージン
        low_f = Parameter.low_f  # CH4密度(kg/m^3)
        low_o = Parameter.low_o  # O2密度(kg/m^3)
        ge = Parameter.ge  # 重力加速度

        for row in reader:
            if index == 0:
                of = float(row[0])          # 混合比
                P_c = float(row[1]) / 10    # 燃焼室圧力
                T_c = float(row[3])         # 断熱火炎温度
                Cpc = float(row[4]) * 1000  # 比熱(J/kg/K)
                gam_c = float(row[5])       # 比熱比
            if index == 1:
                T_t = float(row[3])          # スロート部温度
                Cpt = float(row[4]) * 1000  # 比熱(J/kg/K)
                gamt = float(row[5])        # 比熱比
            if index == 3:
                T_1 = float(row[3])         # スロート上流の適当な出口温度
                Cp1 = float(row[4]) * 1000  # 比熱(J/kg/K)
                gam1 = float(row[5])        # 比熱
                M_1 = float(row[10])        # マッハ数
                aeat1 = float(row[8])       # 開口比
            if index == 4:
                T_2 = float(row[3])          # 温度
                Cp2 = float(row[4]) * 1000  # 比熱(J/kg/K)
                gam2 = float(row[5])        # 比熱比
                M_2 = float(row[10])        # マッハ数
                aeat2 = float(row[8])       # 開口比
            if index == 5:
                T_3 = float(row[3])          # 温度
                Cp3 = float(row[4]) * 1000  # 比熱(J/kg/K)
                gam3 = float(row[5])        # 比熱比
                M_3 = float(row[10])        # マッハ数
                aeat3 = float(row[8])       # 開口比
            if index == 6:
                T_5 = float(row[3])          # 温度
                Cp5 = float(row[4]) * 1000  # 比熱(J/kg/K)
                gam5 = float(row[5])        # 比熱比
                M_5 = float(row[10])        # マッハ数
                aeat5 = float(row[8])       # 開口比
            if index == 2:
                P_e = float(row[1]) * 100    # 出口圧力
                T_e = float(row[3])          # 出口温度
                Cpe = float(row[4]) * 1000  # 比熱(J/kg/K)
                game = float(row[5])        # 比熱比
                isp = float(row[6]) / ge    # 海上比推力
                ivac = float(row[7]) / ge   # 真空比推力
                aeat_e = float(row[8])        # 開口比
                cf = float(row[9])          # 推力係数
                Me = float(row[10])         # 出口マッハ数
                cstar = float(row[11])      # 特性排気速度

                # ペイロード比の計算
                M1_M2 = np.exp(delta_v / ge / isp)
                M1 = (1 - eta_s) / (1 + marge) * Mp * M1_M2 / (1 - (marge + eta_s) / (1 + marge) * M1_M2)  # ロケット全体質量
                payload_ratio = Mp / M1                     # ペイロード比
                Mpr = (1 - eta_s) * (M1 - Mp)               # 推進剤質量
                Ms = eta_s / (1 - eta_s) * Mpr              # 機体質量
                # W = M1 * 3.71 / 1000                        # 火星での重量[kN]
                # T = 2.5 * W                                 # 必要推力[kN]
                # F = T / 7                                   # エンジン推力[kN]
                # 地上での質量(用いない値)
                W = M1 * ge / 1000
                # 推力(エンジン機数は1)
                T = Parameter.Thrust
                F = T
                m_rate_total = T / float(row[6]) * 1000     # 質量流量[kg/s]
                m_rate = F / float(row[6]) * 1000           # エンジン1基当たりの質量流量[kg/s]
                burn_time = Mpr / m_rate_total              # 噴射時間[s]

                m_cool = m_rate / (1 + of)                      # メタンの質量流量[kg/s]
                M_mol = 16 * (1 + of) / 3                       # 平均分子量(g/mol)
                R_mol = 8.3143 / M_mol * 1000                   # 気体定数(J/kg/K)
                low = (1 + of) / (1 / low_f + of / low_o)        # 平均密度(kg/m^3)
                eta_s_cal = 1 / (93 / 5.6 * low / 10 ** 3 + 1)  # 構造質量比

                A_t = m_rate / P_c / 10 ** 6 * np.sqrt(T_c * R_mol / gam_c * ((gam_c + 1) / 2) ** ((gam_c + 1) / (gam_c - 1)))  # 燃焼室断面積(m^2)
                A_e = A_t * aeat_e              # ノズル出口面積(m^2)
                D_t = np.sqrt(4 * A_t / np.pi)  # スロート径(m)
                D_e = np.sqrt(4 * A_e / np.pi)  # 出口径(m)
                cstar_cal = P_c * 10 ** 6 * A_t / m_rate  # 特性排気速度の計算(m/s)

            index += 1
    data = [of, P_c, T_c, Cpc, gam_c, T_t, Cpt, gamt, P_e, T_e, Cpe, game, isp, ivac, aeat_e, cf, Me, cstar, T_1, Cp1, gam1, M_1, aeat1, T_2, Cp2, gam2, M_2, aeat2, T_3, Cp3, gam3, M_3, aeat3, T_5, Cp5, gam5, M_5, aeat5, payload_ratio, M1, Mpr, Ms, W, T, F, m_rate_total, m_rate, burn_time, m_cool, M_mol, R_mol, low, eta_s_cal, A_t, A_e, D_t, D_e, cstar_cal]
    # print(data)

    # 一時ファイルを削除
    os.remove(input_file)
    os.remove(file_name + '.out')
    os.remove(file_name + '.plt')
    os.remove(file_name + '.csv')

    # 結果を出力
    return data


if __name__ == '__main__':
    # 結果を出力するファイル
    data_file = Parameter.cea_optimize_result_name

    plt.ion()
    plt.figure()

    # パラメータの準備
    # 燃焼圧
    pressures_combustion = Parameter.pressures_combustion
    # ノズル出口圧
    p_e = Parameter.P_e * 10 / 10**6
    # 混合比
    f_ratios = Parameter.f_ratios
    # 書き出し用の配列
    # 混合比→燃焼圧のループ時
    # data_array = np.zeros((len(pressures_combustion), 15))
    # save_array = np.zeros((len(f_ratios), len(pressures_combustion), 15))  # CSVに書き込むデータのストック場所
    # 燃焼圧→混合比のループ時
    data_array = np.zeros((len(f_ratios), 58))
    save_array = np.zeros((len(pressures_combustion), len(f_ratios), 58))  # CSVに書き込むデータのストック場所

    # ループ計算
    # 横軸を混合比にしたいので、燃焼圧→混合比でループ計算する
    num = 0
    for p_c in pressures_combustion:
        k = 0
        p_ratio = p_c / p_e
        for f_ratio in f_ratios:
            data = cea(f_ratio, p_c, p_ratio)
            data_array[k] = data
            k += 1
        # 単位換算
        press = pressures_combustion / 10  # 燃焼室圧力 MPa
        of = data_array[:, 0]  # O/F
        Isp = data_array[:, 13]  # 海上比推力 sec
        # Ivac = data_array[:, 11]  # 真空中比推力 sec
        payload_ratio = data_array[:, 38]
        # print(data_array)
        save_array[num] = data_array
        num += 1
        # プロット
        # plt.plot(press, payload_ratio, label='O/F = {0}'.format(str(round(f_ratio, 1))))
        # plt.plot(of, payload_ratio, label='press = {0}'.format(str(pressures_combustion)))
        # 横軸に混合比、縦軸にペイロード比をとる
        # print(data_array[:, 38])
        if(num % 2 == 1):
            plt.plot(data_array[:, 0], data_array[:, 38], label='Pc = {0} MPa'.format(str(round(p_c / 10, 1))))

    # PLOTの設定
    # plt.xlim(0, 20)
    # plt.xlabel('Pressure (MPa)')
    plt.xlabel('MR')
    plt.ylabel('payload_ratio')
    plt.grid()
    plt.title('LOX/LNG 100%  equilibrium  frozen  nfz=2')
    plt.legend(loc='best', fontsize=5)

    # 書き出しディレクトリに書き出し
    output_dir = Parameter.output_dir
    os.chdir(output_dir)
    if(0):
        plt.savefig('cea_result.png', dpi=500)
    # CSVで結果を保存
    with open(data_file, 'w', newline='') as f:
        writer = csv.writer(f, lineterminator='\n')  # 改行コード（\n）を指定しておく
        header = ['MR', 'Pc [MPa]', 'Tc [K]', 'Cpc [J/kg/K]', 'gam_c', 'T_t [K]', 'Cpt [J/kg/K]', 'gamt', 'Pe [kPa]', 'Te [K]', 'Cpe [J/kg/K]', 'gam_e', 'Isp [sec]', 'Ivac [sec]', 'aeat_e', 'Cf', 'Me', 'c* [m/s]', 'T1 [K]', 'Cp1 [J/kg/K]'
                    , 'gam1', 'M1', 'aeat1', 'T2 [K]', 'Cp2 [J/kg/K]', 'gam2', 'M2', 'aeat2', 'T3 [K]', 'Cp3 [J/kg/K]', 'gam3', 'M3', 'aeat3', 'T5 [K]', 'Cp5 [J/kg/K]', 'gam5', 'M5', 'aeat5', 'payload_ratio', 'M1 [kg]'
                    , 'Mpr [kg]', 'Ms [kg]', 'W [kN]', 'T [kN]', 'F [kN]', 'm_rate_total [kg/s]', 'm_rate [kg/s]', 'burn_time [s]', 'm_cool [kg/s]', 'M_mol [g/mol]', 'R_mol [J/kg/K]', 'low [kg/m^3]', 'eta_s_cal', 'A_t [m^2]', 'A_e [m^2]', 'D_t [m]', 'D_e [m]', 'c*_cal [m/s]']
        writer.writerow(header)
        for i in range(len(save_array[:, 0, 0])):
            writer.writerows(save_array[i])
