# -*- coding: utf-8 -*-
# ======
# 多段ロケットの最適質量配分（サイジング）問題の計算
# 必要な軌道速度に空力損失、重力損失、推力損失、制御損失を追加し、
# トータルの⊿Vを事前に算出し、その軌道速度に必要なサイジングを行う。
# 初期検討段階にのみ使用可能。
#
# 入力:
#   各段のIsp[秒]
#   各段の構造比(0.0~1.0)（各段の全備重量と推進剤以外の割合）
#   推力[N](オプション)
#   消費推進剤割合（0.0~1.0）（オプション）(搭載推進剤に対する消費推進剤の割合)
#   投棄物重量[kg]（オプション）
#   エンジン基数[基]（オプション）
#   海面上推力の値を出力・使用するかどうか(bool)（オプション）
#   エンジン1基あたりのノズル出口面積[m2]（オプション）
#
# if __name__ == '__main__':後のUSER INPUT部分を編集すると
# パラメータの変更が可能
#
# cf. 半揚 稔雄(2014) 「ミッション解析と軌道設計の基礎」
# 8.5.1 ロケットの飛翔運動 ― 最適質量配分問題 p.193
#
# Copyright (c) 2016-2017 Takahiro Inagawa
# This code is released under the MIT License.
# http://opensource.org/licenses/mit-license.php
# ======

from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division

import sys
import numpy as np
from scipy import optimize
import pandas as pd

import csv

g0 = 9.80665

class Rocket:
    def __init__(self, Isp, stracture_ratio,
                 thrust=0, propellant_consumption_rate=1.0,
                 jettison=0, number_of_engine=1,
                 use_SeaLevel = False, nozzle_exit_area=0):
        """ロケットクラス
        Args:
            Isp (float) : 真空中比推力[sec]
            stracture_ratio (float) : 構造係数、消費推進剤以外の割合 (0.0~1.0)
            thrust (float, optional) : 推力 [N]
            propellant_consumption_rate (float, optional) : 推進剤消費割合(0.0~1.0)
            jettison (float, optional) : 投棄物質量 [kg]
            number_of_engine (int, optional) : エンジン基数
            use_SeaLevel (bool, optional) : 海面上推力を考慮するかどうか (default : False)
            nozzle_exit_area (float, optional) : エンジン1基あたりのノズル出口面積 [m2/基]
        """
        self.Isp = Isp
        self.s = stracture_ratio
        self.mass = 0
        self.mass_ratio = 0
        self.upper_mass = 0
        self.thrust = thrust
        self.burn_time = 0
        self.acc_liftoff = 0
        self.acc_cutoff = 0
        self.pc_rate = propellant_consumption_rate
        self.residual_propellant = 0
        self.jettison = jettison
        self.num_engine = number_of_engine
        self.use_SeaLevel = use_SeaLevel
        self.nozzle_exit_area = nozzle_exit_area
        self.payload = 0
        self.thrust_SL = 0

    def mass_ratio_calc(self, rambda):
        temp = rambda * g0 * self.Isp
        self.mass_ratio = (1 + temp) / (temp * self.s)

    def mass_calc(self):
        self.mass = (self.mass_ratio - 1) / (1 - self.s * self.mass_ratio) * \
                    self.upper_mass
        self.mass_stracture = self.mass * self.s
        self.mass_prop = self.mass - self.mass_stracture
        self.mass_prop_gross = self.mass_prop / self.pc_rate  # 搭載推進剤
        self.mass_prop_residual = self.mass_prop_gross - self.mass_prop  # 残渣推進剤
        self.mass_stracture_net = self.mass_stracture - self.mass_prop_residual - self.jettison  # 構造重量

    def deltaV_calc(self):
        self.m0mf = (self.mass + self.upper_mass) / (self.mass_stracture + self.upper_mass)
        self.deltaV = self.Isp * g0 * np.log(self.m0mf)

    def burn_time_calc(self):
        if (self.use_SeaLevel):
            self.thrust_SL = self.thrust - self.nozzle_exit_area * self.num_engine * 101300
        self.burn_time = self.mass_prop * self.Isp * g0 / self.thrust
        self.acc_liftoff = self.thrust / g0 / (self.upper_mass + self.mass)
        self.acc_cutoff = self.thrust / g0 / (self.upper_mass + self.mass_stracture)

def sizing_Lagrange_multiplier(x, velocity, rocket_array):
    temp = 1
    Isp_std = rocket_array[0].Isp
    for r in rocket_array:
        temp_Isp = r.Isp / Isp_std
        temp *= ((1.0 / (x * g0 * r.Isp * r.s)) + (1.0/r.s)) ** (temp_Isp)
    temp -= np.exp(velocity / g0 / Isp_std)
    return temp


def main(isp, kouzou):
    # print("==== Optimal Rocket Sizing for Multi-Stage Rocket ====")
    #print("==== 多段ロケットの最適質量配分問題 ====")
    # ==== USER INPUT ====
    rocket_name = "test rocket"
    payload = 1500 # [kg]
    velocity = 10.1680 # [km/s]
    # stage = Rocket(Isp[s], stracture_ratio[-], (optinal:thrust[N]
    #                propellant_consumption_rate [-],
    #                jettison [kg], number_of_engine [-],
    #                use_SeaLevel(bool), nozzle_exit_area [m2])):
    stage1 = Rocket(337.8, 0.10, propellant_consumption_rate=0.95, jettison = 50, number_of_engine=1, use_SeaLevel=False)
    stage2 = Rocket(isp, kouzou, propellant_consumption_rate=0.95, jettison = 0, number_of_engine=1,use_SeaLevel=False)
    r_array = [stage1,stage2]
    # ==== USER INPUT END ====

    total_mass = payload
    total_mass_stracture = 0
    total_mass_prop = 0
    velocity *= 1000 # km/s -> m/s
    min_Isp = np.inf
    for r in r_array:
        min_Isp = min(min_Isp, r.Isp)
    limit = -1.0 / g0 / min_Isp # base of exponential must not be negative
    try:
        sol = optimize.brentq(sizing_Lagrange_multiplier,-100, limit, args = (velocity, r_array)) # solver
    except ValueError:
        print(u"*** NO SOLUTION ***")
        print(u"Please modification Isp and structure ratio")
    temp = payload
    for r in r_array:
        r.mass_ratio_calc(sol)
    for r in reversed(r_array):
        r.upper_mass = temp
        r.mass_calc()
        temp = r.upper_mass + r.mass
        r.deltaV_calc()
        if(r.thrust != 0):
            r.burn_time_calc()
    stage = 0
    data_col = [("項目", "単位"),
                ("質量m0", "kg"), ("質量mf", "kg"),
                ("Isp(vac)", "秒"), ("構造係数", "-"),
                ("質量比", "-"), ("各段質量m0", "kg"),
                ("各段質量mf", "kg"), ("構造重量", "kg"),
                ("残渣推進剤", "kg"), ("投棄物","kg"),
                ("ペイロード", "kg"),
                ("正味の構造効率(構造重量/各段m0)","-"),
                ("消費推進剤重量", "kg"), ("搭載推進剤重量", "kg"), ("推進剤消費率", "%"),
                ("delta V", "m/s"),
                ("推力(vac)", "N"), ("燃焼時間", "秒"),
                ("エンジン数", "基"), ("出口面積", "m2/基"),
                ("推力(S.L.)", "N"),
                ("加速度_ignition", "G"), ("加速度_cutoff", "G")]
    data_col_multi = pd.MultiIndex.from_tuples(data_col, names=['項目', '単位'])
    df = pd.DataFrame(columns=data_col_multi)

    for r in r_array:
        stage += 1
        #print("%d: 質量m0     = %.1f [kg]" % (stage, r.upper_mass + r.mass))
        #print("   質量mf     = %.1f [kg]" % (r.upper_mass + r.mass_stracture))
        #print("   各段質量m0 = %.1f [kg]\t構造比 = %.3f" % (r.mass, r.s))
        #print("   各段質量mf = %.1f [kg]\t質量比 = %.3f" % (r.mass_stracture, r.mass_ratio))
        #print("   Isp(vac)   = %d [s]\t\t消費推進剤質量 = %.1f [kg]" % (r.Isp, r.mass_prop))
        #print("   delta_V    = %d [m/s]" % (r.deltaV))
        #print("   構造重量   = %.1f [kg]" % (r.mass_stracture_net))
        #if(r.pc_rate != 1.0):
            #print("   残渣推進剤 = %.1f [kg]\t推進剤消費率 = %.1f [%%]" % ( r.mass_prop_residual, r.pc_rate * 100))
        #if(r.jettison != 0.0):
            #print("   投棄物     = %.1f [kg]" % (r.jettison))
        #print("   正味の構造効率(構造重量/各段m0) = %.3f" % (1 - r.mass_stracture_net / r.mass))
        #if(r.thrust != 0):
            #print("   ----")
            #print("   推力(vac)   = %d [N]\t燃焼時間 = %d [s]" % (r.thrust, r.burn_time))
            #print("   推力(S.L.)  = %d [N]" % (r.thrust_SL))
            #print("   エンジン基数 = %d [基]\t出口面積 = %.3f [m2/基]" %(r.num_engine, r.nozzle_exit_area))
            #print("   加速度@点火 = %.2f [G]\t加速度@CutOff = %.2f [G]" % (r.acc_liftoff, r.acc_cutoff))
        #print("   ====")
        if stage == 2:
            danm0 = r.upper_mass + r.mass
            dandelV = r.deltaV
        total_mass += r.mass
        total_mass_stracture += r.mass_stracture
        total_mass_prop += r.mass_prop
        # ↓ for output csv file
        output_list = [str(stage)+"段",
                       round(r.upper_mass + r.mass, 1),
                       round(r.upper_mass + r.mass_stracture, 1),
                       r.Isp, r.s, round(r.mass_ratio, 2), round(r.mass, 1),
                       round(r.mass_stracture, 1) , round(r.mass_stracture_net, 1),
                       round(r.mass_prop_residual, 1), round(r.jettison, 1),
                       round(r.payload),
                       round(1 - r.mass_stracture_net / r.mass, 3),
                       round(r.mass_prop, 1), round(r.mass_prop_gross, 1),
                       round(r.pc_rate * 100, 1),
                       round(r.deltaV , 1),
                       round(r.thrust), round(r.burn_time, 1),
                       round(r.num_engine), round(r.nozzle_exit_area, 3),
                       round(r.thrust_SL),
                       round(r.acc_liftoff, 2), round(r.acc_cutoff, 2)]
        df_temp = pd.DataFrame([output_list], columns=data_col)
        df = df.append(df_temp)
    #print("ペイロード: %.1f [kg]" % (payload))
    #print("=====")
    #print("Total: 全備質量   = %d [kg]" % (total_mass))
    #print("       構造質量   = %d [kg]" % (total_mass_stracture))
    #print("       消費推進剤 = %d [kg]" % (total_mass_prop))
    #print("       delta_V    = %.2f [km/s]" % (velocity/1000))
    output_list = ["合計", round(total_mass, 1), "", "", "", "", "", "", "", "", "",
                   round(payload), "", "", "", "", round(velocity), "", "", "", "",
                   "", "", ""]
    df_temp = pd.DataFrame([output_list], columns=data_col)
    df = df.append(df_temp)

    df.T.to_csv("rocket_sizing_"+ rocket_name + ".csv",
                float_format="%.4f", encoding="SHIFT-JIS", header = False)
    
    return [dandelV, total_mass, danm0]




def auto():
    input = np.loadtxt('saitekika_input.csv', delimiter=',',encoding='utf-8-sig')
    sz = len(input)
    output = []
    for koumoku in range(sz):
        out = main(input[koumoku][0],input[koumoku][1])
        print(koumoku)
        output.append(out)
    output = np.array(output)
    np.savetxt('saitekika_output.csv',output, delimiter=',')




if __name__ == '__main__':
    auto()