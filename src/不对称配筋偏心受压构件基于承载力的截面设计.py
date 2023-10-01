from math import floor, sqrt


# 矩形截面钢筋混凝土轴心受压构件的稳定系数 varphi
Varphi_b = [1.0, 0.98, 0.95, 0.92, 0.87, 0.81, 0.75, 0.70, 0.65, 0.60, 0.56, 
            0.52, 0.48, 0.44, 0.40, 0.36, 0.32, 0.29, 0.26, 0.23, 0.21, 0.19]
def varphi_b(num):
    """
        num: 长细比
    """
    if num <= 8:
        return 1.0
    left = Varphi_b[floor(num / 2) - 4]
    right = Varphi_b[floor(num / 2) - 3]
    index = num / 2 - floor(num / 2)
    return left * (1 - index) + right * index


class Calc4:
    def __init__(self, l0, b, h, a_s, a_s_, fc, fy, fy_, Es, Nc, M):
        """
        不对称配筋偏心受压构件基于承载力的截面设计[As_未知]
            l0: 柱长(m)
            b, h: 截面尺寸(mm)
            a_s, a_s_: 远离/靠近轴向力一侧钢筋中心到混凝土边缘的距离(mm)
                e_a: 附加偏心距(mm)
                    《混凝土结构设计规范》(GB 50010)把e_a取为20mm和偏心方向的截面尺寸的1/30二者中的较大值。
                eta_s: 考虑P-Delta效应的弯矩增大系数
                    短柱: l0/h <= 5
                    长柱: 5 < l0/h <= 30
                    细长柱: l0/h > 30
                    对于长柱和细长柱，应考虑此二阶效应。
                    eta_s = 1 + (l0/h)^2 * zeta_c / (1300 * e_i / h0)
                    其中 zeta_c 为截面曲率修正系数: zeta_c = 0.5 * fc * A / Nc 当 zeta_c > 1 时, 取 zeta_c = 1
            fc: 混凝土抗压强度(N/mm^2)
            fy, fy_: 远离/靠近轴向力一侧钢筋抗拉强度(N/mm^2)
            Es: 钢筋弹模(N/mm^2)
            Nc: (kN)
            M: (kN*m)
        """
        self.l0 = l0 * 1e3 # mm
        self.b = b # mm
        self.h = h # mm
        self.a_s = a_s # mm
        self.a_s_ = a_s_ # mm
        self.fc = fc # N/mm^2
        self.fy = fy # N/mm^2
        self.fy_ = fy_ # N/mm^2
        self.Es = Es # N/mm^2
        self.Nc = Nc * 1e3 # N
        self.M = M * 1e6 # N*mm

        self.e_0 = self.M / self.Nc
        self.e_a = max(20, self.b / 30, self.h / 30)
        self.e_i = self.e_0 + self.e_a
        self.h_0 = self.h - self.a_s
        self.h_0_ = self.h - self.a_s_
        if self.l0 / self.h > 5:
            self.zeta_c = min(0.5 * self.fc * self.b * self.h / self.Nc, 1)
            self.eta_s = 1 + (self.l0 / self.h) ** 2 * self.zeta_c / (1300 * self.e_i / self.h_0)
        else:
            self.eta_s = 1
        self.As_min = 0.002 * self.b * self.h
        self.e = self.e_i * self.eta_s + self.h / 2 - self.a_s
        self.solve()
    
    def solve(self):
            # 修正系数
        self.fcu = self.fc * 1.4 / 0.67
        if self.fcu <= 50:
            self.alpha1 = 1.0
            self.beta1 = 0.8
        elif self.fcu >= 80:
            self.alpha1 = 0.94
            self.beta1 = 0.74
        else:
            self.alpha1 = 1 + (0.94 - 1) * (self.fcu - 50) / (80 - 50)
            self.beta1 = 0.8 + (0.74 - 0.8) * (self.fcu - 50) / (80 - 50)
            # 混凝土极限抗压应变
        if self.fcu <= 50:
            self.varepsilon_cu = 0.0033
        else:
            self.varepsilon_cu = 0.0033 - (self.fcu - 50) * 1e-5
            # 界限受压区相对高度
        self.xi_b = self.beta1 / (1 + self.fy / (self.Es * self.varepsilon_cu))
        # 初步判断是大偏心
        if (self.eta_s * self.e_i > 0.3 * self.h_0):
            self.As__todo = (self.Nc * self.e - self.alpha1 * self.fc * self.b * self.h_0**2 * self.xi_b * (1 - self.xi_b / 2)) / (self.fy_ * (self.h_0 - self.a_s_))
            if self.As__todo < self.As_min:
                self.As_ = self.As_min
                self.xi = 1 - sqrt(1 - 2 * (self.Nc * self.e - self.fy_ * self.As_ * (self.h_0 - self.a_s_)) / (self.alpha1 * self.fc * self.b * self.h_0**2))
                """ 我觉得 self.As_ 的增大只会让 self.xi 变小, 小偏心情况不可能出现。
                # 若为小偏心
                if self.xi >= self.xi_b:
                    pass # [TODO]
                """
            else:
                self.As_ = self.As__todo
                self.xi = self.xi_b
            self.As_todo = (self.alpha1 * self.fc * self.b * self.h_0 * self.xi + self.fy_ * self.As_ - self.Nc) / self.fy
            if self.As_todo < self.As_min:
                self.As = self.As_min
            else:
                self.As = self.As_todo
        # 初步判断是小偏心
        else:
            self.As_todo2 = self.As_min
            # 对 self.As 进行补充验算
            self.e_i_ = self.e_0 - self.e_a
            self.e_ = self.h / 2 - self.e_i_ - self.a_s_
            if self.As_todo2 < (self.Nc * self.e_ - self.alpha1 * self.fc * self.b * self.h * (self.h_0_ - self.h / 2)) / (self.fy * (self.h_0_ - self.a_s)):
                # As 过少以至于 Nc 作用下 As 先屈服
                self.As2 = (self.Nc * self.e_ - self.alpha1 * self.fc * self.b * self.h * (self.h_0_ - self.h / 2)) / (self.fy * (self.h_0_ - self.a_s))
            else:
                self.As2 = self.As_todo2
            temp = self.alpha1 * self.fc * self.b * self.h_0**2 * (0.8 - self.xi_b)
            self.B = self.fy * self.As2 * (self.h_0 - self.a_s_) / temp - self.a_s_ / self.h_0
            self.C = (self.Nc * (self.e - self.h_0 + self.a_s_) * (0.8 - self.xi_b) - 0.8 * self.fy * self.As2 * (self.h_0 - self.a_s_)) / temp
            self.xi = -self.B + sqrt(self.B**2 - 2 * self.C)
            # 如果是大偏心
            if self.xi < self.xi_b:
                self.As__todo = (self.Nc * self.e - self.alpha1 * self.fc * self.b * self.h_0**2 * self.xi_b * (1 - self.xi_b / 2)) / (self.fy_ * (self.h_0 - self.a_s_))
                if self.As__todo < self.As_min:
                    self.As_ = self.As_min
                else:
                    self.As_ = self.As__todo
                self.As_todo = (self.alpha1 * self.fc * self.b * self.h_0 * self.xi_b + self.fy_ * self.As_ - self.Nc) / self.fy
                if self.As_todo < self.As_min:
                    self.As = self.As_min
                else:
                    self.As = self.As_todo
            else:
                self.As = self.As2
                self.As__todo = (self.Nc * self.e - self.alpha1 * self.fc * self.b * self.h_0**2 * self.xi * (1 - self.xi / 2)) / (self.fy_ * (self.h_0 - self.a_s_))
                if self.As__todo < self.As_min:
                    self.As_ =  self.As_min
                else:
                    self.As_ = self.As__todo
        # 出平面方向校核
        self.varphi = varphi_b(self.l0 / self.b)
        if (self.As + self.As_) / (self.b * self.h) > 0.03:
            self.Ac = self.b * self.h - self.As - self.As_
        else:
            self.Ac = self.b * self.h
        self.Ncu = self.varphi * (self.Ac * self.fc + self.As * self.fy + self.As_ * self.fy_)
        if self.Ncu < self.Nc:
            pass # 出平面方向校核失败

    def process(self):
        """
        以latex公式形式输出计算过程
        """
        process = f"$$\\begin{{aligned}}"
        if self.fcu > 50:
            process += f"&混凝土强度大于50MPa, 取 \\alpha_1 = {self.alpha1}, \\beta_1 = {self.beta1}, \\varepsilon_{{cu}} = 0.0033 - (f_cu - 50) * 10^{{-5}} = {self.varepsilon_cu}, \\\\"
        process += f"&(1)计算 e_i, \\eta_s, e \\\\"
        process += f"&e_0 = \\frac{{M}}{{N}} = {self.e_0:.2f}mm, e_a = {self.e_a:.2f}mm, e_i = e_0 + e_a = {self.e_i:.2f}mm \\\\"
        if self.l0 / self.h > 5:
            process += f"&\\frac{{l_0}}{{h}} = {(self.l0 / self.h):.2f} > 5, 需要计算 \\eta_s 值 \\\\"
            process += f"&\\zeta_c = \\frac{{0.5 f_c A}}{{N_c}} = {(0.5 * self.fc * self.b * self.h / self.Nc):.3f}"
            if (0.5 * self.fc * self.b * self.h / self.Nc) > 1:
                process += f" > 1.0, 取 \\zeta_c = 1.0"
            else:
                process += f" \\leq 1.0"
            process += f", h_0 = h - a_s = {self.h_0:.2f}mm, \\\\"
            process += f"&\\eta_s = 1 + \\frac{{1}}{{1300 * \\frac{{e_i}}{{h_0}}}} (\\frac{{l_0}}{{h}})^2 \\zeta_c = {self.eta_s:.3f} \\\\"
        else:
            process += f"&\\frac{{l_0}}{{h}} = {(self.l0 / self.h):.2f} \\leq 5, 取 \\eta_s = 1.0 \\\\"
        process += f"&e = \\eta_s e_i + \\frac{{h}}{{2}} - a_s = {self.e:.2f}mm \\\\"
        process += f"&(2)判断大小偏心受压 \\\\"
        if (self.eta_s * self.e_i > 0.3 * self.h_0):
            process += f"&\\eta_s e_i = {(self.eta_s * self.e_i):.2f}mm > 0.3 h_0 = {(0.3 * self.h_0):.2f}mm, 初步判断为大偏心受压。\\\\"
            process += f"&(3)计算 A_s^{{'}}, A_s \\\\"
            process += f"&\\xi_b = \\frac{{\\beta_1}}{{1+\\frac{{f_y}}{{E_s \\varepsilon_{{cu}}}}}} = {self.xi_b:.4g} \\\\"
            process += f"&取 \\xi = \\xi_b = {self.xi_b:.4g} \\\\"
            process += f"&A_s^{{'}} = \\frac{{N_c e - \\alpha_1 f_c b h_0^2 \\xi (1 - 0.5 \\xi)}}{{f_y^{{'}} (h_0 - a_s^{{'}})}} = {self.As__todo:.2f}mm^2 \\\\"
            process += f"&\\rho_{{min}} b h = 0.002 * {self.b} * {self.h} = {self.As_min:.2f}mm^2, "
            if self.As__todo < self.As_min:
                process += f"A_s^{{'}} < \\rho_{{min}} b h, 取 A_s^{{'}} = {self.As_:.2f}mm^2 \\\\"
                process += f"&\\xi = 1 - \\sqrt{{1 - 2 * \\frac{{N_c e - f_y^{{'}} A_s^{{'}} (h_0 - a_s^{{'}})}}{{\\alpha_1 f_c b h_0^2}}}} = {self.xi:.4g} < \\xi_b, 确为大偏心受压。\\\\"
            else:
                process += f"A_s^{{'}} \\geq \\rho_{{min}} b h, 可行 \\\\"
            process += f"&A_s = \\frac{{\\alpha_1 f_c b h_0 \\xi + f_y^{{'}} A_s^{{'}} - N_c}}{{f_y}} = {self.As_todo:.2f}mm^2"
            if self.As_todo < self.As_min:
                process += f" < \\rho_{{min}} b h, 取 A_s = {self.As:.2f}mm^2 \\\\"
            else:
                process += f" \\geq \\rho_{{min}} b h, 可行 \\\\"
        else:
            process += f"&\\eta_s e_i = {(self.eta_s * self.e_i):.2f}mm \\leq 0.3 h_0 = {(0.3 * self.h_0):.2f}mm, 初步判断为小偏心受压。\\\\"
            process += f"&(3)计算 A_s^{{'}}, A_s \\\\"
            process += f"&\\xi_b = \\frac{{\\beta_1}}{{1+\\frac{{f_y}}{{E_s \\varepsilon_{{cu}}}}}} = {self.xi_b:.4g} \\\\"
            process += f"&A_s = \\rho_{{min}} b h = {self.As_min:.2f}mm^2 \\\\"
            process += f"&对 A_s 补充验算 \\\\"
            process += f"&e_i^{{'}} = e_0 - e_a = {self.e_i_:.2f}mm, e^{{'}} = \\frac{{h}}{{2}} - e_i^{{'}} - a_s^{{'}} = {self.e_:.2f}mm, h_0^{{'}} = h - a_s^{{'}} = {self.h_0_:.2f}mm \\\\"
            process += f"&A_s \\geq \\frac{{N_c e^{{'}} - \\alpha_1 f_c b h (h_0^{{'}} - 0.5 h)}}{{f_y (h_0^{{'}} - a_s)}} = {(self.Nc * self.e_ - self.alpha1 * self.fc * self.b * self.h * (self.h_0_ - self.h / 2)) / (self.fy * (self.h_0_ - self.a_s)):.2f}mm^2 \\\\"
            if self.As_todo2 < (self.Nc * self.e_ - self.alpha1 * self.fc * self.b * self.h * (self.h_0_ - self.h / 2)) / (self.fy * (self.h_0_ - self.a_s)):
                process += f"&由最小配筋率确定的 A_s 偏小，取 A_s = {self.As2:.2f}mm^2 \\\\"
            else:
                process += f"&验算通过。\\\\"
            process += f"&联立如下三个方程求解：\\\\"
            process += f"&\\begin{{cases}}"
            process += f"&N_c = \\alpha_1 f_c b h_0 \\xi + f_y^{{'}} A_s^{{'}} - \\sigma_s A_s \\\\"
            process += f"&N_c e = \\alpha_1 f_c b h_0^2 \\xi (1 - 0.5 \\xi) + f_y^{{'}} A_s^{{'}} (h_0 - a_s^{{'}}) \\\\"
            process += f"&\\sigma_s = f_y \\frac{{\\xi - 0.8}}{{\\xi_b - 0.8}} \\\\"
            process += f"\\end{{cases}}\\\\"
            process += f"&解得 \\xi = - B + \\sqrt{{B^2 - 2 C}}, 其中 \\\\"
            process += f"&B = \\frac{{f_y A_s (h_0 - a_s^{{'}})}}{{\\alpha_1 f_c b h_0^2 (0.8 - \\xi_b)}} - \\frac{{a_s^{{'}}}}{{h_0}} = {self.B:.4g} \\\\"
            process += f"&C = \\frac{{N_c (e - h_0 + a_s^{{'}}) (0.8 - \\xi_b) - 0.8 f_y A_s (h_0 - a_s^{{'}})}}{{\\alpha_1 f_c b h_0^2 (0.8 - \\xi_b)}} = {self.C:.4g} \\\\"
            process += f"&\\xi = - {self.B:.4g} + \\sqrt{{{self.B:.4g}^2 - 2 * ({self.C:.4g})}} = {self.xi:.4g}"
            if self.xi < self.xi_b:
                process += f" < \\xi_b, 为大偏心受压，需重新计算 A_s^{{'}}, A_s: \\\\"
                process += f"&取 \\xi = \\xi_b = {self.xi_b:.4g} \\\\"
                process += f"&A_s^{{'}} = \\frac{{N_c e - \\alpha_1 f_c b h_0^2 \\xi (1 - 0.5 \\xi)}}{{f_y^{{'}} (h_0 - a_s^{{'}})}} = {self.As__todo:.2f}mm^2 \\\\"
                process += f"&\\rho_{{min}} b h = 0.002 * {self.b} * {self.h} = {self.As_min:.2f}mm^2, "
                if self.As__todo < self.As_min:
                    process += f"&A_s^{{'}} < \\rho_{{min}} b h, 取 A_s^{{'}} = {self.As_:.2f}mm^2 \\\\"
                else:
                    process += f"&A_s^{{'}} \\geq \\rho_{{min}} b h, 可行"
                process += f"&A_s = \\frac{{\\alpha_1 f_c b h_0 \\xi + f_y^{{'}} A_s^{{'}} - N_c}}{{f_y}} = {self.As_todo:.2f}mm^2"
                if self.As_todo < self.As_min:
                    process += f" < \\rho_{{min}} b h, 取 A_s = {self.As:.2f}mm^2 \\\\"
                else:
                    process += f" > \\rho_{{min}} b h, 可行 \\\\"
            else:
                process += f" > \\xi_b, 确实为小偏心受压，故 \\\\"
                process += f"&A_s^{{'}} = \\frac{{N_c e - \\alpha_1 f_c b h_0^2 \\xi (1 - 0.5 \\xi)}}{{f_y^{{'}} (h_0 - a_s^{{'}})}} = {self.As__todo:.2f}mm^2"
                if self.As__todo < self.As_min:
                    process += f" < \\rho_{{min}} b h, 取 A_s^{{'}} = {self.As_:.2f}mm^2 \\\\"
                else:
                    process += f" \\geq \\rho_{{min}} b h, 可行 \\\\"
        process += f"&(4)出平面方向验算 \\\\"
        process += f"&\\frac{{l_0}}{{b}} = {(self.l0 / self.b):.2f}, 查表插值得 \\varphi = {self.varphi:.3f} \\\\"
        process += f"&\\frac{{A_s^{{'}} + A_s}}{{b h}} = {((self.As + self.As_) / (self.b * self.h) * 100):.3f}\\%"
        if (self.As + self.As_) / (self.b * self.h) > 0.03:
            process += f" \\geq 3\\%, 取 A_c = b h - A_s^{{'}} - A_s = {self.Ac:.2f}mm^2 \\\\"
        else:
            process += f" < 3\\%, 取 A_c = b h \\\\"
        process += f"&N_{{cu}} = \\varphi (A_c f_c + A_s^{{'}} f_y^{{'}} + A_s f_y) = {(self.Ncu / 1e3):.3f}kN"
        if self.Ncu < self.Nc:
            process += f" < N_c, 不满足要求。\\\\"
        else:
            process += f" \\geq N_c, 满足要求。\\\\"
        process += f"\\end{{aligned}}$$"
        return process


class Calc5:
    def __init__(self, l0, b, h, a_s, a_s_, As_, fc, fy, fy_, Es, Nc, M):
        """
        不对称配筋偏心受压构件基于承载力的截面设计[As_已知]
            l0: 柱长(m)
            b, h: 截面尺寸(mm)
            a_s, a_s_: 远离/靠近轴向力一侧钢筋中心到混凝土边缘的距离(mm)
                e_a: 附加偏心距(mm)
                    《混凝土结构设计规范》(GB 50010)把e_a取为20mm和偏心方向的截面尺寸的1/30二者中的较大值。
                eta_s: 考虑P-Delta效应的弯矩增大系数
                    短柱: l0/h <= 5
                    长柱: 5 < l0/h <= 30
                    细长柱: l0/h > 30
                    对于长柱和细长柱，应考虑此二阶效应。
                    eta_s = 1 + (l0/h)^2 * zeta_c / (1300 * e_i / h0)
                    其中 zeta_c 为截面曲率修正系数: zeta_c = 0.5 * fc * A / Nc 当 zeta_c > 1 时, 取 zeta_c = 1
            As_: 靠近轴力一侧钢筋面积(mm^2)
            fc: 混凝土抗压强度(N/mm^2)
            fy, fy_: 远离/靠近轴向力一侧钢筋抗拉强度(N/mm^2)
            Es: 钢筋弹模(N/mm^2)
            Nc: (kN)
            M: (kN*m)
        """
        self.l0 = l0 * 1e3 # mm
        self.b = b # mm
        self.h = h # mm
        self.a_s = a_s # mm
        self.a_s_ = a_s_ # mm
        self.As__todo = As_ # mm^2
        self.fc = fc # N/mm^2
        self.fy = fy # N/mm^2
        self.fy_ = fy_ # N/mm^2
        self.Es = Es # N/mm^2
        self.Nc = Nc * 1e3 # N
        self.M = M * 1e6 # N*mm

        self.e_0 = self.M / self.Nc
        self.e_a = max(20, self.b / 30, self.h / 30)
        self.e_i = self.e_0 + self.e_a
        self.h_0 = self.h - self.a_s
        self.h_0_ = self.h - self.a_s_
        if self.l0 / self.h > 5:
            self.zeta_c = min(0.5 * self.fc * self.b * self.h / self.Nc, 1)
            self.eta_s = 1 + (self.l0 / self.h) ** 2 * self.zeta_c / (1300 * self.e_i / self.h_0)
        else:
            self.eta_s = 1
        self.As_min = 0.002 * self.b * self.h
        self.e = self.e_i * self.eta_s + self.h / 2 - self.a_s
        self.solve()

    def solve(self):
        if self.As__todo < self.As_min:
            self.subCalc = Calc4(self.l0 * 1e-3, self.b, self.h, self.a_s, self.a_s_, self.fc, self.fy, self.fy_, self.Es, self.Nc * 1e-3, self.M * 1e-6)
            return
            # 修正系数
        self.fcu = self.fc * 1.4 / 0.67
        if self.fcu <= 50:
            self.alpha1 = 1.0
            self.beta1 = 0.8
        elif self.fcu >= 80:
            self.alpha1 = 0.94
            self.beta1 = 0.74
        else:
            self.alpha1 = 1 + (0.94 - 1) * (self.fcu - 50) / (80 - 50)
            self.beta1 = 0.8 + (0.74 - 0.8) * (self.fcu - 50) / (80 - 50)
            # 混凝土极限抗压应变
        if self.fcu <= 50:
            self.varepsilon_cu = 0.0033
        else:
            self.varepsilon_cu = 0.0033 - (self.fcu - 50) * 1e-5
            # 界限受压区相对高度
        self.xi_b = self.beta1 / (1 + self.fy / (self.Es * self.varepsilon_cu))
        # 初步判断是大偏心
        if (self.eta_s * self.e_i > 0.3 * self.h_0):
            self.As_ = self.As__todo
            self.xi = 1 - sqrt(1 - 2 * (self.Nc * self.e - self.fy_ * self.As_ * (self.h_0 - self.a_s_)) / (self.alpha1 * self.fc * self.b * self.h_0**2))
            # 确实是大偏心
            if self.xi <= self.xi_b:
                if (self.xi * self.h_0) >= (2 * self.a_s_):
                    self.As_todo = (self.alpha1 * self.fc * self.b * self.h_0 * self.xi + self.fy_ * self.As_ - self.Nc) / self.fy
                # 说明 As_ 不能屈服
                else:
                    self.x = 2 * self.a_s_
                    self.As_todo = self.Nc * (self.e_i * self.eta_s - self.h / 2 + self.a_s_) / (self.fy * (self.h_0 - self.a_s_))
                if self.As_todo < self.As_min:
                    self.As = self.As_min
                else:
                    self.As = self.As_todo
            # 其实是小偏心
            else:
                pass # [TODO]
        # 初步判断是小偏心
        else:
            pass # [TODO]
    
    def process(self):
        pass # [TODO]
