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


class Calc2:
    def __init__(self, l0, b, h, a_s, a_s_, As, As_, fc, fy, fy_, Es, e_0):
        """
        不对称配筋偏心受压构件基于承载力的截面校核[已知e_0,求Ncu]
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
            As, As_: 远离/靠近轴力一侧钢筋面积(mm^2)
            fc: 混凝土抗压强度(N/mm^2)
            fy, fy_: 远离/靠近轴向力一侧钢筋抗拉强度(N/mm^2)
            Es: 钢筋弹模(N/mm^2)
            e_0: 初始偏心距(mm)
        """
        self.l0 = l0 * 1e3 # mm
        self.b = b # mm
        self.h = h # mm
        self.a_s = a_s # mm
        self.a_s_ = a_s_ # mm
        self.As = As # mm^2
        self.As_ = As_ # mm^2
        self.fc = fc # N/mm^2
        self.fy = fy # N/mm^2
        self.fy_ = fy_ # N/mm^2
        self.Es = Es # N/mm^2
        self.e_0 = e_0 # mm

        self.e_a = max(20, self.b / 30, self.h / 30)
        self.e_i = self.e_0 + self.e_a
        self.h_0 = self.h - self.a_s
        self.h_0_ = self.h - self.a_s_
        if self.l0 / self.h > 5:
            # self.zeta_c = min(0.5 * self.fc * self.b * self.h / self.Nc, 1)
            self.zeta_c = 1 # [TODO]
            self.eta_s = 1 + (self.l0 / self.h) ** 2 * self.zeta_c / (1300 * self.e_i / self.h_0)
        else:
            self.eta_s = 1
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
        # 先假定为大偏心
        self.xi_todo = (1 - self.e / self.h_0) + sqrt((1 - self.e / self.h_0)**2 + 2 * (self.fy_ * self.As_ * (self.h_0 - self.a_s_ - self.e) + self.fy * self.As * self.e) / (self.alpha1 * self.fc * self.b * self.h_0**2))
        # 确实为大偏心
        if self.xi_todo <= self.xi_b:
            self.status = "大偏心"
            self.xi = self.xi_todo
            # As_ 能屈服
            if (self.xi * self.h_0) >= (2 * self.a_s_):
                self.Ncu_todo = self.alpha1 * self.fc * self.b * self.h_0 * self.xi + self.fy_ * self.As_ - self.fy * self.As
            # As_ 不能屈服, 取 x = 2 a_s_
            else:
                self.Ncu_todo = self.As * self.fy * (self.h_0 - self.a_s_) / (self.eta_s * self.e_i - self.h / 2 + self.a_s_)
        # 实为小偏心
        else:
            self.status = "小偏心"
            self.A = 0.5 * self.alpha1 * self.fc * self.b * self.h_0**2
            self.B = - self.alpha1 * self.fc * self.b * self.h_0 * (self.h_0 - self.e) - self.fy * self.As * self.e / (self.xi_b - 0.8)
            self.C = self.fy_ * self.As_ * (self.e + self.a_s_ - self.h_0) + self.fy * self.As * self.e * 0.8 / (self.xi_b - 0.8)
            self.xi = (-self.B + sqrt(self.B**2 - 4 * self.A * self.C)) / (2 * self.A)
            self.sigma_s = (0.8 - self.xi) / (0.8 - self.xi_b) * self.fy
            self.Ncu_todo = self.alpha1 * self.fc * self.b * self.h_0 * self.xi + self.fy_ * self.As_ - self.sigma_s * self.As
        # 出平面方向校核
        self.varphi = varphi_b(self.l0 / self.b)
        if (self.As + self.As_) / (self.b * self.h) > 0.03:
            self.Ac = self.b * self.h - self.As - self.As_
        else:
            self.Ac = self.b * self.h
        self.Ncu_todo2 = self.varphi * (self.Ac * self.fc + self.As * self.fy + self.As_ * self.fy_)
        self.Ncu = min(self.Ncu_todo, self.Ncu_todo2)
        self.Mu = self.Ncu * self.e_0
    
    def process(self):
        """
        以latex公式形式输出计算过程
        """
        process = f"$$\\begin{{aligned}}"
        if self.fcu > 50:
            process += f"&混凝土强度大于50MPa, 取 \\alpha_1 = {self.alpha1}, \\beta_1 = {self.beta1}, \\varepsilon_{{cu}} = 0.0033 - (f_cu - 50) * 10^{{-5}} = {self.varepsilon_cu}, \\\\"
        process += f"&e_0 = {self.e_0:.2f}mm, e_a = {self.e_a:.2f}mm, e_i = e_0 + e_a = {self.e_i:.2f}mm \\\\"
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
        process += f"&先假定为大偏心受压，联立如下方程 \\\\"
        process += f"&\\begin{{cases}}"
        process += f"&N_{{cu}} = \\alpha_1 f_c b h_0 \\xi + f_y^{{'}} A_s^{{'}} - f_y A_s \\\\"
        process += f"&N_{{cu}} e = \\alpha_1 f_c b h_0^2 \\xi (1 - \\frac{{\\xi}}{{2}}) + f_y^{{'}} A_s^{{'}} (h_0 - a_s^{{'}}) \\\\"
        process += f"\\end{{cases}}\\\\"
        process += f"&解得 \\xi = {self.xi_todo:.4g} \\\\"
        process += f"&\\xi_b = \\frac{{\\beta_1}}{{1+\\frac{{f_y}}{{E_s \\varepsilon_{{cu}}}}}} = {self.xi_b:.4g} \\\\"
        if self.xi_todo <= self.xi_b:
            process += f"&\\xi \\leq \\xi_b, 确实为大偏心受压 \\\\"
            if (self.xi * self.h_0) >= (2 * self.a_s_):
                process += f"&\\xi h_0 = {(self.xi * self.h_0):.2f}mm \\geq 2 a_s^{{'}} = {(2 * self.a_s_):.2f}mm, 故 \\\\"
                process += f"&N_{{cu}} = \\alpha_1 f_c b h_0 \\xi + f_y^{{'}} A_s^{{'}} - f_y A_s = {(self.Ncu_todo * 1e-3):.3f}kN \\\\"
            else:
                process += f"&\\xi h_0 = {(self.xi * self.h_0):.2f}mm < 2 a_s^{{'}} = {(2 * self.a_s_):.2f}mm, A_s^{{'}} 不能屈服, 取 x = 2 a_s^{{'}}, 则 \\\\"
                process += f"&N_{{cu}} = A_s f_y \\frac{{h_0 - a_s^{{'}}}}{{\\eta_s e_i - \\frac{{h}}{{2}} + a_s^{{'}}}} = {(self.Ncu_todo * 1e-3):.3f}kN \\\\"
        else:
            process += f"&\\xi > \\xi_b, 为小偏心受压，联立如下方程 \\\\"
            process += f"&\\begin{{cases}}"
            process += f"&N_{{cu}} = \\alpha_1 f_c b h_0 \\xi + f_y^{{'}} A_s^{{'}} - \\sigma_s A_s \\\\"
            process += f"&N_{{cu}} e = \\alpha_1 f_c b h_0^2 \\xi (1 - \\frac{{\\xi}}{{2}}) + f_y^{{'}} A_s^{{'}} (h_0 - a_s^{{'}}) \\\\"
            process += f"&\\sigma_s = f_y \\frac{{0.8 - \\xi}}{{0.8 - \\xi_b}} \\\\"
            process += f"\\end{{cases}}\\\\"
            process += f"&代入数值化简得 \\xi^2 {(self.B / self.A):+.4g} \\xi {(self.C / self.A):+.4g} = 0, 解得 \\xi = {self.xi:.4g} \\\\"
            process += f"&\\sigma_s = f_y \\frac{{0.8 - \\xi}}{{0.8 - \\xi_b}} = {self.sigma_s:.2f}N/mm^2 \\\\"
            process += f"&N_{{cu}} = \\alpha_1 f_c b h_0 \\xi + f_y^{{'}} A_s^{{'}} - \\sigma_s A_s = {(self.Ncu_todo * 1e-3):.3f}kN \\\\"
        process += f"&\\frac{{l_0}}{{b}} = {(self.l0 / self.b):.2f}, 查表插值得 \\varphi = {self.varphi:.3f} \\\\"
        process += f"&\\frac{{A_s^{{'}} + A_s}}{{b h}} = {((self.As + self.As_) / (self.b * self.h) * 100):.3f}\\%"
        if (self.As + self.As_) / (self.b * self.h) > 0.03:
            process += f" \\geq 3\\%, 取 A_c = b h - A_s^{{'}} - A_s = {self.Ac:.2f}mm^2 \\\\"
        else:
            process += f" < 3\\%, 取 A_c = b h \\\\"
        process += f"&\\varphi (A_c f_c + A_s^{{'}} f_y^{{'}} + A_s f_y) = {(self.Ncu_todo2 * 1e-3):.3f}kN"
        if self.Ncu_todo2 >= self.Ncu_todo:
            process += f" \\geq {(self.Ncu_todo * 1e-3):.3f}kN, 故 N_{{cu}} = {(self.Ncu * 1e-3):.3f}kN \\\\"
        else:
            process += f" < {(self.Ncu_todo * 1e-3):.3f}kN, 故取 N_{{cu}} = {(self.Ncu * 1e-3):.3f}kN \\\\"
        process += f" \\end{{aligned}}$$"
        return process


class Calc3:
    def __init__(self, l0, b, h, a_s, a_s_, As, As_, fc, fy, fy_, Es, Nc):
        """
        不对称配筋偏心受压构件基于承载力的截面校核[已知Nc,求Mu]
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
            As, As_: 远离/靠近轴力一侧钢筋面积(mm^2)
            fc: 混凝土抗压强度(N/mm^2)
            fy, fy_: 远离/靠近轴向力一侧钢筋抗拉强度(N/mm^2)
            Es: 钢筋弹模(N/mm^2)
            Nc: 轴力(kN)
        """
        self.l0 = l0 * 1e3 # mm
        self.b = b # mm
        self.h = h # mm
        self.a_s = a_s # mm
        self.a_s_ = a_s_ # mm
        self.As = As # mm^2
        self.As_ = As_ # mm^2
        self.fc = fc # N/mm^2
        self.fy = fy # N/mm^2
        self.fy_ = fy_ # N/mm^2
        self.Es = Es # N/mm^2
        self.Nc = Nc * 1e3 # N
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
        # 轴压承载力
        self.varphi = varphi_b(self.l0 / min(self.b, self.h))
        if (self.As + self.As_) / (self.b * self.h) > 0.03:
            self.Ac = self.b * self.h - self.As - self.As_
        else:
            self.Ac = self.b * self.h
        if self.Nc > (self.varphi * (self.Ac * self.fc + self.As * self.fy + self.As_ * self.fy_)):
            self.status = "超过了轴压承载力"
            self.Mu = 0
            return
        self.e_a = max(20, self.b / 30, self.h / 30)
        self.h_0 = self.h - self.a_s
        if self.l0 / self.h > 5:
            # self.zeta_c = min(0.5 * self.fc * self.b * self.h / self.Nc, 1)
            self.zeta_c = 1 # [TODO]
            self.eta_s = 1 + (self.l0 / self.h) ** 2 * self.zeta_c / (1300 * self.e_i / self.h_0)
        else:
            self.eta_s = 1
        # 先假定为大偏心受压
        self.xi_todo = (self.Nc - self.fy_ * self.As_ + self.fy * self.As) / (self.alpha1 * self.fc * self.b * self.h_0)
        # 确实为大偏心
        if self.xi_todo <= self.xi_b:
            self.status = "大偏心"
            self.xi = self.xi_todo
            # As_ 能屈服
            if (self.xi * self.h_0) >= (2 * self.a_s_):
                self.status += ", As'能屈服"
                self.Mu = (-self.Nc * (self.eta_s * self.e_a + self.h / 2 - self.a_s) + self.alpha1 * self.fc * self.b * self.h_0**2 * self.xi * (1 - self.xi / 2) + self.fy_ * self.As_ * (self.h_0 - self.a_s_)) / self.eta_s
            # As_ 不能屈服, 取 x = 2 a_s_
            else:
                self.status += ", As'不能屈服"
                self.Mu = (self.fy * self.As * (self.h_0 - self.a_s_) - self.Nc * (self.eta_s * self.e_a - self.h / 2 + self.a_s_)) / self.eta_s
        # 实为小偏心
        else:
            self.status = "小偏心"
            self.xi = (self.Nc + 0.8 * self.fy * self.As / (0.8 - self.xi_b) - self.fy_ * self.As_) / (self.alpha1 * self.fc * self.b * self.h_0 + self.fy * self.As / (0.8 - self.xi_b))
            self.Mu = (self.alpha1 * self.fc * self.b * self.h_0**2 * self.xi * (1 - self.xi / 2) + self.fy_ * self.As_ * (self.h_0 - self.a_s_) - self.Nc * (self.eta_s * self.e_a + self.h / 2 - self.a_s)) / self.eta_s
    
    def process(self):
        process = f"$$\\begin{{aligned}}"
        if self.fcu > 50:
            process += f"&混凝土强度大于50MPa, 取 \\alpha_1 = {self.alpha1}, \\beta_1 = {self.beta1}, \\varepsilon_{{cu}} = 0.0033 - (f_cu - 50) * 10^{{-5}} = {self.varepsilon_cu}, \\\\"
        process += f"&\\frac{{l_0}}{{min(b, h)}} = {(self.l0 / min(self.b, self.h)):.2f}, 查表插值得 \\varphi = {self.varphi:.3f} \\\\"
        process += f"&\\frac{{A_s^{{'}} + A_s}}{{b h}} = {((self.As + self.As_) / (self.b * self.h) * 100):.3f}\\%"
        if (self.As + self.As_) / (self.b * self.h) > 0.03:
            process += f" \\geq 3\\%, 取 A_c = b h - A_s^{{'}} - A_s = {self.Ac:.2f}mm^2 \\\\"
        else:
            process += f" < 3\\%, 取 A_c = b h \\\\"
        process += f"&\\varphi (A_c f_c + A_s^{{'}} f_y^{{'}} + A_s f_y) = {(self.varphi * (self.Ac * self.fc + self.As * self.fy + self.As_ * self.fy_) * 1e-3):.3f}kN"
        if self.Nc > (self.varphi * (self.Ac * self.fc + self.As * self.fy + self.As_ * self.fy_)):
            process += f" < N_c = {(self.Nc * 1e-3)}kN, 即 N_c 超过了构件轴压承载力, M_u = 0 \\\\"
        else:
            process += f" \\geq N_c \\\\"
            process += f"&e_a = {self.e_a:.2f}mm, h_0 = {self.h_0:.2f}mm \\\\"
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
            process += f"&先假定为大偏心受压，由 N_c = \\alpha_1 h_c b h_0 \\xi + f_y^{{'}} A_s^{{'}} - f_y A_s \\\\"
            process += f"&解得 \\xi = {self.xi_todo:.4g} \\\\"
            process += f"&\\xi_b = \\frac{{\\beta_1}}{{1+\\frac{{f_y}}{{E_s \\varepsilon_{{cu}}}}}} = {self.xi_b:.4g} \\\\"
            if self.xi_todo <= self.xi_b:
                process += f"&\\xi \\leq \\xi_b, 确实为大偏心受压 \\\\"
                if (self.xi * self.h_0) >= (2 * self.a_s_):
                    process += f"&\\xi h_0 = {(self.xi * self.h_0):.2f}mm \\geq 2 a_s^{{'}} = {(2 * self.a_s_):.2f}mm, 故 \\\\"
                    process += f"&M_u = N_c e_0 = \\frac{{-N_c (\\eta_s e_a + h/2 - a_s) + \\alpha_1 f_c b h_0^2 \\xi (1 - \\xi / 2) + f_y^{{'}} A_s^{{'}} (h_0 - a_s^{{'}})}}{{\\eta_s}} = {(self.Mu * 1e-6):.3f}kN*m \\\\"
                else:
                    process += f"&\\xi h_0 = {(self.xi * self.h_0):.2f}mm < 2 a_s^{{'}} = {(2 * self.a_s_):.2f}mm, A_s^{{'}} 不能屈服, 取 x = 2 a_s^{{'}}, 则 \\\\"
                    process += f"&M_u = N_c e_0 = \\frac{{f_y A_s (h_0 - a_s^{{'}}) - N_c (\\eta_s e_a - h/2 + a_s^{{'}})}}{{\\eta_s}} = {(self.Mu * 1e-6):.3f}kN*m \\\\"
            else:
                process += f"&\\xi > \\xi_b, 为小偏心受压，联立如下方程 \\\\"
                process += f"&\\begin{{cases}}"
                process += f"&N_{{cu}} = \\alpha_1 f_c b h_0 \\xi + f_y^{{'}} A_s^{{'}} - \\sigma_s A_s \\\\"
                process += f"&\\sigma_s = f_y \\frac{{0.8 - \\xi}}{{0.8 - \\xi_b}} \\\\"
                process += f"\\end{{cases}}\\\\"
                process += f"&解得 \\xi = {self.xi:.4g} \\\\"
                process += f"&M_u = N_c e_0 = \\frac{{\\alpha_1 f_c b h_0^2 \\xi (1 - \\xi / 2) + f_y^{{'}} A_s^{{'}} (h_0 - a_s^{{'}}) - N_c (\\eta_s e_a + h/2 - a_s)}}{{\\eta_s}} = {(self.Mu * 1e-6):.3f}kN*m \\\\"
        process += f" \\end{{aligned}}$$"
        return process
