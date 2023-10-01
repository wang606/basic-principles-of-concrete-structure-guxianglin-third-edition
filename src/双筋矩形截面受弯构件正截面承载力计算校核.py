class Calc8:
    def __init__(self, b, h, h0, as_, As, As_, fc, ft, fy, fy_, Es):
        """
        双筋矩形截面受弯构件正截面承载力计算校核
            b, h: 截面尺寸(mm)
            h0: 有效高度(mm)
            as_: 受压区钢筋中心线到截面受压区边缘的距离(mm)
            As: 受拉区配筋面积(mm^2)
            As_: 受压区配筋面积(mm^2)
            fc: 混凝土抗压强度(N/mm^2)
            ft: 混凝土抗拉强度(N/mm^2)
            fy: 受拉区钢筋抗拉强度(N/mm^2)
            fy_: 受压区钢筋抗拉强度(N/mm^2)
            Es: 受拉区钢筋弹模(N/mm^2)
        """
        self.b = b
        self.h = h
        self.h0 = h0
        self.as_ = as_
        self.As = As
        self.As_ = As_
        self.fc = fc
        self.ft = ft
        self.fy = fy
        self.fy_ = fy_
        self.Es = Es
        self.solve()
    
    def solve(self):
        """
        计算
            As2: 受拉区钢筋抵消受压区钢筋所需面积
            As1: As - As2
            Mu_: 受压区钢筋贡献抗弯弯矩
                fcu: 混凝土强度等级，即立方体抗压强度(N/mm^2)
                alpha1, beta1: 修正系数
                varepsilon_cu: 混凝土极限抗压应变
                xi_b: 界限受压区相对高度
            x_todo: As1下受压区高度[待校核](mm)
            超筋:
                xi: 受压区相对高度
                x: 受压区高度(mm)
                Mu1: As1或混凝土贡献抗弯弯矩(kN*m)
                Mu: 极限抗弯弯矩(kN*m)
            受压区钢筋不能屈服:
                x: 受压区高度[近似取2*as_](mm)
                Mu: 极限抗弯弯矩(kN*m)
            适筋:
                x: [等于x_todo](mm)
                Mu1: As1或混凝土贡献抗弯弯矩(kN*m)
                Mu: 极限抗弯弯矩(kN*m)
        """
        # 将受拉区钢筋分为两部分
        self.As2 = self.As_ * self.fy_ / self.fy
        self.As1 = self.As - self.As2
        self.Mu_ = self.fy_ * self.As_ * (self.h0 - self.as_) * 1e-6
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
        # As1下受压区高度
        self.x_todo = (self.fy * self.As1) / (self.alpha1 * self.fc * self.b)
        # 超筋
        if self.x_todo > (self.xi_b * self.h0):
            temp = (self.xi_b - 0.8) * (self.alpha1 * self.fc * self.b * self.h0) / (self.fy * self.As1)
            self.xi = 0.8 / (1 - temp)
            self.x = self.xi * self.h0
            self.Mu1 = self.alpha1 * self.fc * self.b * self.x * (self.h0 - self.x / 2) * 1e-6
            self.Mu = self.Mu1 + self.Mu_
        # 受压区钢筋不能屈服
        elif self.x_todo < (2 * self.as_):
            # 近似取 self.x = 2 * self.as_
            self.x = 2 * self.as_
            self.Mu = self.fy * self.As * (self.h0 - self.as_) * 1e-6
        # 适筋
        else:
            self.x = self.x_todo
            self.Mu1 = self.fy * self.As1 * (self.h0 - self.x / 2) * 1e-6
            self.Mu = self.Mu1 + self.Mu_

    def process(self):
        """
        以latex公式形式输出计算过程
        """
        process = f"$$\\begin{{aligned}}"
        process += f"&A_{{s2}} = A_s^{{'}} \\frac{{f_y^{{'}}}}{{f_y}} = {self.As2:.5g}mm^2,\\space A_{{s1}} = A_s - A_{{s2}} = {self.As1:.5g}mm^2 \\\\"
        process += f"&M_u^{{'}} = f_y^{{'}} A_s^{{'}} (h_0 - a_s^{{'}}) = {self.Mu_:.2f}kN*m \\\\"
        if self.fcu > 50:
            process += f"&混凝土强度大于50MPa, 取 \\alpha_1 = {self.alpha1}, \\beta_1 = {self.beta1}, \\varepsilon_{{cu}} = 0.0033 - (f_cu - 50) * 10^{{-5}} = {self.varepsilon_cu}, \\\\"
        process += f"&x = \\frac{{f_y A_{{s1}}}}{{\\alpha_1 f_c b}} = {self.x_todo:.2f}mm \\\\"
        process += f"&\\xi_b = \\frac{{\\beta_1}}{{1+\\frac{{f_y}}{{E_s \\varepsilon_{{cu}}}}}} = {self.xi_b:.4g} \\\\"
        if self.x_todo > (self.xi_b * self.h0):
            process += f"&x > \\xi_b h_0,\\space 为超筋,\\space 所以取 \\\\"
            process += f"&\\sigma_s = f_y \\frac{{\\xi - 0.8}}{{\\xi_b - 0.8}}; 联立\\space \\alpha_1 f_c b h_0 \\xi = \\sigma_s A_s \\space 解得：\\\\"
            process += f"&\\xi = {self.xi:.4g},\\space x = \\xi h_0 = {self.x:.2f}mm \\\\"
            process += f"&故\\space M_{{u1}} = \\alpha_1 f_c b x (h_0 - 0.5 x) = {self.Mu1:.2f}kN*m \\\\"
            process += f"&M_u = M_{{u1}} + M_u^{{'}} = {self.Mu:.2f}kN*m"
        elif self.x_todo < (2 * self.as_):
            process += f"&由于 x < 2 a_s^{{'}}, 受压区钢筋不能屈服, 近似取 x = 2 a_s^{{'}}, 则 \\\\"
            process += f"&M_u = f_y A_s (h_0 - a_s^{{'}}) = {self.Mu:.2f}kN*m"
        else:
            process += f"&2 a_s^{{'}} < x < \\xi_b h_0,\\space 为适筋 \\\\"
            process += f"&M_{{u1}} = f_y A_{{s1}} (h_0 - \\frac{{x}}{{2}}) = {self.Mu1:.2f}kN*m \\\\"
            process += f"&M_u = M_{{u1}} + M_u^{{'}} = {self.Mu:.2f}kN*m"
        process += f"\\end{{aligned}}$$"
        return process
