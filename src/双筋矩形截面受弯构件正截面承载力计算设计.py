from math import sqrt


class Calc9:
    def __init__(self, b, h, h0, as_, fc, ft, fy, fy_, Es, M):
        """
        双筋矩形截面受弯构件正截面承载力计算设计[受压钢筋截面面积未知]
            b, h: 截面尺寸(mm)
            h0: 有效高度(mm)
                对钢筋混凝土梁，采用单排配筋时，取 h0 = h - 35; 双排配筋时，取 h0 = h - 60; 
                对钢筋混凝土板，取 h0 = h - 20; 
            as_: 受压区钢筋中心线到截面受压区边缘的距离(mm)
                受压区钢筋一般都是单排, 故 as_ 常取 梁:35mm
            fc: 混凝土抗压强度(N/mm^2)
            ft: 混凝土抗拉强度(N/mm^2)
            fy: 受拉区钢筋抗拉强度(N/mm^2)
            fy_: 受压区钢筋抗拉强度(N/mm^2)
            Es: 受拉区钢筋弹模(N/mm^2)
            M: 设计弯矩(kN*m)
        """
        self.b = b
        self.h = h
        self.h0 = h0
        self.as_ = as_
        self.fc = fc
        self.ft = ft
        self.fy = fy
        self.fy_ = fy_
        self.Es = Es
        self.M = M
        self.solve()
    
    def solve(self):
        """
        计算
                fcu: 混凝土强度等级，即立方体抗压强度(N/mm^2)
                alpha1, beta1: 修正系数
                varepsilon_cu: 混凝土极限抗压应变
                xi_b: 界限受压区相对高度
            x: 界限受压区高度
            As1
            M1
            M_
            As2
            As_
            As
        """
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
        # 界限受压区高度
        self.x = self.xi_b * self.h0
        self.As1 = self.alpha1 * self.fc * self.b * self.x / self.fy
        self.M1 = self.As1 * self.fy * (self.h0 - 0.5 * self.x) * 1e-6
        self.M_ = self.M - self.M1
        self.As2 = self.M_ * 1e6 / ((self.h0 - self.as_) * self.fy)
        self.As_ = self.As2 * self.fy / self.fy_
        self.As = self.As1 + self.As2
    
    def process(self):
        """
        以latex公式形式输出计算过程
        """
        process = f"$$\\begin{{aligned}}"
        if self.fcu > 50:
            process += f"&混凝土强度大于50MPa, 取 \\alpha_1 = {self.alpha1}, \\beta_1 = {self.beta1}, \\varepsilon_{{cu}} = 0.0033 - (f_cu - 50) * 10^{{-5}} = {self.varepsilon_cu}, \\\\"
        process += f"&\\xi_b = \\frac{{\\beta_1}}{{1+\\frac{{f_y}}{{E_s \\varepsilon_{{cu}}}}}} = {self.xi_b:.4g} \\\\"
        process += f"&取\\space x = \\xi_b h_0 = {self.x:.2f}mm \\\\"
        process += f"&A_{{s1}} = \\alpha_1 f_c b x / f_y = {self.As1:.2f}mm^2,\\space M_1 = A_{{s1}} f_y (h_0 - 0.5 x) = {self.M1:.2f}kN*m \\\\"
        process += f"&M^{{'}} = M - M_1 = {self.M_:.2f}kN*m,\\space A_{{s2}} = \\frac{{M^{{'}}}}{{(h_0 - a_s^{{'}}) f_y}} = {self.As2:.2f}mm^2 \\\\"
        process += f"&A_s^{{'}} = A_{{s2}} \\frac{{f_y}}{{f_y^{{'}}}} = {self.As_:.2f}mm^2,\\space A_s = A_{{s1}} + A_{{s2}} = {self.As:.2f}mm^2"
        process += f"\\end{{aligned}}$$"
        return process


class Calc10:
    def __init__(self, b, h, h0, as_, As_, fc, ft, fy, fy_, Es, M):
        """
        双筋矩形截面受弯构件正截面承载力计算设计[受压钢筋截面面积已知]
            b, h: 截面尺寸(mm)
            h0: 有效高度(mm)
                对钢筋混凝土梁，采用单排配筋时，取 h0 = h - 35; 双排配筋时，取 h0 = h - 60; 
                对钢筋混凝土板，取 h0 = h - 20; 
            as_: 受压区钢筋中心线到截面受压区边缘的距离(mm)
                设受压区钢筋直径为 d, 混凝土保护层厚度为 c(梁为25mm)
                由于受压区钢筋一般为单排，故 as_ = c + d/2
            As_: 受压钢筋截面面积(mm^2)
            fc: 混凝土抗压强度(N/mm^2)
            ft: 混凝土抗拉强度(N/mm^2)
            fy: 受拉区钢筋抗拉强度(N/mm^2)
            fy_: 受压区钢筋抗拉强度(N/mm^2)
            Es: 受拉区钢筋弹模(N/mm^2)
            M: 设计弯矩(kN*m)
        """
        self.b = b
        self.h = h
        self.h0 = h0
        self.as_ = as_
        self.As_ = As_
        self.fc = fc
        self.ft = ft
        self.fy = fy
        self.fy_ = fy_
        self.Es = Es
        self.M = M
        self.solve()
    
    def solve(self):
        """
        计算
            As2
            M_
            M1
                fcu
                alpha1, beta1
                varepsilon_cu
                xi_b
            Mu1_max
            x
            超筋:
                subCalc: 按As_未知从新计算
            适筋:
                As1
                As
                受压钢筋不能屈服:
                    rho
                    rho_min
                    少筋:
                        As
        """
        self.As2 = self.As_ * self.fy_ / self.fy
        self.M_ = self.As2 * self.fy * (self.h0 - self.as_) * 1e-6
        self.M1 = self.M - self.M_
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
        """
        # 适筋构件最大抗弯承载力
        self.Mu1_max = self.alpha1 * self.fc * self.b * self.h0 ** 2 * self.xi_b * (1 - 0.5 * self.xi_b) * 1e-6
        if self.M1 > self.Mu1_max:
            # 适筋构件最大抗弯承载力小于设计弯矩，需增大截面面积
            return 
        """
        # 受压区截面高度
        self.x = self.h0 - sqrt(self.h0 ** 2 - 2 * self.M1 * 1e6 / (self.alpha1 * self.fc * self.b))
        # 超筋
        if self.x > (self.xi_b * self.h0):
            # 此时弃用 As_, 按 As_ 未知来计算, as_ 取 35mm
            self.subCalc = Calc9(self.b, self.h, self.h0, 35, self.fc, self.ft, self.fy, self.fy_, self.Es, self.M)
        else:
        # 适筋
            self.As1 = self.alpha1 * self.fc * self.b * self.x / self.fy
            self.As = self.As1 + self.As2
        # 受压区钢筋不能屈服，需要做最小配筋率校核
            if self.x < (2 * self.as_):
                self.rho = self.As / (self.b * self.h)
                self.rho_min = 0.45 * self.ft / self.fy
                if self.rho < self.rho_min:
                    self.As = self.rho_min * self.b * self.h
    
    def process(self):
        """
        以latex公式形式输出计算过程
        """
        process = f"$$\\begin{{aligned}}"
        process += f"&A_{{s2}} = A_s^{{'}} \\frac{{f_y^{{'}}}}{{f_y}} = {self.As2:.2f}mm^2 \\\\"
        process += f"&M^{{'}} = A_{{s2}} f_y (h_0 - a_s^{{'}}) = {self.M_:.2f}kN*m \\\\"
        process += f"&M_1 = M - M^{{'}} = {self.M1:.2f}kN*m \\\\"
        if self.fcu > 50:
            process += f"&混凝土强度大于50MPa, 取 \\alpha_1 = {self.alpha1}, \\beta_1 = {self.beta1}, \\varepsilon_{{cu}} = 0.0033 - (f_cu - 50) * 10^{{-5}} = {self.varepsilon_cu}, \\\\"
        process += f"&\\xi_b = \\frac{{\\beta_1}}{{1+\\frac{{f_y}}{{E_s \\varepsilon_{{cu}}}}}} = {self.xi_b:.4g} \\\\"
        """
        process += f"&适筋构件最大抗弯承载力\\space M_{{u1,max}} = \\alpha_1 f_c b h_0^2 \\xi_b (1 - 0.5 \\xi_b) = {self.Mu1_max:.2f}kN*m \\\\"
        if self.M1 > self.Mu1_max:
            process += f"&M_{{u1,max}} < M_1, 需加大截面尺寸！\\\\"
            process += f"\\end{{aligned}}$$"
            return process
        """
        process += f"&由\\space M_1 = M_{{u1}} = \\alpha_1 f_c b x (h_0 - \\frac{{x}}{{2}}) \\space 解得:\\space x = {self.x:.2f}mm \\\\"
        if self.x > (self.xi_b * self.h0):
            process += f"&x > \\xi_b h_0,需按 A_s^{{'}} 未知从新计算，取 a_s^{{'}} = 35mm, 如下：\\\\"
            process += self.subCalc.process().strip("$$\\begin{{aligned}}").strip("\\end{{aligned}}$$")
        elif self.x < (2 * self.as_):
            process += f"&x < 2 a_s^{{'}}, 需校核最小配筋率 \\\\"
            process += f"&A_{{s1}} = \\alpha_1 f_c b x / f_y = {self.As1:.2f}mm^2 \\\\"
            process += f"&\\rho = \\frac{{A_{{s1}} + A_{{s2}}}}{{b h}} = {self.rho:.4g},\\space 最小配筋率 \\rho_{{min}} = 0.45 \\frac{{f_t}}{{f_y}} = {self.rho_min:.4g} \\\\"
            if self.rho < self.rho_min:
                process += f"&\\rho < \\rho_{{min}}, 取\\space A_s = \\rho_{{min}} b h = {self.As:.2f}mm^2"
            else:
                process += f"&\\rho \\geq \\rho_{{min}},\\space A_s = A_{{s1}} + A_{{s2}} = {self.As:.2f}mm^2"
        else:
            process += f"&2 a_s^{{'}} < x < \\xi_b h_0,\\space 为适筋 \\\\"
            process += f"&A_{{s1}} = \\alpha_1 f_c b x / f_y = {self.As1:.2f}mm^2 \\\\"
            process += f"&A_s = A_{{s1}} + A_{{s2}} = {self.As:.2f}mm^2"
        process += f"\\end{{aligned}}$$"
        return process