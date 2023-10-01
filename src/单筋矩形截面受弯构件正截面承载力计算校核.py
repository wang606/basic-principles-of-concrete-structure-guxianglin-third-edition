class Calc6:
    def __init__(self, b, h, h0, As, fc, ft, fy, Ec, Es):
        """
        单筋矩形截面受弯构件正截面承载力计算校核
            b, h: 截面尺寸(mm)
            h0: 有效高度(mm)
                设 d 为钢筋直径, c 为混凝土保护层厚度,对梁一般为25mm,对板一般为15mm
                采用单排配筋时，取 h0 = h - c - d/2; 双排配筋时，取 h0 = h - c - d - max(12.5, d/2); 
            As: 配筋面积(mm^2)
            fc: 混凝土抗压强度(N/mm^2)
            ft: 混凝土抗拉强度(N/mm^2)
            fy: 钢筋抗拉强度(N/mm^2)
            Ec: 混凝土抗拉弹模(N/mm^2)
            Es: 钢筋弹模(N/mm^2)
        """
        self.b = b
        self.h = h
        self.h0 = h0
        self.As = As
        self.fc = fc
        self.ft = ft
        self.fy = fy
        self.Ec = Ec
        self.Es = Es
        self.solve()
    
    def solve(self):
        """
        计算
            rho_min: 最小配筋率
                fcu: 混凝土强度等级，即立方体抗压强度(N/mm^2)
                alpha1, beta1: 修正系数
                varepsilon_cu: 混凝土极限抗压应变
                xi_b: 界限受压区相对高度
            rho_max: 最大配筋率
            rho: 配筋率
            Mcr: 开裂弯矩(kN*m)
            少筋: 
                Mu: 极限抗弯弯矩(kN*m)
            超筋:
                xi: 受压区相对高度
                sigma_s: 钢筋应力(N/mm^2)
                x: 受压区高度(mm)
                Mu: 极限抗弯弯矩(kN*m)
            适筋:
                x: 受压区高度(mm)
                Mu: 极限抗弯弯矩(kN*m)
        """
        # 最小配筋率
        self.rho_min = 0.45 * self.ft / self.fy
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
        # 最大配筋率
        self.rho_max = self.xi_b * self.alpha1 * self.fc / self.fy
        # 配筋率
        self.rho = self.As / (self.b * self.h)

        # 开裂弯矩
        self.alpha_E = self.Es / self.Ec
        self.alpha_A = 2 * self.alpha_E * self.As / (self.b * self.h)
        self.Mcr = 0.292 * (1 + 2.5 * self.alpha_A) * self.ft * self.b * self.h ** 2 * 1e-6
        # 如果配筋率小于最小配筋率，少筋
        if self.rho < self.rho_min:
            self.Mu = self.Mcr
        # 如果配筋率大于最大配筋率，超筋
        elif self.rho > self.rho_max:
            temp = (self.xi_b - 0.8) * (self.alpha1 * self.fc * self.b * self.h0) / (self.fy * self.As)
            self.xi = 0.8 / (1 - temp)
            self.sigma_s = self.fy * (self.xi - 0.8) / (self.xi_b - 0.8)
            self.x = self.xi * self.h0
            self.Mu = self.alpha1 * self.fc * self.b * self.x * (self.h0 - self.x / 2) * 1e-6
        # 适筋
        else:
            self.x = (self.fy * self.As) / (self.alpha1 * self.fc * self.b)
            self.Mu = self.fy * self.As * (self.h0 - self.x / 2) * 1e-6

    def process(self):
        """
        以latex公式形式输出计算过程
        """
        process = f"$$\\begin{{aligned}}"
        process += f"&最小配筋率 \\rho_{{min}} = 0.45\\frac{{f_t}}{{f_y}} = 0.45 * \\frac{{{self.ft}}}{{{self.fy}}} = {self.rho_min:.4g} \\\\"
        if self.fcu > 50:
            process += f"&混凝土强度大于50MPa, 取 \\alpha_1 = {self.alpha1}, \\beta_1 = {self.beta1}, \\varepsilon_{{cu}} = 0.0033 - (f_cu - 50) * 10^{{-5}} = {self.varepsilon_cu}, \\\\"
        process += f"&界限受压区相对高度 \\xi_b = \\frac{{\\beta_1}}{{1+\\frac{{f_y}}{{E_s \\varepsilon_{{cu}}}}}} = \\frac{{{self.beta1}}}{{1+\\frac{{{self.fy}}}{{{self.Es} * {self.varepsilon_cu}}}}} = {self.xi_b:.4g} \\\\"
        process += f"&最大配筋率 \\rho_{{max}} = \\xi_b \\frac{{\\alpha_1 f_c}}{{f_y}} = {self.xi_b:.4g} * \\frac{{{self.alpha1} * {self.fc}}}{{{self.fy}}} = {self.rho_max:.4g} \\\\"
        process += f"&配筋率 \\rho = \\frac{{A_s}}{{b h}} = \\frac{{{self.As}}}{{{self.b} * {self.h}}} = {self.rho:.4g} \\\\"
        if self.rho < self.rho_min:
            process += f"&因为\\space \\rho < \\rho_{{min}},\\space 为少筋, \\\\"
            process += f"&\\alpha_E = \\frac{{E_s}}{{E_c}} = {self.alpha_E:.4g}, \\alpha_A = \\frac{{2 \\alpha_E A_s}}{{b h}} = {self.alpha_A:.4g}, \\\\"
            process += f"&故\\space M_u = M_{{cr}} = 0.292(1+2.5\\alpha_A)f_tbh^2 = {self.Mcr:.2f}kN*m"
        elif self.rho > self.rho_max:
            process += f"&因为\\space \\rho > \\rho_{{max}},\\space 为超筋, 所以取 \\\\"
            process += f"&\\sigma_s = f_y \\frac{{\\xi - 0.8}}{{\\xi_b - 0.8}}; 联立\\space \\alpha_1 f_c b h_0 \\xi = \\sigma_s A_s \\space 解得：\\\\"
            process += f"&\\xi = {self.xi:.4g},\\space \\sigma_s = {self.sigma_s:.2f}N/mm^2 \\\\"
            process += f"&故\\space M_u = \\alpha_1 f_c b h_0^2 \\xi (1 - 0.5 \\xi) = {self.Mu:.2f}kN*m"
        else:
            process += f"&\\rho_{{min}} < \\rho < \\rho_{{max}},\\space 为适筋, \\\\"
            process += f"&由\\space \\alpha_1 f_c b x = f_y A_s \\space 解得\\space x = {self.x:.2f}mm \\\\"
            process += f"&故\\space M_u = f_y A_s (h_0 - \\frac{{x}}{{2}}) = {self.Mu:.2f}kN*m"
        process += f"\\end{{aligned}}$$"
        return process
