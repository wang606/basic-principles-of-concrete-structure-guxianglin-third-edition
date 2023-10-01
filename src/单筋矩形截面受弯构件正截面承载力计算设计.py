from math import sqrt


class Calc7:
    def __init__(self, b, h, h0, fc, ft, fy, Es, M):
        """
        单筋矩形截面受弯构件正截面承载力计算设计
            b, h: 截面尺寸(mm)
            h0: 有效高度(mm)
                对钢筋混凝土梁，采用单排配筋时，取 h0 = h - 35; 双排配筋时，取 h0 = h - 60; 
                对钢筋混凝土板，取 h0 = h - 20; 
            fc: 混凝土抗压强度(N/mm^2)
            ft: 混凝土抗拉强度(N/mm^2)
            fy: 钢筋抗拉强度(N/mm^2)
            Es: 钢筋弹模(N/mm^2)
            M: 设计弯矩(kN*m)
        """
        self.b = b
        self.h = h
        self.h0 = h0
        self.fc = fc
        self.ft = ft
        self.fy = fy
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
            Mu_max: 适筋构件最大抗弯承载力(kN*m)
            x: 受压区截面高度(mm)
            As_todo: 待校核配筋率的配筋面积(mm^2)
            rho: 配筋率
            rho_min: 最小配筋率
            As: 最终配筋面积(mm^2)
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
        # 适筋构件最大抗弯承载力
        self.Mu_max = self.alpha1 * self.fc * self.b * self.h0 ** 2 * self.xi_b * (1 - 0.5 * self.xi_b) * 1e-6
        if self.M > self.Mu_max:
            # 适筋构件最大抗弯承载力小于设计弯矩，需增大截面面积
            return 
        # 受压区截面高度
        self.x = self.h0 - sqrt(self.h0 ** 2 - 2 * self.M * 1e6 / (self.alpha1 * self.fc * self.b))
        # 配筋面积
        self.As_todo = self.alpha1 * self.fc * self.b * self.x / self.fy
        # 配筋率
        self.rho = self.As_todo / (self.b * self.h)
        # 最小配筋率
        self.rho_min = 0.45 * self.ft / self.fy
        if self.rho < self.rho_min:
            self.As = self.rho_min * self.b * self.h
        else:
            self.As = self.As_todo
    
    def process(self):
        """
        以latex公式形式输出计算过程
        """
        process = f"$$\\begin{{aligned}}"
        if self.fcu > 50:
            process += f"&混凝土强度大于50MPa, 取 \\alpha_1 = {self.alpha1}, \\beta_1 = {self.beta1}, \\varepsilon_{{cu}} = 0.0033 - (f_cu - 50) * 10^{{-5}} = {self.varepsilon_cu}, \\\\"
        process += f"&界限受压区相对高度 \\xi_b = \\frac{{\\beta_1}}{{1+\\frac{{f_y}}{{E_s \\varepsilon_{{cu}}}}}} = \\frac{{{self.beta1}}}{{1+\\frac{{{self.fy}}}{{{self.Es} * {self.varepsilon_cu}}}}} = {self.xi_b:.4g} \\\\"
        process += f"&适筋构件最大抗弯承载力\\space M_{{u,max}} = \\alpha_1 f_c b h_0^2 \\xi_b (1 - 0.5 \\xi_b) = {self.Mu_max:.2f}kN*m \\\\"
        if self.M > self.Mu_max:
            process += f"&M_{{u,max}} < {self.M:.2f}kN*m, 需加大截面尺寸！\\\\"
            process += f"\\end{{aligned}}$$"
            return process
        process += f"&由\\space M = M_u = \\alpha_1 f_c b x (h_0 - \\frac{{x}}{{2}}) \\space 解得:\\space x = {self.x:.2f}mm \\\\"
        process += f"&由\\space \\alpha_1 f_c b x = f_y A_s \\space 解得:\\space A_s = {self.As_todo:.2f}mm^2 \\\\"
        process += f"&此时，配筋率\\space \\rho = \\frac{{A_s}}{{b h}} = {self.rho:.4g} \\\\"
        process += f"&而最小配筋率\\space \\rho_{{min}} = 0.45 \\frac{{f_t}}{{f_y}} = {self.rho_min:.4g} \\\\"
        if self.rho < self.rho_min:
            process += f"&\\rho < \\rho_{{min}},\\space 取\\space A_s = \\rho_{{min}} b h = {self.As:.2f}mm^2"
        else:
            process += f"&\\rho > \\rho_{{min}},\\space 故\\space A_s = {self.As:.2f}mm^2"
        process += f"\\end{{aligned}}$$"
        return process
