class Calc0:
    def __init__(self, b, h, h0, bf_, hf_, As, fc, ft, fy, Ec, Es):
        """
        T形截面受弯构件正截面承载力计算校核
            b, h, bf_, hf_: 截面尺寸(mm)
                bf_的取值请参阅相关规范, 本程序不进行校核
            h0: 有效高度(mm)
                设 d 为钢筋直径, c 为混凝土保护层厚度,对梁一般为25mm,对板一般为15mm
                采用单排配筋时，取 h0 = h - c - d/2; 双排配筋时，取 h0 = h - c - d - max(12.5, d/2); 
            As: 配筋面积(mm^2)
            fc: 混凝土抗压强度(N/mm^2)
            ft: 混凝土抗拉强度(N/mm^2)
            fy: 钢筋抗拉强度(N/mm^2)
            Ec: 混凝土抗拉弹模(N/mm^2) [少筋时需要]
            Es: 钢筋弹模(N/mm^2)
        """
        self.b = b
        self.h = h
        self.h0 = h0
        self.bf_ = bf_
        self.hf_ = hf_
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
                fcu: 混凝土强度等级，即立方体抗压强度(N/mm^2)
                alpha1, beta1: 修正系数
                varepsilon_cu: 混凝土极限抗压应变
                xi_b: 界限受压区相对高度
            第一类T形截面:
                    rho_min: 最小配筋率
                少筋:
                    alpha_E
                    alpha_A
                    Mu
                适筋:
                    x
                    Mu
            第二类T形截面:
                Muf_
                x_todo
                超筋:
                    xi
                    x
                    Mu1
                    Mu
                适筋:
                    x
                    Mu1
                    Mu
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
        # 第一类T形截面
        if (self.fy * self.As) <= (self.alpha1 * self.fc * self.bf_ * self.hf_):
                # 最小配筋率
            self.rho_min = 0.45 * self.ft / self.fy
            # 少筋
            if self.As < (self.rho_min * self.b * self.h):
                # 按 b, h 矩形截面开裂弯矩计算
                self.alpha_E = self.Es / self.Ec
                self.alpha_A = 2 * self.alpha_E * self.As / (self.b * self.h)
                self.Mu = self.Mcr = 0.292 * (1 + 2.5 * self.alpha_A) * self.ft * self.b * self.h ** 2 * 1e-6
            # 适筋
            else:
                # 按 bf_, h 单筋矩形截面计算
                self.x = (self.fy * self.As) / (self.alpha1 * self.fc * self.bf_)
                self.Mu = self.fy * self.As * (self.h0 - self.x / 2) * 1e-6
        # 第二类T形截面
        else:
            self.Muf_ = self.alpha1 * self.fc * (self.bf_ - self.b) * self.hf_ * (self.h0 - self.hf_ / 2) * 1e-6
            self.x_todo = (self.fy * self.As - self.alpha1 * self.fc * (self.bf_ - self.b) * self.hf_) / (self.alpha1 * self.fc * self.b)
            # 超筋
            if self.x_todo > (self.xi_b * self.h0):
                temp1 = 0.8 + (self.xi_b - 0.8) * (self.alpha1 * self.fc * (self.bf_ - self.b) * self.hf_) / (self.fy * self.As)
                temp2 = 1 - (self.xi_b - 0.8) * (self.alpha1 * self.fc * self.b * self.h0) / (self.fy * self.As)
                self.xi = temp1 / temp2
                self.x = self.xi * self.h0
                self.Mu1 = self.alpha1 * self.fc * self.b * self.x * (self.h0 - self.x / 2) * 1e-6
                self.Mu = self.Mu1 + self.Muf_
            # 适筋
            else:
                self.x = self.x_todo
                self.Mu1 = self.alpha1 * self.fc * self.b * self.x * (self.h0 - self.x / 2) * 1e-6
                self.Mu = self.Mu1 + self.Muf_
    
    def process(self):
        """
        以latex公式形式输出计算过程
        """
        process = f"$$\\begin{{aligned}}"
        # process += f"&b_f^{{'}} < l_0 / 3, b_f^{{'}} \\leq b + 12 h_f^{{'}} = {(self.b + 12 * self.hf_):.2f}mm \\\\"
        if self.fcu > 50:
            process += f"&混凝土强度大于50MPa, 取 \\alpha_1 = {self.alpha1}, \\beta_1 = {self.beta1}, \\varepsilon_{{cu}} = 0.0033 - (f_cu - 50) * 10^{{-5}} = {self.varepsilon_cu}, \\\\"
        if (self.fy * self.As) <= (self.alpha1 * self.fc * self.bf_ * self.hf_):
            process += f"&因为 f_y A_s = {(self.fy * self.As * 1e-3):.2f}kN \\leq \\alpha_1 f_c b_f^{{'}} h_f^{{'}} = {(self.alpha1 * self.fc * self.bf_ * self.hf_ * 1e-3):.2f}kN, 为第一类T形截面。\\\\"
            if self.As < (self.rho_min * self.b * self.h):
                process += f"&A_s < \\rho_{{min}} b h = 0.45 \\frac{{f_t}}{{f_y}} b h = {(self.rho_min * self.b * self.h):.2f}mm^2, 少筋，\\alpha_A = \\frac{{2 E_s A_s}}{{E_c b h}} = {self.alpha_A:.2f} \\\\"
                process += f"&M_u = M_{{cr}} = 0.292 (1 + 2.5 \\alpha_A) f_t b h^2 = {self.Mu:.2f}kN*m \\\\"
            else:
                process += f"&A_s \\geq \\rho_{{min}} b h = 0.45 \\frac{{f_t}}{{f_y}} b h = {(self.rho_min * self.b * self.h):.2f}mm^2, 适筋 \\\\"
                process += f"&x = \\frac{{f_y A_s}}{{\\alpha_1 f_c b_f^{{'}}}} = {self.x:.2f}mm \\\\"
                process += f"&M_u = f_y A_s (h_0 - \\frac{{x}}{{2}}) = {self.Mu:.2f}kN*m \\\\"
        else:
            process += f"&因为 f_y A_s = {(self.fy * self.As * 1e-3):.2f}kN > \\alpha_1 f_c b_f^{{'}} h_f^{{'}} = {(self.alpha1 * self.fc * self.bf_ * self.hf_ * 1e-3):.2f}kN, 为第二类T形截面。\\\\"
            process += f"&M_{{uf}}^{{'}} = \\alpha_1 f_c (b_f^{{'}} - b) h_f^{{'}} (h_0 - \\frac{{h_f^{{'}}}}{{2}}) = {self.Muf_:.2f}kN*m \\\\"
            process += f"&x = \\frac{{f_y A_s - \\alpha_1 f_c (b_f^{{'}} - b) h_f^{{'}}}}{{\\alpha_1 f_c b}} = {self.x_todo:.2f}mm \\\\"
            process += f"&\\xi_b = \\frac{{\\beta_1}}{{1+\\frac{{f_y}}{{E_s \\varepsilon_{{cu}}}}}} = \\frac{{{self.beta1}}}{{1+\\frac{{{self.fy}}}{{{self.Es} * {self.varepsilon_cu}}}}} = {self.xi_b:.4g} \\\\"
            if self.x_todo > (self.xi_b * self.h0):
                process += f"&x > \\xi_b h_0, 超筋, 取 \\sigma_s = f_y \\frac{{\\xi - 0.8}}{{\\xi_b - 0.8}}, 联立 \\alpha_1 f_c b h0 \\xi + \\alpha_1 f_c (b_f^{{'}} - b) h_f^{{'}} = \\sigma_s A_s 解得 \\\\"
                process += f"&\\xi = {self.xi:.4g}, x = \\xi h_0 = {self.x:.2f}mm, M_{{u1}} = \\alpha_1 f_c b x (h_0 - \\frac{{x}}{{2}}) = {self.Mu1:.2f}kN*m \\\\"
                process += f"&M_u = M_{{u1}} + M_{{uf}}^{{'}} = {self.Mu:.2f}kN*m \\\\"
            else:
                process += f"&x \\leq \\xi_b h_0, M_{{u1}} = \\alpha_1 f_c b x (h_0 - \\frac{{x}}{{2}}) = {self.Mu1:.2f}kN*m \\\\"
                process += f"&M_u = M_{{u1}} + M_{{uf}}^{{'}} = {self.Mu:.2f}kN*m \\\\"
        process += f"\\end{{aligned}}$$"
        return process
