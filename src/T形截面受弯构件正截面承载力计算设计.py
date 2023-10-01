from math import sqrt


class Calc1:
    def __init__(self, b, h, h0, bf_, hf_, fc, ft, fy, Es, M):
        """
        T形截面受弯构件正截面承载力计算设计
            b, h, bf_, hf_: 截面尺寸(mm)
                bf_的取值请参阅相关规范, 本程序不进行校核
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
        self.bf_ = bf_
        self.hf_ = hf_
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
            第一类T形截面:
                x: 受压区截面高度
                As_todo: 待校核最小配筋率的配筋面积
                rho_min
                少筋:
                    As: = rho_min * b * h
                适筋:
                    As: = As_todo
            第二类T形截面:
                As2
                M_
                M1
                x
                超筋:
                    重新设计截面或在受压区配置钢筋
                适筋:
                    As1
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
        # 第一类T形截面
        if self.M <= (self.alpha1 * self.fc * self.bf_ * self.hf_ * (self.h0 - self.hf_ / 2) * 1e-6):
            # 受压区截面高度
            self.x = self.h0 - sqrt(self.h0 ** 2 - 2 * self.M * 1e6 / (self.alpha1 * self.fc * self.bf_))
            # 配筋面积
            self.As_todo = self.alpha1 * self.fc * self.bf_ * self.x / self.fy
            # 最小配筋率
            self.rho_min = 0.45 * self.ft / self.fy
            if self.As_todo < (self.rho_min * self.b * self.h):
                self.As = self.rho_min * self.b * self.h
            else:
                self.As = self.As_todo
        # 第二类T形截面
        else:
            self.As2 = self.alpha1 * self.fc * (self.bf_ - self.b) * self.hf_ / self.fy
            self.M_ = self.As2 * self.fy * (self.h0 - self.hf_ / 2) * 1e-6
            self.M1 = self.M - self.M_
            self.x = self.h0 - sqrt(self.h0 ** 2 - 2 * self.M1 * 1e6 / (self.alpha1 * self.fc * self.b))
            if self.x > (self.xi_b * self.h0):
                # 超筋 需重新设计截面或在受压区布置钢筋
                return
            else:
                # 适筋
                self.As1 = self.alpha1 * self.fc * self.b * self.x / self.fy
                self.As = self.As1 + self.As2
    
    def process(self):
        """
        以latex公式形式输出计算过程
        """
        process = f"$$\\begin{{aligned}}"
        if self.fcu > 50:
            process += f"&混凝土强度大于50MPa, 取 \\alpha_1 = {self.alpha1}, \\beta_1 = {self.beta1}, \\varepsilon_{{cu}} = 0.0033 - (f_cu - 50) * 10^{{-5}} = {self.varepsilon_cu}, \\\\"
        if self.M <= (self.alpha1 * self.fc * self.bf_ * self.hf_ * (self.h0 - self.hf_ / 2) * 1e-6):
            process += f"&M \\leq \\alpha_1 f_c b_f^{{'}} h_f^{{'}} (h_0 - \\frac{{h_f^{{'}}}}{{2}}) = {(self.alpha1 * self.fc * self.bf_ * self.hf_ * (self.h0 - self.hf_ / 2) * 1e-6):.2f}kN*m, 为第一类T形截面 \\\\"
            process += f"&由 M = M_u = \\alpha_1 f_c b_f^{{'}} x (h_0 - \\frac{{x}}{{2}}) 解得 x = {self.x:.2f}mm \\\\"
            process += f"&此时配筋面积 A_s = \\alpha_1 f_c b_f^{{'}} x / f_y = {self.As_todo:.2f}mm^2 \\\\"
            if self.As_todo < (self.rho_min * self.b * self.h):
                process += f"&由于 A_s < \\rho_{{min}} b h = 0.45 \\frac{{f_t}}{{f_y}} b h = {(self.rho_min * self.b * self.h):.2f}mm^2 \\\\"
                process += f"&取 A_s = \\rho_{{min}} b h = {self.As:.2f}mm^2 \\\\"
            else:
                process += f"&A_s \\geq \\rho_{{min}} b h = 0.45 \\frac{{f_t}}{{f_y}} b h = {(self.rho_min * self.b * self.h):.2f}mm^2, 校核通过。 \\\\"
        else:
            process += f"&M > \\alpha_1 f_c b_f^{{'}} h_f^{{'}} (h_0 - \\frac{{h_f^{{'}}}}{{2}}) = {(self.alpha1 * self.fc * self.bf_ * self.hf_ * (self.h0 - self.hf_ / 2) * 1e-6):.2f}kN*m, 为第二类T形截面 \\\\"
            process += f"&A_{{s2}} = \\alpha_1 f_c (b_f^{{'}} - b) h_f^{{'}} / f_y = {self.As2:.2f}mm^2 \\\\"
            process += f"&M^{{'}} = A_{{s2}} f_y (h_0 - \\frac{{h_f^{{'}}}}{{2}}) = {self.M_:.2f}kN*m, M_1 = M - M^{{'}} = {self.M1:.2f}kN*m \\\\"
            process += f"&由 M_1 = M_{{u1}} = \\alpha_1 f_c b x (h_0 - \\frac{{x}}{{2}}) 解得: x = {self.x:.2f}mm \\\\"
            process += f"&\\xi_b = \\frac{{\\beta_1}}{{1+\\frac{{f_y}}{{E_s \\varepsilon_{{cu}}}}}} = \\frac{{{self.beta1}}}{{1+\\frac{{{self.fy}}}{{{self.Es} * {self.varepsilon_cu}}}}} = {self.xi_b:.4g} \\\\"
            if self.x > (self.xi_b * self.h0):    
                process += f"&x > \\xi_b h_0, 为超筋, 需重新设计截面或考虑在受压区布置钢筋！\\\\"
            else:
                process += f"&x \\leq \\xi_b h_0, A_{{s1}} = \\alpha_1 f_c b x / f_y = {self.As1:.2f}mm^2 \\\\"
                process += f"&A_s = A_{{s1}} + A_{{s2}} = {self.As:.2f}mm^2 \\\\"
        process += f"\\end{{aligned}}$$"
        return process
