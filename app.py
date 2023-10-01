from viktor import ViktorController
from viktor.parametrization import ViktorParametrization, NumberField, OptionField, Lookup, Text, LineBreak
from viktor.views import DataView, DataGroup, DataItem, DataResult
from src.T形截面受弯构件正截面承载力计算校核 import Calc0
from src.T形截面受弯构件正截面承载力计算设计 import Calc1
from src.不对称配筋偏心受压构件基于承载力的截面校核 import Calc2, Calc3
from src.不对称配筋偏心受压构件基于承载力的截面设计 import Calc4, Calc5
from src.单筋矩形截面受弯构件正截面承载力计算校核 import Calc6
from src.单筋矩形截面受弯构件正截面承载力计算设计 import Calc7
from src.双筋矩形截面受弯构件正截面承载力计算校核 import Calc8
from src.双筋矩形截面受弯构件正截面承载力计算设计 import Calc9, Calc10
from src.对称配筋偏心受压构件基于承载力的截面设计 import Calc11


popts = [
    "T形截面受弯构件正截面承载力计算校核", 
    "T形截面受弯构件正截面承载力计算设计", 
    "不对称配筋偏心受压构件基于承载力的截面校核(已知初始偏心距求极限轴向抗压承载力)", 
    "不对称配筋偏心受压构件基于承载力的截面校核(已知轴向力求极限抗弯弯矩)", 
    "不对称配筋偏心受压构件基于承载力的截面设计(受压区钢筋截面积未知)", 
    "不对称配筋偏心受压构件基于承载力的截面设计(受压区钢筋截面积已知)", 
    "单筋矩形截面受弯构件正截面承载力计算校核", 
    "单筋矩形截面受弯构件正截面承载力计算设计", 
    "双筋矩形截面受弯构件正截面承载力计算校核", 
    "双筋矩形截面受弯构件正截面承载力计算设计(受压区钢筋截面积未知)", 
    "双筋矩形截面受弯构件正截面承载力计算设计(受压区钢筋截面积已知)", 
    "对称配筋偏心受压构件基于承载力的截面设计"
]


"""
l0      2 3 4 5 11
b       0 1 2 3 4 5 6 7 8 9 10 11
h       0 1 2 3 4 5 6 7 8 9 10 11
h0      0 1 6 7 8 9 10
bf_     0 1
hf_     0 1
As      0 2 3 6 8
As_     2 3 5 8 10
a_s     2 3 4 5 11
a_s_    2 3 4 5 8 9 10
e_0     2
fc      0 1 2 3 4 5 6 7 8 9 10 11
ft      0 1 6 7 8 9 10
fy      0 1 2 3 4 5 6 7 8 9 10 11
fy_     2 3 4 5 8 9 10
Ec      0 6
Es      0 1 2 3 4 5 6 7 8 9 10 11
Nc      3 4 5 11
M       1 4 5 7 9 10 11
"""


l0_visible = lambda params, **kwargs: params.problem_type in [popts[2], popts[3], popts[4], popts[5], popts[11]]
b_visible = lambda params, **kwargs: params.problem_type in popts
h_visible = lambda params, **kwargs: params.problem_type in popts
h0_visible = lambda params, **kwargs: params.problem_type in [popts[0], popts[1], popts[6], popts[7], popts[8], popts[9], popts[10]]
bf__visible = lambda params, **kwargs: params.problem_type in [popts[0], popts[1]]
hf__visible = lambda params, **kwargs: params.problem_type in [popts[0], popts[1]]
As_visible = lambda params, **kwargs: params.problem_type in [popts[0], popts[2], popts[3], popts[6], popts[8]]
As__visible = lambda params, **kwargs: params.problem_type in [popts[2], popts[3], popts[5], popts[8], popts[10]]
a_s_visible = lambda params, **kwargs: params.problem_type in [popts[2], popts[3], popts[4], popts[5], popts[11]]
a_s__visible = lambda params, **kwargs: params.problem_type in [popts[2], popts[3], popts[4], popts[5], popts[8], popts[9], popts[10]]
e_0_visible = lambda params, **kwargs: params.problem_type in [popts[2]]
fc_visible = lambda params, **kwargs: params.problem_type in popts
ft_visible = lambda params, **kwargs: params.problem_type in [popts[0], popts[1], popts[6], popts[7], popts[8], popts[9], popts[10]]
fy_visible = lambda params, **kwargs:  params.problem_type in popts
fy__visible = lambda params, **kwargs: params.problem_type in [popts[2], popts[3], popts[4], popts[5], popts[8], popts[9], popts[10]]
Ec_visible = lambda params, **kwargs: params.problem_type in [popts[0], popts[6]]
Es_visible = lambda params, **kwargs: params.problem_type in popts
Nc_visible = lambda params, **kwargs: params.problem_type in [popts[3], popts[4], popts[5], popts[11]]
M_visible = lambda params, **kwargs: params.problem_type in [popts[1], popts[4], popts[5], popts[7], popts[9], popts[10], popts[11]]


class Parametrization(ViktorParametrization):
    problem_type = OptionField("题目类型", options=popts, default=popts[0], flex=100)
    l0 = NumberField("柱长 $$l_0 [m]$$", min=0, visible=l0_visible)
    b = NumberField("截面宽度 $$b [mm]$$", min=0, visible=b_visible)
    h = NumberField("截面高度 $$h [mm]$$", min=0, visible=h_visible)
    h0 = NumberField("有效高度 $$h_0 [mm]$$", min=0, max=Lookup("h"), visible=h0_visible)
    bf_ = NumberField("T型截面受压翼缘计算宽度 $$b_f^` [mm]$$", min=Lookup("b"), visible=bf__visible)
    hf_ = NumberField("T型截面受压翼缘高度 $$h_f^` [mm]$$", min=0, max=Lookup("h"), visible=hf__visible)
    As = NumberField("(受拉区)钢筋面积 $$A_s [mm^2]$$", min=0, visible=As_visible)
    As_ = NumberField("受压区钢筋面积 $$A_s^` [mm^2]$$", min=0, visible=As__visible)
    a_s = NumberField("(受拉区)钢筋中心到混凝土边缘的距离 $$a_s [mm]$$", min=0, visible=a_s_visible)
    a_s_ = NumberField("受压区钢筋中心到混凝土边缘的距离 $$a_s^` [mm]$$", min=0, visible=a_s__visible)
    e_0 = NumberField("初始偏心距 $$e_0 [mm]$$", min=0, visible=e_0_visible)
    fc = NumberField("混凝土抗压强度 $$f_c [MPa]$$", min=0, visible=fc_visible)
    ft = NumberField("混凝土抗拉强度 $$f_t [MPa]$$", min=0, visible=ft_visible)
    fy = NumberField("(受拉区)钢筋抗拉强度 $$f_y [MPa]$$", min=0, visible=fy_visible)
    fy_ = NumberField("受压区钢筋抗拉强度 $$f_y^` [MPa]$$", min=0, visible=fy__visible)
    Ec = NumberField("混凝土抗拉弹模 $$E_c [MPa]$$", min=0, visible=Ec_visible)
    Es = NumberField("钢筋弹模 $$E_s [MPa]$$", min=0, visible=Es_visible)
    Nc = NumberField("轴向力 $$N_c [kN]$$", visible=Nc_visible)
    M = NumberField("弯矩 $$M [kN·m]$$", visible=M_visible)
    _ = LineBreak()
    ps = Text("made by [wang606](https://github.com/wang606/basic-principles-of-concrete-structure-guxianglin-third-edition). ")


class Controller(ViktorController):
    label = '混凝土结构基本原理 顾祥林（第三版）'
    parametrization = Parametrization

    @DataView("计算书", duration_guess=1)
    def result(self, params, **kwargs):
        ptype = params.problem_type

        l0 = params.l0
        b = params.b
        h = params.h
        h0 = params.h0
        bf_ = params.bf_
        hf_ = params.hf_
        As = params.As
        As_ = params.As_
        a_s = params.a_s
        a_s_ = params.a_s_
        e_0 = params.e_0
        fc = params.fc
        ft = params.ft
        fy = params.fy
        fy_ = params.fy_
        Ec = params.Ec
        Es = params.Es
        Nc = params.Nc
        M = params.M

        label = f"$${ptype}\\\\\\begin{{aligned}}已知:&\\\\"
        label += f"&l_0 = {l0} m\\\\" if l0_visible(params) else ""
        label += f"&b = {b} mm\\\\" if b_visible(params) else ""
        label += f"&h = {h} mm\\\\" if h_visible(params) else ""
        label += f"&h_0 = {h0} mm\\\\" if h0_visible(params) else ""
        label += f"&b_f^` = {bf_} mm\\\\" if bf__visible(params) else ""
        label += f"&h_f^` = {hf_} mm\\\\" if hf__visible(params) else ""
        label += f"&A_s = {As} mm^2\\\\" if As_visible(params) else ""
        label += f"&A_s^` = {As_} mm^2\\\\" if As__visible(params) else ""
        label += f"&a_s = {a_s} mm\\\\" if a_s_visible(params) else ""
        label += f"&a_s^` = {a_s_} mm\\\\" if a_s__visible(params) else ""
        label += f"&e_0 = {e_0} mm\\\\" if e_0_visible(params) else ""
        label += f"&f_c = {fc} MPa\\\\" if fc_visible(params) else ""
        label += f"&f_t = {ft} MPa\\\\" if ft_visible(params) else ""
        label += f"&f_y = {fy} MPa\\\\" if fy_visible(params) else ""
        label += f"&f_y^` = {fy_} MPa\\\\" if fy__visible(params) else ""
        label += f"&E_c = {Ec} MPa\\\\" if Ec_visible(params) else ""
        label += f"&E_s = {Es} MPa\\\\" if Es_visible(params) else ""
        label += f"&N_c = {Nc} kN\\\\" if Nc_visible(params) else ""
        label += f"&M = {M} kN·m\\\\" if M_visible(params) else ""
        label += f"\\end{{aligned}}$$"

        if ptype == popts[0]:
            model = Calc0(b, h, h0, bf_, hf_, As, fc, ft, fy, Ec, Es)
        elif ptype == popts[1]:
            model = Calc1(b, h, h0, bf_, hf_, fc, ft, fy, Es, M)
        elif ptype == popts[2]:
            model = Calc2(l0, b, h, a_s, a_s_, As, As_, fc, fy, fy_, Es, e_0)
        elif ptype == popts[3]:
            model = Calc3(l0, b, h, a_s, a_s_, As, As_, fc, fy, fy_, Es, Nc)
        elif ptype == popts[4]:
            model = Calc4(l0, b, h, a_s, a_s_, fc, fy, fy_, Es, Nc, M)
        elif ptype == popts[5]:
            model = Calc5(l0, b, h, a_s, a_s_, As_, fc, fy, fy_, Es, Nc, M)
        elif ptype == popts[6]:
            model = Calc6(b, h, h0, As, fc, ft, fy, Ec, Es)
        elif ptype == popts[7]:
            model = Calc7(b, h, h0, fc, ft, fy, Es, M)
        elif ptype == popts[8]:
            model = Calc8(b, h, h0, a_s_, As, As_, fc, ft, fy, fy_, Es)
        elif ptype == popts[9]:
            model = Calc9(b, h, h0, a_s_, fc, ft, fy, fy_, Es, M)
        elif ptype == popts[10]:
            model = Calc10(b, h, h0, a_s_, As_, fc, ft, fy, fy_, Es, M)
        elif ptype == popts[11]:
            model = Calc11(l0, b,  h, a_s, fc, fy, Es, Nc, M)
        else:
            return DataResult(DataGroup(DataItem("Error", "未知题目类型")))
        
        try:
            model.solve()
            model.process()
            return DataResult(DataGroup(DataItem(label, model.process())))
        except:
            return DataResult(DataGroup(DataItem("Error", "计算出错")))