import os
from pymol import cmd
RECEPTOR = r'/home/xiap/bio_project/for_others/AutoMD/AutoVS/temp/receptor/EGFR_demo/receptor_clean.pdb'
OUT = r'/home/xiap/bio_project/for_others/AutoMD/AutoVS/output_figure/EGFR_demo/smoke'
POSES = [('/home/xiap/bio_project/for_others/AutoMD/AutoVS/temp/viz/EGFR_demo/pose_R01.sdf', 'RISDIPLAM'), ('/home/xiap/bio_project/for_others/AutoMD/AutoVS/temp/viz/EGFR_demo/pose_R02.sdf', 'BUCLIZINE HYDROCHLORIDE'), ('/home/xiap/bio_project/for_others/AutoMD/AutoVS/temp/viz/EGFR_demo/pose_R03.sdf', 'VILAZODONE HYDROCHLORIDE')]
RANKS = [1, 2, 3]
COLORS = ['yellow', 'magenta', 'cyan']
import numpy as np
def base_scene():
    cmd.reinitialize()
    cmd.bg_color('white'); cmd.set('ray_opaque_background', 1)
    cmd.set('ray_shadows', 0); cmd.set('depth_cue', 0)
    cmd.set('orthoscopic', 1)  # 正交投影：标签前置时屏幕位置不漂移
    cmd.set('cartoon_transparency', 0.45); cmd.set('stick_radius', 0.16)
    cmd.set('dash_width', 6); cmd.set('dash_gap', 0.3); cmd.set('dash_length', 0.28)
    cmd.set('label_font_id', 5); cmd.set('label_size', 17); cmd.set('label_color', 'black')
    cmd.set('two_sided_lighting', 1); cmd.set('ambient', 0.35); cmd.set('specular', 0.1)
    cmd.load(RECEPTOR, 'rec')
    cmd.hide('everything'); cmd.show('cartoon', 'rec'); cmd.color('skyblue', 'rec')

def place_labels(sel, ref='lig', tag='R', nb=0):
    """标签直接邻近标注（期刊口袋特写样式，像素级探针自标定）：
    每个残基标签沿“配体中心 -> CA”径向放在 CA 外侧近旁（口袋口方向一般是空地），
    文字压配体投影或与已放文字重叠时沿径向继续外推；仅对被推得偏远的标签
    （文字近缘距 CA 0.55~1.1Å）补短引导虚线（文字近端边缘 -> CA），条数受
    全画面虚线总预算约束：3 - nb（nb=本姿势氢键虚线条数，氢键优先）。
    标签与分子分层渲染后乘法合成（见 _composite_label_layers），文字永远纯黑顶层。
    自标定的必要性（实测）：set field 缩放时视口中心随 zoom 场逐姿势漂移 1-3Å，
    get_view 旋转矩阵行主序 reshape 后须转置——理论换算无法消除全部漂移，
    故每姿势渲染 3 个纯色探针（隐藏分子、平面光照），直接实测该视角的
    仿射映射（中心/px/Å/滚转）与文字锚定偏移，标签按像素目标反推世界坐标。"""
    labs = []
    cmd.iterate_state(-1, sel + ' and name CA', 'L.append([x, y, z, resn, str(resi), chain])', space={'L': labs})
    if not labs:
        return
    view = cmd.get_view()
    R = np.array(view[:9]).reshape(3, 3).T  # 实测：不转置会带 9.6° 假滚转 + x/y 各向异性
    pre = np.array(view[9:12]); post = np.array(view[12:15])
    def to_cam(p): return R @ (np.asarray(p, float) + pre) + post
    def to_world(c): return R.T @ (np.asarray(c, float) - post) - pre
    ext = cmd.get_extent(ref)
    ctr = to_cam((np.asarray(ext[0]) + np.asarray(ext[1])) / 2)  # zoom 目标包围盒中心（相机系）
    z_front = float(ctr[2])  # 包围盒中心深度：必在 near/far clip 内（标签层单独渲染，无遮挡问题）
    # ---- 探针自标定：3 纯色球（平面光照去着色偏移）+ 1 条 6 字符测试文字，同分辨率 ray ----
    W, H = 1600, 1200
    probes = [('prA', 0.0, 0.0, (1.0, 0.0, 0.0)), ('prB', 6.0, 0.0, (0.0, 1.0, 0.0)), ('prC', 0.0, 4.5, (0.0, 0.0, 1.0))]
    def _calibrate():
        old = [cmd.get('ambient'), cmd.get('direct'), cmd.get('specular')]
        cmd.set('ambient', 1); cmd.set('direct', 0); cmd.set('specular', 0)
        _hbs = [_n for _n in cmd.get_names('objects') if _n.startswith('hb')]
        for _o in ['rec', 'lig'] + _hbs:
            try: cmd.disable(_o)
            except Exception: pass
        for nm, dx, dy, rgb in probes:
            w = to_world([ctr[0] + dx, ctr[1] + dy, z_front])
            cmd.pseudoatom(nm, pos=[float(v) for v in w])
            cmd.set('sphere_scale', 0.35, nm); cmd.show('spheres', nm)
            cmd.set_color('c_' + nm, list(rgb)); cmd.color('c_' + nm, nm)
        cmd.label('prA', "'XXXXXX'")
        cmd.hide('nonbonded', 'prA or prB or prC')
        cal_png = os.path.join(OUT, '_calib_%s.png' % tag)
        cmd.ray(W, H); cmd.png(cal_png)
        from PIL import Image
        _im = np.asarray(Image.open(cal_png).convert('RGB')).astype(int)
        _cen = {}
        for nm, dx, dy, rgb in probes:
            d = np.abs(_im - (np.array(rgb) * 255)).sum(axis=2)
            ys, xs = np.where(d < 90)
            _cen[nm] = [float(xs.mean()), float(ys.mean())]
        cx0, cy0 = int(_cen['prA'][0]), int(_cen['prA'][1])
        _win = np.zeros(_im.shape[:2], bool)
        _win[max(0, cy0 - 40):cy0 + 40, max(0, cx0 - 150):cx0 + 150] = True
        _ys, _xs = np.where((_im.max(axis=2) < 60) & _win)  # prA 文字（纯黑）像素
        dx_text = float(_xs.min()) - _cen['prA'][0]                       # 文字左缘 - 锚点
        dy_text = float((_ys.min() + _ys.max()) / 2.0) - _cen['prA'][1]   # 文字竖直中心 - 锚点
        tw = float(_xs.max() - _xs.min())                                 # 6 字符文字宽度 px
        cmd.delete('prA'); cmd.delete('prB'); cmd.delete('prC')
        os.remove(cal_png)
        for _o in ['rec', 'lig'] + _hbs:
            try: cmd.enable(_o)
            except Exception: pass
        cmd.set('ambient', old[0]); cmd.set('direct', old[1]); cmd.set('specular', old[2])
        M = np.array([[dx, dy, 1.0] for _, dx, dy, _ in probes])
        D = np.array([_cen[nm] for nm, _, _, _ in probes])
        # 仿射：[px;py] = sol[[0],[1]]@[相机偏移 Å] + sol[2]（中心）——滚转/各向异性/中心漂移全部实测吸收
        return np.linalg.solve(M, D), dx_text, dy_text, tw, float(_ys.max() - _ys.min())
    sol, dx_text, dy_text, tw, th = _calibrate()
    pxA = abs(float(sol[0, 0]))
    # field 的实际视野随各姿势 zoom 场漂移（实测 21-28Å@field=14）：按实测宽/目标宽闭环校正一次，
    # 保证各姿势取景一致（可视宽 28Å）；二次标定的仿射即为最终精确映射。
    target_pxA = W / 28.0
    if abs(pxA - target_pxA) / target_pxA > 0.02:
        cmd.set('field', float(cmd.get('field')) * pxA / target_pxA)
        sol, dx_text, dy_text, tw, th = _calibrate()
        pxA = abs(float(sol[0, 0]))
    def cam2px(cc):
        return np.array([cc[0] * sol[0, 0] + cc[1] * sol[1, 0] + sol[2, 0], cc[0] * sol[0, 1] + cc[1] * sol[1, 1] + sol[2, 1]])
    def px2cam(pp):
        return np.linalg.inv(sol[:2]) @ (np.asarray(pp, float) - sol[2])
    spc = [to_cam([x, y, z]) for x, y, z, _, _, _ in labs]  # 残基 CA 相机坐标
    capx = [cam2px([c[0] - ctr[0], c[1] - ctr[1]]) for c in spc]  # CA 屏幕像素
    # 配体投影包围盒（extent 8 角点），标签文字不得侵入，防止压配体
    _lp = []
    for _cx in ext[0][0], ext[1][0]:
        for _cy in ext[0][1], ext[1][1]:
            for _cz in ext[0][2], ext[1][2]:
                _c = to_cam([_cx, _cy, _cz])
                _lp.append(cam2px([_c[0] - ctr[0], _c[1] - ctr[1]]))
    lig_l, lig_r = min(p[0] for p in _lp), max(p[0] for p in _lp)
    lig_t, lig_b = min(p[1] for p in _lp), max(p[1] for p in _lp)
    # 标签沿“配体中心 -> CA”径向放在 CA 外侧近旁（无引导线）：口袋残基朝口袋口
    # 方向一般是空地，文字离 CA 最近；文字压配体投影或与已放文字重叠时沿径向外推
    placed = []
    _info = {}
    order2 = sorted(range(len(labs)), key=lambda k: float(np.hypot(capx[k][0] - sol[2, 0], capx[k][1] - sol[2, 1])))
    for k in order2:
        ux, uy = capx[k][0] - sol[2, 0], capx[k][1] - sol[2, 1]
        n = float(np.hypot(ux, uy))
        ux, uy = (ux / n, uy / n) if n >= 1.0 else (-1.0, 0.0)
        cx = capx[k][0] + ux * (tw / 2.0 + 0.3 * pxA)  # 文字中心：文字近缘离 CA 约 pad
        cy = capx[k][1] + uy * (tw / 2.0 + 0.3 * pxA)
        for _ in range(4):
            _hit_lig = lig_l - tw / 2.0 < cx < lig_r + tw / 2.0 and lig_t - th / 2.0 < cy < lig_b + th / 2.0
            _hit_lab = any(abs(cx - pc) < tw and abs(cy - py) < th for pc, py in placed)
            if _hit_lig or _hit_lab:
                cx += ux * 0.5 * pxA
                cy += uy * 0.5 * pxA
            else:
                break
        cx = min(max(cx, 8 + tw / 2.0), W - 8 - tw / 2.0)  # 文字不出画布左右边
        cy = min(max(cy, 8 + th / 2.0), H - 8 - th / 2.0)  # 文字不出画布上下边
        placed.append((cx, cy))
        _info[k] = (cx, cy, ux, uy)
        # 标签 pseudoatom（像素目标反推世界坐标）+ 文字；dy_text 校正文字竖直居中
        off = px2cam([cx, cy - dy_text])
        w = to_world([ctr[0] + off[0], ctr[1] + off[1], z_front])
        nm = 'lab_%d' % k
        cmd.pseudoatom(nm, pos=[float(v) for v in w])
        cmd.label(nm, repr(labs[k][3] + labs[k][4]))
        cmd.hide('nonbonded', nm)
        cmd.set('label_size', 15, nm); cmd.set('label_color', 'black', nm)
        print('LABDBG %s k=%d %s cx=%.0f cy=%.0f px/A=%.2f' % (tag, k, labs[k][3] + labs[k][4], cx, cy, pxA))
    # 短引导虚线（受总预算约束）：仅给被避让外推推得偏远（文字近缘距 CA 0.55~1.1Å）
    # 的标签，从文字近端边缘画到 CA；氢键虚线优先占预算，剩余名额取最需要的前几名
    _budget = max(0, 3 - nb)
    _far = []
    for k in order2:
        cx, cy, ux, uy = _info[k]
        _thw = min(tw / 2.0 / abs(ux) if ux else 1e9, th / 2.0 / abs(uy) if uy else 1e9)  # 射线到文字矩形边界
        _dn = float(np.hypot(cx - capx[k][0], cy - capx[k][1])) - _thw
        if 0.55 * pxA < _dn < 1.1 * pxA:
            _far.append((_dn, k, cx, cy, ux, uy, _thw))
    _far.sort(key=lambda t: -t[0])
    for _dn, k, cx, cy, ux, uy, _thw in _far[:_budget]:
        off2 = px2cam([cx - ux * _thw, cy - uy * _thw])  # 文字近端边缘点
        w2 = to_world([ctr[0] + off2[0], ctr[1] + off2[1], z_front])
        cmd.pseudoatom('lda_%d' % k, pos=[float(v) for v in w2])
        cmd.hide('nonbonded', 'lda_%d' % k)
        lead = 'lead_%d' % k
        cmd.distance(lead, 'lda_%d' % k, 'chain %s and resi %s and name CA' % (labs[k][5], labs[k][4]))
        cmd.hide('labels', lead)
        cmd.set('dash_width', 1.5, lead); cmd.set('dash_color', 'gray50', lead); cmd.set('dash_radius', 0.045, lead)
        print('LABDBG %s lead k=%d len=%.0fpx' % (tag, k, _dn))

# ---- Top N 叠加总览 ----
base_scene()
for i, (sdf, name) in enumerate(POSES):
    cmd.load(sdf, 'p%d' % i); cmd.show('sticks', 'p%d' % i)
    cmd.color(COLORS[i], 'p%d and elem C' % i)
cmd.orient('p0', animate=0); cmd.zoom('p0', 10, animate=0); cmd.ray(1800, 1350)
cmd.png(os.path.join(OUT, 'top_overview.png'), dpi=300)

# ---- 逐姿势特写：紫碳配体 + 黄碳关键残基侧链 + 绿虚线氢键 + 标签 ----
for i, (sdf, name) in enumerate(POSES):
    base_scene()
    cmd.load(sdf, 'lig'); cmd.show('sticks', 'lig')
    cmd.color('violetpurple', 'lig and elem C'); cmd.set('stick_radius', 0.24, 'lig')
    cmd.select('pocket', 'byres rec within 4.5 of lig')
    cmd.set('cartoon_transparency', 0.82, 'rec')  # 全受体极淡 cartoon：无切割边缘，消除口袋聚焦的碎片感
    cmd.set('cartoon_transparency', 0.45, 'pocket')  # 口袋区加实，突出结合位点
    # 氢键逐对画并计数：取距离最近的 3 对（全画面虚线总预算 3 = 氢键 + 引导线）
    _pl = []; _pr = []
    cmd.iterate_state(-1, 'lig and (elem N or elem O or elem S)', 'L.append([x, y, z, index])', space={'L': _pl})
    cmd.iterate_state(-1, 'pocket and (elem N or elem O or elem S)', 'L.append([x, y, z, index])', space={'L': _pr})
    _pairs = []
    for _a in _pl:
        for _b in _pr:
            _d = float((_a[0] - _b[0])**2 + (_a[1] - _b[1])**2 + (_a[2] - _b[2])**2)**0.5
            if _d < 3.5:
                _pairs.append((_d, _a[3], _b[3]))
    _pairs.sort()
    n_hb = 0
    for _d, _ia, _ib in _pairs[:3]:
        cmd.distance('hb_%d' % n_hb, 'lig and index %d' % _ia, 'pocket and index %d' % _ib)
        n_hb += 1
    # 逐对象设置样式（该版本 PyMOL 对 measurement 对象名的 or 链选择解析有 bug，单名正常）
    for _j in range(n_hb):
        _hn = 'hb_%d' % _j
        try:
            cmd.hide('labels', _hn)
            cmd.color('forest', _hn); cmd.set('dash_radius', 0.06, _hn); cmd.set('dash_width', 2, _hn)
        except Exception:
            pass
    cmd.select('ires_R1', 'pocket and ((resi 694 and chain A) or (resi 702 and chain A) or (resi 719 and chain A) or (resi 721 and chain A) or (resi 764 and chain A) or (resi 766 and chain A) or (resi 769 and chain A))')
    cmd.show('sticks', 'ires_R1 and sidechain')
    cmd.color('paleyellow', 'ires_R1 and elem C')
    cmd.select('ires_R2', 'pocket and ((resi 694 and chain A) or (resi 699 and chain A) or (resi 702 and chain A) or (resi 719 and chain A) or (resi 768 and chain A) or (resi 769 and chain A) or (resi 831 and chain A))')
    cmd.show('sticks', 'ires_R2 and sidechain')
    cmd.color('paleyellow', 'ires_R2 and elem C')
    cmd.select('ires_R3', 'pocket and ((resi 699 and chain A) or (resi 702 and chain A) or (resi 721 and chain A) or (resi 766 and chain A) or (resi 820 and chain A) or (resi 831 and chain A))')
    cmd.show('sticks', 'ires_R3 and sidechain')
    cmd.color('paleyellow', 'ires_R3 and elem C')
    cmd.orient('lig', animate=0); cmd.zoom('lig', 5.5, animate=0)
    cmd.set('field', 14)  # ortho 视野：field=可视半宽 -> 可视 28x21Å，配体占幅 ~45%，两侧留标签栏（place_labels 与此耦合）
    place_labels('ires_R%d' % RANKS[i] if 'ires_R%d' % RANKS[i] in cmd.get_names('selections') else 'pocket', tag='R%02d' % RANKS[i], nb=n_hb)
    cmd.ray(1600, 1200)
    cmd.png(os.path.join(OUT, 'pose_R%02d.png' % RANKS[i]), dpi=300)
    # 标签层：隐藏分子对象（rec/lig/hb*），仅留标签文字 + 引导线，与分子层同视角二次 ray
    for _o in ['rec', 'lig'] + [_n for _n in cmd.get_names('objects') if _n.startswith('hb')]:
        try:
            cmd.disable(_o)
        except Exception:
            pass
    cmd.ray(1600, 1200)
    cmd.png(os.path.join(OUT, '_labels_R%02d.png' % RANKS[i]), dpi=300)
    for _o in ['rec', 'lig'] + [_n for _n in cmd.get_names('objects') if _n.startswith('hb')]:
        try:
            cmd.enable(_o)
        except Exception:
            pass
cmd.save(os.path.join(OUT, 'poses_3d_session.pse'), 'all')
