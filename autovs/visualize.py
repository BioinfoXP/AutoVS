"""可视化模块：Top 命中姿势的 3D 渲染（PyMOL）与 2D 相互作用图（几何判据）。

自动化定位：`autovs dock` 汇总后自动调用 `run_visualize()`；也可独立执行
`autovs visualize --target ... --tag ...`。

3D：受体 cartoon 半透明 + 口袋残基 stick + 配体 stick（黄色碳）+ 极性接触虚线，
    逐姿势 PNG（白底透明、300 dpi）+ Top N 叠加总览图。
2D：基于距离/角度几何判据的配体-残基相互作用检测
    （氢键 H···A ≤2.6 Å 且 D-H···A >110°；疏水 C···C ≤4.0 Å；
     盐桥 |q|≥1 原子 ↔ 带电侧链 ≤4.0 Å；π 堆积环心距 ≤5.5 Å），
    绘制以配体 2D 结构为中心的径向示意图（PDF 矢量 + PNG）。
    注意：2D 图用于快速检视结合模式，定量结论需另行复核。
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

import numpy as np

from .config import Config
from .utils import ensure_dir, get_logger

log = get_logger("visualize")

# 相互作用类型 -> 颜色
INTERACTION_COLORS = [
    ("HBDonor", "#1e8e3e"),      # 配体为氢键给体（DS 惯例：氢键绿色虚线）
    ("HBAcceptor", "#2ca02c"),   # 配体为氢键受体
    ("SaltBridge", "#d62728"),
    ("PiStacking", "#c24fb0"),
    ("PiCation", "#e377c2"),
    ("Hydrophobic", "#b0b6ba"),
]
# 残基节点配色（DS 风格：类型 -> (面圆形, 文字色)）
NODE_COLORS = {
    "HBDonor": ("#1e8e3e", "white"), "HBAcceptor": ("#1e8e3e", "white"),
    "SaltBridge": ("#d64545", "white"),
    "PiStacking": ("#f2c4ec", "#5a1e50"), "PiCation": ("#f6d3ee", "#5a1e50"),
    "Hydrophobic": ("#cde9cd", "#1e4d1e"),
}
# 姿势颜色（全部经 cmd.get_color_index 验证有效；violetpink/dashblue 在本 PyMOL 不存在）
PYMOL_POSE_COLORS = ["yellow", "magenta", "cyan", "orange", "salmon",
                     "lime", "warmpink", "skyblue", "wheat", "deepteal"]
# 带电侧链原子（盐桥判据）
POSITIVE_SIDECHAIN = {"ARG": {"NH1", "NH2"}, "LYS": {"NZ"}, "HIS": {"ND1", "NE2"}}
NEGATIVE_SIDECHAIN = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}}
AROMATIC_RESIDUES = {"PHE", "TYR", "TRP", "HIS"}


def _inter_color(name: str) -> str:
    for prefix, color in INTERACTION_COLORS:
        if name.startswith(prefix):
            return color
    return "#666666"


def load_top_poses(ranking_csv: Path, top: int) -> list[dict[str, Any]]:
    """从排行榜取前 top 个有 pose SDF 的命中，返回 [{rank, entry_id, name, score, pose_sdf}]。"""
    with Path(ranking_csv).open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    poses = []
    for row in rows:
        if len(poses) >= int(top):
            break
        sdf = row.get("pose_sdf", "")
        if sdf and Path(sdf).exists():
            poses.append({"rank": int(row["rank"]), "entry_id": row["entry_id"],
                          "name": row.get("name", ""),
                          "score": row.get("top_score_kcal_mol", ""), "pose_sdf": sdf})
    return poses


def _pose_label(p: dict) -> str:
    """统一产物前缀：pose_R01_MIDOSTAURIN——3D/2D/对照图共用，文件名即对应关系。"""
    nm = (p.get("name") or p["entry_id"]).strip().replace(" ", "_")
    return f"pose_R{int(p['rank']):02d}_{nm[:40]}"


def match_entries(ranking_csv: Path, entries: list[str]) -> list[dict[str, Any]]:
    """按 entry_id 精确或药名子串（大小写不敏感）从排行榜取指定化合物（保持去重与出现顺序）。"""
    with Path(ranking_csv).open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    out: list[dict[str, Any]] = []
    for ent in entries:
        e = ent.strip().lower()
        if not e:
            continue
        hit = next((r for r in rows if r["entry_id"].lower() == e), None) \
            or next((r for r in rows if e in (r.get("name") or "").lower()), None)
        if not hit:
            log.warning("排行榜中未找到: %s", ent)
            continue
        if not (hit.get("pose_sdf") and Path(hit["pose_sdf"]).exists()):
            log.warning("%s（%s）无可用 pose SDF，跳过", hit["entry_id"], hit.get("name", ""))
            continue
        if all(d["entry_id"] != hit["entry_id"] for d in out):
            out.append({"rank": int(hit["rank"]), "entry_id": hit["entry_id"],
                        "name": hit.get("name", ""),
                        "score": hit.get("top_score_kcal_mol", ""), "pose_sdf": hit["pose_sdf"]})
    return out


def _pair_image(png3d: Path, png2d: Path, out_png: Path) -> Path:
    """3D|2D 等高横拼对照图：一张图内直接建立口袋姿势与相互作用平面的对应关系。"""
    from PIL import Image
    a, b = Image.open(png3d).convert("RGB"), Image.open(png2d).convert("RGB")
    h = min(a.height, b.height)
    a = a.resize((max(1, int(a.width * h / a.height)), h))
    b = b.resize((max(1, int(b.width * h / b.height)), h))
    canvas = Image.new("RGB", (a.width + b.width + 14, h), "white")
    canvas.paste(a, (0, 0))
    canvas.paste(b, (a.width + 14, 0))
    canvas.save(out_png)
    return out_png


def extract_first_pose(sdf_path: str | Path, out_sdf: Path) -> Path:
    """Vina/meeko 输出的 pose.sdf 往往含多条 pose 记录，只保留打分最高的第一条。

    PyMOL 按多 state 加载这类文件后，zoom/get_extent 会覆盖全部 pose 的总包围盒，
    导致视野估算与标签布局失真（标签被推到画面外），必须渲染前抽成单分子 SDF。
    """
    from rdkit import Chem
    mol = next(m for m in Chem.SDMolSupplier(str(sdf_path), removeHs=False) if m)
    if mol is None:
        raise ValueError(f"pose 不可读: {sdf_path}")
    w = Chem.SDWriter(str(out_sdf))
    w.write(mol)
    w.close()
    return out_sdf


# =====================================================================
# 3D 渲染（PyMOL 批处理：独立进程，避免与主程序耦合）
# =====================================================================
def _pymol_script(receptor_pdb: Path, poses: list[dict], out_dir: Path,
                  resi_map: dict[int, list[str]] | None = None) -> str:
    """生成 PyMOL 批处理脚本。resi_map: {rank: ['95', '98', ...]} 相互作用残基编号（打标签用）。"""
    resi_map = resi_map or {}
    lines = [
        "import os", "from pymol import cmd",
        f"RECEPTOR = r'{receptor_pdb}'", f"OUT = r'{out_dir}'",
        f"POSES = {[(p['pose_sdf'], p['name'] or p['entry_id']) for p in poses]!r}",
        f"RANKS = {[p['rank'] for p in poses]!r}",
        f"COLORS = {PYMOL_POSE_COLORS[:len(poses)]!r}",
        "import numpy as np",
        "def base_scene():",
        "    cmd.reinitialize()",
        "    cmd.bg_color('white'); cmd.set('ray_opaque_background', 1)",
        "    cmd.set('ray_shadows', 0); cmd.set('depth_cue', 0)",
        "    cmd.set('orthoscopic', 1)  # 正交投影：标签前置时屏幕位置不漂移",
        "    cmd.set('cartoon_transparency', 0.45); cmd.set('stick_radius', 0.16)",
        "    cmd.set('dash_width', 6); cmd.set('dash_gap', 0.3); cmd.set('dash_length', 0.28)",
        "    cmd.set('label_font_id', 5); cmd.set('label_size', 17); cmd.set('label_color', 'black')",
        "    cmd.set('two_sided_lighting', 1); cmd.set('ambient', 0.35); cmd.set('specular', 0.1)",
        "    cmd.load(RECEPTOR, 'rec')",
        "    cmd.hide('everything'); cmd.show('cartoon', 'rec'); cmd.color('skyblue', 'rec')",
        "",
        "def place_labels(sel, ref='lig', tag='R', nb=0):",
        "    \"\"\"标签直接邻近标注（期刊口袋特写样式，像素级探针自标定）：",
        "    每个残基标签沿“配体中心 -> CA”径向放在 CA 外侧近旁（口袋口方向一般是空地），",
        "    文字压配体投影或与已放文字重叠时沿径向继续外推；仅对被推得偏远的标签",
        "    （文字近缘距 CA 0.55~1.1Å）补短引导虚线（文字近端边缘 -> CA），条数受",
        "    全画面虚线总预算约束：3 - nb（nb=本姿势氢键虚线条数，氢键优先）。",
        "    标签与分子分层渲染后乘法合成（见 _composite_label_layers），文字永远纯黑顶层。",
        "    自标定的必要性（实测）：set field 缩放时视口中心随 zoom 场逐姿势漂移 1-3Å，",
        "    get_view 旋转矩阵行主序 reshape 后须转置——理论换算无法消除全部漂移，",
        "    故每姿势渲染 3 个纯色探针（隐藏分子、平面光照），直接实测该视角的",
        "    仿射映射（中心/px/Å/滚转）与文字锚定偏移，标签按像素目标反推世界坐标。\"\"\"",
        "    labs = []",
        "    cmd.iterate_state(-1, sel + ' and name CA', 'L.append([x, y, z, resn, str(resi), chain])', space={'L': labs})",
        "    if not labs:",
        "        return",
        "    view = cmd.get_view()",
        "    R = np.array(view[:9]).reshape(3, 3).T  # 实测：不转置会带 9.6° 假滚转 + x/y 各向异性",
        "    pre = np.array(view[9:12]); post = np.array(view[12:15])",
        "    def to_cam(p): return R @ (np.asarray(p, float) + pre) + post",
        "    def to_world(c): return R.T @ (np.asarray(c, float) - post) - pre",
        "    ext = cmd.get_extent(ref)",
        "    ctr = to_cam((np.asarray(ext[0]) + np.asarray(ext[1])) / 2)  # zoom 目标包围盒中心（相机系）",
        "    z_front = float(ctr[2])  # 包围盒中心深度：必在 near/far clip 内（标签层单独渲染，无遮挡问题）",
        "    # ---- 探针自标定：3 纯色球（平面光照去着色偏移）+ 1 条 6 字符测试文字，同分辨率 ray ----",
        "    W, H = 1600, 1200",
        "    probes = [('prA', 0.0, 0.0, (1.0, 0.0, 0.0)), ('prB', 6.0, 0.0, (0.0, 1.0, 0.0)), ('prC', 0.0, 4.5, (0.0, 0.0, 1.0))]",
        "    def _calibrate():",
        "        old = [cmd.get('ambient'), cmd.get('direct'), cmd.get('specular')]",
        "        cmd.set('ambient', 1); cmd.set('direct', 0); cmd.set('specular', 0)",
        "        _hbs = [_n for _n in cmd.get_names('objects') if _n.startswith('hb')]",
        "        for _o in ['rec', 'lig'] + _hbs:",
        "            try: cmd.disable(_o)",
        "            except Exception: pass",
        "        for nm, dx, dy, rgb in probes:",
        "            w = to_world([ctr[0] + dx, ctr[1] + dy, z_front])",
        "            cmd.pseudoatom(nm, pos=[float(v) for v in w])",
        "            cmd.set('sphere_scale', 0.35, nm); cmd.show('spheres', nm)",
        "            cmd.set_color('c_' + nm, list(rgb)); cmd.color('c_' + nm, nm)",
        "        cmd.label('prA', \"'XXXXXX'\")",
        "        cmd.hide('nonbonded', 'prA or prB or prC')",
        "        cal_png = os.path.join(OUT, '_calib_%s.png' % tag)",
        "        cmd.ray(W, H); cmd.png(cal_png)",
        "        from PIL import Image",
        "        _im = np.asarray(Image.open(cal_png).convert('RGB')).astype(int)",
        "        _cen = {}",
        "        for nm, dx, dy, rgb in probes:",
        "            d = np.abs(_im - (np.array(rgb) * 255)).sum(axis=2)",
        "            ys, xs = np.where(d < 90)",
        "            _cen[nm] = [float(xs.mean()), float(ys.mean())]",
        "        cx0, cy0 = int(_cen['prA'][0]), int(_cen['prA'][1])",
        "        _win = np.zeros(_im.shape[:2], bool)",
        "        _win[max(0, cy0 - 40):cy0 + 40, max(0, cx0 - 150):cx0 + 150] = True",
        "        _ys, _xs = np.where((_im.max(axis=2) < 60) & _win)  # prA 文字（纯黑）像素",
        "        dx_text = float(_xs.min()) - _cen['prA'][0]                       # 文字左缘 - 锚点",
        "        dy_text = float((_ys.min() + _ys.max()) / 2.0) - _cen['prA'][1]   # 文字竖直中心 - 锚点",
        "        tw = float(_xs.max() - _xs.min())                                 # 6 字符文字宽度 px",
        "        cmd.delete('prA'); cmd.delete('prB'); cmd.delete('prC')",
        "        os.remove(cal_png)",
        "        for _o in ['rec', 'lig'] + _hbs:",
        "            try: cmd.enable(_o)",
        "            except Exception: pass",
        "        cmd.set('ambient', old[0]); cmd.set('direct', old[1]); cmd.set('specular', old[2])",
        "        M = np.array([[dx, dy, 1.0] for _, dx, dy, _ in probes])",
        "        D = np.array([_cen[nm] for nm, _, _, _ in probes])",
        "        # 仿射：[px;py] = sol[[0],[1]]@[相机偏移 Å] + sol[2]（中心）——滚转/各向异性/中心漂移全部实测吸收",
        "        return np.linalg.solve(M, D), dx_text, dy_text, tw, float(_ys.max() - _ys.min())",
        "    sol, dx_text, dy_text, tw, th = _calibrate()",
        "    pxA = abs(float(sol[0, 0]))",
        "    # field 的实际视野随各姿势 zoom 场漂移（实测 21-28Å@field=14）：按实测宽/目标宽闭环校正一次，",
        "    # 保证各姿势取景一致（可视宽 28Å）；二次标定的仿射即为最终精确映射。",
        "    target_pxA = W / 28.0",
        "    if abs(pxA - target_pxA) / target_pxA > 0.02:",
        "        cmd.set('field', float(cmd.get('field')) * pxA / target_pxA)",
        "        sol, dx_text, dy_text, tw, th = _calibrate()",
        "        pxA = abs(float(sol[0, 0]))",
        "    def cam2px(cc):",
        "        return np.array([cc[0] * sol[0, 0] + cc[1] * sol[1, 0] + sol[2, 0], cc[0] * sol[0, 1] + cc[1] * sol[1, 1] + sol[2, 1]])",
        "    def px2cam(pp):",
        "        return np.linalg.inv(sol[:2]) @ (np.asarray(pp, float) - sol[2])",
        "    spc = [to_cam([x, y, z]) for x, y, z, _, _, _ in labs]  # 残基 CA 相机坐标",
        "    capx = [cam2px([c[0] - ctr[0], c[1] - ctr[1]]) for c in spc]  # CA 屏幕像素",
        "    # 配体投影包围盒（extent 8 角点），标签文字不得侵入，防止压配体",
        "    _lp = []",
        "    for _cx in ext[0][0], ext[1][0]:",
        "        for _cy in ext[0][1], ext[1][1]:",
        "            for _cz in ext[0][2], ext[1][2]:",
        "                _c = to_cam([_cx, _cy, _cz])",
        "                _lp.append(cam2px([_c[0] - ctr[0], _c[1] - ctr[1]]))",
        "    lig_l, lig_r = min(p[0] for p in _lp), max(p[0] for p in _lp)",
        "    lig_t, lig_b = min(p[1] for p in _lp), max(p[1] for p in _lp)",
        "    # 标签沿“配体中心 -> CA”径向放在 CA 外侧近旁（无引导线）：口袋残基朝口袋口",
        "    # 方向一般是空地，文字离 CA 最近；文字压配体投影或与已放文字重叠时沿径向外推",
        "    placed = []",
        "    _info = {}",
        "    order2 = sorted(range(len(labs)), key=lambda k: float(np.hypot(capx[k][0] - sol[2, 0], capx[k][1] - sol[2, 1])))",
        "    for k in order2:",
        "        ux, uy = capx[k][0] - sol[2, 0], capx[k][1] - sol[2, 1]",
        "        n = float(np.hypot(ux, uy))",
        "        ux, uy = (ux / n, uy / n) if n >= 1.0 else (-1.0, 0.0)",
        "        cx = capx[k][0] + ux * (tw / 2.0 + 0.3 * pxA)  # 文字中心：文字近缘离 CA 约 pad",
        "        cy = capx[k][1] + uy * (tw / 2.0 + 0.3 * pxA)",
        "        for _ in range(4):",
        "            _hit_lig = lig_l - tw / 2.0 < cx < lig_r + tw / 2.0 and lig_t - th / 2.0 < cy < lig_b + th / 2.0",
        "            _hit_lab = any(abs(cx - pc) < tw and abs(cy - py) < th for pc, py in placed)",
        "            if _hit_lig or _hit_lab:",
        "                cx += ux * 0.5 * pxA",
        "                cy += uy * 0.5 * pxA",
        "            else:",
        "                break",
        "        cx = min(max(cx, 8 + tw / 2.0), W - 8 - tw / 2.0)  # 文字不出画布左右边",
        "        cy = min(max(cy, 8 + th / 2.0), H - 8 - th / 2.0)  # 文字不出画布上下边",
        "        placed.append((cx, cy))",
        "        _info[k] = (cx, cy, ux, uy)",
        "        # 标签 pseudoatom（像素目标反推世界坐标）+ 文字；dy_text 校正文字竖直居中",
        "        off = px2cam([cx, cy - dy_text])",
        "        w = to_world([ctr[0] + off[0], ctr[1] + off[1], z_front])",
        "        nm = 'lab_%d' % k",
        "        cmd.pseudoatom(nm, pos=[float(v) for v in w])",
        "        cmd.label(nm, repr(labs[k][3] + labs[k][4]))",
        "        cmd.hide('nonbonded', nm)",
        "        cmd.set('label_size', 15, nm); cmd.set('label_color', 'black', nm)",
        "        print('LABDBG %s k=%d %s cx=%.0f cy=%.0f px/A=%.2f' % (tag, k, labs[k][3] + labs[k][4], cx, cy, pxA))",
        "    # 短引导虚线（受总预算约束）：仅给被避让外推推得偏远（文字近缘距 CA 0.55~1.1Å）",
        "    # 的标签，从文字近端边缘画到 CA；氢键虚线优先占预算，剩余名额取最需要的前几名",
        "    _budget = max(0, 3 - nb)",
        "    _far = []",
        "    for k in order2:",
        "        cx, cy, ux, uy = _info[k]",
        "        _thw = min(tw / 2.0 / abs(ux) if ux else 1e9, th / 2.0 / abs(uy) if uy else 1e9)  # 射线到文字矩形边界",
        "        _dn = float(np.hypot(cx - capx[k][0], cy - capx[k][1])) - _thw",
        "        if 0.55 * pxA < _dn < 1.1 * pxA:",
        "            _far.append((_dn, k, cx, cy, ux, uy, _thw))",
        "    _far.sort(key=lambda t: -t[0])",
        "    for _dn, k, cx, cy, ux, uy, _thw in _far[:_budget]:",
        "        off2 = px2cam([cx - ux * _thw, cy - uy * _thw])  # 文字近端边缘点",
        "        w2 = to_world([ctr[0] + off2[0], ctr[1] + off2[1], z_front])",
        "        cmd.pseudoatom('lda_%d' % k, pos=[float(v) for v in w2])",
        "        cmd.hide('nonbonded', 'lda_%d' % k)",
        "        lead = 'lead_%d' % k",
        "        cmd.distance(lead, 'lda_%d' % k, 'chain %s and resi %s and name CA' % (labs[k][5], labs[k][4]))",
        "        cmd.hide('labels', lead)",
        "        cmd.set('dash_width', 1.5, lead); cmd.set('dash_color', 'gray50', lead); cmd.set('dash_radius', 0.045, lead)",
        "        print('LABDBG %s lead k=%d len=%.0fpx' % (tag, k, _dn))",
        "",
        "# ---- Top N 叠加总览 ----",
        "base_scene()",
        "for i, (sdf, name) in enumerate(POSES):",
        "    cmd.load(sdf, 'p%d' % i); cmd.show('sticks', 'p%d' % i)",
        "    cmd.color(COLORS[i], 'p%d and elem C' % i)",
        "cmd.orient('p0', animate=0); cmd.zoom('p0', 10, animate=0); cmd.ray(1800, 1350)",
        "cmd.png(os.path.join(OUT, 'top_overview.png'), dpi=300)",
        "",
        "# ---- 逐姿势特写：紫碳配体 + 黄碳关键残基侧链 + 绿虚线氢键 + 标签 ----",
        "for i, (sdf, name) in enumerate(POSES):",
        "    base_scene()",
        "    cmd.load(sdf, 'lig'); cmd.show('sticks', 'lig')",
        "    cmd.color('violetpurple', 'lig and elem C'); cmd.set('stick_radius', 0.24, 'lig')",
        "    cmd.select('pocket', 'byres rec within 4.5 of lig')",
        "    cmd.set('cartoon_transparency', 0.82, 'rec')  # 全受体极淡 cartoon：无切割边缘，消除口袋聚焦的碎片感",
        "    cmd.set('cartoon_transparency', 0.45, 'pocket')  # 口袋区加实，突出结合位点",
        "    # 氢键逐对画并计数：取距离最近的 3 对（全画面虚线总预算 3 = 氢键 + 引导线）",
        "    _pl = []; _pr = []",
        "    cmd.iterate_state(-1, 'lig and (elem N or elem O or elem S)', 'L.append([x, y, z, index])', space={'L': _pl})",
        "    cmd.iterate_state(-1, 'pocket and (elem N or elem O or elem S)', 'L.append([x, y, z, index])', space={'L': _pr})",
        "    _pairs = []",
        "    for _a in _pl:",
        "        for _b in _pr:",
        "            _d = float((_a[0] - _b[0])**2 + (_a[1] - _b[1])**2 + (_a[2] - _b[2])**2)**0.5",
        "            if _d < 3.5:",
        "                _pairs.append((_d, _a[3], _b[3]))",
        "    _pairs.sort()",
        "    n_hb = 0",
        "    for _d, _ia, _ib in _pairs[:3]:",
        "        cmd.distance('hb_%d' % n_hb, 'lig and index %d' % _ia, 'pocket and index %d' % _ib)",
        "        n_hb += 1",
        "    # 逐对象设置样式（该版本 PyMOL 对 measurement 对象名的 or 链选择解析有 bug，单名正常）",
        "    for _j in range(n_hb):",
        "        _hn = 'hb_%d' % _j",
        "        try:",
        "            cmd.hide('labels', _hn)",
        "            cmd.color('forest', _hn); cmd.set('dash_radius', 0.06, _hn); cmd.set('dash_width', 2, _hn)",
        "        except Exception:",
        "            pass",
    ]
    # 只显示非疏水相互作用残基的侧链（画面聚焦，参考教科书图做法）；无则兜底显示全口袋
    has_ires = False
    for p in poses:
        sel = " or ".join(f"(resi {rid} and chain {chain})"
                          for chain, rid in resi_map.get(p["rank"], []))
        if sel:
            has_ires = True
            lines += [
                f"    cmd.select('ires_R{p['rank']}', 'pocket and ({sel})')",
                f"    cmd.show('sticks', 'ires_R{p['rank']} and sidechain')",
                f"    cmd.color('paleyellow', 'ires_R{p['rank']} and elem C')",
            ]
    if not has_ires:
        lines.append("    cmd.show('sticks', 'pocket and sidechain'); cmd.color('paleyellow', 'pocket and elem C')")
    lines += [
        "    cmd.orient('lig', animate=0); cmd.zoom('lig', 5.5, animate=0)",
        "    cmd.set('field', 14)  # ortho 视野：field=可视半宽 -> 可视 28x21Å，配体占幅 ~45%，两侧留标签栏（place_labels 与此耦合）",
        "    place_labels('ires_R%d' % RANKS[i] if 'ires_R%d' % RANKS[i] in cmd.get_names('selections') else 'pocket', tag='R%02d' % RANKS[i], nb=n_hb)",
        "    cmd.ray(1600, 1200)",
        "    cmd.png(os.path.join(OUT, 'pose_R%02d.png' % RANKS[i]), dpi=300)",
        "    # 标签层：隐藏分子对象（rec/lig/hb*），仅留标签文字 + 引导线，与分子层同视角二次 ray",
        "    for _o in ['rec', 'lig'] + [_n for _n in cmd.get_names('objects') if _n.startswith('hb')]:",
        "        try:",
        "            cmd.disable(_o)",
        "        except Exception:",
        "            pass",
        "    cmd.ray(1600, 1200)",
        "    cmd.png(os.path.join(OUT, '_labels_R%02d.png' % RANKS[i]), dpi=300)",
        "    for _o in ['rec', 'lig'] + [_n for _n in cmd.get_names('objects') if _n.startswith('hb')]:",
        "        try:",
        "            cmd.enable(_o)",
        "        except Exception:",
        "            pass",
        "cmd.save(os.path.join(OUT, 'poses_3d_session.pse'), 'all')",
        "",
    ]
    return "\n".join(lines)


def _composite_label_layers(out_dir: Path) -> None:
    """标签层（白底黑字）× 分子层 乘法合成：文字永远纯黑顶层，不被前景 cartoon 罩灰。

    乘法混合：标签层白底(255) 处分子原样保留，黑字(0)/灰引导线处按比例压暗，
    视觉上等价于把文字印刷在分子图最上层。合成后删除标签层中间文件。
    """
    from PIL import Image

    for lab_png in sorted(out_dir.glob("_labels_R*.png")):
        mol_png = lab_png.with_name(lab_png.name.replace("_labels_", "pose_"))
        if not mol_png.exists():
            continue
        mol = np.asarray(Image.open(mol_png).convert("RGB"), dtype=float)
        lab = np.asarray(Image.open(lab_png).convert("RGB"), dtype=float)
        if mol.shape != lab.shape:
            log.warning("标签层尺寸不一致，跳过合成: %s", lab_png.name)
            continue
        out = np.clip(mol * lab / 255.0, 0, 255).astype(np.uint8)
        Image.fromarray(out).save(mol_png, dpi=(300, 300))
        lab_png.unlink()


def render_poses_3d(cfg: Config, poses: list[dict], out_dir: Path,
                    resi_map: dict[int, list[str]] | None = None) -> list[Path]:
    """PyMOL 渲染 Top 姿势特写 + 总览，返回生成的 PNG 列表。"""
    if not poses:
        return []
    pymol_bin = cfg.tools.pymol
    if not pymol_bin or not Path(pymol_bin).exists():
        log.warning("pymol 不可用，跳过 3D 渲染（framework.yaml tools.pymol）")
        return []
    ensure_dir(out_dir)
    script = out_dir / "_pymol_render.py"
    script.write_text(_pymol_script(cfg.paths.temp_dir / "receptor" / cfg.target.name / "receptor_clean.pdb",
                                    poses, out_dir, resi_map), encoding="utf-8")
    from .utils import run_cmd
    run_cmd([str(pymol_bin), "-cq", str(script)], cwd=out_dir, log=out_dir / "pymol.log", timeout=600)
    _composite_label_layers(out_dir)  # 标签层 x 分子层合成，文字置顶防罩灰
    outputs = sorted(out_dir.glob("pose_R*.png")) + [out_dir / "top_overview.png"]
    log.info("3D 渲染完成: %d 张 PNG -> %s", len(outputs), out_dir)
    return [p for p in outputs if p.exists()]


# =====================================================================
# 2D 相互作用：几何判据检测（自研，无外部转换器依赖）
# =====================================================================
def _parse_pdb_atoms(pdb_path: Path) -> dict[str, Any]:
    """解析受体 PDB 的 ATOM 记录：坐标/元素/残基/氢-重原子归属。"""
    recs = []
    for line in Path(pdb_path).read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        element = line[76:78].strip().upper() or line[12:16].strip()[0]
        recs.append({
            "name": line[12:16].strip(),
            "resname": line[17:20].strip().upper(),
            "chain": line[21:22].strip() or "A",
            "resid": int(line[22:26]),
            "element": element,
            "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
        })
    heavy_idx = [i for i, r in enumerate(recs) if r["element"] != "H"]
    h_idx = [i for i, r in enumerate(recs) if r["element"] == "H"]
    P_heavy = np.array([recs[i]["xyz"] for i in heavy_idx])
    # 氢归属：最近重原子（<1.3 Å）
    h_parent = {}
    for i in h_idx:
        d = np.linalg.norm(P_heavy - recs[i]["xyz"], axis=1)
        j = int(np.argmin(d))
        if d[j] < 1.3:
            h_parent[i] = heavy_idx[j]
    return {"recs": recs, "heavy_idx": heavy_idx, "h_idx": h_idx, "h_parent": h_parent}


def _ring_centroids(mol, conf_xyz: np.ndarray) -> list[tuple[np.ndarray, list[int]]]:
    """配体芳香环：(环心坐标, 环重原子 idx 列表)。"""
    ri = mol.GetRingInfo()
    out = []
    for ring in ri.AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            out.append((conf_xyz[list(ring)].mean(axis=0), [int(i) for i in ring]))
    return out


def _protein_rings(recs: dict[str, Any]) -> list[tuple[str, np.ndarray]]:
    """芳香残基环心：PHE/TYR/TRP/HIS 的六元/五元环（取环原子坐标均值）。"""
    rings = []
    cur_res_atoms: dict[tuple, list[int]] = {}
    for i, r in enumerate(recs["recs"]):
        if r["resname"] in AROMATIC_RESIDUES:
            cur_res_atoms.setdefault((r["resname"], r["chain"], r["resid"]), []).append(i)
    for (resname, chain, resid), idx in cur_res_atoms.items():
        xyz = np.array([recs["recs"][i]["xyz"] for i in idx])
        rings.append((f"{resname}{resid}.{chain}", xyz.mean(axis=0)))
    return rings


def detect_interactions(pose_sdf: Path, receptor_pdb: Path) -> list[dict[str, Any]]:
    """几何判据检测配体-残基相互作用。

    返回去重列表，每项 {residue, type, detail, lig_atoms}——lig_atoms 为参与的
    配体重原子 idx（2D 图把作用线连到配体具体原子用）。
    """
    from rdkit import Chem

    lig = next((m for m in Chem.SDMolSupplier(str(pose_sdf), removeHs=False) if m), None)
    if lig is None:
        raise ValueError(f"pose 不可读: {pose_sdf}")
    L = np.array(lig.GetConformer().GetPositions())

    recs = _parse_pdb_atoms(Path(receptor_pdb))
    R = np.array([recs["recs"][i]["xyz"] for i in recs["heavy_idx"]])
    heavy_meta = [recs["recs"][i] for i in recs["heavy_idx"]]
    out: list[dict[str, Any]] = []

    def res_label(meta) -> str:
        return f"{meta['resname']}{meta['resid']}.{meta['chain']}"

    def add(label, itype, detail, lig_atoms=None):
        if not any(e["residue"] == label and e["type"] == itype for e in out):
            out.append({"residue": label, "type": itype, "detail": detail,
                        "lig_atoms": sorted({int(a) for a in (lig_atoms or [])})})

    # --- 配体原子角色 ---
    l_heavy, l_h = [], []
    for a in lig.GetAtoms():
        if a.GetAtomicNum() == 1:
            l_h.append(a.GetIdx())
        else:
            l_heavy.append(a.GetIdx())
    L_heavy = L[l_heavy]
    # 氢键给体（配体侧）：H 连接 N/O
    lig_donor_h = []
    for hi in l_h:
        parent = next((n for n in lig.GetAtomWithIdx(hi).GetNeighbors() if n.GetAtomicNum() in (7, 8)), None)
        if parent is not None:
            lig_donor_h.append(hi)
    # 氢键受体（配体侧）：N/O/F 重原子
    lig_acceptors = [idx for idx in l_heavy
                     if lig.GetAtomWithIdx(idx).GetAtomicNum() in (7, 8, 9)
                     and lig.GetAtomWithIdx(idx).GetFormalCharge() <= 0]
    # 盐桥带电原子（配体侧）
    lig_charged = [idx for idx in l_heavy if abs(lig.GetAtomWithIdx(idx).GetFormalCharge()) >= 1]

    # --- 蛋白侧角色 ---
    prot_pos = [k for k, m in enumerate(heavy_meta)
                if m["name"] in POSITIVE_SIDECHAIN.get(m["resname"], set())]
    prot_neg = [k for k, m in enumerate(heavy_meta)
                if m["name"] in NEGATIVE_SIDECHAIN.get(m["resname"], set())]

    HB_DIST = 2.6  # H···受体距离上限（Å）

    # 1) 配体给体 H ↔ 蛋白受体（HBDonor：配体提供氢键）
    for hi in lig_donor_h:
        parent = next(n.GetIdx() for n in lig.GetAtomWithIdx(hi).GetNeighbors() if n.GetAtomicNum() != 1)
        d = np.linalg.norm(R - L[hi], axis=1)
        for k in np.argsort(d)[:30]:
            if d[k] > HB_DIST or heavy_meta[k]["element"] not in ("N", "O", "S"):
                break
            # 角度 D-H···A > 110°
            v1 = L[parent] - L[hi]; v2 = R[k] - L[hi]
            cosang = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            if np.degrees(np.arccos(np.clip(cosang, -1, 1))) > 110:
                add(res_label(heavy_meta[k]), "HBDonor",
                    f"H···A {d[k]:.2f} Å", [parent])
    # 2) 蛋白给体 H ↔ 配体受体（HBAcceptor：配体接受氢键）
    acc_set = set(lig_acceptors)
    for pi in recs["h_idx"]:
        parent_i = recs["h_parent"].get(pi)
        if parent_i is None or recs["recs"][parent_i]["element"] not in ("N", "O", "S"):
            continue
        d = np.linalg.norm(L_heavy - recs["recs"][pi]["xyz"], axis=1)
        for li in np.argsort(d)[:10]:
            atom_idx = l_heavy[li]  # L_heavy 下标 -> 配体原子 idx
            if d[li] > HB_DIST or atom_idx not in acc_set:
                continue
            v1 = recs["recs"][parent_i]["xyz"] - recs["recs"][pi]["xyz"]
            v2 = L_heavy[li] - recs["recs"][pi]["xyz"]
            cosang = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
            if np.degrees(np.arccos(np.clip(cosang, -1, 1))) > 110:
                add(res_label(recs["recs"][parent_i]), "HBAcceptor",
                    f"H···A {d[li]:.2f} Å", [atom_idx])
    # 3) 疏水接触：C···C ≤ 4.0 Å
    lig_c = [idx for idx in l_heavy if lig.GetAtomWithIdx(idx).GetAtomicNum() == 6]
    prot_c = [k for k, m in enumerate(heavy_meta) if m["element"] == "C"]
    if lig_c and prot_c:
        D = np.linalg.norm(L[lig_c][:, None, :] - R[prot_c][None, :, :], axis=2)
        ii, kk = np.where(D <= 4.0)
        for li, k in zip(ii, kk):
            add(res_label(heavy_meta[prot_c[k]]), "Hydrophobic",
                f"C···C {D[li, k]:.2f} Å", [lig_c[li]])
    # 4) 盐桥
    for li in lig_charged:
        q = lig.GetAtomWithIdx(li).GetFormalCharge()
        targets = prot_neg if q > 0 else prot_pos
        for k in targets:
            dist = np.linalg.norm(L[li] - R[k])
            if dist <= 4.0:
                add(res_label(heavy_meta[k]), "SaltBridge", f"{dist:.2f} Å", [li])
    # 5) π 堆积 / 阳离子-π
    centroids = _ring_centroids(lig, L)
    for label, cxyz in _protein_rings(recs):
        for c, ring_atoms in centroids:
            dist = float(np.linalg.norm(c - cxyz))
            if dist <= 5.5:
                add(label, "PiStacking", f"centroid {dist:.2f} Å", ring_atoms)
        if centroids and lig_charged:
            for li in lig_charged:
                c_near, atoms_near = min(centroids, key=lambda t: np.linalg.norm(L[li] - t[0]))
                if np.linalg.norm(L[li] - c_near) <= 5.0:
                    add(label, "PiCation", "cation-pi", atoms_near)
    return out


def _parse_res_label(label: str) -> tuple[str, str, str]:
    """'THR95.A' -> ('THR', 'A', '95')（残基名/链/编号）。"""
    resname = label[:3]
    resid, chain = label[3:].rsplit(".", 1)
    return resname, chain, resid


def plot_interactions_2d(pose_sdf: Path, receptor_h_pdb: Path, out_base: Path,
                         title: str = "", edges: list[dict] | None = None) -> dict[str, Path]:
    """Discovery Studio 风格 2D 相互作用图：配体键线式居中，残基彩色圆形节点环绕，
    彩色虚线从残基连到配体具体原子（PDF 矢量 + PNG）。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from rdkit import Chem
    from rdkit.Chem import AllChem

    out_base = Path(out_base)
    ensure_dir(out_base.parent)
    if edges is None:
        edges = detect_interactions(pose_sdf, receptor_h_pdb)

    lig = next((m for m in Chem.SDMolSupplier(str(pose_sdf), removeHs=False) if m), None)
    try:
        Chem.Kekulize(lig, clearAromaticFlags=True)  # 芳香键转单双交替，便于画键线式
    except Exception:
        pass
    AllChem.Compute2DCoords(lig)
    P = np.array(lig.GetConformer().GetPositions())[:, :2]
    span = max(P.max(0) - P.min(0)) / 2 or 1.0
    P = (P - (P.max(0) + P.min(0)) / 2) / span * 0.55  # 归一化到 [-0.55, 0.55]

    heavy = [a.GetIdx() for a in lig.GetAtoms() if a.GetAtomicNum() > 1]
    pos = {i: P[i] for i in heavy}
    ELEM_COLOR = {"N": "#2f4fd0", "O": "#d02020", "S": "#c8970a", "P": "#e07b00",
                  "F": "#3a9a3a", "Cl": "#3a9a3a", "Br": "#8a4a2a", "I": "#7a3aa0"}

    fig, ax = plt.subplots(figsize=(7.4, 5.6))
    ax.set_aspect("equal"); ax.axis("off")

    # ---- 键线式骨架（重原子，双键画平行双线）----
    for b in lig.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i not in pos or j not in pos:
            continue
        x1, y1 = pos[i]; x2, y2 = pos[j]
        d = np.array([x2 - x1, y2 - y1])
        n = d / (np.linalg.norm(d) or 1.0)
        t = np.array([-n[1], n[0]]) * 0.024
        if b.GetBondTypeAsDouble() == 2:
            inn = 0.14  # 内线两端各缩短 14%
            ax.plot([x1 + d[0] * inn, x2 - d[0] * inn], [y1 + d[1] * inn, y2 - d[1] * inn],
                    color="#2b2b2b", lw=1.1, zorder=2)
            ax.plot([x1 - t[0], x2 - t[0]], [y1 - t[1], y2 - t[1]], color="#2b2b2b", lw=1.1, zorder=2)
        else:
            ax.plot([x1, x2], [y1, y2], color="#2b2b2b", lw=1.3, zorder=2)
    # ---- 杂原子标签（碳不标）----
    for i in heavy:
        sym = lig.GetAtomWithIdx(i).GetSymbol()
        if sym != "C":
            ax.scatter([pos[i][0]], [pos[i][1]], s=110, color="white", zorder=3.4, edgecolors="none")
            ax.text(pos[i][0], pos[i][1], sym, ha="center", va="center", fontsize=9,
                    color=ELEM_COLOR.get(sym, "#333333"), zorder=4, fontweight="bold")

    # ---- 残基节点：按相互作用原子方向环绕排布，环向去重 ----
    mol_c = P[heavy].mean(axis=0)
    by_res: dict[str, list[dict]] = {}
    for e in edges:
        by_res.setdefault(e["residue"], []).append(e)
    items = []
    for res, es in by_res.items():
        idxs = [a for e in es for a in e.get("lig_atoms", []) if a in pos]
        v = (P[idxs].mean(axis=0) - mol_c) if idxs else np.array([0.0, -1.0])
        items.append([res, es, float(np.degrees(np.arctan2(v[1], v[0])))])
    items.sort(key=lambda t: t[2])
    # 配体包围盒半宽/半高：椭圆环绕与节点避让都用
    HW = float(P[heavy][:, 0].max()) + 0.06
    HH = float(P[heavy][:, 1].max()) + 0.06

    def _ring_r(ang):
        t = np.deg2rad(ang)
        c, s = abs(np.cos(t)), abs(np.sin(t))
        te = min(HW / c if c > 1e-9 else 1e9, HH / s if s > 1e-9 else 1e9)
        return te + 0.36                          # 节点半径 ~0.16 + 虚线呼吸段

    MIN_GAP = 360.0 / max(len(items), 18)  # 环向最小角间隔
    for a_i in range(len(items)):
        for b_i in range(a_i + 1, len(items)):
            if abs(items[b_i][2] - items[a_i][2]) < MIN_GAP:
                items[b_i][2] = items[a_i][2] + MIN_GAP

    # 椭圆上节点直径级避让：按环序逐个向前找第一个不与已放节点冲突的角度；
    # 一圈放不下时半径外扩 0.55 进下一圈（保证必然有解，绝不在失败位置硬放）
    def _node_xy(ang, lap=0):
        t = np.deg2rad(ang)
        rr = _ring_r(ang) + 0.55 * lap
        return np.array([np.cos(t) * rr, np.sin(t) * rr])

    placed_ang: list[tuple[float, int]] = []
    for it in items:
        done = False
        for lap in range(3):
            a = it[2]
            for _ in range(130):
                xy = _node_xy(a, lap)
                if all(float(np.hypot(*(_node_xy(pa, pl) - xy))) >= 0.38
                       for pa, pl in placed_ang):
                    done = True
                    break
                a += 3.0
            if done:
                break
        it[2] = a
        it.append(lap if done else 3)   # 圈数写入 items，摆放时用（否则第二圈仍画在第一圈半径上）
        placed_ang.append((a, lap if done else 3))

    PRIO = {"SaltBridge": 0, "HBDonor": 1, "HBAcceptor": 1,
            "PiStacking": 2, "PiCation": 2, "Hydrophobic": 3}
    # 节点椭圆环绕（DS 风格）：沿方向角求配体包围盒的射线交点 + 固定间隙，
    # 节点紧贴配体轮廓，虚线极短；比固定大圆半径大幅减少长虚线横穿画面
    _nx, _ny = [], []
    for res, es, ang, lap in items:
        rr = _ring_r(ang) + 0.55 * lap
        theta = np.deg2rad(ang)
        nx, ny = np.cos(theta) * rr, np.sin(theta) * rr
        _nx.append(nx); _ny.append(ny)
        # 同残基同类型合并为一条虚线（连到参与原子质心），避免多线汇聚标注重叠
        by_type: dict[str, list[int]] = {}
        for e in es:
            by_type.setdefault(e["type"], []).extend(
                a for a in e.get("lig_atoms", []) if a in pos)
        for typ, atoms in by_type.items():
            tx, ty = pos[atoms[0]] if len(atoms) == 1 else P[atoms].mean(axis=0)
            u = np.array([tx - nx, ty - ny]); u = u / (np.linalg.norm(u) or 1.0)
            ax.plot([nx + u[0] * 0.150, tx - u[0] * 0.05], [ny + u[1] * 0.150, ty - u[1] * 0.05],
                    ls=(0, (4, 2.4)), color=_inter_color(typ), lw=1.5, alpha=0.95, zorder=1.8)
        main = min((e["type"] for e in es), key=lambda t: PRIO.get(t, 9))
        face, tcol = NODE_COLORS.get(main, (_inter_color(main), "white"))
        ax.scatter([nx], [ny], s=720, color=face, zorder=4.6,
                   edgecolors="white", linewidths=1.4, alpha=0.98)
        resname, chain, resid = _parse_res_label(res)
        ax.text(nx, ny, f"{resname}\n{chain}:{resid}", ha="center", va="center",
                fontsize=7.2, color=tcol, fontweight="bold", zorder=5, linespacing=0.95)

    # ---- 画幅收紧到节点范围（消除大片空白）----
    pad = 0.24
    ax.set_xlim(min(_nx) - pad, max(_nx) + pad)
    ax.set_ylim(min(_ny) - pad, max(_ny) + pad)

    # ---- 图例（仅显示出现的类型）----
    seen = [c for c in INTERACTION_COLORS if any(e["type"] == c[0] for e in edges)]
    if seen:
        handles = [plt.Line2D([], [], color=c, lw=2.2, ls=(0, (4, 2.4)), label=k)
                   for k, c in seen]
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.01),
                  ncol=min(3, len(handles)), frameon=False, fontsize=9.6)
    if title:
        ax.set_title(title, fontsize=12, fontweight="bold", pad=10)
    fig.tight_layout()
    pdf = out_base.with_suffix(".pdf"); png = out_base.with_suffix(".png")
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {"pdf": pdf, "png": png}


# =====================================================================
# 受体加氢（供相互作用几何判据的氢键角度检测）
# =====================================================================
def prepare_receptor_h(receptor_pdb: Path, out_pdb: Path, ph: float = 7.4) -> Path:
    """pdbfixer + OpenMM 为受体加氢（pH 7.4）；失败退化为重原子（给体检测缺失）。"""
    out_pdb = Path(out_pdb)
    ensure_dir(out_pdb.parent)
    if out_pdb.exists():
        return out_pdb
    try:
        from openmm.app import PDBFile, Modeller
        from pdbfixer import PDBFixer

        fixer = PDBFixer(filename=str(receptor_pdb))
        fixer.findMissingResidues()
        fixer.missingResidues = {}          # 保守：不补建缺失环，只加氢
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(False)
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(float(ph))
        modeller = Modeller(fixer.topology, fixer.positions)
        with open(out_pdb, "w", encoding="utf-8") as fh:
            PDBFile.writeFile(modeller.topology, modeller.positions, fh, keepIds=True)
        log.info("受体加氢完成（pdbfixer pH %.1f）: %s", ph, out_pdb)
    except Exception as exc:
        log.warning("pdbfixer 加氢失败（%s: %s），退化为重原子受体", type(exc).__name__, exc)
        import MDAnalysis as mda

        u = mda.Universe(str(receptor_pdb))
        u.select_atoms("not name H*").write(str(out_pdb))
    return out_pdb


# =====================================================================
# 编排入口
# =====================================================================
def run_visualize(cfg: Config, tag: str, top: int, entries: list[str] | None = None) -> dict[str, Any]:
    """可视化指定批次：Top N 或用户指定化合物（entries），返回产物清单。

    输出为任务自包含目录 output_figure/{target}/{tag}/：
      pose_R{rank}_{药名}_3d.png   3D 口袋特写
      pose_R{rank}_{药名}_2d.pdf/.png  2D 相互作用图
      pose_R{rank}_{药名}_pair.png 3D|2D 横拼对照（文件名即对应关系）
      top_overview.png / MANIFEST.csv / {tag}_visualize_summary.json
    """
    assert cfg.target is not None
    ranking_csv = cfg.paths.output_data_dir / cfg.target.name / f"{tag}_ranking.csv"
    if not ranking_csv.exists():
        raise FileNotFoundError(f"排行榜不存在: {ranking_csv}（先执行 autovs dock）")
    if entries:
        poses = match_entries(ranking_csv, entries)
        if not poses:
            raise FileNotFoundError(f"指定的化合物在排行榜中无可用结果: {entries}")
    else:
        poses = load_top_poses(ranking_csv, top)
    if not poses:
        raise FileNotFoundError("排行榜中没有可用的 pose SDF")
    # 任务自包含目录：一个 target+tag 一个文件夹，ranking 图谱清单齐备，整包可拷走复用
    out_fig = ensure_dir(cfg.paths.output_figure_dir / cfg.target.name / tag)

    # 受体加氢（氢键角度判据）
    receptor_h = prepare_receptor_h(cfg.paths.temp_dir / "receptor" / cfg.target.name / "receptor_clean.pdb",
                                    ensure_dir(cfg.paths.temp_dir / "receptor" / cfg.target.name) / "receptor_h.pdb")
    # 相互作用检测一次，3D 标签与 2D 绘图共用；渲染前先抽单 pose SDF（防 PyMOL 多 state 加载失真）
    single_dir = ensure_dir(cfg.paths.temp_dir / "viz" / cfg.target.name)
    for p in poses:
        p["pose_sdf"] = str(extract_first_pose(p["pose_sdf"], single_dir / f"pose_R{p['rank']:02d}.sdf"))
        p["edges"] = detect_interactions(Path(p["pose_sdf"]), receptor_h)
    # 3D 标签集合：非疏水优先、疏水按最近距离补充，总数上限 7（相互作用多的配体防标签爆炸）
    resi_map = {}
    for p in poses:
        non_hyd = [e["residue"] for e in p["edges"] if e["type"] != "Hydrophobic"]
        hyd = sorted((float(e["detail"].split()[-2]), e["residue"])
                     for e in p["edges"] if e["type"] == "Hydrophobic")
        picked: list[str] = []
        for r in non_hyd + [r for _, r in hyd]:
            if len(picked) >= 7:
                break
            if r not in picked:
                picked.append(r)
        resi_map[p["rank"]] = sorted({_parse_res_label(r)[1:3] for r in picked})

    # 3D（PyMOL 批处理，内部用 pose_RXX.png 命名，渲染后重命名统一前缀）
    render_poses_3d(cfg, poses, out_fig, resi_map)
    for p in poses:
        src = out_fig / f"pose_R{p['rank']:02d}.png"
        if src.exists():
            src.replace(out_fig / f"{_pose_label(p)}_3d.png")
    pngs = sorted(out_fig.glob("pose_R*_3d.png")) + ([out_fig / "top_overview.png"]
                                                     if (out_fig / "top_overview.png").exists() else [])

    # 2D（Discovery Studio 风格）+ 3D|2D 横拼对照
    items = []
    for p in poses:
        base = out_fig / f"{_pose_label(p)}_2d"
        png2d = out_fig / f"{_pose_label(p)}_2d.png"
        png3d = out_fig / f"{_pose_label(p)}_3d.png"
        paths = plot_interactions_2d(Path(p["pose_sdf"]), receptor_h, base,
                                     title=f"{p['name'] or p['entry_id']}  (rank {p['rank']})",
                                     edges=p["edges"])
        pair = _pair_image(png3d, png2d, out_fig / f"{_pose_label(p)}_pair.png") if png3d.exists() else ""
        items.append({"rank": p["rank"], "entry_id": p["entry_id"], "name": p["name"],
                      "score": p.get("score", ""), "pair": str(pair),
                      "interactions": [f"{e['residue']}:{e['type']}" for e in p["edges"]],
                      **{k: str(v) for k, v in paths.items()}})

    # MANIFEST：任务产物索引（合并式：--entry 增量可视化时更新对应行而非覆盖历史）
    mf_path = out_fig / "MANIFEST.csv"
    merged: dict[str, dict[str, str]] = {}
    if mf_path.exists():
        with mf_path.open(encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                merged[row["entry_id"]] = row
    for it in items:
        merged[it["entry_id"]] = {
            "rank": it["rank"], "entry_id": it["entry_id"], "name": it["name"],
            "score_kcal_mol": it["score"],
            "fig_3d": str(out_fig / f"{_pose_label(it)}_3d.png"),
            "fig_2d_png": it["png"], "fig_2d_pdf": it["pdf"], "fig_pair": it["pair"],
            "interactions": "; ".join(it["interactions"]),
        }
    with mf_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "entry_id", "name", "score_kcal_mol", "fig_3d", "fig_2d_png",
                    "fig_2d_pdf", "fig_pair", "interactions"])
        for row in sorted(merged.values(), key=lambda r: int(r["rank"])):
            w.writerow([row[c] for c in ["rank", "entry_id", "name", "score_kcal_mol",
                                         "fig_3d", "fig_2d_png", "fig_2d_pdf", "fig_pair", "interactions"]])

    summary = {"target": cfg.target.name, "tag": tag, "n_poses": len(poses),
               "3d_png": [str(p) for p in pngs], "2d_items": items}
    from .utils import write_json
    write_json(out_fig / f"{tag}_visualize_summary.json", summary)
    log.info("可视化完成: %d 组（3D+2D+对照）-> %s", len(items), out_fig)
    return summary
