# -*- coding: utf-8 -*-
"""厳密整数演算による `hform`（第 1 跳びの損失）の検算。

★NormalizedTraceDescent.lean:65-66 の表:
    測り方              使う跳び  指数        値          axDecay 3 2 = 3^{1/6} との比較
    跡(上の層)          i = 8     2*8/18=8/9  3^{0.889}   超える
    sharp(上の層)       i = 8     8/18 = 4/9  3^{0.444}   超える
    ★sharp(第 1 跳び)  i₁ = 2    2/18 = 1/9  3^{0.111}   ★収まる
  ★この表は本日まで再現していない。

★★測る 2 点:
  (1) 表の跳び i = 8 と i₁ = 2 が正しいか（zeta-tower-check.py で既に確認済み。再掲）。
  (2) ★★**「収まる」は予算の話であって、降下がその損失を達成できるかは別**である。
      前々波で回復した x（WildDescentDistanceOnly の反証候補）について、
      第 1 跳びの損失 3^{i₁/e_L} = 3^{2/18} を達成できるかを測る:
        要求は v_M(x − x') >= v(eps) − 2 = 11 − 2 = 9
        実際は d(x, E₁) = 7
      ⇒ 達成できるか?
"""

import importlib.util
import sys

Z = r"D:/Math_ABC3/tools/zeta27-distance-check.py"
spec = importlib.util.spec_from_file_location("z27c", Z)
z = importlib.util.module_from_spec(spec)
sys.modules["z27c"] = z
spec.loader.exec_module(z)

T = r"D:/Math_ABC3/tools/zeta-tower-check.py"
spec2 = importlib.util.spec_from_file_location("ztow2", T)
zt = importlib.util.module_from_spec(spec2)
sys.modules["ztow2"] = zt
spec2.loader.exec_module(zt)

HK = [a for a in range(1, 27) if a % 3 == 1 and a != 1]


def main():
    # --- (1) 跳びの再掲 ---
    F = zt.Field(3, 3)
    iv = {a: F.v(F.sub(F.sigma_pi(a), F.PI)) for a in F.units()}
    subK = [a for a in F.units() if a % 3 == 1 and a != 1]
    subE = [a for a in F.units() if a % 9 == 1 and a != 1]
    i1 = min(iv[a] for a in subK) - 1
    i_top = min(iv[a] for a in subE) - 1
    eL = F.N
    print("=== (1) NormalizedTraceDescent.lean:65-66 の表の跳び ===")
    print(f"  e_L = {eL}")
    print(f"  ★第 1 跳び i1      = {i1}   (表の字面 2)")
    print(f"  上の層の跳び i      = {i_top}   (表の字面 8)")
    print(f"  跡(上の層) 指数     = (3-1)*{i_top}/{eL} = {2*i_top}/{eL}"
          f"  (表の字面 8/9 = {16}/18)")
    print(f"  sharp(上の層) 指数  = {i_top}/{eL}       (表の字面 4/9 = {8}/18)")
    print(f"  ★sharp(第1跳び)    = {i1}/{eL}        (表の字面 1/9 = {2}/18)")
    print(f"  axDecay 3 2         = 3^(1/6) = 3^({3}/18)")
    print(f"  ⇒ 跡 {2*i_top} > 3 : {2*i_top > 3} /  sharp(上) {i_top} > 3 : {i_top > 3}"
          f" /  ★sharp(第1) {i1} <= 3 : {i1 <= 3}")
    print()

    # --- (2) 達成できるか ---
    x = [71, 78, 9, 50, 30, 54, 67, 43, 40, 3, 25, 11, 40, 44, 74, 63, 38, 24]
    eps = min(z.vM(z.sub(z.apply_sigma(a, x), x)) for a in HK)
    dE = z.dist_to_F(x)
    dK = z.dist_to_K(x)
    print("=== (2) 第 1 跳びの損失を降下が達成できるか（回復した x で） ===")
    print(f"  eps = {eps},  d(x, E1) = {dE},  d(x, K) = {dK}")
    for name, jump in (("第 1 跳び i1", i1), ("axDecay 3 2", 3), ("上の層 i", i_top)):
        need = eps - jump
        print(f"  {name:14s}: 損失 3^({jump}/18) ⇒ 要求 v(x-x') >= {need}"
              f" ;  実際 d(x,E1) = {dE}  ⇒ 達成 {dE >= need}")
    print()
    print(f"  ★実際に達成できる損失の指数 = eps - d(x,E1) = {eps} - {dE} = {eps - dE}")
    print(f"  ★第 1 跳びの上界 {i1} を {eps - dE - i1} だけ超えている")
    print(f"  ★axDecay 3 2 の {3} を {eps - dE - 3} だけ超えている")
    print()
    print("=== (3) それでも対の予算は通る（BudgetFiniteExcess と同じ理由） ===")
    print(f"  対の予算 12 に対し 必要量 = eps - d(x,K) = {eps} - {dK} = {eps - dK}"
          f"  ⇒ slack {12 - (eps - dK)}")


if __name__ == "__main__":
    main()
