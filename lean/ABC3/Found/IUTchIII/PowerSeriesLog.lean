import Mathlib.RingTheory.PowerSeries.Log

/-!
# 形式的べき級数の対数の加法公式

mathlib の `PowerSeries.log` / `PowerSeries.logOf`
(`Mathlib/RingTheory/PowerSeries/Log.lean`)には
**加法公式 `logOf (f * g) = logOf f + logOf g` が無い**(2026-08-14 実測)。
p進対数の乗法性(`Found/IUTchIII/PadicLog.lean`)を出すのに要るので、ここで補う。

## 証明の筋(mathlib の `exp_unique_of_derivative_eq_self` と同じ形)

1. `deriv_log : d⁄dX (log A) = mk fun n ↦ (-1)^n` ——**mathlib にある**。
2. その幾何級数は `1 + X` の逆元(`geom_mul_one_add`)。
3. 代入は環準同型なので、`constantCoeff f = 1` のとき
   `(d⁄dX (log A)).subst (f - 1) * f = 1`(`substDerivLog_mul`)。
4. ゆえに `d⁄dX (logOf f) * f = d⁄dX f`(`deriv_logOf_mul_self`)——対数微分。
5. `f`, `g` の定数項が 1 なら `f * g` は**単元**なので、両辺に掛けて比較でき、
   `derivative.ext`(微分と定数項が一致すれば等しい)で結論。

## ★★測定の訂正(2026-08-14)

前回「**`PowerSeries.log` は mathlib に無い**」と報告したが、**誤りだった**。
`Mathlib/RingTheory/PowerSeries/Log.lean`(101 行、2026 年追加)に
`log` / `coeff_log` / `deriv_log` / `logOf` / `HasSubst.log` などが**ある**。
無いのは**加法公式だけ**だった。

誤りの原因は探索手順: `WellKnown.lean` の中しか見ずに「`PowerSeries.exp` はあるが
`log` は無い」と結論した。`ls Mathlib/RingTheory/PowerSeries/` で
**ディレクトリを列挙していれば `Log.lean` が目に入った**。
`tools/check.mjs` 冒頭 A2 の (S2)「全出現を列挙する」を、
**ファイル名の列挙にも適用すべきだった**。

## mathlib への貢献としての位置づけ

本ファイルの内容は p進とは独立で、`Mathlib/RingTheory/PowerSeries/Log.lean` に
そのまま足せる形にしてある。

## ★隔たりは埋まった(2026-08-15)。ただし道は違った。

`logOf_mul` は**形式的べき級数の等式**であり、これを ℚ_p の
`Found/IUTchIII/PadicLog.lean` の `logOneAdd` に移すには**評価の橋**が要る。
当初は次の2つを想定していたが、**どちらも使わなかった**:

```
-- (i) mathlib の eval₂ での評価と我々の tsum 定義の一致  → 下記の理由で不可
-- (ii) 2変数版の評価(ℚ⟦X,Y⟧ で立てて (x,y) で評価)  → MvPowerSeries に微分が無く不可
```

実際に通した道は `Found/IUTchIII/PadicLogMul.lean` にある。
**帳簿用の変数を 1 つだけ入れる**:

- `ℚ_p⟦X⟧` の中で `(1 + C x·X)(1 + C y·X) = 1 + C(x+y)·X + C(xy)·X²` とし、
  本ファイルの `logOf_mul` を当てる。`x`, `y` は**係数**であって
  代入する点では無いので、形式的な側は最後まで `ℚ_p⟦X⟧`(X進位相=線形位相)に留まる。
- ℚ_p 側では、係数の総和 `∑' n, coeff n (logOf (1 + P))` が
  `logOneAdd (P.eval 1)` に等しいことを**手で**示す
  (`summable_and_tsum_coeff_logOf`)。二重級数の並べ替えであり、
  分母は `1/d` ただ一つなので `‖coeff d (log)‖ ≤ d` で押さえられる。

結果として `ABC3.Found.IUTchIII.logOneAdd_mul` が `sorry` 無しである。

### ★★実測: mathlib の評価 API は **ℚ_p には使えない**

`MvPowerSeries.eval₂Hom` / `aeval` は

```
[IsTopologicalSemiring R] [IsUniformAddGroup R] [IsUniformAddGroup S]
[CompleteSpace S] [T2Space S] [IsTopologicalRing S] [IsLinearTopology S S]
```

を要求する。**`IsLinearTopology ℚ_[p] ℚ_[p]` はインスタンスとして存在しない**
(2026-08-15 実測: `infer_instance` が失敗)。数学的にも当然で、線形位相は
**イデアル**からなる 0 の近傍基を要求するが、ℚ_p は体なのでイデアルは `0` と全体しかなく、
ノルム位相は離散でも密着でもない。mathlib の評価 API は
**`ℚ⟦Y⟧` のような進位相の環のために作られており、ノルム体には適用できない**。

`Mathlib/RingTheory/PowerSeries/Restricted.lean` は係数が有界なべき級数
(`IsRestricted`)を扱うが、**評価・総和の理論は無い**(全宣言を列挙して確認:
`zero` / `one` / `monomial` / `C` / `add` / `neg` / `smul` / `convergenceSet` / `mul` のみ)。

### ★★実測: `MvPowerSeries` には微分が無い

`ls Mathlib/RingTheory/MvPowerSeries/` は全 14 ファイル(Basic, Equiv, Evaluation, Expand,
GaussNorm, Inverse, LexOrder, LinearTopology, NoZeroDivisors, Order, PiTopology, Rename,
Substitution, Trunc)。**`Derivative.lean` は無く**、
`grep -rniE "derivativ|derivation|d⁄d"` は全ファイルで **0 件**。
ゆえに `deriv_log` / `derivative.ext` に相当するものも無く、
`Mathlib/RingTheory/PowerSeries/Log.lean` は丸ごと `namespace PowerSeries` なので
`MvPowerSeries.log` も無い。**2変数の形式的加法公式を取るには微分から作る必要がある**。
上記の1変数の道は、それを回避するために採った。

### 検討して却下した案(2026-08-15)

`ℚ⟦Y⟧` を係数環とする 1 変数版 `ℚ⟦Y⟧⟦X⟧` で `f := 1 + X`, `g := 1 + C Y` と置く案は、
**`constantCoeff g = 1 + Y ≠ 1`** なので `logOf_mul` の仮定を満たさない
(Lean で確認済み: `constantCoeff (1 + C X) = 1` は `False` に帰着する)。
`C a` の定数項は `a` であって 0 ではないため、`HasSubst (g - 1)` も成り立たない。
-/

namespace ABC3.Found.IUTchIII

open PowerSeries

variable {A : Type*} [CommRing A] [Algebra ℚ A]

/-- 幾何級数 `Σ (-1)ⁿ Xⁿ` は `1 + X` の逆元。 -/
theorem geom_mul_one_add :
    (mk fun n => algebraMap ℚ A ((-1 : ℚ) ^ n)) * (1 + X) = 1 := by
  ext n
  rw [mul_add, mul_one, map_add]
  cases n with
  | zero => simp
  | succ m =>
    rw [coeff_succ_mul_X, coeff_mk, coeff_mk, ← map_add, coeff_one]
    simp [pow_succ]

omit [Algebra ℚ A] in
/-- `constantCoeff f = 1` なら `f - 1` は代入可能。 -/
theorem hasSubst_sub_one {f : A⟦X⟧} (hf : constantCoeff f = 1) : HasSubst (f - 1) :=
  HasSubst.of_constantCoeff_zero' (by rw [map_sub, map_one, hf, sub_self])

/-- `log(1+X)` の微分を `f - 1` で代入したものは `f` の逆元。 -/
theorem substDerivLog_mul {f : A⟦X⟧} (hf : constantCoeff f = 1) :
    (d⁄dX A (log A)).subst (f - 1) * f = 1 := by
  have hsub : HasSubst (f - 1) := hasSubst_sub_one hf
  have h1 : (d⁄dX A (log A)) * (1 + X) = 1 := by
    rw [deriv_log]; exact geom_mul_one_add
  have h2 := congrArg (substAlgHom hsub) h1
  simp only [map_mul, map_one, map_add, coe_substAlgHom, subst_X hsub] at h2
  rwa [show (1 : A⟦X⟧) + (f - 1) = f by ring] at h2

/-- ★**対数微分**: `d⁄dX (logOf f) * f = d⁄dX f`。 -/
theorem deriv_logOf_mul_self {f : A⟦X⟧} (hf : constantCoeff f = 1) :
    d⁄dX A (logOf f) * f = d⁄dX A f := by
  have hsub : HasSubst (f - 1) := hasSubst_sub_one hf
  rw [logOf_eq, derivative_subst A hsub, map_sub, Derivation.map_one_eq_zero, sub_zero]
  calc (d⁄dX A (log A)).subst (f - 1) * d⁄dX A f * f
      = ((d⁄dX A (log A)).subst (f - 1) * f) * d⁄dX A f := by ring
    _ = d⁄dX A f := by rw [substDerivLog_mul hf, one_mul]

omit [Algebra ℚ A] in
/-- 定数項が 1 のべき級数は単元。 -/
theorem isUnit_of_constantCoeff_one {f : A⟦X⟧} (hf : constantCoeff f = 1) : IsUnit f := by
  rw [PowerSeries.isUnit_iff_constantCoeff, hf]
  exact isUnit_one

/-- ★★**加法公式** `logOf (f * g) = logOf f + logOf g`。

mathlib の `PowerSeries.log` に欠けていたもの。 -/
theorem logOf_mul [IsAddTorsionFree A] {f g : A⟦X⟧}
    (hf : constantCoeff f = 1) (hg : constantCoeff g = 1) :
    logOf (f * g) = logOf f + logOf g := by
  have hfg : constantCoeff (f * g) = 1 := by rw [map_mul, hf, hg, one_mul]
  refine derivative.ext ?_ ?_
  · -- 両辺の微分が一致することを、単元 `f * g` を掛けて比較する
    have hunit : IsUnit (f * g) := isUnit_of_constantCoeff_one hfg
    refine (IsUnit.mul_left_inj hunit).mp ?_
    · have e1 := deriv_logOf_mul_self hfg
      have e2 := deriv_logOf_mul_self hf
      have e3 := deriv_logOf_mul_self hg
      rw [Derivation.leibniz, smul_eq_mul, smul_eq_mul] at e1
      calc d⁄dX A (logOf (f * g)) * (f * g)
          = f * d⁄dX A g + g * d⁄dX A f := e1
        _ = f * (d⁄dX A (logOf g) * g) + g * (d⁄dX A (logOf f) * f) := by rw [e2, e3]
        _ = (d⁄dX A (logOf f) + d⁄dX A (logOf g)) * (f * g) := by ring
        _ = d⁄dX A (logOf f + logOf g) * (f * g) := by rw [map_add]
  · rw [constantCoeff_logOf hfg, map_add, constantCoeff_logOf hf, constantCoeff_logOf hg, add_zero]

end ABC3.Found.IUTchIII
