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

## ★残っている隔たり(未証明。型のレベルで書く)

`logOf_mul` は**形式的べき級数の等式**であり、これを ℚ_p の
`Found/IUTchIII/PadicLog.lean` の `logOneAdd` に移すには、**評価の橋**が要る。
必要な補題は次の2つで、いずれも**未証明**である:

```
-- (i) 収束域での「べき級数の評価」と我々の tsum 定義の一致
theorem logOneAdd_eq_eval {x : ℚ_[p]} (hx : ‖x‖ < 1) :
    logOneAdd x = PowerSeries.eval₂ (algebraMap ℚ ℚ_[p]) x (PowerSeries.log ℚ)

-- (ii) 2変数版の評価。logOf_mul を ℚ⟦X,Y⟧ で立て、(x, y) で評価する
theorem logOneAdd_mul {x y : ℚ_[p]} (hx : ‖x‖ < 1) (hy : ‖y‖ < 1) :
    logOneAdd (x + y + x * y) = logOneAdd x + logOneAdd y
```

道具は mathlib にある——`PowerSeries.eval₂_eq_tsum (hφ : Continuous φ) (ha : HasEval a) :
eval₂ φ a f = ∑' d, φ (coeff d f) * a ^ d` と、代数準同型 `PowerSeries.aeval`。
★ただし (ii) には**1変数では届かない**: `logOf_mul` の `f`, `g` は定数項 1 の
**べき級数**であって定数ではないので、`f := C (1+x)` とは置けない
(`constantCoeff (C a) = a ≠ 1`)。`MvPowerSeries` の2変数版か、
`ℚ⟦Y⟧` を係数環とする1変数版を経由する必要がある。**そこは書いていない。**
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
