import ABC3.Found.GenEll.HeightConstruction
import ABC3.Found.GenEll.Conductor
import ABC3.Found.GenEll.BDClass

/-!
# [GenEll] Proposition 1.4, (ii) —— 構成した高さは非負である(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

## ★★★原文の仮定は、因子で表すと**アルキメデス側にしか効かない**

原文 (ii) は

> (ii) If some positive tensor power of the line bundle LQ on XQ is generated
> by global sections [for instance, if the line bundle LQ is ample], then

を仮定して `ht_L̄ ≳ 0` を結論する。

★★**因子で表すと事情が変わる。** `ArithCartier` の有限素点側は
**イデアル層**(`O_X` の部分層)であり、引き戻しても
`𝓞_F` の**イデアル**であって分数イデアルにはならない。
したがって `idealADiv` の係数は**常に非負**である。

★★★**つまり原文の「大域切断で生成される」という仮定は、
因子表示では「有効因子で代表できる」ことに吸収されており、
残る仮定はグリーン関数の非負性だけである。**

## ★★★しかも結論は `≳` より強く `≥` が出る

`Definition 1.2, (ii)` の `≳` は印字の向きが食い違っている
(`BDClass.lean` / `definition-1-2.html` に記録した)。
★★本ファイルが出すのは **`0 ≤ ht`(定数 `C` を要しない真の不等式)**であり、
`BDle` と `BDge` のどちらの読みよりも強い。★向きの争点が消える。

## ★ここには `PullbackMul` が要らない

`Proposition 1.4, (i)`(加法性)は引き戻しの積の保存を要したが、
**(ii) は 1 つの因子しか見ないので要らない**。★無条件である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F : Type) [Field F] [NumberField F]

/-! ## ★有効な算術因子の次数は非負 -/

/-- ★**有効な算術因子の次数は非負**。

有限素点側は係数 `≥ 0` かつ `log q_v > 0`、アルキメデス側は係数そのもの。 -/
theorem deg_nonneg_of_isEffective (a : ADiv F) (ha : a.IsEffective) : 0 ≤ deg a := by
  have h1 : 0 ≤ a.fin.sum (fun v n => (n : ℝ) * Real.log (residueCard v)) := by
    rw [Finsupp.sum]
    refine Finset.sum_nonneg fun v _ => ?_
    exact mul_nonneg (by exact_mod_cast ha.1 v) (log_residueCard_pos v).le
  have h2 : 0 ≤ a.arc.sum (fun _ r => r) := by
    rw [Finsupp.sum]
    exact Finset.sum_nonneg fun v _ => ha.2 v
  rw [deg]
  linarith

/-- ★**正規化しても非負**(`[F : ℚ] > 0` なので符号は変わらない)。 -/
theorem degNormalized_nonneg_of_isEffective (a : ADiv F) (ha : a.IsEffective) :
    0 ≤ degNormalized a := by
  rw [degNormalized]
  exact div_nonneg (deg_nonneg_of_isEffective F a ha) (Nat.cast_nonneg _)

/-! ## ★★引き戻した算術因子の有効性 -/

/-- ★★**引き戻した算術因子は有効である** —— グリーン関数が非負でありさえすれば。

★★**有限素点側には仮定が要らない。**
`pullbackIdeal` は `𝓞_F` の**イデアル**であり、`idealADiv` の係数は
素因子の重複度(自然数)だからである。 -/
theorem pullbackADiv_isEffective {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) (hg : ∀ p, 0 ≤ D.green p) :
    (pullbackADiv F D xF).IsEffective := by
  constructor
  · intro v
    rw [pullbackADiv_fin]
    exact (idealADiv_isEffective F (pullbackIdeal F D.divisor xF)).1 v
  · intro v
    rw [pullbackADiv_arc, archADiv_apply]
    exact hg _

/-! ## ★★★`Proposition 1.4, (ii)` -/

/-- ★★★**構成した高さは非負である**。

原文 (GenEll p.6):
> Proposition 1.4. (Basic Properties of Heights) In the notation of the above

★★★**原文の `≳ 0` より強い。** 定数 `C` を要しない真の不等式 `0 ≤ ht` である。
`Definition 1.2, (ii)` の `≳` の向きの食い違い(`BDClass.lean` に記録)は、
ここでは争点にならない。

★★**仮定はグリーン関数の非負性だけ**である——
有限素点側は構成から自動的に非負になる。
★`PullbackMul` は要らない((i) と違って因子を 1 つしか見ないため)。 -/
theorem htArith_nonneg {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) (hg : ∀ p, 0 ≤ D.green p) :
    0 ≤ htArith F D xF :=
  degNormalized_nonneg_of_isEffective F _ (pullbackADiv_isEffective F D xF hg)

/-- ★**BD-class の水準でも非負** —— 原文の `≳ 0` の形。

★`BDge α β ≝ ∃ C, ∀ x, α x − β x ≤ C`(印字どおりの向き)。
`0 ≤ ht` から `C = 0` で出る。 -/
theorem htArith_bdge_zero {X : Scheme.{0}} (D : ArithCartier X)
    (hg : ∀ p, 0 ≤ D.green p) :
    BDge (fun _ : specRingOfIntegers F ⟶ X => (0 : ℝ)) (fun xF => htArith F D xF) :=
  ⟨0, fun xF => by simpa using htArith_nonneg F D xF hg⟩

/-! ## ★出典の紐付け(`.src`)

★条つきである。`Proposition 1.4` 全体には (i)〜(iv) の 4 条すべてが要る。 -/

def htArith_nonneg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(構成した高さの非負性——因子表示では有限素点側は自動)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
