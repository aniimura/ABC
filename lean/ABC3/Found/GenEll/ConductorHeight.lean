import ABC3.Found.GenEll.HeightNonneg

/-!
# [GenEll] Proposition 1.6 —— 導手は**構成した高さ**で抑えられる(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.6. (Conductor Bounded by the Height) Let D ⊆ X be an effective Cartier divisor,

## ★★★`CartierPullback.lean` が名指しした穴を塞ぐ

`CartierPullback.lean` の `logCond_le_degNormalized_pullback` の docstring はこう書いた:

> ★★原文の右辺は「高さ」だが、その手前の段——「導手 ≤ 引き戻した因子の次数」——が
> これである。**高さとの接続には `X` 上の算術直線束(`Definition 1.1`、
> 複素解析空間が要るので未着手)が要る。混同しない。**

★★★**本日 `htArith` を構成したので、その接続がつく。**

## ★★機構は 2 段

1. `logCond ≤ deg_F(D_x)` —— 被約化で次数は減る(`logCond_le_degNormalized_pullback`)
2. `deg_F(D_x) ≤ ht` —— ★**アルキメデス側の寄与が非負**だから

★★★**2 段目が本ファイルの内容である。**
`ht` は有限素点側(= `deg_F(D_x)`)とアルキメデス側(= Green 関数の値の和)の和なので、
Green 関数が非負なら有限素点側だけを取ったものより大きい。

## ★原文の `≲` より強い形が出る

原文は BD-class の不等式 `log-cond_D ≲ ht` を主張する。
★★構成の側では **`log-cond_D(x) ≤ ht(x)`(定数を要しない真の不等式)**が出る。
★これは `Proposition 1.4, (ii)` のときと同じ現象である
(`HeightNonneg.lean` を参照)。

## ★★これは `Proposition 1.6` 全体ではない

原文の `Proposition 1.6` は `X` の正規性・固有性や
`L = O_X(D)` との対応を含む。★本ファイルが取るのは
**「導手 ≤ 構成した高さ」という不等式そのもの**である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

variable (F : Type) [Field F] [NumberField F]

/-! ## ★引き戻しの次数は 2 成分に分かれる -/

/-- ★**引き戻した算術因子の次数** = 有限素点側 + アルキメデス側。

★`idealADiv` はアルキメデス側が `0` なので、その次数は有限素点側だけである。 -/
theorem deg_pullbackADiv {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) :
    deg (pullbackADiv F D xF)
      = deg (idealADiv F (pullbackIdeal F D.divisor xF))
        + (archADiv F D.green xF).sum (fun _ r => r) := by
  have harc : (idealADiv F (pullbackIdeal F D.divisor xF)).arc = 0 :=
    idealADiv_arc F _
  rw [deg, deg, pullbackADiv_fin, pullbackADiv_arc, harc, Finsupp.sum_zero_index]
  ring

/-- ★**アルキメデス側の寄与は非負**(Green 関数が非負なら)。 -/
theorem archADiv_sum_nonneg {X : Scheme.{0}} (g : GreenFn X)
    (xF : specRingOfIntegers F ⟶ X) (hg : ∀ p, 0 ≤ g p) :
    0 ≤ (archADiv F g xF).sum (fun _ r => r) := by
  rw [Finsupp.sum]
  refine Finset.sum_nonneg fun v _ => ?_
  rw [archADiv_apply]
  exact mul_nonneg (Nat.cast_nonneg _) (hg _)

/-! ## ★★★`Proposition 1.6` の不等式 -/

/-- ★★**有限素点側だけを取ったものは、高さ以下である**。

★★アルキメデス側の寄与が非負だからである。 -/
theorem degNormalized_pullbackIdeal_le_htArith {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X) (hg : ∀ p, 0 ≤ D.green p) :
    degNormalized (idealADiv F (pullbackIdeal F D.divisor xF)) ≤ htArith F D xF := by
  have hpos : (0 : ℝ) < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  rw [htArith, degNormalized, degNormalized, deg_pullbackADiv]
  exact div_le_div_of_nonneg_right
    (by linarith [archADiv_sum_nonneg F D.green xF hg]) hpos.le

/-- ★★★**[GenEll] Proposition 1.6 の不等式** —— 導手は構成した高さで抑えられる。

原文 (GenEll p.9):
> Proposition 1.6. (Conductor Bounded by the Height) Let D ⊆ X be an effective Cartier divisor,

★★★**原文の `≲`(BD-class、定数差を許す)より強い。**
定数を要しない真の不等式 `log-cond_D(x) ≤ ht(x)` が出る。

★★機構は 2 段:
1. 被約化で次数は減る(`logCond_le_degNormalized_pullback`)
2. アルキメデス側の寄与が非負(本ファイル)

★`CartierPullback.lean` が「高さとの接続には `Definition 1.1` が要る」と
名指しした穴を、本日構成した `htArith` で塞いだものである。 -/
theorem logCond_le_htArith {X : Scheme.{0}} (D : ArithCartier X)
    (xF : specRingOfIntegers F ⟶ X)
    (h : pullbackIdeal F D.divisor xF ≠ 0) (hg : ∀ p, 0 ≤ D.green p) :
    logCond F D.divisor xF ≤ htArith F D xF :=
  le_trans (logCond_le_degNormalized_pullback F D.divisor xF h)
    (degNormalized_pullbackIdeal_le_htArith F D xF hg)

/-- ★**BD-class の水準でも成り立つ**(原文の `≲` の形)。

★定数 `C = 0` で足りる——上で真の不等式が出ているからである。 -/
theorem logCond_bdle_htArith {X : Scheme.{0}} (D : ArithCartier X)
    (h : ∀ xF : specRingOfIntegers F ⟶ X, pullbackIdeal F D.divisor xF ≠ 0)
    (hg : ∀ p, 0 ≤ D.green p) :
    BDge (fun xF : specRingOfIntegers F ⟶ X => logCond F D.divisor xF)
      (fun xF => htArith F D xF) :=
  ⟨0, fun xF => by
    show logCond F D.divisor xF - htArith F D xF ≤ 0
    linarith [logCond_le_htArith F D xF (h xF) hg]⟩

/-! ## ★出典の紐付け(`.src`)

★条つきである。原文の `Proposition 1.6` は `X` の正規性・固有性や
`L = O_X(D)` との対応を含み、本ファイルが取るのは不等式そのものである。 -/

def logCond_le_htArith.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.6(構成した高さとの不等式のみ——L = O_X(D) の対応は含まない)",
    sectionId := "genell-prop-1-6" }

end ABC3.Found.GenEll
