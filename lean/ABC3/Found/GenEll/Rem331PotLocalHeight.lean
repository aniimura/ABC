/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.LocalHeightDelta
import ABC3.Found.GenEll.LocalHeightRamified

/-!
# [GenEll] Remark 3.3.1 —— **潜在的局所高さは `L` の取り方に依らない**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.16。

原文 (GenEll p.16):
> that this definition is independent of the choice of L].

## ★★★詰まりは「何で測るか」だった

`LocalHeightRamified.lean` は `ord_w(x)/e = ord_v(x)`（`ordAt_div_ramificationIdx`）を
**`x : Kˣ` について**証明していた。★しかし局所高さは `v_K(q_E)` であり、
**潜在的乗法還元では `q_E` は `K` に無い**（拡大 `L` の上にしかない）。
★★だからこの補題は当たらない——そう見えていた。

★★★**当たる。`j`-不変量で測ればよい。**

原文 `Definition 3.3` の局所高さは極小判別式の付値である
（`LocalHeightDelta.lean` の `localHeight_eq_vAdd_Delta`）。そして

    j = c₄³/Δ,   乗法的還元では v(c₄) = 0

なので `v(j) = −v(Δ) = −（局所高さ）`。★★★★**`j` は `K` に在る**
（曲線の不変量であってモデルの取り方に依らない）ので、
`ordAt_div_ramificationIdx_indep` が **`x = j` でそのまま当たる**。

## ★★★★機構

| 段 | 宣言 |
|---|---|
| 乗法的還元なら `v(c₄) = 0` | ★`vAdd_c4_eq_zero`（mathlib の `HasMultiplicativeReduction` の欄そのもの） |
| 局所高さ `= −v(j)` | ★★`localHeight_eq_neg_vAdd_j` |
| `ordAt` と `vAdd(tateDvrVal)` は同じもの | ★`ordAt_eq_vAdd` |
| `ord_w(x)/e` は `L` に依らない | `ordAt_div_ramificationIdx_indep`（`LocalHeightRamified.lean`） |
| ★**まとめ** | ★★★`potLocalHeight_indep`（本ファイル） |

★★`v(c₄) = 0` は mathlib の `WeierstrassCurve.HasMultiplicativeReduction` が
**クラスの欄として持っている**（`multiplicativeReduction : v(c₄) = 1`、乗法的表記）。
★探しに行く必要が無かった——**仮定の中に既に在った**。

## ★★逸脱の記録（CLAUDE.md の「逸脱」）

★**原文は `E ⊗_K L` の局所高さを `e(L/K)` で割る**と書くが、本ファイルは
「`L` 上の極小 Weierstrass モデル `W_L` で、その `j` が `K` の `j` の像であるもの」
として受ける。★★`E ⊗_K L` の極小モデルが実際にそうなることは
`j` がモデルに依らない不変量であることから従うが、**その段は含んでいない**。

★★★**`A, B` は離散付値環**として受ける（原文の局所体 `K` とその有限拡大 `L`）。
`ordAt` は Dedekind 環の高さ 1 素点で定義されているが、DVR ではそれが極大イデアル
ただ 1 つなので一致する。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve ABC3.Found.GaloisRep

/-! ## ★★★★★乗法的還元では `c₄` は単元 -/

/-- ★★★★★**乗法的還元なら `v(c₄) = 0`**。

原文 (GenEll p.16):
> that this definition is independent of the choice of L].

★mathlib の `WeierstrassCurve.HasMultiplicativeReduction` が
`multiplicativeReduction : v(c₄) = 1`（乗法的表記）を**クラスの欄として持っている**。
★★加法的表記に直すだけである。 -/
theorem vAdd_c4_eq_zero {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [hm : W.HasMultiplicativeReduction R]
    (hc4 : W.c₄ ≠ 0) :
    vAdd (tateDvrVal R K) (Units.mk0 W.c₄ hc4) = 0 := by
  have h := hm.multiplicativeReduction
  have hv := valuation_eq_ofAdd_neg_vAdd (R := R) (K := K) (Units.mk0 W.c₄ hc4)
  simp only [Units.val_mk0] at hv
  rw [h] at hv
  have hc : ((Multiplicative.ofAdd (-(vAdd (tateDvrVal R K) (Units.mk0 W.c₄ hc4))) :
      Multiplicative ℤ) : WithZero (Multiplicative ℤ)) = (1 : WithZero (Multiplicative ℤ)) := hv.symm
  have h2 : (Multiplicative.ofAdd (-(vAdd (tateDvrVal R K) (Units.mk0 W.c₄ hc4))) :
      Multiplicative ℤ) = 1 := by exact_mod_cast hc
  have h3 := Multiplicative.ofAdd.injective h2
  simpa using h3

/-! ## ★★★★★★★局所高さは `j` の付値で決まる -/

/-- ★★★★★★★**局所高さ `= −v(j)`**。

原文 (GenEll p.16):
> that this definition is independent of the choice of L].

★`j = c₄³/Δ` と `v(c₄) = 0`、そして `localHeight = v(Δ)`
（`LocalHeightDelta.lean` の `localHeight_eq_vAdd_Delta`）から。

★★★**これが `Remark 3.3.1` の要である**——`j` は `K` に在る不変量なので、
潜在的乗法還元でも拡大 `L` を取らずに測れる。 -/
theorem localHeight_eq_neg_vAdd_j {R : Type} [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [hell : W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) (hc4 : W.c₄ ≠ 0) (hj : W.j ≠ 0) :
    (localHeightOf W h : ℤ) = - vAdd (tateDvrVal R K) (Units.mk0 W.j hj) := by
  have hΔ : W.Δ ≠ 0 := hell.isUnit.ne_zero
  have hjeq : (Units.mk0 W.j hj) = (Units.mk0 W.Δ hΔ)⁻¹ * (Units.mk0 W.c₄ hc4) ^ 3 := by
    apply Units.ext
    push_cast
    unfold WeierstrassCurve.j
    simp
  have hval : vAdd (tateDvrVal R K) (Units.mk0 W.j hj)
      = - vAdd (tateDvrVal R K) (Units.mk0 W.Δ hΔ)
        + 3 * vAdd (tateDvrVal R K) (Units.mk0 W.c₄ hc4) := by
    rw [hjeq]
    show Multiplicative.toAdd (tateDvrVal R K ((Units.mk0 W.Δ hΔ)⁻¹ * (Units.mk0 W.c₄ hc4) ^ 3)) = _
    rw [map_mul, map_inv, map_pow]
    simp [vAdd]
  rw [hval, vAdd_c4_eq_zero W hc4, localHeight_eq_vAdd_Delta W h]
  ring

/-- ★★★`ordAt` と `vAdd (tateDvrVal)` は同じもの（DVR では極大イデアルが唯一の高さ 1 素点）。 -/
theorem ordAt_eq_vAdd {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K] (x : Kˣ) :
    ordAt (IsDiscreteValuationRing.maximalIdeal R) x = vAdd (tateDvrVal R K) x := by
  have h1 := valuation_eq_ofAdd_neg_ordAt (IsDiscreteValuationRing.maximalIdeal R) x
  have h2 := valuation_eq_ofAdd_neg_vAdd (R := R) (K := K) x
  rw [h1] at h2
  have h3 : (Multiplicative.ofAdd (-(ordAt (IsDiscreteValuationRing.maximalIdeal R) x)) :
      Multiplicative ℤ) = Multiplicative.ofAdd (-(vAdd (tateDvrVal R K) x)) := by
    exact_mod_cast h2
  have h4 := Multiplicative.ofAdd.injective h3
  omega

/-- ★★★★★★★★**局所高さ `= −ord(j)`** —— `ordAt` の語彙で。 -/
theorem localHeight_eq_neg_ordAt_j {R : Type} [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R] [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [WeierstrassCurve.IsElliptic W] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) (hc4 : W.c₄ ≠ 0) (hj : W.j ≠ 0) :
    (localHeightOf W h : ℤ)
      = - ordAt (IsDiscreteValuationRing.maximalIdeal R) (Units.mk0 W.j hj) := by
  rw [ordAt_eq_vAdd]
  exact localHeight_eq_neg_vAdd_j W h hc4 hj

/-! ## ★★★★★★★★`Remark 3.3.1` -/

/-- ★★★★★★★★**[GenEll] Remark 3.3.1** —— 潜在的局所高さは `L` の取り方に依らない。

原文 (GenEll p.16):
> that this definition is independent of the choice of L].

★★★**原文が「one verifies immediately」で畳んだ 1 文**である。★中身は 2 つ:

1. 局所高さは `−ord(j)` である（`localHeight_eq_neg_ordAt_j`）——★`j` は `K` に在る
2. `ord_w(x)/e` は `L` に依らない（`ordAt_div_ramificationIdx_indep`）——★`x = j` で当たる

★★**「`q_E` は `K` に無いからこの補題は当たらない」という見立てが誤りだった。**
測る対象を `q_E` から `j` に取り替えると当たる。 -/
theorem potLocalHeight_indep
    {A K : Type} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A]
    [Field K] [Algebra A K] [IsFractionRing A K]
    (L L' : Type) {B B' : Type}
    [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    [Algebra A B] [Module.IsTorsionFree A B]
    [Field L] [Algebra K L] [Algebra A L] [IsScalarTower A K L]
    [Algebra B L] [IsFractionRing B L] [IsScalarTower A B L]
    [IsAdicComplete (IsLocalRing.maximalIdeal B) B]
    [CommRing B'] [IsDomain B'] [IsDiscreteValuationRing B']
    [Algebra A B'] [Module.IsTorsionFree A B']
    [Field L'] [Algebra K L'] [Algebra A L'] [IsScalarTower A K L']
    [Algebra B' L'] [IsFractionRing B' L'] [IsScalarTower A B' L']
    [IsAdicComplete (IsLocalRing.maximalIdeal B') B']
    [(IsDiscreteValuationRing.maximalIdeal B).asIdeal.LiesOver
      (IsDiscreteValuationRing.maximalIdeal A).asIdeal]
    [(IsDiscreteValuationRing.maximalIdeal B').asIdeal.LiesOver
      (IsDiscreteValuationRing.maximalIdeal A).asIdeal]
    (j : Kˣ)
    (WL : WeierstrassCurve L) [WeierstrassCurve.IsElliptic WL] [WL.IsMinimal B]
    (hL : WL.HasSplitMultiplicativeReduction B) (hc4L : WL.c₄ ≠ 0) (hjL0 : WL.j ≠ 0)
    (hjL : WL.j = algebraMap K L (j : K))
    (WL' : WeierstrassCurve L') [WeierstrassCurve.IsElliptic WL'] [WL'.IsMinimal B']
    (hL' : WL'.HasSplitMultiplicativeReduction B') (hc4L' : WL'.c₄ ≠ 0) (hjL0' : WL'.j ≠ 0)
    (hjL' : WL'.j = algebraMap K L' (j : K)) :
    ((localHeightOf WL hL : ℤ) : ℚ)
        / (((IsDiscreteValuationRing.maximalIdeal A).asIdeal.ramificationIdx
            (IsDiscreteValuationRing.maximalIdeal B).asIdeal : ℕ) : ℚ)
      = ((localHeightOf WL' hL' : ℤ) : ℚ)
        / (((IsDiscreteValuationRing.maximalIdeal A).asIdeal.ramificationIdx
            (IsDiscreteValuationRing.maximalIdeal B').asIdeal : ℕ) : ℚ) := by
  have hu : Units.map (algebraMap K L).toMonoidHom j = Units.mk0 WL.j hjL0 :=
    Units.ext (by simp [hjL])
  have hu' : Units.map (algebraMap K L').toMonoidHom j = Units.mk0 WL'.j hjL0' :=
    Units.ext (by simp [hjL'])
  have hkey := ordAt_div_ramificationIdx_indep (A := A) (K := K) L L'
    (IsDiscreteValuationRing.maximalIdeal A) (IsDiscreteValuationRing.maximalIdeal B)
    (IsDiscreteValuationRing.maximalIdeal B') j
  rw [hu, hu'] at hkey
  rw [localHeight_eq_neg_ordAt_j WL hL hc4L hjL0, localHeight_eq_neg_ordAt_j WL' hL' hc4L' hjL0']
  push_cast
  rw [neg_div, neg_div, hkey]

/-! ### ★★★★★★★★項目全体の `.src`

★`.src` は「その原典項目を**完全に**実装した」という主張である
（`tools/genell-progress.mjs` の規則）。 -/

/-- ★★★★★★★★**[GenEll] Remark 3.3.1** —— 実装された。

原文 (GenEll p.16):
> that this definition is independent of the choice of L].

## ★主張

| 原文 | 宣言 |
|---|---|
| 潜在的乗法還元での局所高さ `= v_L(q_{E⊗L})/e(L/K) ∈ ℚ` | ★`localHeightOf WL hL / e`（本ファイル） |
| ★**`L` の取り方に依らない** | ★★★`potLocalHeight_indep`（本ファイル） |
| その実質: `v_L(q) = e(L/K)·v_K(j) の符号違い` | ★`localHeight_eq_neg_ordAt_j` ＋ `ordAt_div_ramificationIdx_indep` |

## ★★★★逸脱の記録（CLAUDE.md の「逸脱」）

### 1. `E ⊗_K L` の極小モデルを仮定として受ける

原文は `E ⊗_K L` の局所高さを `e(L/K)` で割ると書く。★本実装は
「`L` 上の極小 Weierstrass モデル `W_L` で、その `j` が `K` の `j` の像であるもの」
として受ける（`hjL : WL.j = algebraMap K L j`）。
★★`E ⊗_K L` の極小モデルが実際にそうなることは `j` がモデルに依らない不変量である
ことから従うが、**その段は含んでいない**。

### 2. 局所体を離散付値環で受ける

原文の `K` は `ℚ_p` の有限拡大、`L` はその有限拡大である。
★本実装は `A, B` を離散付値環として受ける。`ordAt` は Dedekind 環の高さ 1 素点で
定義されているが、DVR ではそれが極大イデアルただ 1 つなので一致する
（`ordAt_eq_vAdd`）。

### 3. ★★★測る対象を `q_E` から `j` に取り替えた

★**これは逸脱ではなく発見である。** `LocalHeightRamified.lean` の
`ordAt_div_ramificationIdx_indep` は `x : Kˣ` を要求するので、
「`q_E` は `K` に無いから当たらない」と見えていた。
★★★**`j` は `K` に在る**（曲線の不変量）ので、`localHeight = −v(j)` を経由すると当たる。 -/
def remark_3_3_1.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 16, item := "Remark 3.3.1",
    sectionId := "genell-rem-3-3-1" }

def remark_3_3_1.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "potLocalHeight_indep(L の取り方に依らないこと——主張そのもの)"
      (.inProject "ABC3" "ABC3.Found.GenEll.potLocalHeight_indep") 16,
    .citation "[ABC3]" "localHeight_eq_neg_ordAt_j(局所高さ = −ord(j))"
      (.inProject "ABC3" "ABC3.Found.GenEll.localHeight_eq_neg_ordAt_j") 16,
    .citation "[ABC3]" "ordAt_div_ramificationIdx_indep(ord_w(x)/e は L に依らない)"
      (.inProject "ABC3" "ABC3.Found.GenEll.ordAt_div_ramificationIdx_indep") 16,
    .citation "[ABC3]" "localHeight_eq_vAdd_Delta(局所高さ = 極小判別式の付値)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.localHeight_eq_vAdd_Delta") 15,
    .citation "[mathlib]" "WeierstrassCurve.HasMultiplicativeReduction(v(c₄) = 0 を欄として持つ)"
      (.inMathlib "WeierstrassCurve.HasMultiplicativeReduction") 16,
    .implicitStep
      ("★逸脱 1: E ⊗_K L の極小モデルを仮定として受ける(hjL : WL.j = algebraMap K L j)。" ++
       "実際にそうなることは j がモデルに依らない不変量であることから従うが、" ++
       "その段は含んでいない") 16,
    .implicitStep
      ("★逸脱 2: 局所体を離散付値環で受ける。ordAt は Dedekind 環の高さ 1 素点で" ++
       "定義されているが、DVR では極大イデアルただ 1 つなので一致する(ordAt_eq_vAdd)") 16,
    .implicitStep
      ("★★★測る対象を q_E から j に取り替えた。ordAt_div_ramificationIdx_indep は" ++
       "x : Kˣ を要求するので「q_E は K に無いから当たらない」と見えていたが、" ++
       "j は K に在る(曲線の不変量)ので localHeight = −v(j) を経由すると当たる") 16 ]

/-! ### ★出典の紐付け(`.src`) -/

def vAdd_c4_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 16,
    item := "Remark 3.3.1(乗法的還元なら v(c₄) = 0)",
    sectionId := "genell-rem-3-3-1" }

def localHeight_eq_neg_vAdd_j.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 16,
    item := "Remark 3.3.1(局所高さ = −v(j))",
    sectionId := "genell-rem-3-3-1" }

def ordAt_eq_vAdd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 16,
    item := "Remark 3.3.1(ordAt と vAdd(tateDvrVal) は同じもの)",
    sectionId := "genell-rem-3-3-1" }

def localHeight_eq_neg_ordAt_j.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 16,
    item := "Remark 3.3.1(局所高さ = −ord(j)——ordAt の語彙で)",
    sectionId := "genell-rem-3-3-1" }

def potLocalHeight_indep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 16,
    item := "Remark 3.3.1(潜在的局所高さは L の取り方に依らない)",
    sectionId := "genell-rem-3-3-1" }

end ABC3.Found.GenEll
