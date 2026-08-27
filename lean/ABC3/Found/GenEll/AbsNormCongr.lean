/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Conductor
import ABC3.Found.GenEll.CartierPullback
import ABC3.Found.GenEll.LogDiffCongr
import ABC3.Found.GenEll.LogCondSigma

/-!
# [GenEll] Definition 1.5, (iv) —— **絶対ノルムは環同型で保たれる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> determines a well-defined log-conductor function log-condD on UX(Q) ⊆X(Q).

## ★★何のために要るか

`log-cond_D` を `U_X(ℚ̄)` の関数として well-defined にするには、
`log-diff` と同じ 2 つの側が要る（`LogDiffCongr.lean` の構図）:

| 側 | 主張 | 状態 |
|---|---|---|
| **幾何** | 同じ点の最小定義体は互いに対応する | ★`minField_baseChange`（取得済み） |
| **代数** | 対応する体は同じ `log-cond` を与える | ★**その核が本ファイル** |

★`log-cond_D(x) = deg_F(f^D_x)` で `f^D_x = (D_x)_red` は**イデアル**なので、
代数の側の核は「**イデアルの絶対ノルムが環同型で保たれる**」ことである。

## ★★★mathlib に無かった（2026-08-27 実測）

`Ideal.absNorm` を環同型で移す補題は見つからなかった
（`exact?` が空振り）。★ただし材料はすべてある:

* `Ideal.absNorm_apply` —— `absNorm I = Submodule.cardQuot I`
* `Ideal.quotientEquiv` —— 環同型は剰余環の同型を誘導する

★★**商の濃度が等しいことに落ちる**ので、3 行で出る。
`absNorm` を「ℤ-加群としての指数」ではなく**剰余環の濃度**として読むのが要点である。
-/

namespace ABC3.Found.GenEll

open Ideal

/-- ★★**イデアルの絶対ノルムは環同型で保たれる**。

原文 (GenEll p.8):
> determines a well-defined log-conductor function log-condD on UX(Q) ⊆X(Q).

★★★中身は「環同型が剰余環の同型を誘導する」だけである
（`Ideal.quotientEquiv`）。`absNorm` を剰余環の濃度として読む
（`Ideal.absNorm_apply`）と、そこに落ちる。

★これが `log-cond` の同型不変性（`Definition 1.5, (iv)` の代数の側）の核である。 -/
theorem absNorm_map_ringEquiv {R S : Type} [CommRing R] [IsDedekindDomain R] [Module.Free ℤ R]
    [CommRing S] [IsDedekindDomain S] [Module.Free ℤ S]
    (e : R ≃+* S) (I : Ideal R) :
    Ideal.absNorm (I.map (e : R →+* S)) = Ideal.absNorm I := by
  rw [Ideal.absNorm_apply, Ideal.absNorm_apply, Submodule.cardQuot, Submodule.cardQuot]
  have h : (R ⧸ I) ≃+* (S ⧸ I.map (e : R →+* S)) :=
    Ideal.quotientEquiv I (I.map (e : R →+* S)) e rfl
  exact congrArg _ (Cardinal.mk_congr h.toEquiv).symm

/-- ★根基も環同型で移る（`(−)_red` が同型で保たれること）。 -/
theorem radical_map_ringEquiv {R S : Type} [CommRing R] [CommRing S]
    (e : R ≃+* S) (I : Ideal R) :
    (I.map (e : R →+* S)).radical = I.radical.map (e : R →+* S) :=
  (Ideal.map_radical_of_surjective e.surjective
    (by rw [RingHom.ker_coe_equiv]; exact bot_le)).symm

/-- ★★★**導手の絶対ノルムは環同型で保たれる** —— `(−)_red` を挟んでも同じ。

★`log-cond_D(x) = deg_F((D_x)_red)` の中身がこれである。 -/
theorem absNorm_radical_map_ringEquiv {R S : Type} [CommRing R] [IsDedekindDomain R]
    [Module.Free ℤ R] [CommRing S] [IsDedekindDomain S] [Module.Free ℤ S]
    (e : R ≃+* S) (I : Ideal R) :
    Ideal.absNorm ((I.map (e : R →+* S)).radical) = Ideal.absNorm I.radical := by
  rw [radical_map_ringEquiv e I, absNorm_map_ringEquiv e I.radical]

/-! ## ★★★★★`log-cond` の同型不変性 -/

theorem map_ringEquiv_ne_zero {R S : Type} [CommRing R] [CommRing S]
    (e : R ≃+* S) {I : Ideal R} (h : I ≠ 0) : I.map (e : R →+* S) ≠ 0 := by
  intro hm
  apply h
  rw [Ideal.zero_eq_bot] at hm ⊢
  have hc : Ideal.comap (e : R →+* S) (Ideal.map (e : R →+* S) I)
      = Ideal.comap (e : R →+* S) (⊥ : Ideal S) := congrArg _ hm
  rw [Ideal.comap_map_of_bijective (e : R →+* S) e.bijective,
    Ideal.comap_bot_of_injective (e : R →+* S) e.injective] at hc
  exact hc

open AlgebraicGeometry CategoryTheory NumberField in
/-- ★★★★★**イデアルが対応すれば `log-cond` は一致する** ——
`Definition 1.5, (iv)` の**代数の側**。

原文 (GenEll p.8):
> determines a well-defined log-conductor function log-condD on UX(Q) ⊆X(Q).

★★**欠けている 1 点を仮定 `h` として分離してある** ——
「引き戻したイデアルが同型で対応する」ことである。
★それは「アフィンでイデアル層の `comap` が `Ideal.map` に対応する」から出るはずだが、
**mathlib に無い**（`LogDiffPoint.lean` の `comap_baseChange` に実測を記録）。

★★★**それ以外はすべて本定理で取ってある**:
`log-cond = log(absNorm (D_x)_red) / [F:ℚ]` の
分子は `absNorm_radical_map_ringEquiv`、分母は `finrank_congr` で一致する。

★`h` が埋まれば、`(iii)` の `logDiff_minField_baseChange` と同じ形で
`log-cond_D` が `U_X(ℚ̄)` の関数として well-defined になる。 -/
theorem logCond_congr_of_ideal_map (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] {X : Scheme.{0}} (D : X.IdealSheafData)
    (xF : specRingOfIntegers F ⟶ X) (xK : specRingOfIntegers K ⟶ X)
    (eF : F ≃+* K) (eO : (𝓞 F) ≃+* (𝓞 K))
    (hne : pullbackIdeal F D xF ≠ 0)
    (h : pullbackIdeal K D xK = (pullbackIdeal F D xF).map (eO : (𝓞 F) →+* (𝓞 K))) :
    logCond K D xK = logCond F D xF := by
  have hrad := radical_ne_zero hne
  have hradK : (pullbackIdeal K D xK).radical ≠ 0 :=
    radical_ne_zero (h ▸ map_ringEquiv_ne_zero eO hne)
  rw [logCond, logCond, conductorADiv, conductorADiv, degNormalized, degNormalized,
    deg_idealADiv K _ hradK, deg_idealADiv F _ hrad, h,
    absNorm_radical_map_ringEquiv eO (pullbackIdeal F D xF), finrank_congr eF]

/-! ### ★出典の紐付け(`.src`) -/

def logCond_congr_of_ideal_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iv)(イデアルが対応すれば log-cond は一致する——代数の側)",
    sectionId := "genell-def-1-5" }

def logCond_congr_of_ideal_map.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "absNorm_radical_map_ringEquiv(分子——導手の絶対ノルムは同型で保たれる)"
      (.inProject "ABC3" "ABC3.Found.GenEll.absNorm_radical_map_ringEquiv") 8,
    .citation "[ABC3]" "finrank_congr(分母——同型な数体は同じ次数)"
      (.inProject "ABC3" "ABC3.Found.GenEll.finrank_congr") 8,
    .implicitStep
      ("★★欠けている 1 点を仮定 h として分離してある ——「引き戻したイデアルが同型で対応する」。" ++
       "それは「アフィンでイデアル層の comap が Ideal.map に対応する」から出るはずだが" ++
       "mathlib に無い(2026-08-27 実測、LogDiffPoint.lean の comap_baseChange に記録)。" ++
       "★h が埋まれば (iii) の logDiff_minField_baseChange と同じ形で " ++
       "log-cond_D が U_X(ℚ̄) の関数として well-defined になる") 8 ]

def absNorm_map_ringEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iv)(イデアルの絶対ノルムは環同型で保たれる)",
    sectionId := "genell-def-1-5" }

def absNorm_map_ringEquiv.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Ideal.absNorm_apply(absNorm = 剰余環の濃度)"
      (.inMathlib "Ideal.absNorm_apply") 8,
    .citation "[mathlib]" "Ideal.quotientEquiv(環同型は剰余環の同型を誘導する)"
      (.inMathlib "Ideal.quotientEquiv") 8,
    .implicitStep
      ("★★mathlib に absNorm を環同型で移す補題は無い(2026-08-27 実測、exact? が空振り)。" ++
       "★absNorm を「ℤ-加群としての指数」ではなく剰余環の濃度として読むと 3 行で出る") 8 ]

def absNorm_radical_map_ringEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iv)(導手の絶対ノルムは環同型で保たれる)",
    sectionId := "genell-def-1-5" }

end ABC3.Found.GenEll
