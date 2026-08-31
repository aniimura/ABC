/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.MetricRatio
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★`Spec ℤ` から来た算術直線束の高さは**定数**である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z)

## ★★★★★★★★★★これは何か —— 原文が 1 行で使う段

原文は `Proposition 1.4, (ii)` の証明でこう書く:

> Now observe that if M̄ is an arithmetic line bundle that arises [by pull-back to X]
> from an arithmetic line bundle on Spec(ℤ), then `ht_{L̄⊗M̄} ≈_{X(ℚ̄)} ht_L̄` [cf. assertion (i)].

★★本ファイルはそれを**定数まで込めて等号で**示す:

    `ht_{L̄ ⊗ g^*N̄}(x) = ht_L̄(x) + deg_ℚ(N̄)`   （`htMetricAlg_mul_pullback`）

★★★原文の `≈` より強い（差が**厳密に定数**である）。

## ★★★★★★★★機構 —— 底変換不変性がそのまま効く

    `ht_{g^*N̄}(x_F) = deg_F((x_F ≫ g)^* N̄)`   （引き戻しの関手性、`§9-747`）

★`Spec 𝓞_ℚ` は**終対象**なので `x_F ≫ g` は
`Spec.map (algebraMap 𝓞_ℚ 𝓞_F)` **そのもの**である（`eq_algebraMap_of_terminal`）。

★★したがって `§9-799` の底変換不変性

    `deg_K(L̄|_{Spec 𝓞_K}) = deg_F(L̄)`

を `F = ℚ`, `K = F` で当てるだけで `= deg_ℚ(N̄)` が出る
——**`F` にも `x_F` にも依らない定数**である。

## ★★★★★終対象性の配管

`𝓞_ℚ` は `ℤ` そのものではないので、`Rat.ringOfIntegersEquiv : 𝓞 ℚ ≃+* ℤ` で
`Spec 𝓞_ℚ ≅ Spec ℤ` を作り、mathlib の `specZIsTerminal` を運ぶ
（`specRatIsTerminal`）。★これだけで「任意の射が `algebraMap` である」が出る。

## ★★★★★★★これが原文の (ii) の証明の骨格である

`exists_htMetricU_ge_of_twist`——捻った束 `L̄ ⊗ g^*N̄` の側で
`§9-805` の下界を取れば、定数分ずらして元の高さの下界になる。

★残るのは原文の「`L_ℚ` のある正冪が大域切断で生成される ⟹ 捻った束が
`|t| ≤ 1` の切断を持つ」——**分母を払う段**である。本ファイルには含まない。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite NumberField
open ABC3.Found.GenEll

/-! ## ★★★(1) 引き戻しの関手性を高さで読む -/

/-- ★★**`ht_{g^*N̄}(x_F) = ht_N̄(x_F ≫ g)`**。 -/
theorem htMetricU_comp {X Y : Scheme.{0}} (F : Type) [Field F] [NumberField F]
    (g : X ⟶ Y) (N : AInv Y) (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) :
    htMetricU F (AInv.pullback g N) xF = htMetricU F N (xF ≫ g) := by
  show degAPicM F (APicMPullback xF (APicM.mk (AInv.pullback g N))) = _
  have h : APicM.mk (AInv.pullback g N) = APicMPullback g (APicM.mk N) := rfl
  rw [h, APicMPullback_comp]
  rfl

/-! ## ★★★★★(2) `Spec 𝓞_ℚ` は終対象 -/

/-- ★★**`Spec (𝓞 ℚ)` は終対象である**。

★`Rat.ringOfIntegersEquiv : 𝓞 ℚ ≃+* ℤ` で mathlib の `specZIsTerminal` を運ぶ。 -/
noncomputable def specRatIsTerminal : Limits.IsTerminal (Spec (CommRingCat.of (𝓞 ℚ))) :=
  Limits.IsTerminal.ofIso specZIsTerminal
    (Scheme.Spec.mapIso (Iso.op (Rat.ringOfIntegersEquiv.toCommRingCatIso)))

/-- ★★★**`Spec 𝓞_F ⟶ Spec 𝓞_ℚ` の射は `algebraMap` そのものである**。 -/
theorem eq_algebraMap_of_terminal (F : Type) [Field F] [NumberField F]
    (u : Spec (CommRingCat.of (𝓞 F)) ⟶ Spec (CommRingCat.of (𝓞 ℚ))) :
    u = Spec.map (CommRingCat.ofHom (algebraMap (𝓞 ℚ) (𝓞 F))) :=
  specRatIsTerminal.hom_ext u _

/-! ## ★★★★★★★★★(3) 高さは定数 -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★**`Spec ℤ` から来た算術直線束の高さは定数 `deg_ℚ(N̄)` である**。

原文 (GenEll p.6):
> that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z)

★機構は「終対象性で射が決まる」＋「`§9-799` の底変換不変性」だけである。 -/
theorem htMetricU_pullback_const (F : Type) [Field F] [NumberField F] {X : Scheme.{0}}
    (g : X ⟶ Spec (CommRingCat.of (𝓞 ℚ))) (N : AInv (Spec (CommRingCat.of (𝓞 ℚ))))
    (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X) :
    htMetricU F (AInv.pullback g N) xF = degAPicM ℚ (APicM.mk N) := by
  rw [htMetricU_comp F g N xF, eq_algebraMap_of_terminal F (xF ≫ g)]
  exact degAPicM_baseChange ℚ F (APicM.mk N)

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★`X(ℚ̄)` の上で定数である。 -/
theorem htMetricAlg_pullback_const {X : Scheme.{0}}
    (g : X ⟶ Spec (CommRingCat.of (𝓞 ℚ))) (N : AInv (Spec (CommRingCat.of (𝓞 ℚ))))
    (x : AlgPointAnyClass X) :
    htMetricAlg (AInv.pullback g N) x = degAPicM ℚ (APicM.mk N) := by
  induction x using AlgPointAnyClass.ind with
  | _ p => exact @htMetricU_pullback_const p.fld p.instField p.instNF X g N p.map

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★**`ht_{L̄ ⊗ g^*N̄} = ht_L̄ + deg_ℚ(N̄)`**。

原文 (GenEll p.6):
> that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z)

★★原文は `≈`（BD-同値）と書くが、**差はきっかり定数**である。 -/
theorem htMetricAlg_mul_pullback {X : Scheme.{0}} (L : AInv X)
    (g : X ⟶ Spec (CommRingCat.of (𝓞 ℚ))) (N : AInv (Spec (CommRingCat.of (𝓞 ℚ))))
    (x : AlgPointAnyClass X) :
    htMetricAlg (L.mul (AInv.pullback g N)) x
      = htMetricAlg L x + degAPicM ℚ (APicM.mk N) := by
  rw [htMetricAlg_mul, htMetricAlg_pullback_const]

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★**BD-類は `Spec ℤ` からの捻りで変わらない**——原文の `≈` そのもの。 -/
theorem heightBDClass_mul_pullback {X : Scheme.{0}} (L : AInv X)
    (g : X ⟶ Spec (CommRingCat.of (𝓞 ℚ))) (N : AInv (Spec (CommRingCat.of (𝓞 ℚ))))
    (S : Set (AlgPointAnyClass X)) :
    heightBDClass (L.mul (AInv.pullback g N)) S = heightBDClass L S := by
  refine (BDClass.mk_eq_mk _ _).2 ⟨|degAPicM ℚ (APicM.mk N)|, fun x => ?_⟩
  show |heightOn (L.mul (AInv.pullback g N)) S x - heightOn L S x| ≤ _
  show |htMetricAlg (L.mul (AInv.pullback g N)) (x : AlgPointAnyClass X)
    - htMetricAlg L (x : AlgPointAnyClass X)| ≤ _
  rw [htMetricAlg_mul_pullback L g N]
  simp

/-! ## ★★★★★★★★★★(4) 原文の (ii) の証明の骨格 -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★**捻った束の側で下界を取れば元の高さの下界が出る**。

原文 (GenEll p.6):
> that arises [by pull-back to X] from an arithmetic line bundle on Spec(Z)

★これが原文の `Proposition 1.4, (ii)` の証明の骨格である
——`§9-805` の下界を `L̄ ⊗ g^*N̄` に当て、定数 `deg_ℚ(N̄)` の分だけずらす。

★★残るのは「`L_ℚ` のある正冪が大域切断で生成される ⟹ 捻った束が切断を持つ」
（**分母を払う段**）であり、本ファイルには含まない。 -/
theorem exists_htMetricU_ge_of_twist {X : Scheme.{0}} [CompactSpace X]
    (hval : ValuativeCriterion (specZIsTerminal.from X))
    [Nonempty (Spec (CommRingCat.of ℂ) ⟶ X)]
    (L : AInv X) (g : X ⟶ Spec (CommRingCat.of (𝓞 ℚ)))
    (N : AInv (Spec (CommRingCat.of (𝓞 ℚ))))
    (hc : (L.mul (AInv.pullback g N)).carrier.IsContinuous)
    (s : ((L.mul (AInv.pullback g N)).carrier.sheaf.obj (op ⊤) : Type)) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ (F : Type) [Field F] [NumberField F]
      (xF : Spec (CommRingCat.of (𝓞 F)) ⟶ X),
      pullSecT xF (L.mul (AInv.pullback g N)).carrier.sheaf s ≠ 0 →
      -C ≤ htMetricU F L xF := by
  obtain ⟨C, hC, h⟩ := exists_htMetricU_ge hval (L.mul (AInv.pullback g N)) hc s
  refine ⟨C + |degAPicM ℚ (APicM.mk N)|, by positivity, fun F _ _ xF hne => ?_⟩
  have h1 := h F xF hne
  have h2 : htMetricU F (L.mul (AInv.pullback g N)) xF
      = htMetricU F L xF + degAPicM ℚ (APicM.mk N) := by
    rw [htMetricU_mul, htMetricU_pullback_const F g N xF]
  have h3 : degAPicM ℚ (APicM.mk N) ≤ |degAPicM ℚ (APicM.mk N)| := le_abs_self _
  rw [h2] at h1
  linarith

/-! ## ★出典の紐付け(`.src`) -/

def specRatIsTerminal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(Spec 𝓞_ℚ は終対象——Spec ℤ から来る束の同定に使う)",
    sectionId := "genell-prop-1-4" }

def htMetricU_pullback_const.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(Spec ℤ から来た算術直線束の高さは定数 deg_ℚ(N̄))",
    sectionId := "genell-prop-1-4" }

def htMetricAlg_mul_pullback.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(ht_{L̄⊗M̄} = ht_L̄ + deg_ℚ(N̄)——原文の ≈ を等号で)",
    sectionId := "genell-prop-1-4" }

def exists_htMetricU_ge_of_twist.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (ii)(捻った束の側で下界を取る段——分母を払う段は含まない)",
    sectionId := "genell-prop-1-4" }

def htMetricAlg_mul_pullback.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "APicMPullback_comp(引き戻しの関手性、§9-747)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.APicMPullback_comp") 3,
    .citation "[ABC3]" "degAPicM_baseChange(deg_K(L|Spec 𝓞_K) = deg_F(L)、§9-799)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.degAPicM_baseChange") 4,
    .citation "[mathlib]" "specZIsTerminal / Rat.ringOfIntegersEquiv"
      (.inMathlib "AlgebraicGeometry.specZIsTerminal") 6,
    .implicitStep
      ("★残るのは原文の「L_ℚ のある正冪が大域切断で生成される ⟹ 捻った束が " ++
       "|t| ≤ 1 の切断を持つ」——**分母を払う段**である。" ++
       "★★因子表示では Found/GenEll/VerticalTwist.lean が同種の段を持つ") 6 ]

end ABC3.Found.Arakelov
