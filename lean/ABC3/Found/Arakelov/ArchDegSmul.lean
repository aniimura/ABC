/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.ArchDeg
import ABC3.Found.Arakelov.AMetricHom
import ABC3.Found.Arakelov.ArcEvalCont
import ABC3.Meta.Claim

/-!
# 段 D の**アルキメデス側** —— `archDeg (c · s) = archDeg s − log |N(c)| / [F:Q]`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★★★これは何か

`deg_F(L̄)` の古典的な定義

    `deg_F(L̄) = ( log #(Gamma(L)/O_F·s) − Sum_sigma log |s|_sigma ) / [F:Q]`

は**切断 `s ≠ 0` の取り方に依らない**——それが台帳の**段 D** である。
`s` を `c·s`（`c ∈ O_F`、`c ≠ 0`）に取り替えると

| 側 | ずれ |
|---|---|
| 有限側（`degFin_smul_norm`、§9-774） | `+ log |N(c)|` |
| ★**アルキメデス側（本ファイル）** | `− log |N(c)| / [F:Q]` |

★★本ファイルは**アルキメデス側のずれを計算し切る**。
★★★両者が相殺するのが**積公式**であり、正規化（`/[F:Q]`）を揃えれば
`degFin_smul_norm` の右辺とちょうど打ち消す。

## ★★★★★機構（3 段）

    (1) evalOn_specMap        : `Spec.map f` での点の値は `f` そのもの（Gamma-Spec 随伴の自然性）
    (2) sum_log_norm_emb      : `Sum_sigma log |sigma(u)| = log |N(u)|`（norm_eq_prod_embeddings）
    (3) archDeg_smul          : `AMetric.norm_smul` に (1)(2) を差し込む

★(1) は「複素点で関数を評価する」ことの**環論的な正体**である——
`embSpecPoint F sigma` は `Spec.map (sigma ∘ algebraMap)` なので、
そこでの評価は `sigma ∘ algebraMap` の適用そのものになる。
★★これは在庫の `evalOn_appLE`（`rfl`）と mathlib の `Scheme.ΓSpecIso_naturality` の合成で出る。

★★★(2) は mathlib の `Algebra.norm_eq_prod_embeddings`（`N(u) = prod_sigma sigma(u)`）に
`RingHom.equivRatAlgHom`（`F →+* C` と `F →ₐ[Q] C` の対応）を噛ませたものである。
`Real.log_prod` で積が和になる。

## ★残っている段（明示）

★★段 D の**残り**は、有限側とアルキメデス側を**同じ切断 `s` の上で足す**こと
——`degFin` は `AInv.toInvSheaf` / `moduleSpecΓFunctor` の世界、
`archDeg` は `AMetric` / `L.sheaf.obj (op ⊤)` の世界にあり、その同一視が要る。
★★★段 E（加法性）は別。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite NumberField

/-! ## ★★★★★★(1) `Spec.map f` での点の値 -/

/-- ★★★★★★**`Spec.map f` が定める点での関数の値は `f` の適用そのものである**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★在庫の `evalOn_appLE`（`evalOn` は `appLE` で書ける、`rfl`）と
mathlib の `Scheme.ΓSpecIso_naturality` の合成である。 -/
theorem evalOn_specMap (R : CommRingCat.{0}) (f : R ⟶ CommRingCat.of ℂ)
    (h : (Spec.map f) ⁻¹ᵁ (⊤ : (Spec R).Opens) = ⊤)
    (c : (Γ(Spec R, (⊤ : (Spec R).Opens)) : Type)) :
    evalOn (Spec.map f) ⊤ h c = f.hom ((Scheme.ΓSpecIso R).hom.hom c) := by
  rw [evalOn_appLE]
  have hEq : ((Spec.map f).appLE ⊤ ⊤ (le_of_eq h.symm)) = (Spec.map f).appTop :=
    (Scheme.Hom.app_eq_appLE (Spec.map f)).symm
  rw [hEq]
  exact congrArg (fun (m : _ ⟶ _) => CommRingCat.Hom.hom m c) (Scheme.ΓSpecIso_naturality f)

/-- ★★★★★**埋め込みが定める複素点での値は `σ` の適用である**。 -/
theorem evalOn_embSpecPoint (F : Type) [Field F] [NumberField F] (σ : F →+* ℂ)
    (h : (embSpecPoint F σ) ⁻¹ᵁ (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens) = ⊤)
    (c : (Γ(Spec (CommRingCat.of (𝓞 F)),
      (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)) :
    evalOn (embSpecPoint F σ) ⊤ h c
      = σ (algebraMap (𝓞 F) F ((Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom.hom c)) :=
  evalOn_specMap (CommRingCat.of (𝓞 F)) _ h c

/-! ## ★★★★★★★(2) 埋め込みの絶対値の積はノルム -/

/-- ★整数の絶対値の実数への像。 -/
theorem natAbs_cast_real (n : ℤ) : ((n.natAbs : ℕ) : ℝ) = |(n : ℝ)| := by
  rw [Nat.cast_natAbs]; simp

/-- ★★★★★★★**`Σ_σ log |σ(u)| = log |N_{F/Q}(u)|`**（体の側）。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★mathlib の `Algebra.norm_eq_prod_embeddings`（`N(u) = Π_σ σ(u)`、`σ : F →ₐ[ℚ] ℂ`）に
`RingHom.equivRatAlgHom` を噛ませて `F →+* ℂ` の側に移し、`Real.log_prod` で和にした。 -/
theorem sum_log_norm_emb_field (F : Type) [Field F] [NumberField F] (u : F) (hu : u ≠ 0) :
    ∑ σ : F →+* ℂ, Real.log ‖σ u‖ = Real.log |(Algebra.norm ℚ u : ℝ)| := by
  classical
  have hprod : ∏ σ : F →+* ℂ, ‖σ u‖ = |(Algebra.norm ℚ u : ℝ)| := by
    have h := congrArg (‖·‖) (Algebra.norm_eq_prod_embeddings (K := ℚ) (L := F) (E := ℂ) u)
    simp only [norm_prod] at h
    rw [Fintype.prod_equiv RingHom.equivRatAlgHom (fun f : F →+* ℂ => ‖f u‖)
      (fun φ : F →ₐ[ℚ] ℂ => ‖φ u‖) (fun _ => by simp [RingHom.equivRatAlgHom_apply])]
    rw [← h]
    simp [Complex.norm_ratCast]
  rw [← hprod, Real.log_prod]
  intro σ _
  simpa using (map_ne_zero σ).mpr hu

/-- ★★★★★★★★**`Σ_σ log |σ(u)| = log |N(u)|`**（整数環の側）。

★有限側の `degFin_smul_norm`（§9-774）が `log |N(u)|` の形で書かれているので、
★★**両者はこの補題を蝶番にして直接つながる**。 -/
theorem sum_log_norm_emb (F : Type) [Field F] [NumberField F] (u : 𝓞 F) (hu : u ≠ 0) :
    ∑ σ : F →+* ℂ, Real.log ‖σ (algebraMap (𝓞 F) F u)‖
      = Real.log ((((Algebra.norm ℤ) u).natAbs : ℕ) : ℝ) := by
  have hu' : (algebraMap (𝓞 F) F u) ≠ 0 := by
    have h := RingOfIntegers.coe_eq_zero_iff (K := F) (x := u)
    simpa [h] using hu
  rw [sum_log_norm_emb_field F _ hu', natAbs_cast_real]
  congr 2
  rw [← Algebra.coe_norm_int]
  push_cast
  rfl

/-! ## ★★★★★★★★★★(3) 段 D のアルキメデス側 -/

/-- ★★★★★★★★★★**段 D のアルキメデス側**——切断を `c` 倍すると
アルキメデス次数は `log |N(c)| / [F:ℚ]` だけ下がる。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★★有限側の `degFin_smul_norm`（§9-774）は

    `degFin (c·s) = degFin s + log |N(c)|`

であった。★★★★正規化（`/[F:ℚ]`）を揃えれば**ちょうど打ち消す**——
それが `deg_F` が切断の取り方に依らないことであり、**積公式**である。 -/
theorem archDeg_smul (F : Type) [Field F] [NumberField F]
    (L : AMetric (Spec (CommRingCat.of (𝓞 F))))
    (c : (Γ(Spec (CommRingCat.of (𝓞 F)),
      (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type))
    (s : (L.sheaf.obj (op ⊤) : Type))
    (hc : (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom.hom c ≠ 0)
    (hs : ∀ σ : F →+* ℂ, L.norm s (embSpecPoint F σ) ≠ 0) :
    archDeg F L (c • s)
      = archDeg F L s
        - Real.log ((((Algebra.norm ℤ)
              ((Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom.hom c)).natAbs : ℕ) : ℝ)
            / (Module.finrank ℚ F : ℝ) := by
  set u := (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom.hom c with hu
  have htop : ∀ σ : F →+* ℂ,
      (embSpecPoint F σ) ⁻¹ᵁ (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens) = ⊤ := by
    intro σ; simp
  have hterm : ∀ σ : F →+* ℂ,
      Real.log (L.norm (c • s) (embSpecPoint F σ))
        = Real.log ‖σ (algebraMap (𝓞 F) F u)‖ + Real.log (L.norm s (embSpecPoint F σ)) := by
    intro σ
    rw [AMetric.norm_smul L c s (embSpecPoint F σ) (htop σ),
      evalOn_embSpecPoint F σ (htop σ) c]
    refine Real.log_mul ?_ (hs σ)
    simpa using (map_ne_zero σ).mpr
      (by simpa [RingOfIntegers.coe_eq_zero_iff] using hc :
        (algebraMap (𝓞 F) F u) ≠ 0)
  show -(∑ σ : F →+* ℂ, Real.log (L.norm (c • s) (embSpecPoint F σ))) / _ = _
  simp only [hterm, Finset.sum_add_distrib]
  rw [sum_log_norm_emb F u (by simpa [hu] using hc)]
  show _ = -(∑ σ : F →+* ℂ, Real.log (L.norm s (embSpecPoint F σ))) / _ - _
  ring

/-! ### ★出典の紐付け(`.src`) -/

def evalOn_specMap.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(Spec.map f が定める点での関数の値は f そのもの)",
    sectionId := "genell-def-1-1-ii" }

def sum_log_norm_emb.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(Σ_σ log |σ(u)| = log |N(u)|——積公式のアルキメデス側)",
    sectionId := "genell-def-1-1-ii" }

def archDeg_smul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(段 D のアルキメデス側——archDeg(c·s) = archDeg(s) − log |N(c)|/[F:Q])",
    sectionId := "genell-def-1-1-ii" }

def archDeg_smul.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "AMetric.norm_smul(|c·s|(p) = |c(p)|·|s|(p))"
      (.inProject "ABC3" "ABC3.Found.Arakelov.AMetric.norm_smul") 3,
    .citation "[ABC3]" "evalOn_appLE(evalOn は appLE で書ける、rfl)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.evalOn_appLE") 3,
    .citation "[mathlib]" "Scheme.ΓSpecIso_naturality(Γ-Spec 随伴の自然性)"
      (.inMathlib "AlgebraicGeometry.Scheme.ΓSpecIso_naturality") 4,
    .citation "[mathlib]" "Algebra.norm_eq_prod_embeddings(N(u) = Π_σ σ(u))"
      (.inMathlib "Algebra.norm_eq_prod_embeddings") 4,
    .citation "[mathlib]" "Algebra.coe_norm_int(整数環のノルムと体のノルムの一致)"
      (.inMathlib "Algebra.coe_norm_int") 4,
    .citation "[ABC3]" "degFin_smul_norm(段 D の有限側、§9-774)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.degFin_smul_norm") 4,
    .implicitStep
      ("★★段 D の残りは、有限側とアルキメデス側を**同じ切断 s の上で足す**ことである" ++
       "——degFin は AInv.toInvSheaf / moduleSpecΓFunctor の世界、" ++
       "archDeg は AMetric / L.sheaf.obj (op ⊤) の世界にあり、その同一視が要る") 4 ]

end ABC3.Found.Arakelov
