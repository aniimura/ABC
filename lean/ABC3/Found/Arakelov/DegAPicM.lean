/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.DegArithIsometry
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★**`deg_F : APic(Spec 𝓞_F) → ℝ` を無条件に作る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★★★これは何か

原文は `deg_F : APic(Spec(O_F)) → R` を
「`ADiv(F)/APrc(F) ≅ APic(Spec(O_F))`（[Szp] Prop 1.1）」で移して定義する。
★本ファイルは**その同型を経由せず**、古典的な定義から直接作る:

    `deg_F([L̄]) = ( log #(Γ(L)/Γ(X,⊤)·s) − Σ_σ log |s|_σ ) / [F:ℚ]`   （`s ≠ 0` は任意）

★★これで `MetricHeight.lean` の `SzpData`（[Szp] の同型を担ぐ仮定）が**不要になる**。

## ★★★★★★★機構

| 段 | 内容 | 出典 |
|---|---|---|
| 切断の存在 | `Γ(L)` は可逆なので非零元をもつ | `§9-779` ＋ `nontrivial_of_invertible` |
| ★**ノルムの非零性** | `s ≠ 0` ならどの複素点でもノルムは非零 | ★本ファイル |
| 切断に依らない | 段 D | `degArithPre_congr`（`§9-780`） |
| 加法性 | 段 E | `degArithPre_mul`（`§9-780`） |
| 等長不変 | | `degArithPre_of_isIsometry`（`§9-781`） |

★★★**ノルムの非零性**の証明が本ファイルの中身である:
`s ∈ Γ(L)` と `t ∈ Γ(L⁻¹)`（非零）を取ると、
`L ⊗ L⁻¹ ≅ Ō_X` の下で `s ⊗ t` は `𝓞_F` の**非零元** `c` に対応する
（`mk_injective_of_ne_zero`、`§9-778`）。
★★★★ノルムは掛け算になる（`AMetric.norm_mul`）ので

    `|s|_σ · |t|_σ = |c|_σ = |σ(c)| ≠ 0`

であり、したがって `|s|_σ ≠ 0` である。
★★★★★`σ` は体の埋め込みなので単射、`c ≠ 0` から `σ(c) ≠ 0` が出る。

## ★残っている段（明示）

★★`X(ℚ̄)` の上に降ろすには底変換 `deg_K(L|_{Spec 𝓞_K}) = deg_F(L)` の計量版が要る。
★★★因子の側ではそれは `degNormalized_baseChange` として証明済みである。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite NumberField
open ABC3.Found.GenEll
open scoped TensorProduct

/-! ## ★★★★可逆加群は非自明である -/

/-- ★**可逆加群は非自明である**（係数環が非自明なら）。

★`Mᵛ ⊗ M ≃ₗ R` なので、`M` が自明なら `R` も自明になってしまう。 -/
theorem nontrivial_of_invertible (R : Type) [CommRing R] [Nontrivial R] (M : Type)
    [AddCommGroup M] [Module R M] [Module.Invertible R M] : Nontrivial M := by
  rcases subsingleton_or_nontrivial M with h | h
  · exfalso
    haveI := h
    have e := Module.Invertible.linearEquiv R M
    haveI : Subsingleton (Module.Dual R M ⊗[R] M) := by infer_instance
    exact (not_subsingleton R) (e.symm.injective.subsingleton)
  · exact h

/-- ★★可逆層の大域切断には非零元がある。 -/
theorem exists_ne_zero_gammaPre {X : Scheme.{0}} [Nontrivial (Γ(X, (⊤ : X.Opens)) : Type)]
    (L : AInv X) : ∃ s : (L.carrier.sheaf.obj (op ⊤) : Type), s ≠ 0 := by
  haveI := invertible_gammaPre L
  haveI : Nontrivial ((L.carrier.sheaf.obj (op ⊤)) : Type) :=
    nontrivial_of_invertible (Γ(X, (⊤ : X.Opens)) : Type) ((L.carrier.sheaf.obj (op ⊤)) : Type)
  exact exists_ne 0

/-! ## ★★★★★★★★★★非零切断のノルムは非零である -/

/-- ★★★★★★★★★★**非零切断のノルムはどの複素点でも非零である**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★★機構: `t ∈ Γ(L⁻¹)` を非零に取ると、`L ⊗ L⁻¹ ≅ Ō_X` の下で
`s ⊗ t` は `𝓞_F` の**非零元** `c` に対応する（`mk_injective_of_ne_zero`）。
★★★★ノルムは掛け算になる（`AMetric.norm_mul`）ので

    `|s|_σ · |t|_σ = |c|_σ = |σ(c)| ≠ 0`

であり、したがって `|s|_σ ≠ 0` である。 -/
theorem norm_ne_zero_of_ne_zero (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) (hs : s ≠ 0)
    (σ : F →+* ℂ) :
    L.carrier.norm s (embSpecPoint F σ) ≠ 0 := by
  have htop : (embSpecPoint F σ) ⁻¹ᵁ
      (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens) = ⊤ := by simp
  obtain ⟨φ, hφ⟩ := L.isInv
  haveI := invertible_gammaPre L
  haveI hinvL : Module.Invertible
      (Γ(Spec (CommRingCat.of (𝓞 F)),
        (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)
      ((L.inv.sheaf.obj (op ⊤)) : Type) := invertible_gammaPre L.symm
  haveI : Nontrivial ((L.inv.sheaf.obj (op ⊤)) : Type) :=
    nontrivial_of_invertible
      (Γ(Spec (CommRingCat.of (𝓞 F)),
        (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)
      ((L.inv.sheaf.obj (op ⊤)) : Type)
  obtain ⟨t, ht⟩ := exists_ne (0 : ((L.inv.sheaf.obj (op ⊤)) : Type))
  have hinj := mk_injective_of_ne_zero
    (Γ(Spec (CommRingCat.of (𝓞 F)), (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)
    ((L.carrier.sheaf.obj (op ⊤)) : Type) ((L.inv.sheaf.obj (op ⊤)) : Type) s hs
  have hst : (s ⊗ₜ[(Γ(Spec (CommRingCat.of (𝓞 F)),
      (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)] t) ≠ 0 := by
    intro h
    exact ht (hinj (by simpa using h))
  have hcinj : Function.Injective (fun x => (φ.hom.app (op ⊤)) x) :=
    ((PresheafOfModules.evaluation (Spec (CommRingCat.of (𝓞 F))).ringCatSheaf.obj
      (op ⊤)).mapIso φ).toLinearEquiv.injective
  have hc : (φ.hom.app (op ⊤)) (s ⊗ₜ[(Γ(Spec (CommRingCat.of (𝓞 F)),
      (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)] t) ≠ 0 := by
    intro h
    exact hst (hcinj (h.trans (map_zero _).symm))
  have h1 := norm_of_isIsometry hφ (s ⊗ₜ[(Γ(Spec (CommRingCat.of (𝓞 F)),
      (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)] t) (embSpecPoint F σ)
  have h2 := AMetric.norm_mul L.carrier L.inv s t (embSpecPoint F σ)
  have h3 := AMetric.one_norm ((φ.hom.app (op ⊤)) (s ⊗ₜ[(Γ(Spec (CommRingCat.of (𝓞 F)),
      (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)] t)) (embSpecPoint F σ) htop
  have h4 : evalOn (embSpecPoint F σ) ⊤ htop ((φ.hom.app (op ⊤))
      (s ⊗ₜ[(Γ(Spec (CommRingCat.of (𝓞 F)),
        (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)] t)) ≠ 0 := by
    rw [evalOn_embSpecPoint F σ htop]
    intro h
    refine gammaSpec_hom_ne_zero (CommRingCat.of (𝓞 F)) _ hc ?_
    have hσ : (algebraMap (𝓞 F) F)
        ((Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).hom.hom
          ((φ.hom.app (op ⊤)) (s ⊗ₜ[(Γ(Spec (CommRingCat.of (𝓞 F)),
            (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)] t))) = 0 :=
      (map_eq_zero (f := σ)).mp h
    simpa [RingOfIntegers.coe_eq_zero_iff] using hσ
  intro hzero
  rw [hzero, zero_mul] at h2
  rw [h2] at h1
  rw [h3] at h1
  exact h4 (norm_eq_zero.mp h1)

/-! ## ★★★★★★★★★★`deg_F` を作る -/

/-- ★選んだ非零切断。 -/
noncomputable def chosenSec (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F)))) : (L.carrier.sheaf.obj (op ⊤) : Type) :=
  (exists_ne_zero_gammaPre L).choose

theorem chosenSec_ne_zero (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F)))) : chosenSec F L ≠ 0 :=
  (exists_ne_zero_gammaPre L).choose_spec

/-- ★★★★★★★★★★**`deg_F(L̄)`**（切断を選ばない形）。 -/
noncomputable def degAInv (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F)))) : ℝ :=
  degArithPre F L (chosenSec F L)

/-- ★★★**どの非零切断で測っても同じ**（段 D）。 -/
theorem degAInv_eq (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F))))
    (s : (L.carrier.sheaf.obj (op ⊤) : Type)) (hs : s ≠ 0) :
    degAInv F L = degArithPre F L s :=
  degArithPre_congr F L _ s (chosenSec_ne_zero F L) hs
    (fun σ => norm_ne_zero_of_ne_zero F L _ (chosenSec_ne_zero F L) σ)
    (fun σ => norm_ne_zero_of_ne_zero F L s hs σ)

/-- ★★★★★★★★**加法性**（段 E）。 -/
theorem degAInv_mul (F : Type) [Field F] [NumberField F]
    (L M : AInv (Spec (CommRingCat.of (𝓞 F)))) :
    degAInv F (L.mul M) = degAInv F L + degAInv F M := by
  haveI := invertible_gammaPre L
  haveI := invertible_gammaPre M
  have hs := chosenSec_ne_zero F L
  have ht := chosenSec_ne_zero F M
  have hinj := mk_injective_of_ne_zero
    (Γ(Spec (CommRingCat.of (𝓞 F)), (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)
    ((L.carrier.sheaf.obj (op ⊤)) : Type) ((M.carrier.sheaf.obj (op ⊤)) : Type)
    (chosenSec F L) hs
  have hst : ((chosenSec F L) ⊗ₜ[(Γ(Spec (CommRingCat.of (𝓞 F)),
      (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)] (chosenSec F M)) ≠ 0 := by
    intro h
    exact ht (hinj (by simpa using h))
  rw [degAInv_eq F (L.mul M) _ hst, degAInv_eq F L _ hs, degAInv_eq F M _ ht]
  exact degArithPre_mul F L M _ _ hs ht
    (fun σ => norm_ne_zero_of_ne_zero F L _ hs σ)
    (fun σ => norm_ne_zero_of_ne_zero F M _ ht σ)

/-- ★★★★★**等長同型類だけで決まる**。 -/
theorem degAInv_congr_isometric (F : Type) [Field F] [NumberField F]
    (L M : AInv (Spec (CommRingCat.of (𝓞 F))))
    (h : Isometric L.carrier M.carrier) : degAInv F L = degAInv F M := by
  obtain ⟨φ, hφ⟩ := h
  have hs := chosenSec_ne_zero F L
  have hinjφ : Function.Injective (fun x => (φ.hom.app (op ⊤)) x) :=
    ((PresheafOfModules.evaluation (Spec (CommRingCat.of (𝓞 F))).ringCatSheaf.obj
      (op ⊤)).mapIso φ).toLinearEquiv.injective
  have hne : (φ.hom.app (op ⊤)) (chosenSec F L) ≠ 0 := by
    intro hz
    exact hs (hinjφ (hz.trans (map_zero _).symm))
  rw [degAInv_eq F M _ hne, degArithPre_of_isIsometry F L M hφ (chosenSec F L)]
  rfl

/-- ★★★★★★★★★★**`deg_F : APic(Spec 𝓞_F) → ℝ`**——[Szp] の同型を経由しない形。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★★★これが `MetricHeight.lean` の `SzpData` を不要にする本体である。 -/
noncomputable def degAPicM (F : Type) [Field F] [NumberField F] :
    APicM (Spec (CommRingCat.of (𝓞 F))) → ℝ :=
  Quotient.lift (degAInv F) (fun L M h => degAInv_congr_isometric F L M h)

@[simp] theorem degAPicM_mk (F : Type) [Field F] [NumberField F]
    (L : AInv (Spec (CommRingCat.of (𝓞 F)))) :
    degAPicM F (APicM.mk L) = degAInv F L := rfl

/-- ★★★★★★★★★★**`deg_F` は準同型である**。 -/
theorem degAPicM_mul (F : Type) [Field F] [NumberField F]
    (a b : APicM (Spec (CommRingCat.of (𝓞 F)))) :
    degAPicM F (a * b) = degAPicM F a + degAPicM F b := by
  induction a using Quotient.ind with
  | _ L =>
    induction b using Quotient.ind with
    | _ M => exact degAInv_mul F L M

@[simp] theorem degAPicM_one (F : Type) [Field F] [NumberField F] :
    degAPicM F (1 : APicM (Spec (CommRingCat.of (𝓞 F)))) = 0 := by
  have h := degAPicM_mul F (1 : APicM (Spec (CommRingCat.of (𝓞 F)))) 1
  rw [one_mul] at h
  linarith

theorem degAPicM_inv (F : Type) [Field F] [NumberField F]
    (a : APicM (Spec (CommRingCat.of (𝓞 F)))) :
    degAPicM F a⁻¹ = - degAPicM F a := by
  have h := degAPicM_mul F a a⁻¹
  rw [mul_inv_cancel] at h
  rw [degAPicM_one] at h
  linarith

/-! ### ★出典の紐付け(`.src`) -/

def norm_ne_zero_of_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(非零切断のノルムはどの複素点でも非零)",
    sectionId := "genell-def-1-1-ii" }

def degAPicM.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(deg_F : APic(Spec 𝓞_F) → ℝ——[Szp] を経由しない構成)",
    sectionId := "genell-def-1-1-ii" }

def degAPicM_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(deg_F は準同型である)",
    sectionId := "genell-def-1-1-ii" }

def degAPicM_mul.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "degArithPre_mul(段 E、§9-780)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.degArithPre_mul") 4,
    .citation "[ABC3]" "degArithPre_congr(段 D、§9-780)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.degArithPre_congr") 4,
    .citation "[ABC3]" "degArithPre_of_isIsometry(等長不変、§9-781)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.degArithPre_of_isIsometry") 4,
    .implicitStep
      ("★★★原文は deg_F を ADiv(F)/APrc(F) ≅ APic(Spec 𝓞_F)([Szp] Prop 1.1)で移して" ++
       "定義するが、本構成はその同型を**経由しない**(逸脱の記録)。" ++
       "★消費側(Definition 1.2 の ht)が要るのは deg_F が準同型で底変換で不変なことだけなので、" ++
       "後続の証明に影響は出ない") 4,
    .implicitStep
      ("★★X(ℚ̄) の上に降ろすには底変換 deg_K(L|_{Spec 𝓞_K}) = deg_F(L) の計量版が要る。" ++
       "因子の側ではそれは degNormalized_baseChange として証明済み") 4 ]

end ABC3.Found.Arakelov
