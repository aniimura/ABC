/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.NumberTheory.NumberField.InfinitePlace.Embeddings
import Mathlib.FieldTheory.Galois.Basic
import ABC3.Meta.Claim

/-!
# 埋め込みの**数え上げ** —— `σ : F ↪ ℂ` の延長はちょうど `[K:F]` 個（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

## ★★★★★★★★これは何か——底変換の**アルキメデス側**

`deg_K(L|_{Spec 𝓞_K}) = deg_F(L)` のアルキメデス側は、`§9-784` の

    `|pullSec f s|_τ = |s|_{τ|_F}`

と合わせると

    `Σ_{τ : K ↪ ℂ} h(τ|_F) = [K:F] · Σ_{σ : F ↪ ℂ} h(σ)`

に帰着する。★本ファイルはそれを入れる。

## ★★★★★機構（2 段）

    card_ext_eq         : `σ` の延長の個数は `[K:F]`（`AlgHom.card`）
    sum_over_extensions : ファイバーごとにまとめる（`Finset.sum_fiberwise`）

★`σ : F →+* ℂ` を使って `ℂ` を `F`-代数と見ると（`σ.toAlgebra`）、
`σ` の延長の集合は `K →ₐ[F] ℂ` と**同じもの**である。
★★mathlib の `AlgHom.card`（`#(E →ₐ[F] K) = [E:F]`、`K` は代数閉体）がそのまま効く。

## ★★★正規化との関係

`archDeg` は `1/[F:ℚ]` で割っている。`[K:ℚ] = [K:F]·[F:ℚ]` なので、
アルキメデス側は `[K:F]` 倍になり、正規化がそれを打ち消す。
-/

namespace ABC3.Found.Arakelov

open NumberField

/-- ★★★★**`σ` の延長は `F`-代数準同型と同じものである**。

★`σ : F →+* ℂ` で `ℂ` を `F`-代数と見ると、
`τ` が `σ` の延長であることと `τ` が `F`-代数準同型であることは同値。 -/
noncomputable def extEquivAlgHom (F K : Type) [Field F] [Field K] [Algebra F K]
    [Algebra F ℂ] (σ : F →+* ℂ) (hσ : algebraMap F ℂ = σ) :
    {τ : K →+* ℂ // τ.comp (algebraMap F K) = σ} ≃ (K →ₐ[F] ℂ) where
  toFun t := AlgHom.mk t.1 (fun x => by
      have h := congrFun (congrArg (fun (m : F →+* ℂ) => (m : F → ℂ)) t.2) x
      rw [hσ]
      exact h)
  invFun φ := ⟨φ.toRingHom, by
    ext x
    have h2 := φ.commutes x
    rw [hσ] at h2
    exact h2⟩
  left_inv t := by ext x; rfl
  right_inv φ := by ext x; rfl

/-- ★★★★★★**`σ` の延長はちょうど `[K:F]` 個**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

★mathlib の `AlgHom.card`（`#(E →ₐ[F] K) = [E:F]`）そのままである。 -/
theorem card_ext_eq (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] (σ : F →+* ℂ) :
    Nat.card {τ : K →+* ℂ // τ.comp (algebraMap F K) = σ} = Module.finrank F K := by
  letI : Algebra F ℂ := σ.toAlgebra
  have hσ : algebraMap F ℂ = σ := rfl
  rw [Nat.card_congr (extEquivAlgHom F K σ hσ), Nat.card_eq_fintype_card]
  exact AlgHom.card F K ℂ

/-- ★★★★★★★★★★**埋め込みの和はファイバーごとにまとまる**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) →X is any morphism that gives rise to x.

    `Σ_{τ : K ↪ ℂ} h(τ|_F) = [K:F] · Σ_{σ : F ↪ ℂ} h(σ)`

★★これが底変換の**アルキメデス側**の核である
——`§9-784` の `norm_pullSec` と合わせると `archDeg` が不変になる。 -/
theorem sum_over_extensions (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] (g : (F →+* ℂ) → ℝ) :
    ∑ τ : K →+* ℂ, g (τ.comp (algebraMap F K))
      = (Module.finrank F K : ℝ) * ∑ σ : F →+* ℂ, g σ := by
  classical
  rw [← Finset.sum_fiberwise Finset.univ (fun τ : K →+* ℂ => τ.comp (algebraMap F K))
       (fun τ => g (τ.comp (algebraMap F K)))]
  have key : ∀ σ : F →+* ℂ,
      (∑ τ ∈ Finset.univ.filter (fun τ : K →+* ℂ => τ.comp (algebraMap F K) = σ),
        g (τ.comp (algebraMap F K))) = (Module.finrank F K : ℝ) * g σ := by
    intro σ
    rw [Finset.sum_congr rfl (fun τ hτ => by rw [(Finset.mem_filter.mp hτ).2])]
    rw [Finset.sum_const, nsmul_eq_mul]
    congr 1
    have h1 : (Finset.univ.filter (fun τ : K →+* ℂ => τ.comp (algebraMap F K) = σ)).card
        = Fintype.card {τ : K →+* ℂ // τ.comp (algebraMap F K) = σ} :=
      (Fintype.card_subtype _).symm
    rw [h1, ← Nat.card_eq_fintype_card, card_ext_eq F K σ]
  rw [Finset.sum_congr rfl (fun σ _ => key σ), ← Finset.mul_sum]

/-! ### ★出典の紐付け(`.src`) -/

def card_ext_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(σ : F ↪ ℂ の延長はちょうど [K:F] 個)",
    sectionId := "genell-def-1-1-ii" }

def sum_over_extensions.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(埋め込みの和はファイバーごとにまとまる——底変換のアルキメデス側)",
    sectionId := "genell-def-1-1-ii" }

def sum_over_extensions.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "AlgHom.card(#(E →ₐ[F] K) = [E:F]、K は代数閉体)"
      (.inMathlib "AlgHom.card") 4,
    .citation "[mathlib]" "Finset.sum_fiberwise"
      (.inMathlib "Finset.sum_fiberwise") 4,
    .citation "[ABC3]" "norm_pullSec(|f^*s|(q) = |s|(q ≫ f)、§9-784)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.norm_pullSec") 3,
    .implicitStep
      ("★★残るのは配管である: embSpecPoint K τ ≫ Spec.map(algebraMap 𝓞_F 𝓞_K) " ++
       "= embSpecPoint F (τ|_F) を示し、§9-784 と本定理を組む") 4 ]

end ABC3.Found.Arakelov
