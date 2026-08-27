/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.AMetricNorm
import Mathlib.NumberTheory.NumberField.InfinitePlace.Embeddings
import ABC3.Meta.Claim

/-!
# `deg_F` の**アルキメデス部分**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

## ★★★★★★★★台帳 `arakelov-degF-finite-places` の**段 C**

古典的な定義は

    `deg_F(L̄) = (1/[F:ℚ])·( log #(Γ(L)/𝒪_F·s) − Σ_σ log |s|(σ) )`

である（`s ≠ 0` は任意）。★本ファイルは**第 2 項**を作る。

## ★★★★★★複素点は埋め込みである

`Spec ℂ ⟶ Spec 𝒪_F` は環準同型 `𝒪_F → ℂ` に対応し、
`𝒪_F` の分数体が `F` なのでそれは埋め込み `σ : F ↪ ℂ` から来る。
★`Fintype.card (F →+* ℂ) = [F:ℚ]`（mathlib の `Embeddings.card`）が
正規化の分母をちょうど消す。

## ★★★★★加法性

    `archDeg (L̄ ⊗ M̄) (s ⊗ t) = archDeg L̄ s + archDeg M̄ t`

★機構は `AMetric.norm_mul`（`|s ⊗ t| = |s| · |t|`、2026-08-28）と `Real.log_mul`。
★★仮定は「`s`, `t` がどの埋め込みでも消えない」——
それは `s ≠ 0`（かつ `Γ(L)` が階数 1）から出るはずだが、
その段は**アフィンの橋（段 A）**を要するので本ファイルには含めない。

## ★残っている段（明示）

★★段 A（アフィンの橋）・段 B（有限部分 `log #(Γ(L)/𝒪_F·s)`）・
段 D（`s` の取り方に依らないこと＝積公式）・段 E（全体の加法性）。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite NumberField

variable (F : Type) [Field F] [NumberField F]

/-- ★**埋め込み `σ : F →+* ℂ` が定める `Spec 𝒪_F` の複素点**。 -/
noncomputable def embSpecPoint (σ : F →+* ℂ) :
    Spec (CommRingCat.of ℂ) ⟶ Spec (CommRingCat.of (𝓞 F)) :=
  Spec.map (CommRingCat.ofHom (σ.comp (algebraMap (𝓞 F) F)))

/-- ★★★★★★**`deg_F` のアルキメデス部分**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

    `archDeg L̄ s = −(Σ_σ log |s|(σ)) / [F:ℚ]` -/
noncomputable def archDeg (L : AMetric (Spec (CommRingCat.of (𝓞 F))))
    (s : L.sheaf.obj (op ⊤)) : ℝ :=
  -(∑ σ : F →+* ℂ, Real.log (L.norm s (embSpecPoint F σ))) / (Module.finrank ℚ F : ℝ)

/-- ★★★★★**アルキメデス部分は加法的である**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

★機構は `AMetric.norm_mul`（`|s ⊗ t| = |s| · |t|`）と `Real.log_mul`。
★★仮定は「`s`, `t` がどの埋め込みでも消えない」である。 -/
theorem archDeg_mul (L M : AMetric (Spec (CommRingCat.of (𝓞 F))))
    (s : L.sheaf.obj (op ⊤)) (t : M.sheaf.obj (op ⊤))
    (hs : ∀ σ : F →+* ℂ, L.norm s (embSpecPoint F σ) ≠ 0)
    (ht : ∀ σ : F →+* ℂ, M.norm t (embSpecPoint F σ) ≠ 0) :
    archDeg F (L * M)
        (s ⊗ₜ[(Γ(Spec (CommRingCat.of (𝓞 F)), (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens))
          : Type)] t)
      = archDeg F L s + archDeg F M t := by
  have hsum : ∀ σ : F →+* ℂ,
      Real.log ((L * M).norm (s ⊗ₜ[(Γ(Spec (CommRingCat.of (𝓞 F)),
          (⊤ : (Spec (CommRingCat.of (𝓞 F))).Opens)) : Type)] t) (embSpecPoint F σ))
        = Real.log (L.norm s (embSpecPoint F σ))
          + Real.log (M.norm t (embSpecPoint F σ)) :=
    fun σ => AMetric.log_norm_mul L M s t (embSpecPoint F σ) (hs σ) (ht σ)
  show -(∑ σ : F →+* ℂ, Real.log ((L * M).norm _ (embSpecPoint F σ))) / _ = _
  simp only [hsum, Finset.sum_add_distrib, archDeg]
  ring

/-- ★★**計量を `exp (-c)` 倍する類の話**——ノルムを一律 `r > 0` 倍すると
アルキメデス次数は `log r` だけ下がる。

★★★これが「`archDeg ≡ 0`」の退化を殺す欄である。 -/
theorem archDeg_of_norm_smul (L L' : AMetric (Spec (CommRingCat.of (𝓞 F))))
    (s : L.sheaf.obj (op ⊤)) (s' : L'.sheaf.obj (op ⊤)) (r : ℝ) (hr : 0 < r)
    (hnz : ∀ σ : F →+* ℂ, L.norm s (embSpecPoint F σ) ≠ 0)
    (h : ∀ σ : F →+* ℂ,
      L'.norm s' (embSpecPoint F σ) = r * L.norm s (embSpecPoint F σ)) :
    archDeg F L' s' = archDeg F L s - Real.log r := by
  have hcard : (Fintype.card (F →+* ℂ) : ℝ) = (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast congrArg (Nat.cast (R := ℝ)) (Embeddings.card F ℂ)
  have hpos : (Module.finrank ℚ F : ℝ) ≠ 0 := by
    have := Module.finrank_pos (R := ℚ) (M := F)
    positivity
  have hsum : ∀ σ : F →+* ℂ, Real.log (L'.norm s' (embSpecPoint F σ))
      = Real.log r + Real.log (L.norm s (embSpecPoint F σ)) :=
    fun σ => by rw [h σ, Real.log_mul hr.ne' (hnz σ)]
  show -(∑ σ : F →+* ℂ, Real.log (L'.norm s' (embSpecPoint F σ))) / _ = _
  simp only [hsum, Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ,
    nsmul_eq_mul, hcard, archDeg]
  field_simp
  ring

/-! ### ★出典の紐付け(`.src`) -/

def embSpecPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(埋め込みが定める Spec 𝒪_F の複素点)",
    sectionId := "genell-def-1-1-ii" }

def archDeg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(deg_F のアルキメデス部分。有限部分は含まない)",
    sectionId := "genell-def-1-1-ii" }

def archDeg_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4,
    item := "Definition 1.1, (ii)(アルキメデス部分の加法性。非消失を仮定した形)",
    sectionId := "genell-def-1-1-ii" }

def archDeg_mul.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "AMetric.norm_mul(|s ⊗ t| = |s| · |t|)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.AMetric.norm_mul") 4,
    .citation "[mathlib]" "NumberField.Embeddings.card(#(F →+* ℂ) = [F:ℚ])"
      (.inMathlib "NumberField.Embeddings.card") 4,
    .implicitStep
      ("★非消失の仮定は s ≠ 0(かつ Γ(L) が階数 1)から出るはずだが、" ++
       "その段は**アフィンの橋(台帳の段 A)**を要するので本ファイルには含めない") 4,
    .implicitStep
      ("★★残っている段: 段 A(アフィンの橋)・段 B(有限部分 log #(Γ(L)/𝒪_F·s))・" ++
       "段 D(s の取り方に依らないこと＝積公式)・段 E(全体の加法性)") 4 ]

end ABC3.Found.Arakelov
