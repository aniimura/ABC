/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm62Slim
import ABC3.Found.FrdI.Prop55RlfRefl

/-!
# [FrdI] Theorem 6.4, (i) —— Div-slim は `Φ` から `Φ^pf`・`Φ^rlf` へ移る

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

原文 (FrdI p.114):
> and rationally standard type, but not of group-like type; D is Frobenius-slim

原文 (FrdI p.115):
> that D is Div-slim [relative to Φ, hence also relative to Φpf, Φrlf].

## ★★角括弧の「hence also」の中身

`IsDivSlim Φ` は **`Aut(𝒟_A → 𝒟) → Aut(𝒟_A → Mon)` が単射**という条件で、
`𝒞` を一切見ない(`𝒟` と `Φ` だけで決まる)。

★★`Φ^pf` / `Φ^rlf` へ移るのは、**単位射が単射**だからである:

| 移送先 | 単位射 | 単射性の根拠 |
|---|---|---|
| `Φ^pf` | `Pf.of : M → M^pf` | `Pf.of_injective_of_divisorial`(divisorial) |
| `Φ^rlf` | `toSc : M → ℝ≥0 ⊗_ℕ M` | `toSc_injective`(perf-factorial ＋ divisorial) |

★どちらも `Φ^X.map σ = Φ^X.map σ'` を単位射で引き戻して
`Φ.map σ = Φ.map σ'` に落とすだけである。

## ★★★原文の 3 文と、この節点の位置

原文 `Theorem 6.4, (i)` の第 2・3 文は 3 つを言う:

1. `𝒟` は **Frobenius-slim** —— [Mzk7] `Corollary 1.1.6` を入力に
   `isFrobeniusSlim_of_embeds_profinite`(在庫、`Thm62Slim.lean`)。
2. `𝒟` は `Φ`(したがって `Φ^pf`, `Φ^rlf`)に関して **Div-slim** ——
   ★`Φ` の側は「数体の自己同型で全素点を固定するものは恒等」
   (`eq_one_of_fixes_prime` ＋ 完全分解する素点の存在)、
   **`Φ^pf` / `Φ^rlf` への伝播が本ファイル**。
3. `𝒟` が **slim ⟺ `Z = {1}`** —— `isSlimCat_iff_of_mulEquiv`(在庫、`Thm62Slim.lean`)。

★1 と 3 は `Theorem 6.2, (iv)` と同じ一般形がそのまま使え、
原文自身も「entirely similar to the proof given for Theorem 6.2, (iv)」と書く。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal

universe v u w

/-! ## ★0. `Pf.map` は単位射と可換 -/

/-- ★`Pf.map f (Pf.of a) = Pf.of (f a)`。 -/
theorem Pf.map_of {M N : Type w} [AddCommMonoid M] [AddCommMonoid N] (f : M →+ N) (a : M) :
    Pf.map f (Pf.of a) = Pf.of (f a) := by
  rw [Pf.of_apply, Pf.of_apply, Pf.map_mk]

section DivSlimTransport

variable {D : Type u} [Category.{v} D] {Φ : MonoidOn.{v, u, w} D}

/-- ★★★★**Div-slim は `Φ` から `Φ^pf` へ移る**。

原文 (FrdI p.115):
> that D is Div-slim [relative to Φ, hence also relative to Φpf, Φrlf].

★中身は `Pf.of` が単射(`Φ` が divisorial)というだけである ——
`Pf.map (Φ.map σ) = Pf.map (Φ.map σ')` を `Pf.of a` に当てれば
`Pf.of (Φ.map σ a) = Pf.of (Φ.map σ' a)`、単射性で `Φ.map σ = Φ.map σ'`。 -/
theorem isDivSlim_pfOn (hsh : ∀ A : D, IsSharp (Φ.val A))
    (hdiv : ∀ A : D, IsDivisorial (Φ.val A))
    (h : IsDivSlim Φ) : IsDivSlim (Φ.pfOn hsh) := by
  intro A η₁ η₂ hη
  refine h A ?_
  apply Iso.ext
  apply NatTrans.ext
  funext X
  have hX : (overPhiAut (Φ.pfOn hsh) A η₁).hom.app X
      = (overPhiAut (Φ.pfOn hsh) A η₂).hom.app X := by rw [hη]
  have h1 : Pf.map (Φ.map (η₁.hom.app X.unop)) = Pf.map (Φ.map (η₂.hom.app X.unop)) :=
    congrArg AddCommMonCat.Hom.hom hX
  show Φ.functor.map ((η₁.hom.app X.unop).op) = Φ.functor.map ((η₂.hom.app X.unop).op)
  refine AddCommMonCat.ext (fun a => ?_)
  refine Pf.of_injective_of_divisorial (hdiv _) ?_
  have h2 := congrArg (fun t : Pf (Φ.val (X.unop).left) →+ Pf (Φ.val (X.unop).left) =>
    t (Pf.of a)) h1
  exact (Pf.map_of _ a).symm.trans (h2.trans (Pf.map_of _ a))

/-- ★★★★**Div-slim は `Φ` から `Φ^rlf` へ移る**。

★`Pf` のときと同じ形で、単位射が `toSc : M → ℝ≥0 ⊗_ℕ M` になるだけ。
その単射性は `toSc_injective`(perf-factorial ＋ divisorial、`Prop55RlfRefl.lean`)。 -/
theorem isDivSlim_phiScOn_nnreal
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := ℝ≥0) (Φ.map α)))
    (ι : ∀ A : D, Prime (Φ.val A) → Pf (Φ.val A) → ℝ≥0)
    (hpf : ∀ A : D, IsPerfFactorialWith (Φ.val A) (ι A))
    (hdiv : ∀ A : D, IsDivisorial (Φ.val A))
    (h : IsDivSlim Φ) : IsDivSlim (phiScOn ℝ≥0 Φ hcharInj) := by
  intro A η₁ η₂ hη
  refine h A ?_
  apply Iso.ext
  apply NatTrans.ext
  funext X
  have hX : (overPhiAut (phiScOn ℝ≥0 Φ hcharInj) A η₁).hom.app X
      = (overPhiAut (phiScOn ℝ≥0 Φ hcharInj) A η₂).hom.app X := by rw [hη]
  have h1 : scMap (S := ℝ≥0) (Φ.map (η₁.hom.app X.unop))
      = scMap (S := ℝ≥0) (Φ.map (η₂.hom.app X.unop)) :=
    congrArg AddCommMonCat.Hom.hom hX
  show Φ.functor.map ((η₁.hom.app X.unop).op) = Φ.functor.map ((η₂.hom.app X.unop).op)
  refine AddCommMonCat.ext (fun a => ?_)
  refine toSc_injective (hpf _) (hdiv _) ?_
  have h2 := congrArg (fun t : ScT ℝ≥0 (Φ.val (X.unop).left)
      →+ ScT ℝ≥0 (Φ.val (X.unop).left) => t (toSc a)) h1
  exact (scMap_toSc _ a).symm.trans (h2.trans (scMap_toSc _ a))

end DivSlimTransport

/-! ### ★出典の紐付け -/

def isDivSlim_pfOn.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — Div-slim は Φ から Φ^pf へ移る",
    sectionId := "frdi-thm-6-4" }

def isDivSlim_phiScOn_nnreal.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — Div-slim は Φ から Φ^rlf へ移る",
    sectionId := "frdi-thm-6-4" }

def isDivSlim_pfOn.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "Pf.of_injective_of_divisorial(単位射の単射性)"
      (.inProject "ABC3" "ABC3.Found.FrdI.Pf.of_injective_of_divisorial") 114,
    .citation "[ABC3]" "toSc_injective(実化の側の単位射の単射性)"
      (.inProject "ABC3" "ABC3.Found.FrdI.toSc_injective") 114,
    .derivation "原文の角括弧「hence also relative to Φ^pf, Φ^rlf」の中身" 115 ]

end ABC3.Found.FrdI
