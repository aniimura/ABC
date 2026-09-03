/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Profinite
import ABC3.Found.FrdI.Remark312
import ABC3.Found.FrdI.Prop113
import ABC3.Found.FrdI.Def45

/-!
# Thm62Slim —— `[FrdI] Theorem 6.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.FrdI

open CategoryTheory ABC3.Meta
universe v u w w2 v2 u2

/-! ## ★1. Frobenius-slim —— 原文 2 段目

原文 (FrdI p.112):
> finite, it follows formally that Z, ZG(H) are also residually finite, hence that D is
-/

/-- ★★★★★★**`Theorem 6.2, (iv)` 前半** ——
`Aut(𝒟_A → 𝒟)` がどれも副有限群 `G` に埋め込めるなら `𝒟` は **Frobenius-slim**。

★仮定 `h` が [Mzk7] `Corollary 1.1.6`(`Z_G(H) ≃ Aut(𝒟_{Spec L} → 𝒟)`)の
**使う分だけ**の形である —— 同型まで要らず、**単射準同型で足りる**。

★中身は
「副有限 ⟹ residually finite」(`Profinite.residuallyFinite_of_profinite`)
＋「単射で下りる」(`Profinite.residuallyFinite_of_injective`)
＋ `Remark 3.1.2`(`isFrobeniusSlim_of_residuallyFinite`)。 -/
theorem isFrobeniusSlim_of_embeds_profinite {E : Type u} [Category.{v} E]
    {G : Type w} [Group G] [TopologicalSpace G] [IsTopologicalGroup G] [CompactSpace G]
    [T2Space G] [TotallyDisconnectedSpace G]
    (h : ∀ A : E, ∃ f : Aut (Over.forget A) →* G, Function.Injective f) :
    IsFrobeniusSlim E := by
  haveI : Group.ResiduallyFinite G := Profinite.residuallyFinite_of_profinite
  refine isFrobeniusSlim_of_residuallyFinite E (fun A => ?_)
  obtain ⟨f, hf⟩ := h A
  exact Profinite.residuallyFinite_of_injective f hf

/-- ★★**同型版** —— `Z_G(H) ≃ Aut(𝒟_A → 𝒟)` の形をそのまま受ける。 -/
theorem isFrobeniusSlim_of_mulEquiv_subgroup {E : Type u} [Category.{v} E]
    {G : Type w} [Group G] [TopologicalSpace G] [IsTopologicalGroup G] [CompactSpace G]
    [T2Space G] [TotallyDisconnectedSpace G]
    (Hsub : E → Subgroup G) (e : ∀ A : E, Aut (Over.forget A) ≃* (Hsub A)) :
    IsFrobeniusSlim E :=
  isFrobeniusSlim_of_embeds_profinite
    (fun A => ⟨((Hsub A).subtype).comp ((e A) : Aut (Over.forget A) →* (Hsub A)),
      Subtype.val_injective.comp (e A).injective⟩)

/-! ## ★2. slim ⟺ `Z = {1}` —— 原文 3 段目の前半

原文 (FrdI p.112):
> the form "ZG(H)", it follows formally that D is slim if and only if Z = {1}, and
-/

/-- ★★★★★**`Theorem 6.2, (iv)` 中盤** —— `Aut(𝒟_A → 𝒟) ≃ Z_G(H_A)` のもとで
`𝒟` が slim ⟺ すべての `A` で `Z_G(H_A) = {1}`。

★原文の `Z = ⋃_H Z_G(H)` は、添字ごとに読むとこの形になる。 -/
theorem isSlimCat_iff_of_mulEquiv {E : Type u} [Category.{v} E]
    {G : Type w} [Group G] (Hsub : E → Subgroup G)
    (e : ∀ A : E, Aut (Over.forget A) ≃* (Hsub A)) :
    IsSlimCat E ↔ ∀ A : E, Hsub A = ⊥ := by
  constructor
  · intro hslim A
    rw [eq_bot_iff]
    intro z hz
    have h1 : (e A).symm ⟨z, hz⟩ = 1 := hslim A _
    have h2 : (⟨z, hz⟩ : (Hsub A)) = 1 := by
      have := congrArg (e A) h1
      simpa using this
    simpa using congrArg Subtype.val h2
  · intro hbot A
    -- ★`rw [hbot A]` は motive が壊れる(`Hsub A` が項の型に出る)ので、
    --   ★**書き換えず** `Hsub A ≤ ⊥` として当てる。
    have hle : Hsub A ≤ ⊥ := le_of_eq (hbot A)
    -- ★`IsSlimCat` の結論は `Iso.refl` の形、こちらは `Aut` の `1`。
    --   定義的に等しいので、`Aut` 側で作ってから `exact` で渡す。
    have key : ∀ η : Aut (Over.forget A), η = 1 := by
      intro η
      have h1 : ((e A) η : G) = 1 := Subgroup.mem_bot.mp (hle ((e A) η).2)
      have h2 : (e A) η = 1 := Subtype.ext (by simpa using h1)
      simpa using congrArg (e A).symm h2
    intro η
    exact key η

/-! ## ★3. Div-slim —— 原文 3 段目の後半

原文 (FrdI p.112):
> that D is Div-slim [relative to Φ] if and only if, for every 1 = z ∈ Z, there exists
-/

variable {D : Type u} [Category.{v} D] {Φ : MonoidOn.{v, u, w} D}

/-- ★★★★★**`Theorem 6.2, (iv)` 後半** —— `Aut(𝒟_A → 𝒟) ≃ Z_G(H_A)` のもとで
`𝒟` が `Φ` に関して Div-slim ⟺ `Z_G(H_A)` の側で見た作用が単射。

★原文の「for every `1 ≠ z ∈ Z`, there exists a finite Galois extension `L`
such that `z` acts nontrivially on `Φ(L)`」は、群準同型の核が自明という
言い方をしているが、我々は `IsDivSlim` の定義(単射性)のまま移送する。

★★**この同値は純粋な移送**である —— `(e A).symm` が全単射なので、
`overPhiAut Φ A` の単射性と、それを `Z_G(H_A)` に引き戻したものの単射性は同値。 -/
theorem isDivSlim_iff_of_mulEquiv {G : Type w2} [Group G] (Hsub : D → Subgroup G)
    (e : ∀ A : D, Aut (Over.forget A) ≃* (Hsub A)) :
    IsDivSlim Φ ↔ ∀ A : D, Function.Injective
      (fun z : (Hsub A) => overPhiAut Φ A ((e A).symm z)) := by
  constructor
  · intro hdiv A
    exact (hdiv A).comp (e A).symm.injective
  · intro hinj A η₁ η₂ hη
    have h := hinj A (a₁ := (e A) η₁) (a₂ := (e A) η₂) (by simpa using hη)
    simpa using congrArg (e A).symm h

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def isFrobeniusSlim_of_embeds_profinite.src : Source :=
  { paper := "FrdI", pdfPage := 112,
    item := "Theorem 6.2, (iv) — D は Frobenius-slim",
    sectionId := "frdi-thm-6-2" }

def isFrobeniusSlim_of_embeds_profinite.needs : List ProofObligation :=
  [ .otherPaper "[Mzk7]" "Corollary 1.1.6 — Z_G(H) ≃ Aut(𝒟_{Spec L} → 𝒟)。★仮定 h として受ける" 14,
    .citation "[ABC3]" "副有限群は Group.ResiduallyFinite"
      (.inProject "ABC3" "ABC3.Found.FrdI.Profinite.residuallyFinite_of_profinite") 112,
    .citation "[ABC3]" "単射準同型に沿って ResiduallyFinite は下りる"
      (.inProject "ABC3" "ABC3.Found.FrdI.Profinite.residuallyFinite_of_injective") 112,
    .citation "[ABC3]" "Remark 3.1.2 — Aut が residually finite なら Frobenius-slim"
      (.inProject "ABC3" "ABC3.Found.FrdI.isFrobeniusSlim_of_residuallyFinite") 112 ]

def isSlimCat_iff_of_mulEquiv.src : Source :=
  { paper := "FrdI", pdfPage := 112,
    item := "Theorem 6.2, (iv) — D is slim if and only if Z = {1}",
    sectionId := "frdi-thm-6-2" }

def isSlimCat_iff_of_mulEquiv.needs : List ProofObligation :=
  [ .otherPaper "[Mzk7]" "Corollary 1.1.6 — Z_G(H) ≃ Aut(𝒟_{Spec L} → 𝒟)。★仮定 e として受ける" 14,
    .citation "[ABC3]" "IsSlimCat"
      (.inProject "ABC3" "ABC3.Found.FrdI.IsSlimCat") 112,
    .derivation "Z = ⋃_H Z_G(H) を添字ごとに読み替える(原文の「Z is the union of subgroups of the form Z_G(H)」)" 112 ]

def isDivSlim_iff_of_mulEquiv.src : Source :=
  { paper := "FrdI", pdfPage := 112,
    item := "Theorem 6.2, (iv) — D is Div-slim iff z acts nontrivially on Φ(L)",
    sectionId := "frdi-thm-6-2" }

def isDivSlim_iff_of_mulEquiv.needs : List ProofObligation :=
  [ .otherPaper "[Mzk7]" "Corollary 1.1.6 — Z_G(H) ≃ Aut(𝒟_{Spec L} → 𝒟)。★仮定 e として受ける" 14,
    .citation "[ABC3]" "IsDivSlim"
      (.inProject "ABC3" "ABC3.Found.FrdI.IsDivSlim") 112,
    .derivation "全単射との合成は単射性を保つ(両向き)" 112 ]

end ABC3.Found.FrdI
